"""
Wraps the vaccine discovery pipeline as a REST API.

Endpoints:
  POST /api/pipeline/run          - Start a pipeline run (async)
  GET  /api/pipeline/status/{id}  - Get run status and progress
  GET  /api/pipeline/results/{id} - Get completed results
  GET  /api/runs                  - List all runs
  GET  /api/health                - Health check
  WS   /ws/pipeline/{id}          - Real-time progress via WebSocket

Run:
  uvicorn api.main:app --reload --port 8000
"""

import requests as http_requests
from src.agents.predictors.coverage_agent import CoverageAgent
from src.agents.predictors.safety_filter import SafetyFilterAgent
from src.agents.predictors.bcell_predictor import BCellPredictorAgent
from src.agents.predictors.tcell_predictor import TCellPredictorAgent
from src.models.candidate import (
    CandidateProtein, CandidateStatus, PipelineRun,
    EpitopeResult, EpitopeType, ConfidenceTier,
)
import os
import sys
import uuid
import json
import time
import logging
import asyncio
from datetime import datetime
from typing import List, Dict, Any, Optional
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

from fastapi import FastAPI, HTTPException, BackgroundTasks, WebSocket, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("kozi.api")

# --APP ─────────────────────────────────────────────────────────────────

app = FastAPI(
    title="TOP_DEEP - Vaccine Discovery API",
    description="Vaccine target discovery pipeline",
    version="2.0.0",
)

@app.get("/api/test-db")
def test_db():
    from src.storage.supabase_client import db
    return {"connected": db.test_connection()}

@app.get("/")
def root():
    return {"message": "Welcome to the Kozi AI Vaccine Discovery API. Visit /docs for API documentation."}

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://localhost:3001",
                   "https://kozi-ai.com", "https://playground-kozi-ai.netlify.app", "https://playground-kozi-ai-old.netlify.app", "https://playground.kozi-ai.com"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Thread pool for running pipeline (CPU-bound IEDB calls)
executor = ThreadPoolExecutor(max_workers=2)

# In-memory run tracking (production: use Supabase)
active_runs: Dict[str, Dict] = {}


# --REQUEST/RESPONSE MODELS ────────────────────────────────────────────

class PipelineRequest(BaseModel):
    """Request to start a pipeline run."""
    input_type: str = Field(
        description="'pathogen', 'uniprot_id', or 'sequence'")
    input_value: str = Field(
        description="Pathogen name, UniProt ID, or amino acid sequence")
    protein_name: Optional[str] = Field(
        None, description="Optional protein name")
    run_safety: bool = Field(True, description="Run safety screening (N6)")
    run_coverage: bool = Field(
        True, description="Run population coverage (N7)")
    max_proteins: int = Field(
        3, description="Max proteins to analyze (for pathogen search)")


class PipelineStatus(BaseModel):
    run_id: str
    status: str  # "pending", "running", "completed", "failed"
    current_node: Optional[str] = None
    progress: float = 0.0  # 0.0 to 1.0
    message: Optional[str] = None
    started_at: Optional[str] = None
    completed_at: Optional[str] = None


class EpitopeResponse(BaseModel):
    sequence: str
    epitope_type: str
    hla_allele: Optional[str]
    ic50_nm: Optional[float]
    percentile_rank: Optional[float]
    confidence: str
    allergenicity_safe: Optional[bool]
    toxicity_safe: Optional[bool]


class CandidateResponse(BaseModel):
    protein_id: str
    protein_name: str
    sequence_length: int
    ctl_count: int
    ctl_strong: int
    htl_count: int
    bcell_count: int
    global_coverage_pct: float
    african_coverage_pct: float
    epitopes: List[EpitopeResponse]
    decisions: List[Dict[str, Any]]
    coverage_detail: Optional[Dict[str, Any]] = None


class PipelineResult(BaseModel):
    run_id: str
    status: str
    timing: Dict[str, float]
    candidates: List[CandidateResponse]


# --HELPER: FETCH FROM UNIPROT ──────────────────────────────────────────


def fetch_protein_by_id(uniprot_id: str) -> Optional[Dict]:
    try:
        resp = http_requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta", timeout=15)
        resp.raise_for_status()
        lines = resp.text.strip().split("\n")
        header = lines[0][1:]
        sequence = "".join(l.strip() for l in lines[1:])
        resp2 = http_requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json", timeout=15)
        name = header.split()[0]
        if resp2.status_code == 200:
            data = resp2.json()
            name = data.get("proteinDescription", {}).get(
                "recommendedName", {}).get("fullName", {}).get("value", name)
        return {"protein_id": uniprot_id, "protein_name": name, "sequence": sequence}
    except Exception as e:
        logger.error(f"UniProt fetch failed: {e}")
        return None


def search_pathogen_proteins(pathogen_name: str, max_results: int = 5) -> List[Dict]:
    try:
        query = f'(organism_name:"{pathogen_name}") AND (reviewed:true)'
        resp = http_requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={"query": query, "format": "json", "size": max_results * 3,
                    "fields": "accession,protein_name,organism_name,length,cc_subcellular_location"},
            timeout=30,
        )
        resp.raise_for_status()
        proteins = []
        for entry in resp.json().get("results", []):
            accession = entry.get("primaryAccession", "")
            name_obj = entry.get("proteinDescription", {}
                                 ).get("recommendedName", {})
            protein_name = name_obj.get("fullName", {}).get("value", "Unknown")
            location = ""
            for comment in entry.get("comments", []):
                if comment.get("commentType") == "SUBCELLULAR LOCATION":
                    for loc in comment.get("subcellularLocations", []):
                        location += loc.get("location", {}
                                            ).get("value", "") + " "
            is_surface = any(kw in location.lower() for kw in [
                "membrane", "secreted", "cell surface", "extracellular", "outer membrane", "cell wall"])
            proteins.append({
                "protein_id": accession, "protein_name": protein_name,
                "is_surface": is_surface, "location": location.strip(),
                "length": entry.get("sequence", {}).get("length", 0),
            })
        proteins.sort(key=lambda x: (not x["is_surface"],))
        return proteins[:max_results]
    except Exception as e:
        logger.error(f"UniProt search failed: {e}")
        return []


def fetch_sequence(protein_id: str) -> Optional[str]:
    try:
        resp = http_requests.get(
            f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta", timeout=15)
        resp.raise_for_status()
        lines = resp.text.strip().split("\n")
        return "".join(l.strip() for l in lines[1:])
    except Exception:
        return None


# --PIPELINE EXECUTION ──────────────────────────────────────────────────

def run_pipeline_sync(run_id: str, candidates: List[CandidateProtein], run_safety: bool, run_coverage: bool):
    """Run the full pipeline synchronously (called in background thread)."""
    run = active_runs[run_id]
    run["status"] = "running"
    run["started_at"] = datetime.now().isoformat()
    start = time.time()

    try:
        # N2: Antigen screening (VaxiJen + Phobius)
        run["current_node"] = "N2"
        run["message"] = "Screening antigens (VaxiJen, Phobius)..."
        run["progress"] = 0.05
        n2_start = time.time()
        try:
            from src.tools.vaxijen_client import vaxijen
            from src.tools.phobius_client import phobius
            for c in candidates:
                # VaxiJen antigenicity
                c.vaxijen_score = vaxijen.predict_antigenicity(c.sequence, "bacteria")
                # Phobius transmembrane prediction
                phobius_result = phobius.predict_transmembrane(c.sequence, c.protein_id)
                c.tmhmm_helices = phobius_result.get("num_tm_helices", 0)
                if phobius_result.get("has_signal_peptide"):
                    c.psortb_localization = "secreted"
                elif phobius_result.get("num_tm_helices", 0) == 0:
                    c.psortb_localization = "cytoplasmic"
                elif phobius_result.get("num_tm_helices", 0) == 1:
                    c.psortb_localization = "single_pass_membrane"
                else:
                    c.psortb_localization = "multi_pass_membrane"
                c.add_decision(
                    stage="antigen_screening",
                    decision="screened",
                    reasoning=f"VaxiJen={c.vaxijen_score:.2f}, TM_helices={c.tmhmm_helices}, localization={c.psortb_localization}",
                )
                logger.info(f"N2: {c.protein_name} — VaxiJen={c.vaxijen_score:.2f}, TM={c.tmhmm_helices}, loc={c.psortb_localization}")
        except Exception as e:
            logger.warning(f"N2 screening failed (continuing): {e}")
        n2_time = time.time() - n2_start

        # N3: T-cell
        run["current_node"] = "N3"
        run["message"] = "Predicting T-cell epitopes..."
        run["progress"] = 0.1
        n3_start = time.time()
        n3 = TCellPredictorAgent()
        candidates = n3.run(candidates)
        n3_time = time.time() - n3_start
        run["progress"] = 0.35

        # N4: B-cell
        run["current_node"] = "N4"
        run["message"] = "Predicting B-cell epitopes..."
        n4_start = time.time()
        n4 = BCellPredictorAgent()
        candidates = n4.run(candidates)
        n4_time = time.time() - n4_start
        run["progress"] = 0.55

        # N6: Safety
        n6_time = 0
        if run_safety:
            run["current_node"] = "N6"
            run["message"] = "Safety screening (AllerTOP, ToxinPred)..."
            n6_start = time.time()
            n6 = SafetyFilterAgent()
            candidates = n6.run(candidates)
            n6_time = time.time() - n6_start
        run["progress"] = 0.8

        # N7: Coverage
        n7_time = 0
        if run_coverage:
            run["current_node"] = "N7"
            run["message"] = "Calculating population coverage..."
            n7_start = time.time()
            n7 = CoverageAgent()
            candidates = n7.run(candidates)
            n7_time = time.time() - n7_start
        run["progress"] = 0.95

        # Save to Supabase
        try:
            from src.storage.supabase_client import db
            from src.models.candidate import PipelineRun as PipelineRunModel
            pipeline_run = PipelineRunModel(
                run_id=run_id,
                pathogen_name=candidates[0].protein_name if candidates else "unknown",
                input_type=run.get("input_type", "unknown"),
                raw_input=candidates[0].sequence[:100] if candidates else "",
                current_stage="completed",
                status="completed",
            )
            db.create_run(pipeline_run)
            db.update_run_stage(run_id, "completed", "completed")
            for candidate in candidates:
                db.save_candidate(run_id, candidate)
            logger.info(f"Saved to Supabase: run {run_id}")
        except Exception as e:
            logger.warning(f"Supabase save failed: {e}")

        total_time = time.time() - start

        # Build results
        run["status"] = "completed"
        run["progress"] = 1.0
        run["current_node"] = None
        run["message"] = "Complete"
        run["completed_at"] = datetime.now().isoformat()
        run["timing"] = {
            "total_seconds": round(total_time, 1),
            "n2_screening": round(n2_time, 1),
            "n3_tcell": round(n3_time, 1),
            "n4_bcell": round(n4_time, 1),
            "n6_safety": round(n6_time, 1),
            "n7_coverage": round(n7_time, 1),
        }
        run["candidates"] = candidates

    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        run["status"] = "failed"
        run["message"] = str(e)
        run["progress"] = 0


def candidate_to_response(c: CandidateProtein) -> CandidateResponse:
    """Convert CandidateProtein to API response."""
    epitopes = []
    for ep in c.ctl_epitopes:
        epitopes.append(EpitopeResponse(
            sequence=ep.sequence, epitope_type="CTL", hla_allele=ep.hla_allele,
            ic50_nm=ep.ic50_nm, percentile_rank=ep.percentile_rank,
            confidence=ep.confidence_tier.value,
            allergenicity_safe=ep.allergenicity_safe, toxicity_safe=ep.toxicity_safe,
        ))
    for ep in c.htl_epitopes:
        epitopes.append(EpitopeResponse(
            sequence=ep.sequence, epitope_type="HTL", hla_allele=ep.hla_allele,
            ic50_nm=ep.ic50_nm, percentile_rank=ep.percentile_rank,
            confidence=ep.confidence_tier.value,
            allergenicity_safe=ep.allergenicity_safe, toxicity_safe=ep.toxicity_safe,
        ))
    for ep in c.bcell_epitopes:
        epitopes.append(EpitopeResponse(
            sequence=ep.sequence, epitope_type="B-cell", hla_allele=None,
            ic50_nm=None, percentile_rank=None,
            confidence=ep.confidence_tier.value,
            allergenicity_safe=ep.allergenicity_safe, toxicity_safe=ep.toxicity_safe,
        ))

    # Extract coverage detail from decisions
    coverage_detail = None
    for d in c.decisions:
        if d.get("stage") == "coverage_analysis" and "per_population" in d:
            coverage_detail = d["per_population"]

    return CandidateResponse(
        protein_id=c.protein_id,
        protein_name=c.protein_name,
        sequence_length=len(c.sequence),
        ctl_count=len(c.ctl_epitopes),
        ctl_strong=len(
            [e for e in c.ctl_epitopes if e.confidence_tier == ConfidenceTier.HIGH]),
        htl_count=len(c.htl_epitopes),
        bcell_count=len(c.bcell_epitopes),
        global_coverage_pct=round((c.hla_coverage_global or 0) * 100, 1),
        african_coverage_pct=round((c.hla_coverage_africa or 0) * 100, 1),
        epitopes=epitopes,
        decisions=c.decisions,
        coverage_detail=coverage_detail,
    )


# --API ENDPOINTS ───────────────────────────────────────────────────────

@app.get("/api/health")
async def health():
    return {"status": "ok", "service": "kozi-pipeline", "version": "2.0.0"}


@app.post("/api/pipeline/run", response_model=PipelineStatus)
async def start_pipeline(req: PipelineRequest, background_tasks: BackgroundTasks):
    """Start a new pipeline run. Returns immediately with run_id for polling."""
    run_id = str(uuid.uuid4())

    # Build candidate list based on input type
    candidates = []

    if req.input_type == "uniprot_id":
        prot = fetch_protein_by_id(req.input_value)
        if not prot:
            raise HTTPException(
                status_code=404, detail=f"Protein {req.input_value} not found on UniProt")
        candidates.append(CandidateProtein(
            protein_id=prot["protein_id"],
            protein_name=prot["protein_name"],
            sequence=prot["sequence"],
            source="uniprot",
            stage="antigen_screening",
            status=CandidateStatus.ACTIVE,
        ))

    elif req.input_type == "sequence":
        seq = req.input_value.upper().replace(" ", "").replace("\n", "")
        if len(seq) < 10:
            raise HTTPException(
                status_code=400, detail="Sequence too short (minimum 10 amino acids)")
        candidates.append(CandidateProtein(
            protein_id="user_input",
            protein_name=req.protein_name or "Custom protein",
            sequence=seq,
            source="user_input",
            stage="antigen_screening",
            status=CandidateStatus.ACTIVE,
        ))

    elif req.input_type == "pathogen":
        proteins = search_pathogen_proteins(
            req.input_value, max_results=req.max_proteins)
        if not proteins:
            raise HTTPException(
                status_code=404, detail=f"No proteins found for '{req.input_value}'")
        for p in proteins:
            seq = fetch_sequence(p["protein_id"])
            if seq and len(seq) >= 20:
                candidates.append(CandidateProtein(
                    protein_id=p["protein_id"],
                    protein_name=f"{p['protein_name']} ({p['protein_id']})",
                    sequence=seq,
                    source="uniprot",
                    stage="antigen_screening",
                    status=CandidateStatus.ACTIVE,
                ))
    else:
        raise HTTPException(
            status_code=400, detail=f"Invalid input_type: {req.input_type}")

    if not candidates:
        raise HTTPException(
            status_code=400, detail="No valid proteins to analyze")

    # Register run
    active_runs[run_id] = {
        "run_id": run_id,
        "status": "pending",
        "current_node": None,
        "progress": 0.0,
        "message": f"Starting analysis of {len(candidates)} protein(s)...",
        "input_type": req.input_type,
        "input_value": req.input_value,
        "protein_count": len(candidates),
        "started_at": None,
        "completed_at": None,
        "candidates": None,
        "timing": None,
    }

    # Run pipeline in background thread
    executor.submit(run_pipeline_sync, run_id, candidates,
                    req.run_safety, req.run_coverage)

    return PipelineStatus(
        run_id=run_id,
        status="pending",
        progress=0.0,
        message=f"Queued: {len(candidates)} protein(s) to analyze",
    )


@app.get("/api/pipeline/status/{run_id}", response_model=PipelineStatus)
async def get_pipeline_status(run_id: str):
    """Poll for pipeline progress."""
    run = active_runs.get(run_id)
    if not run:
        raise HTTPException(status_code=404, detail="Run not found")

    return PipelineStatus(
        run_id=run["run_id"],
        status=run["status"],
        current_node=run.get("current_node"),
        progress=run.get("progress", 0),
        message=run.get("message"),
        started_at=run.get("started_at"),
        completed_at=run.get("completed_at"),
    )


@app.get("/api/pipeline/results/{run_id}", response_model=PipelineResult)
async def get_pipeline_results(run_id: str):
    """Get completed pipeline results."""
    run = active_runs.get(run_id)
    if not run:
        raise HTTPException(status_code=404, detail="Run not found")

    if run["status"] != "completed":
        raise HTTPException(
            status_code=409, detail=f"Run is {run['status']}, not completed")

    candidates_response = []
    for c in run.get("candidates", []):
        candidates_response.append(candidate_to_response(c))

    return PipelineResult(
        run_id=run_id,
        status="completed",
        timing=run.get("timing", {}),
        candidates=candidates_response,
    )


@app.get("/api/runs")
async def list_runs():
    """List all runs (most recent first)."""
    runs = []
    for run_id, run in active_runs.items():
        runs.append({
            "run_id": run_id,
            "status": run["status"],
            "input_type": run.get("input_type"),
            "input_value": run.get("input_value"),
            "protein_count": run.get("protein_count"),
            "progress": run.get("progress", 0),
            "started_at": run.get("started_at"),
            "completed_at": run.get("completed_at"),
        })
    return sorted(runs, key=lambda x: x.get("started_at") or "", reverse=True)


# --WEBSOCKET FOR REAL-TIME PROGRESS

@app.websocket("/ws/pipeline/{run_id}")
async def pipeline_websocket(websocket: WebSocket, run_id: str):
    """WebSocket for real-time pipeline progress updates."""
    await websocket.accept()

    try:
        last_progress = -1
        while True:
            run = active_runs.get(run_id)
            if not run:
                await websocket.send_json({"error": "Run not found"})
                break

            progress = run.get("progress", 0)
            if progress != last_progress:
                await websocket.send_json({
                    "run_id": run_id,
                    "status": run["status"],
                    "current_node": run.get("current_node"),
                    "progress": progress,
                    "message": run.get("message"),
                })
                last_progress = progress

            if run["status"] in ("completed", "failed"):
                # Send final update
                await websocket.send_json({
                    "run_id": run_id,
                    "status": run["status"],
                    "progress": 1.0 if run["status"] == "completed" else progress,
                    "message": run.get("message"),
                    "timing": run.get("timing"),
                })
                break

            await asyncio.sleep(1)  # Poll every second

    except WebSocketDisconnect:
        logger.info(f"WebSocket disconnected for run {run_id}")
    except Exception as e:
        logger.error(f"WebSocket error: {e}")



# --STARTUP

@app.on_event("startup")
async def startup():
    logger.info("Kozi AI API starting...")
    logger.info(
        "Endpoints: POST /api/pipeline/run, GET /api/pipeline/status/{id}, GET /api/pipeline/results/{id}")
