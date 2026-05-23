"""
Wraps the vaccine discovery pipeline as a REST API.

Endpoints:
  POST /api/pipeline/run          - Start a pipeline run (requires auth)
  GET  /api/pipeline/status/{id}  - Get run status and progress (requires auth)
  GET  /api/pipeline/results/{id} - Get completed results — falls back to Supabase
  GET  /api/runs                  - List runs for authenticated user (requires auth)
  GET  /api/health                - Extended health: tool status, N7 method (public)
  WS   /ws/pipeline/{id}          - Real-time progress (token query param)

Pipeline nodes:
  N1  Data curation    provenance audit, protein load record
  N2  Antigen screen   VaxiJen 2.0 + Phobius (organism-aware VaxiJen)
  N3  T-cell predict   NetMHCpan 4.1 / MHCflurry 2.0 fallback
  N4  B-cell predict   IEDB BepiPred 2.0
  N5  Structure        AlphaFold DB
  N6  Safety filter    AllerTOP + AllergenFP + ToxinPred + BLAST (unscored != safe)
  N7  Coverage         IEDB tool v3.0.1 / AFND 2020 fallback
  N8  Construct        ProtParam (Biopython), RS09 adjuvant, linker assembly
"""

import requests as http_requests
from src.agents.predictors.coverage_agent import CoverageAgent
from src.agents.predictors.safety_filter import SafetyFilterAgent, safety_filter
from src.agents.predictors.bcell_predictor import BCellPredictorAgent
from src.agents.predictors.tcell_predictor import TCellPredictorAgent
from src.agents.predictors.structure_agent import StructureAgent
from src.agents.predictors.construct_designer import ConstructDesignerAgent
from src.models.candidate import (
    CandidateProtein, CandidateStatus, PipelineRun,
    EpitopeResult, EpitopeType, ConfidenceTier,
)
import os
import sys
import uuid
import time
import logging
import asyncio
from datetime import datetime
from typing import List, Dict, Any, Optional
from concurrent.futures import ThreadPoolExecutor

from fastapi import FastAPI, HTTPException, BackgroundTasks, WebSocket, WebSocketDisconnect, Depends, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from api.auth import require_user, optional_user, UserClaims

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("kozi.api")

# -- APP

app = FastAPI(
    title="TOP_DEEP — Vaccine Discovery API",
    description="TOPE_DEEP epitope-based vaccine target discovery pipeline",
    version="2.1.0",
)

@app.get("/api/test-db")
def test_db():
    from src.storage.supabase_client import db
    return {"connected": db.test_connection()}

@app.get("/")
def root():
    return {"message": "Kozi AI TOPE_DEEP API. See /docs."}

_CORS_ORIGINS: list[str] = [
    o.strip()
    for o in os.getenv("CORS_ORIGINS", "http://localhost:3000").split(",")
    if o.strip()
]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

executor = ThreadPoolExecutor(max_workers=2)
active_runs: Dict[str, Dict] = {}


# -- MODELS

class PipelineRequest(BaseModel):
    input_type: str  = Field(description="'pathogen', 'uniprot_id', or 'sequence'")
    input_value: str = Field(description="Pathogen name, UniProt ID, or amino acid sequence")
    protein_name: Optional[str] = Field(None)
    run_safety:   bool = Field(True)
    run_coverage: bool = Field(True)
    max_proteins: int  = Field(3)

class PipelineStatus(BaseModel):
    run_id: str
    status: str
    current_node: Optional[str] = None
    progress: float = 0.0
    message: Optional[str] = None
    started_at: Optional[str] = None
    completed_at: Optional[str] = None

class EpitopeResponse(BaseModel):
    sequence: str
    epitope_type: str
    hla_allele: Optional[str]
    ic50_nm: Optional[float]
    ic50_note: Optional[str] = "approximated_from_percentile_rank"
    percentile_rank: Optional[float]
    confidence: str
    allergenicity_safe: Optional[bool]
    toxicity_safe: Optional[bool]
    safety_verdict: Optional[str] = None
    method_used: Optional[str] = None

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
    structure_source: Optional[str] = None
    structure_pdb_url: Optional[str] = None
    phobius_localization: Optional[str] = None
    vaxijen_score: Optional[float] = None
    vaxijen_method: Optional[str] = None
    epitopes: List[EpitopeResponse]
    decisions: List[Dict[str, Any]]
    coverage_detail: Optional[Dict[str, Any]] = None

class PipelineResult(BaseModel):
    run_id: str
    status: str
    timing: Dict[str, float]
    candidates: List[CandidateResponse]
    construct_report: Optional[Dict[str, Any]] = None


# -- UNIPROT HELPERS

def _infer_organism_class(input_value: str) -> str:
    v = input_value.lower()
    if any(w in v for w in ["virus","viral","sars","covid","influenza","hiv",
                              "ebola","dengue","zika","rabies","hepatitis","rsv","hpv"]):
        return "virus"
    if any(w in v for w in ["plasmodium","malaria","leishmania","trypanosoma",
                              "toxoplasma","schistosoma","parasite","helminth"]):
        return "parasite"
    return "bacteria"

def fetch_protein_by_id(uniprot_id: str) -> Optional[Dict]:
    try:
        r = http_requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta", timeout=15)
        r.raise_for_status()
        lines = r.text.strip().split("\n")
        seq = "".join(l.strip() for l in lines[1:])
        r2 = http_requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json", timeout=15)
        name = lines[0][1:].split()[0]
        if r2.status_code == 200:
            d = r2.json()
            name = d.get("proteinDescription",{}).get("recommendedName",{}).get(
                "fullName",{}).get("value", name)
        return {"protein_id": uniprot_id, "protein_name": name, "sequence": seq}
    except Exception as e:
        logger.error(f"UniProt fetch failed: {e}")
        return None

def search_pathogen_proteins(pathogen_name: str, max_results: int = 5) -> List[Dict]:
    try:
        r = http_requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={
                "query": f'(organism_name:"{pathogen_name}") AND (reviewed:true)',
                "format": "json", "size": max_results * 3,
                "fields": "accession,protein_name,organism_name,length,cc_subcellular_location",
            }, timeout=30,
        )
        r.raise_for_status()
        proteins = []
        for entry in r.json().get("results", []):
            acc = entry.get("primaryAccession","")
            name = entry.get("proteinDescription",{}).get("recommendedName",{}).get(
                "fullName",{}).get("value","Unknown")
            location = ""
            for comment in entry.get("comments",[]):
                if comment.get("commentType") == "SUBCELLULAR LOCATION":
                    for loc in comment.get("subcellularLocations",[]):
                        location += loc.get("location",{}).get("value","") + " "
            is_surface = any(kw in location.lower() for kw in [
                "membrane","secreted","cell surface","extracellular","outer membrane","cell wall"])
            proteins.append({
                "protein_id": acc, "protein_name": name,
                "is_surface": is_surface, "location": location.strip(),
                "length": entry.get("sequence",{}).get("length",0),
            })
        proteins.sort(key=lambda x: (not x["is_surface"],))
        return proteins[:max_results]
    except Exception as e:
        logger.error(f"UniProt search failed: {e}")
        return []

def fetch_sequence(protein_id: str) -> Optional[str]:
    try:
        r = http_requests.get(
            f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta", timeout=15)
        r.raise_for_status()
        lines = r.text.strip().split("\n")
        return "".join(l.strip() for l in lines[1:])
    except Exception:
        return None


# -- SUPABASE RECONSTRUCTION

def reconstruct_from_supabase(run_id: str, user_id: str) -> Optional[Dict]:
    try:
        from src.storage.supabase_client import db
        run_row = db.get_run(run_id)
        if not run_row:
            return None
        if run_row.get("user_id") and run_row["user_id"] != user_id:
            return None

        candidates_rows = db.get_candidates_for_run(run_id)
        if not candidates_rows:
            return None

        candidates_response = []
        for c_row in candidates_rows:
            ep_rows = db.get_epitopes_for_candidate(c_row["id"])
            epitopes = [
                EpitopeResponse(
                    sequence=ep["sequence"],
                    epitope_type=ep["epitope_type"],
                    hla_allele=ep.get("hla_allele"),
                    ic50_nm=ep.get("ic50_nm"),
                    ic50_note="approximated_from_percentile_rank",
                    percentile_rank=ep.get("percentile_rank"),
                    confidence=ep.get("confidence_tier","uncertain"),
                    allergenicity_safe=ep.get("allergenicity_safe"),
                    toxicity_safe=ep.get("toxicity_safe"),
                    safety_verdict=(ep.get("tool_outputs") or {}).get("safety_verdict"),
                    method_used=(ep.get("tool_outputs") or {}).get("method_used"),
                )
                for ep in ep_rows
            ]
            coverage_detail = None
            for d in (c_row.get("decisions") or []):
                if d.get("stage") == "coverage_analysis" and "per_population" in d:
                    coverage_detail = d["per_population"]

            ctl   = [e for e in epitopes if e.epitope_type == "CTL"]
            htl   = [e for e in epitopes if e.epitope_type == "HTL"]
            bcell = [e for e in epitopes if "B-cell" in e.epitope_type]

            # Extract vaxijen info from decisions
            vaxijen_score = c_row.get("vaxijen_score")
            vaxijen_method = None
            phobius_loc = c_row.get("structure_source")
            for d in (c_row.get("decisions") or []):
                if d.get("stage") == "antigen_screening":
                    vaxijen_method = d.get("vaxijen_method")
                    phobius_loc = d.get("phobius_localization") or phobius_loc

            candidates_response.append(CandidateResponse(
                protein_id=c_row["protein_id"],
                protein_name=c_row["protein_name"],
                sequence_length=len(c_row.get("sequence") or ""),
                ctl_count=len(ctl),
                ctl_strong=len([e for e in ctl if e.confidence == "high"]),
                htl_count=len(htl),
                bcell_count=len(bcell),
                global_coverage_pct=round((c_row.get("hla_coverage_global") or 0)*100, 1),
                african_coverage_pct=round((c_row.get("hla_coverage_africa") or 0)*100, 1),
                structure_source=c_row.get("structure_source"),
                structure_pdb_url=c_row.get("structure_pdb_path"),
                phobius_localization=phobius_loc,
                vaxijen_score=vaxijen_score,
                vaxijen_method=vaxijen_method,
                epitopes=epitopes,
                decisions=c_row.get("decisions") or [],
                coverage_detail=coverage_detail,
            ))

        return {
            "run_id": run_id,
            "status": run_row.get("status","completed"),
            "timing": {"total_seconds": 0},
            "candidates": candidates_response,
            "construct_report": None,
        }
    except Exception as e:
        logger.error(f"Supabase reconstruction failed for {run_id}: {e}")
        return None


# -- PIPELINE EXECUTION

def run_pipeline_sync(
    run_id: str,
    candidates: List[CandidateProtein],
    run_safety: bool,
    run_coverage: bool,
    organism_class: str = "bacteria",
):
    run = active_runs[run_id]
    run["status"] = "running"
    run["started_at"] = datetime.now().isoformat()
    start = time.time()
    n1_time = n2_time = n3_time = n4_time = n5_time = 0
    n6_time = n7_time = n8_time = 0

    try:
        # N1: Data curation audit
        run["current_node"] = "N1"
        run["message"] = "Recording data provenance..."
        run["progress"] = 0.02
        n1_start = time.time()
        for c in candidates:
            c.add_decision(
                stage="data_curation",
                decision="protein_loaded",
                reasoning=(
                    f"Protein loaded from {c.source}. "
                    f"UniProt ID: {c.protein_id}. "
                    f"Sequence length: {len(c.sequence)} aa. "
                    f"Input type: {run.get('input_type','unknown')}. "
                    f"Organism class inferred: {organism_class}."
                ),
                source=c.source,
                sequence_length=len(c.sequence),
                input_type=run.get("input_type"),
                organism_class=organism_class,
            )
            logger.info(
                f"N1: {c.protein_name} ({c.protein_id}) "
                f"loaded from {c.source}, {len(c.sequence)} aa"
            )
        n1_time = time.time() - n1_start
        run["progress"] = 0.04

        # N2: Antigen screening
        run["current_node"] = "N2"
        run["message"] = "Screening antigens (VaxiJen, Phobius)..."
        run["progress"] = 0.05
        n2_start = time.time()
        try:
            from src.tools.vaxijen_client import vaxijen
            from src.tools.phobius_client import phobius
            for c in candidates:
                # VaxiJen — organism-aware
                score = vaxijen.predict_antigenicity(c.sequence, organism_class)
                c.vaxijen_score = score
                # Detect if fallback was used (real VaxiJen returns score from server)
                vaxijen_method = (
                    "VaxiJen_2.0_real"
                    if vaxijen.is_server_available()
                    else "VaxiJen_2.0_ACC_fallback"
                )

                phobius_result = phobius.predict_transmembrane(c.sequence, c.protein_id)
                tm_helices = phobius_result.get("num_tm_helices", 0)
                c.tmhmm_helices = tm_helices

                # Use phobius_localization — not psortb (PSORTb not available on Railway)
                if phobius_result.get("has_signal_peptide"):
                    c.psortb_localization = "secreted"
                elif tm_helices == 0:
                    c.psortb_localization = "cytoplasmic"
                elif tm_helices == 1:
                    c.psortb_localization = "single_pass_membrane"
                else:
                    c.psortb_localization = "multi_pass_membrane"

                c.add_decision(
                    stage="antigen_screening",
                    decision="screened",
                    reasoning=(
                        f"VaxiJen={c.vaxijen_score:.3f} "
                        f"[organism={organism_class}, method={vaxijen_method}]. "
                        f"Phobius TM_helices={tm_helices}, "
                        f"localization={c.psortb_localization}. "
                        f"Note: localization predicted by Phobius "
                        f"(Kall et al. 2004), not PSORTb. "
                        f"PSORTb requires Docker-in-Docker, unavailable on Railway."
                    ),
                    vaxijen_score=c.vaxijen_score,
                    vaxijen_method=vaxijen_method,
                    organism_class=organism_class,
                    tm_helices=tm_helices,
                    phobius_localization=c.psortb_localization,
                    localization_tool="Phobius_2.0",
                )
                logger.info(
                    f"N2: {c.protein_name} VaxiJen={c.vaxijen_score:.3f} "
                    f"[{vaxijen_method}] TM={tm_helices} loc={c.psortb_localization}"
                )
        except Exception as e:
            logger.warning(f"N2 screening failed (continuing): {e}")
        n2_time = time.time() - n2_start

        # N3: T-cell
        run["current_node"] = "N3"
        run["message"] = "Predicting T-cell epitopes (NetMHCpan, NetMHCIIpan)..."
        run["progress"] = 0.10
        n3_start = time.time()
        n3 = TCellPredictorAgent()
        candidates = n3.run(candidates)
        n3_time = time.time() - n3_start
        run["progress"] = 0.35

        # N4: B-cell
        run["current_node"] = "N4"
        run["message"] = "Predicting B-cell epitopes (IEDB BepiPred 2.0)..."
        n4_start = time.time()
        n4 = BCellPredictorAgent()
        candidates = n4.run(candidates)
        n4_time = time.time() - n4_start
        run["progress"] = 0.55

        # N5: Structure
        run["current_node"] = "N5"
        run["message"] = "Retrieving AlphaFold structures..."
        n5_start = time.time()
        n5 = StructureAgent()
        candidates = n5.run(candidates)
        n5_time = time.time() - n5_start
        run["progress"] = 0.65

        # N6: Safety
        if run_safety:
            run["current_node"] = "N6"
            run["message"] = "Safety screening (AllerTOP, AllergenFP, ToxinPred, BLAST)..."
            n6_start = time.time()
            n6 = SafetyFilterAgent()
            candidates = n6.run(candidates)
            n6_time = time.time() - n6_start
        run["progress"] = 0.80

        # N7: Coverage
        if run_coverage:
            run["current_node"] = "N7"
            run["message"] = "Calculating population coverage (IEDB / AFND 2020)..."
            n7_start = time.time()
            n7 = CoverageAgent()
            candidates = n7.run(candidates)
            n7_time = time.time() - n7_start
        run["progress"] = 0.92

        # N8: Construct
        run["current_node"] = "N8"
        run["message"] = "Assembling multi-epitope construct (ProtParam, RS09 adjuvant)..."
        n8_start = time.time()
        n8 = ConstructDesignerAgent()
        candidates, construct_report = n8.run(candidates)
        n8_time = time.time() - n8_start
        run["construct_report"] = construct_report
        run["progress"] = 0.98

        # Save to Supabase
        try:
            from src.storage.supabase_client import db
            from src.models.candidate import PipelineRun as PipelineRunModel
            user_id = run.get("user_id")
            pipeline_run = PipelineRunModel(
                run_id=run_id,
                pathogen_name=candidates[0].protein_name if candidates else "unknown",
                input_type=run.get("input_type","unknown"),
                raw_input=candidates[0].sequence[:100] if candidates else "",
                current_stage="completed",
                status="completed",
            )
            db.create_run(pipeline_run, user_id=user_id)
            db.update_run_stage(run_id, "completed", "completed")
            for candidate in candidates:
                db.save_candidate(run_id, candidate, user_id=user_id)
            logger.info(f"Saved to Supabase: run {run_id} (user={user_id})")
        except Exception as e:
            logger.warning(f"Supabase save failed: {e}")

        total_time = time.time() - start
        run["status"]       = "completed"
        run["progress"]     = 1.0
        run["current_node"] = None
        run["message"]      = "Complete"
        run["completed_at"] = datetime.now().isoformat()
        run["timing"] = {
            "total_seconds":  round(total_time, 1),
            "n1_curation":    round(n1_time,    1),
            "n2_screening":   round(n2_time,    1),
            "n3_tcell":       round(n3_time,    1),
            "n4_bcell":       round(n4_time,    1),
            "n5_structure":   round(n5_time,    1),
            "n6_safety":      round(n6_time,    1),
            "n7_coverage":    round(n7_time,    1),
            "n8_construct":   round(n8_time,    1),
        }
        run["candidates"] = candidates

    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        run["status"]  = "failed"
        run["message"] = str(e)
        run["progress"] = 0


def candidate_to_response(c: CandidateProtein) -> CandidateResponse:
    epitopes = []
    for ep in c.ctl_epitopes:
        to = ep.tool_outputs or {}
        epitopes.append(EpitopeResponse(
            sequence=ep.sequence, epitope_type="CTL", hla_allele=ep.hla_allele,
            ic50_nm=ep.ic50_nm, ic50_note=to.get("ic50_note","approximated_from_percentile_rank"),
            percentile_rank=ep.percentile_rank,
            confidence=ep.confidence_tier.value,
            allergenicity_safe=ep.allergenicity_safe, toxicity_safe=ep.toxicity_safe,
            safety_verdict=to.get("safety_verdict"),
            method_used=to.get("method_used"),
        ))
    for ep in c.htl_epitopes:
        to = ep.tool_outputs or {}
        epitopes.append(EpitopeResponse(
            sequence=ep.sequence, epitope_type="HTL", hla_allele=ep.hla_allele,
            ic50_nm=ep.ic50_nm, ic50_note=to.get("ic50_note","approximated_from_percentile_rank"),
            percentile_rank=ep.percentile_rank,
            confidence=ep.confidence_tier.value,
            allergenicity_safe=ep.allergenicity_safe, toxicity_safe=ep.toxicity_safe,
            safety_verdict=to.get("safety_verdict"),
            method_used=to.get("method_used"),
        ))
    for ep in c.bcell_epitopes:
        to = ep.tool_outputs or {}
        epitopes.append(EpitopeResponse(
            sequence=ep.sequence, epitope_type="B-cell", hla_allele=None,
            ic50_nm=None, ic50_note=None, percentile_rank=None,
            confidence=ep.confidence_tier.value,
            allergenicity_safe=ep.allergenicity_safe, toxicity_safe=ep.toxicity_safe,
            safety_verdict=to.get("safety_verdict"),
            method_used=to.get("method_used"),
        ))
    coverage_detail = None
    for d in c.decisions:
        if d.get("stage") == "coverage_analysis" and "per_population" in d:
            coverage_detail = d["per_population"]

    # Extract vaxijen method from decisions
    vaxijen_method = None
    for d in c.decisions:
        if d.get("stage") == "antigen_screening":
            vaxijen_method = d.get("vaxijen_method")
            break

    return CandidateResponse(
        protein_id=c.protein_id, protein_name=c.protein_name,
        sequence_length=len(c.sequence),
        ctl_count=len(c.ctl_epitopes),
        ctl_strong=len([e for e in c.ctl_epitopes if e.confidence_tier == ConfidenceTier.HIGH]),
        htl_count=len(c.htl_epitopes), bcell_count=len(c.bcell_epitopes),
        global_coverage_pct=round((c.hla_coverage_global or 0)*100, 1),
        african_coverage_pct=round((c.hla_coverage_africa or 0)*100, 1),
        structure_source=c.structure_source, structure_pdb_url=c.structure_pdb_path,
        phobius_localization=c.psortb_localization,
        vaxijen_score=c.vaxijen_score, vaxijen_method=vaxijen_method,
        epitopes=epitopes, decisions=c.decisions, coverage_detail=coverage_detail,
    )


# -- API ENDPOINTS

@app.get("/api/health")
async def health():
    """
    Extended health check — public endpoint.
    Returns N6 circuit breaker state, N7 coverage method,
    and pipeline dependency status.
    """
    n6_status = safety_filter.get_tool_status()
    n6_degraded = all(v == "open" for v in n6_status.values())

    try:
        from src.agents.predictors.coverage_agent import _load_iedb_tool
        n7_method = "iedb_tool_v3.0.1" if _load_iedb_tool() else "afnd_2020_fallback"
    except Exception:
        n7_method = "afnd_2020_fallback"

    return {
        "status": "degraded" if n6_degraded else "ok",
        "service": "kozi-pipeline",
        "version": "2.1.0",
        "pipeline": {
            "N1_data_curation":   {"status": "operational", "method": "UniProt_REST"},
            "N2_antigen_screen":  {"status": "operational", "method": "VaxiJen_2.0+Phobius",
                                   "note": "PSORTb replaced by Phobius (Docker-in-Docker not available on Railway)"},
            "N3_tcell":           {"status": "operational", "method": "IEDB_NetMHCpan4.1+MHCflurry_fallback"},
            "N4_bcell":           {"status": "operational", "method": "IEDB_BepiPred2.0"},
            "N5_structure":       {"status": "operational", "method": "AlphaFold_DB_REST"},
            "N6_safety": {
                "status": "degraded" if n6_degraded else "operational",
                "tools": n6_status,
                "note": (
                    "All safety tools unavailable. Epitopes marked unscored (None). "
                    "Do not use results for wet-lab work without manual screening."
                    if n6_degraded else None
                ),
            },
            "N7_coverage": {
                "status": "operational",
                "method": n7_method,
                "note": (
                    "Using AFND 2020 static frequency tables. "
                    "Commit src/tools/population_coverage/ to enable IEDB tool v3.0.1."
                    if n7_method == "afnd_2020_fallback" else None
                ),
            },
            "N8_construct": {"status": "operational", "method": "ProtParam_Biopython+RS09"},
        },
    }

@app.post("/api/pipeline/run", response_model=PipelineStatus)
async def start_pipeline(
    req: PipelineRequest,
    background_tasks: BackgroundTasks,
    user: UserClaims = Depends(require_user),
):
    from fastapi import BackgroundTasks
    run_id  = str(uuid.uuid4())
    user_id = user.sub
    organism_class = _infer_organism_class(req.input_value)
    candidates = []

    if req.input_type == "uniprot_id":
        prot = fetch_protein_by_id(req.input_value)
        if not prot:
            raise HTTPException(404, f"Protein {req.input_value} not found on UniProt")
        candidates.append(CandidateProtein(
            protein_id=prot["protein_id"], protein_name=prot["protein_name"],
            sequence=prot["sequence"], source="uniprot",
            stage="antigen_screening", status=CandidateStatus.ACTIVE,
        ))
    elif req.input_type == "sequence":
        seq = req.input_value.upper().replace(" ","").replace("\n","")
        if len(seq) < 10:
            raise HTTPException(400, "Sequence too short (minimum 10 amino acids)")
        candidates.append(CandidateProtein(
            protein_id="user_input", protein_name=req.protein_name or "Custom protein",
            sequence=seq, source="user_input",
            stage="antigen_screening", status=CandidateStatus.ACTIVE,
        ))
    elif req.input_type == "pathogen":
        proteins = search_pathogen_proteins(req.input_value, max_results=req.max_proteins)
        if not proteins:
            raise HTTPException(404, f"No proteins found for '{req.input_value}'")
        for p in proteins:
            seq = fetch_sequence(p["protein_id"])
            if seq and len(seq) >= 20:
                candidates.append(CandidateProtein(
                    protein_id=p["protein_id"],
                    protein_name=f"{p['protein_name']} ({p['protein_id']})",
                    sequence=seq, source="uniprot",
                    stage="antigen_screening", status=CandidateStatus.ACTIVE,
                ))
    else:
        raise HTTPException(400, f"Invalid input_type: {req.input_type}")

    if not candidates:
        raise HTTPException(400, "No valid proteins to analyze")

    active_runs[run_id] = {
        "run_id": run_id, "user_id": user_id,
        "status": "pending", "current_node": None, "progress": 0.0,
        "message": f"Starting analysis of {len(candidates)} protein(s)...",
        "input_type": req.input_type, "input_value": req.input_value,
        "protein_count": len(candidates), "organism_class": organism_class,
        "started_at": None, "completed_at": None,
        "candidates": None, "construct_report": None, "timing": None,
    }

    executor.submit(
        run_pipeline_sync, run_id, candidates,
        req.run_safety, req.run_coverage, organism_class
    )
    return PipelineStatus(
        run_id=run_id, status="pending", progress=0.0,
        message=f"Queued: {len(candidates)} protein(s) [{organism_class}]",
    )

@app.get("/api/pipeline/status/{run_id}", response_model=PipelineStatus)
async def get_pipeline_status(run_id: str, user: UserClaims = Depends(require_user)):
    run = active_runs.get(run_id)
    if not run:
        raise HTTPException(404, "Run not found")
    if run.get("user_id") and run["user_id"] != user.sub:
        raise HTTPException(403, "Access denied")
    return PipelineStatus(
        run_id=run["run_id"], status=run["status"],
        current_node=run.get("current_node"), progress=run.get("progress",0),
        message=run.get("message"), started_at=run.get("started_at"),
        completed_at=run.get("completed_at"),
    )

@app.get("/api/pipeline/results/{run_id}", response_model=PipelineResult)
async def get_pipeline_results(run_id: str, user: UserClaims = Depends(require_user)):
    """
    Results endpoint — checks active_runs first, falls back to Supabase.
    Prevents 404 after Railway restarts wipe in-memory state.
    """
    run = active_runs.get(run_id)

    if run:
        if run.get("user_id") and run["user_id"] != user.sub:
            raise HTTPException(403, "Access denied")
        if run["status"] != "completed":
            raise HTTPException(409, f"Run is {run['status']}, not completed")
        return PipelineResult(
            run_id=run_id, status="completed",
            timing=run.get("timing",{}),
            candidates=[candidate_to_response(c) for c in run.get("candidates",[])],
            construct_report=run.get("construct_report"),
        )

    logger.info(f"Run {run_id} not in active_runs — checking Supabase")
    reconstructed = reconstruct_from_supabase(run_id, user.sub)

    if reconstructed is None:
        try:
            from src.storage.supabase_client import db
            row = db.get_run(run_id)
            if row and row.get("user_id") and row["user_id"] != user.sub:
                raise HTTPException(403, "Access denied")
        except HTTPException:
            raise
        except Exception:
            pass
        raise HTTPException(404, "Run not found")

    return PipelineResult(
        run_id=run_id, status=reconstructed["status"],
        timing=reconstructed["timing"],
        candidates=reconstructed["candidates"],
        construct_report=reconstructed.get("construct_report"),
    )

@app.get("/api/runs")
async def list_runs(user: UserClaims = Depends(require_user)):
    runs = []
    for run_id, run in active_runs.items():
        if run.get("user_id") and run["user_id"] != user.sub:
            continue
        global_coverage = None
        if run.get("candidates"):
            covs = [
                round((c.hla_coverage_global or 0)*100, 1)
                for c in run["candidates"]
                if hasattr(c,"hla_coverage_global") and c.hla_coverage_global
            ]
            if covs: global_coverage = max(covs)
        runs.append({
            "run_id": run_id, "status": run["status"],
            "input_type": run.get("input_type"), "input_value": run.get("input_value"),
            "protein_count": run.get("protein_count"),
            "progress": run.get("progress",0),
            "started_at": run.get("started_at"), "completed_at": run.get("completed_at"),
            "global_coverage_pct": global_coverage,
            "has_construct": run.get("construct_report") is not None,
        })
    try:
        from src.storage.supabase_client import db
        db_runs = (
            db.client.table("runs").select("*")
            .or_(f"user_id.eq.{user.sub},user_id.is.null")
            .order("created_at", desc=True).limit(50).execute()
        )
        existing_ids = {r["run_id"] for r in runs}
        for db_run in db_runs.data or []:
            if db_run["id"] not in existing_ids:
                cands = (
                    db.client.table("candidates")
                    .select("hla_coverage_global").eq("run_id", db_run["id"]).execute()
                )
                cov = None
                if cands.data:
                    covs = [round((c["hla_coverage_global"] or 0)*100,1)
                            for c in cands.data if c.get("hla_coverage_global")]
                    cov = max(covs) if covs else None
                runs.append({
                    "run_id": db_run["id"], "status": db_run.get("status","unknown"),
                    "input_type": db_run.get("input_type"),
                    "input_value": db_run.get("pathogen_name"),
                    "protein_count": None,
                    "progress": 1.0 if db_run.get("status") == "completed" else 0,
                    "started_at": db_run.get("created_at"),
                    "completed_at": db_run.get("completed_at"),
                    "global_coverage_pct": cov, "has_construct": False,
                })
    except Exception as e:
        logger.warning(f"Supabase list_runs failed: {e}")
    return sorted(runs, key=lambda x: x.get("started_at") or "", reverse=True)


@app.websocket("/ws/pipeline/{run_id}")
async def pipeline_websocket(
    websocket: WebSocket,
    run_id: str,
    token: Optional[str] = Query(None),
):
    await websocket.accept()
    if not token:
        await websocket.send_json({"error": "Missing token", "code": 401})
        await websocket.close(code=4001)
        return
    try:
        from api.auth import _verify_token
        user = _verify_token(token)
    except Exception:
        await websocket.send_json({"error": "Invalid or expired token", "code": 401})
        await websocket.close(code=4001)
        return
    try:
        last_progress = -1
        while True:
            run = active_runs.get(run_id)
            if not run:
                await websocket.send_json({"error": "Run not found"})
                break
            if run.get("user_id") and run["user_id"] != user.sub:
                await websocket.send_json({"error": "Access denied", "code": 403})
                break
            progress = run.get("progress", 0)
            if progress != last_progress:
                await websocket.send_json({
                    "run_id": run_id, "status": run["status"],
                    "current_node": run.get("current_node"),
                    "progress": progress, "message": run.get("message"),
                })
                last_progress = progress
            if run["status"] in ("completed","failed"):
                await websocket.send_json({
                    "run_id": run_id, "status": run["status"],
                    "progress": 1.0 if run["status"] == "completed" else progress,
                    "message": run.get("message"), "timing": run.get("timing"),
                })
                break
            await asyncio.sleep(1)
    except WebSocketDisconnect:
        logger.info(f"WebSocket disconnected for run {run_id}")
    except Exception as e:
        logger.error(f"WebSocket error: {e}")


@app.on_event("startup")
async def startup():
    logger.info("Kozi AI TOPE_DEEP API starting (v2.1.0)...")
    if not os.getenv("SUPABASE_JWT_SECRET"):
        logger.warning("SUPABASE_JWT_SECRET not set — protected endpoints will return 500")
    # Log N7 method at startup
    try:
        from src.agents.predictors.coverage_agent import _load_iedb_tool
        n7 = "IEDB tool v3.0.1" if _load_iedb_tool() else "AFND 2020 fallback"
        logger.info(f"N7 coverage method: {n7}")
    except Exception:
        logger.warning("N7 coverage method check failed")