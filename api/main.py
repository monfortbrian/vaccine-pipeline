"""
TOPE_DEEP REST API v1.0.0

Pure HTTP wrapper. All pipeline logic lives in src/orchestrator/.

Endpoints:
  POST /api/v1/pipeline/run           Start a run
  GET  /api/v1/pipeline/status/{id}   Poll run progress
  GET  /api/v1/pipeline/results/{id}  Fetch completed results
  GET  /api/v1/runs                   List runs (paginated)
  GET  /api/v1/health                 System and agent health (public)
  WS   /ws/pipeline/{id}              Real-time progress stream

Legacy routes /api/pipeline/* and /api/runs kept for frontend compatibility.
"""

from dotenv import load_dotenv
load_dotenv()

from src.utils.logger import configure_root
configure_root()

import os
import sys
import uuid
import time
import json
import logging
import asyncio
from datetime import datetime
from typing import List, Dict, Any, Optional
from concurrent.futures import ThreadPoolExecutor

import requests as http_requests
from fastapi import (
    FastAPI, HTTPException, Depends, Query,
    Request, WebSocket, WebSocketDisconnect, APIRouter
)
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from api.auth import require_user, UserClaims

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("tope_deep.api")


# ── App ───────────────────────────────────────────────────────────────────────

app = FastAPI(
    title="TOPE_DEEP REST API",
    description=(
        "A computational vaccine discovery pipeline. 10 Agents, one workflow, full audit trail."
    ),
    version="1.0.0",
)

limiter = Limiter(key_func=get_remote_address)
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        o.strip()
        for o in os.getenv("CORS_ORIGINS", "http://localhost:3000").split(",")
        if o.strip()
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

executor = ThreadPoolExecutor(max_workers=4)


# ── Job state: Redis with in-memory fallback ──────────────────────────────────

_memory_store: Dict[str, Dict] = {}


def _redis_client():
    url = os.getenv("REDIS_URL", "")
    if not url:
        return None
    try:
        import redis
        return redis.from_url(url, decode_responses=True, socket_timeout=2)
    except Exception:
        return None


def _set_run(run_id: str, data: Dict):
    r = _redis_client()
    if r:
        try:
            r.setex(f"tope:run:{run_id}", 86400, json.dumps(data, default=str))
            return
        except Exception as e:
            logger.warning(f"Redis write failed, using memory: {e}")
    _memory_store[run_id] = data


def _get_run(run_id: str) -> Optional[Dict]:
    r = _redis_client()
    if r:
        try:
            v = r.get(f"tope:run:{run_id}")
            if v:
                return json.loads(v)
        except Exception as e:
            logger.warning(f"Redis read failed, checking memory: {e}")
    return _memory_store.get(run_id)


def _patch_run(run_id: str, patch: Dict):
    data = _get_run(run_id) or {}
    data.update(patch)
    _set_run(run_id, data)


# ── Pydantic models ───────────────────────────────────────────────────────────

class RunRequest(BaseModel):
    input_type:      str           = Field(description="pathogen | uniprot_id | sequence")
    input_value:     str           = Field(description="Pathogen name, UniProt accession, or amino acid sequence")
    protein_name:    Optional[str] = None
    run_safety:      bool          = True
    run_coverage:    bool          = True
    run_literature:  bool          = True
    run_experiment:  bool          = True
    max_proteins:    int           = Field(3, ge=1, le=10)
    lab_constraints: str           = "standard academic lab"


class RunStatus(BaseModel):
    run_id:        str
    status:        str
    current_agent: Optional[str] = None
    progress:      float         = 0.0
    message:       Optional[str] = None
    started_at:    Optional[str] = None
    completed_at:  Optional[str] = None


class EpitopeOut(BaseModel):
    model_config = {"protected_namespaces": ()}
    sequence:           str
    epitope_type:       str
    hla_allele:         Optional[str]
    ic50_nm:            Optional[float]
    percentile_rank:    Optional[float]
    confidence:         str
    allergenicity_safe: Optional[bool]
    toxicity_safe:      Optional[bool]
    model_categories:   List[str] = []
    safety_verdict:     Optional[str] = None
    method_used:        Optional[str] = None


class CandidateOut(BaseModel):
    protein_id:           str
    protein_name:         str
    sequence_length:      int
    ctl_count:            int
    ctl_strong:           int
    htl_count:            int
    bcell_count:          int
    global_coverage_pct:  float
    african_coverage_pct: float
    structure_source:     Optional[str]   = None
    structure_pdb_url:    Optional[str]   = None
    phobius_localization: Optional[str]   = None
    vaxijen_score:        Optional[float] = None
    vaxijen_method:       Optional[str]   = None
    epitopes:             List[EpitopeOut]
    decisions:            List[Dict[str, Any]]
    coverage_detail:      Optional[Dict[str, Any]] = None


class RunResult(BaseModel):
    run_id:           str
    status:           str
    timing:           Dict[str, Any]
    candidates:       List[CandidateOut]
    construct_report: Optional[Dict[str, Any]] = None


class RunListItem(BaseModel):
    run_id:              str
    status:              str
    input_type:          Optional[str]
    input_value:         Optional[str]
    protein_count:       Optional[int]
    progress:            float
    started_at:          Optional[str]
    completed_at:        Optional[str]
    global_coverage_pct: Optional[float]
    has_construct:       bool


class RunListResponse(BaseModel):
    runs:     List[RunListItem]
    total:    int
    page:     int
    per_page: int
    has_more: bool


# ── UniProt helpers ───────────────────────────────────────────────────────────

def _infer_organism(v: str) -> str:
    v = v.lower()
    if any(w in v for w in ["virus","viral","sars","covid","influenza","hiv","ebola","dengue","hpv","zika"]):
        return "virus"
    if any(w in v for w in ["plasmodium","malaria","leishmania","trypanosoma","parasite","helminth","schistosoma"]):
        return "parasite"
    return "bacteria"


def _fetch_by_id(uid: str) -> Optional[Dict]:
    try:
        r = http_requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=15)
        r.raise_for_status()
        lines = r.text.strip().split("\n")
        seq   = "".join(l for l in lines[1:])
        r2    = http_requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uid}.json", timeout=15)
        name  = lines[0][1:].split()[0]
        if r2.status_code == 200:
            name = (r2.json()
                    .get("proteinDescription", {})
                    .get("recommendedName", {})
                    .get("fullName", {})
                    .get("value", name))
        return {"protein_id": uid, "protein_name": name, "sequence": seq}
    except Exception as e:
        logger.error(f"UniProt fetch failed for {uid}: {e}")
        return None


def _search_pathogen(name: str, n: int = 5) -> List[Dict]:
    try:
        r = http_requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={
                "query":  f'(organism_name:"{name}") AND (reviewed:true)',
                "format": "json",
                "size":   n * 3,
                "fields": "accession,protein_name,length,cc_subcellular_location",
            },
            timeout=30,
        )
        r.raise_for_status()
        out = []
        for e in r.json().get("results", []):
            acc   = e.get("primaryAccession", "")
            pname = (e.get("proteinDescription", {})
                     .get("recommendedName", {})
                     .get("fullName", {})
                     .get("value", "Unknown"))
            loc   = " ".join(
                loc_item.get("location", {}).get("value", "")
                for c in e.get("comments", [])
                if c.get("commentType") == "SUBCELLULAR LOCATION"
                for loc_item in c.get("subcellularLocations", [])
            )
            surface = any(
                k in loc.lower()
                for k in ["membrane", "secreted", "cell surface", "extracellular", "outer membrane"]
            )
            out.append({
                "protein_id":   acc,
                "protein_name": pname,
                "is_surface":   surface,
                "length":       e.get("sequence", {}).get("length", 0),
            })
        out.sort(key=lambda x: not x["is_surface"])
        return out[:n]
    except Exception as e:
        logger.error(f"UniProt pathogen search failed: {e}")
        return []


def _fetch_seq(uid: str) -> Optional[str]:
    try:
        r = http_requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=15)
        r.raise_for_status()
        lines = r.text.strip().split("\n")
        return "".join(l for l in lines[1:])
    except Exception:
        return None


# ── Supabase reconstruction (results fallback after restart) ──────────────────

def _reconstruct_from_supabase(run_id: str, user_id: str) -> Optional[Dict]:
    try:
        from src.storage.supabase_client import db

        row = db.get_run(run_id)
        if not row:
            logger.warning(f"Reconstruction: run {run_id} not found in Supabase")
            return None

        # FIX: allow reconstruction for runs with NULL user_id (saved before auth was set)
        # or runs belonging to this user
        row_user = row.get("user_id")
        if row_user and row_user != user_id:
            logger.warning(f"Reconstruction: run {run_id} belongs to different user")
            return None

        candidate_rows = db.get_candidates_for_run(run_id)
        if not candidate_rows:
            logger.warning(f"Reconstruction: no candidates for run {run_id}")
            return None

        candidates_out = []
        for cr in candidate_rows:
            ep_rows = db.get_epitopes_for_candidate(cr["id"])
            eps = [
                EpitopeOut(
                    sequence=e["sequence"],
                    epitope_type=e["epitope_type"],
                    hla_allele=e.get("hla_allele"),
                    ic50_nm=e.get("ic50_nm"),
                    percentile_rank=e.get("percentile_rank"),
                    confidence=e.get("confidence_tier", "uncertain"),
                    allergenicity_safe=e.get("allergenicity_safe"),
                    toxicity_safe=e.get("toxicity_safe"),
                    model_categories=(e.get("tool_outputs") or {}).get("model_categories", []),
                    safety_verdict=(e.get("tool_outputs") or {}).get("safety_verdict"),
                    method_used=(e.get("tool_outputs") or {}).get("method_used"),
                )
                for e in ep_rows
            ]
            ctl   = [e for e in eps if e.epitope_type == "CTL"]
            htl   = [e for e in eps if e.epitope_type == "HTL"]
            bcell = [e for e in eps if "B-cell" in e.epitope_type]

            decisions = cr.get("decisions") or []
            cov_detail = next(
                (d.get("per_population") for d in decisions
                 if d.get("stage") == "coverage_analysis" and d.get("per_population")),
                None
            )
            vaxijen_method = next(
                (d.get("vaxijen_method") for d in decisions
                 if d.get("stage") == "antigen_screening"),
                None
            )
            phobius_loc = next(
                (d.get("phobius_localization") for d in decisions
                 if d.get("stage") == "antigen_screening"),
                None
            )

            candidates_out.append(CandidateOut(
                protein_id=cr["protein_id"],
                protein_name=cr["protein_name"],
                sequence_length=len(cr.get("sequence") or ""),
                ctl_count=len(ctl),
                ctl_strong=sum(1 for e in ctl if e.confidence == "high"),
                htl_count=len(htl),
                bcell_count=len(bcell),
                global_coverage_pct=round((cr.get("hla_coverage_global") or 0) * 100, 1),
                african_coverage_pct=round((cr.get("hla_coverage_africa") or 0) * 100, 1),
                structure_source=cr.get("structure_source"),
                structure_pdb_url=cr.get("structure_pdb_path"),
                phobius_localization=phobius_loc,
                vaxijen_score=cr.get("vaxijen_score"),
                vaxijen_method=vaxijen_method,
                epitopes=eps,
                decisions=decisions,
                coverage_detail=cov_detail,
            ))

        # Reconstruct timing from candidates decisions if available
        timing = {"total_seconds": 0}
        try:
            timing_row = row.get("timing") or {}
            if timing_row:
                timing = timing_row
        except Exception:
            pass

        return {
            "run_id":           run_id,
            "status":           row.get("status", "completed"),
            "timing":           timing,
            "candidates":       candidates_out,
            "construct_report": row.get("construct_report"),
        }
    except Exception as e:
        logger.error(f"Supabase reconstruction failed for {run_id}: {e}", exc_info=True)
        return None


# ── Candidate serializer ──────────────────────────────────────────────────────

def _serialize_candidate(c) -> CandidateOut:
    from src.models.candidate import ConfidenceTier

    def _ep(ep_, type_str: str) -> EpitopeOut:
        to = ep_.tool_outputs or {}
        return EpitopeOut(
            sequence=ep_.sequence,
            epitope_type=type_str,
            hla_allele=ep_.hla_allele,
            ic50_nm=ep_.ic50_nm,
            percentile_rank=ep_.percentile_rank,
            confidence=ep_.confidence_tier.value,
            allergenicity_safe=ep_.allergenicity_safe,
            toxicity_safe=ep_.toxicity_safe,
            model_categories=to.get("model_categories", ["HUMAN"]),
            safety_verdict=to.get("safety_verdict"),
            method_used=to.get("method_used"),
        )

    eps = (
        [_ep(e, "CTL")    for e in c.ctl_epitopes] +
        [_ep(e, "HTL")    for e in c.htl_epitopes] +
        [_ep(e, "B-cell") for e in c.bcell_epitopes]
    )

    cov_detail = next(
        (d.get("per_population") for d in c.decisions
         if d.get("stage") == "coverage_analysis" and d.get("per_population")),
        None
    )
    vaxijen_method = next(
        (d.get("vaxijen_method") for d in c.decisions
         if d.get("stage") == "antigen_screening"),
        None
    )

    return CandidateOut(
        protein_id=c.protein_id,
        protein_name=c.protein_name,
        sequence_length=len(c.sequence),
        ctl_count=len(c.ctl_epitopes),
        ctl_strong=sum(1 for e in c.ctl_epitopes if e.confidence_tier == ConfidenceTier.HIGH),
        htl_count=len(c.htl_epitopes),
        bcell_count=len(c.bcell_epitopes),
        global_coverage_pct=round((c.hla_coverage_global or 0) * 100, 1),
        african_coverage_pct=round((c.hla_coverage_africa or 0) * 100, 1),
        structure_source=c.structure_source,
        structure_pdb_url=c.structure_pdb_path,
        phobius_localization=c.psortb_localization,
        vaxijen_score=c.vaxijen_score,
        vaxijen_method=vaxijen_method,
        epitopes=eps,
        decisions=c.decisions,
        coverage_detail=cov_detail,
    )


# ── Pipeline runner (runs in thread pool) ─────────────────────────────────────

def _execute_pipeline(run_id: str, candidates: list, config: Dict):
   
    from src.orchestrator import PipelineOrchestrator

    def _progress(agent: str, pct: float, msg: str):
        _patch_run(run_id, {
            "current_agent": agent,
            "progress":      pct,
            "message":       msg,
        })

    # Get user_id BEFORE pipeline runs (Redis has it at this point)
    user_id = (_get_run(run_id) or {}).get("user_id")

    result       = {}
    final_status = "failed"

    try:
        _patch_run(run_id, {"status": "running", "started_at": datetime.now().isoformat()})

        orch   = PipelineOrchestrator()
        result = orch.run(run_id, candidates, config, progress_callback=_progress)

        final_status = "completed"
        _patch_run(run_id, {
            "status":                "completed",
            "progress":              1.0,
            "current_agent":         None,
            "message":               "Complete",
            "completed_at":          datetime.now().isoformat(),
            "timing":                result.get("timing", {}),
            "construct_report":      result.get("construct_report"),
            "candidates_serialized": [
                _serialize_candidate(c).model_dump()
                for c in result.get("candidates", [])
            ],
        })

    except Exception as e:
        logger.error(f"Pipeline failed [{run_id}]: {e}", exc_info=True)
        _patch_run(run_id, {
            "status":   "failed",
            "message":  str(e),
            "progress": 0,
        })

    finally:
        # ── ALWAYS write to Supabase ──────────────────────────────────────────
        # Even partial / failed runs are persisted so history page shows them.
        # This block runs whether pipeline succeeded, partially failed, or crashed.
        try:
            from src.storage.supabase_client import db
            from src.models.candidate import PipelineRun as PR

            completed_at = datetime.now().isoformat() if final_status == "completed" else None
            timing       = result.get("timing", {}) if result else {}

            # Upsert the run row (create or update)
            db.upsert_run(
                run_id       = run_id,
                user_id      = user_id,
                input_type   = config.get("input_type", "unknown"),
                pathogen_name= (candidates[0].protein_name if candidates else "unknown"),
                status       = final_status,
                completed_at = completed_at,
                timing       = timing,
            )

            # Save candidates (only if pipeline produced them)
            for c in result.get("candidates", []):
                try:
                    db.save_candidate(run_id, c, user_id=user_id)
                except Exception as ce:
                    logger.warning(f"Candidate save failed [{run_id}]: {ce}")

            logger.info(f"Supabase persisted [{run_id}] status={final_status}")

        except Exception as e:
            logger.error(f"Supabase finally-block save failed [{run_id}]: {e}", exc_info=True)


# ── v1 router ─────────────────────────────────────────────────────────────────

v1 = APIRouter(prefix="/api/v1")


@v1.get("/health", tags=["system"])
async def health():
    from src.agents.safety_filter import safety_filter

    n6_status   = safety_filter.get_tool_status()
    n6_degraded = all(
        v in ("open", "degraded_no_db")
        for k, v in n6_status.items()
        if k not in ("allergenonline_version", "human_swissprot_version")
    )

    n7_method = "afnd_2020_fallback"
    try:
        from src.agents.coverage_agent import _load_iedb_tool
        if _load_iedb_tool():
            n7_method = "iedb_tool_v3.0.1"
    except Exception:
        pass

    qdrant_ok = False
    try:
        from qdrant_client import QdrantClient
        QdrantClient(":memory:")
        qdrant_ok = True
    except ImportError:
        pass

    redis_ok = False
    try:
        r = _redis_client()
        if r:
            r.ping()
            redis_ok = True
    except Exception:
        pass

    return {
        "status":  "degraded" if n6_degraded else "ok",
        "version": "1.0.0",
        "infrastructure": {
            "redis":    "connected" if redis_ok else "in_memory_fallback",
            "supabase": "configured" if os.getenv("SUPABASE_URL") else "not_configured",
        },
        "agents": {
            "Data curator": {
                "status": "operational",
                "tool":   "UniProt REST + NCBI",
            },
            "Antigen screener": {
                "status": "operational",
                "tool":   "VaxiJen 2.0 ACC local + Phobius 2.0",
                "note":   "PSORTb unavailable on Railway, Phobius used as proxy",
            },
            "T-cell predictor": {
                "status": "operational",
                "tool":   "IEDB NetMHCpan 4.1 EL + NetMHCIIpan 4.3 + MHCflurry 2.0 fallback",
            },
            "B-cell predictor": {
                "status": "operational",
                "tool":   "IEDB BepiPred 2.0",
            },
            "Structure agent": {
                "status": "operational",
                "tool":   "AlphaFold DB REST (EBI)",
            },
            "Safety filter": {
                "status": "degraded" if n6_degraded else "operational",
                "tool":   "FAO/WHO 2001 + AllerTOP v2.0 + HemoPI + FDA/EMA 8-mer",
                "detail": n6_status,
            },
            "Coverage agent": {
                "status":      "operational",
                "tool":        n7_method,
            },
            "Construct designer": {
                "status": "operational",
                "tool":   "ProtParam (Biopython)",
            },
            "Literature agent": {
                "status":          "operational" if qdrant_ok else "degraded",
                "tool":            "PubMed E-utilities + Qdrant + sentence-transformers",
                "claude_synthesis": bool(os.getenv("ANTHROPIC_API_KEY")),
            },
            "Experiment planner": {
                "status":           "operational",
                "tool":             "Claude API + template fallback",
                "claude_available": bool(os.getenv("ANTHROPIC_API_KEY")),
            },
        },
    }


@v1.post("/pipeline/run", response_model=RunStatus, tags=["pipeline"])
@limiter.limit("5/minute")
async def start_run(
    request: Request,
    body: RunRequest,
    user: UserClaims = Depends(require_user),
):
    from src.models.candidate import CandidateProtein, CandidateStatus
    from src.validation.candidate_validator import validate_input

    val = validate_input(body.input_type, body.input_value)
    if not val.valid:
        raise HTTPException(
            status_code=400,
            detail={"error": val.error, "suggestions": val.suggestions},
        )

    run_id   = str(uuid.uuid4())
    organism = _infer_organism(body.input_value)
    candidates = []

    if body.input_type == "uniprot_id":
        p = _fetch_by_id(body.input_value)
        if not p:
            raise HTTPException(404, f"UniProt accession '{body.input_value}' not found or unreachable.")
        candidates.append(CandidateProtein(
            protein_id=p["protein_id"], protein_name=p["protein_name"],
            sequence=p["sequence"], source="uniprot",
            stage="antigen_screening", status=CandidateStatus.ACTIVE,
        ))

    elif body.input_type == "sequence":
        seq = body.input_value.upper().replace(" ", "").replace("\n", "")
        if len(seq) < 20:
            raise HTTPException(400, "Sequence too short, minimum 20 amino acids.")
        candidates.append(CandidateProtein(
            protein_id="user_input",
            protein_name=body.protein_name or "Custom protein",
            sequence=seq, source="user_input",
            stage="antigen_screening", status=CandidateStatus.ACTIVE,
        ))

    elif body.input_type == "pathogen":
        proteins = _search_pathogen(body.input_value, n=body.max_proteins)
        if not proteins:
            raise HTTPException(
                404,
                f"No reviewed UniProt entries found for '{body.input_value}'. "
                f"Check spelling or try the full scientific name."
            )
        for p in proteins:
            seq = _fetch_seq(p["protein_id"])
            if seq and len(seq) >= 20:
                candidates.append(CandidateProtein(
                    protein_id=p["protein_id"],
                    protein_name=f"{p['protein_name']} ({p['protein_id']})",
                    sequence=seq, source="uniprot",
                    stage="antigen_screening", status=CandidateStatus.ACTIVE,
                ))
    else:
        raise HTTPException(400, "input_type must be 'pathogen', 'uniprot_id', or 'sequence'.")

    if not candidates:
        raise HTTPException(400, "No valid protein sequences to analyse.")

    config = {
        "organism_class":  organism,
        "input_type":      body.input_type,
        "run_safety":      body.run_safety,
        "run_coverage":    body.run_coverage,
        "run_literature":  body.run_literature,
        "run_experiment":  body.run_experiment,
        "lab_constraints": body.lab_constraints,
    }

    _set_run(run_id, {
        "run_id":        run_id,
        "user_id":       user.sub,
        "status":        "pending",
        "current_agent": None,
        "progress":      0.0,
        "message":       f"Queued {len(candidates)} protein(s) [{organism}]",
        "input_type":    body.input_type,
        "input_value":   body.input_value,
        "protein_count": len(candidates),
        "started_at":    None,
        "completed_at":  None,
        "input_warning": val.warning,
    })

    executor.submit(_execute_pipeline, run_id, candidates, config)

    return RunStatus(
        run_id=run_id,
        status="pending",
        progress=0.0,
        message=f"Queue {len(candidates)} protein(s)" + (
            f" | Warning: {val.warning}" if val.warning else ""
        ),
    )


@v1.get("/pipeline/status/{run_id}", response_model=RunStatus, tags=["pipeline"])
async def get_status(run_id: str, user: UserClaims = Depends(require_user)):
    run = _get_run(run_id)
    if not run:
        raise HTTPException(404, "Run not found.")
    if run.get("user_id") and run["user_id"] != user.sub:
        raise HTTPException(403, "Access denied.")
    return RunStatus(
        run_id=run_id,
        status=run["status"],
        current_agent=run.get("current_agent"),
        progress=run.get("progress", 0),
        message=run.get("message"),
        started_at=run.get("started_at"),
        completed_at=run.get("completed_at"),
    )


@v1.get("/pipeline/results/{run_id}", response_model=RunResult, tags=["pipeline"])
async def get_results(run_id: str, user: UserClaims = Depends(require_user)):
    run = _get_run(run_id)

    if run:
        if run.get("user_id") and run["user_id"] != user.sub:
            raise HTTPException(403, "Access denied.")
        if run["status"] not in ("completed", "failed"):
            raise HTTPException(409, f"Run status is '{run['status']}', not completed yet.")
        if run["status"] == "failed":
            # Try Supabase for partial results even on failed runs
            rec = _reconstruct_from_supabase(run_id, user.sub)
            if rec and rec.get("candidates"):
                return RunResult(
                    run_id=run_id,
                    status="completed",  # partial but usable
                    timing=rec["timing"],
                    candidates=rec["candidates"],
                    construct_report=rec.get("construct_report"),
                )
            raise HTTPException(409, "Run failed with no partial results.")

        candidates = [CandidateOut(**c) for c in run.get("candidates_serialized", [])]
        return RunResult(
            run_id=run_id,
            status="completed",
            timing=run.get("timing", {}),
            candidates=candidates,
            construct_report=run.get("construct_report"),
        )

    # Not in job store, try Supabase (covers Railway restarts, expired Redis)
    logger.info(f"Run {run_id} not in job store, checking Supabase")
    rec = _reconstruct_from_supabase(run_id, user.sub)
    if not rec:
        raise HTTPException(
            404,
            "Run not found. It may have expired or belong to another account."
        )
    return RunResult(
        run_id=run_id,
        status=rec["status"],
        timing=rec["timing"],
        candidates=rec["candidates"],
        construct_report=rec.get("construct_report"),
    )


@v1.get("/runs", response_model=RunListResponse, tags=["runs"])
async def list_runs(
    user:     UserClaims = Depends(require_user),
    page:     int = Query(1, ge=1, description="Page number"),
    per_page: int = Query(20, ge=1, le=100, description="Results per page"),
):
    runs  = []
    total = 0

    try:
        from src.storage.supabase_client import db
        offset = (page - 1) * per_page

        count_res = (
            db.client.table("runs")
            .select("id", count="exact")
            .or_(f"user_id.eq.{user.sub},user_id.is.null")
            .execute()
        )
        total = count_res.count or 0

        db_runs = (
            db.client.table("runs")
            .select("*")
            .or_(f"user_id.eq.{user.sub},user_id.is.null")
            .order("created_at", desc=True)
            .range(offset, offset + per_page - 1)
            .execute()
        )

        for r in db_runs.data or []:
            cands = (
                db.client.table("candidates")
                .select("hla_coverage_global")
                .eq("run_id", r["id"])
                .execute()
            )
            covs = [
                round((c["hla_coverage_global"] or 0) * 100, 1)
                for c in (cands.data or [])
                if c.get("hla_coverage_global")
            ]
            runs.append(RunListItem(
                run_id=r["id"],
                status=r.get("status", "unknown"),
                input_type=r.get("input_type"),
                input_value=r.get("pathogen_name"),
                protein_count=None,
                progress=1.0 if r.get("status") == "completed" else 0.0,
                started_at=r.get("created_at"),
                completed_at=r.get("completed_at"),
                global_coverage_pct=max(covs) if covs else None,
                has_construct=False,
            ))
    except Exception as e:
        logger.error(f"Supabase list_runs failed: {e}", exc_info=True)

    return RunListResponse(
        runs=runs,
        total=total,
        page=page,
        per_page=per_page,
        has_more=(page * per_page) < total,
    )


# ── WebSocket ─────────────────────────────────────────────────────────────────

@app.websocket("/ws/pipeline/{run_id}")
async def ws_pipeline(
    websocket: WebSocket,
    run_id:    str,
    token:     Optional[str] = Query(None),
):
    await websocket.accept()

    if not token:
        await websocket.send_json({"error": "token required", "code": 401})
        await websocket.close(code=4001)
        return

    try:
        from api.auth import _verify_token
        ws_user = _verify_token(token)
    except Exception:
        await websocket.send_json({"error": "invalid or expired token", "code": 401})
        await websocket.close(code=4001)
        return

    try:
        last_progress = -1.0
        while True:
            run = _get_run(run_id)
            if not run:
                await websocket.send_json({"error": "run not found"})
                break
            if run.get("user_id") and run["user_id"] != ws_user.sub:
                await websocket.send_json({"error": "access denied"})
                break

            progress = run.get("progress", 0)
            if progress != last_progress:
                await websocket.send_json({
                    "run_id":        run_id,
                    "status":        run["status"],
                    "current_agent": run.get("current_agent"),
                    "progress":      progress,
                    "message":       run.get("message"),
                })
                last_progress = progress

            if run["status"] in ("completed", "failed"):
                await websocket.send_json({
                    "run_id":   run_id,
                    "status":   run["status"],
                    "progress": 1.0 if run["status"] == "completed" else progress,
                    "timing":   run.get("timing"),
                })
                break

            await asyncio.sleep(1)

    except WebSocketDisconnect:
        logger.info(f"WebSocket disconnected: {run_id}")
    except Exception as e:
        logger.error(f"WebSocket error [{run_id}]: {e}")


# ── Register routers ──────────────────────────────────────────────────────────

app.include_router(v1)

legacy = APIRouter(prefix="/api")

@legacy.get("/health")
async def _legacy_health():
    return await health()

@legacy.post("/pipeline/run")
@limiter.limit("5/minute")
async def _legacy_run(request: Request, body: RunRequest, user: UserClaims = Depends(require_user)):
    return await start_run(request, body, user)

@legacy.get("/pipeline/status/{run_id}")
async def _legacy_status(run_id: str, user: UserClaims = Depends(require_user)):
    return await get_status(run_id, user)

@legacy.get("/pipeline/results/{run_id}")
async def _legacy_results(run_id: str, user: UserClaims = Depends(require_user)):
    return await get_results(run_id, user)

@legacy.get("/runs")
async def _legacy_runs(request: Request, user: UserClaims = Depends(require_user)):
    page     = int(request.query_params.get("page", 1))
    per_page = int(request.query_params.get("per_page", 20))
    return await list_runs(user, page, per_page)

app.include_router(legacy)


# ── Root ──────────────────────────────────────────────────────────────────────

@app.get("/")
def root():
    return {
        "service": "TOPE_DEEP API",
        "version": "1.0.0",
        "docs":    "/docs",
        "health":  "/api/v1/health",
    }


@app.get("/api/test-db")
def test_db():
    from src.storage.supabase_client import db
    return {"connected": db.test_connection()}


# ── Startup ───────────────────────────────────────────────────────────────────

@app.on_event("startup")
async def startup():
    logger.info("TOPE_DEEP API v1.0.0 starting...")

    r = _redis_client()
    if r:
        try:
            r.ping()
            logger.info("Redis: connected")
        except Exception:
            logger.warning("Redis: unreachable, in-memory job store active")
    else:
        logger.warning("Redis: REDIS_URL not set in-memory only. Railway deploys will lose run state.")

    if not os.getenv("SUPABASE_JWT_SECRET"):
        logger.warning("SUPABASE_JWT_SECRET not set, auth will fail")

    if not os.getenv("ANTHROPIC_API_KEY"):
        logger.info("ANTHROPIC_API_KEY not set, It will use template fallback")

    try:
        from src.agents.coverage_agent import _load_iedb_tool
        method = "IEDB tool v3.0.1" if _load_iedb_tool() else "AFND 2020 fallback"
        logger.info(f"Coverage: {method}")
    except Exception:
        logger.warning("Coverage agent check failed")

    try:
        from qdrant_client import QdrantClient
        QdrantClient(":memory:")
        logger.info("Qdrant: in-memory ready")
    except ImportError:
        logger.warning("Qdrant: not installed, pip install qdrant-client")

    try:
        from sentence_transformers import SentenceTransformer
        logger.info("SentenceTransformers: ready")
    except ImportError:
        logger.warning("SentenceTransformers: not installed")