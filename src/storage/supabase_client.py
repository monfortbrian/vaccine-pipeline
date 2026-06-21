"""
Persistence layer.
Service role key used for all server-side writes (bypasses RLS).
"""

import os
import json
import logging
from typing import List, Optional, Dict, Any
from datetime import datetime
from supabase import create_client, Client
from src.models.candidate import PipelineRun, CandidateProtein, EpitopeResult

if os.getenv("RAILWAY_ENVIRONMENT") is None:
    from dotenv import load_dotenv
    load_dotenv()

logger = logging.getLogger(__name__)


class SupabaseClient:
    def __init__(self):
        self.url = os.getenv("SUPABASE_URL")
        self.key = os.getenv("SUPABASE_SERVICE_KEY") or os.getenv("SUPABASE_KEY")
        if not self.url or not self.key:
            raise ValueError("SUPABASE_URL and SUPABASE_SERVICE_KEY must be set")
        self.client: Client = create_client(self.url, self.key)

    # ── RUNS ──────────────────────────────────────────────────────────────────

    def upsert_run(
        self,
        run_id:        str,
        user_id:       Optional[str],
        input_type:    str,
        pathogen_name: str,
        status:        str,
        completed_at:  Optional[str] = None,
        timing:        Optional[Dict] = None,
        construct_report: Optional[Dict] = None,
    ) -> str:
        """
        Upsert a run row INSERT if new, UPDATE if exists.
        Called from main.py's finally block runs always regardless of
        whether pipeline succeeded or partially failed.
        """
        data: Dict[str, Any] = {
            "id":           run_id,
            "input_type":   input_type,
            "pathogen_name": pathogen_name,
            "status":       status,
        }
        if user_id:
            data["user_id"] = user_id
        if completed_at:
            data["completed_at"] = completed_at
        if timing:
            data["timing"] = timing
        if construct_report:
            data["construct_report"] = construct_report

        try:
            self.client.table("runs").upsert(data, on_conflict="id").execute()
            logger.info(f"Upserted run {run_id} status={status} user={user_id}")
        except Exception as e:
            # Fallback: try plain insert if upsert not supported on this schema version
            logger.warning(f"Upsert failed, trying insert: {e}")
            try:
                self.client.table("runs").insert(data).execute()
            except Exception as e2:
                # Last resort: try update
                logger.warning(f"Insert failed, trying update: {e2}")
                self.client.table("runs").update({
                    k: v for k, v in data.items() if k != "id"
                }).eq("id", run_id).execute()

        return run_id

    def create_run(self, run: PipelineRun, user_id: Optional[str] = None) -> str:
        """Legacy kept for compatibility. Prefer upsert_run()."""
        data = {
            "id":            run.run_id,
            "pathogen_name": run.pathogen_name,
            "input_type":    run.input_type,
            "raw_input":     run.raw_input,
            "current_stage": run.current_stage,
            "status":        run.status,
        }
        if user_id:
            data["user_id"] = user_id
        try:
            self.client.table("runs").upsert(data, on_conflict="id").execute()
        except Exception:
            self.client.table("runs").insert(data).execute()
        logger.info(f"Created run {run.run_id} (user={user_id})")
        return run.run_id

    def update_run_stage(self, run_id: str, stage: str, status: str = None):
        update_data: Dict[str, Any] = {"current_stage": stage}
        if status:
            update_data["status"] = status
            if status in ["completed", "failed"]:
                update_data["completed_at"] = datetime.now().isoformat()
        self.client.table("runs").update(update_data).eq("id", run_id).execute()

    def update_run_construct(self, run_id: str, construct_report: Dict, timing: Dict):
        try:
            self.client.table("runs").update({
                "construct_report": construct_report,
                "timing":           timing,
            }).eq("id", run_id).execute()
            logger.info(f"Persisted construct + timing for run {run_id}")
        except Exception as e:
            logger.warning(f"Could not persist construct/timing: {e}")

    def get_run(self, run_id: str) -> Optional[Dict]:
        result = self.client.table("runs").select("*").eq("id", run_id).execute()
        return result.data[0] if result.data else None

    # ── CANDIDATES ────────────────────────────────────────────────────────────

    def save_candidate(
        self, run_id: str, candidate: CandidateProtein, user_id: Optional[str] = None
    ) -> str:
        try:
            existing = (
                self.client.table("candidates")
                .select("id")
                .eq("run_id", run_id)
                .eq("protein_id", candidate.protein_id)
                .execute()
            )
            data: Dict[str, Any] = {
                "run_id":                 run_id,
                "protein_id":             candidate.protein_id,
                "protein_name":           candidate.protein_name,
                "sequence":               candidate.sequence,
                "source":                 candidate.source,
                "stage":                  candidate.stage,
                "status":                 candidate.status.value,
                "confidence_tier":        candidate.confidence_tier.value,
                "flags":                  candidate.flags,
                "psortb_localization":    candidate.psortb_localization,
                "tmhmm_helices":          candidate.tmhmm_helices,
                "vaxijen_score":          candidate.vaxijen_score,
                "blast_human_identity":   candidate.blast_human_identity,
                "structure_source":       candidate.structure_source,
                "structure_pdb_path":     candidate.structure_pdb_path,
                "conservation_profile":   candidate.conservation_profile,
                "hla_coverage_global":    candidate.hla_coverage_global,
                "hla_coverage_africa":    candidate.hla_coverage_africa,
                "decisions":              candidate.decisions,
                "updated_at":             datetime.now().isoformat(),
            }
            if user_id:
                data["user_id"] = user_id

            if existing.data:
                self.client.table("candidates").update(data).eq("id", existing.data[0]["id"]).execute()
                candidate_db_id = existing.data[0]["id"]
            else:
                result = self.client.table("candidates").insert(data).execute()
                candidate_db_id = result.data[0]["id"]

            self._save_epitopes(candidate_db_id, run_id, candidate, user_id=user_id)
            logger.info(f"Saved candidate {candidate.protein_id}")
            return candidate_db_id

        except Exception as e:
            logger.error(f"Failed to save candidate {candidate.protein_id}: {e}")
            raise

    def get_candidates_for_run(self, run_id: str) -> List[Dict]:
        result = self.client.table("candidates").select("*").eq("run_id", run_id).execute()
        return result.data or []

    # ── EPITOPES ──────────────────────────────────────────────────────────────

    def _save_epitopes(
        self,
        candidate_id: str,
        run_id:       str,
        candidate:    CandidateProtein,
        user_id:      Optional[str] = None,
    ):
        self.client.table("epitopes").delete().eq("candidate_id", candidate_id).execute()

        all_epitopes = candidate.ctl_epitopes + candidate.htl_epitopes + candidate.bcell_epitopes
        if not all_epitopes:
            return

        rows = []
        for ep in all_epitopes:
            tool_outputs = ep.tool_outputs or {}
            row: Dict[str, Any] = {
                "candidate_id":     candidate_id,
                "run_id":           run_id,
                "sequence":         ep.sequence,
                "epitope_type":     ep.epitope_type.value,
                "hla_allele":       ep.hla_allele,
                "ic50_nm":          ep.ic50_nm,
                "percentile_rank":  ep.percentile_rank,
                "conservation_score": ep.conservation_score,
                "allergenicity_safe": ep.allergenicity_safe,
                "toxicity_safe":    ep.toxicity_safe,
                "confidence_tier":  ep.confidence_tier.value,
                "tool_outputs":     tool_outputs,
            }
            if user_id:
                row["user_id"] = user_id
            rows.append(row)

        # Batch insert in chunks to avoid payload limits
        chunk_size = 100
        for i in range(0, len(rows), chunk_size):
            self.client.table("epitopes").insert(rows[i:i+chunk_size]).execute()

        logger.info(f"Saved {len(rows)} epitopes for candidate {candidate_id}")

    def get_epitopes_for_candidate(self, candidate_id: str) -> List[Dict]:
        result = (
            self.client.table("epitopes")
            .select("*")
            .eq("candidate_id", candidate_id)
            .execute()
        )
        return result.data or []

    # ── DECISIONS ─────────────────────────────────────────────────────────────

    def log_decision(
        self, run_id: str, candidate_id: str, stage: str,
        decision: str, reasoning: str,
        input_data: Dict[str, Any] = None, user_id: Optional[str] = None,
    ):
        try:
            row: Dict[str, Any] = {
                "run_id":       run_id,
                "candidate_id": candidate_id,
                "stage":        stage,
                "decision":     decision,
                "reasoning":    reasoning,
                "input_data":   input_data or {},
            }
            if user_id:
                row["user_id"] = user_id
            self.client.table("decisions").insert(row).execute()
        except Exception as e:
            logger.error(f"Failed to log decision: {e}")

    # ── DIAGNOSTICS ───────────────────────────────────────────────────────────

    def test_connection(self) -> bool:
        try:
            self.client.table("runs").select("id").limit(1).execute()
            return True
        except Exception as e:
            logger.error(f"DB connection failed: {e}")
            return False


db = SupabaseClient()