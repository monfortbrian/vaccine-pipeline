import os
from typing import List, Optional, Dict, Any
from supabase import create_client, Client
from src.models.candidate import PipelineRun, CandidateProtein, EpitopeResult
import json
from datetime import datetime
import logging

if os.getenv("RAILWAY_ENVIRONMENT") is None:
    from dotenv import load_dotenv
    load_dotenv()

logger = logging.getLogger(__name__)

class SupabaseClient:
    def __init__(self):
        self.url = os.getenv("SUPABASE_URL")
        self.key = os.getenv("SUPABASE_KEY")

        if not self.url or not self.key:
            raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in environment")

        self.client: Client = create_client(self.url, self.key)

    def create_run(self, run: PipelineRun) -> str:
        """Create a new pipeline run in the database."""
        try:
            result = self.client.table("runs").insert({
                "id": run.run_id,
                "pathogen_name": run.pathogen_name,
                "input_type": run.input_type,
                "raw_input": run.raw_input,
                "target_populations": run.target_populations,
                "coverage_threshold": run.coverage_threshold,
                "max_candidates_output": run.max_candidates_output,
                "current_stage": run.current_stage,
                "status": run.status,
                "config": {
                    "coverage_loop_count": run.coverage_loop_count,
                    "coverage_met": run.coverage_met
                }
            }).execute()

            logger.info(f"Created run {run.run_id} for pathogen {run.pathogen_name}")
            return run.run_id

        except Exception as e:
            logger.error(f"Failed to create run: {e}")
            raise

    def update_run_stage(self, run_id: str, stage: str, status: str = None):
        """Update the current stage of a pipeline run."""
        update_data = {"current_stage": stage}
        if status:
            update_data["status"] = status
            if status in ["completed", "failed"]:
                update_data["completed_at"] = datetime.now().isoformat()

        self.client.table("runs").update(update_data).eq("id", run_id).execute()
        logger.info(f"Updated run {run_id} to stage {stage}")

    def save_candidate(self, run_id: str, candidate: CandidateProtein) -> str:
        """Save or update a candidate protein."""
        try:
            # Check if candidate exists
            existing = self.client.table("candidates").select("id").eq("run_id", run_id).eq("protein_id", candidate.protein_id).execute()

            candidate_data = {
                "run_id": run_id,
                "protein_id": candidate.protein_id,
                "protein_name": candidate.protein_name,
                "sequence": candidate.sequence,
                "source": candidate.source,
                "stage": candidate.stage,
                "status": candidate.status.value,
                "confidence_tier": candidate.confidence_tier.value,
                "flags": candidate.flags,
                "psortb_localization": candidate.psortb_localization,
                "tmhmm_helices": candidate.tmhmm_helices,
                "vaxijen_score": candidate.vaxijen_score,
                "blast_human_identity": candidate.blast_human_identity,
                "structure_source": candidate.structure_source,
                "structure_pdb_path": candidate.structure_pdb_path,
                "conservation_profile": candidate.conservation_profile,
                "hla_coverage_global": candidate.hla_coverage_global,
                "hla_coverage_africa": candidate.hla_coverage_africa,
                "decisions": candidate.decisions,
                "updated_at": datetime.now().isoformat()
            }

            if existing.data:
                # Update existing
                result = self.client.table("candidates").update(candidate_data).eq("id", existing.data[0]["id"]).execute()
                candidate_db_id = existing.data[0]["id"]
            else:
                # Insert new
                result = self.client.table("candidates").insert(candidate_data).execute()
                candidate_db_id = result.data[0]["id"]

            # Save epitopes
            self._save_epitopes(candidate_db_id, run_id, candidate)

            logger.info(f"Saved candidate {candidate.protein_id}")
            return candidate_db_id

        except Exception as e:
            logger.error(f"Failed to save candidate {candidate.protein_id}: {e}")
            raise

    def _save_epitopes(self, candidate_id: str, run_id: str, candidate: CandidateProtein):
        """Save epitopes for a candidate."""
        # Delete existing epitopes for this candidate
        self.client.table("epitopes").delete().eq("candidate_id", candidate_id).execute()

        all_epitopes = candidate.ctl_epitopes + candidate.htl_epitopes + candidate.bcell_epitopes

        if not all_epitopes:
            return

        epitope_data = []
        for epitope in all_epitopes:
            epitope_data.append({
                "candidate_id": candidate_id,
                "run_id": run_id,
                "sequence": epitope.sequence,
                "epitope_type": epitope.epitope_type.value,
                "hla_allele": epitope.hla_allele,
                "ic50_nm": epitope.ic50_nm,
                "percentile_rank": epitope.percentile_rank,
                "conservation_score": epitope.conservation_score,
                "allergenicity_safe": epitope.allergenicity_safe,
                "toxicity_safe": epitope.toxicity_safe,
                "confidence_tier": epitope.confidence_tier.value,
                "tool_outputs": epitope.tool_outputs
            })

        self.client.table("epitopes").insert(epitope_data).execute()
        logger.info(f"Saved {len(epitope_data)} epitopes for candidate {candidate_id}")

    def log_decision(self, run_id: str, candidate_id: str, stage: str, decision: str, reasoning: str, input_data: Dict[str, Any] = None):
        """Log an AI decision for audit trail."""
        try:
            self.client.table("decisions").insert({
                "run_id": run_id,
                "candidate_id": candidate_id,
                "stage": stage,
                "decision": decision,
                "reasoning": reasoning,
                "input_data": input_data or {}
            }).execute()

            logger.info(f"Logged decision: {stage} -> {decision}")

        except Exception as e:
            logger.error(f"Failed to log decision: {e}")
            raise

    def get_run(self, run_id: str) -> Optional[Dict]:
        """Get a pipeline run by ID."""
        result = self.client.table("runs").select("*").eq("id", run_id).execute()
        return result.data[0] if result.data else None

    def get_candidates_for_run(self, run_id: str) -> List[Dict]:
        """Get all candidates for a run."""
        result = self.client.table("candidates").select("*").eq("run_id", run_id).execute()
        return result.data

    def get_epitopes_for_candidate(self, candidate_id: str) -> List[Dict]:
        """Get all epitopes for a candidate."""
        result = self.client.table("epitopes").select("*").eq("candidate_id", candidate_id).execute()
        return result.data

    def test_connection(self) -> bool:
        """Test database connection."""
        try:
            result = self.client.table("runs").select("id").limit(1).execute()
            logger.info("Database connection successful")
            return True
        except Exception as e:
            logger.error(f"Database connection failed: {e}")
            return False

# Global instance
db = SupabaseClient()