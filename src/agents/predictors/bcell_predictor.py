"""
B-CELL PREDICTOR AGENT - MVP-2 NODE N4
Predicts antibody epitopes using IEDB BepiPred 2.0.
"""

import logging
from typing import List, Dict, Any
from src.models.candidate import CandidateProtein, EpitopeResult, EpitopeType, ConfidenceTier

logger = logging.getLogger(__name__)


class BCellPredictorAgent:
    def __init__(self):
        self.stage_name = "bcell_prediction"
        self._iedb_bcell = None

    @property
    def iedb_bcell(self):
        if self._iedb_bcell is None:
            from src.tools.iedb_bcell_client import iedb_bcell
            self._iedb_bcell = iedb_bcell
        return self._iedb_bcell

    def run(self, candidates: List[CandidateProtein]) -> List[CandidateProtein]:
        """Run B-cell epitope prediction on candidates with T-cell epitopes."""
        logger.info("Starting B-cell epitope prediction")

        # Only predict for candidates that have T-cell epitopes
        tcell_candidates = [c for c in candidates
                            if c.status.value == "active" and
                            (len(c.ctl_epitopes) > 0 or len(c.htl_epitopes) > 0)]

        logger.info(
            f"   Processing {len(tcell_candidates)} candidates with T-cell epitopes")

        for i, candidate in enumerate(tcell_candidates):
            logger.info(
                f"   Predicting B-cell epitopes {i+1}/{len(tcell_candidates)}: {candidate.protein_name}")

            try:
                # Use consensus method for better accuracy
                iedb_predictions = self.iedb_bcell.predict_consensus_epitopes(
                    candidate.sequence)

                # Convert to EpitopeResult objects
                bcell_epitopes = self._process_iedb_predictions(
                    iedb_predictions)

                candidate.bcell_epitopes = bcell_epitopes
                candidate.stage = self.stage_name

                # Add decision record
                epitope_count = len(bcell_epitopes)
                high_conf_count = len(
                    [e for e in bcell_epitopes if e.confidence_tier == ConfidenceTier.HIGH])

                candidate.add_decision(
                    stage=self.stage_name,
                    decision="epitopes_predicted",
                    reasoning=f"Predicted {epitope_count} B-cell epitopes ({high_conf_count} high confidence)",
                    bcell_epitope_count=epitope_count,
                    high_confidence_count=high_conf_count
                )

                logger.info(
                    f"       B-cell epitopes: {epitope_count} ({high_conf_count} high conf)")

            except Exception as e:
                logger.error(f"       B-cell prediction failed: {e}")
                candidate.flags.append("bcell_prediction_failed")
                candidate.add_decision(
                    stage=self.stage_name,
                    decision="prediction_failed",
                    reasoning=f"Error: {str(e)}"
                )

        total_bcell = sum(len(c.bcell_epitopes) for c in tcell_candidates)
        logger.info(
            f"B-cell prediction complete: {total_bcell} total epitopes")

        return candidates

    def _process_iedb_predictions(self, predictions: List[Dict[str, Any]]) -> List[EpitopeResult]:
        """Process IEDB B-cell predictions into EpitopeResult objects."""
        epitopes = []

        for pred in predictions:
            try:
                confidence_str = pred.get('confidence', 'uncertain')
                confidence_tier = self._map_confidence(confidence_str)

                epitope = EpitopeResult(
                    sequence=pred['sequence'],
                    epitope_type=EpitopeType.B_CELL_LINEAR,
                    hla_allele=None,  # B-cell epitopes don't use HLA
                    ic50_nm=None,     # Different binding mechanism
                    percentile_rank=None,
                    confidence_tier=confidence_tier,
                    tool_outputs=pred
                )

                epitopes.append(epitope)

            except Exception as e:
                logger.warning(f"Failed to process B-cell prediction: {e}")

        return epitopes

    def _map_confidence(self, confidence_str: str) -> ConfidenceTier:
        """Map string confidence to ConfidenceTier."""
        confidence_map = {
            'high': ConfidenceTier.HIGH,
            'medium': ConfidenceTier.MEDIUM,
            'low': ConfidenceTier.LOW,
            'uncertain': ConfidenceTier.UNCERTAIN
        }
        return confidence_map.get(confidence_str.lower(), ConfidenceTier.UNCERTAIN)


# Global instance
bcell_predictor = BCellPredictorAgent()
