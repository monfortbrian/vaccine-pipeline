"""
T-Cell Predictor Agent (Node N3)
Predicts CTL and HTL epitopes using IEDB APIs.
"""

import logging
from typing import List, Dict, Any
from src.models.candidate import CandidateProtein, EpitopeResult, EpitopeType, ConfidenceTier

logger = logging.getLogger(__name__)


class TCellPredictorAgent:
    def __init__(self):
        self.stage_name = "tcell_prediction"
        self._iedb = None

    @property
    def iedb(self):
        if self._iedb is None:
            from src.tools.iedb_client import iedb
            self._iedb = iedb
        return self._iedb

    def run(self, candidates: List[CandidateProtein]) -> List[CandidateProtein]:
        """Run T-cell epitope prediction."""
        logger.info("Starting T-cell epitope prediction")

        active_candidates = [c for c in candidates if c.status.value == "active"]

        for i, candidate in enumerate(active_candidates):
            logger.info(f"Predicting epitopes {i+1}/{len(active_candidates)}: {candidate.protein_name}")

            try:
                # Predict CTL epitopes
                ctl_predictions = self.iedb.predict_mhc_i_binding(candidate.sequence)
                candidate.ctl_epitopes = self._process_ctl_predictions(ctl_predictions)

                # Predict HTL epitopes
                htl_predictions = self.iedb.predict_mhc_ii_binding(candidate.sequence)
                candidate.htl_epitopes = self._process_htl_predictions(htl_predictions)

                candidate.stage = self.stage_name

                ctl_count = len([e for e in candidate.ctl_epitopes if e.confidence_tier != ConfidenceTier.UNCERTAIN])
                htl_count = len([e for e in candidate.htl_epitopes if e.confidence_tier != ConfidenceTier.UNCERTAIN])

                candidate.add_decision(
                    stage=self.stage_name,
                    decision="epitopes_predicted",
                    reasoning=f"Predicted {ctl_count} CTL and {htl_count} HTL epitopes",
                    ctl_epitope_count=ctl_count,
                    htl_epitope_count=htl_count
                )

                logger.info(f" CTL: {ctl_count}, HTL: {htl_count}")

            except Exception as e:
                logger.error(f" T-cell prediction failed: {e}")
                candidate.flags.append("tcell_prediction_failed")

        return candidates

    def _process_ctl_predictions(self, predictions: List[Dict[str, Any]]) -> List[EpitopeResult]:
        epitopes = []
        for pred in predictions:
            try:
                ic50 = pred.get('ic50_nm', 50000)
                if ic50 > 5000:
                    continue

                epitope = EpitopeResult(
                    sequence=pred['sequence'],
                    epitope_type=EpitopeType.CTL,
                    hla_allele=pred['allele'],
                    ic50_nm=ic50,
                    percentile_rank=pred.get('percentile_rank'),
                    confidence_tier=self._score_ctl_confidence(pred),
                    tool_outputs=pred
                )
                epitopes.append(epitope)
            except Exception as e:
                logger.warning(f"Failed to process CTL: {e}")

        epitopes.sort(key=lambda x: x.ic50_nm or 50000)
        return epitopes[:20]

    def _process_htl_predictions(self, predictions: List[Dict[str, Any]]) -> List[EpitopeResult]:
        epitopes = []
        for pred in predictions:
            try:
                ic50 = pred.get('ic50_nm', 50000)
                if ic50 > 10000:
                    continue

                epitope = EpitopeResult(
                    sequence=pred['sequence'],
                    epitope_type=EpitopeType.HTL,
                    hla_allele=pred['allele'],
                    ic50_nm=ic50,
                    percentile_rank=pred.get('percentile_rank'),
                    confidence_tier=self._score_htl_confidence(pred),
                    tool_outputs=pred
                )
                epitopes.append(epitope)
            except Exception as e:
                logger.warning(f"Failed to process HTL: {e}")

        epitopes.sort(key=lambda x: x.ic50_nm or 50000)
        return epitopes[:15]

    def _score_ctl_confidence(self, prediction: Dict[str, Any]) -> ConfidenceTier:
        """CTL confidence based on percentile rank (preferred) or IC50."""
        rank = prediction.get('percentile_rank')
        if rank is not None:
            if rank < 0.5:
                return ConfidenceTier.HIGH
            elif rank < 2.0:
                return ConfidenceTier.MEDIUM
            elif rank < 10.0:
                return ConfidenceTier.LOW
            else:
                return ConfidenceTier.UNCERTAIN

        ic50 = prediction.get('ic50_nm', 50000)
        if ic50 < 50:
            return ConfidenceTier.HIGH
        elif ic50 < 500:
            return ConfidenceTier.MEDIUM
        elif ic50 < 5000:
            return ConfidenceTier.LOW
        else:
            return ConfidenceTier.UNCERTAIN

    def _score_htl_confidence(self, prediction: Dict[str, Any]) -> ConfidenceTier:
        """
        HTL confidence based on percentile rank.

        MHC-II binding thresholds are DIFFERENT from MHC-I:
        - MHC-II percentile ranks are naturally higher
        - rank < 2.0 = strong binder (HIGH)
        - rank < 5.0 = good binder (MEDIUM)
        - rank < 10.0 = weak binder (LOW)

        This is the standard IEDB recommended threshold for MHC-II.
        """
        rank = prediction.get('percentile_rank')
        if rank is not None:
            if rank < 2.0:
                return ConfidenceTier.HIGH
            elif rank < 5.0:
                return ConfidenceTier.MEDIUM
            elif rank < 10.0:
                return ConfidenceTier.LOW
            else:
                return ConfidenceTier.UNCERTAIN

        ic50 = prediction.get('ic50_nm', 50000)
        if ic50 < 500:
            return ConfidenceTier.HIGH
        elif ic50 < 2000:
            return ConfidenceTier.MEDIUM
        elif ic50 < 10000:
            return ConfidenceTier.LOW
        else:
            return ConfidenceTier.UNCERTAIN


# Global instance
tcell_predictor = TCellPredictorAgent()