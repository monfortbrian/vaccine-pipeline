"""
BCELL PREDICTOR AGENT

Tool: IEDB BepiPred 2.0 (linear B-cell epitope prediction)
Reference: Jespersen et al. (2017) Nucleic Acids Res 45:W24-W29

B-cell epitopes have no HLA restriction, they are recognised by antibodies
directly on the protein surface. No allele field, no IC50.

Animal model note:
  Rabbit is the standard model for B-cell epitope validation.
  Rabbit polyclonal antisera against linear BepiPred epitopes is a
  standard first-pass validation protocol in vaccine immunology.
  All B-cell epitopes carry model_categories: ["RABBIT"] to signal
  this validation path to the wet-lab team.

Confidence scoring:
  HIGH   : BepiPred score >= 0.9 (very strong linear epitope signal)
  MEDIUM : BepiPred score >= 0.7
  LOW    : BepiPred score >= 0.5 (default BepiPred threshold)
  UNCERTAIN : score < 0.5 or unavailable
"""

import logging
from typing import List, Dict, Any
from src.models.candidate import CandidateProtein, EpitopeResult, EpitopeType, ConfidenceTier

logger = logging.getLogger("tope_deep.agents.N4")

BEPIPRED_THRESHOLD = 0.5


class BCellPredictorAgent:
    def __init__(self):
        self.stage_name = "bcell_prediction"
        self._iedb_bcell = None

    @property
    def iedb_bcell(self):
        if not self._iedb_bcell:
            from src.tools.iedb_bcell_client import iedb_bcell
            self._iedb_bcell = iedb_bcell
        return self._iedb_bcell

    def run(self, candidates: List[CandidateProtein]) -> List[CandidateProtein]:
        logger.info("N4: B-Cell Predictor, IEDB BepiPred 2.0")
        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"   Processing {len(active)} candidates with T-cell epitopes")

        for i, candidate in enumerate(active):
            logger.info(
                f"   [{i+1}/{len(active)}] {candidate.protein_name}"
            )
            try:
                raw = self.iedb_bcell.predict_bcell_epitopes(candidate.sequence)
                candidate.bcell_epitopes = self._process(raw)
                candidate.stage = self.stage_name

                high = sum(1 for e in candidate.bcell_epitopes if e.confidence_tier == ConfidenceTier.HIGH)

                candidate.add_decision(
                    stage=self.stage_name,
                    decision="bcell_predicted",
                    reasoning=(
                        f"B-cell epitopes: {len(candidate.bcell_epitopes)} total "
                        f"({high} high confidence, BepiPred score >= 0.9). "
                        f"Method: IEDB BepiPred 2.0 (Jespersen et al. 2017). "
                        f"All B-cell epitopes flagged for rabbit polyclonal antisera validation "
                        f"standard first-pass protocol for linear epitope confirmation. "
                        f"No HLA restriction applies to B-cell epitopes."
                    ),
                    bcell_count=len(candidate.bcell_epitopes),
                    bcell_high_confidence=high,
                    method="IEDB_BepiPred_2.0",
                    rabbit_validation_path=True,
                )
                logger.info(
                    f"      B-cell: {len(candidate.bcell_epitopes)} "
                    f"({high} high conf) [IEDB_BepiPred_2.0]"
                )
            except Exception as e:
                logger.error(f"      N4 failed for {candidate.protein_name}: {e}")
                candidate.add_decision(
                    stage=self.stage_name,
                    decision="prediction_failed",
                    reasoning=f"BepiPred 2.0 call failed: {str(e)}. No B-cell epitopes predicted.",
                )

        logger.info("N4 complete")
        return candidates

    def _process(self, predictions: List[Dict[str, Any]]) -> List[EpitopeResult]:
        epitopes = []
        seen: set = set()

        for pred in predictions:
            try:
                seq   = pred.get("sequence", "")
                score = pred.get("score", 0.0)

                if score < BEPIPRED_THRESHOLD or seq in seen or len(seq) < 8:
                    continue

                seen.add(seq)
                epitopes.append(EpitopeResult(
                    sequence=seq,
                    epitope_type=EpitopeType.BCELL,
                    hla_allele=None,
                    ic50_nm=None,
                    percentile_rank=None,
                    confidence_tier=self._score(score),
                    tool_outputs={
                        **pred,
                        "method_used":       "IEDB_BepiPred_2.0",
                        "bepipred_score":    score,
                        # Rabbit model flag, B-cell validation standard
                        "model_categories":  ["RABBIT"],
                        "rabbit_validation": True,
                        "validation_note":   (
                            "Rabbit polyclonal antisera standard protocol. "
                            "Jespersen et al. (2017) Nucleic Acids Res 45:W24-W29"
                        ),
                    },
                ))
            except Exception as e:
                logger.warning(f"      B-cell process error: {e}")

        epitopes.sort(
            key=lambda e: e.tool_outputs.get("bepipred_score", 0),
            reverse=True,
        )
        return epitopes[:10]

    @staticmethod
    def _score(score: float) -> ConfidenceTier:
        if score >= 0.9:  return ConfidenceTier.HIGH
        if score >= 0.7:  return ConfidenceTier.MEDIUM
        if score >= 0.5:  return ConfidenceTier.LOW
        return ConfidenceTier.UNCERTAIN


bcell_predictor = BCellPredictorAgent()