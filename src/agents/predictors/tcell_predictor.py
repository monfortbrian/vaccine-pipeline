"""
T-CELL PREDICTOR AGENT — TOPE_DEEP NODE N3
Predicts CTL (MHC-I) and HTL (MHC-II) epitopes.

Primary:  IEDB tools-cluster-interface (NetMHCpan 4.1 EL, NetMHCIIpan 4.3)
Fallback: MHCflurry 2.0 for MHC-I when IEDB is down (O'Brien et al. 2019)
          No fallback for MHC-II — logged explicitly in audit trail.

IC50 note:
  Values are APPROXIMATED from percentile rank using the IEDB standard
  rank-to-IC50 mapping (Sette & Sidney 1999). They are not measured IC50 values.
  The field ic50_nm should be read as "binding affinity estimate (nM)".
  method_used field records this on every epitope.

Confidence thresholds:
  CTL (MHC-I):  rank < 0.5 = HIGH, < 2.0 = MEDIUM, < 10.0 = LOW
  HTL (MHC-II): rank < 2.0 = HIGH, < 5.0 = MEDIUM, < 10.0 = LOW
  MHC-II thresholds differ — naturally higher ranks due to longer binding groove.
  Reference: IEDB recommended thresholds (tools.iedb.org/mhcii/help/).
"""

import logging
from typing import List, Dict, Any, Optional
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
        logger.info("N3: Starting T-cell epitope prediction")
        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"   {len(active)} candidates")

        for i, candidate in enumerate(active):
            logger.info(
                f"   [{i+1}/{len(active)}] {candidate.protein_name} "
                f"({len(candidate.sequence)} aa)"
            )
            try:
                # CTL — MHC-I
                ctl_raw = self.iedb.predict_mhc_i_binding(candidate.sequence)
                candidate.ctl_epitopes = self._process_ctl(ctl_raw)

                # HTL — MHC-II
                htl_raw = self.iedb.predict_mhc_ii_binding(candidate.sequence)
                candidate.htl_epitopes = self._process_htl(htl_raw)

                candidate.stage = self.stage_name

                ctl_high = len([e for e in candidate.ctl_epitopes
                                if e.confidence_tier == ConfidenceTier.HIGH])
                htl_high = len([e for e in candidate.htl_epitopes
                                if e.confidence_tier == ConfidenceTier.HIGH])

                # Detect which prediction method was actually used
                ctl_method = _infer_method(ctl_raw, "CTL")
                htl_method = _infer_method(htl_raw, "HTL")
                htl_failed = len(htl_raw) == 0

                candidate.add_decision(
                    stage=self.stage_name,
                    decision="epitopes_predicted",
                    reasoning=(
                        f"CTL: {len(candidate.ctl_epitopes)} epitopes "
                        f"({ctl_high} high confidence, IC50 <500 nM). "
                        f"HTL: {len(candidate.htl_epitopes)} epitopes "
                        f"({htl_high} high confidence). "
                        + ("HTL prediction returned no results — "
                           "NetMHCIIpan may be unavailable. "
                           "No MHC-II fallback exists. "
                           "Population coverage will be MHC-I only for this candidate. "
                           if htl_failed else "") +
                        f"CTL method: {ctl_method}. "
                        f"HTL method: {htl_method}. "
                        "IC50 values are approximated from percentile rank "
                        "(IEDB rank-to-IC50 mapping, Sette & Sidney 1999). "
                        "Not equivalent to measured IC50."
                    ),
                    ctl_count=len(candidate.ctl_epitopes),
                    ctl_high_confidence=ctl_high,
                    htl_count=len(candidate.htl_epitopes),
                    htl_high_confidence=htl_high,
                    ctl_method=ctl_method,
                    htl_method=htl_method,
                    htl_failed=htl_failed,
                )

                logger.info(
                    f"      CTL: {len(candidate.ctl_epitopes)} "
                    f"({ctl_high} high) [{ctl_method}]"
                )
                logger.info(
                    f"      HTL: {len(candidate.htl_epitopes)} "
                    f"({htl_high} high) [{htl_method}]"
                    + (" — WARNING: no HTL predictions" if htl_failed else "")
                )

            except Exception as e:
                logger.error(f"      N3 failed for {candidate.protein_name}: {e}")
                candidate.flags.append("tcell_prediction_failed")
                candidate.add_decision(
                    stage=self.stage_name,
                    decision="prediction_failed",
                    reasoning=f"N3 exception: {str(e)}. No epitopes predicted.",
                )

        logger.info("N3: T-cell prediction complete")
        return candidates

    # ── PROCESSORS ────────────────────────────────────────────────────────────

    def _process_ctl(self, predictions: List[Dict[str, Any]]) -> List[EpitopeResult]:
        epitopes = []
        for pred in predictions:
            try:
                ic50 = pred.get("ic50_nm", 50000)
                if ic50 > 5000:
                    continue
                epitopes.append(EpitopeResult(
                    sequence=pred["sequence"],
                    epitope_type=EpitopeType.CTL,
                    hla_allele=pred["allele"],
                    ic50_nm=ic50,
                    percentile_rank=pred.get("percentile_rank"),
                    confidence_tier=self._score_ctl(pred),
                    tool_outputs={
                        **pred,
                        "ic50_note": "approximated_from_percentile_rank",
                        "method_used": pred.get(
                            "prediction_method", "IEDB_NetMHCpan4.1"
                        ),
                    },
                ))
            except Exception as e:
                logger.warning(f"CTL process error: {e}")
        epitopes.sort(key=lambda x: x.ic50_nm or 50000)
        return epitopes[:20]

    def _process_htl(self, predictions: List[Dict[str, Any]]) -> List[EpitopeResult]:
        epitopes = []
        for pred in predictions:
            try:
                ic50 = pred.get("ic50_nm", 50000)
                if ic50 > 10000:
                    continue
                epitopes.append(EpitopeResult(
                    sequence=pred["sequence"],
                    epitope_type=EpitopeType.HTL,
                    hla_allele=pred["allele"],
                    ic50_nm=ic50,
                    percentile_rank=pred.get("percentile_rank"),
                    confidence_tier=self._score_htl(pred),
                    tool_outputs={
                        **pred,
                        "ic50_note": "approximated_from_percentile_rank",
                        "method_used": pred.get(
                            "prediction_method", "IEDB_NetMHCIIpan4.3"
                        ),
                    },
                ))
            except Exception as e:
                logger.warning(f"HTL process error: {e}")
        epitopes.sort(key=lambda x: x.ic50_nm or 50000)
        return epitopes[:15]

    # ── CONFIDENCE SCORING ─────────────────────────────────────────────────────

    def _score_ctl(self, pred: Dict[str, Any]) -> ConfidenceTier:
        """
        MHC-I thresholds (NetMHCpan 4.1 EL):
          rank < 0.5 = strong binder = HIGH
          rank < 2.0 = weak binder   = MEDIUM
          rank < 10  = moderate      = LOW
        Reference: Reynisson et al., Nucleic Acids Research 2020.
        """
        rank = pred.get("percentile_rank")
        if rank is not None:
            if rank < 0.5:  return ConfidenceTier.HIGH
            if rank < 2.0:  return ConfidenceTier.MEDIUM
            if rank < 10.0: return ConfidenceTier.LOW
            return ConfidenceTier.UNCERTAIN
        ic50 = pred.get("ic50_nm", 50000)
        if ic50 < 50:    return ConfidenceTier.HIGH
        if ic50 < 500:   return ConfidenceTier.MEDIUM
        if ic50 < 5000:  return ConfidenceTier.LOW
        return ConfidenceTier.UNCERTAIN

    def _score_htl(self, pred: Dict[str, Any]) -> ConfidenceTier:
        """
        MHC-II thresholds (NetMHCIIpan 4.3):
          rank < 2.0 = strong = HIGH
          rank < 5.0 = good   = MEDIUM
          rank < 10  = weak   = LOW
        MHC-II ranks are naturally higher than MHC-I.
        Reference: IEDB MHC-II help (tools.iedb.org/mhcii/help/).
        """
        rank = pred.get("percentile_rank")
        if rank is not None:
            if rank < 2.0:  return ConfidenceTier.HIGH
            if rank < 5.0:  return ConfidenceTier.MEDIUM
            if rank < 10.0: return ConfidenceTier.LOW
            return ConfidenceTier.UNCERTAIN
        ic50 = pred.get("ic50_nm", 50000)
        if ic50 < 500:   return ConfidenceTier.HIGH
        if ic50 < 2000:  return ConfidenceTier.MEDIUM
        if ic50 < 10000: return ConfidenceTier.LOW
        return ConfidenceTier.UNCERTAIN


def _infer_method(predictions: List[Dict], epitope_class: str) -> str:
    """Infer which prediction method was used from the first prediction dict."""
    if not predictions:
        return f"none_{epitope_class.lower()}_unavailable"
    method = predictions[0].get("prediction_method", "")
    if "MHCflurry" in method:
        return "MHCflurry_2.0_affinity_fallback"
    if "NetMHCpan" in method or "IEDB" in method:
        return "IEDB_NetMHCpan4.1_EL"
    if "NetMHCIIpan" in method:
        return "IEDB_NetMHCIIpan4.3"
    return method or "unknown"


tcell_predictor = TCellPredictorAgent()