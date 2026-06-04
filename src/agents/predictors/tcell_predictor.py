"""
TOPE_DEEP NODE N3
NetMHCpan 4.1 (MHC-I) + NetMHCIIpan 4.3 (MHC-II) via IEDB tools cluster.
MHCflurry 2.0 fallback for MHC-I when IEDB is unavailable.

Animal model alleles (Phase 1):
  Mouse H-2:    H-2-Kb, H-2-Db, H-2-Kd, H-2-Dd (C57BL/6 and BALB/c strains)
  Macaque Mamu: Mamu-A*01, Mamu-A*02, Mamu-B*17 (NHP model for TB/HIV/malaria)

Epitopes binding both human HLA and animal alleles are flagged as
wet-lab validation candidates testable in preclinical models.
"""

import logging
from typing import List, Dict, Any, Optional
from src.models.candidate import CandidateProtein, EpitopeResult, EpitopeType, ConfidenceTier

logger = logging.getLogger(__name__)

# Human HLA alleles used in N3 prediction
HUMAN_HLA_SUPERTYPES = [
    "HLA-A*02:01", "HLA-A*24:02", "HLA-A*03:01", "HLA-A*01:01",
    "HLA-A*11:01", "HLA-B*07:02", "HLA-B*44:02", "HLA-B*35:01",
]

# Animal model MHC alleles; Phase 1
MOUSE_H2_ALLELES = ["H-2-Kb", "H-2-Db", "H-2-Kd", "H-2-Dd"]
MACAQUE_MAMU_ALLELES = ["Mamu-A*01", "Mamu-A*02", "Mamu-B*17"]

# All alleles queried in one IEDB call
ALL_MHC_I_ALLELES = HUMAN_HLA_SUPERTYPES + MOUSE_H2_ALLELES + MACAQUE_MAMU_ALLELES


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
        logger.info("N3: T-cell epitope prediction (human HLA + mouse H-2 + macaque Mamu)")
        active = [c for c in candidates if c.status.value == "active"]

        for i, candidate in enumerate(active):
            logger.info(f"   [{i+1}/{len(active)}] {candidate.protein_name} ({len(candidate.sequence)} aa)")
            try:
                ctl_raw = self.iedb.predict_mhc_i_binding(candidate.sequence)
                candidate.ctl_epitopes = self._process_ctl(ctl_raw)

                htl_raw = self.iedb.predict_mhc_ii_binding(candidate.sequence)
                candidate.htl_epitopes = self._process_htl(htl_raw)

                candidate.stage = self.stage_name

                ctl_high = len([e for e in candidate.ctl_epitopes if e.confidence_tier == ConfidenceTier.HIGH])
                htl_high = len([e for e in candidate.htl_epitopes if e.confidence_tier == ConfidenceTier.HIGH])
                ctl_method = _infer_ctl_method(ctl_raw)
                htl_method = _infer_htl_method(htl_raw)
                htl_failed = len(htl_raw) == 0

                # Count animal model cross-reactive epitopes
                mouse_reactive = len([e for e in candidate.ctl_epitopes
                                     if e.tool_outputs.get("animal_model_alleles")])
                mamu_reactive  = len([e for e in candidate.ctl_epitopes
                                     if e.tool_outputs.get("mamu_alleles")])

                candidate.add_decision(
                    stage=self.stage_name,
                    decision="epitopes_predicted",
                    reasoning=(
                        f"CTL: {len(candidate.ctl_epitopes)} epitopes "
                        f"({ctl_high} high confidence). "
                        f"HTL: {len(candidate.htl_epitopes)} epitopes "
                        f"({htl_high} high confidence). "
                        f"Animal model: {mouse_reactive} mouse H-2 cross-reactive, "
                        f"{mamu_reactive} macaque Mamu cross-reactive. "
                        + ("HTL prediction unavailable no MHC-II fallback. " if htl_failed else "") +
                        f"CTL method: {ctl_method}. HTL method: {htl_method}. "
                        "IC50 approximated from percentile rank (Sette & Sidney 1999)."
                    ),
                    ctl_count=len(candidate.ctl_epitopes),
                    ctl_high_confidence=ctl_high,
                    htl_count=len(candidate.htl_epitopes),
                    htl_high_confidence=htl_high,
                    ctl_method=ctl_method,
                    htl_method=htl_method,
                    htl_failed=htl_failed,
                    mouse_h2_reactive=mouse_reactive,
                    mamu_reactive=mamu_reactive,
                )

                logger.info(f"      CTL: {len(candidate.ctl_epitopes)} ({ctl_high} high) [{ctl_method}]")
                logger.info(f"      HTL: {len(candidate.htl_epitopes)} ({htl_high} high) [{htl_method}]")
                logger.info(f"      Animal: {mouse_reactive} mouse H-2 | {mamu_reactive} Mamu cross-reactive")

            except Exception as e:
                logger.error(f"      N3 failed: {e}")
                candidate.flags.append("tcell_prediction_failed")
                candidate.add_decision(
                    stage=self.stage_name,
                    decision="prediction_failed",
                    reasoning=f"N3 exception: {str(e)}.",
                )

        return candidates

    def _process_ctl(self, predictions: List[Dict[str, Any]]) -> List[EpitopeResult]:
        epitopes = []
        seen_sequences = set()

        for pred in predictions:
            try:
                ic50 = pred.get("ic50_nm", 50000)
                if ic50 > 5000:
                    continue
                seq = pred["sequence"]

                # Classify allele species
                allele = pred.get("allele", "")
                is_human  = any(allele.startswith(h) for h in ["HLA-A", "HLA-B", "HLA-C"])
                is_mouse  = allele.startswith("H-2")
                is_mamu   = allele.startswith("Mamu")

                # Group predictions by sequence
                if seq not in seen_sequences:
                    seen_sequences.add(seq)
                    ep = EpitopeResult(
                        sequence=seq,
                        epitope_type=EpitopeType.CTL,
                        hla_allele=allele if is_human else None,
                        ic50_nm=ic50,
                        percentile_rank=pred.get("percentile_rank"),
                        confidence_tier=self._score_ctl(pred),
                        tool_outputs={
                            **pred,
                            "ic50_note": "approximated_from_percentile_rank",
                            "method_used": _infer_ctl_method([pred]),
                            "animal_model_alleles": [allele] if (is_mouse or is_mamu) else [],
                            "mamu_alleles": [allele] if is_mamu else [],
                            "human_hla_alleles": [allele] if is_human else [],
                        },
                    )
                    epitopes.append(ep)
                else:
                    # Add allele to existing epitope
                    for ep in epitopes:
                        if ep.sequence == seq:
                            if is_human and allele not in ep.tool_outputs.get("human_hla_alleles", []):
                                ep.tool_outputs["human_hla_alleles"].append(allele)
                                if ep.hla_allele is None:
                                    ep.hla_allele = allele
                            if is_mouse and allele not in ep.tool_outputs.get("animal_model_alleles", []):
                                ep.tool_outputs["animal_model_alleles"].append(allele)
                            if is_mamu and allele not in ep.tool_outputs.get("mamu_alleles", []):
                                ep.tool_outputs["mamu_alleles"].append(allele)
                                ep.tool_outputs["animal_model_alleles"].append(allele)
                            break

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
                    hla_allele=pred.get("allele"),
                    ic50_nm=ic50,
                    percentile_rank=pred.get("percentile_rank"),
                    confidence_tier=self._score_htl(pred),
                    tool_outputs={
                        **pred,
                        "ic50_note": "approximated_from_percentile_rank",
                        "method_used": "IEDB_NetMHCIIpan4.3",
                    },
                ))
            except Exception as e:
                logger.warning(f"HTL process error: {e}")
        epitopes.sort(key=lambda x: x.ic50_nm or 50000)
        return epitopes[:15]

    def _score_ctl(self, pred: Dict[str, Any]) -> ConfidenceTier:
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


def _infer_ctl_method(predictions: List[Dict]) -> str:
    if not predictions:
        return "none_CTL_unavailable"
    method = predictions[0].get("prediction_method", "")
    if "MHCflurry" in method:
        return "MHCflurry_2.0_affinity_fallback"
    return "IEDB_NetMHCpan4.1_EL"


def _infer_htl_method(predictions: List[Dict]) -> str:
    if not predictions:
        return "none_HTL_unavailable"
    return "IEDB_NetMHCIIpan4.3"


tcell_predictor = TCellPredictorAgent()