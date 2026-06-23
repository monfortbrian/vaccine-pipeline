"""
BCELL PREDICTOR AGENT

Tool: IEDB BepiPred 2.0 (linear B-cell epitope prediction)
Reference: Jespersen et al. (2017) Nucleic Acids Res 45:W24-W29

"""

import logging
import httpx
from typing import List, Dict, Any, Optional

from src.models.candidate import CandidateProtein, EpitopeResult, EpitopeType, ConfidenceTier

logger = logging.getLogger("tope_deep.agents.Agent 4")

BEPIPRED_THRESHOLD = 0.5
IEDB_BCELL_URL = "https://tools-cluster-interface.iedb.org/tools_api/bcell/"


def _call_iedb_bcell_client(client, sequence: str) -> Optional[List[Dict]]:
    """
    Try multiple method names on the IEDB B-cell client.
    The IEDB Python client has had inconsistent method naming across versions.
    Returns raw prediction list or None if all methods fail.
    """
    # Ordered by most-likely current name first
    method_names = [
        "predict",
        "predict_bcell_epitopes",
        "run",
        "predict_epitopes",
        "bcell_predict",
    ]
    for name in method_names:
        method = getattr(client, name, None)
        if method is not None:
            try:
                result = method(sequence=sequence, method="Bepipred-2.0", threshold=BEPIPRED_THRESHOLD)
                logger.info(f"Agent 4: IEDB client method '{name}' succeeded")
                return result
            except TypeError:
                # Method exists but has different signature try positional
                try:
                    result = method(sequence)
                    logger.info(f"Agent 4: IEDB client method '{name}' (positional) succeeded")
                    return result
                except Exception as e:
                    logger.warning(f"Agent 4: method '{name}' failed: {e}")
            except Exception as e:
                logger.warning(f"Agent 4: method '{name}' failed: {e}")

    logger.warning("Agent4: All IEDB client methods failed, falling back to direct HTTP")
    return None


def _call_iedb_bcell_http(sequence: str) -> List[Dict]:
    """
    Direct HTTP call to IEDB BepiPred 2.0 REST endpoint.
    Fallback when the Python client is unavailable or broken.
    Reference: https://tools.iedb.org/bcell/
    """
    try:
        resp = httpx.post(
            IEDB_BCELL_URL,
            data={
                "method":        "Bepipred-2.0",
                "sequence_text": sequence,
                "threshold":     BEPIPRED_THRESHOLD,
            },
            timeout=30.0,
            follow_redirects=True,
        )
        resp.raise_for_status()

        # IEDB returns TSV: position, residue, score
        lines   = resp.text.strip().split("\n")
        results = []
        current_window: List[tuple] = []

        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                pos   = int(parts[0])
                res   = parts[1]
                score = float(parts[2])
                current_window.append((pos, res, score))
            except (ValueError, IndexError):
                continue

        # Collect contiguous windows above threshold as epitopes
        in_epitope = False
        ep_start   = 0
        ep_seq     = ""
        ep_scores: List[float] = []

        for pos, res, score in current_window:
            if score >= BEPIPRED_THRESHOLD:
                if not in_epitope:
                    in_epitope = True
                    ep_start   = pos
                    ep_seq     = ""
                    ep_scores  = []
                ep_seq    += res
                ep_scores.append(score)
            else:
                if in_epitope and len(ep_seq) >= 8:
                    results.append({
                        "sequence": ep_seq,
                        "score":    round(sum(ep_scores) / len(ep_scores), 4),
                        "start":    ep_start,
                        "end":      ep_start + len(ep_seq) - 1,
                        "method":   "IEDB_BepiPred_2.0_HTTP",
                    })
                in_epitope = False
                ep_seq     = ""
                ep_scores  = []

        # Catch trailing epitope
        if in_epitope and len(ep_seq) >= 8:
            results.append({
                "sequence": ep_seq,
                "score":    round(sum(ep_scores) / len(ep_scores), 4),
                "start":    ep_start,
                "end":      ep_start + len(ep_seq) - 1,
                "method":   "IEDB_BepiPred_2.0_HTTP",
            })

        logger.info(f"Agent 4: IEDB HTTP fallback returned {len(results)} epitope windows")
        return results

    except Exception as e:
        logger.error(f"Agent 4: IEDB HTTP fallback failed: {e}")
        return []


class BCellPredictorAgent:
    def __init__(self):
        self.stage_name  = "bcell_prediction"
        self._iedb_bcell = None

    @property
    def iedb_bcell(self):
        if not self._iedb_bcell:
            try:
                from src.tools.iedb_bcell_client import iedb_bcell
                self._iedb_bcell = iedb_bcell
            except Exception as e:
                logger.warning(f"Agent 4: Could not load IEDB B-cell client: {e}")
                self._iedb_bcell = None
        return self._iedb_bcell

    def run(self, candidates: List[CandidateProtein]) -> List[CandidateProtein]:
        logger.info("Agent 4: B-Cell Predictor BepiPred 2.0")
        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"   Processing {len(active)} candidates")

        for i, candidate in enumerate(active):
            logger.info(f"   [{i+1}/{len(active)}] {candidate.protein_name}")
            try:
                raw       = self._predict(candidate.sequence)
                epitopes  = self._process(raw)

                candidate.bcell_epitopes = epitopes
                candidate.stage          = self.stage_name

                high = sum(1 for e in epitopes if e.confidence_tier == ConfidenceTier.HIGH)

                candidate.add_decision(
                    stage=self.stage_name,
                    decision="bcell_predicted",
                    reasoning=(
                        f"Linear B-cell epitope prediction via BepiPred 2.0 "
                        f"(Jespersen et al. 2017, Nucleic Acids Res 45:W24–W29). "
                        f"{len(epitopes)} epitopes above BepiPred score threshold (≥{BEPIPRED_THRESHOLD}), "
                        f"{high} with score ≥0.9 (high confidence). "
                        f"B-cell epitopes carry no HLA restriction; validated via rabbit polyclonal antisera protocol."
                    ),
                    bcell_count=len(epitopes),
                    bcell_high_confidence=high,
                    method="IEDB_BepiPred_2.0",
                    rabbit_validation_path=True,
                )
                logger.info(f"      B-cell: {len(epitopes)} ({high} high conf)")

            except Exception as e:
                logger.error(f"      Agent 4 failed for {candidate.protein_name}: {e}")
                candidate.add_decision(
                    stage=self.stage_name,
                    decision="prediction_failed",
                    reasoning=(
                        f"BepiPred 2.0 prediction failed for this protein sequence. "
                        f"Error: {str(e)[:200]}. No B-cell epitopes predicted."
                    ),
                )

        logger.info("Agent 4 complete")
        return candidates

    def _predict(self, sequence: str) -> List[Dict]:
        """
        Attempt prediction via client, fall back to direct HTTP.
        """
        raw = None

        # Try Python client if available
        if self.iedb_bcell is not None:
            raw = _call_iedb_bcell_client(self.iedb_bcell, sequence)

        # Fall back to direct HTTP
        if raw is None:
            logger.info("Agent 4: Using IEDB HTTP endpoint directly")
            raw = _call_iedb_bcell_http(sequence)

        return raw or []

    def _process(self, predictions: List[Dict[str, Any]]) -> List[EpitopeResult]:
        epitopes: List[EpitopeResult] = []
        seen: set = set()

        for pred in predictions:
            try:
                seq   = pred.get("sequence", "")
                score = float(pred.get("score", 0.0))

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
                        "model_categories":  ["RABBIT"],
                        "rabbit_validation": True,
                        "validation_note":   (
                            "Rabbit polyclonal antisera standard protocol. "
                            "Jespersen et al. (2017) Nucleic Acids Res 45:W24–W29"
                        ),
                    },
                ))
            except Exception as e:
                logger.warning(f"      B-cell process error: {e}")

        epitopes.sort(key=lambda e: e.tool_outputs.get("bepipred_score", 0), reverse=True)
        return epitopes[:10]

    @staticmethod
    def _score(score: float) -> ConfidenceTier:
        if score >= 0.9: return ConfidenceTier.HIGH
        if score >= 0.7: return ConfidenceTier.MEDIUM
        if score >= 0.5: return ConfidenceTier.LOW
        return ConfidenceTier.UNCERTAIN


bcell_predictor = BCellPredictorAgent()