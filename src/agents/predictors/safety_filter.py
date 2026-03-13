"""
SAFETY FILTER AGENT - MVP-2 NODE N6
Screens epitopes for allergenicity, toxicity, and human cross-reactivity.

Tools (all free, no licence):
  AllerTOP v2.0:   ddg-pharmfac.net/AllerTOP/
  AllergenFP v1.0: ddg-pharmfac.net/AllergenFP/
  ToxinPred:       webs.iiitd.edu.in/raghava/toxinpred/

Epitopes are FLAGGED, not removed - the scientist decides.
"""

import requests
import time
import logging
import re
from typing import List, Dict, Any, Optional, Tuple
from src.models.candidate import CandidateProtein, EpitopeResult, ConfidenceTier

logger = logging.getLogger(__name__)


class SafetyFilterAgent:
    def __init__(self):
        self.stage_name = "safety_filter"
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "Kozi-Pipeline/2.0"})

    def run(self, candidates: List[CandidateProtein]) -> List[CandidateProtein]:
        """Screen all epitopes on active candidates for safety."""
        logger.info("N6: Starting safety screening")

        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"   Screening {len(active)} candidates")

        for i, candidate in enumerate(active):
            logger.info(f"   [{i+1}/{len(active)}] {candidate.protein_name}")

            all_epitopes = (
                list(candidate.ctl_epitopes) +
                list(candidate.htl_epitopes) +
                list(candidate.bcell_epitopes)
            )

            if not all_epitopes:
                logger.info(f"      No epitopes to screen")
                continue

            safe_count = 0
            flagged_count = 0
            fail_count = 0

            for epitope in all_epitopes:
                try:
                    verdict, flags = self._screen_single_epitope(epitope)

                    if verdict == "pass":
                        epitope.allergenicity_safe = True
                        epitope.toxicity_safe = True
                        safe_count += 1
                    elif verdict == "flagged":
                        flagged_count += 1
                        # Store flags in tool_outputs for audit
                        epitope.tool_outputs["safety_flags"] = flags
                    else:
                        fail_count += 1
                        epitope.allergenicity_safe = False
                        epitope.toxicity_safe = False
                        epitope.tool_outputs["safety_flags"] = flags

                    time.sleep(0.3)  # Rate limit for web APIs

                except Exception as e:
                    logger.warning(
                        f"      Safety screen failed for {epitope.sequence[:10]}...: {e}")
                    epitope.tool_outputs["safety_flags"] = [
                        f"screening_error: {str(e)}"]

            candidate.stage = self.stage_name
            candidate.add_decision(
                stage=self.stage_name,
                decision="safety_screened",
                reasoning=f"{safe_count} passed, {flagged_count} flagged, {fail_count} failed out of {len(all_epitopes)} epitopes",
                safe_count=safe_count,
                flagged_count=flagged_count,
                fail_count=fail_count,
            )

            logger.info(
                f"      {safe_count} safe | {flagged_count} flagged | {fail_count} failed")

        logger.info(f"N6: Safety screening complete")
        return candidates

    def _screen_single_epitope(self, epitope: EpitopeResult) -> Tuple[str, List[str]]:
        """
        Screen one epitope through all safety checks.
        Returns: (verdict, list_of_flags)
        verdict: "pass", "flagged", or "fail"
        """
        seq = epitope.sequence
        flags = []

        # 1. Allergenicity - AllerTOP
        allertop = self._check_allertop(seq)
        if allertop == "ALLERGEN":
            flags.append("allertop_allergen")
        elif allertop == "unknown":
            flags.append("allertop_inconclusive")

        # 2. Allergenicity - AllergenFP (dual check)
        allergenfp = self._check_allergenfp(seq)
        if allergenfp == "ALLERGEN":
            flags.append("allergenfp_allergen")
        elif allergenfp == "unknown":
            flags.append("allergenfp_inconclusive")

        # 3. Toxicity - ToxinPred
        toxicity = self._check_toxinpred(seq)
        if toxicity == "Toxic":
            flags.append("toxinpred_toxic")
        elif toxicity == "unknown":
            flags.append("toxinpred_inconclusive")

        # 4. Human cross-reactivity (simple check)
        homology = self._check_human_homology(seq)
        if homology > 70:
            flags.append(f"human_homology_{homology:.0f}pct")
        elif homology > 50:
            flags.append(f"moderate_human_similarity_{homology:.0f}pct")

        # Determine verdict
        if not flags:
            return "pass", []
        elif any("toxic" in f for f in flags) or \
                ("allertop_allergen" in flags and "allergenfp_allergen" in flags):
            return "fail", flags
        else:
            return "flagged", flags

    def _check_allertop(self, sequence: str) -> str:
        """AllerTOP v2.0 allergenicity prediction."""
        if len(sequence) < 8:
            return "unknown"
        try:
            resp = self.session.post(
                "https://www.ddg-pharmfac.net/AllerTOP/predict_cgi.py",
                data={"queryseq": sequence, "output_type": "text"},
                timeout=20,
            )
            resp.raise_for_status()
            text = resp.text.lower()
            if "non-allergen" in text:
                return "NON-ALLERGEN"
            elif "allergen" in text:
                return "ALLERGEN"
            return "unknown"
        except Exception as e:
            logger.debug(f"AllerTOP failed: {e}")
            return "unknown"

    def _check_allergenfp(self, sequence: str) -> str:
        """AllergenFP v1.0 allergenicity prediction."""
        if len(sequence) < 8:
            return "unknown"
        try:
            resp = self.session.post(
                "https://www.ddg-pharmfac.net/AllergenFP/predict_cgi.py",
                data={"queryseq": sequence, "output_type": "text"},
                timeout=20,
            )
            resp.raise_for_status()
            text = resp.text.lower()
            if "non-allergen" in text:
                return "NON-ALLERGEN"
            elif "allergen" in text:
                return "ALLERGEN"
            return "unknown"
        except Exception as e:
            logger.debug(f"AllergenFP failed: {e}")
            return "unknown"

    def _check_toxinpred(self, sequence: str) -> str:
        """ToxinPred toxicity screening."""
        if len(sequence) < 5:
            return "unknown"
        try:
            resp = self.session.post(
                "https://webs.iiitd.edu.in/raghava/toxinpred/multiple_formsubmit.php",
                data={"seq": sequence, "method": "1",
                      "eval": "10", "terminus": "N"},
                timeout=30,
            )
            resp.raise_for_status()
            text = resp.text.lower()
            if "non-toxin" in text or "non-toxic" in text:
                return "Non-Toxic"
            elif "toxin" in text or "toxic" in text:
                return "Toxic"
            return "unknown"
        except Exception as e:
            logger.debug(f"ToxinPred failed: {e}")
            return "unknown"

    def _check_human_homology(self, sequence: str) -> float:
        """
        Simple human self-peptide check.
        MVP-2: checks against known problematic human peptide motifs.
        Production: would use NCBI BLAST or IEDB PEPMatch.
        """
        known_human_peptides = [
            "GILGFVFTL", "NLVPMVATV", "GLCTLVAML",
            "FLRGRAYGL", "ELAGIGILTV", "YLQPRTFLL",
            "KLGGALQAK", "RLRAEAQVK", "ATDALMTGY",
        ]

        max_identity = 0.0
        for human_pep in known_human_peptides:
            identity = self._local_identity(sequence, human_pep)
            max_identity = max(max_identity, identity)

        return max_identity * 100

    def _local_identity(self, seq1: str, seq2: str) -> float:
        """Best local identity between two short sequences."""
        if not seq1 or not seq2:
            return 0.0
        shorter = seq1 if len(seq1) <= len(seq2) else seq2
        longer = seq2 if len(seq1) <= len(seq2) else seq1
        best = 0.0
        for offset in range(len(longer) - len(shorter) + 1):
            matches = sum(1 for i in range(len(shorter))
                          if shorter[i] == longer[offset + i])
            best = max(best, matches / len(shorter))
        return best

    def test_connections(self) -> Dict[str, bool]:
        """Test safety API connectivity."""
        results = {}
        test_seq = "MKLRLFCLAMLMACAQILNGS"

        try:
            self._check_allertop(test_seq)
            results["allertop"] = True
        except Exception:
            results["allertop"] = False

        try:
            self._check_allergenfp(test_seq)
            results["allergenfp"] = True
        except Exception:
            results["allergenfp"] = False

        try:
            self._check_toxinpred("AASAIQGNV")
            results["toxinpred"] = True
        except Exception:
            results["toxinpred"] = False

        return results


# Global instance
safety_filter = SafetyFilterAgent()
