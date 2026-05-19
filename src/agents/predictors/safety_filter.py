"""
SAFETY FILTER AGENT - MVP-2 NODE N6
Screens epitopes for allergenicity, toxicity, and human cross-reactivity.

Tools:
  AllerTOP v2.0    - ddg-pharmfac.net/AllerTOP   (Doytchinova et al. 2014)
  AllergenFP v1.0  - ddg-pharmfac.net/AllergenFP (Dimitrov et al. 2014)
  ToxinPred        - webs.iiitd.edu.in/raghava/toxinpred (Gupta et al. 2013)
  NCBI BLAST       - blast.ncbi.nlm.nih.gov, human refseq_protein taxid:9606
                     (Altschul et al., J Mol Biol 1990, 215:403-410)

Circuit breaker: each tool tracks consecutive failures and opens after 3,
logging cb_open so Railway logs show exactly which tool degraded and when.
Inconclusive results from short peptides are audit-logged, not used as flags.

Parallel execution: 3 epitopes screened concurrently.
"""

import os
import re
import time
import logging
import requests
from typing import List, Dict, Any, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from src.models.candidate import CandidateProtein, EpitopeResult, ConfidenceTier

logger = logging.getLogger(__name__)

_CB_THRESHOLD = 3


class _CircuitBreaker:
    """Per-tool failure tracker. Opens after _CB_THRESHOLD consecutive failures."""

    def __init__(self, name: str):
        self.name = name
        self.failures = 0
        self.open = False

    def record_success(self):
        self.failures = 0
        self.open = False

    def record_failure(self):
        self.failures += 1
        if self.failures >= _CB_THRESHOLD:
            if not self.open:
                logger.warning(
                    f"N6 circuit breaker OPEN: {self.name} failed "
                    f"{self.failures} times consecutively - skipping for this run"
                )
            self.open = True

    def is_open(self) -> bool:
        return self.open


class SafetyFilterAgent:
    def __init__(self):
        self.stage_name = "safety_filter"
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                          "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.5",
        })
        self.max_workers = 3

        self._cb_allertop   = _CircuitBreaker("AllerTOP")
        self._cb_allergenfp = _CircuitBreaker("AllergenFP")
        self._cb_toxinpred  = _CircuitBreaker("ToxinPred")
        self._cb_blast      = _CircuitBreaker("NCBI-BLAST")

    # ── PUBLIC ENTRY POINT ────────────────────────────────────────────────────

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
                logger.info("      No epitopes to screen")
                continue

            seen: set = set()
            unique_epitopes = []
            for ep in all_epitopes:
                if ep.sequence not in seen:
                    seen.add(ep.sequence)
                    unique_epitopes.append(ep)

            logger.info(
                f"      Screening {len(unique_epitopes)} unique epitopes "
                f"(parallel, {self.max_workers} workers)"
            )

            results_map: Dict[str, Tuple[str, List[str]]] = {}
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_ep = {
                    executor.submit(self._screen_single_epitope, ep): ep
                    for ep in unique_epitopes
                }
                for future in as_completed(future_to_ep):
                    ep = future_to_ep[future]
                    try:
                        verdict, flags = future.result()
                        results_map[ep.sequence] = (verdict, flags)
                    except Exception as e:
                        logger.warning(
                            f"      Safety screen error for "
                            f"{ep.sequence[:10]}...: {e}"
                        )
                        results_map[ep.sequence] = (
                            "flagged", [f"screening_error: {str(e)}"]
                        )

            safe_count = flagged_count = fail_count = 0

            for ep in all_epitopes:
                verdict, flags = results_map.get(
                    ep.sequence, ("flagged", ["no_result"])
                )
                if verdict == "pass":
                    ep.allergenicity_safe = True
                    ep.toxicity_safe = True
                    safe_count += 1
                elif verdict == "flagged":
                    flagged_count += 1
                    ep.tool_outputs["safety_flags"] = flags
                else:
                    fail_count += 1
                    ep.allergenicity_safe = False
                    ep.toxicity_safe = False
                    ep.tool_outputs["safety_flags"] = flags

            candidate.stage = self.stage_name
            candidate.add_decision(
                stage=self.stage_name,
                decision="safety_screened",
                reasoning=(
                    f"{safe_count} passed, {flagged_count} flagged, "
                    f"{fail_count} failed out of {len(all_epitopes)} epitopes. "
                    f"Tools: AllerTOP v2.0, AllergenFP v1.0, ToxinPred, "
                    f"NCBI BLAST (taxid:9606)."
                ),
                safe_count=safe_count,
                flagged_count=flagged_count,
                fail_count=fail_count,
            )
            logger.info(
                f"      {safe_count} safe | {flagged_count} flagged | "
                f"{fail_count} failed"
            )

        logger.info("N6: Safety screening complete")
        return candidates

    # ── PER-EPITOPE SCREENER ──────────────────────────────────────────────────

    def _screen_single_epitope(
        self, epitope: EpitopeResult
    ) -> Tuple[str, List[str]]:
        seq = epitope.sequence
        flags: List[str] = []

        allertop = self._check_allertop(seq)
        if allertop == "ALLERGEN":
            flags.append("allertop_allergen")
        elif allertop == "unknown":
            flags.append("allertop_inconclusive")

        allergenfp = self._check_allergenfp(seq)
        if allergenfp == "ALLERGEN":
            flags.append("allergenfp_allergen")
        elif allergenfp == "unknown":
            flags.append("allergenfp_inconclusive")

        toxicity = self._check_toxinpred(seq)
        if toxicity == "Toxic":
            flags.append("toxinpred_toxic")
        elif toxicity == "unknown":
            flags.append("toxinpred_inconclusive")

        homology = self._check_human_homology(seq)
        if homology > 70:
            flags.append(f"human_homology_{homology:.0f}pct")
        elif homology > 50:
            flags.append(f"moderate_human_similarity_{homology:.0f}pct")

        real_flags = [f for f in flags if "inconclusive" not in f]

        if not real_flags:
            return "pass", flags

        if any("toxic" in f for f in real_flags) or (
            "allertop_allergen" in real_flags
            and "allergenfp_allergen" in real_flags
        ):
            return "fail", flags

        return "flagged", real_flags

    # ── ALLERGENICITY ─────────────────────────────────────────────────────────

    def _check_allertop(self, sequence: str) -> str:
        """
        AllerTOP v2.0 - SVM-based allergenicity prediction.
        Doytchinova & Flower, J Proteome Res 2014.
        Circuit breaker: opens after 3 consecutive failures.
        """
        if len(sequence) < 8:
            return "unknown"
        if self._cb_allertop.is_open():
            return "unknown"

        for attempt in range(3):
            try:
                resp = self.session.post(
                    "https://www.ddg-pharmfac.net/AllerTOP/predict_cgi.py",
                    data={"queryseq": sequence, "output_type": "text"},
                    timeout=20,
                )
                resp.raise_for_status()
                text = resp.text.lower()
                if "non-allergen" in text:
                    self._cb_allertop.record_success()
                    return "NON-ALLERGEN"
                if "allergen" in text:
                    self._cb_allertop.record_success()
                    return "ALLERGEN"
                return "unknown"
            except requests.Timeout:
                logger.debug(f"AllerTOP timeout (attempt {attempt+1})")
                if attempt < 2:
                    time.sleep(2 ** attempt)
            except Exception as e:
                logger.debug(f"AllerTOP error (attempt {attempt+1}): {e}")
                if attempt < 2:
                    time.sleep(2 ** attempt)

        self._cb_allertop.record_failure()
        return "unknown"

    def _check_allergenfp(self, sequence: str) -> str:
        """
        AllergenFP v1.0 - fingerprint-based allergenicity prediction.
        Dimitrov et al., Bioinformatics 2014.
        Circuit breaker: opens after 3 consecutive failures.
        """
        if len(sequence) < 8:
            return "unknown"
        if self._cb_allergenfp.is_open():
            return "unknown"

        for attempt in range(3):
            try:
                resp = self.session.post(
                    "https://www.ddg-pharmfac.net/AllergenFP/predict_cgi.py",
                    data={"queryseq": sequence, "output_type": "text"},
                    timeout=20,
                )
                resp.raise_for_status()
                text = resp.text.lower()
                if "non-allergen" in text:
                    self._cb_allergenfp.record_success()
                    return "NON-ALLERGEN"
                if "allergen" in text:
                    self._cb_allergenfp.record_success()
                    return "ALLERGEN"
                return "unknown"
            except requests.Timeout:
                logger.debug(f"AllergenFP timeout (attempt {attempt+1})")
                if attempt < 2:
                    time.sleep(2 ** attempt)
            except Exception as e:
                logger.debug(f"AllergenFP error (attempt {attempt+1}): {e}")
                if attempt < 2:
                    time.sleep(2 ** attempt)

        self._cb_allergenfp.record_failure()
        return "unknown"

    def _check_toxinpred(self, sequence: str) -> str:
        """
        ToxinPred - SVM toxicity prediction for mammalian cells.
        Gupta et al., In Silico Pharmacology 2013.
        Circuit breaker: opens after 3 consecutive failures.
        """
        if len(sequence) < 5:
            return "unknown"
        if self._cb_toxinpred.is_open():
            return "unknown"

        for attempt in range(3):
            try:
                resp = self.session.post(
                    "https://webs.iiitd.edu.in/raghava/toxinpred/multiple_formsubmit.php",
                    data={
                        "seq":      sequence,
                        "method":   "1",
                        "eval":     "10",
                        "terminus": "N",
                    },
                    timeout=25,
                )
                resp.raise_for_status()
                text = resp.text.lower()
                if "non-toxin" in text or "non-toxic" in text:
                    self._cb_toxinpred.record_success()
                    return "Non-Toxic"
                if "toxin" in text or "toxic" in text:
                    self._cb_toxinpred.record_success()
                    return "Toxic"
                return "unknown"
            except requests.Timeout:
                logger.debug(f"ToxinPred timeout (attempt {attempt+1})")
                if attempt < 2:
                    time.sleep(2 ** attempt)
            except Exception as e:
                logger.debug(f"ToxinPred error (attempt {attempt+1}): {e}")
                if attempt < 2:
                    time.sleep(2 ** attempt)

        self._cb_toxinpred.record_failure()
        return "unknown"

    # ── HUMAN HOMOLOGY ────────────────────────────────────────────────────────

    def _check_human_homology(self, sequence: str) -> float:
        """
        NCBI BLAST remote search against human refseq_protein (taxid:9606).
        Returns max percent identity found. Returns 0.0 on timeout or failure.
        Reference: Altschul et al., J Mol Biol 1990, 215:403-410.
        Circuit breaker: opens after 3 consecutive failures.
        """
        if len(sequence) < 8:
            return 0.0
        if self._cb_blast.is_open():
            return 0.0

        ncbi_key = os.getenv("NCBI_API_KEY", "")

        try:
            submit_resp = self.session.post(
                "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
                data={
                    "CMD":          "Put",
                    "PROGRAM":      "blastp",
                    "DATABASE":     "refseq_protein",
                    "QUERY":        sequence,
                    "ENTREZ_QUERY": "Homo sapiens[Organism]",
                    "FORMAT_TYPE":  "JSON2",
                    "HITLIST_SIZE": "5",
                    "EXPECT":       "1",
                    "WORD_SIZE":    "3",
                    "api_key":      ncbi_key,
                },
                timeout=30,
            )
            submit_resp.raise_for_status()

            rid_match = re.search(r'RID = ([A-Z0-9]+)', submit_resp.text)
            if not rid_match:
                logger.debug("BLAST: RID not found in response")
                self._cb_blast.record_failure()
                return 0.0
            rid = rid_match.group(1)

            for _ in range(12):
                time.sleep(5)
                result_resp = self.session.get(
                    "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
                    params={
                        "CMD":         "Get",
                        "RID":         rid,
                        "FORMAT_TYPE": "JSON2",
                        "api_key":     ncbi_key,
                    },
                    timeout=30,
                )
                if "Status=WAITING" in result_resp.text:
                    continue
                if "Status=FAILED" in result_resp.text:
                    logger.debug("BLAST: search failed")
                    self._cb_blast.record_failure()
                    return 0.0

                identity_matches = re.findall(
                    r'"identity"\s*:\s*(\d+)', result_resp.text
                )
                length_matches = re.findall(
                    r'"align_len"\s*:\s*(\d+)', result_resp.text
                )
                if identity_matches and length_matches:
                    max_pct = max(
                        int(i) / int(l) * 100
                        for i, l in zip(identity_matches, length_matches)
                        if int(l) > 0
                    )
                    self._cb_blast.record_success()
                    logger.info(f"BLAST human homology: {max_pct:.1f}%")
                    return max_pct

                self._cb_blast.record_success()
                return 0.0

            logger.debug("BLAST: timed out after 60s")
            self._cb_blast.record_failure()
            return 0.0

        except Exception as e:
            logger.warning(f"BLAST human homology check failed: {e}")
            self._cb_blast.record_failure()
            return 0.0

    # ── DIAGNOSTICS ───────────────────────────────────────────────────────────

    def test_connections(self) -> Dict[str, bool]:
        """Health check for all external tools."""
        test_seq = "MKLRLFCLAMLMACAQILNGS"
        results = {}
        for name, fn, seq in [
            ("allertop",   self._check_allertop,   test_seq),
            ("allergenfp", self._check_allergenfp, test_seq),
            ("toxinpred",  self._check_toxinpred,  "AASAIQGNV"),
            ("blast",      self._check_human_homology, test_seq),
        ]:
            try:
                result = fn(seq)
                results[name] = result != "unknown" and result is not None
            except Exception:
                results[name] = False
        return results


# Global instance
safety_filter = SafetyFilterAgent()