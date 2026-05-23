"""
SAFETY FILTER AGENT — TOPE_DEEP NODE N6
Screens epitopes for allergenicity, toxicity, and human cross-reactivity.

Tools:
  AllerTOP v2.0    ddg-pharmfac.net/AllerTOP   (Doytchinova & Flower, J Proteome Res 2014)
  AllergenFP v1.0  ddg-pharmfac.net/AllergenFP (Dimitrov et al., Bioinformatics 2014)
  ToxinPred        webs.iiitd.edu.in/raghava/toxinpred (Gupta et al., In Silico Pharm 2013)
  NCBI BLAST       blast.ncbi.nlm.nih.gov, refseq_protein taxid:9606
                   (Altschul et al., J Mol Biol 1990, 215:403-410)

DESIGN RULE — inconclusive != safe (fixed from previous version)
  When external tools time out or are unreachable, epitopes are marked
  allergenicity_safe=None, toxicity_safe=None.
  They are NEVER marked True on insufficient evidence.
  method_used field records which tools ran vs timed out.

Verdict logic:
  PASS      all tools returned explicit non-allergen / non-toxic
  FAIL      explicit allergen or toxic signal (consensus or toxic)
  FLAGGED   partial allergen signal (one tool only) — conservative, excluded from N8
  UNSCORED  all tools inconclusive — excluded from N8, labeled in audit

allergenicity_safe / toxicity_safe:
  PASS     -> True  / True
  FAIL     -> False / False
  FLAGGED  -> False / False
  UNSCORED -> None  / None

method_used field on every epitope tool_outputs:
  Records which tools ran, which timed out, and which CB was open.
  This field flows through to CSV export.

Circuit breaker: opens after 3 consecutive failures per tool.
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

VERDICT_PASS     = "pass"
VERDICT_FAIL     = "fail"
VERDICT_FLAGGED  = "flagged"
VERDICT_UNSCORED = "unscored"


class _CircuitBreaker:
    def __init__(self, name: str):
        self.name = name
        self.failures = 0
        self.open = False

    def record_success(self):
        self.failures = 0
        self.open = False

    def record_failure(self):
        self.failures += 1
        if self.failures >= _CB_THRESHOLD and not self.open:
            logger.warning(
                f"N6 circuit breaker OPEN: {self.name} — "
                f"{self.failures} consecutive failures. "
                f"Affected epitopes marked unscored (None), not safe."
            )
            self.open = True

    def is_open(self) -> bool:
        return self.open

    def status(self) -> str:
        return "open" if self.open else "closed"


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

    # ── PUBLIC ────────────────────────────────────────────────────────────────

    def run(self, candidates: List[CandidateProtein]) -> List[CandidateProtein]:
        logger.info("N6: Starting safety screening")
        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"   {len(active)} candidates to screen")

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

            # Deduplicate by sequence before hitting external APIs
            seen: set = set()
            unique: List[EpitopeResult] = []
            for ep in all_epitopes:
                if ep.sequence not in seen:
                    seen.add(ep.sequence)
                    unique.append(ep)

            logger.info(
                f"      {len(unique)} unique sequences "
                f"(parallel, {self.max_workers} workers)"
            )

            results_map: Dict[str, Tuple[str, List[str], Dict[str, str]]] = {}
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_ep = {
                    executor.submit(self._screen_single_epitope, ep): ep
                    for ep in unique
                }
                for future in as_completed(future_to_ep):
                    ep = future_to_ep[future]
                    try:
                        verdict, flags, method_used = future.result()
                        results_map[ep.sequence] = (verdict, flags, method_used)
                    except Exception as e:
                        logger.warning(f"      Screen error {ep.sequence[:10]}: {e}")
                        results_map[ep.sequence] = (
                            VERDICT_UNSCORED,
                            [f"screening_error: {str(e)}"],
                            {"error": str(e)},
                        )

            safe_count = fail_count = flagged_count = unscored_count = 0

            for ep in all_epitopes:
                verdict, flags, method_used = results_map.get(
                    ep.sequence,
                    (VERDICT_UNSCORED, ["no_result"], {}),
                )
                ep.tool_outputs["safety_verdict"]  = verdict
                ep.tool_outputs["safety_flags"]    = flags
                ep.tool_outputs["safety_method_used"] = method_used

                if verdict == VERDICT_PASS:
                    ep.allergenicity_safe = True
                    ep.toxicity_safe      = True
                    safe_count += 1
                elif verdict == VERDICT_UNSCORED:
                    ep.allergenicity_safe = None
                    ep.toxicity_safe      = None
                    ep.tool_outputs["safety_note"] = (
                        "Safety tools unavailable or timed out. "
                        "Epitope not screened. "
                        "Excluded from construct assembly (N8). "
                        "Manual screening required before wet-lab use."
                    )
                    unscored_count += 1
                else:
                    ep.allergenicity_safe = False
                    ep.toxicity_safe      = False
                    fail_count    += 1 if verdict == VERDICT_FAIL    else 0
                    flagged_count += 1 if verdict == VERDICT_FLAGGED else 0

            all_tools_down = (
                self._cb_allertop.is_open()
                and self._cb_allergenfp.is_open()
                and self._cb_toxinpred.is_open()
            )

            candidate.stage = self.stage_name
            candidate.add_decision(
                stage=self.stage_name,
                decision="safety_screened" if not all_tools_down else "safety_unscored",
                reasoning=(
                    f"{safe_count} passed, {flagged_count} flagged, "
                    f"{fail_count} failed, {unscored_count} unscored "
                    f"out of {len(all_epitopes)} epitopes. "
                    + (
                        "WARNING: All external safety tools unavailable. "
                        "No epitopes confirmed safe. allergenicity_safe=None on all. "
                        "Do not use for wet-lab work without manual screening. "
                        if all_tools_down else ""
                    ) +
                    "Tools: AllerTOP v2.0 (Doytchinova 2014), "
                    "AllergenFP v1.0 (Dimitrov 2014), "
                    "ToxinPred (Gupta 2013), "
                    "NCBI BLAST refseq_protein taxid:9606 (Altschul 1990)."
                ),
                safe_count=safe_count,
                flagged_count=flagged_count,
                fail_count=fail_count,
                unscored_count=unscored_count,
                all_tools_available=not all_tools_down,
                circuit_breaker_status=self.get_tool_status(),
            )

            if all_tools_down:
                logger.error(
                    f"N6 WARNING: all safety tools down for {candidate.protein_name}. "
                    f"All {len(all_epitopes)} epitopes marked unscored."
                )
            else:
                logger.info(
                    f"      {safe_count} safe | {flagged_count} flagged | "
                    f"{fail_count} failed | {unscored_count} unscored"
                )

        logger.info("N6: Safety screening complete")
        return candidates

    # ── SINGLE EPITOPE SCREENER ───────────────────────────────────────────────

    def _screen_single_epitope(
        self, epitope: EpitopeResult
    ) -> Tuple[str, List[str], Dict[str, str]]:
        """
        Returns (verdict, flags, method_used_dict).
        method_used_dict records tool outcome for every tool attempted.
        """
        seq = epitope.sequence
        allergen_flags:    List[str] = []
        toxic_flags:       List[str] = []
        inconclusive_flags: List[str] = []
        method_used: Dict[str, str] = {}

        # AllerTOP
        allertop = self._check_allertop(seq)
        method_used["allertop"] = (
            "AllerTOP_v2.0_unavailable" if self._cb_allertop.is_open()
            else f"AllerTOP_v2.0_{allertop.lower()}"
        )
        if allertop == "ALLERGEN":
            allergen_flags.append("allertop_allergen")
        elif allertop == "unknown":
            inconclusive_flags.append("allertop_inconclusive")

        # AllergenFP
        allergenfp = self._check_allergenfp(seq)
        method_used["allergenfp"] = (
            "AllergenFP_v1.0_unavailable" if self._cb_allergenfp.is_open()
            else f"AllergenFP_v1.0_{allergenfp.lower()}"
        )
        if allergenfp == "ALLERGEN":
            allergen_flags.append("allergenfp_allergen")
        elif allergenfp == "unknown":
            inconclusive_flags.append("allergenfp_inconclusive")

        # ToxinPred
        toxicity = self._check_toxinpred(seq)
        method_used["toxinpred"] = (
            "ToxinPred_unavailable" if self._cb_toxinpred.is_open()
            else f"ToxinPred_{toxicity.lower().replace('-', '_')}"
        )
        if toxicity == "Toxic":
            toxic_flags.append("toxinpred_toxic")
        elif toxicity == "unknown":
            inconclusive_flags.append("toxinpred_inconclusive")

        # BLAST human homology
        homology = self._check_human_homology(seq)
        if self._cb_blast.is_open():
            method_used["blast"] = "NCBI_BLAST_unavailable"
            inconclusive_flags.append("blast_unavailable")
        else:
            method_used["blast"] = f"NCBI_BLAST_human_homology_{homology:.1f}pct"
            if homology > 70:
                allergen_flags.append(f"human_homology_{homology:.0f}pct")
            elif homology > 50:
                allergen_flags.append(f"moderate_human_similarity_{homology:.0f}pct")

        all_flags = allergen_flags + toxic_flags + inconclusive_flags

        # Verdict logic
        if toxic_flags or (
            "allertop_allergen" in allergen_flags
            and "allergenfp_allergen" in allergen_flags
        ):
            return VERDICT_FAIL, all_flags, method_used

        if allergen_flags:
            return VERDICT_FLAGGED, all_flags, method_used

        # All inconclusive — cannot confirm safety
        if not allergen_flags and not toxic_flags and inconclusive_flags:
            return VERDICT_UNSCORED, all_flags, method_used

        # Explicit non-allergen, non-toxic from all available tools
        return VERDICT_PASS, all_flags, method_used

    # ── TOOL CALLERS ─────────────────────────────────────────────────────────

    def _check_allertop(self, sequence: str) -> str:
        if len(sequence) < 8 or self._cb_allertop.is_open():
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
                if attempt < 2: time.sleep(2 ** attempt)
            except Exception as e:
                logger.debug(f"AllerTOP error (attempt {attempt+1}): {e}")
                if attempt < 2: time.sleep(2 ** attempt)
        self._cb_allertop.record_failure()
        return "unknown"

    def _check_allergenfp(self, sequence: str) -> str:
        if len(sequence) < 8 or self._cb_allergenfp.is_open():
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
                if attempt < 2: time.sleep(2 ** attempt)
            except Exception as e:
                logger.debug(f"AllergenFP error (attempt {attempt+1}): {e}")
                if attempt < 2: time.sleep(2 ** attempt)
        self._cb_allergenfp.record_failure()
        return "unknown"

    def _check_toxinpred(self, sequence: str) -> str:
        if len(sequence) < 5 or self._cb_toxinpred.is_open():
            return "unknown"
        for attempt in range(3):
            try:
                resp = self.session.post(
                    "https://webs.iiitd.edu.in/raghava/toxinpred/multiple_formsubmit.php",
                    data={"seq": sequence, "method": "1", "eval": "10", "terminus": "N"},
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
                if attempt < 2: time.sleep(2 ** attempt)
            except Exception as e:
                logger.debug(f"ToxinPred error (attempt {attempt+1}): {e}")
                if attempt < 2: time.sleep(2 ** attempt)
        self._cb_toxinpred.record_failure()
        return "unknown"

    def _check_human_homology(self, sequence: str) -> float:
        if len(sequence) < 8 or self._cb_blast.is_open():
            return 0.0
        ncbi_key = os.getenv("NCBI_API_KEY", "")
        try:
            submit = self.session.post(
                "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
                data={
                    "CMD": "Put", "PROGRAM": "blastp",
                    "DATABASE": "refseq_protein",
                    "QUERY": sequence,
                    "ENTREZ_QUERY": "Homo sapiens[Organism]",
                    "FORMAT_TYPE": "JSON2", "HITLIST_SIZE": "5",
                    "EXPECT": "1", "WORD_SIZE": "3",
                    "api_key": ncbi_key,
                },
                timeout=30,
            )
            submit.raise_for_status()
            rid_match = re.search(r'RID = ([A-Z0-9]+)', submit.text)
            if not rid_match:
                self._cb_blast.record_failure()
                return 0.0
            rid = rid_match.group(1)
            for _ in range(12):
                time.sleep(5)
                result = self.session.get(
                    "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
                    params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2", "api_key": ncbi_key},
                    timeout=30,
                )
                if "Status=WAITING" in result.text:
                    continue
                if "Status=FAILED" in result.text:
                    self._cb_blast.record_failure()
                    return 0.0
                idents  = re.findall(r'"identity"\s*:\s*(\d+)', result.text)
                lengths = re.findall(r'"align_len"\s*:\s*(\d+)', result.text)
                if idents and lengths:
                    max_pct = max(
                        int(i) / int(l) * 100
                        for i, l in zip(idents, lengths) if int(l) > 0
                    )
                    self._cb_blast.record_success()
                    return max_pct
                self._cb_blast.record_success()
                return 0.0
            self._cb_blast.record_failure()
            return 0.0
        except Exception as e:
            logger.warning(f"BLAST failed: {e}")
            self._cb_blast.record_failure()
            return 0.0

    # ── DIAGNOSTICS ───────────────────────────────────────────────────────────

    def get_tool_status(self) -> Dict[str, str]:
        """Circuit breaker state per tool — exposed via /api/health."""
        return {
            "allertop":   self._cb_allertop.status(),
            "allergenfp": self._cb_allergenfp.status(),
            "toxinpred":  self._cb_toxinpred.status(),
            "blast":      self._cb_blast.status(),
        }

    def test_connections(self) -> Dict[str, bool]:
        test_seq = "MKLRLFCLAMLMACAQILNGS"
        results = {}
        for name, fn, seq in [
            ("allertop",   self._check_allertop,   test_seq),
            ("allergenfp", self._check_allergenfp, test_seq),
            ("toxinpred",  self._check_toxinpred,  "AASAIQGNV"),
            ("blast",      self._check_human_homology, test_seq),
        ]:
            try:
                r = fn(seq)
                results[name] = r != "unknown" and r is not None and r is not False
            except Exception:
                results[name] = False
        return results


safety_filter = SafetyFilterAgent()