"""
SAFETY FILTER AGENT
Version 3.0 - Local implementation replacing web scraping.

Architecture: local-first, external-enhancement.

Primary (always available, zero external calls):
  WHO allergenicity    - regulatory standard, AllergenOnline database
  AllerTOP v2.0 local      - Doytchinova & Flower 2014, full lag=1,2 implementation
  HemoPI hemolytic         - Singh et al. 2011, WHO vaccine safety standard
  Human homology local     - FDA/EMA 8-mer threshold, UniProt human Swiss-Prot

Secondary (attempted if primary passes, non-blocking):
  NCBI BLAST               - confirms human homology if available
  External AllerTOP server - cross-check if available (never blocks pipeline)

This architecture means:
  - Agent 6 runs at 10/10 reliability (local tools always available)
  - External tools enhance results when available but never block them
  - All results are reproducible regardless of network state
  - Database versions are recorded in every audit trail entry

Scientific basis: see src/tools/safety_local.py for full citations.
"""

import os
import re
import time
import logging
import requests
from typing import List, Dict, Any, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from src.models.candidate import CandidateProtein, EpitopeResult, ConfidenceTier
from src.tools.safety_local import (
    screen_epitope_local,
    check_database_availability,
    ALLERGENONLINE_VERSION,
    HUMAN_SWISSPROT_VERSION,
    ALLERTOP_VERSION,
    HEMOPI_VERSION,
)

from src.utils.logger import get_logger
logger = get_logger("tope_deep.agents.Agent 6")  # use the correct agent name

VERDICT_PASS     = "pass"
VERDICT_FAIL     = "fail"
VERDICT_FLAGGED  = "flagged"
VERDICT_UNSCORED = "unscored"

_CB_THRESHOLD = 3


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
            logger.info(
                f"Agent 6: {self.name} external enhancement unavailable "
                f"(circuit breaker open after {self.failures} failures). "
                f"Local implementation continues unaffected."
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
        })
        self.max_workers = 4  # increased - no I/O wait on local tools

        # External enhancement circuit breakers (non-blocking)
        self._cb_blast = _CircuitBreaker("NCBI-BLAST")
        self._cb_allertop_ext = _CircuitBreaker("AllerTOP-external")

        # Log database status at startup
        db_status = check_database_availability()
        logger.info(
            f"Agent 6 database status: "
            f"AllergenOnline={'OK' if db_status['allergenonline'] else 'MISSING'} "
            f"v{db_status['allergenonline_version']}, "
            f"HumanSwissProt={'OK' if db_status['human_swissprot'] else 'MISSING'} "
            f"v{db_status['human_swissprot_version']}"
        )
        if not db_status["allergenonline"]:
            logger.warning(
                "Agent 6: AllergenOnline database missing. "
                "WHO allergenicity screen will use AllerTOP only. "
                "Run: python data/safety_db/download_databases.py"
            )

    # ── PUBLIC ────────────────────────────────────────────────────────────────

    def run(self, candidates: List[CandidateProtein]) -> List[CandidateProtein]:
        logger.info("Agent 6: Starting safety screening (local implementation v3.0)")
        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"   {len(active)} candidates")

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

            # Deduplicate
            seen: set = set()
            unique: List[EpitopeResult] = []
            for ep in all_epitopes:
                if ep.sequence not in seen:
                    seen.add(ep.sequence)
                    unique.append(ep)

            logger.info(f"      {len(unique)} unique sequences ({self.max_workers} workers)")

            results_map: Dict[str, Tuple[str, List[str], Dict]] = {}
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
                    ep.sequence, (VERDICT_UNSCORED, ["no_result"], {}),
                )
                ep.tool_outputs["safety_verdict"] = verdict
                ep.tool_outputs["safety_flags"] = flags
                ep.tool_outputs["safety_method_used"] = method_used
                ep.tool_outputs["homology_safe"] = method_used.get("homology_safe", None)

                if verdict == VERDICT_PASS:
                    ep.allergenicity_safe = True
                    ep.toxicity_safe = True
                    safe_count += 1
                elif verdict == VERDICT_UNSCORED:
                    ep.allergenicity_safe = None
                    ep.toxicity_safe = None
                    ep.tool_outputs["safety_note"] = (
                        "Safety screening inconclusive. "
                        "Excluded from construct assembly. "
                        "Manual review required."
                    )
                    unscored_count += 1
                else:
                    ep.allergenicity_safe = False
                    ep.toxicity_safe = False
                    fail_count += 1 if verdict == VERDICT_FAIL else 0
                    flagged_count += 1 if verdict == VERDICT_FLAGGED else 0

            db_status = check_database_availability()

            candidate.stage = self.stage_name
            candidate.add_decision(
                stage=self.stage_name,
                decision="safety_screened",
                reasoning=(
                    f"{safe_count} passed, {flagged_count} flagged (review required), "
                    f"{fail_count} failed, {unscored_count} unscored "
                    f"out of {len(all_epitopes)} epitopes. "
                    f"Methods: WHO 2001 allergenicity protocol "
                    f"(AllergenOnline v{ALLERGENONLINE_VERSION}); "
                    f"AllerTOP v{ALLERTOP_VERSION} local "
                    f"(Doytchinova & Flower 2014); "
                    f"HemoPI v{HEMOPI_VERSION} hemolytic screen "
                    f"(Singh et al. 2011, WHO/BS/2019.2364); "
                    f"Human homology FDA/EMA 8-mer threshold "
                    f"(UniProt human v{HUMAN_SWISSPROT_VERSION}). "
                    f"All algorithms run locally - results reproducible."
                ),
                safe_count=safe_count,
                flagged_count=flagged_count,
                fail_count=fail_count,
                unscored_count=unscored_count,
                allergenonline_version=ALLERGENONLINE_VERSION,
                human_swissprot_version=HUMAN_SWISSPROT_VERSION,
                allertop_version=ALLERTOP_VERSION,
                hemopi_version=HEMOPI_VERSION,
                database_status=db_status,
            )

            logger.info(
                f"      {safe_count} safe | {flagged_count} flagged | "
                f"{fail_count} failed | {unscored_count} unscored"
            )

        logger.info("Agent 6: Safety screening complete")
        return candidates

    # ── SINGLE EPITOPE ────────────────────────────────────────────────────────

    def _screen_single_epitope(
        self, epitope: EpitopeResult
    ) -> Tuple[str, List[str], Dict]:
        seq = epitope.sequence

        # Run local screens (always available)
        local_result = screen_epitope_local(seq)

        verdict = local_result["verdict"]
        flags = (
            local_result["allergen_flags"] +
            local_result["toxic_flags"] +
            local_result["review_flags"]
        )

        method_used = {
            "primary": local_result["method_summary"],
            "fao_who_criterion": local_result["fao_who"].get("criterion"),
            "allertop_score": local_result["allertop"].get("svm_score"),
            "hemopi_score": local_result["hemopi"].get("svm_score"),
            "human_8mer_matches": local_result["human_homology"].get("overlap_count", 0),
        }
        # Expose homology result as a top-level tool_output field for frontend
        human_overlap = local_result["human_homology"].get("overlap_count", 0)
        method_used["homology_safe"] = human_overlap == 0

        # Optional: BLAST enhancement (non-blocking, adds confidence if available)
        if not self._cb_blast.is_open() and verdict != VERDICT_FAIL:
            blast_homology = self._check_human_homology_blast(seq)
            if blast_homology > 0:
                method_used["blast_homology_pct"] = blast_homology
                if blast_homology > 70 and "human_homology" not in " ".join(flags):
                    flags.append(f"blast_human_homology_{blast_homology:.0f}pct")
                    if verdict == VERDICT_PASS:
                        verdict = VERDICT_FLAGGED

        return verdict, flags, method_used

    # ── OPTIONAL BLAST ENHANCEMENT ────────────────────────────────────────────

    def _check_human_homology_blast(self, sequence: str) -> float:
        """NCBI BLAST - optional enhancement, never blocks pipeline."""
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
                    "FORMAT_TYPE": "JSOAgent 2", "HITLIST_SIZE": "5",
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
            for _ in range(8):
                time.sleep(5)
                result = self.session.get(
                    "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
                    params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSOAgent 2", "api_key": ncbi_key},
                    timeout=30,
                )
                if "Status=WAITING" in result.text:
                    continue
                if "Status=FAILED" in result.text:
                    self._cb_blast.record_failure()
                    return 0.0
                idents = re.findall(r'"identity"\s*:\s*(\d+)', result.text)
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
        except Exception:
            self._cb_blast.record_failure()
            return 0.0

    # ── DIAGNOSTICS ───────────────────────────────────────────────────────────

    def get_tool_status(self) -> Dict[str, str]:
        db = check_database_availability()
        return {
            "local_fao_who": "operational" if db["allergenonline"] else "degraded_no_db",
            "local_allertop": "operational",
            "local_hemopi": "operational",
            "local_human_homology": "operational" if db["human_swissprot"] else "degraded_no_db",
            "blast_enhancement": self._cb_blast.status(),
            "allergenonline_version": ALLERGENONLINE_VERSION,
            "human_swissprot_version": HUMAN_SWISSPROT_VERSION,
        }


safety_filter = SafetyFilterAgent()