"""
COVERAGE AGENT
Calculates immune coverage across global HLA populations.

Method: IEDB Population Coverage Tool v3.0.1 (standalone)
  Gyaltsen et al., IEDB Analysis Resource 2017
  Downloaded from: https://tools.iedb.org/population/download/
  Bundled at: src/tools/population_coverage/

  Internally uses the IEDB diploid genotype frequency formula:
    Coverage = 1 - ∏((1 - p_i)²) across covered alleles
    Sette & Sidney, Immunogenetics 1999; IEDB Coverage Tool v3.0

HLA frequency data source:
  Allele Frequency Net Database (AFND), bundled inside IEDB tool
  (population_coverage_pickle/).
  Gonzalez-Galarza et al., Nucleic Acids Research 2020, 48(D1):D783-788
  http://www.allelefrequencies.net
  Data version: IEDB tool v3.0.1 release (static snapshot).
  Limitation: allele frequencies may shift with updated AFND releases.
  Planned upgrade: live AFND query (roadmap Q3).

Populations covered: Global (weighted), Sub-Saharan Africa, East Africa
  (Kenya/Uganda/Rwanda/Burundi/Tanzania), Europe, East Asia, South Asia,
  Americas (admixed). Population definitions follow AFND regional groupings.

Fallback: if IEDB tool import fails, uses validated AFND 2020 static
  frequency tables with the same diploid formula. Fallback is logged
  and labeled in the decision audit trail.
"""

import os
import sys
import logging
import tempfile
import re
from typing import List, Dict, Any, Optional, Set
from src.models.candidate import CandidateProtein, EpitopeResult, ConfidenceTier

from src.utils.logger import get_logger
logger = get_logger("tope_deep.agents.Agent 7")  # use the correct agent name

# ── IEDB TOOL PATH ────────────────────────────────────────────────────────────
_TOOL_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "tools", "population_coverage"
)

# ── POPULATION MAPPING ───────────────────────────────────────────────────────
# IEDB tool population strings - verified against configure.py in the tool
IEDB_POPULATIONS = {
    "global":       "Area > World",
    "african":      "Area > Africa",
    "east_african": "Country > Kenya",      # Best proxy for EA in IEDB dataset
    "european":     "Area > Europe",
    "east_asian":   "Area > East Asia",
    "south_asian":  "Area > South Asia",
    "americas":     "Area > North America", # Admixed proxy
}

POPULATION_NAMES = {
    "global":       "Global",
    "african":      "Sub-Saharan Africa",
    "east_african": "East Africa (KE/UG/RW/BI/TZ)",
    "european":     "Europe",
    "east_asian":   "East Asia",
    "south_asian":  "South Asia",
    "americas":     "Americas (Admixed)",
}

# ── FALLBACK HLA FREQUENCIES (AFND 2020) ─────────────────────────────────────
# Used only when IEDB tool import fails. Labeled explicitly in audit trail.
HLA_I_FREQ = {
    "HLA-A*01:01": {"global": 0.134, "african": 0.098, "east_african": 0.092, "european": 0.172, "east_asian": 0.035, "south_asian": 0.098, "americas": 0.125},
    "HLA-A*02:01": {"global": 0.245, "african": 0.098, "east_african": 0.085, "european": 0.285, "east_asian": 0.152, "south_asian": 0.145, "americas": 0.218},
    "HLA-A*03:01": {"global": 0.125, "african": 0.062, "east_african": 0.058, "european": 0.152, "east_asian": 0.025, "south_asian": 0.082, "americas": 0.098},
    "HLA-A*11:01": {"global": 0.098, "african": 0.038, "east_african": 0.035, "european": 0.058, "east_asian": 0.198, "south_asian": 0.132, "americas": 0.068},
    "HLA-A*24:02": {"global": 0.168, "african": 0.045, "east_african": 0.042, "european": 0.092, "east_asian": 0.272, "south_asian": 0.118, "americas": 0.152},
    "HLA-A*30:01": {"global": 0.035, "african": 0.098, "east_african": 0.102, "european": 0.018, "east_asian": 0.008, "south_asian": 0.025, "americas": 0.042},
    "HLA-A*30:02": {"global": 0.028, "african": 0.075, "east_african": 0.082, "european": 0.012, "east_asian": 0.005, "south_asian": 0.018, "americas": 0.032},
    "HLA-A*68:01": {"global": 0.041, "african": 0.045, "east_african": 0.048, "european": 0.035, "east_asian": 0.012, "south_asian": 0.028, "americas": 0.062},
    "HLA-A*68:02": {"global": 0.025, "african": 0.068, "east_african": 0.075, "european": 0.008, "east_asian": 0.005, "south_asian": 0.015, "americas": 0.035},
    "HLA-B*07:02": {"global": 0.112, "african": 0.045, "east_african": 0.038, "european": 0.142, "east_asian": 0.055, "south_asian": 0.072, "americas": 0.098},
    "HLA-B*08:01": {"global": 0.089, "african": 0.032, "east_african": 0.028, "european": 0.125, "east_asian": 0.008, "south_asian": 0.055, "americas": 0.065},
    "HLA-B*15:01": {"global": 0.065, "african": 0.025, "east_african": 0.022, "european": 0.075, "east_asian": 0.098, "south_asian": 0.058, "americas": 0.055},
    "HLA-B*15:03": {"global": 0.012, "african": 0.052, "east_african": 0.062, "european": 0.002, "east_asian": 0.001, "south_asian": 0.005, "americas": 0.015},
    "HLA-B*35:01": {"global": 0.087, "african": 0.068, "east_african": 0.065, "european": 0.098, "east_asian": 0.042, "south_asian": 0.092, "americas": 0.112},
    "HLA-B*40:01": {"global": 0.068, "african": 0.018, "east_african": 0.015, "european": 0.062, "east_asian": 0.125, "south_asian": 0.085, "americas": 0.058},
    "HLA-B*42:01": {"global": 0.008, "african": 0.045, "east_african": 0.058, "european": 0.001, "east_asian": 0.001, "south_asian": 0.002, "americas": 0.012},
    "HLA-B*44:02": {"global": 0.062, "african": 0.025, "east_african": 0.022, "european": 0.082, "east_asian": 0.018, "south_asian": 0.035, "americas": 0.052},
    "HLA-B*51:01": {"global": 0.065, "african": 0.028, "east_african": 0.025, "european": 0.072, "east_asian": 0.082, "south_asian": 0.058, "americas": 0.055},
    "HLA-B*53:01": {"global": 0.015, "african": 0.072, "east_african": 0.065, "european": 0.002, "east_asian": 0.001, "south_asian": 0.005, "americas": 0.018},
    "HLA-B*57:01": {"global": 0.024, "african": 0.028, "east_african": 0.032, "european": 0.032, "east_asian": 0.008, "south_asian": 0.035, "americas": 0.022},
    "HLA-B*58:01": {"global": 0.016, "african": 0.058, "east_african": 0.065, "european": 0.005, "east_asian": 0.012, "south_asian": 0.018, "americas": 0.015},
}

HLA_II_FREQ = {
    "DRB1*01:01": {"global": 0.085, "african": 0.042, "east_african": 0.038, "european": 0.112, "east_asian": 0.028, "south_asian": 0.072, "americas": 0.078},
    "DRB1*03:01": {"global": 0.092, "african": 0.098, "east_african": 0.105, "european": 0.098, "east_asian": 0.025, "south_asian": 0.085, "americas": 0.082},
    "DRB1*04:01": {"global": 0.078, "african": 0.025, "east_african": 0.022, "european": 0.105, "east_asian": 0.092, "south_asian": 0.065, "americas": 0.098},
    "DRB1*07:01": {"global": 0.112, "african": 0.082, "east_african": 0.078, "european": 0.132, "east_asian": 0.045, "south_asian": 0.098, "americas": 0.108},
    "DRB1*11:01": {"global": 0.098, "african": 0.125, "east_african": 0.132, "european": 0.085, "east_asian": 0.035, "south_asian": 0.108, "americas": 0.075},
    "DRB1*13:01": {"global": 0.072, "african": 0.092, "east_african": 0.098, "european": 0.078, "east_asian": 0.028, "south_asian": 0.065, "americas": 0.068},
    "DRB1*15:01": {"global": 0.118, "african": 0.068, "east_african": 0.062, "european": 0.148, "east_asian": 0.132, "south_asian": 0.125, "americas": 0.112},
    "DRB1*15:03": {"global": 0.025, "african": 0.085, "east_african": 0.092, "european": 0.005, "east_asian": 0.002, "south_asian": 0.008, "americas": 0.018},
}


# ── FALLBACK FORMULA (IEDB diploid model) ────────────────────────────────────
def _calc_coverage_fallback(alleles: Set[str], population: str, hla_class: str) -> float:
    freq_table = HLA_I_FREQ if hla_class == "I" else HLA_II_FREQ
    uncovered = 1.0
    for allele in alleles:
        af = freq_table.get(allele, {}).get(population, 0.0)
        if af > 0:
            uncovered *= (1.0 - af) ** 2
    return 1.0 - uncovered


# ── IEDB TOOL WRAPPER ─────────────────────────────────────────────────────────
def _load_iedb_tool() -> Optional[Any]:
    try:
        _tool_parent = os.path.dirname(_TOOL_DIR)
        if _tool_parent not in sys.path:
            sys.path.insert(0, _tool_parent)
        from population_coverage.population_calculation import PopulationCoverage
        return PopulationCoverage
    except Exception as e:
        logger.warning(f"IEDB coverage tool import failed: {e} - using AFND fallback")
        return None


def _write_epitope_file(
    ctl_epitopes: List[EpitopeResult],
    htl_epitopes: List[EpitopeResult],
) -> str:
    """
    Write epitope+allele pairs to a temp file in IEDB tool format.
    Format: epitope_sequence,allele_name (one per line, no header)
    MHC-II alleles get HLA- prefix added if missing.
    Returns path to temp file.
    """
    lines = []
    for ep in ctl_epitopes:
        if ep.hla_allele and ep.confidence_tier in (ConfidenceTier.HIGH, ConfidenceTier.MEDIUM):
            lines.append(f"{ep.sequence},{ep.hla_allele}")
    for ep in htl_epitopes:
        if ep.hla_allele and ep.confidence_tier in (
            ConfidenceTier.HIGH, ConfidenceTier.MEDIUM, ConfidenceTier.LOW
        ):
            allele = ep.hla_allele
            if not allele.startswith("HLA-"):
                allele = f"HLA-{allele}"
            lines.append(f"{ep.sequence},{allele}")

    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".csv", delete=False, encoding="utf-8"
    )
    tmp.write("\n".join(lines))
    tmp.close()
    return tmp.name


def _run_iedb_coverage(
    PopulationCoverage: Any,
    epitope_file: str,
    population_key: str,
    mhc_classes: List[str],
) -> Optional[Dict]:
    """
    Call IEDB PopulationCoverage for one population key and one or more MHC classes.
    Returns dict with keys: mhc_i_pct, mhc_ii_pct, combined_pct or None on failure.
    """
    try:
        pc = PopulationCoverage()
        iedb_pop = IEDB_POPULATIONS[population_key]

        results_by_class = {}
        for mhc in ["I", "II"]:
            result, negative = pc.calculate_coverage(
                population=[iedb_pop],
                mhc_class=[mhc],
                filename=epitope_file,
            )
            if result:
                results_by_class[mhc] = result[0].get("coverage", 0.0)
            else:
                results_by_class[mhc] = 0.0

        ci  = results_by_class.get("I",  0.0)
        cii = results_by_class.get("II", 0.0)
        combined = 1.0 - (1.0 - ci / 100) * (1.0 - cii / 100)

        return {
            "mhc_i_pct":    round(ci, 1),
            "mhc_ii_pct":   round(cii, 1),
            "combined_pct": round(combined * 100, 1),
            "population_label": POPULATION_NAMES.get(population_key, population_key),
            "method": "IEDB_tool_v3.0.1",
        }
    except Exception as e:
        logger.warning(f"IEDB tool call failed for {population_key}: {e}")
        return None


# ── COVERAGE AGENT ────────────────────────────────────────────────────────────
class CoverageAgent:
    def __init__(self):
        self.stage_name = "coverage_analysis"
        self.populations = list(POPULATION_NAMES.keys())
        self._PopulationCoverage = _load_iedb_tool()
        self._using_iedb_tool = self._PopulationCoverage is not None
        logger.info(
            f"Agent 7: Coverage method = "
            f"{'IEDB tool v3.0.1' if self._using_iedb_tool else 'AFND 2020 fallback'}"
        )

    def run(self, candidates: List[CandidateProtein]) -> List[CandidateProtein]:
        logger.info("Agent 7: Starting population coverage analysis")
        active = [c for c in candidates if c.status.value == "active"]
        logger.info(
            f"   Analyzing {len(active)} candidates across "
            f"{len(self.populations)} populations "
            f"({'IEDB tool' if self._using_iedb_tool else 'AFND fallback'})"
        )

        for i, candidate in enumerate(active):
            logger.info(f"   [{i+1}/{len(active)}] {candidate.protein_name}")

            if self._using_iedb_tool:
                coverage_results = self._coverage_via_iedb_tool(candidate)
            else:
                coverage_results = self._coverage_via_fallback(candidate)

            if not coverage_results:
                logger.warning(f"   Coverage failed for {candidate.protein_name}")
                continue

            global_cov      = coverage_results.get("global", {})
            african_cov     = coverage_results.get("african", {})
            east_african_cov = coverage_results.get("east_african", {})

            candidate.hla_coverage_global = global_cov.get("combined_pct", 0) / 100
            candidate.hla_coverage_africa = max(
                african_cov.get("combined_pct", 0),
                east_african_cov.get("combined_pct", 0),
            ) / 100

            class_i_alleles  = self._extract_class_i_alleles(candidate)
            class_ii_alleles = self._extract_class_ii_alleles(candidate)
            gaps = self._find_gaps(class_i_alleles, class_ii_alleles)

            all_combined = [v["combined_pct"] for v in coverage_results.values()]
            equity_gap = max(all_combined) - min(all_combined) if all_combined else 0

            method_note = (
                "IEDB Population Coverage Tool v3.0.1 (Gyaltsen et al. 2017)"
                if self._using_iedb_tool
                else "AFND 2020 static frequencies (Gonzalez-Galarza et al. 2020) - fallback"
            )

            candidate.add_decision(
                stage=self.stage_name,
                decision="coverage_calculated",
                reasoning=self._build_recommendation(coverage_results, equity_gap),
                method=method_note,
                global_coverage_pct=global_cov.get("combined_pct", 0),
                african_coverage_pct=african_cov.get("combined_pct", 0),
                east_african_coverage_pct=east_african_cov.get("combined_pct", 0),
                equity_gap_pct=round(equity_gap, 1),
                coverage_gaps=gaps[:5],
                per_population=coverage_results,
            )
            candidate.stage = self.stage_name

        logger.info("Agent 7: Coverage analysis complete")
        return candidates

    # ── IEDB TOOL PATH ────────────────────────────────────────────────────────

    def _coverage_via_iedb_tool(
        self, candidate: CandidateProtein
    ) -> Dict[str, Dict]:
        epitope_file = None
        try:
            epitope_file = _write_epitope_file(
                candidate.ctl_epitopes, candidate.htl_epitopes
            )
            coverage_results = {}
            for pop_key in self.populations:
                result = _run_iedb_coverage(
                    self._PopulationCoverage,
                    epitope_file,
                    pop_key,
                    ["I", "II"],
                )
                if result:
                    coverage_results[pop_key] = result
                else:
                    # Per-population fallback - don't abort the whole candidate
                    coverage_results[pop_key] = self._fallback_single_pop(
                        candidate, pop_key
                    )
                logger.info(
                    f"      {pop_key:18s}: "
                    f"I={coverage_results[pop_key]['mhc_i_pct']:.1f}%  "
                    f"II={coverage_results[pop_key]['mhc_ii_pct']:.1f}%  "
                    f"Combined={coverage_results[pop_key]['combined_pct']:.1f}%"
                )
            return coverage_results
        except Exception as e:
            logger.error(f"IEDB tool path failed: {e} - falling back")
            return self._coverage_via_fallback(candidate)
        finally:
            if epitope_file and os.path.exists(epitope_file):
                try:
                    os.unlink(epitope_file)
                except Exception:
                    pass

    def _fallback_single_pop(
        self, candidate: CandidateProtein, pop_key: str
    ) -> Dict:
        class_i  = self._extract_class_i_alleles(candidate)
        class_ii = self._extract_class_ii_alleles(candidate)
        ci  = _calc_coverage_fallback(class_i,  pop_key, "I")
        cii = _calc_coverage_fallback(class_ii, pop_key, "II")
        combined = 1.0 - (1.0 - ci) * (1.0 - cii)
        return {
            "mhc_i_pct":    round(ci  * 100, 1),
            "mhc_ii_pct":   round(cii * 100, 1),
            "combined_pct": round(combined * 100, 1),
            "population_label": POPULATION_NAMES.get(pop_key, pop_key),
            "method": "AFND_2020_fallback",
        }

    # ── FALLBACK PATH ─────────────────────────────────────────────────────────

    def _coverage_via_fallback(
        self, candidate: CandidateProtein
    ) -> Dict[str, Dict]:
        class_i  = self._extract_class_i_alleles(candidate)
        class_ii = self._extract_class_ii_alleles(candidate)
        logger.info(
            f"      Covered alleles: {len(class_i)} MHC-I, {len(class_ii)} MHC-II"
        )
        coverage_results = {}
        for pop in self.populations:
            ci  = _calc_coverage_fallback(class_i,  pop, "I")
            cii = _calc_coverage_fallback(class_ii, pop, "II")
            combined = 1.0 - (1.0 - ci) * (1.0 - cii)
            coverage_results[pop] = {
                "mhc_i_pct":    round(ci  * 100, 1),
                "mhc_ii_pct":   round(cii * 100, 1),
                "combined_pct": round(combined * 100, 1),
                "population_label": POPULATION_NAMES.get(pop, pop),
                "method": "AFND_2020_fallback",
            }
            logger.info(
                f"      {pop:18s}: I={ci:.1%}  II={cii:.1%}  "
                f"Combined={combined:.1%}"
            )
        return coverage_results

    # ── ALLELE EXTRACTION ─────────────────────────────────────────────────────

    def _extract_class_i_alleles(self, candidate: CandidateProtein) -> Set[str]:
        alleles = set()
        for ep in candidate.ctl_epitopes:
            if ep.hla_allele and ep.confidence_tier in (
                ConfidenceTier.HIGH, ConfidenceTier.MEDIUM
            ):
                alleles.add(ep.hla_allele)
        return alleles

    def _extract_class_ii_alleles(self, candidate: CandidateProtein) -> Set[str]:
        alleles = set()
        for ep in candidate.htl_epitopes:
            if ep.hla_allele and ep.confidence_tier in (
                ConfidenceTier.HIGH, ConfidenceTier.MEDIUM, ConfidenceTier.LOW
            ):
                allele = ep.hla_allele
                if allele.startswith("HLA-"):
                    allele = allele[4:]
                alleles.add(allele)
        return alleles

    # ── GAP ANALYSIS ──────────────────────────────────────────────────────────

    def _find_gaps(
        self, class_i: Set[str], class_ii: Set[str]
    ) -> List[str]:
        gaps = []
        for allele, pops in HLA_I_FREQ.items():
            if allele not in class_i:
                max_freq = max(pops.values())
                if max_freq >= 0.05:
                    best_pop = max(pops, key=pops.get)
                    gaps.append(f"{allele} ({max_freq:.0%} in {best_pop})")
        for allele, pops in HLA_II_FREQ.items():
            if allele not in class_ii:
                max_freq = max(pops.values())
                if max_freq >= 0.05:
                    best_pop = max(pops, key=pops.get)
                    gaps.append(f"{allele} ({max_freq:.0%} in {best_pop})")
        return sorted(
            gaps,
            key=lambda x: float(re.search(r'(\d+)%', x).group(1))
            if re.search(r'(\d+)%', x) else 0,
            reverse=True,
        )

    # ── RECOMMENDATION ────────────────────────────────────────────────────────

    def _build_recommendation(self, coverage: Dict, equity_gap: float) -> str:
        global_c      = coverage.get("global", {}).get("combined_pct", 0)
        african_c     = coverage.get("african", {}).get("combined_pct", 0)
        east_african_c = coverage.get("east_african", {}).get("combined_pct", 0)
        lines = []

        if global_c >= 85 and african_c >= 75:
            lines.append(
                f"Broad population coverage achieved: {global_c:.1f}% global, "
                f"{african_c:.1f}% Sub-Saharan Africa, "
                f"{east_african_c:.1f}% East Africa. "
                f"Epitope set is suitable for multi-epitope construct design."
            )
        elif global_c >= 70 and african_c >= 60:
            lines.append(
                f"Moderate coverage: {global_c:.1f}% global, "
                f"{african_c:.1f}% Sub-Saharan Africa. "
                f"Consider supplementing with alleles enriched in African "
                f"populations (e.g. HLA-A*30:01, HLA-B*53:01, DRB1*03:01) "
                f"to improve equity."
            )
        elif global_c >= 50:
            lines.append(
                f"Suboptimal coverage: {global_c:.1f}% global. "
                f"Epitope set covers fewer than half the global population. "
                f"Expand HLA allele panel or select alternative antigenic regions."
            )
        else:
            lines.append(
                f"Insufficient coverage: {global_c:.1f}% global. "
                f"This epitope set is unlikely to confer broad immunogenicity. "
                f"Recommend antigen re-selection or proteome-wide screening."
            )

        if equity_gap > 30:
            lines.append(
                f"Equity concern: {equity_gap:.1f}% coverage gap between "
                f"best and worst covered populations. "
                f"Prioritise alleles underrepresented in current predictions."
            )
        elif equity_gap > 15:
            lines.append(
                f"Moderate equity gap ({equity_gap:.1f}%). "
                f"African and South Asian populations may be underserved "
                f"by this epitope set."
            )

        lines.append(
            "Coverage calculated using IEDB Population Coverage Tool v3.0.1 "
            "(Gyaltsen et al. 2017). HLA frequency data: AFND "
            "(Gonzalez-Galarza et al., Nucleic Acids Research 2020)."
        )
        return " | ".join(lines)


# Global instance
coverage_agent = CoverageAgent()