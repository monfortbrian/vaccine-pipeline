"""
COVERAGE AGENT - MVP-2 NODE N7
Calculates immune coverage across global HLA populations.

Uses the standard IEDB diploid genotype frequency formula:
  Coverage = 1 - product((1 - allele_freq)^2) for each covered allele

Data: HLA frequencies from Allele Frequency Net Database (AFND)
      Gonzalez-Galarza et al., Nucleic Acids Research, 2020

Populations: Global, Sub-Saharan Africa, East Africa, Europe,
             East Asia, South Asia, Americas
"""

import logging
import re
from typing import List, Dict, Any, Optional, Set
from src.models.candidate import CandidateProtein, EpitopeResult, ConfidenceTier

logger = logging.getLogger(__name__)


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
    "DRB1*01:01":  {"global": 0.085, "african": 0.042, "east_african": 0.038, "european": 0.112, "east_asian": 0.028, "south_asian": 0.072, "americas": 0.078},
    "DRB1*03:01":  {"global": 0.092, "african": 0.098, "east_african": 0.105, "european": 0.098, "east_asian": 0.025, "south_asian": 0.085, "americas": 0.082},
    "DRB1*04:01":  {"global": 0.078, "african": 0.025, "east_african": 0.022, "european": 0.105, "east_asian": 0.092, "south_asian": 0.065, "americas": 0.098},
    "DRB1*07:01":  {"global": 0.112, "african": 0.082, "east_african": 0.078, "european": 0.132, "east_asian": 0.045, "south_asian": 0.098, "americas": 0.108},
    "DRB1*11:01":  {"global": 0.098, "african": 0.125, "east_african": 0.132, "european": 0.085, "east_asian": 0.035, "south_asian": 0.108, "americas": 0.075},
    "DRB1*13:01":  {"global": 0.072, "african": 0.092, "east_african": 0.098, "european": 0.078, "east_asian": 0.028, "south_asian": 0.065, "americas": 0.068},
    "DRB1*15:01":  {"global": 0.118, "african": 0.068, "east_african": 0.062, "european": 0.148, "east_asian": 0.132, "south_asian": 0.125, "americas": 0.112},
    "DRB1*15:03":  {"global": 0.025, "african": 0.085, "east_african": 0.092, "european": 0.005, "east_asian": 0.002, "south_asian": 0.008, "americas": 0.018},
}

POPULATION_NAMES = {
    "global": "Global",
    "african": "Sub-Saharan Africa",
    "east_african": "East Africa (KE/UG/RW/TZ)",
    "european": "Europe",
    "east_asian": "East Asia",
    "south_asian": "South Asia",
    "americas": "Americas (Admixed)",
}


def calc_coverage(alleles: Set[str], population: str, hla_class: str = "I") -> float:
    """
    IEDB standard diploid coverage formula.
    Coverage = 1 - product( (1 - af)^2 ) for each covered allele
    """
    freq_table = HLA_I_FREQ if hla_class == "I" else HLA_II_FREQ
    uncovered = 1.0
    for allele in alleles:
        af = freq_table.get(allele, {}).get(population, 0.0)
        if af > 0:
            phenotype_freq = 1.0 - (1.0 - af) ** 2
            uncovered *= (1.0 - phenotype_freq)
    return 1.0 - uncovered


class CoverageAgent:
    def __init__(self):
        self.stage_name = "coverage_analysis"
        self.populations = list(POPULATION_NAMES.keys())

    def run(self, candidates: List[CandidateProtein]) -> List[CandidateProtein]:
        """Calculate population coverage for each candidate."""
        logger.info("N7: Starting population coverage analysis")

        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"   Analyzing {len(active)} candidates across {len(self.populations)} populations")

        for i, candidate in enumerate(active):
            logger.info(f"   [{i+1}/{len(active)}] {candidate.protein_name}")

            # Extract MHC-I alleles from CTL epitopes (HIGH and MEDIUM confidence)
            class_i_alleles = set()
            for ep in candidate.ctl_epitopes:
                if ep.hla_allele and ep.confidence_tier in (ConfidenceTier.HIGH, ConfidenceTier.MEDIUM):
                    class_i_alleles.add(ep.hla_allele)

            # Extract MHC-II alleles from HTL epitopes (include LOW - MHC-II ranks are naturally higher)
            class_ii_alleles = set()
            for ep in candidate.htl_epitopes:
                if ep.hla_allele and ep.confidence_tier in (ConfidenceTier.HIGH, ConfidenceTier.MEDIUM, ConfidenceTier.LOW):
                    allele = ep.hla_allele
                    if allele.startswith("HLA-"):
                        allele = allele[4:]
                    class_ii_alleles.add(allele)

            logger.info(f"      Covered alleles: {len(class_i_alleles)} MHC-I, {len(class_ii_alleles)} MHC-II")

            # Calculate per-population coverage
            coverage_results = {}
            for pop in self.populations:
                ci = calc_coverage(class_i_alleles, pop, "I")
                cii = calc_coverage(class_ii_alleles, pop, "II")
                combined = 1.0 - (1.0 - ci) * (1.0 - cii)

                coverage_results[pop] = {
                    "mhc_i_pct": round(ci * 100, 1),
                    "mhc_ii_pct": round(cii * 100, 1),
                    "combined_pct": round(combined * 100, 1),
                    "population_label": POPULATION_NAMES.get(pop, pop),
                }

                logger.info(f"      {pop:18s}: I={ci:.1%}  II={cii:.1%}  Combined={combined:.1%}")

            # Set candidate-level coverage fields
            global_cov = coverage_results.get("global", {})
            african_cov = coverage_results.get("african", {})
            east_african_cov = coverage_results.get("east_african", {})

            candidate.hla_coverage_global = global_cov.get("combined_pct", 0) / 100
            candidate.hla_coverage_africa = max(
                african_cov.get("combined_pct", 0),
                east_african_cov.get("combined_pct", 0),
            ) / 100

            # Find coverage gaps
            gaps = self._find_gaps(class_i_alleles, class_ii_alleles)

            # Equity analysis
            all_combined = [v["combined_pct"] for v in coverage_results.values()]
            equity_gap = max(all_combined) - min(all_combined) if all_combined else 0

            candidate.add_decision(
                stage=self.stage_name,
                decision="coverage_calculated",
                reasoning=self._build_recommendation(coverage_results, equity_gap),
                global_coverage_pct=global_cov.get("combined_pct", 0),
                african_coverage_pct=african_cov.get("combined_pct", 0),
                east_african_coverage_pct=east_african_cov.get("combined_pct", 0),
                equity_gap_pct=round(equity_gap, 1),
                coverage_gaps=gaps[:5],
                per_population=coverage_results,
            )

            candidate.stage = self.stage_name

        logger.info("N7: Coverage analysis complete")
        return candidates

    def _find_gaps(self, class_i: Set[str], class_ii: Set[str]) -> List[str]:
        """Find high-frequency alleles that are NOT covered."""
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
            key=lambda x: float(re.search(r'(\d+)%', x).group(1)) if re.search(r'(\d+)%', x) else 0,
            reverse=True,
        )

    def _build_recommendation(self, coverage: Dict, equity_gap: float) -> str:
        """Generate actionable recommendation."""
        global_c = coverage.get("global", {}).get("combined_pct", 0)
        african_c = coverage.get("african", {}).get("combined_pct", 0)

        if global_c >= 90 and african_c >= 80:
            return f"EXCELLENT: {global_c:.0f}% global + {african_c:.0f}% African. Ready for construct design."
        elif global_c >= 70 and african_c >= 60:
            return f"GOOD: {global_c:.0f}% global, {african_c:.0f}% African. Consider adding African-enriched alleles."
        elif global_c >= 50:
            return f"MODERATE: {global_c:.0f}% global. Needs more allele coverage for robust design."
        else:
            return f"LOW: {global_c:.0f}% global. Consider alternative antigens or expanded HLA panel."


# Global instance
coverage_agent = CoverageAgent()