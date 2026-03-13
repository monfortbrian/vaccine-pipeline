"""
KOZI MVP-2 PIPELINE ORCHESTRATOR
Chains: N3 (T-cell) -> N4 (B-cell) -> N6 (Safety) -> N7 (Coverage)

Usage:
    python run_pipeline.py                           # Demo: M. tuberculosis ESAT-6
    python run_pipeline.py --protein P9WNK7          # UniProt ID
    python run_pipeline.py --sequence MTEQQWNF...    # Raw sequence
    python run_pipeline.py --name "My Protein" --sequence MTEQ...

Output:
    results/<timestamp>_results.json
    results/<timestamp>_epitopes.csv
"""

from src.agents.predictors.coverage_agent import CoverageAgent
from src.agents.predictors.safety_filter import SafetyFilterAgent
from src.agents.predictors.bcell_predictor import BCellPredictorAgent
from src.agents.predictors.tcell_predictor import TCellPredictorAgent
from src.models.candidate import (
    CandidateProtein, CandidateStatus, PipelineRun,
    EpitopeResult, EpitopeType, ConfidenceTier,
)
import json
import csv
import logging
import time
import sys
import os
import argparse
import requests
from datetime import datetime
from typing import List, Dict, Any, Optional
from pathlib import Path

# Ensure project root is on path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# Logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("kozi.pipeline")


# - DEMO PROTEIN-──────────────────────────────────────────────────────

DEMO_PROTEINS = [
    CandidateProtein(
        protein_id="P9WNK7",
        protein_name="ESAT-6 (EsxA) - M. tuberculosis",
        sequence=(
            "MTEQQWNFAGIEAAASAIQGNVTSIHSLLDEGKQSLTKLAAAWGGSGSEAYQGVQQKWD"
            "ATATELNNALQNLARTISEAGQAMASTEGNVTGMFA"
        ),
        source="uniprot",
        stage="antigen_screening",
        status=CandidateStatus.ACTIVE,
        vaxijen_score=0.87,
        psortb_localization="secreted",
        tmhmm_helices=0,
    ),
]


# - PIPELINE-──────────────────────────────────────────────────────────

class KoziMVP2Pipeline:
    """
    Orchestrates: N3 T-Cell -> N4 B-Cell -> N6 Safety -> N7 Coverage
    """

    def __init__(self, output_dir: str = "results"):
        self.n3 = TCellPredictorAgent()
        self.n4 = BCellPredictorAgent()
        self.n6 = SafetyFilterAgent()
        self.n7 = CoverageAgent()
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run(self, candidates: List[CandidateProtein]) -> Dict[str, Any]:
        """Execute the full MVP-2 pipeline."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        start = time.time()

        print()
        print("=" * 62)
        print("   KOZI AI - MVP-2 Vaccine Discovery Pipeline")
        print(f"   {len(candidates)} candidate protein(s)")
        print(f"   Nodes: N3 (T-cell) → N4 (B-cell) → N6 (Safety) → N7 (Coverage)")
        print("=" * 62)

        # - N3: T-Cell Prediction-
        print("\nNODE N3: T-Cell Epitope Prediction-")
        n3_start = time.time()
        candidates = self.n3.run(candidates)
        n3_time = time.time() - n3_start
        print(f"   N3 completed in {n3_time:.1f}s")

        # - N4: B-Cell Prediction-
        print("\nNODE N4: B-Cell Epitope Prediction-")
        n4_start = time.time()
        candidates = self.n4.run(candidates)
        n4_time = time.time() - n4_start
        print(f"   N4 completed in {n4_time:.1f}s")

        # - N6: Safety Filter-
        print("\nNODE N6: Safety Screening-")
        n6_start = time.time()
        candidates = self.n6.run(candidates)
        n6_time = time.time() - n6_start
        print(f"   N6 completed in {n6_time:.1f}s")

        # - N7: Coverage Analysis-
        print("\nNODE N7: Population Coverage-")
        n7_start = time.time()
        candidates = self.n7.run(candidates)
        n7_time = time.time() - n7_start
        print(f"   N7 completed in {n7_time:.1f}s")

        total_time = time.time() - start

        # - Build results-
        results = self._build_results(candidates, timestamp, total_time,
                                      n3_time, n4_time, n6_time, n7_time)

        # - Export-
        json_path = self.output_dir / f"{timestamp}_results.json"
        csv_path = self.output_dir / f"{timestamp}_epitopes.csv"

        with open(json_path, "w") as f:
            json.dump(results, f, indent=2, default=str)

        self._export_csv(candidates, csv_path)

        # - Summary-
        self._print_summary(results, json_path, csv_path)

        return results

    def _build_results(self, candidates, timestamp, total, n3, n4, n6, n7):
        """Build complete results dict."""
        summary_candidates = []

        for c in candidates:
            ctl_count = len(c.ctl_epitopes)
            htl_count = len(c.htl_epitopes)
            bcell_count = len(c.bcell_epitopes)

            ctl_strong = len(
                [e for e in c.ctl_epitopes if e.confidence_tier == ConfidenceTier.HIGH])
            htl_strong = len(
                [e for e in c.htl_epitopes if e.confidence_tier == ConfidenceTier.HIGH])
            bcell_high = len(
                [e for e in c.bcell_epitopes if e.confidence_tier == ConfidenceTier.HIGH])

            # Count safety results
            all_eps = list(c.ctl_epitopes) + \
                list(c.htl_epitopes) + list(c.bcell_epitopes)
            safe = len([e for e in all_eps if e.allergenicity_safe is True])
            flagged = len(
                [e for e in all_eps if e.tool_outputs.get("safety_flags")])

            summary_candidates.append({
                "protein_id": c.protein_id,
                "protein_name": c.protein_name,
                "sequence_length": len(c.sequence),
                "epitopes": {
                    "ctl": {"total": ctl_count, "strong_binders": ctl_strong},
                    "htl": {"total": htl_count, "strong_binders": htl_strong},
                    "bcell": {"total": bcell_count, "high_confidence": bcell_high},
                },
                "safety": {"safe": safe, "flagged": flagged},
                "coverage": {
                    "global_pct": round((c.hla_coverage_global or 0) * 100, 1),
                    "african_pct": round((c.hla_coverage_africa or 0) * 100, 1),
                },
                "decisions": c.decisions,
                "all_epitopes": {
                    "ctl": [self._ep_to_dict(e) for e in c.ctl_epitopes],
                    "htl": [self._ep_to_dict(e) for e in c.htl_epitopes],
                    "bcell": [self._ep_to_dict(e) for e in c.bcell_epitopes],
                },
            })

        return {
            "pipeline": "Kozi AI MVP-2",
            "version": "2.0.0",
            "timestamp": timestamp,
            "timing": {
                "total_seconds": round(total, 1),
                "n3_tcell": round(n3, 1),
                "n4_bcell": round(n4, 1),
                "n6_safety": round(n6, 1),
                "n7_coverage": round(n7, 1),
            },
            "candidates": summary_candidates,
        }

    def _ep_to_dict(self, ep: EpitopeResult) -> Dict:
        return {
            "sequence": ep.sequence,
            "type": ep.epitope_type.value,
            "hla_allele": ep.hla_allele,
            "ic50_nm": ep.ic50_nm,
            "percentile_rank": ep.percentile_rank,
            "confidence": ep.confidence_tier.value,
            "allergenicity_safe": ep.allergenicity_safe,
            "toxicity_safe": ep.toxicity_safe,
        }

    def _export_csv(self, candidates: List[CandidateProtein], path: Path):
        """Export all epitopes to CSV for lab teams."""
        rows = []
        for c in candidates:
            for ep in c.ctl_epitopes:
                rows.append(self._csv_row(c, ep, "CTL"))
            for ep in c.htl_epitopes:
                rows.append(self._csv_row(c, ep, "HTL"))
            for ep in c.bcell_epitopes:
                rows.append(self._csv_row(c, ep, "B-cell"))

        if rows:
            fieldnames = ["protein", "protein_id", "epitope_sequence", "type",
                          "hla_allele", "ic50_nm", "confidence",
                          "allergenicity_safe", "toxicity_safe"]
            with open(path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)

    def _csv_row(self, cand: CandidateProtein, ep: EpitopeResult, ep_type: str) -> Dict:
        return {
            "protein": cand.protein_name,
            "protein_id": cand.protein_id,
            "epitope_sequence": ep.sequence,
            "type": ep_type,
            "hla_allele": ep.hla_allele or "N/A",
            "ic50_nm": ep.ic50_nm or "",
            "confidence": ep.confidence_tier.value,
            "allergenicity_safe": ep.allergenicity_safe,
            "toxicity_safe": ep.toxicity_safe,
        }

    def _print_summary(self, results, json_path, csv_path):
        """Print final summary to terminal."""
        print()
        print("=" * 62)
        print("   PIPELINE COMPLETE")
        print("=" * 62)

        for c in results["candidates"]:
            ep = c["epitopes"]
            cov = c["coverage"]
            saf = c["safety"]
            print(f"\n   {c['protein_name']} ({c['sequence_length']}aa)")
            print(f"   {'─' * 50}")
            print(
                f"   CTL epitopes:     {ep['ctl']['total']:>3}  ({ep['ctl']['strong_binders']} strong binders)")
            print(
                f"   HTL epitopes:     {ep['htl']['total']:>3}  ({ep['htl']['strong_binders']} strong binders)")
            print(
                f"   B-cell epitopes:  {ep['bcell']['total']:>3}  ({ep['bcell']['high_confidence']} high confidence)")
            print(
                f"   Safety:           {saf['safe']} safe, {saf['flagged']} flagged")
            print(f"   Global coverage:  {cov['global_pct']}%")
            print(f"   African coverage: {cov['african_pct']}%")

        t = results["timing"]
        print(f"\n   Time: {t['total_seconds']}s total "
              f"(N3:{t['n3_tcell']}s N4:{t['n4_bcell']}s N6:{t['n6_safety']}s N7:{t['n7_coverage']}s)")
        print(f"\n   Results: {json_path}")
        print(f"   CSV:     {csv_path}")
        print("=" * 62)
        print()


# - CLI-───────────────────────────────────────────────────────────────

def fetch_uniprot(uniprot_id: str) -> Optional[str]:
    """Fetch sequence from UniProt REST API."""
    try:
        resp = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta",
            timeout=15,
        )
        resp.raise_for_status()
        lines = resp.text.strip().split("\n")
        return "".join(l.strip() for l in lines if not l.startswith(">"))
    except Exception as e:
        logger.error(f"UniProt fetch failed: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Kozi AI MVP-2 Pipeline")
    parser.add_argument("--demo", action="store_true",
                        help="Run demo with M. tuberculosis ESAT-6")
    parser.add_argument("--protein", type=str,
                        help="UniProt ID to analyze")
    parser.add_argument("--name", type=str, default="input_protein",
                        help="Protein name (used with --sequence or --protein)")
    parser.add_argument("--sequence", type=str,
                        help="Raw amino acid sequence")
    parser.add_argument("--output-dir", type=str, default="results",
                        help="Output directory")

    args = parser.parse_args()
    candidates = []

    if args.demo or (not args.protein and not args.sequence):
        candidates = DEMO_PROTEINS
        print("Running DEMO: M. tuberculosis ESAT-6")

    elif args.sequence:
        candidates = [CandidateProtein(
            protein_id=args.name,
            protein_name=args.name,
            sequence=args.sequence.upper().replace(" ", ""),
            source="user_input",
            stage="antigen_screening",
            status=CandidateStatus.ACTIVE,
        )]

    elif args.protein:
        seq = fetch_uniprot(args.protein)
        if not seq:
            print(f"ERROR: Could not fetch {args.protein} from UniProt")
            sys.exit(1)
        candidates = [CandidateProtein(
            protein_id=args.protein,
            protein_name=args.name or args.protein,
            sequence=seq,
            source="uniprot",
            stage="antigen_screening",
            status=CandidateStatus.ACTIVE,
        )]

    pipeline = KoziMVP2Pipeline(output_dir=args.output_dir)
    pipeline.run(candidates)


if __name__ == "__main__":
    main()
