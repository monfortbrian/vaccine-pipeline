"""
TOPE_DEEP CLI Runner

Direct pipeline execution without the HTTP API.
Useful for local testing, batch processing, and CI benchmarks.

Usage:
  python run_pipeline.py --demo
  python run_pipeline.py --protein P9WNK7
  python run_pipeline.py --sequence MTEQQWNFAG... --name "ESAT-6"
  python run_pipeline.py --protein P9WNK7 --no-literature --no-experiment
"""

import sys, os, json, csv, uuid, time, logging, argparse, requests
from datetime import datetime
from pathlib import Path
from typing import List, Optional

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(message)s", datefmt="%H:%M:%S")
logger = logging.getLogger("tope_deep.cli")

from src.models.candidate import CandidateProtein, CandidateStatus, ConfidenceTier, EpitopeResult
from src.orchestrator import PipelineOrchestrator


DEMO = CandidateProtein(
    protein_id="P9WNK7",
    protein_name="ESAT-6 Mycobacterium tuberculosis",
    sequence="MTEQQWNFAGIEAAASAIQGNVTSIHSLLDEGKQSLTKLAAAWGGSGSEAYQGVQQKWDATATELNNALQNLARTISEAGQAMASTEGNVTGMFA",
    source="uniprot",
    stage="antigen_screening",
    status=CandidateStatus.ACTIVE,
    vaxijen_score=0.87,
    psortb_localization="secreted",
)


def fetch_uniprot(uid: str) -> Optional[str]:
    try:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=15)
        r.raise_for_status()
        lines = r.text.strip().split("\n")
        return "".join(l for l in lines if not l.startswith(">"))
    except Exception as e:
        logger.error(f"UniProt fetch failed: {e}")
        return None


def ep_row(c: CandidateProtein, ep: EpitopeResult, t: str) -> dict:
    cats = ep.tool_outputs.get("model_categories", ["HUMAN"])
    return {
        "protein_id":       c.protein_id,
        "protein_name":     c.protein_name,
        "sequence":         ep.sequence,
        "type":             t,
        "hla_allele":       ep.hla_allele or "",
        "ic50_nm":          ep.ic50_nm or "",
        "percentile_rank":  ep.percentile_rank or "",
        "confidence":       ep.confidence_tier.value,
        "model_categories": ",".join(cats),
        "allergenicity_safe": ep.allergenicity_safe,
        "toxicity_safe":      ep.toxicity_safe,
    }


def main():
    p = argparse.ArgumentParser(description="TOPE_DEEP - Agentic Vaccine Discovery Pipeline")
    p.add_argument("--demo",           action="store_true",  help="Run ESAT-6 demo (P9WNK7)")
    p.add_argument("--protein",        type=str,             help="UniProt accession")
    p.add_argument("--sequence",       type=str,             help="Amino acid sequence")
    p.add_argument("--name",           type=str,             help="Protein name (used with --sequence)")
    p.add_argument("--no-safety",      action="store_true",  help="Skip Agent 6 safety screening")
    p.add_argument("--no-coverage",    action="store_true",  help="Skip Agent 7 population coverage")
    p.add_argument("--no-literature",  action="store_true",  help="Skip Agent 9 literature agent")
    p.add_argument("--no-experiment",  action="store_true",  help="Skip Agent 10 experiment planner")
    p.add_argument("--output",         type=str, default="results", help="Output directory")
    args = p.parse_args()

    candidates = []

    if args.demo or (not args.protein and not args.sequence):
        candidates = [DEMO]
        print("\nTOPE_DEEP Demo run: ESAT-6 (P9WNK7), Mycobacterium tuberculosis")

    elif args.sequence:
        seq = args.sequence.upper().replace(" ", "").replace("\n", "")
        candidates = [CandidateProtein(
            protein_id=args.name or "user_input",
            protein_name=args.name or "Custom protein",
            sequence=seq, source="user_input",
            stage="antigen_screening", status=CandidateStatus.ACTIVE,
        )]

    elif args.protein:
        seq = fetch_uniprot(args.protein)
        if not seq:
            print(f"ERROR: Could not fetch {args.protein} from UniProt"); sys.exit(1)
        candidates = [CandidateProtein(
            protein_id=args.protein,
            protein_name=args.name or args.protein,
            sequence=seq, source="uniprot",
            stage="antigen_screening", status=CandidateStatus.ACTIVE,
        )]

    run_id = str(uuid.uuid4())
    config = {
        "organism_class":  "bacteria",
        "input_type":      "demo" if args.demo else ("uniprot_id" if args.protein else "sequence"),
        "run_safety":      not args.no_safety,
        "run_coverage":    not args.no_coverage,
        "run_literature":  not args.no_literature,
        "run_experiment":  not args.no_experiment,
        "lab_constraints": "standard academic lab",
    }

    print(f"\nRun ID: {run_id}")
    print(f"Agents: Agent 1-Agent 10 | Safety={'yes' if config['run_safety'] else 'no'} | Coverage={'yes' if config['run_coverage'] else 'no'}")
    print("-" * 60)

    o = PipelineOrchestrator()

    def cb(node, pct, msg):
        bar = "█" * int(pct * 20) + "░" * (20 - int(pct * 20))
        print(f"\r  [{bar}] {int(pct*100):3d}%  {node}: {msg[:50]:<50}", end="", flush=True)

    result = o.run(run_id, candidates, config, progress_callback=cb)
    print()

    # Output
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")

    json_path = out / f"TD_{ts}_{run_id[:8]}_results.json"
    csv_path  = out / f"TD_{ts}_{run_id[:8]}_epitopes.csv"

    with open(json_path, "w") as f:
        json.dump(result, f, indent=2, default=str)

    rows = []
    for c in result["candidates"]:
        for ep in c.ctl_epitopes:   rows.append(ep_row(c, ep, "CTL"))
        for ep in c.htl_epitopes:   rows.append(ep_row(c, ep, "HTL"))
        for ep in c.bcell_epitopes: rows.append(ep_row(c, ep, "B-cell"))

    if rows:
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)

    t = result["timing"]
    print(f"\n{'─'*60}")
    print(f"  TOPE_DEEP Run complete")
    print(f"  Total: {t['total_seconds']}s")
    for k, v in t.items():
        if k != "total_seconds" and v is not None:
            print(f"    {k:<18} {v:.2f}s")
    print(f"\n  Results: {json_path}")
    print(f"  CSV:     {csv_path}")
    print(f"{'─'*60}\n")


if __name__ == "__main__":
    main()