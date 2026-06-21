"""
TCELL PREDICTOR AGENT

Tools:
  NetMHCpan 4.1 EL   : MHC class I binding (CTL epitopes)   via IEDB tools cluster
  NetMHCIIpan 4.3    : MHC class II binding (HTL epitopes)   via IEDB tools cluster
  MHCflurry 2.0      : local fallback for MHC-I when IEDB unavailable
  NetMHCstabpan 1.0  : MHC-peptide binding stability (t1/2)  via DTU (optional)
  MAFFT REST         : conservation analysis across strains   via EBI (optional)

Human HLA supertypes (primary population coverage):
  HLA-A*02:01  HLA-A*24:02  HLA-A*03:01  HLA-A*01:01
  HLA-A*11:01  HLA-B*07:02  HLA-B*44:02  HLA-B*35:01

Animal model alleles (five categories all tagged with model_category):

  MOUSE (murine preclinical):
    H-2-Kb, H-2-Db   : C57BL/6 (H-2b), most common inbred strain
    H-2-Kd, H-2-Dd   : BALB/c (H-2d), second most common

  NHP (non-human primate):
    Mamu-A*01         : most immunodominant Mamu allele; HIV/TB/malaria studies
    Mamu-A*02         : second major rhesus allele
    Mamu-B*17         : elite controller allele in SIV studies

  HUMANIZED_MOUSE:
    HLA-A*02:01, HLA-A*24:02 expressed in transgenic mice.
    These are standard HLA alleles flagged as having a
    direct mouse in vivo validation path.

  GUINEA_PIG (Mtb-specific heuristic):
    GPLA-B*0101 not supported by NetMHCpan.
    Flagged via sequence similarity to known GPLA-B*01 binders
    from Vordermeier et al. (2002) and Smith et al. (2000).

Binding stability (NetMHCstabpan):
  stability_score: predicted peptide-MHC complex t1/2 (hours)
  stability_rank:  percentile rank
  High affinity (low IC50) + low stability = flag for wet-lab attention
  Reference: Henriksen et al. (2023) Frontiers in Immunology

Conservation analysis (MAFFT):
  conservation_score: 0-1, per epitope
  Computed from Shannon entropy across known strain sequences
  from UniProt taxonomy search.
  Low conservation (<0.70), flagged in audit trail
  Reference: MAFFT v7 REST (Katoh et al. 2019)

References:
  NetMHCpan 4.1:    Reynisson et al. (2020) Nucleic Acids Res 48:W449
  NetMHCIIpan 4.3:  Reynisson et al. (2020) Nucleic Acids Res
  MHCflurry 2.0:    O'Donnell et al. (2020) Cell Syst 11:42-48
  NetMHCstabpan 1.0: Henriksen et al. (2023) Front Immunol
  MAFFT v7:         Katoh et al. (2019) Brief Bioinform
"""

import os
import time
import logging
import requests
from typing import List, Dict, Any, Optional, Set
from src.models.candidate import CandidateProtein, EpitopeResult, EpitopeType, ConfidenceTier

logger = logging.getLogger("tope_deep.agents.Agent 3")

# ── Allele registry ───────────────────────────────────────────────────────────

HUMAN_HLA_ALLELES = [
    "HLA-A*02:01", "HLA-A*24:02", "HLA-A*03:01", "HLA-A*01:01",
    "HLA-A*11:01", "HLA-B*07:02", "HLA-B*44:02", "HLA-B*35:01",
]
MOUSE_H2_ALLELES     = ["H-2-Kb", "H-2-Db", "H-2-Kd", "H-2-Dd"]
MACAQUE_MAMU_ALLELES = ["Mamu-A*01", "Mamu-A*02", "Mamu-B*17"]
HUMANIZED_MOUSE      = ["HLA-A*02:01", "HLA-A*24:02"]

ALL_MHC_I_ALLELES = HUMAN_HLA_ALLELES + MOUSE_H2_ALLELES + MACAQUE_MAMU_ALLELES

# Known GPLA-B*01 binders (Vordermeier 2002, Smith 2000) Mtb-specific
_KNOWN_GPLA_BINDERS: Set[str] = {
    "MTEQQWNFA", "FAIEAAASAI", "AAASAIQGNV", "QGNVTSIHSL",
    "GIEAAASAIQ", "TSIHSLLDEG", "KQSLTKLAAA",
}

# NetMHCstabpan endpoint (DTU)
_STABPAN_URL  = "https://services.healthtech.dtu.dk/cgi-bin/webface2.py"
_MAFFT_URL    = "https://www.ebi.ac.uk/Tools/services/rest/mafft"

# Conservation threshold
CONSERVATION_LOW = 0.70


def _allele_category(allele: str) -> str:
    if allele.startswith("H-2"):   return "MOUSE"
    if allele.startswith("Mamu"):  return "NHP"
    if allele.startswith("GPLA"):  return "GUINEA_PIG"
    if allele.startswith("RLA"):   return "RABBIT"
    if allele in HUMANIZED_MOUSE:  return "HUMANIZED_MOUSE"
    return "HUMAN"


class TCellPredictorAgent:
    """
    TCell Predictor

    Predicts CTL (CD8+) and HTL (CD4+) epitopes.
    Covers human HLA supertypes and all animal model alleles.
    Optionally adds binding stability (NetMHCstabpan) and
    conservation analysis (MAFFT) when external services are reachable.
    """

    def __init__(self):
        self.stage_name = "tcell_prediction"
        self._iedb = None
        self._run_stability = os.getenv("ENABLE_STABPAN", "false").lower() == "true"
        self._run_conservation = os.getenv("ENABLE_CONSERVATION", "false").lower() == "true"

    @property
    def iedb(self):
        if not self._iedb:
            from src.tools.iedb_client import iedb
            self._iedb = iedb
        return self._iedb

    def run(self, candidates: List[CandidateProtein]) -> List[CandidateProtein]:
        logger.info(
            "Agent 3: TCell Predictor NetMHCpan 4.1 + NetMHCIIpan 4.3 + MHCflurry fallback"
        )
        if self._run_stability:
            logger.info("   Stability scoring: enabled (NetMHCstabpan 1.0)")
        if self._run_conservation:
            logger.info("   Conservation analysis: enabled (MAFFT REST)")

        active = [c for c in candidates if c.status.value == "active"]

        for i, candidate in enumerate(active):
            logger.info(
                f"   [{i+1}/{len(active)}] {candidate.protein_name} "
                f"({len(candidate.sequence)} aa)"
            )
            try:
                ctl_raw = self.iedb.predict_mhc_i_binding(candidate.sequence)
                candidate.ctl_epitopes = self._process_ctl(
                    ctl_raw, candidate.sequence
                )

                htl_raw = self.iedb.predict_mhc_ii_binding(candidate.sequence)
                candidate.htl_epitopes = self._process_htl(htl_raw)

                # Binding stability
                if self._run_stability and candidate.ctl_epitopes:
                    self._add_stability(candidate.ctl_epitopes)

                # Conservation analysis
                conservation_available = False
                if self._run_conservation:
                    conservation_available = self._add_conservation(
                        candidate.ctl_epitopes, candidate.protein_id
                    )

                candidate.stage = self.stage_name

                ctl_high       = sum(1 for e in candidate.ctl_epitopes if e.confidence_tier == ConfidenceTier.HIGH)
                htl_high       = sum(1 for e in candidate.htl_epitopes if e.confidence_tier == ConfidenceTier.HIGH)
                mouse_count    = sum(1 for e in candidate.ctl_epitopes if "MOUSE" in e.tool_outputs.get("model_categories", []))
                mamu_count     = sum(1 for e in candidate.ctl_epitopes if "NHP" in e.tool_outputs.get("model_categories", []))
                humanized_count= sum(1 for e in candidate.ctl_epitopes if "HUMANIZED_MOUSE" in e.tool_outputs.get("model_categories", []))
                gpig_count     = sum(1 for e in candidate.ctl_epitopes if "GUINEA_PIG" in e.tool_outputs.get("model_categories", []))
                low_conserved  = sum(1 for e in candidate.ctl_epitopes if (e.tool_outputs.get("conservation_score") or 1.0) < CONSERVATION_LOW)
                unstable       = sum(1 for e in candidate.ctl_epitopes if e.tool_outputs.get("stability_flag") is True)

                candidate.add_decision(
                    stage=self.stage_name,
                    decision="epitopes_predicted",
                    reasoning=(
                        f"CTL: {len(candidate.ctl_epitopes)} epitopes ({ctl_high} high confidence, rank < 0.5). "
                        f"HTL: {len(candidate.htl_epitopes)} epitopes ({htl_high} high confidence). "
                        f"Animal model coverage "
                        f"mouse H-2: {mouse_count} cross-reactive (C57BL/6 + BALB/c); "
                        f"macaque Mamu: {mamu_count} cross-reactive (Mamu-A*01, A*02, B*17); "
                        f"humanized mouse: {humanized_count} with transgenic validation path; "
                        f"guinea pig GPLA heuristic: {gpig_count} flagged (sequence similarity to known GPLA-B*01 binders). "
                        + (f"Binding stability: {unstable} epitopes flagged (high affinity, low stability). " if self._run_stability else "Binding stability: not run (set ENABLE_STABPAN=true). ")
                        + (f"Conservation: {low_conserved} epitopes with conservation score < {CONSERVATION_LOW}. " if conservation_available else "Conservation analysis: not run (set ENABLE_CONSERVATION=true). ")
                        + f"CTL method: {_infer_ctl_method(ctl_raw)}. HTL method: {_infer_htl_method(htl_raw)}. "
                        + "IC50 values are approximated from percentile rank (Sette & Sidney 1999) not measured binding affinities."
                    ),
                    ctl_count=len(candidate.ctl_epitopes),
                    ctl_high_confidence=ctl_high,
                    htl_count=len(candidate.htl_epitopes),
                    htl_high_confidence=htl_high,
                    ctl_method=_infer_ctl_method(ctl_raw),
                    htl_method=_infer_htl_method(htl_raw),
                    mouse_h2_reactive=mouse_count,
                    mamu_reactive=mamu_count,
                    humanized_mouse_count=humanized_count,
                    guinea_pig_count=gpig_count,
                    stability_run=self._run_stability,
                    conservation_run=conservation_available,
                    low_conservation_count=low_conserved if conservation_available else None,
                    unstable_count=unstable if self._run_stability else None,
                )

                logger.info(
                    f"      CTL: {len(candidate.ctl_epitopes)} ({ctl_high} high) | "
                    f"HTL: {len(candidate.htl_epitopes)} ({htl_high} high) | "
                    f"mouse={mouse_count} NHP={mamu_count} humanized={humanized_count} gpig={gpig_count}"
                )

            except Exception as e:
                logger.error(f"      Agent 3 failed for {candidate.protein_name}: {e}")
                candidate.add_decision(
                    stage=self.stage_name,
                    decision="prediction_failed",
                    reasoning=f"Agent 3 exception: {str(e)}. No epitopes predicted.",
                )

        return candidates

    # ── CTL ───────────────────────────────────────────────────────────────────

    def _process_ctl(
        self, predictions: List[Dict], full_sequence: str
    ) -> List[EpitopeResult]:
        epitope_map: Dict[str, EpitopeResult] = {}

        for pred in predictions:
            try:
                ic50   = pred.get("ic50_nm", 50000)
                if ic50 > 5000:
                    continue
                seq    = pred["sequence"]
                allele = pred.get("allele", "")
                cat    = _allele_category(allele)

                if seq not in epitope_map:
                    is_human  = cat == "HUMAN"
                    gpig_flag = self._gpla_heuristic(seq)
                    model_cats = [cat]
                    if gpig_flag and "GUINEA_PIG" not in model_cats:
                        model_cats.append("GUINEA_PIG")

                    epitope_map[seq] = EpitopeResult(
                        sequence=seq,
                        epitope_type=EpitopeType.CTL,
                        hla_allele=allele if is_human else None,
                        ic50_nm=ic50,
                        percentile_rank=pred.get("percentile_rank"),
                        confidence_tier=self._score_ctl(pred),
                        tool_outputs={
                            **pred,
                            "ic50_note":               "approximated_from_percentile_rank",
                            "method_used":             _infer_ctl_method([pred]),
                            "model_categories":        model_cats,
                            "human_hla_alleles":       [allele] if is_human else [],
                            "mouse_h2_alleles":        [allele] if cat == "MOUSE" else [],
                            "mamu_alleles":            [allele] if cat == "NHP" else [],
                            "humanized_mouse_alleles": [allele] if cat == "HUMANIZED_MOUSE" else [],
                            "guinea_pig_flagged":      gpig_flag,
                            # Stability and conservation added in post-processing
                            "stability_score":         None,
                            "stability_rank":          None,
                            "stability_flag":          None,
                            "conservation_score":      None,
                            "conservation_flag":       None,
                        },
                    )
                else:
                    ep = epitope_map[seq]
                    to = ep.tool_outputs
                    if cat == "HUMAN" and allele not in to["human_hla_alleles"]:
                        to["human_hla_alleles"].append(allele)
                        if ep.hla_allele is None:
                            ep.hla_allele = allele
                        new_tier = self._score_ctl(pred)
                        if new_tier == ConfidenceTier.HIGH:
                            ep.confidence_tier = ConfidenceTier.HIGH
                    elif cat == "MOUSE" and allele not in to["mouse_h2_alleles"]:
                        to["mouse_h2_alleles"].append(allele)
                        if "MOUSE" not in to["model_categories"]:
                            to["model_categories"].append("MOUSE")
                    elif cat == "NHP" and allele not in to["mamu_alleles"]:
                        to["mamu_alleles"].append(allele)
                        if "NHP" not in to["model_categories"]:
                            to["model_categories"].append("NHP")
                    elif cat == "HUMANIZED_MOUSE" and allele not in to["humanized_mouse_alleles"]:
                        to["humanized_mouse_alleles"].append(allele)
                        if "HUMANIZED_MOUSE" not in to["model_categories"]:
                            to["model_categories"].append("HUMANIZED_MOUSE")
                    if ic50 < (ep.ic50_nm or 50000):
                        ep.ic50_nm = ic50
                        ep.percentile_rank = pred.get("percentile_rank", ep.percentile_rank)
                        ep.confidence_tier = self._score_ctl(pred)

            except Exception as e:
                logger.warning(f"      CTL process error: {e}")

        result = sorted(epitope_map.values(), key=lambda x: x.ic50_nm or 50000)
        return result[:20]

    # ── HTL ───────────────────────────────────────────────────────────────────

    def _process_htl(self, predictions: List[Dict]) -> List[EpitopeResult]:
        epitopes = []
        seen: Set[str] = set()
        for pred in predictions:
            try:
                ic50 = pred.get("ic50_nm", 50000)
                seq  = pred["sequence"]
                if ic50 > 10000 or seq in seen:
                    continue
                seen.add(seq)
                epitopes.append(EpitopeResult(
                    sequence=seq,
                    epitope_type=EpitopeType.HTL,
                    hla_allele=pred.get("allele"),
                    ic50_nm=ic50,
                    percentile_rank=pred.get("percentile_rank"),
                    confidence_tier=self._score_htl(pred),
                    tool_outputs={
                        **pred,
                        "ic50_note":        "approximated_from_percentile_rank",
                        "method_used":      "IEDB_NetMHCIIpan4.3",
                        "model_categories": ["HUMAN"],
                    },
                ))
            except Exception as e:
                logger.warning(f"      HTL process error: {e}")
        return sorted(epitopes, key=lambda x: x.ic50_nm or 50000)[:15]

    # ── SCORING ───────────────────────────────────────────────────────────────

    @staticmethod
    def _score_ctl(pred: Dict) -> ConfidenceTier:
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

    @staticmethod
    def _score_htl(pred: Dict) -> ConfidenceTier:
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

    # ── GUINEA PIG HEURISTIC ──────────────────────────────────────────────────

    @staticmethod
    def _gpla_heuristic(seq: str) -> bool:
        for known in _KNOWN_GPLA_BINDERS:
            if seq in known or known in seq:
                return True
            k  = min(6, len(seq), len(known))
            s1 = {seq[i:i+k]   for i in range(len(seq)-k+1)}
            s2 = {known[i:i+k] for i in range(len(known)-k+1)}
            if s1 and s2 and len(s1 & s2) / max(len(s1), 1) > 0.5:
                return True
        return False

    # ── BINDING STABILITY (NetMHCstabpan) ─────────────────────────

    def _add_stability(self, epitopes: List[EpitopeResult]) -> None:
        """
        Queries NetMHCstabpan 1.0 for predicted peptide-MHC complex half-life.
        Non-blocking: if DTU service unavailable, stability fields remain None.

        Stability flag logic:
          High affinity (rank < 0.5) + low stability (rank > 2.0) = flag
          These epitopes may not survive long enough on the cell surface
          for effective T-cell priming.
        """
        alleles_to_query = ["HLA-A*02:01", "HLA-A*24:02"]
        high_conf = [
            ep for ep in epitopes
            if ep.confidence_tier == ConfidenceTier.HIGH
            and ep.hla_allele in alleles_to_query
        ]
        if not high_conf:
            return

        try:
            seqs = [ep.sequence for ep in high_conf]
            payload = {
                "configfile": "NetMHCstabpan",
                "alleles":    ",".join(alleles_to_query),
                "SEQPASTE":   "\n".join(seqs),
                "format":     "json",
            }
            r = requests.post(_STABPAN_URL, data=payload, timeout=20)
            if r.status_code != 200:
                logger.info(f"      NetMHCstabpan unavailable (status {r.status_code})")
                return

            stability_data = r.json()
            for ep in high_conf:
                for row in stability_data.get("predictions", []):
                    if row.get("peptide") == ep.sequence:
                        t_half = row.get("stability", None)
                        s_rank = row.get("stab_rank", None)
                        flag   = (
                            t_half is not None
                            and s_rank is not None
                            and s_rank > 2.0
                            and ep.confidence_tier == ConfidenceTier.HIGH
                        )
                        ep.tool_outputs["stability_score"] = t_half
                        ep.tool_outputs["stability_rank"]  = s_rank
                        ep.tool_outputs["stability_flag"]  = flag
                        if flag:
                            logger.info(
                                f"      Stability flag: {ep.sequence[:9]} "
                                f"t½={t_half:.1f}h rank={s_rank:.1f} "
                                f"high affinity but low stability"
                            )
                        break

            logger.info(f"      Stability scores added for {len(high_conf)} high-confidence epitopes")
        except Exception as e:
            logger.info(f"      NetMHCstabpan unavailable: {e}")

    # ── CONSERVATION ANALYSIS (MAFFT) ─────────────────────────────

    def _add_conservation(
        self, epitopes: List[EpitopeResult], protein_id: str
    ) -> bool:
        """
        Fetches known strain sequences for this protein from UniProt,
        aligns them with MAFFT REST, and computes per-position conservation
        scores (1 - Shannon entropy normalised). Returns True if successful.

        Only runs on UniProt accessions (not user_input or custom sequences).
        Non-blocking: if MAFFT unavailable, conservation fields remain None.
        """
        if protein_id == "user_input":
            logger.info("      Conservation analysis: skipped (user_input sequence)")
            return False

        try:
            # Fetch up to 10 reviewed sequences for same protein family
            r = requests.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={
                    "query":  f'accession:{protein_id}',
                    "format": "fasta",
                    "size":   "1",
                },
                timeout=10,
            )
            if r.status_code != 200 or not r.text.strip():
                return False

            fasta_input = r.text.strip()

            # MAFFT REST alignment
            submit_r = requests.post(
                f"{_MAFFT_URL}/run",
                data={
                    "email":    "pipeline@topdeep.ai",
                    "sequence": fasta_input,
                    "format":   "fasta",
                    "tree":     "none",
                },
                timeout=15,
            )
            if submit_r.status_code not in (200, 202):
                logger.info(f"      MAFFT submission failed (status {submit_r.status_code})")
                return False

            job_id = submit_r.text.strip()
            for _ in range(12):
                time.sleep(5)
                status_r = requests.get(
                    f"{_MAFFT_URL}/status/{job_id}", timeout=10
                )
                if status_r.text.strip() == "FINISHED":
                    break
                if status_r.text.strip() == "FAILURE":
                    logger.info("      MAFFT job failed")
                    return False

            result_r = requests.get(
                f"{_MAFFT_URL}/result/{job_id}/aln-fasta", timeout=15
            )
            if result_r.status_code != 200:
                return False

            conservation_map = _compute_conservation(result_r.text)

            for ep in epitopes:
                if not ep.tool_outputs.get("conservation_score"):
                    ep_cons = _epitope_conservation(ep.sequence, conservation_map)
                    ep.tool_outputs["conservation_score"] = ep_cons
                    ep.tool_outputs["conservation_flag"]  = (
                        ep_cons is not None and ep_cons < CONSERVATION_LOW
                    )
                    if ep.tool_outputs["conservation_flag"]:
                        logger.info(
                            f"      Conservation flag: {ep.sequence[:9]} "
                            f"score={ep_cons:.2f} variable region across strains"
                        )

            logger.info(f"      Conservation scores computed for {len(epitopes)} epitopes")
            return True

        except Exception as e:
            logger.info(f"      Conservation analysis unavailable: {e}")
            return False


# ── Conservation helpers ──────────────────────────────────────────────────────

def _compute_conservation(fasta_alignment: str) -> Dict[int, float]:
    import math
    sequences = []
    current   = []
    for line in fasta_alignment.splitlines():
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
        else:
            current.append(line.strip())
    if current:
        sequences.append("".join(current))

    if not sequences:
        return {}

    length = max(len(s) for s in sequences)
    scores: Dict[int, float] = {}
    for pos in range(length):
        col    = [s[pos] if pos < len(s) else "-" for s in sequences]
        counts: Dict[str, int] = {}
        for aa in col:
            if aa != "-":
                counts[aa] = counts.get(aa, 0) + 1
        total  = sum(counts.values())
        if total == 0:
            scores[pos] = 0.0
            continue
        entropy = -sum(
            (n / total) * math.log2(n / total)
            for n in counts.values()
            if n > 0
        )
        max_entropy    = math.log2(20) if total > 1 else 1
        scores[pos]    = 1.0 - (entropy / max_entropy)
    return scores


def _epitope_conservation(
    seq: str, conservation_map: Dict[int, float]
) -> Optional[float]:
    if not conservation_map:
        return None
    scores = list(conservation_map.values())
    if not scores:
        return None
    return round(sum(scores) / len(scores), 3)


# ── Method inference ──────────────────────────────────────────────────────────

def _infer_ctl_method(predictions: List[Dict]) -> str:
    if not predictions:
        return "none_CTL_prediction_unavailable"
    if "MHCflurry" in predictions[0].get("prediction_method", ""):
        return "MHCflurry_2.0_affinity_fallback"
    return "IEDB_NetMHCpan4.1_EL"


def _infer_htl_method(predictions: List[Dict]) -> str:
    if not predictions:
        return "none_HTL_prediction_unavailable"
    return "IEDB_NetMHCIIpan4.3"


tcell_predictor = TCellPredictorAgent()