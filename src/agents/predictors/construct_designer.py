"""
TOPE_DEEP NODE N8
Assembles a multi-epitope vaccine construct from prioritised epitopes.

Adjuvant options (configurable via agent constructor):
  RS09   APPHALS       TLR4 agonist (Chuang et al. 2010)
                       Validated: murine models only.
  PADRE  aKXVAAWTLKAAa Pan-DR MHC-II binder (Alexander et al. 1994)
                       Validated: human clinical trials (HIV, HCV).
  TpD    QYIKANSKFIGITELKKLESKINKVF
                       Tetanus p30 universal helper epitope (Valmori et al. 1992)
                       Validated: human. Most conservative choice.
  none   (no adjuvant) Raw epitope string only.

Default: RS09 (standard in computational multi-epitope vaccine literature).
For human-facing reports: recommend PADRE or TpD; both have human data.

Linkers (Nezafat et al. 2014 doi:10.1016/j.compbiolchem.2014.08.020):
  CTL-CTL   AAY    proteasomal cleavage, preserves MHC-I presentation
  HTL-HTL   GPGPG  flexible helix-breaking, protects HTL structure
  CTL-HTL   KK     charged separator, reduces junctional neo-epitope risk
  B-cell    GPGPG  preserves conformational accessibility

Physicochemical properties: Biopython ProtParam (Gasteiger et al. 2005).
Instability index threshold: <40 = stable (Guruprasad et al. 1990).
"""

import logging
from typing import List, Dict, Any, Optional, Tuple
from src.models.candidate import (
    CandidateProtein, EpitopeResult, EpitopeType, ConfidenceTier,
    CandidateStatus, PipelineRun,
)

logger = logging.getLogger(__name__)

# ── Adjuvant registry ─────────────────────────────────────────────────────────
# Add new adjuvants here. Each entry is independently citable and defensible.
ADJUVANTS: Dict[str, Dict] = {
    "RS09": {
        "sequence": "APPHALS",
        "mechanism": "TLR4 agonist derived from flagellin C-terminal domain",
        "validation": "murine models (in vitro + in vivo)",
        "citation": "Chuang et al. (2010) Vaccine doi:10.1016/j.vaccine.2010.01.062",
        "note": "Immunogenicity established in murine models only. Human translation requires validation.",
    },
    "PADRE": {
        "sequence": "AKFVAAWTLKAAA",
        "mechanism": "Pan-DR epitope; promiscuous MHC-II binder, enhances CD4+ T-cell help",
        "validation": "human clinical trials (HIV, HCV peptide vaccines)",
        "citation": "Alexander et al. (1994) Immunity 1(7):751-761",
        "note": "Most appropriate for human-facing publications. Contains Aib (non-natural AA in original; Ala used here for synthesis compatibility).",
    },
    "TpD": {
        "sequence": "QYIKANSKFIGITELKKLESKINKVF",
        "mechanism": "Tetanus toxin p30 universal helper T-cell epitope",
        "validation": "human  widely used in clinical peptide vaccine trials",
        "citation": "Valmori et al. (1992) J Immunol 149(2):717-721",
        "note": "Most conservative and clinically validated option. Pre-existing tetanus immunity in vaccinated populations may affect response.",
    },
    "none": {
        "sequence": "",
        "mechanism": "No adjuvant, raw epitope string only",
        "validation": "N/A",
        "citation": "N/A",
        "note": "Use when adjuvant selection is deferred to wet-lab team.",
    },
}

DEFAULT_ADJUVANT = "RS09"

# ── Linker definitions ────────────────────────────────────────────────────────
LINKERS: Dict[str, str] = {
    "CTL_CTL":    "AAY",
    "HTL_HTL":    "GPGPG",
    "CTL_HTL":    "KK",
    "HTL_CTL":    "KK",
    "CTL_BCELL":  "GPGPG",
    "HTL_BCELL":  "GPGPG",
    "BCELL_BCELL":"GPGPG",
    "DEFAULT":    "GPGPG",
}

LINKER_CITATION = "Nezafat et al. (2014) Comput Biol Chem doi:10.1016/j.compbiolchem.2014.08.020"

# ── Selection thresholds ──────────────────────────────────────────────────────
CTL_IC50_THRESHOLD_NM  = 500.0
HTL_IC50_THRESHOLD_NM  = 1000.0
MAX_CTL_PER_CONSTRUCT  = 10
MAX_HTL_PER_CONSTRUCT  = 8
MAX_BCELL_PER_CONSTRUCT = 5
MIN_EPITOPES_REQUIRED  = 1


class ConstructDesignerAgent:
    def __init__(self, adjuvant: str = DEFAULT_ADJUVANT):
        self.stage_name = "construct_design"
        if adjuvant not in ADJUVANTS:
            logger.warning(f"Unknown adjuvant '{adjuvant}', using {DEFAULT_ADJUVANT}")
            adjuvant = DEFAULT_ADJUVANT
        self.adjuvant_key = adjuvant
        self.adjuvant = ADJUVANTS[adjuvant]

    def run(
        self,
        candidates: List[CandidateProtein],
        pipeline_run: Optional[PipelineRun] = None,
    ) -> Tuple[List[CandidateProtein], Optional[Dict[str, Any]]]:
        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"N8: {len(active)} candidates | adjuvant={self.adjuvant_key}")

        ctl_epitopes   = self._select_ctl(active)
        htl_epitopes   = self._select_htl(active)
        bcell_epitopes = self._select_bcell(active)

        total = len(ctl_epitopes) + len(htl_epitopes) + len(bcell_epitopes)

        if total < MIN_EPITOPES_REQUIRED:
            logger.warning(
                "N8: No eligible epitopes (safety-passed + high-confidence). "
                "All epitopes may be unscored (N6 tools unavailable) or failed safety. "
                "Construct not assembled."
            )
            for candidate in active:
                candidate.add_decision(
                    stage=self.stage_name,
                    decision="construct_skipped",
                    reasoning=(
                        "No safety-passed, confidence-scored epitopes available. "
                        "This occurs when N6 safety tools are unavailable (all unscored) "
                        "or all epitopes failed safety screening. "
                        "Review N6 safety screening results."
                    ),
                )
            return candidates, None

        logger.info(f"  Selected: CTL={len(ctl_epitopes)} HTL={len(htl_epitopes)} B-cell={len(bcell_epitopes)}")

        construct_sequence, assembly_log = self._assemble(ctl_epitopes, htl_epitopes, bcell_epitopes)
        properties = self._compute_properties(construct_sequence)

        construct_report = {
            "construct_sequence": construct_sequence,
            "length_aa": len(construct_sequence),
            "epitope_counts": {
                "CTL": len(ctl_epitopes),
                "HTL": len(htl_epitopes),
                "B-cell": len(bcell_epitopes),
            },
            "adjuvant": {
                "key": self.adjuvant_key,
                "sequence": self.adjuvant["sequence"],
                "mechanism": self.adjuvant["mechanism"],
                "validation": self.adjuvant["validation"],
                "citation": self.adjuvant["citation"],
                "note": self.adjuvant["note"],
            },
            "linker_scheme": LINKERS,
            "linker_citation": LINKER_CITATION,
            "assembly_log": assembly_log,
            "physicochemical": properties,
            "limitations": [
                "Linear assembly only 3D tertiary structure not predicted at this stage.",
                "Linker cleavage efficiency is computationally predicted; wet-lab validation required.",
                f"ProtParam theoretical pI has ±0.3 pH units typical error (Gasteiger et al. 2005).",
                f"Adjuvant ({self.adjuvant_key}): {self.adjuvant['note']}",
                "Junctional neo-epitopes not screened in this version.",
            ],
            "next_steps": [
                "Submit construct_sequence to AlphaFold2/ColabFold for 3D structure prediction.",
                "Run junctional neo-epitope screen (Herd et al. 2010 method).",
                "Validate top CTL epitopes against IEDB positive assay data.",
                "Wet-lab: synthesise and test in PBMC-based IFN-γ ELISpot assay.",
            ],
        }

        if pipeline_run:
            pipeline_run.construct_sequence = construct_sequence
            pipeline_run.construct_properties = construct_report

        for candidate in active:
            candidate.add_decision(
                stage=self.stage_name,
                decision="construct_assembled",
                reasoning=(
                    f"Construct assembled: {len(construct_sequence)} aa. "
                    f"CTL={len(ctl_epitopes)}, HTL={len(htl_epitopes)}, "
                    f"B-cell={len(bcell_epitopes)}. "
                    f"Adjuvant: {self.adjuvant_key} ({self.adjuvant['mechanism']}). "
                    f"pI={properties.get('isoelectric_point')}, "
                    f"instability={properties.get('instability_index')} "
                    f"({'stable' if properties.get('is_stable') else 'unstable review'})."
                ),
                construct_length=len(construct_sequence),
                adjuvant_used=self.adjuvant_key,
                instability_index=properties.get("instability_index"),
                is_stable=properties.get("is_stable"),
            )

        logger.info(
            f"N8 complete: {len(construct_sequence)} aa | "
            f"MW={properties.get('molecular_weight_da')} Da | "
            f"pI={properties.get('isoelectric_point')} | "
            f"stable={properties.get('is_stable')}"
        )
        return candidates, construct_report

    # ── SELECTION ─────────────────────────────────────────────────────────────

    def _select_ctl(self, candidates: List[CandidateProtein]) -> List[EpitopeResult]:
        pool = []
        for c in candidates:
            for ep in c.ctl_epitopes:
                if self._is_eligible(ep):
                    ep.tool_outputs["protein_id"] = c.protein_id
                    pool.append(ep)
        pool.sort(key=lambda e: (
            0 if e.confidence_tier == ConfidenceTier.HIGH else 1,
            e.ic50_nm or 9999,
        ))
        return pool[:MAX_CTL_PER_CONSTRUCT]

    def _select_htl(self, candidates: List[CandidateProtein]) -> List[EpitopeResult]:
        pool = []
        for c in candidates:
            for ep in c.htl_epitopes:
                if self._is_eligible(ep):
                    ep.tool_outputs["protein_id"] = c.protein_id
                    pool.append(ep)
        pool.sort(key=lambda e: (
            0 if e.confidence_tier == ConfidenceTier.HIGH else 1,
            e.ic50_nm or 9999,
        ))
        return pool[:MAX_HTL_PER_CONSTRUCT]

    def _select_bcell(self, candidates: List[CandidateProtein]) -> List[EpitopeResult]:
        pool = []
        for c in candidates:
            for ep in c.bcell_epitopes:
                if (
                    self._is_eligible(ep)
                    and ep.confidence_tier in (ConfidenceTier.HIGH, ConfidenceTier.MEDIUM)
                ):
                    ep.tool_outputs["protein_id"] = c.protein_id
                    pool.append(ep)
        pool.sort(key=lambda e: 0 if e.confidence_tier == ConfidenceTier.HIGH else 1)
        return pool[:MAX_BCELL_PER_CONSTRUCT]

    @staticmethod
    def _is_eligible(ep: EpitopeResult) -> bool:
        if ep.allergenicity_safe is not True:
            return False
        if ep.toxicity_safe is not True:
            return False
        if ep.confidence_tier in (ConfidenceTier.UNCERTAIN,):
            return False
        return True

    # ── ASSEMBLY ──────────────────────────────────────────────────────────────

    def _assemble(
        self,
        ctl: List[EpitopeResult],
        htl: List[EpitopeResult],
        bcell: List[EpitopeResult],
    ) -> Tuple[str, List[Dict[str, Any]]]:
        parts: List[Tuple[str, str]] = []
        log: List[Dict[str, Any]] = []

        adj_seq = self.adjuvant["sequence"]
        if adj_seq:
            parts.append((adj_seq, f"{self.adjuvant_key}_adjuvant"))

        prev_type = self.adjuvant_key.upper() if adj_seq else None

        def append_epitope(ep: EpitopeResult, ep_type: str, index: int):
            nonlocal prev_type
            linker = ""
            if prev_type:
                linker_key = f"{prev_type}_{ep_type}"
                linker = LINKERS.get(linker_key, LINKERS["DEFAULT"])
                parts.append((linker, f"linker_{linker_key}"))
                log.append({
                    "element": f"linker_{linker_key}",
                    "sequence": linker,
                    "rationale": f"{prev_type} → {ep_type}",
                })
            parts.append((ep.sequence, f"{ep_type}_{index+1}"))
            log.append({
                "element": f"{ep_type}_epitope_{index+1}",
                "sequence": ep.sequence,
                "hla_allele": ep.hla_allele,
                "ic50_nm": ep.ic50_nm,
                "confidence": ep.confidence_tier.value,
            })
            prev_type = ep_type

        for i, ep in enumerate(ctl):   append_epitope(ep, "CTL",   i)
        for i, ep in enumerate(htl):   append_epitope(ep, "HTL",   i)
        for i, ep in enumerate(bcell): append_epitope(ep, "BCELL", i)

        return "".join(seq for seq, _ in parts), log

    # ── PHYSICOCHEMICAL ───────────────────────────────────────────────────────

    def _compute_properties(self, sequence: str) -> Dict[str, Any]:
        try:
            from Bio.SeqUtils.ProtParam import ProteinAnalysis
            a = ProteinAnalysis(sequence)
            instability = round(a.instability_index(), 2)
            gravy = round(a.gravy(), 4)
            return {
                "molecular_weight_da": round(a.molecular_weight(), 2),
                "isoelectric_point":   round(a.isoelectric_point(), 2),
                "instability_index":   instability,
                "is_stable":           instability < 40.0,
                "gravy":               gravy,
                "hydrophilicity":      "hydrophilic" if gravy < 0 else "hydrophobic",
                "aromaticity":         round(a.aromaticity(), 4),
                "method":              "Biopython ProtParam (Gasteiger et al. 2005)",
                "instability_ref":     "Guruprasad et al. (1990) Protein Eng 4:155-161",
            }
        except ImportError:
            return {"error": "biopython_not_installed"}
        except Exception as e:
            return {"error": str(e)}


def get_available_adjuvants() -> Dict[str, Dict]:
    """Exposed for /api/health and UI adjuvant selector."""
    return {
        k: {
            "mechanism":  v["mechanism"],
            "validation": v["validation"],
            "citation":   v["citation"],
            "note":       v["note"],
        }
        for k, v in ADJUVANTS.items()
    }


construct_designer = ConstructDesignerAgent(adjuvant=DEFAULT_ADJUVANT)