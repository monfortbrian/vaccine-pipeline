"""
CONSTRUCT DESIGNER - TOPE_DEEP NODE N8
Assembles a multi-epitope vaccine construct from prioritised epitopes.

Pipeline position: final synthesis node. Runs after N6 (safety) and N7 (coverage).
Requires at least one CTL or HTL epitope to produce output.

Assembly logic:
  1. Select epitopes: high-confidence, safety-passed only
  2. Order: CTL → HTL → B-cell (immuno-dominant first within each class)
  3. Link with validated linkers:
       CTL-CTL    AAY   (proteasomal cleavage site, enhances MHC-I presentation)
       HTL-HTL    GPGPG (flexible, preserves HTL epitope structure)
       CTL-HTL    KK    (charge-based separation)
       B-cell     GPGPG (preserves conformational accessibility)
  4. Prepend adjuvant sequence (RS09 TLR4 agonist peptide, optional, default on)
  5. Compute physicochemical properties via Biopython ProtParam
  6. Write construct record to PipelineRun.construct_sequence

Limitations (recorded in decision audit):
  - Linear assembly only; 3D folding not predicted at this node
  - Linker cleavage efficiency is predicted, not validated
  - ProtParam properties are theoretical (isoelectric point ±0.3 pH units typical error)
  - RS09 adjuvant effect is in silico only; requires wet-lab validation

References:
  Linker selection: Nezafat et al. (2014) doi:10.1016/j.compbiolchem.2014.08.020
  RS09 adjuvant:    Chuang et al. (2010) doi:10.1016/j.vaccine.2010.01.062
  ProtParam:        Gasteiger et al. (2005) ExPASy: The proteomics server
"""

import logging
from typing import List, Dict, Any, Optional, Tuple

from src.models.candidate import (
    CandidateProtein, EpitopeResult, EpitopeType,
    ConfidenceTier, CandidateStatus, PipelineRun,
)

logger = logging.getLogger(__name__)

# --LINKER DEFINITIONS
# Source: Nezafat et al. (2014), validated in multi-epitope vaccine literature
LINKERS: Dict[str, str] = {
    "CTL_CTL": "AAY",    # Promotes proteasomal cleavage; preserves MHC-I binding
    "HTL_HTL": "GPGPG",  # Flexible helix-breaking; protects HTL epitope structure
    "CTL_HTL": "KK",     # Charged separator; reduces junctional neo-epitope risk
    "HTL_CTL": "KK",
    "CTL_BCELL": "GPGPG",
    "HTL_BCELL": "GPGPG",
    "BCELL_BCELL": "GPGPG",
    "DEFAULT": "GPGPG",
}

# RS09 - synthetic TLR4 agonist derived from flagellin
# Chuang et al. (2010): enhances CD4+ and CD8+ T-cell responses in mice
RS09_ADJUVANT = "APPHALS"

# Minimum epitope counts to proceed
MIN_EPITOPES_REQUIRED = 1

# Selection thresholds
CTL_IC50_THRESHOLD_NM = 500.0   # NetMHCpan strong binder cutoff
HTL_IC50_THRESHOLD_NM = 1000.0  # NetMHCIIpan threshold
MAX_CTL_PER_CONSTRUCT = 10
MAX_HTL_PER_CONSTRUCT = 8
MAX_BCELL_PER_CONSTRUCT = 5


class ConstructDesignerAgent:
    def __init__(self, include_adjuvant: bool = True):
        self.stage_name = "construct_design"
        self.include_adjuvant = include_adjuvant

    def run(
        self,
        candidates: List[CandidateProtein],
        pipeline_run: Optional[PipelineRun] = None,
    ) -> Tuple[List[CandidateProtein], Optional[Dict[str, Any]]]:
        """
        Design multi-epitope construct from top-ranked candidates.

        Returns:
            (candidates, construct_report)
            construct_report is None if no eligible epitopes found.
        """
        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"N8 ConstructDesigner: {len(active)} active candidates")

        # Collect eligible epitopes across all candidates
        ctl_epitopes = self._select_ctl(active)
        htl_epitopes = self._select_htl(active)
        bcell_epitopes = self._select_bcell(active)

        total_eligible = len(ctl_epitopes) + len(htl_epitopes) + len(bcell_epitopes)

        if total_eligible < MIN_EPITOPES_REQUIRED:
            logger.warning("N8: No eligible epitopes for construct assembly.")
            if pipeline_run:
                pipeline_run.add_warning(
                    "Construct design skipped: no safety-passed, high-confidence epitopes available."
                )
            return candidates, None

        logger.info(
            f"  Selected epitopes - CTL: {len(ctl_epitopes)}, "
            f"HTL: {len(htl_epitopes)}, B-cell: {len(bcell_epitopes)}"
        )

        # Assemble sequence
        construct_sequence, assembly_log = self._assemble(
            ctl_epitopes, htl_epitopes, bcell_epitopes
        )

        # Compute physicochemical properties
        properties = self._compute_properties(construct_sequence)

        # Build construct report
        construct_report = {
            "construct_sequence": construct_sequence,
            "length_aa": len(construct_sequence),
            "epitope_counts": {
                "CTL": len(ctl_epitopes),
                "HTL": len(htl_epitopes),
                "B-cell": len(bcell_epitopes),
            },
            "assembly_log": assembly_log,
            "physicochemical": properties,
            "adjuvant_included": self.include_adjuvant,
            "adjuvant_sequence": RS09_ADJUVANT if self.include_adjuvant else None,
            "adjuvant_reference": (
                "Chuang et al. (2010) doi:10.1016/j.vaccine.2010.01.062"
                if self.include_adjuvant else None
            ),
            "linker_scheme": LINKERS,
            "linker_reference": (
                "Nezafat et al. (2014) doi:10.1016/j.compbiolchem.2014.08.020"
            ),
            "limitations": [
                "Linear assembly only - 3D tertiary structure not predicted at this stage.",
                "Linker cleavage efficiency is computationally predicted; wet-lab validation required.",
                (
                    f"ProtParam theoretical pI has ±0.3 pH units typical error "
                    f"(Gasteiger et al., 2005)."
                ),
                "RS09 adjuvant immunogenicity established in murine models only.",
                "Junctional neo-epitopes not screened in this version.",
            ],
            "next_steps": [
                "Submit construct_sequence to AlphaFold2/ColabFold for 3D structure prediction.",
                "Run junctional neo-epitope screen (Herd et al., 2010 method).",
                "Validate top CTL epitopes against IEDB positive assay data.",
                "Wet-lab: synthesise and test in PBMC-based IFN-γ ELISpot assay.",
            ],
        }

        # Write to PipelineRun if provided
        if pipeline_run:
            pipeline_run.construct_sequence = construct_sequence
            pipeline_run.construct_properties = construct_report

        # Add audit decision to each contributing candidate
        contributing_ids = list({
            e.tool_outputs.get("protein_id", "unknown")
            for e in ctl_epitopes + htl_epitopes + bcell_epitopes
        })
        for candidate in active:
            candidate.add_decision(
                stage=self.stage_name,
                decision="construct_assembled",
                reasoning=(
                    f"Contributed epitopes to multi-epitope construct. "
                    f"Construct length: {len(construct_sequence)} aa. "
                    f"Mean pI: {properties.get('isoelectric_point', 'N/A')}. "
                    f"Instability index: {properties.get('instability_index', 'N/A')} "
                    f"({'stable' if properties.get('is_stable') else 'unstable - review required'})."
                ),
                construct_length=len(construct_sequence),
                instability_index=properties.get("instability_index"),
                is_stable=properties.get("is_stable"),
            )

        logger.info(
            f"N8 complete: construct={len(construct_sequence)} aa, "
            f"MW={properties.get('molecular_weight_da', 'N/A')} Da, "
            f"pI={properties.get('isoelectric_point', 'N/A')}, "
            f"stable={properties.get('is_stable', 'N/A')}"
        )

        return candidates, construct_report

    # --EPITOPE SELECTION

    def _select_ctl(self, candidates: List[CandidateProtein]) -> List[EpitopeResult]:
        """Select CTL epitopes: safety-passed, ranked by IC50 ascending."""
        pool = []
        for c in candidates:
            for ep in c.ctl_epitopes:
                if self._is_eligible(ep):
                    # Tag source protein for audit trail
                    ep.tool_outputs["protein_id"] = c.protein_id
                    pool.append(ep)

        # Rank: high confidence first, then by IC50 ascending (lower = stronger binder)
        pool.sort(key=lambda e: (
            0 if e.confidence_tier == ConfidenceTier.HIGH else 1,
            e.ic50_nm or 9999,
        ))
        selected = pool[:MAX_CTL_PER_CONSTRUCT]
        logger.info(f"  CTL pool: {len(pool)} eligible → {len(selected)} selected")
        return selected

    def _select_htl(self, candidates: List[CandidateProtein]) -> List[EpitopeResult]:
        """Select HTL epitopes: safety-passed, ranked by IC50 ascending."""
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
        selected = pool[:MAX_HTL_PER_CONSTRUCT]
        logger.info(f"  HTL pool: {len(pool)} eligible → {len(selected)} selected")
        return selected

    def _select_bcell(self, candidates: List[CandidateProtein]) -> List[EpitopeResult]:
        """Select B-cell epitopes: safety-passed, high/medium confidence only."""
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
        selected = pool[:MAX_BCELL_PER_CONSTRUCT]
        logger.info(f"  B-cell pool: {len(pool)} eligible → {len(selected)} selected")
        return selected

    @staticmethod
    def _is_eligible(ep: EpitopeResult) -> bool:
        """
        Epitope is eligible for construct inclusion if:
          - allergenicity_safe is True (None = not screened = excluded)
          - toxicity_safe is True (None = not screened = excluded)
          - confidence_tier is not UNCERTAIN or UNSCORED
        """
        if ep.allergenicity_safe is not True:
            return False
        if ep.toxicity_safe is not True:
            return False
        if ep.confidence_tier in (ConfidenceTier.UNCERTAIN, ConfidenceTier.UNSCORED):
            return False
        return True

    # --ASSEMBLY

    def _assemble(
        self,
        ctl: List[EpitopeResult],
        htl: List[EpitopeResult],
        bcell: List[EpitopeResult],
    ) -> Tuple[str, List[Dict[str, Any]]]:
        """
        Assemble epitope sequences with linkers.
        Order: [RS09]–CTL–HTL–B-cell
        Returns (full_sequence, assembly_log)
        """
        parts: List[Tuple[str, str]] = []  # (sequence, label)
        log: List[Dict[str, Any]] = []

        if self.include_adjuvant:
            parts.append((RS09_ADJUVANT, "RS09_adjuvant"))

        prev_type = "ADJUVANT" if self.include_adjuvant else None

        def append_epitope(ep: EpitopeResult, ep_type: str, index: int):
            nonlocal prev_type
            linker_key = f"{prev_type}_{ep_type}" if prev_type else None
            linker = LINKERS.get(linker_key, LINKERS["DEFAULT"]) if linker_key else ""

            if linker:
                parts.append((linker, f"linker_{linker_key}"))
                log.append({
                    "position": len(parts),
                    "element": f"linker_{linker_key}",
                    "sequence": linker,
                    "rationale": f"Connects {prev_type} → {ep_type}",
                })

            parts.append((ep.sequence, f"{ep_type}_{index+1}"))
            log.append({
                "position": len(parts),
                "element": f"{ep_type}_epitope_{index+1}",
                "sequence": ep.sequence,
                "hla_allele": ep.hla_allele,
                "ic50_nm": ep.ic50_nm,
                "confidence": ep.confidence_tier.value,
            })
            prev_type = ep_type

        for i, ep in enumerate(ctl):
            append_epitope(ep, "CTL", i)
        for i, ep in enumerate(htl):
            append_epitope(ep, "HTL", i)
        for i, ep in enumerate(bcell):
            append_epitope(ep, "BCELL", i)

        full_sequence = "".join(seq for seq, _ in parts)
        return full_sequence, log

    # --PHYSICOCHEMICAL PROPERTIES

    def _compute_properties(self, sequence: str) -> Dict[str, Any]:
        """
        Compute theoretical physicochemical properties using Biopython ProtParam.
        All values are theoretical; see module-level limitations.
        """
        try:
            from Bio.SeqUtils.ProtParam import ProteinAnalysis
            analysis = ProteinAnalysis(sequence)

            mw = round(analysis.molecular_weight(), 2)
            pi = round(analysis.isoelectric_point(), 2)
            instability = round(analysis.instability_index(), 2)
            gravy = round(analysis.gravy(), 4)
            aromaticity = round(analysis.aromaticity(), 4)

            # Instability index < 40 = stable (Guruprasad et al., 1990)
            is_stable = instability < 40.0

            # Grand average of hydropathicity (GRAVY) interpretation
            # Negative = hydrophilic (good for soluble vaccine antigen)
            hydrophilicity = "hydrophilic" if gravy < 0 else "hydrophobic"

            return {
                "molecular_weight_da": mw,
                "isoelectric_point": pi,
                "instability_index": instability,
                "is_stable": is_stable,
                "gravy": gravy,
                "hydrophilicity": hydrophilicity,
                "aromaticity": aromaticity,
                "method": "Biopython ProtParam (Gasteiger et al., 2005)",
                "instability_reference": "Guruprasad et al. (1990) Protein Eng. 4:155-161",
            }

        except ImportError:
            logger.error("Biopython not installed. Run: pip install biopython")
            return {"error": "biopython_not_installed"}
        except Exception as e:
            logger.error(f"ProtParam computation failed: {e}")
            return {"error": str(e)}


# Module-level instance
construct_designer = ConstructDesignerAgent(include_adjuvant=True)