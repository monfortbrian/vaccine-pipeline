"""
STRUCTURE AGENT - TOPE_DEEP NODE N5
Version 2.0 - Fixed pLDDT field name resolution.

Change from v1: AlphaFold DB API field name for mean pLDDT has changed
across API versions. Now tries multiple field names in order:
  meanPlddt, pLDDT, globalMetricValue, confidenceScore

Also computes mean from pLDDT array if scalar field is 0 or missing.
"""

import logging
import requests
from typing import List, Optional, Dict, Any
from src.models.candidate import CandidateProtein

logger = logging.getLogger(__name__)

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
REQUEST_TIMEOUT = 20


class StructureAgent:
    def __init__(self):
        self.stage_name = "structure_retrieval"

    def run(self, candidates: List[CandidateProtein]) -> List[CandidateProtein]:
        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"N5 StructureAgent: {len(active)} active candidates")

        for i, candidate in enumerate(active):
            logger.info(
                f"  [{i+1}/{len(active)}] {candidate.protein_name} "
                f"({candidate.protein_id})"
            )

            if not self._is_uniprot_id(candidate.protein_id):
                candidate.structure_source = "unavailable"
                candidate.add_decision(
                    stage=self.stage_name,
                    decision="skipped",
                    reasoning=(
                        f"'{candidate.protein_id}' is not a UniProt accession. "
                        f"AlphaFold DB requires UniProt ID. "
                        f"Submit to ColabFold manually."
                    ),
                    structure_source="unavailable",
                )
                continue

            result = self._fetch_alphafold(candidate.protein_id)

            if result:
                candidate.structure_source = "alphafold_db"
                candidate.structure_pdb_path = result["cif_url"]

                plddt = result["mean_plddt"]
                if plddt >= 90:
                    plddt_label = "very high confidence (≥90)"
                elif plddt >= 70:
                    plddt_label = "confident (70–89)"
                elif plddt >= 50:
                    plddt_label = "low confidence (50–69)"
                else:
                    plddt_label = "very low confidence (<50) - backbone unreliable"

                candidate.add_decision(
                    stage=self.stage_name,
                    decision="structure_retrieved",
                    reasoning=(
                        f"AlphaFold DB entry found. "
                        f"Model version: {result['model_version']}. "
                        f"Mean pLDDT: {plddt:.1f}/100 ({plddt_label}). "
                        f"Coverage: {result['fragment_coverage']}. "
                        f"pLDDT interpretation: Jumper et al. (2021) "
                        f"doi:10.1038/s41586-021-03819-2."
                    ),
                    structure_source="alphafold_db",
                    alphafold_entry_id=result["entry_id"],
                    model_version=result["model_version"],
                    mean_plddt=plddt,
                    plddt_field_used=result["plddt_field"],
                    pdb_url=result["pdb_url"],
                    cif_url=result["cif_url"],
                    fragment_coverage=result["fragment_coverage"],
                )
                logger.info(
                    f"    AlphaFold: pLDDT={plddt:.1f} "
                    f"[field={result['plddt_field']}] "
                    f"version={result['model_version']}"
                )
            else:
                candidate.structure_source = "unavailable"
                candidate.structure_pdb_path = None
                candidate.add_decision(
                    stage=self.stage_name,
                    decision="structure_unavailable",
                    reasoning=(
                        f"No AlphaFold DB entry for '{candidate.protein_id}'. "
                        f"B-cell conformational predictions use linear sequence only "
                        f"(reduced accuracy). "
                        f"Submit to ColabFold for de novo structure prediction."
                    ),
                    structure_source="unavailable",
                    colabfold_hint=(
                        "https://colab.research.google.com/github/sokrypton/ColabFold/"
                        "blob/main/AlphaFold2.ipynb"
                    ),
                )

        total = sum(1 for c in candidates if c.structure_source == "alphafold_db")
        logger.info(f"N5 complete: {total}/{len(active)} have AlphaFold structures")
        return candidates

    def _fetch_alphafold(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        url = f"{ALPHAFOLD_API}/{uniprot_id}"
        try:
            resp = requests.get(
                url,
                headers={"Accept": "application/json"},
                timeout=REQUEST_TIMEOUT,
            )
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            data = resp.json()
            if not data:
                return None

            entry = data[0]

            # ── pLDDT field resolution - tries multiple field names ──────────
            # AlphaFold DB has changed field names across API versions.
            # Tries each in priority order; falls back to computing from array.
            mean_plddt = 0.0
            plddt_field = "not_found"

            for field in ["meanPlddt", "pLDDT", "globalMetricValue", "confidenceScore"]:
                val = entry.get(field)
                if val and float(val) > 0:
                    mean_plddt = round(float(val), 2)
                    plddt_field = field
                    break

            # If all scalar fields are 0/missing, try computing from pLDDT array
            if mean_plddt == 0.0:
                plddt_array = entry.get("plddt") or entry.get("confidenceValues")
                if plddt_array and isinstance(plddt_array, list) and len(plddt_array) > 0:
                    mean_plddt = round(sum(plddt_array) / len(plddt_array), 2)
                    plddt_field = "computed_from_array"

            if mean_plddt == 0.0:
                logger.warning(
                    f"pLDDT could not be resolved for {uniprot_id}. "
                    f"Available fields: {list(entry.keys())}"
                )
                plddt_field = "unavailable"

            return {
                "entry_id":          entry.get("entryId", uniprot_id),
                "model_version":     entry.get("latestVersion", "unknown"),
                "pdb_url":           entry.get("pdbUrl", ""),
                "cif_url":           entry.get("cifUrl", entry.get("pdbUrl", "")),
                "mean_plddt":        mean_plddt,
                "plddt_field":       plddt_field,
                "fragment_coverage": (
                    f"{entry.get('uniprotStart', '?')}-{entry.get('uniprotEnd', '?')}"
                ),
            }

        except requests.exceptions.Timeout:
            logger.warning(f"AlphaFold DB timeout for {uniprot_id}")
            return None
        except Exception as e:
            logger.warning(f"AlphaFold DB error for {uniprot_id}: {e}")
            return None

    @staticmethod
    def _is_uniprot_id(protein_id: str) -> bool:
        import re
        pattern = (
            r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$"
            r"|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
        )
        return bool(re.match(pattern, protein_id.strip().upper()))


structure_agent = StructureAgent()