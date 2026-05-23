"""
STRUCTURE AGENT - TOPE_DEEP NODE N5
Retrieves 3D protein structure metadata from AlphaFold DB using UniProt accession.

Scope:
  - Fetches AlphaFold DB entry for each candidate by UniProt ID
  - Stores structure URL, mean pLDDT confidence, and per-residue pLDDT summary
  - Falls back to ColabFold URL hint if AlphaFold DB returns no entry
  - Never raises - pipeline continues regardless of structure availability
  - Does NOT download PDB files locally (no disk I/O in this node)

Data sources:
  AlphaFold DB REST API  https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}
  ColabFold server hint  https://colabfold.com (manual fallback only)

Output fields written to CandidateProtein:
  structure_source      "alphafold_db" | "colabfold_hint" | "unavailable"
  structure_pdb_path    AlphaFold CIF/PDB URL (not a local path in this node)
"""

import logging
import requests
from typing import List, Optional, Dict, Any

from src.models.candidate import CandidateProtein, ConfidenceTier

logger = logging.getLogger(__name__)

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
REQUEST_TIMEOUT = 15  # seconds


class StructureAgent:
    def __init__(self):
        self.stage_name = "structure_retrieval"

    def run(self, candidates: List[CandidateProtein]) -> List[CandidateProtein]:
        """
        Retrieve AlphaFold structure metadata for active candidates.
        Candidates with non-UniProt IDs (e.g. 'user_input') are skipped
        with a logged warning - not failed.
        """
        active = [c for c in candidates if c.status.value == "active"]
        logger.info(f"N5 StructureAgent: processing {len(active)} active candidates")

        for i, candidate in enumerate(active):
            logger.info(
                f"  [{i+1}/{len(active)}] Fetching structure: "
                f"{candidate.protein_name} ({candidate.protein_id})"
            )

            if not self._is_uniprot_id(candidate.protein_id):
                logger.warning(
                    f"  Skipping {candidate.protein_id} - not a UniProt accession. "
                    f"ColabFold manual submission required."
                )
                candidate.structure_source = "unavailable"
                candidate.add_decision(
                    stage=self.stage_name,
                    decision="skipped",
                    reasoning=(
                        f"protein_id '{candidate.protein_id}' is not a UniProt accession. "
                        f"AlphaFold DB lookup requires a valid UniProt ID. "
                        f"Submit sequence manually to ColabFold for structure prediction."
                    ),
                    structure_source="unavailable",
                )
                continue

            result = self._fetch_alphafold(candidate.protein_id)

            if result:
                candidate.structure_source = "alphafold_db"
                # structure_pdb_path stores the canonical CIF URL
                # (field name is legacy from local-path era - we store the remote URL)
                candidate.structure_pdb_path = result["cif_url"]
                candidate.add_decision(
                    stage=self.stage_name,
                    decision="structure_retrieved",
                    reasoning=(
                        f"AlphaFold DB entry found. "
                        f"Model version: {result['model_version']}. "
                        f"Mean pLDDT: {result['mean_plddt']:.1f}/100 "
                        f"({'high' if result['mean_plddt'] >= 70 else 'low'} confidence). "
                        f"Fragment coverage: {result['fragment_coverage']}."
                    ),
                    structure_source="alphafold_db",
                    alphafold_entry_id=result["entry_id"],
                    model_version=result["model_version"],
                    mean_plddt=result["mean_plddt"],
                    pdb_url=result["pdb_url"],
                    cif_url=result["cif_url"],
                    fragment_coverage=result["fragment_coverage"],
                )
                logger.info(
                    f"    AlphaFold hit: mean_pLDDT={result['mean_plddt']:.1f}, "
                    f"version={result['model_version']}"
                )

            else:
                candidate.structure_source = "unavailable"
                candidate.structure_pdb_path = None
                candidate.add_decision(
                    stage=self.stage_name,
                    decision="structure_unavailable",
                    reasoning=(
                        f"No AlphaFold DB entry for UniProt ID '{candidate.protein_id}'. "
                        f"Possible reasons: protein not in reviewed Swiss-Prot set, "
                        f"sequence too short, or not yet modelled. "
                        f"Downstream conformational B-cell predictions will use "
                        f"linear sequence only (reduced accuracy)."
                    ),
                    structure_source="unavailable",
                    colabfold_hint=(
                        f"https://colab.research.google.com/github/sokrypton/ColabFold/"
                        f"blob/main/AlphaFold2.ipynb"
                    ),
                )
                logger.warning(
                    f"    No AlphaFold entry for {candidate.protein_id}. "
                    f"Conformational epitope accuracy is reduced."
                )

        total_with_structure = sum(
            1 for c in candidates if c.structure_source == "alphafold_db"
        )
        logger.info(
            f"N5 complete: {total_with_structure}/{len(active)} candidates "
            f"have AlphaFold structures"
        )
        return candidates

    # --PRIVATE

    def _fetch_alphafold(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Query AlphaFold DB REST API.
        Returns a normalised dict or None on any error.

        API returns a list; we take the first entry (highest-version model).
        """
        url = f"{ALPHAFOLD_API}/{uniprot_id}"
        try:
            resp = requests.get(url, timeout=REQUEST_TIMEOUT)

            if resp.status_code == 404:
                return None  # legitimate miss, not an error

            resp.raise_for_status()
            data = resp.json()

            if not data:
                return None

            entry = data[0]  # first = latest model version

            # pLDDT summary - API returns meanPlddt directly on the entry
            mean_plddt = entry.get("meanPlddt") or 0.0

            return {
                "entry_id": entry.get("entryId", uniprot_id),
                "model_version": entry.get("latestVersion", "unknown"),
                "pdb_url": entry.get("pdbUrl", ""),
                "cif_url": entry.get("cifUrl", entry.get("pdbUrl", "")),
                "mean_plddt": round(float(mean_plddt), 2),
                # fraction of sequence covered by this model (1 fragment = full length)
                "fragment_coverage": (
                    f"{entry.get('uniprotStart', '?')}-{entry.get('uniprotEnd', '?')}"
                ),
            }

        except requests.exceptions.Timeout:
            logger.warning(f"AlphaFold DB timeout for {uniprot_id}")
            return None
        except requests.exceptions.RequestException as e:
            logger.warning(f"AlphaFold DB request failed for {uniprot_id}: {e}")
            return None
        except (KeyError, IndexError, ValueError) as e:
            logger.warning(f"AlphaFold DB response parse error for {uniprot_id}: {e}")
            return None

    @staticmethod
    def _is_uniprot_id(protein_id: str) -> bool:
        """
        Validates UniProt accession format.
        Accepted formats: P12345, A0A000, O00001 etc.
        Rejects: 'user_input', NCBI GI numbers, custom strings.
        """
        import re
        # UniProt accession: 6 or 10 alphanumeric chars, specific pattern
        pattern = r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
        return bool(re.match(pattern, protein_id.strip().upper()))


# Module-level instance - matches pattern of bcell_predictor.py
structure_agent = StructureAgent()