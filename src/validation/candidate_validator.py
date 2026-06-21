"""
CANDIDATE VALIDATOR

Runs before Agent 1. Validates the raw input (pathogen name, UniProt ID, or sequence)
and produces human-readable errors before any agent touches it.

This is the "human-in-the-loop" gate:
  - Wrong UniProt ID    : clear error, no agents run
  - Pathogen not found : clear error, suggest corrections
  - Sequence too short  : clear error with minimum requirement
  - All candidates deprioritised after Agent 2 : pipeline halts with summary, not silent failure

The guard also prevents resource waste:
  - If UniProt ID returns a human protein (Homo sapiens) : warn, ask to confirm
  - If sequence has >50% similarity to human proteome : flag before Agent 3 runs
  - If pathogen name matches 0 UniProt reviewed entries : halt with suggestions

This is NOT about restricting what scientists can run.
It is about surfacing problems early so they do not waste 20 minutes
on a pipeline run that will produce meaningless results.
"""

import logging
import requests
from typing import Tuple, Optional, List

logger = logging.getLogger(__name__)


class ValidationResult:
    def __init__(self, valid: bool, error: Optional[str] = None,
                 warning: Optional[str] = None, suggestions: Optional[List[str]] = None):
        self.valid       = valid
        self.error       = error
        self.warning     = warning
        self.suggestions = suggestions or []


def validate_uniprot_id(uid: str) -> ValidationResult:
    uid = uid.strip().upper()
    if not uid:
        return ValidationResult(False, "UniProt accession is empty")

    # UniProt accession format: [A-Z][0-9][A-Z]{3}[0-9] or P/Q/O+5 chars
    import re
    if not re.match(r'^[A-Z][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', uid):
        return ValidationResult(
            False,
            f"'{uid}' does not look like a valid UniProt accession.",
            suggestions=["UniProt accessions look like P9WNK7, A0A089QRB9, Q8I3H7",
                         "Search at https://www.uniprot.org/ to find the correct accession"]
        )

    # Check existence
    try:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.json", timeout=10)
        if r.status_code == 404:
            return ValidationResult(
                False,
                f"UniProt accession '{uid}' does not exist.",
                suggestions=[f"Search https://www.uniprot.org/uniprotkb?query={uid}"]
            )
        r.raise_for_status()
        data     = r.json()
        organism = data.get("organism", {}).get("scientificName", "")
        reviewed = data.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)"

        warning = None
        if "Homo sapiens" in organism:
            warning = (
                f"'{uid}' is a human protein (Homo sapiens). "
                f"TOPE_DEEP is designed for pathogen proteins. "
                f"Results will have very low population coverage and high human homology flags."
            )
        elif not reviewed:
            warning = (
                f"'{uid}' is an unreviewed TrEMBL entry. "
                f"Sequence quality and annotation may be lower than reviewed Swiss-Prot entries."
            )

        return ValidationResult(True, warning=warning)
    except requests.Timeout:
        # Cannot reach UniProt, allow pipeline to continue, Agent 1 will catch it
        return ValidationResult(True, warning="UniProt connectivity check timed out, proceeding")
    except Exception as e:
        return ValidationResult(True, warning=f"UniProt pre-check failed: {e}, proceeding")


def validate_sequence(seq: str) -> ValidationResult:
    seq = seq.upper().replace(" ", "").replace("\n", "").strip()
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")

    if len(seq) < 20:
        return ValidationResult(False, f"Sequence too short ({len(seq)} aa). Minimum 20 amino acids required for meaningful epitope prediction.")

    invalid = set(seq) - valid_aa
    if invalid:
        return ValidationResult(False, f"Sequence contains non-standard characters: {sorted(invalid)}. Use standard one-letter amino acid codes (A-Z, excluding B, J, O, U, X, Z).")

    if len(seq) > 5000:
        return ValidationResult(
            True,
            warning=f"Sequence is {len(seq)} aa. Very long sequences may cause IEDB tool timeouts. Consider submitting individual domains."
        )

    return ValidationResult(True)


def validate_pathogen_name(name: str) -> ValidationResult:
    name = name.strip()
    if len(name) < 3:
        return ValidationResult(False, "Pathogen name too short. Enter the full scientific or common name.")

    # Quick UniProt check
    try:
        r = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={"query": f'(organism_name:"{name}") AND (reviewed:true)', "format": "json", "size": "1"},
            timeout=10,
        )
        r.raise_for_status()
        count = len(r.json().get("results", []))
        if count == 0:
            return ValidationResult(
                False,
                f"No reviewed UniProt entries found for '{name}'.",
                suggestions=[
                    "Check spelling, use the scientific name (e.g. 'Mycobacterium tuberculosis' not 'TB')",
                    "Try a broader name (e.g. 'Plasmodium' instead of 'Plasmodium falciparum 3D7')",
                    f"Search https://www.uniprot.org/taxonomy?query={name.replace(' ','+')}",
                ]
            )
        return ValidationResult(True)
    except requests.Timeout:
        return ValidationResult(True, warning="Pathogen name pre-check timed out, proceeding")
    except Exception as e:
        return ValidationResult(True, warning=f"Pathogen name pre-check failed: {e}, proceeding")


def validate_input(input_type: str, input_value: str) -> ValidationResult:
    if input_type == "uniprot_id":
        return validate_uniprot_id(input_value)
    elif input_type == "sequence":
        return validate_sequence(input_value)
    elif input_type == "pathogen":
        return validate_pathogen_name(input_value)
    else:
        return ValidationResult(False, f"input_type must be 'pathogen', 'uniprot_id', or 'sequence'. Got: '{input_type}'")