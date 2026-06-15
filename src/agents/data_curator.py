"""
DATA CURATOR AGENT

Responsibility: record data provenance for each candidate protein.
The actual protein fetching happens in main.py (before the pipeline starts)
because it requires HTTP calls with different logic per input type.

N1's job within the pipeline:
  - Validate sequences
  - Record provenance decision (source, length, organism class)
  - Flag sequences that are too short, invalid, or highly repetitive
  - Mark candidates ACTIVE or DISCARDED based on basic quality gates

N1 does NOT:
  - Call PSORTb (unavailable on Railway, N2 uses Phobius instead)
  - Call Claude for antigen decisions (non-reproducible, expensive)
  - Run VaxiJen (that is N2's scope)

Consistent with all other agents: run(candidates) -> candidates
"""

import logging
from typing import List
from src.models.candidate import CandidateProtein, CandidateStatus, ConfidenceTier

from src.utils.logger import get_logger
logger = get_logger("tope_deep.agents.N1")  # use the correct agent name

# Quality gates
MIN_SEQUENCE_LENGTH = 20    # below this: no meaningful epitope prediction possible
MAX_SEQUENCE_LENGTH = 5000  # above this: IEDB call may time out; warn and continue
VALID_AMINO_ACIDS   = set("ACDEFGHIKLMNPQRSTVWY")
MAX_REPETITIVE_RUN  = 10    # consecutive identical amino acids = flag


class DataCuratorAgent:
    """
    Data Curator

    Records provenance and validates sequences for all candidate proteins.
    Marks candidates DISCARDED if they fail minimum quality gates.
    All other candidates proceed ACTIVE to N2.
    """

    def __init__(self):
        self.stage_name = "data_curation"

    def run(
        self,
        candidates: List[CandidateProtein],
        organism_class: str = "bacteria",
        input_type: str = "unknown",
    ) -> List[CandidateProtein]:
        logger.info(f"N1: Data Curator {len(candidates)} candidates")

        for i, c in enumerate(candidates):
            logger.info(
                f"   [{i+1}/{len(candidates)}] {c.protein_name} "
                f"({c.protein_id}) from {c.source}, {len(c.sequence)} aa"
            )

            flags, discard_reason = self._validate(c)

            if discard_reason:
                c.status = ConfidenceTier.UNCERTAIN  # type: ignore
                c.status = CandidateStatus.DISCARDED
                decision = "discarded"
                reasoning = (
                    f"Candidate discarded at N1: {discard_reason}. "
                    f"Flags: {', '.join(flags)}. "
                    f"No further processing."
                )
            else:
                c.status = CandidateStatus.ACTIVE
                decision = "protein_loaded"
                reasoning = (
                    f"Protein loaded from {c.source}. "
                    f"UniProt ID: {c.protein_id}. "
                    f"Sequence length: {len(c.sequence)} aa. "
                    f"Input type: {input_type}. "
                    f"Organism class inferred: {organism_class}. "
                    + (f"Flags (non-blocking): {', '.join(flags)}. " if flags else "No quality flags. ")
                    + "Advancing to N2 Antigen Screener."
                )
                if len(c.sequence) > MAX_SEQUENCE_LENGTH:
                    logger.warning(
                        f"   N1 WARNING: {c.protein_name} is {len(c.sequence)} aa "
                        f"IEDB calls may be slow or time out."
                    )

            c.add_decision(
                stage=self.stage_name,
                decision=decision,
                reasoning=reasoning,
                source=c.source,
                sequence_length=len(c.sequence),
                input_type=input_type,
                organism_class=organism_class,
                quality_flags=flags,
            )

        active  = sum(1 for c in candidates if c.status == CandidateStatus.ACTIVE)
        discard = sum(1 for c in candidates if c.status == CandidateStatus.DISCARDED)
        logger.info(f"N1 complete: {active} active, {discard} discarded")
        return candidates

    # ── VALIDATION ────────────────────────────────────────────────────────────

    def _validate(self, candidate: CandidateProtein):
        """
        Returns (flags, discard_reason).
        discard_reason is None if candidate should proceed.
        flags are non-blocking notes attached to the audit trail.
        """
        seq    = candidate.sequence.upper().strip()
        flags  = []

        # Hard failures, discard
        if len(seq) < MIN_SEQUENCE_LENGTH:
            return (["too_short"], f"Sequence too short ({len(seq)} aa, minimum {MIN_SEQUENCE_LENGTH})")

        invalid = set(seq) - VALID_AMINO_ACIDS
        if invalid:
            return (["invalid_amino_acids"], f"Non-standard amino acids found: {sorted(invalid)}")

        if not seq:
            return (["empty_sequence"], "Empty sequence after stripping")

        # Soft warnings, flag but proceed
        if len(seq) > MAX_SEQUENCE_LENGTH:
            flags.append(f"very_long_{len(seq)}aa")

        if self._is_repetitive(seq):
            flags.append("repetitive_sequence")

        if self._has_signal_peptide_signature(seq):
            flags.append("possible_signal_peptide")

        return (flags, None)

    @staticmethod
    def _is_repetitive(seq: str) -> bool:
        if len(seq) < 20:
            return False
        run = max_run = 1
        for i in range(1, len(seq)):
            run = run + 1 if seq[i] == seq[i-1] else 1
            max_run = max(max_run, run)
        return max_run > MAX_REPETITIVE_RUN or max_run > len(seq) * 0.3

    @staticmethod
    def _has_signal_peptide_signature(seq: str) -> bool:
        if len(seq) < 25:
            return False
        n_terminal   = seq[:25]
        hydrophobic  = sum(1 for aa in n_terminal if aa in "AILMFWYV")
        return hydrophobic / len(n_terminal) > 0.5


data_curator = DataCuratorAgent()