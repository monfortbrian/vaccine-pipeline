"""
ANTIGEN SCREENER AGENT

Tools:
  VaxiJen 2.0  - antigenicity score via local ACC implementation
                 (cloudflare blocks cloud IPs; ACC model runs locally)
  Phobius      - transmembrane topology and signal peptide prediction
                 (replaces PSORTb which requires Docker-in-Docker)

Early exit gate: if VaxiJen score < 0.3 AND sequence is cytoplasmic with
no signal peptide, the candidate is deprioritised with a clear audit record.
All other candidates proceed to Agent 3 regardless of score; the score is
recorded for scientific context, not used as a hard filter by default.

References:
  VaxiJen 2.0: Doytchinova & Flower (2007) BMC Bioinformatics 8:4
  Phobius: Kall et al. (2004) J Mol Biol 338:1027-1036
"""

import logging
from typing import List
from src.models.candidate import CandidateProtein, CandidateStatus

from src.utils.logger import get_logger
logger = get_logger("tope_deep.agents.Agent 2")  # use the correct agent name

VAXIJEN_LOW_THRESHOLD = 0.3   # below this + cytoplasmic = deprioritise


class AntigenScreenerAgent:
    """
    Antigen Screener

    Scores antigenicity and predicts subcellular localisation for all active candidates.
    Non-blocking: if either tool fails, candidate proceeds with missing score flagged.
    Deprioritises (does not discard) candidates with very low scores.
    """

    def __init__(self):
        self.stage_name = "antigen_screening"
        self._vaxijen = None
        self._phobius = None

    @property
    def vaxijen(self):
        if not self._vaxijen:
            from src.tools.vaxijen_client import vaxijen
            self._vaxijen = vaxijen
        return self._vaxijen

    @property
    def phobius(self):
        if not self._phobius:
            from src.tools.phobius_client import phobius
            self._phobius = phobius
        return self._phobius

    def run(
        self,
        candidates: List[CandidateProtein],
        organism_class: str = "bacteria",
    ) -> List[CandidateProtein]:
        logger.info("Agent 2: Antigen Screener - VaxiJen 2.0 ACC + Phobius")
        active = [c for c in candidates if c.status == CandidateStatus.ACTIVE]

        for i, c in enumerate(active):
            logger.info(f"   [{i+1}/{len(active)}] {c.protein_name}")
            vaxijen_score  = None
            vaxijen_method = "unavailable"
            localization   = "unknown"
            tm_helices     = 0
            failed_tools   = []

            try:
                vaxijen_score  = self.vaxijen.predict_antigenicity(c.sequence, organism_class)
                c.vaxijen_score = vaxijen_score
                vaxijen_method = (
                    "VaxiJen_2.0_real"
                    if self.vaxijen.is_server_available()
                    else "VaxiJen_2.0_ACC_local"
                )
            except Exception as e:
                logger.warning(f"      VaxiJen failed: {e}")
                failed_tools.append("vaxijen")

            try:
                result     = self.phobius.predict_transmembrane(c.sequence, c.protein_id)
                tm_helices = result.get("num_tm_helices", 0)
                c.tmhmm_helices = tm_helices
                if result.get("has_signal_peptide"):
                    localization = "secreted"
                elif tm_helices == 0:
                    localization = "cytoplasmic"
                elif tm_helices == 1:
                    localization = "single_pass_membrane"
                else:
                    localization = "multi_pass_membrane"
                c.psortb_localization = localization
            except Exception as e:
                logger.warning(f"      Phobius failed: {e}")
                failed_tools.append("phobius")

            # Deprioritise if clearly poor antigen and cytoplasmic
            deprioritised = (
                vaxijen_score is not None
                and vaxijen_score < VAXIJEN_LOW_THRESHOLD
                and localization == "cytoplasmic"
            )
            if deprioritised:
                c.status = CandidateStatus.DEPRIORITIZED

            c.add_decision(
                stage=self.stage_name,
                decision="deprioritised" if deprioritised else "screened",
                reasoning=(
                    f"VaxiJen={vaxijen_score:.3f} [organism={organism_class}, method={vaxijen_method}]. "
                    f"Phobius: TM_helices={tm_helices}, localization={localization}. "
                    + (f"Deprioritised: VaxiJen score below {VAXIJEN_LOW_THRESHOLD} with cytoplasmic localisation. "
                       f"Candidate will not proceed to Agent 3. " if deprioritised else "Advancing to Agent 3. ")
                    + (f"Tool failures (non-blocking): {', '.join(failed_tools)}. " if failed_tools else "")
                    + "Phobius used as PSORTb proxy (PSORTb requires Docker-in-Docker, unavailable on Railway)."
                ),
                vaxijen_score=vaxijen_score,
                vaxijen_method=vaxijen_method,
                organism_class=organism_class,
                tm_helices=tm_helices,
                phobius_localization=localization,
                localization_tool="Phobius_2.0",
                deprioritised=deprioritised,
            )
            logger.info(
                f"      VaxiJen={vaxijen_score:.3f if vaxijen_score else 'n/a'} "
                f"loc={localization} TM={tm_helices} "
                f"{'[DEPRIORITISED]' if deprioritised else ''}"
            )

        active_count = sum(1 for c in candidates if c.status == CandidateStatus.ACTIVE)
        logger.info(f"Agent 2 complete: {active_count} active candidates proceed to Agent 3")
        return candidates


antigen_screener = AntigenScreenerAgent()