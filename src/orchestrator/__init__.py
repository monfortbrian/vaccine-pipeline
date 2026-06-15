"""
TOPE_DEEP Pipeline Orchestrator

Chains all 10 agents in sequence. This is the single place
where pipeline execution logic lives. main.py calls this.

Agent sequence:
  DataCuratorAgent       : validate, record provenance
  AntigenScreenerAgent   : VaxiJen + Phobius (inline in orchestrator, extracted below)
  TCellPredictorAgent    : NetMHCpan 4.1 + NetMHCIIpan 4.3 + MHCflurry fallback
  BCellPredictorAgent    : IEDB BepiPred 2.0
  StructureAgent         : AlphaFold DB
  SafetyFilterAgent      : FAO/WHO + AllerTOP + HemoPI + human homology
  CoverageAgent          : IEDB population coverage tool / AFND 2020 fallback
  ConstructDesignerAgent : ProtParam, configurable adjuvant, linker assembly
  LiteratureAgent        : PubMed + Qdrant + Claude synthesis
  ExperimentPlannerAgent : Claude API + template fallback
"""

import time
import logging
from typing import List, Dict, Any, Optional
from src.models.candidate import CandidateProtein, CandidateStatus

logger = logging.getLogger(__name__)


class PipelineOrchestrator:
    def __init__(self):
        self._agents_loaded = False

    def _load_agents(self):
        if self._agents_loaded:
            return
        from src.agents.data_curator import DataCuratorAgent
        from src.agents.tcell_predictor import TCellPredictorAgent
        from src.agents.bcell_predictor import BCellPredictorAgent
        from src.agents.structure_agent import StructureAgent
        from src.agents.safety_filter import SafetyFilterAgent
        from src.agents.coverage_agent import CoverageAgent
        from src.agents.construct_designer import ConstructDesignerAgent
        from src.agents.literature_agent import LiteratureAgent
        from src.agents.experiment_planner import ExperimentPlannerAgent

        self.n1 = DataCuratorAgent()
        self.n3 = TCellPredictorAgent()
        self.n4 = BCellPredictorAgent()
        self.n5 = StructureAgent()
        self.n6 = SafetyFilterAgent()
        self.n7 = CoverageAgent()
        self.n8 = ConstructDesignerAgent()
        self.n9 = LiteratureAgent()
        self.n10 = ExperimentPlannerAgent()
        self._agents_loaded = True

    def run(
        self,
        run_id: str,
        candidates: List[CandidateProtein],
        run_config: Dict[str, Any],
        progress_callback=None,
    ) -> Dict[str, Any]:
        """
        Execute the full 10-agent pipeline.

        run_config keys:
          organism_class   str   "bacteria" | "virus" | "parasite"
          input_type       str   "pathogen" | "uniprot_id" | "sequence"
          run_safety       bool
          run_coverage     bool
          run_literature   bool
          run_experiment   bool
          lab_constraints  str

        progress_callback(node, progress, message) optional, called after each agent.
        """
        self._load_agents()

        organism_class  = run_config.get("organism_class",  "bacteria")
        input_type      = run_config.get("input_type",      "unknown")
        run_safety      = run_config.get("run_safety",      True)
        run_coverage    = run_config.get("run_coverage",    True)
        run_literature  = run_config.get("run_literature",  True)
        run_experiment  = run_config.get("run_experiment",  True)
        lab_constraints = run_config.get("lab_constraints", "standard academic lab")

        def _cb(node, progress, message):
            if progress_callback:
                progress_callback(node, progress, message)

        timings: Dict[str, Optional[float]] = {
            "n1_curation":    None,
            "n2_screening":   None,
            "n3_tcell":       None,
            "n4_bcell":       None,
            "n5_structure":   None,
            "n6_safety":      None,
            "n7_coverage":    None,
            "n8_construct":   None,
            "n9_literature":  None,
            "n10_experiment": None,
        }

        start_total = time.time()

        # Data Curator
        _cb("N1", 0.02, "Recording data provenance...")
        t = time.time()
        candidates = self.n1.run(candidates, organism_class=organism_class, input_type=input_type)
        timings["n1_curation"] = round(time.time() - t, 2)
        _cb("N1", 0.05, f"N1 complete: {sum(1 for c in candidates if c.status == CandidateStatus.ACTIVE)} candidates active")

        # Antigen Screener (VaxiJen + Phobius)
        _cb("N2", 0.06, "Antigen screening: VaxiJen 2.0 + Phobius...")
        t = time.time()
        candidates = self._run_n2(candidates, organism_class)
        timings["n2_screening"] = round(time.time() - t, 2)
        _cb("N2", 0.10, "N2 complete")

        # TCell Predictor
        _cb("N3", 0.11, "Predicting T-cell epitopes: NetMHCpan 4.1, NetMHCIIpan 4.3...")
        t = time.time()
        candidates = self.n3.run(candidates)
        timings["n3_tcell"] = round(time.time() - t, 2)
        _cb("N3", 0.35, "N3 complete")

        # BCell Predictor
        _cb("N4", 0.36, "Predicting B-cell epitopes: BepiPred 2.0...")
        t = time.time()
        candidates = self.n4.run(candidates)
        timings["n4_bcell"] = round(time.time() - t, 2)
        _cb("N4", 0.50, "N4 complete")

        # Structure Agent
        _cb("N5", 0.51, "Retrieving 3D structures: AlphaFold DB...")
        t = time.time()
        candidates = self.n5.run(candidates)
        timings["n5_structure"] = round(time.time() - t, 2)
        _cb("N5", 0.58, "N5 complete")

        # Safety Filter
        if run_safety:
            _cb("N6", 0.59, "Safety screening: FAO/WHO, AllerTOP, HemoPI, human homology...")
            t = time.time()
            candidates = self.n6.run(candidates)
            timings["n6_safety"] = round(time.time() - t, 2)
        _cb("N6", 0.72, "N6 complete" if run_safety else "N6 skipped")

        # Coverage Agent
        if run_coverage:
            _cb("N7", 0.73, "Calculating HLA population coverage: IEDB / AFND 2020...")
            t = time.time()
            candidates = self.n7.run(candidates)
            timings["n7_coverage"] = round(time.time() - t, 2)
        _cb("N7", 0.82, "N7 complete" if run_coverage else "N7 skipped")

        # Construct Designer
        _cb("N8", 0.83, "Assembling multi-epitope construct: ProtParam, adjuvant, linkers...")
        t = time.time()
        candidates, construct_report = self.n8.run(candidates)
        timings["n8_construct"] = round(time.time() - t, 2)
        _cb("N8", 0.88, "N8 complete")

        # Literature Agent
        if run_literature:
            _cb("N9", 0.89, "Searching published literature: PubMed, Qdrant...")
            t = time.time()
            try:
                candidates = self.n9.run(candidates, run_id=run_id)
            except Exception as e:
                logger.warning(f"N9 failed (non-blocking): {e}")
            timings["n9_literature"] = round(time.time() - t, 2)
        _cb("N9", 0.94, "N9 complete" if run_literature else "N9 skipped")

        # Experiment Planner
        if run_experiment:
            _cb("N10", 0.95, "Generating wet-lab validation roadmap: Claude API...")
            t = time.time()
            try:
                candidates = self.n10.run(candidates, lab_constraints=lab_constraints)
            except Exception as e:
                logger.warning(f"N10 failed (non-blocking): {e}")
            timings["n10_experiment"] = round(time.time() - t, 2)
        _cb("N10", 0.99, "N10 complete" if run_experiment else "N10 skipped")

        total = round(time.time() - start_total, 1)
        timings["total_seconds"] = total

        logger.info(
            f"Pipeline complete: run={run_id} total={total}s "
            f"N3={timings['n3_tcell']}s N4={timings['n4_bcell']}s "
            f"N6={timings['n6_safety']}s N9={timings['n9_literature']}s"
        )

        return {
            "candidates":       candidates,
            "construct_report": construct_report,
            "timing":           timings,
        }

    def _run_n2(
        self, candidates: List[CandidateProtein], organism_class: str
    ) -> List[CandidateProtein]:
        """
        Antigen Screener.
        VaxiJen 2.0 (local ACC fallback) + Phobius transmembrane topology.
        Non-blocking: if tools fail, candidates proceed with missing scores flagged.
        """
        try:
            from src.tools.vaxijen_client import vaxijen
            from src.tools.phobius_client import phobius
        except ImportError as e:
            logger.warning(f"N2 tools unavailable: {e}")
            return candidates

        active = [c for c in candidates if c.status == CandidateStatus.ACTIVE]
        for c in active:
            try:
                score = vaxijen.predict_antigenicity(c.sequence, organism_class)
                c.vaxijen_score = score
                method = "VaxiJen_2.0_real" if vaxijen.is_server_available() else "VaxiJen_2.0_ACC_local"

                result = phobius.predict_transmembrane(c.sequence, c.protein_id)
                tm = result.get("num_tm_helices", 0)
                c.tmhmm_helices = tm

                if result.get("has_signal_peptide"):
                    c.psortb_localization = "secreted"
                elif tm == 0:
                    c.psortb_localization = "cytoplasmic"
                elif tm == 1:
                    c.psortb_localization = "single_pass_membrane"
                else:
                    c.psortb_localization = "multi_pass_membrane"

                c.add_decision(
                    stage="antigen_screening",
                    decision="screened",
                    reasoning=(
                        f"VaxiJen={c.vaxijen_score:.3f} [organism={organism_class}, method={method}]. "
                        f"Phobius: TM_helices={tm}, localization={c.psortb_localization} "
                        f"(Kall et al. 2004). PSORTb unavailable on Railway, Phobius used as proxy."
                    ),
                    vaxijen_score=c.vaxijen_score,
                    vaxijen_method=method,
                    organism_class=organism_class,
                    tm_helices=tm,
                    phobius_localization=c.psortb_localization,
                    localization_tool="Phobius_2.0",
                )
                logger.info(
                    f"   N2: {c.protein_name} VaxiJen={c.vaxijen_score:.3f} "
                    f"[{method}] TM={tm} loc={c.psortb_localization}"
                )
            except Exception as e:
                logger.warning(f"   N2 failed for {c.protein_name} (continuing): {e}")
                c.add_decision(
                    stage="antigen_screening",
                    decision="screening_failed",
                    reasoning=f"N2 antigen screening failed: {e}. Candidate proceeds unscored.",
                )
        return candidates


orchestrator = PipelineOrchestrator()
