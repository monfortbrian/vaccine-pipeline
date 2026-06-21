"""
EXPERIMENT PLANNER AGENT
Generates a prioritized wet-lab validation roadmap for each candidate.


"""

import os
import re
import json
import time
import logging
import requests
from typing import List, Dict, Any, Optional

from src.models.candidate import CandidateProtein, ConfidenceTier
from src.utils.logger import get_logger

logger = get_logger("tope_deep.agents.Agent 10")


# ── Template fallback (always returns non-None) ───────────────────────────────

def _template_plan(candidate: CandidateProtein, context: Dict) -> Dict:
    """
    Structured wet-lab plan without Claude.
    GUARANTEE: always returns a non-None, non-empty dict.
    Reference: Sette & Rappuoli (2010) Immunity 33(4):530-541.
    """
    high_ctl = [
        ep for ep in candidate.ctl_epitopes
        if ep.confidence_tier == ConfidenceTier.HIGH
        and ep.allergenicity_safe is True
    ][:5]

    mouse_reactive = [
        ep for ep in candidate.ctl_epitopes
        if ep.tool_outputs.get("animal_model_alleles")
        and ep.allergenicity_safe is True
    ][:3]

    safe_count   = context.get("safe_count", 0)
    failure_sigs = context.get("failure_signals", [])
    prior_valid  = context.get("prior_validated", False)

    return {
        "priority_reasoning": (
            f"{len(high_ctl)} high-confidence CTL epitopes identified for {candidate.protein_name}, "
            f"all safety-cleared (allergenicity and haemolytic screens passed). "
            f"{'Prior experimental evidence in literature supports this target.' if prior_valid else 'No prior experimental validation found in published literature.'}"
        ),
        "priority_epitopes": [
            {
                "sequence":   ep.sequence,
                "hla_allele": ep.hla_allele,
                "ic50_nm":    ep.ic50_nm,
                "confidence": ep.confidence_tier.value,
                "rationale":  "High-confidence MHC-I binder, allergenicity and toxicity cleared",
            }
            for ep in high_ctl
        ],
        "phase_1": {
            "timeline":     "Week 1–4",
            "key_assay":    "PBMC IFN-γ ELISpot",
            "protocol_notes": (
                f"Synthesise top {min(5, len(high_ctl))} CTL peptides (9–11mer, >95% purity). "
                f"Stimulate PBMCs from HLA-typed donors at 10 μg/mL. "
                f"Readout: IFN-γ spot-forming cells per 10⁶ PBMCs. "
                f"Threshold: ≥2× background AND ≥50 SFU/10⁶ = positive."
            ),
            "expected_cost_usd": "$3,000–$8,000",
            "go_criteria":  "≥2 of top 5 epitopes elicit IFN-γ in ≥2 HLA-matched donors",
            "reference":    "Janetzki et al. (2015) Cancer Immunol Immunother 64:1695–1703",
        },
        "phase_2": {
            "timeline":     "Month 1–3",
            "key_assay":    (
                "Murine immunisation (C57BL/6, BALB/c)"
                if mouse_reactive else
                "Skip, no H-2 cross-reactive epitopes predicted"
            ),
            "protocol_notes": (
                f"Immunise C57BL/6 (H-2b) mice with peptide + CFA. "
                f"Boost at day 14, harvest splenocytes at day 21. "
                f"ICS: CD8+/IFN-γ+/TNF-α+. "
                f"Cross-reactive epitopes: {', '.join(ep.sequence for ep in mouse_reactive) or 'none'}."
            ) if mouse_reactive else "No H-2 cross-reactive epitopes, proceed directly to phase 3.",
            "expected_cost_usd": "$5,000–$15,000" if mouse_reactive else "$0 (skipped)",
            "go_criteria":  "CD8+ T-cell response in ≥3/5 mice, ≥0.5% antigen-specific cells",
        },
        "phase_3": {
            "timeline":     "Month 3–6",
            "key_assay":    "MHC-I tetramer staining",
            "protocol_notes": (
                "Produce HLA-A*02:01 tetramers for top CTL epitopes. "
                "Stain PBMCs from HLA-A*02:01+ donors. "
                "Confirms antigen-specific T-cell frequency directly."
            ),
            "expected_cost_usd": "$10,000–$30,000",
            "go_criteria":  "Tetramer+ CD8+ T cells >0.1% of total CD8+ in ≥2 donors",
            "reference":    "Klenerman et al. (2002) Nat Rev Immunol 2:263–272",
        },
        "critical_risks": failure_sigs + (
            ["Prior experimental validation exists, review published HLA restrictions before synthesis"]
            if prior_valid else []
        ) + (
            [f"Only {safe_count} CTL epitopes safety-cleared, limited synthesis options"]
            if safe_count < 3 else []
        ),
        "elispot_protocol": (
            f"1. Synthesise top {min(5, len(high_ctl))} CTL epitopes as 9–11mer peptides (>95% purity, lyophilised).\n"
            f"2. Coat PVDF plates (Millipore MSIP) with anti-IFN-γ antibody overnight at 4°C.\n"
            f"3. Rest PBMCs 4h at 37°C, then plate 2×10⁵/well.\n"
            f"4. Stimulate with peptide (10 μg/mL) + brefeldin A for 18h.\n"
            f"5. Detect with biotinylated anti-IFN-γ → streptavidin-ALP → BCIP substrate.\n"
            f"6. Count spots: ≥2× background AND ≥50 SFU/10⁶ = positive response."
        ),
        "immunisation_schedule": (
            "Day 0:  Prime - peptide (100 μg) + CFA, subcutaneous.\n"
            "Day 14: Boost - peptide (50 μg) + IFA, subcutaneous.\n"
            "Day 21: Harvest -  terminal bleed, splenocyte isolation.\n"
            "Day 22: ICS and ELISpot readout."
        ),
        "nhp_plan": (
            "NHP validation indicated if murine model (phase 2) meets go criteria.\n"
            "Macaque (Macaca mulatta): use Mamu-A*01 and Mamu-A*02 tetramers.\n"
            "Dose: 200 μg peptide + Montanide ISA 51, weeks 0, 4, 8.\n"
            "PBMC isolation at week 10 via Ficoll gradient.\n"
            "Primary readout: ICS (CD8+/IFN-γ+) and tetramer staining."
        ) if mouse_reactive else (
            "NHP model recommended directly (no H-2 cross-reactive epitopes for murine screening).\n"
            "Proceed to Mamu allele-based ELISpot first."
        ),
        "estimated_cost_usd": {
            "phase_1": "$3,000–$8,000",
            "phase_2": "$5,000–$15,000" if mouse_reactive else "$0 (skipped)",
            "phase_3": "$10,000–$30,000",
            "total":   "$18,000–$53,000 for full validation",
        },
        "references": [
            "Sette & Rappuoli (2010) Immunity 33(4):530–541",
            "Janetzki et al. (2015) Cancer Immunol Immunother 64:1695–1703",
            "Klenerman et al. (2002) Nat Rev Immunol 2:263–272",
            "Seder et al. (2008) Nat Immunol 9:239–245",
        ],
        "generated_by": "TOPE_DEEP template",
    }


# ── Claude-generated plan ─────────────────────────────────────────────────────

def _claude_plan(candidate: CandidateProtein, context: Dict) -> Optional[Dict]:
    """
    Generate structured wet-lab plan via Claude API.
    Returns None on any failure, caller must fall back to template.
    FIX: explicit JSON parse with fallback, no silent None return on bad JSON.
    """
    api_key = os.getenv("ANTHROPIC_API_KEY", "")
    if not api_key:
        logger.info("Agent 10: ANTHROPIC_API_KEY not set, using template plan")
        return None

    high_ctl = [
        ep for ep in candidate.ctl_epitopes
        if ep.confidence_tier == ConfidenceTier.HIGH
        and ep.allergenicity_safe is True
    ][:8]

    mouse_reactive = [
        ep.sequence for ep in candidate.ctl_epitopes
        if ep.tool_outputs.get("animal_model_alleles") and ep.allergenicity_safe is True
    ][:5]

    mamu_reactive = [
        ep.sequence for ep in candidate.ctl_epitopes
        if ep.tool_outputs.get("mamu_alleles") and ep.allergenicity_safe is True
    ][:3]

    prompt = f"""You are an expert computational immunologist at a vaccine research institute.
Generate a structured wet-lab validation plan for this vaccine candidate.

CANDIDATE: {candidate.protein_name} ({candidate.protein_id})
Sequence length: {getattr(candidate, 'sequence_length', len(candidate.sequence or ''))} aa
VaxiJen antigenicity: {candidate.vaxijen_score}

TOP CTL EPITOPES (safety-cleared, high confidence):
{chr(10).join(f"  {ep.sequence} | {ep.hla_allele} | IC50={ep.ic50_nm}nM | rank={ep.percentile_rank}" for ep in high_ctl)}

ANIMAL MODEL CROSS-REACTIVE:
  Mouse H-2: {', '.join(mouse_reactive) if mouse_reactive else 'none'}
  Macaque Mamu: {', '.join(mamu_reactive) if mamu_reactive else 'none'}

POPULATION COVERAGE:
  Global: {context.get('global_coverage_pct')}% | African: {context.get('african_coverage_pct')}%

SAFETY: {context.get('safe_count')} epitopes cleared, {context.get('fail_count')} failed

LITERATURE:
  Prior validation: {context.get('prior_validated')}
  Failure signals: {', '.join(context.get('failure_signals', [])) or 'none'}

LAB CONSTRAINTS: {context.get('lab_constraints', 'standard academic lab')}

Return ONLY a JSON object with this exact structure (no markdown, no preamble):
{{
  "priority_reasoning": "2-3 sentences on why these epitopes are prioritised",
  "phase_1": {{
    "timeline": "Week 1-4",
    "key_assay": "PBMC IFN-γ ELISpot",
    "protocol_notes": "specific adaptations for this candidate",
    "expected_cost_usd": "range",
    "go_criteria": "measurable threshold"
  }},
  "phase_2": {{
    "timeline": "Month 1-3",
    "key_assay": "murine immunisation or skip",
    "protocol_notes": "specific notes",
    "expected_cost_usd": "range",
    "go_criteria": "measurable threshold"
  }},
  "phase_3": {{
    "timeline": "Month 3-6",
    "key_assay": "MHC tetramer or NHP",
    "protocol_notes": "specific notes",
    "expected_cost_usd": "range",
    "go_criteria": "measurable threshold"
  }},
  "elispot_protocol": "step-by-step ELISpot protocol for this candidate",
  "immunisation_schedule": "immunisation timeline with doses",
  "nhp_plan": "NHP validation plan or note if not needed",
  "critical_risks": ["list of specific risks"],
  "regulatory_note": "one sentence on regulatory pathway"
}}"""

    try:
        resp = requests.post(
            "https://api.anthropic.com/v1/messages",
            headers={
                "x-api-key":          api_key,
                "anthropic-version":  "2023-06-01",
                "content-type":       "application/json",
            },
            json={
                "model":    "claude-sonnet-4-6",
                "max_tokens": 1200,
                "messages": [{"role": "user", "content": prompt}],
            },
            timeout=30,
        )
        resp.raise_for_status()
        text = resp.json()["content"][0]["text"].strip()

        # Strip markdown fences if Claude adds them despite instructions
        text = re.sub(r'^```json\s*', '', text)
        text = re.sub(r'\s*```$',     '', text)

        parsed = json.loads(text)
        parsed["generated_by"] = "claude-sonnet-4-6"
        return parsed

    except json.JSONDecodeError as e:
        logger.warning(f"Agent 10: Claude returned invalid JSON: {e} using template")
        return None
    except Exception as e:
        logger.warning(f"Agent 10: Claude plan failed: {e} using template")
        return None


# ── Agent ─────────────────────────────────────────────────────────────────────

class ExperimentPlannerAgent:
    """
    Experiment Planner Agent
    Uses Claude when ANTHROPIC_API_KEY is set.
    Always falls back to template, never produces empty output.
    """

    def __init__(self):
        self.stage_name = "experiment_planning"

    def run(
        self,
        candidates: List[CandidateProtein],
        lab_constraints: str = "standard academic lab",
    ) -> List[CandidateProtein]:
        logger.info("Agent 10: Starting experiment planning")
        active     = [c for c in candidates if c.status.value == "active"]
        has_claude = bool(os.getenv("ANTHROPIC_API_KEY", ""))
        logger.info(f"   {len(active)} candidates | Claude={'yes' if has_claude else 'no (template fallback)'}")

        for i, candidate in enumerate(active):
            logger.info(f"   [{i+1}/{len(active)}] {candidate.protein_name}")
            start = time.time()

            # Gather context from previous agent decisions
            safety_decision = next(
                (d for d in candidate.decisions if d.get("stage") == "safety_filter"), {}
            )
            lit_decision = next(
                (d for d in candidate.decisions if d.get("stage") == "literature_search"), {}
            )

            context = {
                "safe_count":           safety_decision.get("safe_count", 0),
                "unscored_count":       safety_decision.get("unscored_count", 0),
                "fail_count":           safety_decision.get("fail_count", 0),
                "global_coverage_pct":  round((candidate.hla_coverage_global or 0) * 100, 1),
                "african_coverage_pct": round((candidate.hla_coverage_africa or 0) * 100, 1),
                "prior_validated":      lit_decision.get("prior_validated", False),
                "failure_signals":      lit_decision.get("failure_signals", []),
                "literature_summary":   lit_decision.get("literature_summary"),
                "lab_constraints":      lab_constraints,
            }

            # Try Claude first, template is guaranteed fallback
            plan = _claude_plan(candidate, context)
            if plan is None:
                plan = _template_plan(candidate, context)
                generation_method = "template"
            else:
                generation_method = "claude-sonnet-4-6"

            # Guarantee plan is never None or empty
            if not plan:
                plan = _template_plan(candidate, context)
                generation_method = "template"

            elapsed = round(time.time() - start, 1)

            # Scientific reasoning only no internal method names
            priority_count = len([
                ep for ep in candidate.ctl_epitopes
                if ep.confidence_tier == ConfidenceTier.HIGH and ep.allergenicity_safe is True
            ])

            candidate.add_decision(
                stage=self.stage_name,
                decision="plan_generated",
                reasoning=(
                    f"Wet-lab validation plan generated for {candidate.protein_name}. "
                    f"{priority_count} high-confidence, safety-cleared CTL epitopes prioritised for synthesis. "
                    f"Plan covers ELISpot (phase 1), murine model (phase 2), MHC tetramer (phase 3). "
                    f"Lab constraints: {lab_constraints}."
                ),
                # FIX: use 'plan' key, this is what the frontend reads
                # from decisions[].plan in results page ExperimentSection
                plan=plan,
                lab_constraints=lab_constraints,
                plan_time_s=elapsed,
            )

            # Safely set attribute, CandidateProtein may not define it
            try:
                candidate.experiment_plan = plan
            except AttributeError:
                pass

            logger.info(f"      Plan generated via {generation_method} in {elapsed}s")

        logger.info("Agent 10: Experiment planning complete")
        return candidates


experiment_planner = ExperimentPlannerAgent()