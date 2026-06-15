"""
EXPERIMENT PLANNER AGENT
Generates a prioritized wet-lab validation roadmap for each candidate.

This is the final agent in the pipeline. It synthesizes outputs from all
previous agents and produces a structured, actionable wet-lab plan.

The plan is tailored to:
  - The specific epitopes predicted (CTL vs HTL, confidence tiers)
  - The animal models available (mouse H-2 and macaque Mamu alleles from N3)
  - Population coverage profile (from N7)
  - Safety screening results (from N6)
  - Literature evidence and failure signals (from N9)
  - Client lab constraints (species, equipment, budget from request params)


All generation is done via Claude API (claude Sonnet 4.6).

References:
  PBMC ELISpot: Janetzki et al. (2015) Cancer Immunol Immunother
  MHC tetramer: Klenerman et al. (2002) Nat Rev Immunol
  BCG boost protocol: Goonetilleke et al. (2006) J Exp Med
"""

import os
import time
import logging
import requests
from typing import List, Dict, Any, Optional

from src.models.candidate import CandidateProtein, ConfidenceTier

from src.utils.logger import get_logger
logger = get_logger("tope_deep.agents.N9")  # use the correct agent name


# ── Template fallback (no Claude) ─────────────────────────────────────────────

def _template_plan(candidate: CandidateProtein, context: Dict) -> Dict:
    """
    Structured wet-lab plan without Claude.
    Based on standard computational immunology validation workflow.
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

    plan = {
        "priority_epitopes": [
            {
                "sequence":    ep.sequence,
                "hla_allele":  ep.hla_allele,
                "ic50_nm":     ep.ic50_nm,
                "confidence":  ep.confidence_tier.value,
                "rationale":   "High confidence MHC-I binder, safety-cleared",
            }
            for ep in high_ctl
        ],

        "phase_1_assays": {
            "timeline":    "Week 1-4",
            "description": "Peptide synthesis and in vitro immunogenicity screening",
            "assays": [
                {
                    "name":      "Peptide synthesis",
                    "protocol":  "Synthesize top 5 CTL epitopes as 9-11mer peptides (>95% purity, lyophilized)",
                    "vendor":    "GenScript, Bachem, or Mimotopes",
                    "cost_usd":  "~$200-400 per peptide",
                    "reference": "Standard solid-phase peptide synthesis",
                },
                {
                    "name":      "PBMC IFN-γ ELISpot",
                    "protocol":  (
                        "Stimulate PBMCs from HLA-typed healthy donors with each peptide (10μg/mL). "
                        "Readout: IFN-γ spot-forming cells per 10^6 PBMCs. "
                        "Threshold: ≥2x background AND ≥50 SFU/10^6 = positive."
                    ),
                    "hla_donors": [ep.hla_allele for ep in high_ctl if ep.hla_allele],
                    "cost_usd":  "~$500-1000 per donor per peptide panel",
                    "reference": "Janetzki et al. (2015) Cancer Immunol Immunother 64:1695-1703",
                },
            ],
            "go_criteria": "≥2 of top 5 epitopes elicit IFN-γ response in ≥2 donors",
        },

        "phase_2_assays": {
            "timeline":    "Month 1-3",
            "description": "Mouse model validation (if mouse H-2 cross-reactive epitopes exist)",
            "assays": [
                {
                    "name":      "Murine immunization",
                    "protocol":  (
                        f"Immunize C57BL/6 (H-2b) or BALB/c (H-2d) mice with construct peptide "
                        f"+ CFA/IFA adjuvant. "
                        f"Boost at day 14. Harvest splenocytes at day 21."
                    ),
                    "mouse_epitopes": [ep.sequence for ep in mouse_reactive],
                    "strains":   ["C57BL/6 (H-2b)", "BALB/c (H-2d)"],
                    "cost_usd":  "~$2000-5000 per cohort",
                    "reference": "Standard murine immunization protocol",
                    "skip_if":   len(mouse_reactive) == 0,
                },
                {
                    "name":      "Intracellular cytokine staining (ICS)",
                    "protocol":  "Stimulate splenocytes with peptide pool. Stain for CD8+/IFN-γ+/TNF-α+.",
                    "cost_usd":  "~$500-1000",
                    "reference": "Seder et al. (2008) Nat Immunol 9:239-245",
                },
            ],
            "go_criteria": "CD8+ T-cell response in ≥3/5 mice, ≥0.5% antigen-specific cells",
        },

        "phase_3_assays": {
            "timeline":    "Month 3-6",
            "description": "Advanced validation NHP model or human T-cell line confirmation",
            "assays": [
                {
                    "name":       "MHC-I tetramer staining",
                    "protocol":   "Produce HLA-A*02:01 or HLA-A*24:02 tetramers for top CTL epitopes. Stain PBMCs.",
                    "cost_usd":   "~$3000-8000 per tetramer",
                    "reference":  "Klenerman et al. (2002) Nat Rev Immunol 2:263-272",
                    "notes":      "Confirms antigen-specific T-cell frequency directly",
                },
            ],
            "go_criteria": "Tetramer+ CD8+ T cells >0.1% of total CD8+ in ≥2 donors",
        },

        "failure_risk_flags": failure_sigs + (
            ["prior experimental validation exists review known HLA restrictions before synthesis"]
            if prior_valid else []
        ) + (
            [f"only {safe_count} of {len(candidate.ctl_epitopes)} CTL epitopes safety-cleared limited synthesis options"]
            if safe_count < 3 else []
        ),

        "estimated_cost_usd": {
            "phase_1": "$3,000 - $8,000",
            "phase_2": "$5,000 - $15,000 (if mouse-reactive epitopes present)",
            "phase_3": "$10,000 - $30,000",
            "total":   "$18,000 - $53,000 for full validation",
        },

        "references": [
            "Sette & Rappuoli (2010) Immunity 33(4):530-541 reverse vaccinology framework",
            "Janetzki et al. (2015) Cancer Immunol Immunother 64:1695-1703 ELISpot harmonization",
            "Klenerman et al. (2002) Nat Rev Immunol 2:263-272 MHC tetramer",
            "Seder et al. (2008) Nat Immunol 9:239-245 T-cell quality assessment",
        ],

        "generated_by": "TOPE_DEEP template (set ANTHROPIC_API_KEY for AI-generated plan)",
    }

    return plan


# ── Claude-generated plan ─────────────────────────────────────────────────────

def _claude_plan(candidate: CandidateProtein, context: Dict) -> Optional[Dict]:
    """Generate structured wet-lab plan via Claude API."""
    api_key = os.getenv("ANTHROPIC_API_KEY", "")
    if not api_key:
        return None

    high_ctl = [
        ep for ep in candidate.ctl_epitopes
        if ep.confidence_tier == ConfidenceTier.HIGH
        and ep.allergenicity_safe is True
    ][:8]

    mouse_reactive = [
        ep.sequence for ep in candidate.ctl_epitopes
        if ep.tool_outputs.get("animal_model_alleles")
        and ep.allergenicity_safe is True
    ][:5]

    mamu_reactive = [
        ep.sequence for ep in candidate.ctl_epitopes
        if ep.tool_outputs.get("mamu_alleles")
        and ep.allergenicity_safe is True
    ][:3]

    prompt = f"""You are an expert computational immunologist at a vaccine research institute.
Generate a structured wet-lab validation plan for this vaccine candidate.

CANDIDATE: {candidate.protein_name} ({candidate.protein_id})
Sequence length: {candidate.sequence_length} aa
VaxiJen antigenicity: {candidate.vaxijen_score}

TOP CTL EPITOPES (safety-cleared, high confidence):
{chr(10).join(f"  {ep.sequence} | {ep.hla_allele} | IC50={ep.ic50_nm}nM | rank={ep.percentile_rank}" for ep in high_ctl)}

ANIMAL MODEL CROSS-REACTIVE:
  Mouse H-2: {', '.join(mouse_reactive) if mouse_reactive else 'none'}
  Macaque Mamu: {', '.join(mamu_reactive) if mamu_reactive else 'none'}

POPULATION COVERAGE:
  Global: {context.get('global_coverage_pct')}% | African: {context.get('african_coverage_pct')}%

SAFETY: {context.get('safe_count')} epitopes cleared, {context.get('unscored_count')} unscored, {context.get('fail_count')} failed

LITERATURE:
  Prior validation: {context.get('prior_validated')}
  Failure signals: {', '.join(context.get('failure_signals', [])) or 'none'}
  Summary: {context.get('literature_summary', 'not available')}

LAB CONSTRAINTS: {context.get('lab_constraints', 'standard academic lab')}

Generate a JSON response with this exact structure:
{{
  "priority_reasoning": "2-3 sentences explaining why these epitopes are prioritized",
  "phase_1": {{
    "timeline": "Week 1-4",
    "key_assay": "PBMC IFN-γ ELISpot",
    "protocol_notes": "specific protocol adaptations for this candidate",
    "expected_cost_usd": "range",
    "go_criteria": "measurable threshold"
  }},
  "phase_2": {{
    "timeline": "Month 1-3",
    "key_assay": "murine immunization" or "skip no H-2 cross-reactive epitopes",
    "protocol_notes": "specific notes",
    "expected_cost_usd": "range",
    "go_criteria": "measurable threshold"
  }},
  "phase_3": {{
    "timeline": "Month 3-6",
    "key_assay": "MHC tetramer or NHP model",
    "protocol_notes": "specific notes",
    "expected_cost_usd": "range",
    "go_criteria": "measurable threshold"
  }},
  "critical_risks": ["list of specific risks for this candidate"],
  "regulatory_note": "one sentence on regulatory pathway implications"
}}
Return only valid JSON, no markdown."""

    try:
        resp = requests.post(
            "https://api.anthropic.com/v1/messages",
            headers={
                "x-api-key": api_key,
                "anthropic-version": "2023-06-01",
                "content-type": "application/json",
            },
            json={
                "model": "claude-sonnet-4-6",
                "max_tokens": 1000,
                "messages": [{"role": "user", "content": prompt}],
            },
            timeout=30,
        )
        resp.raise_for_status()
        import json
        text = resp.json()["content"][0]["text"].strip()
        parsed = json.loads(text)
        parsed["generated_by"] = "Claude claude-sonnet-4-6 via TOPE_DEEP N10"
        return parsed
    except Exception as e:
        logger.warning(f"N10: Claude plan generation failed: {e}")
        return None


class ExperimentPlannerAgent:
    """
    Experiment Planner Agent
    Generates a prioritized wet-lab validation roadmap per candidate.
    Uses Claude API when available, falls back to template-based plan.
    """

    def __init__(self):
        self.stage_name = "experiment_planning"

    def run(
        self,
        candidates: List[CandidateProtein],
        lab_constraints: str = "standard academic lab",
    ) -> List[CandidateProtein]:
        logger.info("N10: Starting experiment planning")
        active = [c for c in candidates if c.status.value == "active"]
        has_claude = bool(os.getenv("ANTHROPIC_API_KEY", ""))
        logger.info(f"   {len(active)} candidates | Claude={'yes' if has_claude else 'no (template fallback)'}")

        for i, candidate in enumerate(active):
            logger.info(f"   [{i+1}/{len(active)}] {candidate.protein_name}")
            start = time.time()

            # Gather context from previous agents
            safety_decision = next(
                (d for d in candidate.decisions if d.get("stage") == "safety_filter"), {}
            )
            lit_decision = next(
                (d for d in candidate.decisions if d.get("stage") == "literature_search"), {}
            )

            context = {
                "safe_count":          safety_decision.get("safe_count", 0),
                "unscored_count":      safety_decision.get("unscored_count", 0),
                "fail_count":          safety_decision.get("fail_count", 0),
                "global_coverage_pct": round((candidate.hla_coverage_global or 0) * 100, 1),
                "african_coverage_pct":round((candidate.hla_coverage_africa or 0) * 100, 1),
                "prior_validated":     lit_decision.get("prior_validated", False),
                "failure_signals":     lit_decision.get("failure_signals", []),
                "literature_summary":  lit_decision.get("literature_summary"),
                "lab_constraints":     lab_constraints,
            }

            # Try Claude first, fall back to template
            plan = _claude_plan(candidate, context)
            if plan is None:
                plan = _template_plan(candidate, context)
                generation_method = "template"
            else:
                generation_method = "claude-sonnet-4-6"

            elapsed = time.time() - start

            candidate.add_decision(
                stage=self.stage_name,
                decision="plan_generated",
                reasoning=(
                    f"Wet-lab validation plan generated for {candidate.protein_name}. "
                    f"Method: {generation_method}. "
                    f"Priority epitopes: {len(plan.get('priority_epitopes', plan.get('phase_1', {})))}. "
                    f"Lab constraints: {lab_constraints}. "
                    f"Generation time: {elapsed:.1f}s."
                ),
                experiment_plan=plan,
                generation_method=generation_method,
                lab_constraints=lab_constraints,
                plan_time_s=round(elapsed, 1),
            )

            # Store plan on candidate for easy access
            candidate.experiment_plan = plan

            logger.info(
                f"      Plan generated via {generation_method} in {elapsed:.1f}s"
            )

        logger.info("N10: Experiment planning complete")
        return candidates


experiment_planner = ExperimentPlannerAgent()