"""
KOZI AI - Vaccine Discovery Platform
Streamlit interface combining MVP-1 (data acquisition) and MVP-2 (analysis pipeline).

Run: streamlit run streamlit_app.py

Flow:
  1. Researcher enters pathogen name, UniProt ID, or raw sequence
  2. MVP-1: Fetch proteins from UniProt, filter by surface localization
  3. MVP-2: N3 (T-cell) → N4 (B-cell) → N6 (Safety) → N7 (Coverage)
  4. Display ranked results with charts, tables, CSV download
"""

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import requests
import time
import json
import logging
from datetime import datetime
from typing import List, Dict, Any, Optional

# Pipeline imports
from src.models.candidate import (
    CandidateProtein, CandidateStatus, EpitopeResult,
    EpitopeType, ConfidenceTier,
)
from src.agents.predictors.tcell_predictor import TCellPredictorAgent
from src.agents.predictors.bcell_predictor import BCellPredictorAgent
from src.agents.predictors.safety_filter import SafetyFilterAgent
from src.agents.predictors.coverage_agent import CoverageAgent

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("kozi.streamlit")

# --PAGE CONFIG ─────────────────────────────────────────────────────────

st.set_page_config(
    page_title="Kozi AI - Vaccine Discovery",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --CUSTOM CSS ──────────────────────────────────────────────────────────

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');



/* Global background */
.stApp {
    background-color: #000000;
}

.block-container {
    max-width: 1050px;
    margin: auto;
    padding-top: 2rem;
}
hr {
    border: none;
    border-top: 1px solid #e5e5e5;
}
/* Header */
.kozi-header {
    text-align: center;
    padding: 1.5rem 0 1rem;
}

.kozi-header h1 {
    font-size: 2rem;
    font-weight: 700;
    color: #ffffff;
    margin-bottom: 0.3rem;
}

.kozi-header p {
    color: rgba(f,f,f,0.8);
    font-size: 1rem;
    margin: 0;
}

/* Buttons */
.stButton>button {
    background-color: #000000;
    color: #ffffff;
    border-radius: 6px;
    border: none;
    padding: 10px 18px;
    font-weight: 500;
}

.stButton>button:hover {
    background-color: #000000;
    color: #ffffff;
    opacity: 0.9;
}

/* Metric cards */
.metric-row {
    display: flex;
    gap: 12px;
    margin: 1rem 0;
}



.metric-card .value {
    font-size: 1.8rem;
    font-weight: 700;
    color: #000000;
    line-height: 1.2;
}

.metric-card .label {
    font-size: 0.85rem;
    color: rgba(0,0,0,0.8);
    margin-top: 2px;
}

.metric-card .sub {
    font-size: 0.75rem;
    color: rgba(0,0,0,0.6);
}

/* Pipeline steps */
.step-indicator {
    display: flex;
    gap: 8px;
    margin: 1rem 0;
    align-items: center;
}

.step {
    flex: 1;
    text-align: center;
    padding: 8px 12px;
    border-radius: 6px;
    font-size: 0.8rem;
    font-weight: 500;
    background: #ffffff;
    border: 1px solid #e5e5e5;
    color: #000000;
}

.step-arrow {
    color: #000000;
    font-size: 1.2rem;
}

.metric-card {
    flex: 1;
    background: #ffffff;
    border-radius: 8px;
    padding: 16px 20px;
    border: 1px solid #e5e5e5;
    transition: all 0.15s ease;
}

.metric-card:hover {
    transform: translateY(-2px);
    border-color: #000000;
}




.epitope-card:hover {
    border-color: #000000;
}

.epitope-card {
    background: #ffffff;
    border-radius: 6px;
    padding: 10px 14px;
    margin: 6px 0;
    border: 1px solid #e5e5e5;
    transition: all 0.15s ease;
}

.epitope-card .seq {
    font-family: 'Courier New', monospace;
    font-size: 0.95rem;
    font-weight: 600;
    color: #000000;
    letter-spacing: 1px;
}

.epitope-card .meta {
    font-size: 0.78rem;
    color: rgba(0,0,0,0.7);
    margin-top: 3px;
}

/* Tables */
.stDataFrame {
    border: 1px solid #e5e5e5;
}

/* Expanders */
.streamlit-expanderHeader {
    color: #000000;
}

/* Download buttons */
.stDownloadButton>button {
    background-color: #000000;
    color: #ffffff;
    border-radius: 6px;
    border: none;
    padding: 10px 18px;
}

</style>
""", unsafe_allow_html=True)


# --UNIPROT FUNCTIONS (MVP-1 simplified) ────────────────────────────────

def fetch_protein_by_id(uniprot_id: str) -> Optional[Dict]:
    """Fetch a single protein from UniProt by ID."""
    try:
        # Get FASTA sequence
        resp = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta",
            timeout=15,
        )
        resp.raise_for_status()
        lines = resp.text.strip().split("\n")
        header = lines[0][1:]  # Remove >
        sequence = "".join(l.strip() for l in lines[1:])

        # Get metadata
        resp2 = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json",
            timeout=15,
        )
        metadata = {}
        if resp2.status_code == 200:
            data = resp2.json()
            metadata = {
                "organism": data.get("organism", {}).get("scientificName", "Unknown"),
                "protein_name": data.get("proteinDescription", {}).get(
                    "recommendedName", {}).get("fullName", {}).get("value", header.split()[0]),
            }

        return {
            "protein_id": uniprot_id,
            "protein_name": metadata.get("protein_name", header.split()[0]),
            "organism": metadata.get("organism", "Unknown"),
            "sequence": sequence,
            "length": len(sequence),
        }
    except Exception as e:
        logger.error(f"UniProt fetch failed for {uniprot_id}: {e}")
        return None


def search_pathogen_proteins(pathogen_name: str, max_results: int = 10) -> List[Dict]:
    """
    Search UniProt for surface proteins of a pathogen.
    This is a simplified MVP-1 N1+N2: fetch + basic surface filter.
    """
    try:
        # Search UniProt for reviewed proteins from this organism
        query = f'(organism_name:"{pathogen_name}") AND (reviewed:true)'
        params = {
            "query": query,
            "format": "json",
            "size": min(max_results * 3, 50),  # Fetch extra to filter
            "fields": "accession,protein_name,organism_name,length,cc_subcellular_location",
        }

        resp = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params=params,
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()

        proteins = []
        for entry in data.get("results", []):
            accession = entry.get("primaryAccession", "")
            name_obj = entry.get("proteinDescription", {})
            rec_name = name_obj.get("recommendedName", {})
            protein_name = rec_name.get("fullName", {}).get("value", "Unknown")

            # Check subcellular location for surface/secreted proteins
            location = ""
            for comment in entry.get("comments", []):
                if comment.get("commentType") == "SUBCELLULAR LOCATION":
                    for loc in comment.get("subcellularLocations", []):
                        loc_val = loc.get("location", {}).get("value", "")
                        location += loc_val + " "

            # Basic surface filter (MVP-1 N2 simplified)
            is_surface = any(kw in location.lower() for kw in [
                "membrane", "secreted", "cell surface", "extracellular",
                "outer membrane", "cell wall", "exported",
            ])

            organism = entry.get("organism", {}).get(
                "scientificName", "Unknown")
            length = entry.get("sequence", {}).get("length", 0)

            proteins.append({
                "protein_id": accession,
                "protein_name": protein_name,
                "organism": organism,
                "length": length,
                "location": location.strip() or "Unknown",
                "is_surface": is_surface,
            })

        # Sort: surface proteins first, then by length (larger often more antigenic)
        proteins.sort(key=lambda x: (not x["is_surface"], -x["length"]))
        return proteins[:max_results]

    except Exception as e:
        logger.error(f"UniProt search failed: {e}")
        return []


def fetch_sequence(protein_id: str) -> Optional[str]:
    """Fetch just the sequence for a protein ID."""
    try:
        resp = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta",
            timeout=15,
        )
        resp.raise_for_status()
        lines = resp.text.strip().split("\n")
        return "".join(l.strip() for l in lines[1:])
    except Exception:
        return None


# --PIPELINE RUNNER ─────────────────────────────────────────────────────

def run_pipeline_on_candidate(candidate: CandidateProtein, progress_callback=None) -> CandidateProtein:
    """Run MVP-2 pipeline (N3→N4→N6→N7) on a single candidate."""

    # N3: T-cell
    if progress_callback:
        progress_callback("N3: Predicting T-cell epitopes...", 0.2)
    n3 = TCellPredictorAgent()
    candidates = n3.run([candidate])
    candidate = candidates[0]

    # N4: B-cell
    if progress_callback:
        progress_callback("N4: Predicting B-cell epitopes...", 0.45)
    n4 = BCellPredictorAgent()
    candidates = n4.run([candidate])
    candidate = candidates[0]

    # N6: Safety
    if progress_callback:
        progress_callback("N6: Safety screening...", 0.65)
    n6 = SafetyFilterAgent()
    candidates = n6.run([candidate])
    candidate = candidates[0]

    # N7: Coverage
    if progress_callback:
        progress_callback("N7: Calculating population coverage...", 0.85)
    n7 = CoverageAgent()
    candidates = n7.run([candidate])
    candidate = candidates[0]

    if progress_callback:
        progress_callback("Complete!", 1.0)

    return candidate


# --UI COMPONENTS ───────────────────────────────────────────────────────

def render_header():
    st.markdown("""
    <div class="kozi-header">
        <h1>🧬 Kozi AI</h1>
        <p>Automated vaccine target discovery - from pathogen to epitope candidates in minutes</p>
    </div>
    """, unsafe_allow_html=True)


def render_metrics(candidate: CandidateProtein):
    ctl_strong = len(
        [e for e in candidate.ctl_epitopes if e.confidence_tier == ConfidenceTier.HIGH])
    htl_count = len(candidate.htl_epitopes)
    bcell_count = len(candidate.bcell_epitopes)
    global_cov = round((candidate.hla_coverage_global or 0) * 100, 1)
    african_cov = round((candidate.hla_coverage_africa or 0) * 100, 1)

    st.markdown(f"""
    <div class="metric-row">
        <div class="metric-card purple">
            <div class="value">{len(candidate.ctl_epitopes)}</div>
            <div class="label">CTL epitopes</div>
            <div class="sub">{ctl_strong} strong binders</div>
        </div>
        <div class="metric-card blue">
            <div class="value">{htl_count}</div>
            <div class="label">HTL epitopes</div>
            <div class="sub">Helper T-cell targets</div>
        </div>
        <div class="metric-card">
            <div class="value">{bcell_count}</div>
            <div class="label">B-cell epitopes</div>
            <div class="sub">Antibody targets</div>
        </div>
        <div class="metric-card green">
            <div class="value">{global_cov}%</div>
            <div class="label">Global coverage</div>
            <div class="sub">MHC-I + MHC-II combined</div>
        </div>
        <div class="metric-card amber">
            <div class="value">{african_cov}%</div>
            <div class="label">African coverage</div>
            <div class="sub">Priority population</div>
        </div>
    </div>
    """, unsafe_allow_html=True)


def render_coverage_chart(candidate: CandidateProtein):
    """Population coverage bar chart from decision audit trail."""
    # Find coverage decision
    cov_decision = None
    for d in candidate.decisions:
        if d.get("stage") == "coverage_analysis":
            cov_decision = d
            break

    if not cov_decision or "per_population" not in cov_decision:
        st.info("No coverage data available")
        return

    pop_data = cov_decision["per_population"]
    pops = []
    mhc_i = []
    mhc_ii = []
    combined = []

    for key, val in pop_data.items():
        pops.append(val.get("population_label", key))
        mhc_i.append(val.get("mhc_i_pct", 0))
        mhc_ii.append(val.get("mhc_ii_pct", 0))
        combined.append(val.get("combined_pct", 0))

    fig = go.Figure()
    fig.add_trace(go.Bar(name="MHC-I (killer T-cells)", x=pops, y=mhc_i,
                         marker_color="#7C3AED", opacity=0.85))
    fig.add_trace(go.Bar(name="MHC-II (helper T-cells)", x=pops, y=mhc_ii,
                         marker_color="#0D9488", opacity=0.85))
    fig.add_trace(go.Scatter(name="Combined", x=pops, y=combined,
                             mode="lines+markers", line=dict(color="#1E293B", width=2),
                             marker=dict(size=8)))

    fig.update_layout(
        barmode="group",
        title="Population coverage by immune response type",
        yaxis_title="Coverage (%)",
        yaxis=dict(range=[0, 100]),
        height=400,
        template="plotly_white",
        font=dict(family="Inter, sans-serif"),
        legend=dict(orientation="h", yanchor="bottom",
                    y=1.02, xanchor="right", x=1),
    )

    fig.add_hline(y=80, line_dash="dash", line_color="#DC2626", opacity=0.5,
                  annotation_text="80% target", annotation_position="top right")

    st.plotly_chart(fig, use_container_width=True)


def render_epitope_table(candidate: CandidateProtein):
    """Epitope results as a sortable table."""
    rows = []

    for ep in candidate.ctl_epitopes:
        rows.append({
            "Sequence": ep.sequence,
            "Type": "CTL",
            "HLA allele": ep.hla_allele or "",
            "IC50 (nM)": round(ep.ic50_nm, 1) if ep.ic50_nm else "",
            "Rank": round(ep.percentile_rank, 2) if ep.percentile_rank else "",
            "Confidence": ep.confidence_tier.value,
            "Safe": "Yes" if ep.allergenicity_safe else ("Flagged" if ep.tool_outputs.get("safety_flags") else "-"),
        })

    for ep in candidate.htl_epitopes:
        rows.append({
            "Sequence": ep.sequence,
            "Type": "HTL",
            "HLA allele": ep.hla_allele or "",
            "IC50 (nM)": round(ep.ic50_nm, 1) if ep.ic50_nm else "",
            "Rank": round(ep.percentile_rank, 2) if ep.percentile_rank else "",
            "Confidence": ep.confidence_tier.value,
            "Safe": "Yes" if ep.allergenicity_safe else ("Flagged" if ep.tool_outputs.get("safety_flags") else "-"),
        })

    for ep in candidate.bcell_epitopes:
        rows.append({
            "Sequence": ep.sequence,
            "Type": "B-cell",
            "HLA allele": "N/A",
            "IC50 (nM)": "",
            "Rank": "",
            "Confidence": ep.confidence_tier.value,
            "Safe": "Yes" if ep.allergenicity_safe else ("Flagged" if ep.tool_outputs.get("safety_flags") else "-"),
        })

    if rows:
        df = pd.DataFrame(rows)
        st.dataframe(df, use_container_width=True, hide_index=True)
        return df
    else:
        st.info("No epitopes found")
        return None


def render_top_epitopes(candidate: CandidateProtein):
    """Show top CTL epitopes as styled cards."""
    top_ctl = sorted(candidate.ctl_epitopes,
                     key=lambda e: e.ic50_nm or 99999)[:6]

    if not top_ctl:
        return

    cols = st.columns(3)
    for i, ep in enumerate(top_ctl):
        conf_class = ep.confidence_tier.value
        with cols[i % 3]:
            st.markdown(f"""
            <div class="epitope-card {conf_class}">
                <div class="seq">{ep.sequence}</div>
                <div class="meta">{ep.hla_allele} &nbsp;·&nbsp; IC50 {ep.ic50_nm:.0f} nM &nbsp;·&nbsp; rank {ep.percentile_rank:.2f}</div>
            </div>
            """, unsafe_allow_html=True)


def build_csv_export(candidate: CandidateProtein) -> str:
    """Build CSV string for download."""
    rows = []
    for ep in candidate.ctl_epitopes:
        rows.append(f"{candidate.protein_name},{candidate.protein_id},{ep.sequence},CTL,{ep.hla_allele},{ep.ic50_nm},{ep.confidence_tier.value},{ep.allergenicity_safe},{ep.toxicity_safe}")
    for ep in candidate.htl_epitopes:
        rows.append(f"{candidate.protein_name},{candidate.protein_id},{ep.sequence},HTL,{ep.hla_allele},{ep.ic50_nm},{ep.confidence_tier.value},{ep.allergenicity_safe},{ep.toxicity_safe}")
    for ep in candidate.bcell_epitopes:
        rows.append(
            f"{candidate.protein_name},{candidate.protein_id},{ep.sequence},B-cell,N/A,,{ep.confidence_tier.value},{ep.allergenicity_safe},{ep.toxicity_safe}")

    header = "protein,protein_id,epitope_sequence,type,hla_allele,ic50_nm,confidence,allergenicity_safe,toxicity_safe"
    return header + "\n" + "\n".join(rows)


# --MAIN APP ────────────────────────────────────────────────────────────

def main():
    render_header()

    # --SIDEBAR ──
    with st.sidebar:
        st.markdown("### Discovery options")

        input_mode = st.radio(
            "Input method",
            ["Pathogen name", "UniProt protein ID", "Paste sequence"],
            index=0,
        )

        if input_mode == "Pathogen name":
            pathogen = st.text_input(
                "Pathogen name",
                value="Mycobacterium tuberculosis",
                help="Scientific name of the pathogen",
            )
            max_proteins = st.slider("Max proteins to analyze", 1, 10, 3)

        elif input_mode == "UniProt protein ID":
            protein_id = st.text_input(
                "UniProt ID",
                value="P9WNK7",
                help="Example: P9WNK7 (M. tuberculosis ESAT-6)",
            )

        else:
            protein_name = st.text_input(
                "Protein name", value="Custom protein")
            raw_sequence = st.text_area(
                "Amino acid sequence",
                value="MTEQQWNFAGIEAAASAIQGNVTSIHSLLDEGKQSLTKLAAAWGGSGSEAYQGVQQKWDATATELNNALQNLARTISEAGQAMASTEGNVTGMFA",
                height=120,
            )

        st.markdown("---")
        st.markdown("### Pipeline settings")
        run_safety = st.checkbox("Run safety screening (N6)", value=True,
                                 help="AllerTOP + AllergenFP + ToxinPred - adds ~2 min per protein")
        run_coverage = st.checkbox("Run population coverage (N7)", value=True)

        st.markdown("---")
        run_btn = st.button("Run discovery pipeline",
                            type="primary", use_container_width=True)

    # --MAIN CONTENT ──
    if not run_btn:
        # Welcome screen
        st.markdown("---")
        col1, col2, col3 = st.columns(3)

        with col1:
            st.markdown("#### T-cell prediction")
            st.markdown(
                "Find peptides that bind HLA molecules to activate killer and helper T-cells. Uses IEDB NetMHCpan across 15 HLA-I and 7 HLA-II alleles.")

        with col2:
            st.markdown("#### Safety screening")
            st.markdown(
                "Screen every epitope for allergenicity (AllerTOP + AllergenFP) and toxicity (ToxinPred). Dual redundancy on allergenicity.")

        with col3:
            st.markdown("#### Population coverage")
            st.markdown(
                "Calculate what percentage of 7 global populations would be protected. Built-in equity analysis for African populations.")

        st.markdown("---")
        st.markdown("##### How it works")
        st.markdown(
            "Enter a pathogen name → Kozi fetches proteins from UniProt → filters for surface antigens → "
            "predicts T-cell and B-cell epitopes → screens for safety → calculates global population coverage → "
            "exports ranked results as CSV for your lab team."
        )

        st.markdown("---")
        st.caption(
            "Kozi AI - umukozi means 'worker' in Kinyarwanda. We build agents that work so scientists can discover.")
        return

    # --RUN PIPELINE ──
    st.markdown("---")

    # Step 1: Get protein(s) based on input mode
    candidates_to_run = []

    if input_mode == "Pathogen name":
        with st.status(f"Searching UniProt for {pathogen} proteins...", expanded=True) as status:
            st.write(
                f"Querying UniProt for reviewed proteins from *{pathogen}*...")
            proteins = search_pathogen_proteins(
                pathogen, max_results=max_proteins)

            if not proteins:
                st.error(
                    f"No proteins found for '{pathogen}'. Try the exact scientific name.")
                return

            st.write(f"Found {len(proteins)} proteins. Fetching sequences...")

            for p in proteins:
                seq = fetch_sequence(p["protein_id"])
                if seq and len(seq) >= 20:
                    candidates_to_run.append(CandidateProtein(
                        protein_id=p["protein_id"],
                        protein_name=f"{p['protein_name']} ({p['protein_id']})",
                        sequence=seq,
                        source="uniprot",
                        stage="antigen_screening",
                        status=CandidateStatus.ACTIVE,
                    ))

            status.update(
                label=f"Found {len(candidates_to_run)} proteins to analyze", state="complete")

    elif input_mode == "UniProt protein ID":
        with st.spinner(f"Fetching {protein_id} from UniProt..."):
            prot = fetch_protein_by_id(protein_id)
            if not prot:
                st.error(f"Could not fetch {protein_id}. Check the ID.")
                return
            candidates_to_run.append(CandidateProtein(
                protein_id=prot["protein_id"],
                protein_name=f"{prot['protein_name']} ({prot['protein_id']})",
                sequence=prot["sequence"],
                source="uniprot",
                stage="antigen_screening",
                status=CandidateStatus.ACTIVE,
            ))

    else:
        seq = raw_sequence.upper().replace(" ", "").replace("\n", "")
        if len(seq) < 10:
            st.error("Sequence too short (minimum 10 amino acids)")
            return
        candidates_to_run.append(CandidateProtein(
            protein_id="user_input",
            protein_name=protein_name,
            sequence=seq,
            source="user_input",
            stage="antigen_screening",
            status=CandidateStatus.ACTIVE,
        ))

    if not candidates_to_run:
        st.error("No valid proteins to analyze")
        return

    # Step 2: Run pipeline on each candidate
    results = []
    progress_bar = st.progress(0)
    status_text = st.empty()

    total = len(candidates_to_run)
    for idx, candidate in enumerate(candidates_to_run):
        st.markdown(
            f"#### Analyzing: {candidate.protein_name} ({len(candidate.sequence)} aa)")

        def update_progress(msg, pct):
            overall = (idx + pct) / total
            progress_bar.progress(min(overall, 1.0))
            status_text.text(f"[{idx+1}/{total}] {msg}")

        try:
            result = run_pipeline_on_candidate(
                candidate, progress_callback=update_progress)
            results.append(result)
        except Exception as e:
            st.error(f"Pipeline failed for {candidate.protein_name}: {e}")
            logger.error(f"Pipeline error: {e}", exc_info=True)

    progress_bar.progress(1.0)
    status_text.text("Pipeline complete!")
    time.sleep(0.5)
    status_text.empty()
    progress_bar.empty()

    # Step 3: Display results
    if not results:
        st.error("No results generated")
        return

    st.markdown("---")
    st.markdown("## Results")

    # If multiple proteins, show comparison first
    if len(results) > 1:
        st.markdown("### Protein ranking")
        ranking_data = []
        for r in results:
            ranking_data.append({
                "Protein": r.protein_name,
                "CTL": len(r.ctl_epitopes),
                "HTL": len(r.htl_epitopes),
                "B-cell": len(r.bcell_epitopes),
                "Total": r.get_total_epitopes(),
                "Global %": round((r.hla_coverage_global or 0) * 100, 1),
                "African %": round((r.hla_coverage_africa or 0) * 100, 1),
            })
        ranking_df = pd.DataFrame(ranking_data)
        ranking_df = ranking_df.sort_values("Total", ascending=False)
        st.dataframe(ranking_df, use_container_width=True, hide_index=True)
        st.markdown("---")

    # Show detailed results per protein
    for result in results:
        if len(results) > 1:
            st.markdown(f"### {result.protein_name}")

        # Metrics row
        render_metrics(result)

        # Top epitopes
        st.markdown("#### Top CTL epitopes")
        render_top_epitopes(result)

        # Coverage chart
        if run_coverage:
            st.markdown("#### Population coverage")
            render_coverage_chart(result)

        # Full epitope table
        with st.expander("View all epitopes", expanded=False):
            df = render_epitope_table(result)

        # Decision audit trail
        with st.expander("Decision audit trail", expanded=False):
            for d in result.decisions:
                st.markdown(f"**{d['stage']}**: {d['reasoning']}")

        st.markdown("---")

    # Step 4: Export
    st.markdown("### Export results")

    col1, col2, col3 = st.columns(3)

    with col1:
        # CSV download
        all_csv_parts = []
        for r in results:
            all_csv_parts.append(build_csv_export(r))

        # Combine (skip duplicate headers)
        header = all_csv_parts[0].split("\n")[0] if all_csv_parts else ""
        csv_body = "\n".join(
            "\n".join(part.split("\n")[1:]) for part in all_csv_parts if part
        )
        full_csv = header + "\n" + csv_body

        st.download_button(
            "Download epitopes (CSV)",
            data=full_csv,
            file_name=f"kozi_epitopes_{datetime.now().strftime('%Y%m%d_%H%M')}.csv",
            mime="text/csv",
            use_container_width=True,
        )

    with col2:
        # JSON download
        json_data = {
            "pipeline": "Kozi AI MVP-2",
            "timestamp": datetime.now().isoformat(),
            "candidates": [
                {
                    "protein_id": r.protein_id,
                    "protein_name": r.protein_name,
                    "sequence_length": len(r.sequence),
                    "ctl_epitopes": len(r.ctl_epitopes),
                    "htl_epitopes": len(r.htl_epitopes),
                    "bcell_epitopes": len(r.bcell_epitopes),
                    "global_coverage": round((r.hla_coverage_global or 0) * 100, 1),
                    "african_coverage": round((r.hla_coverage_africa or 0) * 100, 1),
                    "decisions": r.decisions,
                }
                for r in results
            ],
        }
        st.download_button(
            "Download audit trail (JSON)",
            data=json.dumps(json_data, indent=2, default=str),
            file_name=f"kozi_results_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
            mime="application/json",
            use_container_width=True,
        )

    with col3:
        total_epitopes = sum(r.get_total_epitopes() for r in results)
        st.info(f"{total_epitopes} epitopes across {len(results)} protein(s)")

    st.caption(
        "Kozi AI - 3 minutes instead of 3 weeks. 6 APIs. 7 populations. Full audit trail.")


if __name__ == "__main__":
    main()
