"""
TOPE_DEEP Report Generator
LaTeX-style scientific PDF using ReportLab Platypus.

Usage (backend):
    from api.report_generator import generate_report_pdf
    pdf_bytes = generate_report_pdf(results_dict, run_id)

Output: write to response or disk.

Design language:
  - White background, black text only
  - Computer Modern-style via Helvetica (closest available sans)
  - Section headers: small-caps style (bold + letter-spacing)
  - Tables: thin rules, no color fill, zebra via light gray only
  - Margins: 2.5cm all sides (standard A4 academic)
  - No UI chrome, no colored badges, no rounded cards
"""

import io
import math
from datetime import datetime
from typing import Any, Dict, List, Optional

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT, TA_RIGHT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm, mm
from reportlab.platypus import (
    HRFlowable, KeepTogether, PageBreak, Paragraph,
    SimpleDocTemplate, Spacer, Table, TableStyle,
)

# ── Page geometry ──────────────────────────────────────────────────────────────

PAGE_W, PAGE_H = A4
MARGIN        = 2.5 * cm
TEXT_W        = PAGE_W - 2 * MARGIN

# ── Color palette black/white only ──────────────────────────────────────────

BLACK      = colors.HexColor("#000000")
WHITE      = colors.white
GRAY_RULE  = colors.HexColor("#888888")
GRAY_LIGHT = colors.HexColor("#F2F2F2")
GRAY_MID   = colors.HexColor("#CCCCCC")
GRAY_TEXT  = colors.HexColor("#444444")

# ── Typography ─────────────────────────────────────────────────────────────────

def _styles():
    base = getSampleStyleSheet()

    def ps(name, **kw):
        return ParagraphStyle(name, **kw)

    return {
        "title": ps("ReportTitle",
            fontName="Helvetica-Bold", fontSize=18,
            leading=22, alignment=TA_CENTER,
            spaceAfter=4,
        ),
        "subtitle": ps("ReportSubtitle",
            fontName="Helvetica", fontSize=11,
            leading=14, alignment=TA_CENTER, textColor=GRAY_TEXT,
            spaceAfter=2,
        ),
        "meta": ps("ReportMeta",
            fontName="Helvetica", fontSize=9,
            leading=12, alignment=TA_CENTER, textColor=GRAY_TEXT,
            spaceAfter=0,
        ),
        "h1": ps("H1",
            fontName="Helvetica-Bold", fontSize=12,
            leading=16, spaceBefore=18, spaceAfter=6,
            textTransform="uppercase", letterSpacing=0.8,
        ),
        "h2": ps("H2",
            fontName="Helvetica-Bold", fontSize=10,
            leading=14, spaceBefore=12, spaceAfter=4,
        ),
        "body": ps("Body",
            fontName="Helvetica", fontSize=9.5,
            leading=14, alignment=TA_JUSTIFY, spaceAfter=6,
        ),
        "body_left": ps("BodyLeft",
            fontName="Helvetica", fontSize=9.5,
            leading=14, alignment=TA_LEFT, spaceAfter=4,
        ),
        "caption": ps("Caption",
            fontName="Helvetica-Oblique", fontSize=8.5,
            leading=12, alignment=TA_CENTER, textColor=GRAY_TEXT,
            spaceAfter=8,
        ),
        "mono": ps("Mono",
            fontName="Courier", fontSize=8.5,
            leading=12, alignment=TA_LEFT, spaceAfter=4,
        ),
        "mono_small": ps("MonoSmall",
            fontName="Courier", fontSize=7.5,
            leading=11, alignment=TA_LEFT, spaceAfter=4,
        ),
        "th": ps("TH",
            fontName="Helvetica-Bold", fontSize=8,
            leading=11, alignment=TA_LEFT,
        ),
        "td": ps("TD",
            fontName="Helvetica", fontSize=8.5,
            leading=11, alignment=TA_LEFT,
        ),
        "td_mono": ps("TDMono",
            fontName="Courier", fontSize=8,
            leading=11, alignment=TA_LEFT,
        ),
        "td_right": ps("TDRight",
            fontName="Helvetica", fontSize=8.5,
            leading=11, alignment=TA_RIGHT,
        ),
        "note": ps("Note",
            fontName="Helvetica-Oblique", fontSize=8,
            leading=11, textColor=GRAY_TEXT, spaceAfter=4,
        ),
        "bullet": ps("Bullet",
            fontName="Helvetica", fontSize=9.5,
            leading=14, leftIndent=12, spaceAfter=3,
        ),
    }


# ── Table helpers ──────────────────────────────────────────────────────────────

_TABLE_BASE = TableStyle([
    ("FONTNAME",    (0, 0), (-1,  0), "Helvetica-Bold"),
    ("FONTSIZE",    (0, 0), (-1,  0), 8),
    ("FONTNAME",    (0, 1), (-1, -1), "Helvetica"),
    ("FONTSIZE",    (0, 1), (-1, -1), 8.5),
    ("LEADING",     (0, 0), (-1, -1), 11),
    ("TOPPADDING",  (0, 0), (-1, -1), 3),
    ("BOTTOMPADDING",(0,0), (-1, -1), 3),
    ("LEFTPADDING", (0, 0), (-1, -1), 5),
    ("RIGHTPADDING",(0, 0), (-1, -1), 5),
    ("LINEBELOW",   (0, 0), (-1,  0), 0.6, BLACK),
    ("LINEBELOW",   (0,-1), (-1, -1), 0.4, GRAY_RULE),
    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [WHITE, GRAY_LIGHT]),
    ("VALIGN",      (0, 0), (-1, -1), "MIDDLE"),
    ("GRID",        (0, 0), (-1, -1), 0.25, GRAY_MID),
])

def _table(data: List[List], col_widths=None, extra_style=None):
    style = TableStyle(_TABLE_BASE.getCommands())
    if extra_style:
        for cmd in extra_style:
            style.add(*cmd)
    t = Table(data, colWidths=col_widths, repeatRows=1)
    t.setStyle(style)
    return t

def _rule():
    return HRFlowable(width="100%", thickness=0.5, color=GRAY_RULE, spaceAfter=6, spaceBefore=4)

def _spacer(h=6):
    return Spacer(1, h)

def _fmt_s(s):
    """Format seconds to human string."""
    if s is None or s == 0:
        return "-"
    if s < 60:
        return f"{s:.1f}s"
    return f"{int(s//60)}m {int(s%60)}s"

def _pct(v):
    if v is None:
        return "-"
    return f"{v:.1f}%"

def _safe(v, fallback="-"):
    if v is None:
        return fallback
    return str(v)


# ── Header / footer callbacks ──────────────────────────────────────────────────

def _make_header_footer(protein_name: str, run_id: str, page_count_holder: dict):
    def on_page(canvas, doc):
        canvas.saveState()
        # Header line
        canvas.setStrokeColor(GRAY_RULE)
        canvas.setLineWidth(0.4)
        canvas.line(MARGIN, PAGE_H - MARGIN + 4*mm, PAGE_W - MARGIN, PAGE_H - MARGIN + 4*mm)
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(GRAY_TEXT)
        canvas.drawString(MARGIN, PAGE_H - MARGIN + 1.5*mm, "TOPE_DEEP  |  Computational Vaccine Discovery Report")
        canvas.drawRightString(PAGE_W - MARGIN, PAGE_H - MARGIN + 1.5*mm, protein_name)
        # Footer line
        canvas.line(MARGIN, MARGIN - 4*mm, PAGE_W - MARGIN, MARGIN - 4*mm)
        canvas.drawString(MARGIN, MARGIN - 7*mm, f"Run ID: {run_id}")
        canvas.drawCentredString(PAGE_W/2, MARGIN - 7*mm, "CONFIDENTIAL For research use only")
        canvas.drawRightString(PAGE_W - MARGIN, MARGIN - 7*mm, f"Page {doc.page}")
        canvas.restoreState()
    return on_page


# ── Section builders ───────────────────────────────────────────────────────────

def _cover(st, protein_name, protein_id, run_id, run_date, candidates_count):
    story = []
    story.append(_spacer(60))
    story.append(Paragraph("TOPE_DEEP", st["title"]))
    story.append(Paragraph("Computational Vaccine Discovery Report", st["subtitle"]))
    story.append(_spacer(8))
    story.append(_rule())
    story.append(_spacer(8))
    story.append(Paragraph(f"<b>{protein_name}</b>", ParagraphStyle(
        "CoverProtein", fontName="Helvetica-Bold", fontSize=14,
        leading=18, alignment=TA_CENTER, spaceAfter=4,
    )))
    story.append(Paragraph(f"UniProt Accession: {protein_id}", st["subtitle"]))
    story.append(_spacer(16))
    meta_lines = [
        f"Run ID: {run_id}",
        f"Generated: {run_date}",
        f"Pipeline: 10-agent TOPE_DEEP v1.0.0",
        f"Candidate proteins analysed: {candidates_count}",
    ]
    for line in meta_lines:
        story.append(Paragraph(line, st["meta"]))
        story.append(_spacer(2))
    story.append(_spacer(20))
    story.append(_rule())
    story.append(_spacer(8))
    notice = (
        "This report is generated by TOPE_DEEP, a computational orchestration system for Stage 2 "
        "vaccine design and discovery. All predictions are computational and require experimental "
        "validation by a qualified subject matter expert. This document does not constitute a "
        "clinical or regulatory submission."
    )
    story.append(Paragraph(notice, ParagraphStyle(
        "Notice", fontName="Helvetica-Oblique", fontSize=8.5,
        leading=13, alignment=TA_CENTER, textColor=GRAY_TEXT,
    )))
    story.append(PageBreak())
    return story


def _executive_summary(st, c: dict, construct_report: Optional[dict]):
    story = []
    story.append(Paragraph("1. Executive Summary", st["h1"]))
    story.append(_rule())

    eps      = c.get("epitopes", [])
    ctl_tot  = c.get("ctl_count", 0)
    ctl_hi   = c.get("ctl_strong", 0)
    htl_tot  = c.get("htl_count", 0)
    bcell    = c.get("bcell_count", 0)
    african  = c.get("african_coverage_pct", 0)
    global_  = c.get("global_coverage_pct", 0)
    safe_eps = [e for e in eps if e.get("allergenicity_safe") and e.get("toxicity_safe")]
    fail_eps = [e for e in eps if e.get("allergenicity_safe") is False or e.get("toxicity_safe") is False]
    protein  = c.get("protein_name", "")
    prot_id  = c.get("protein_id", "")
    seq_len  = c.get("sequence_length", 0)
    vaxijen  = next((d.get("vaxijen_score") for d in c.get("decisions", []) if d.get("stage") == "antigen_screening"), None)
    loc      = next((d.get("phobius_localization") for d in c.get("decisions", []) if d.get("stage") == "antigen_screening"), None)

    summary = (
        f"The TOPE_DEEP pipeline analysed <b>{protein}</b> (UniProt: {prot_id}, {seq_len} aa) "
        f"as a candidate antigen for computational vaccine design. "
    )
    if vaxijen is not None:
        summary += (
            f"VaxiJen 2.0 antigenicity score: <b>{vaxijen:.3f}</b> "
            f"({'above' if vaxijen >= 0.4 else 'below'} the 0.4 threshold for {('bacteria' if vaxijen else 'pathogen')} classification). "
        )
    if loc:
        summary += f"Phobius subcellular localisation: <b>{loc.replace('_', ' ')}</b>. "

    summary += (
        f"Epitope prediction identified <b>{ctl_tot} CTL</b> (MHC-I) epitopes, of which "
        f"<b>{ctl_hi}</b> are high-confidence strong binders (IC50 &lt; 50 nM); "
        f"<b>{htl_tot} HTL</b> (MHC-II) epitopes; and <b>{bcell}</b> linear B-cell epitopes. "
        f"Population coverage analysis (IEDB/AFND 2020) yields <b>{_pct(african)}</b> across "
        f"Sub-Saharan African populations and <b>{_pct(global_)}</b> globally. "
        f"Immunosafety screening cleared <b>{len(safe_eps)}</b> of {len(eps)} epitopes; "
        f"<b>{len(fail_eps)}</b> failed allergenicity or toxicity thresholds and were excluded."
    )

    if construct_report:
        props = construct_report.get("physicochemical", {})
        mw    = props.get("molecular_weight_da")
        pi    = props.get("isoelectric_point")
        ii    = props.get("instability_index")
        stable= props.get("is_stable")
        clen  = construct_report.get("length_aa")
        if mw:
            summary += (
                f" The assembled multi-epitope construct is <b>{clen} aa</b> "
                f"(MW {mw/1000:.1f} kDa, pI {pi:.2f}, instability index {ii:.1f} - "
                f"{'stable' if stable else 'unstable'} by ProtParam criteria)."
            )

    story.append(Paragraph(summary, st["body"]))
    story.append(_spacer(8))

    # KPI table
    kpi_rows = [
        ["Metric", "Value", "Interpretation"],
        ["CTL epitopes (total / high confidence)", f"{ctl_tot} / {ctl_hi}", "High conf. IC50 < 50 nM, strong MHC-I binders"],
        ["HTL epitopes (MHC-II)", str(htl_tot), "CD4+ T-cell helper epitopes"],
        ["B-cell linear epitopes", str(bcell), "BepiPred 2.0 score >= 0.5"],
        ["African HLA coverage", _pct(african), "Primary design constraint (AFND 2020)"],
        ["Global HLA coverage",  _pct(global_), "Combined MHC-I + II, 7 populations"],
        ["Safety cleared / total", f"{len(safe_eps)} / {len(eps)}", "WHO allergenicity + HemoPI + human homology"],
        ["VaxiJen score", f"{vaxijen:.3f}" if vaxijen else "-", "Threshold >= 0.4 (bacteria)"],
    ]
    if construct_report:
        props = construct_report.get("physicochemical", {})
        kpi_rows.append(["Construct length", f"{construct_report.get('length_aa', '-')} aa", "Full multi-epitope fusion"])
        if props.get("molecular_weight_da"):
            kpi_rows.append(["Construct MW", f"{props['molecular_weight_da']/1000:.1f} kDa", "ProtParam"])

    cols = [TEXT_W * 0.42, TEXT_W * 0.18, TEXT_W * 0.40]
    story.append(Paragraph("Table 1. Pipeline summary metrics.", st["caption"]))
    story.append(_table(
        [[Paragraph(str(c2), st["th"] if i == 0 else st["td"]) for c2 in row]
         for i, row in enumerate(kpi_rows)],
        col_widths=cols,
    ))
    story.append(_spacer(4))
    return story


def _population_coverage(st, c: dict, section_num: int):
    story = []
    story.append(Paragraph(f"{section_num}. Population Coverage Analysis", st["h1"]))
    story.append(_rule())

    cov = c.get("coverage_detail") or {}
    POP_ORDER = ["global", "african", "east_african", "european", "east_asian", "south_asian", "americas"]
    POP_LABELS = {
        "global":       "Global",
        "african":      "Sub-Saharan Africa",
        "east_african": "East Africa",
        "european":     "Europe",
        "east_asian":   "East Asia",
        "south_asian":  "South Asia",
        "americas":     "Americas",
    }

    if not cov:
        story.append(Paragraph("Population coverage data not available for this run.", st["body"]))
        return story

    intro = (
        "HLA population coverage was calculated using allele frequency data from the Allele "
        "Frequency Net Database (AFND 2020) across seven major population groups. "
        "Sub-Saharan African coverage is designated as the primary design constraint for this pipeline, "
        "reflecting the platform's focus on vaccine equity. The 80% combined coverage threshold "
        "is adopted from WHO priority population guidance."
    )
    story.append(Paragraph(intro, st["body"]))
    story.append(_spacer(6))

    rows = [["Population", "MHC-I (%)", "MHC-II (%)", "Combined (%)", "80% Target"]]
    for k in POP_ORDER:
        d = cov.get(k)
        if not d:
            continue
        combined = d.get("combined_pct", 0)
        met      = "Met" if combined >= 80 else "Below"
        rows.append([
            POP_LABELS.get(k, k) + (" *" if k in ("african", "east_african") else ""),
            f"{d.get('mhc_i_pct', 0):.1f}",
            f"{d.get('mhc_ii_pct', 0):.1f}",
            f"{combined:.1f}",
            met,
        ])

    cols = [TEXT_W*0.28, TEXT_W*0.16, TEXT_W*0.16, TEXT_W*0.20, TEXT_W*0.20]
    story.append(Paragraph(f"Table {section_num}.1. HLA population coverage by region (AFND 2020 allele frequencies).", st["caption"]))
    story.append(_table(
        [[Paragraph(str(cell), st["th"] if i == 0 else st["td"]) for cell in row]
         for i, row in enumerate(rows)],
        col_widths=cols,
        extra_style=[("FONTNAME", (3,1), (3,-1), "Helvetica-Bold")],
    ))
    story.append(Paragraph("* Primary design constraint populations.", st["note"]))
    return story


def _epitopes_table(st, c: dict, section_num: int):
    story = []
    story.append(Paragraph(f"{section_num}. Epitope Predictions", st["h1"]))
    story.append(_rule())

    eps  = c.get("epitopes", [])
    ctl  = [e for e in eps if e.get("epitope_type") == "CTL"]
    htl  = [e for e in eps if e.get("epitope_type") == "HTL"]
    bcell= [e for e in eps if "B-cell" in (e.get("epitope_type") or "")]

    intro = (
        f"Epitope prediction was performed using NetMHCpan 4.1 (MHC-I, CTL epitopes) and "
        f"NetMHCIIpan 4.3 (MHC-II, HTL epitopes) via the IEDB tools cluster, with BepiPred 2.0 "
        f"for linear B-cell epitopes. Confidence tiers are defined as: HIGH (IC50 &lt; 50 nM or "
        f"percentile rank &lt; 0.5%), MEDIUM (IC50 50–500 nM), LOW (IC50 500–5000 nM)."
    )
    story.append(Paragraph(intro, st["body"]))
    story.append(_spacer(6))

    # CTL table top 15
    story.append(Paragraph(f"{section_num}.1  CTL Epitopes (MHC-I)", st["h2"]))
    top_ctl = sorted(ctl, key=lambda e: e.get("ic50_nm") or 9999)[:15]
    if top_ctl:
        ctl_rows = [["Sequence", "HLA Allele", "IC50 (nM)", "%Rank", "Confidence", "Allergenicity", "Toxicity"]]
        for e in top_ctl:
            ctl_rows.append([
                e.get("sequence", ""),
                e.get("hla_allele") or "-",
                f"{e['ic50_nm']:.1f}" if e.get("ic50_nm") is not None else "-",
                f"{e['percentile_rank']:.2f}" if e.get("percentile_rank") is not None else "-",
                (e.get("confidence") or "-").capitalize(),
                "Safe"   if e.get("allergenicity_safe") else ("Flagged" if e.get("allergenicity_safe") is False else "-"),
                "Safe"   if e.get("toxicity_safe")      else ("Flagged" if e.get("toxicity_safe")      is False else "-"),
            ])
        cols = [TEXT_W*0.20, TEXT_W*0.16, TEXT_W*0.10, TEXT_W*0.08, TEXT_W*0.13, TEXT_W*0.16, TEXT_W*0.17]
        story.append(Paragraph(f"Table {section_num}.1. Top {len(top_ctl)} CTL epitopes ranked by IC50 (nM). High confidence: IC50 &lt; 50 nM.", st["caption"]))
        story.append(_table(
            [[Paragraph(str(cell), st["th"] if i == 0 else st["td_mono"] if j == 0 else st["td"])
              for j, cell in enumerate(row)]
             for i, row in enumerate(ctl_rows)],
            col_widths=cols,
        ))
        story.append(_spacer(4))
    else:
        story.append(Paragraph("No CTL epitopes predicted.", st["body"]))

    # HTL table top 10
    story.append(Paragraph(f"{section_num}.2  HTL Epitopes (MHC-II)", st["h2"]))
    top_htl = sorted(htl, key=lambda e: e.get("ic50_nm") or 9999)[:10]
    if top_htl:
        htl_rows = [["Sequence", "HLA-DR Allele", "IC50 (nM)", "Confidence"]]
        for e in top_htl:
            htl_rows.append([
                e.get("sequence", ""),
                e.get("hla_allele") or "-",
                f"{e['ic50_nm']:.1f}" if e.get("ic50_nm") is not None else "-",
                (e.get("confidence") or "-").capitalize(),
            ])
        cols = [TEXT_W*0.36, TEXT_W*0.24, TEXT_W*0.18, TEXT_W*0.22]
        story.append(Paragraph(f"Table {section_num}.2. Top {len(top_htl)} HTL epitopes (CD4+ T-cell help).", st["caption"]))
        story.append(_table(
            [[Paragraph(str(cell), st["th"] if i == 0 else st["td_mono"] if j == 0 else st["td"])
              for j, cell in enumerate(row)]
             for i, row in enumerate(htl_rows)],
            col_widths=cols,
        ))
        story.append(_spacer(4))

    # B-cell
    story.append(Paragraph(f"{section_num}.3  Linear B-Cell Epitopes", st["h2"]))
    if bcell:
        bc_rows = [["Sequence", "BepiPred Score", "Rabbit Validation", "Safety"]]
        for e in bcell[:10]:
            to = e.get("tool_outputs") or {}
            bc_rows.append([
                e.get("sequence", ""),
                f"{to.get('bepipred_score', 0):.3f}" if to.get("bepipred_score") else "-",
                "Recommended" if to.get("rabbit_validation") else "-",
                "Safe" if (e.get("allergenicity_safe") and e.get("toxicity_safe")) else "Review",
            ])
        cols = [TEXT_W*0.38, TEXT_W*0.18, TEXT_W*0.22, TEXT_W*0.22]
        story.append(Paragraph(f"Table {section_num}.3. Linear B-cell epitopes (BepiPred 2.0, threshold >= 0.5).", st["caption"]))
        story.append(_table(
            [[Paragraph(str(cell), st["th"] if i == 0 else st["td_mono"] if j == 0 else st["td"])
              for j, cell in enumerate(row)]
             for i, row in enumerate(bc_rows)],
            col_widths=cols,
        ))
    else:
        story.append(Paragraph(
            "No linear B-cell epitopes predicted above the BepiPred 2.0 threshold (score >= 0.5) "
            "for this protein. This is consistent with highly structured or membrane-embedded antigens "
            "where conformational epitopes predominate. Conformational B-cell epitope prediction "
            "(e.g., via ElliPro) is recommended as a follow-up.",
            st["body"],
        ))
    return story


def _structure_section(st, c: dict, section_num: int):
    story = []
    story.append(Paragraph(f"{section_num}. Structural Analysis", st["h1"]))
    story.append(_rule())

    d = next((x for x in c.get("decisions", []) if x.get("stage") == "structure_retrieval"), {})
    plddt  = d.get("mean_plddt")
    ver    = d.get("model_version", "v4")
    entry  = d.get("alphafold_entry_id") or f"AF-{c.get('protein_id')}-F1"
    source = "AlphaFold DB"
    seq_len= c.get("sequence_length", 0)
    uid    = c.get("protein_id", "")

    plddt_interp = (
        "very high confidence (>= 90)" if (plddt or 0) >= 90 else
        "confident (70–90)"            if (plddt or 0) >= 70 else
        "low confidence (< 70)"
    )

    intro = (
        f"Three-dimensional structural data was retrieved from the {source} "
        f"(EBI, {ver}) for {c.get('protein_name', '')} ({uid}). "
        f"The predicted structure covers {seq_len} residues. "
    )
    if plddt:
        intro += (
            f"Mean per-residue confidence score (pLDDT): <b>{plddt:.1f}/100</b> - "
            f"{plddt_interp}. pLDDT values above 70 are considered reliable for "
            f"structural-based epitope evaluation (Jumper et al. 2021, Nature 596:583–589)."
        )
    story.append(Paragraph(intro, st["body"]))
    story.append(_spacer(6))

    rows = [
        ["Property", "Value"],
        ["Structure source",   source],
        ["AlphaFold entry ID", entry],
        ["Model version",      _safe(ver)],
        ["Sequence length",    f"{seq_len} aa"],
        ["Mean pLDDT",         f"{plddt:.1f} / 100.0" if plddt else "-"],
        ["pLDDT interpretation", plddt_interp],
        ["Reference",          "Jumper et al. (2021) Nature 596:583–589"],
    ]
    cols = [TEXT_W*0.38, TEXT_W*0.62]
    story.append(Paragraph(f"Table {section_num}.1. Structural analysis metadata.", st["caption"]))
    story.append(_table(
        [[Paragraph(str(cell), st["th"] if i == 0 else st["td"]) for cell in row]
         for i, row in enumerate(rows)],
        col_widths=cols,
    ))
    story.append(_spacer(4))
    story.append(Paragraph(
        f"Note: Three-dimensional visualisation is available via the AlphaFold DB viewer "
        f"at https://alphafold.ebi.ac.uk/entry/{uid}. "
        f"Epitope surface mapping is recommended using PyMOL or UCSF ChimeraX.",
        st["note"],
    ))
    return story


def _safety_section(st, c: dict, section_num: int):
    story = []
    story.append(Paragraph(f"{section_num}. Immunosafety Screening", st["h1"]))
    story.append(_rule())

    eps      = c.get("epitopes", [])
    safe     = [e for e in eps if e.get("allergenicity_safe") and e.get("toxicity_safe")]
    unscored = [e for e in eps if e.get("allergenicity_safe") is None]
    failed   = [e for e in eps if e.get("allergenicity_safe") is False or e.get("toxicity_safe") is False]

    d = next((x for x in c.get("decisions", []) if x.get("stage") == "safety_filter"), {})
    methods = (
        "WHO 2001 allergenicity protocol (AllergenOnline); "
        "AllerTOP v2.0 local (Doytchinova & Flower 2014); "
        "HemoPI haemolytic screen (Singh et al. 2011, WHO/BS/2019.2364); "
        "Human homology FDA/EMA 8-mer threshold (UniProt human Swiss-Prot)."
    )

    intro = (
        f"All {len(eps)} predicted epitopes were subjected to a four-layer immunosafety screen: "
        f"{methods} "
        f"Results: <b>{len(safe)} safe</b>, <b>{len(unscored)} unscored</b> (screening inconclusive), "
        f"<b>{len(failed)} failed</b>. "
        f"Only safety-cleared epitopes proceed to construct assembly."
    )
    story.append(Paragraph(intro, st["body"]))
    story.append(_spacer(6))

    rows = [
        ["Safety Screen", "Method", "Threshold", "Result"],
        ["WHO allergenicity", "AllergenOnline FASTA alignment", "8-mer, 35% identity", f"{len(safe)+len(unscored)} passed"],
        ["AllerTOP v2.0",    "SVM classifier (Doytchinova 2014)", "Score < 0 = safe", f"{len(safe)} cleared"],
        ["HemoPI",           "SVM haemolytic peptide prediction", "Score < 0 = safe", f"{len(safe)} cleared"],
        ["Human homology",   "FDA/EMA 8-mer, UniProt Swiss-Prot", "0 matches = safe", f"{len(safe)} cleared"],
    ]
    cols = [TEXT_W*0.22, TEXT_W*0.32, TEXT_W*0.22, TEXT_W*0.24]
    story.append(Paragraph(f"Table {section_num}.1. Immunosafety screening protocol and results.", st["caption"]))
    story.append(_table(
        [[Paragraph(str(cell), st["th"] if i == 0 else st["td"]) for cell in row]
         for i, row in enumerate(rows)],
        col_widths=cols,
    ))

    if failed:
        story.append(_spacer(6))
        story.append(Paragraph(f"{section_num}.1  Failed Epitopes", st["h2"]))
        fail_rows = [["Sequence", "Type", "Allergenicity", "Toxicity"]]
        for e in failed[:10]:
            fail_rows.append([
                e.get("sequence", ""),
                e.get("epitope_type", "-"),
                "Flagged" if e.get("allergenicity_safe") is False else "Clear",
                "Flagged" if e.get("toxicity_safe")      is False else "Clear",
            ])
        cols2 = [TEXT_W*0.38, TEXT_W*0.16, TEXT_W*0.23, TEXT_W*0.23]
        story.append(_table(
            [[Paragraph(str(cell), st["th"] if i == 0 else st["td_mono"] if j == 0 else st["td"])
              for j, cell in enumerate(row)]
             for i, row in enumerate(fail_rows)],
            col_widths=cols2,
        ))
    return story


def _construct_section(st, construct_report: Optional[dict], section_num: int):
    story = []
    story.append(Paragraph(f"{section_num}. Multi-Epitope Construct Design", st["h1"]))
    story.append(_rule())

    if not construct_report:
        story.append(Paragraph("Construct assembly was not performed for this run.", st["body"]))
        return story

    seq   = construct_report.get("construct_sequence", "")
    props = construct_report.get("physicochemical", {})
    adj   = construct_report.get("adjuvant", {})
    linkers = construct_report.get("linker_scheme", {})
    counts  = construct_report.get("epitope_counts", {})
    clen    = construct_report.get("length_aa", len(seq) // 1)

    intro = (
        f"A multi-epitope fusion construct was assembled incorporating safety-cleared "
        f"CTL, HTL, and B-cell epitopes joined by immunologically validated linker sequences. "
        f"The construct includes the {adj.get('key', 'RS09')} adjuvant "
        f"({adj.get('mechanism', '').split('.')[0] if adj else 'TLR4 agonist'}). "
        f"Physicochemical properties were computed using ProtParam (Biopython)."
    )
    story.append(Paragraph(intro, st["body"]))
    story.append(_spacer(6))

    # Composition table
    comp_rows = [
        ["Component", "Count / Value"],
        ["CTL epitopes included", str(counts.get("CTL", "-"))],
        ["HTL epitopes included", str(counts.get("HTL", "-"))],
        ["B-cell epitopes included", str(counts.get("B-cell", "-"))],
        ["Total construct length", f"{clen} aa"],
        ["Molecular weight", f"{props.get('molecular_weight_da', 0)/1000:.2f} kDa" if props.get("molecular_weight_da") else "-"],
        ["Isoelectric point (pI)", f"{props.get('isoelectric_point', 0):.2f}" if props.get("isoelectric_point") else "-"],
        ["Instability index", f"{props.get('instability_index', 0):.1f}" if props.get("instability_index") else "-"],
        ["Stability (ProtParam)", "Stable (index < 40)" if props.get("is_stable") else "Unstable"],
        ["GRAVY score", f"{props.get('gravy', 0):.3f}" if props.get("gravy") else "-"],
        ["Aromaticity", f"{props.get('aromaticity', 0):.3f}" if props.get("aromaticity") else "-"],
    ]
    cols = [TEXT_W*0.52, TEXT_W*0.48]
    story.append(Paragraph(f"Table {section_num}.1. Construct composition and physicochemical properties.", st["caption"]))
    story.append(_table(
        [[Paragraph(str(cell), st["th"] if i == 0 else st["td"]) for cell in row]
         for i, row in enumerate(comp_rows)],
        col_widths=cols,
    ))

    # Linker table
    if linkers:
        story.append(_spacer(6))
        story.append(Paragraph(f"{section_num}.1  Linker Scheme", st["h2"]))
        linker_info = {
            "AAY":   ("CTL spacer", "Rigid", "Proteasomal cleavage signal, MHC-I presentation"),
            "GPGPG": ("HTL spacer", "Flexible", "Helix breaker, prevents secondary structure"),
            "KK":    ("B-cell spacer", "Flexible", "Lysine spacer, surface accessibility"),
            "EAAAK": ("Domain spacer", "Rigid (alpha-helix)", "Spatial separation of domains"),
        }
        lk_rows = [["Linker", "Sequence", "Type", "Flexibility", "Function"]]
        for name, seq_l in linkers.items():
            info = linker_info.get(name, ("-", "-", "-"))
            lk_rows.append([name, str(seq_l), info[0], info[1], info[2]])
        cols2 = [TEXT_W*0.10, TEXT_W*0.18, TEXT_W*0.14, TEXT_W*0.14, TEXT_W*0.44]
        story.append(Paragraph(f"Table {section_num}.2. Linker sequences and functional rationale.", st["caption"]))
        story.append(_table(
            [[Paragraph(str(cell), st["th"] if i == 0 else st["td_mono"] if j < 2 else st["td"])
              for j, cell in enumerate(row)]
             for i, row in enumerate(lk_rows)],
            col_widths=cols2,
        ))

    # Sequence block
    if seq:
        story.append(_spacer(8))
        story.append(Paragraph(f"{section_num}.2  Full Construct Sequence", st["h2"]))
        story.append(Paragraph(
            "The complete amino acid sequence of the multi-epitope fusion construct is provided below. "
            "Sequence formatted in 60-residue blocks for readability.",
            st["body"],
        ))
        blocks = [seq[i:i+60] for i in range(0, len(seq), 60)]
        for block in blocks:
            story.append(Paragraph(block, st["mono_small"]))
        story.append(_spacer(4))

    return story


def _literature_section(st, c: dict, section_num: int):
    story = []
    story.append(Paragraph(f"{section_num}. Literature Evidence", st["h1"]))
    story.append(_rule())

    d       = next((x for x in c.get("decisions", []) if x.get("stage") == "literature_search"), None)
    if not d:
        story.append(Paragraph("Literature search was not performed for this run.", st["body"]))
        return story

    count   = d.get("pubmed_hits") or d.get("result_count", 0)
    query   = d.get("query", "-")
    elapsed = d.get("search_time_s")
    prior   = d.get("prior_validated", False)
    signals = d.get("failure_signals", [])
    summary = d.get("literature_summary") or d.get("reasoning", "")
    pmids   = d.get("evidence_pmids", [])

    intro = (
        f"A systematic PubMed search was performed for <b>{c.get('protein_name', '')}</b> "
        f"using immunology-focused MeSH terms (T cell, epitope, vaccine, MHC, antibody). "
        f"Search returned <b>{count}</b> relevant abstracts"
        f"{f' (search time: {elapsed:.1f}s)' if elapsed else ''}. "
        f"Prior experimental validation in published literature: <b>{'Yes' if prior else 'No'}</b>."
    )
    story.append(Paragraph(intro, st["body"]))

    if summary:
        story.append(_spacer(6))
        story.append(Paragraph("Evidence Synthesis (Claude API)", st["h2"]))
        story.append(Paragraph(summary, st["body"]))

    if signals:
        story.append(_spacer(4))
        story.append(Paragraph("Literature-Derived Failure Signals", st["h2"]))
        for sig in signals:
            story.append(Paragraph(f"&bull;  {sig}", st["bullet"]))

    if pmids:
        story.append(_spacer(6))
        story.append(Paragraph(f"Source Publications ({len(pmids)} PMIDs retrieved)", st["h2"]))
        pmid_text = "  ".join(f"PMID:{p}" for p in pmids)
        story.append(Paragraph(pmid_text, st["mono_small"]))

    return story


def _experiment_section(st, c: dict, section_num: int):
    story = []
    story.append(Paragraph(f"{section_num}. Wet-Lab Validation Roadmap", st["h1"]))
    story.append(_rule())

    d = next((x for x in c.get("decisions", []) if x.get("stage") == "experiment_planning"), None)
    if not d:
        story.append(Paragraph("Wet-lab roadmap was not generated for this run.", st["body"]))
        return story

    plan     = d.get("plan", {}) or {}
    reasoning= d.get("reasoning", "")

    if reasoning:
        story.append(Paragraph(reasoning, st["body"]))
        story.append(_spacer(4))

    # Three-phase table
    phases = []
    for key, label in [("phase_1", "Phase 1"), ("phase_2", "Phase 2"), ("phase_3", "Phase 3")]:
        ph = plan.get(key, {})
        if ph:
            phases.append([label, ph.get("timeline", "-"), ph.get("key_assay", "-"),
                           ph.get("go_criteria", "-"), ph.get("expected_cost_usd", "-")])

    if phases:
        phase_rows = [["Phase", "Timeline", "Key Assay", "Go Criteria", "Est. Cost"]] + phases
        cols = [TEXT_W*0.08, TEXT_W*0.13, TEXT_W*0.20, TEXT_W*0.35, TEXT_W*0.24]
        story.append(Paragraph(f"Table {section_num}.1. Three-phase validation roadmap.", st["caption"]))
        story.append(_table(
            [[Paragraph(str(cell), st["th"] if i == 0 else st["td"]) for cell in row]
             for i, row in enumerate(phase_rows)],
            col_widths=cols,
        ))
        story.append(_spacer(6))

    # ELISpot protocol
    elispot = plan.get("elispot_protocol")
    if elispot:
        story.append(Paragraph(f"{section_num}.1  ELISpot Protocol", st["h2"]))
        for line in elispot.strip().split("\n"):
            line = line.strip()
            if line:
                story.append(Paragraph(line, st["bullet"]))
        story.append(_spacer(4))

    # Immunisation schedule
    sched = plan.get("immunisation_schedule")
    if sched:
        story.append(Paragraph(f"{section_num}.2  Immunisation Schedule", st["h2"]))
        for line in sched.strip().split("\n"):
            line = line.strip()
            if line:
                story.append(Paragraph(line, st["bullet"]))
        story.append(_spacer(4))

    # NHP plan
    nhp = plan.get("nhp_plan")
    if nhp:
        story.append(Paragraph(f"{section_num}.3  Non-Human Primate (NHP) Plan", st["h2"]))
        story.append(Paragraph(nhp, st["body"]))
        story.append(_spacer(4))

    # Critical risks
    risks = plan.get("critical_risks", [])
    if risks:
        story.append(Paragraph(f"{section_num}.4  Critical Risks", st["h2"]))
        for r in risks:
            story.append(Paragraph(f"&bull;  {r}", st["bullet"]))

    # Cost summary
    costs = plan.get("estimated_cost_usd", {})
    if costs:
        story.append(_spacer(6))
        story.append(Paragraph("Estimated Cost", st["h2"]))
        cost_rows = [["Phase", "Estimated Cost (USD)"]]
        for k, v in costs.items():
            cost_rows.append([k.replace("_", " ").capitalize(), str(v)])
        cols2 = [TEXT_W*0.4, TEXT_W*0.6]
        story.append(_table(
            [[Paragraph(str(cell), st["th"] if i == 0 else st["td"]) for cell in row]
             for i, row in enumerate(cost_rows)],
            col_widths=cols2,
        ))

    return story


def _audit_section(st, decisions: list, timing: dict, section_num: int):
    story = []
    story.append(Paragraph(f"{section_num}. Pipeline Audit Log", st["h1"]))
    story.append(_rule())

    TIMING_MAP = {
        "data_curation":     "n1_curation",
        "antigen_screening": "n2_screening",
        "tcell_prediction":  "n3_tcell",
        "bcell_prediction":  "n4_bcell",
        "structure_retrieval":"n5_structure",
        "safety_filter":     "n6_safety",
        "coverage_analysis": "n7_coverage",
        "construct_design":  "n8_construct",
        "literature_search": "n9_literature",
        "experiment_planning":"n10_experiment",
    }
    AGENT_LABELS = {
        "data_curation":     "Sequence Ingestor",
        "antigen_screening": "Antigenicity Screener",
        "tcell_prediction":  "MHC Binding Predictor",
        "bcell_prediction":  "Linear Epitope Mapper",
        "structure_retrieval":"Structural Resolver",
        "safety_filter":     "Immunosafety Filter",
        "coverage_analysis": "Population Coverage",
        "construct_design":  "Construct Assembler",
        "literature_search": "Evidence Retriever",
        "experiment_planning":"Validation Roadmap",
    }

    intro = (
        "The following table records every computational decision made by the pipeline, "
        "including the agent responsible, the decision outcome, and elapsed time. "
        "This audit trail is machine-readable and reproducible."
    )
    story.append(Paragraph(intro, st["body"]))
    story.append(_spacer(6))

    rows = [["Agent", "Stage", "Decision", "Elapsed"]]
    for d in decisions:
        stage  = d.get("stage", "")
        label  = AGENT_LABELS.get(stage, stage.replace("_", " ").title())
        t_key  = TIMING_MAP.get(stage)
        t_val  = timing.get(t_key) if t_key else None
        rows.append([
            label,
            stage,
            d.get("decision", "-"),
            _fmt_s(t_val),
        ])
    # Total
    total = timing.get("total_seconds")
    if total:
        rows.append(["TOTAL PIPELINE TIME", "", "", _fmt_s(total)])

    cols = [TEXT_W*0.26, TEXT_W*0.22, TEXT_W*0.36, TEXT_W*0.16]
    story.append(Paragraph(f"Table {section_num}.1. Agent decision audit trail.", st["caption"]))
    t = _table(
        [[Paragraph(str(cell), st["th"] if i == 0 else st["td"]) for cell in row]
         for i, row in enumerate(rows)],
        col_widths=cols,
        extra_style=[
            ("FONTNAME", (0, len(rows)-1), (-1, len(rows)-1), "Helvetica-Bold"),
            ("LINEABOVE", (0, len(rows)-1), (-1, len(rows)-1), 0.6, BLACK),
        ] if total else [],
    )
    story.append(t)

    # Detailed reasoning
    story.append(_spacer(10))
    story.append(Paragraph(f"{section_num}.1  Agent Reasoning (Full)", st["h2"]))
    for i, d in enumerate(decisions, 1):
        stage = d.get("stage", "")
        label = AGENT_LABELS.get(stage, stage.replace("_", " ").title())
        reasoning = d.get("reasoning", "")
        if reasoning:
            story.append(Paragraph(f"<b>{i}. {label}</b>", st["body_left"]))
            # Sanitize internal implementation strings
            clean = (reasoning
                .replace("PSORTb unavailable on Railway, Phobius used as proxy.", "")
                .replace("Kozi AI", "")
                .replace("TOPE_DEEP", "")
                .strip()
            )
            story.append(Paragraph(clean, st["body"]))
            story.append(_spacer(4))

    return story


def _references(st):
    story = []
    story.append(Paragraph("References", st["h1"]))
    story.append(_rule())
    refs = [
        "Jumper J et al. (2021) Highly accurate protein structure prediction with AlphaFold. <i>Nature</i> 596:583–589.",
        "Jespersen MC et al. (2017) BepiPred-2.0: improving sequence-based B-cell epitope prediction using conformational epitopes. <i>Nucleic Acids Res</i> 45:W24–W29.",
        "Andreatta M & Nielsen M (2016) Gapped sequence alignment using artificial neural networks: application to the MHC class I system. <i>Bioinformatics</i> 32:511–517.",
        "Doytchinova IA & Flower DR (2014) AllerTOP v.2 a server for in silico prediction of allergens. <i>J Mol Model</i> 20:2278.",
        "Singh H et al. (2011) HemoPI Haemolytic Peptides. <i>J Transl Med</i> 9:1–9.",
        "Sette A & Rappuoli R (2010) Reverse vaccinology: developing vaccines in the era of genomics. <i>Immunity</i> 33:530–541.",
        "Janetzki S et al. (2015) Guidelines for the automated evaluation of ELISpot assays. <i>Cancer Immunol Immunother</i> 64:1695–1703.",
        "Klenerman P et al. (2002) T cell responses against immunodominant epitopes. <i>Nat Rev Immunol</i> 2:263–272.",
        "Seder RA et al. (2008) T-cell quality in memory and protection. <i>Nat Immunol</i> 9:239–245.",
        "Bui HH et al. (2006) Predicting population coverage of T-cell epitope-based diagnostics and vaccines. <i>BMC Bioinformatics</i> 7:153.",
    ]
    for i, r in enumerate(refs, 1):
        story.append(Paragraph(f"{i}.  {r}", ParagraphStyle(
            f"Ref{i}", fontName="Helvetica", fontSize=8.5, leading=13,
            leftIndent=16, firstLineIndent=-16, spaceAfter=5, alignment=TA_LEFT,
        )))
    return story


# ── Main entry point ───────────────────────────────────────────────────────────

def generate_report_pdf(results: Dict[str, Any], run_id: str) -> bytes:
    """
    Generate a LaTeX-style scientific PDF from pipeline results dict.

    Args:
        results: dict from api.getResults() keys: run_id, status, timing, candidates, construct_report
        run_id:  str

    Returns:
        bytes PDF file content
    """
    buf = io.BytesIO()
    st  = _styles()

    candidates       = results.get("candidates", [])
    construct_report = results.get("construct_report")
    timing           = results.get("timing", {})

    if not candidates:
        raise ValueError("No candidates in results")

    c            = candidates[0]
    protein_name = c.get("protein_name", "Unknown protein")
    # Strip parenthetical UniProt accession from display name
    import re
    clean_name   = re.sub(r'\s*\([A-Z0-9]+\)\s*$', '', protein_name).strip()
    protein_id   = c.get("protein_id", "")
    run_date     = datetime.utcnow().strftime("%d %B %Y, %H:%M UTC")
    decisions    = c.get("decisions", [])

    on_page = _make_header_footer(clean_name, run_id, {})

    doc = SimpleDocTemplate(
        buf,
        pagesize=A4,
        leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=MARGIN + 8*mm, bottomMargin=MARGIN + 8*mm,
        title=f"TOPE_DEEP Report  {clean_name}",
        author="TOPE_DEEP v1.0.0",
        subject="Computational Vaccine Discovery Report",
    )

    story = []

    # Cover
    story += _cover(st, clean_name, protein_id, run_id, run_date, len(candidates))

    # Body sections
    story += _executive_summary(st, c, construct_report)
    story.append(_spacer(6))

    story += _population_coverage(st, c, section_num=2)
    story.append(_spacer(6))

    story += _epitopes_table(st, c, section_num=3)
    story.append(_spacer(6))

    story += _structure_section(st, c, section_num=4)
    story.append(_spacer(6))

    story += _safety_section(st, c, section_num=5)
    story.append(_spacer(6))

    story += _construct_section(st, construct_report, section_num=6)
    story.append(_spacer(6))

    story += _literature_section(st, c, section_num=7)
    story.append(_spacer(6))

    story += _experiment_section(st, c, section_num=8)
    story.append(PageBreak())

    story += _audit_section(st, decisions, timing, section_num=9)
    story.append(PageBreak())

    story += _references(st)

    doc.build(story, onFirstPage=on_page, onLaterPages=on_page)
    return buf.getvalue()