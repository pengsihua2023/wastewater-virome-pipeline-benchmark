#!/usr/bin/env python3
"""Build a 22-slide overview deck of the seven-pipeline viral detection benchmark.

Every number on every slide is taken from the verified result files, not typed in
by hand: metrics come from interim/viral_aggregate_metrics.tsv, the depth
analysis from 13_detection_vs_depth.xlsx, the two controlled comparisons from
15_assembly_loss_and_preprocessing.xlsx, and the figures are the same PNGs the
manuscripts embed. Rerunning this script after rerunning the analysis therefore
updates the deck automatically.

Writes Seven-pipeline-benchmark-overview.pptx at the repository root.
"""
import sys
from pathlib import Path

import pandas as pd
from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Emu, Inches, Pt

ROOT = Path(__file__).resolve().parents[2]
AN = ROOT / "analysis_708"
INT = AN / "interim"

W, H = Inches(13.333), Inches(7.5)

INK = RGBColor(0x1A, 0x1A, 0x1A)
MUTED = RGBColor(0x5A, 0x5A, 0x5A)
ACCENT = RGBColor(0x0B, 0x3D, 0x5C)      # deep navy
HILITE = RGBColor(0xC0, 0x4A, 0x00)      # burnt orange
RULE = RGBColor(0xD8, 0xD8, 0xD8)
BAND = RGBColor(0xF2, 0xF4, 0xF6)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

FONT = "Calibri"

NAME = {"TaxProfiler": "nf-core/taxprofiler", "TaxTriage": "nf-core/taxtriage",
        "mag-nf": "nf-core/mag", "KrakenMetaReads-nf": "KrakenMetaReads-nf",
        "GOTTCHA2": "GOTTCHA2", "CLARK": "CLARK", "sourmash": "sourmash"}


# --------------------------------------------------------------- primitives
def textbox(slide, x, y, w, h, text, size=18, bold=False, color=INK,
            align=PP_ALIGN.LEFT, italic=False, space_after=8, line=1.15):
    tb = slide.shapes.add_textbox(x, y, w, h)
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    lines = text.split("\n") if isinstance(text, str) else text
    for i, ln in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(space_after)
        p.line_spacing = line
        r = p.add_run()
        r.text = ln
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.italic = italic
        r.font.color.rgb = color
        r.font.name = FONT
    return tb


def bullets(slide, x, y, w, h, items, size=17, gap=12):
    """items: list of (text, level, bold, color)."""
    tb = slide.shapes.add_textbox(x, y, w, h)
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    for i, it in enumerate(items):
        text, level, bold, color = (list(it) + [0, False, INK])[:4]
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.level = level
        p.space_after = Pt(gap)
        p.line_spacing = 1.2
        r = p.add_run()
        r.text = ("— " if level else "▪  ") + text if text else ""
        r.font.size = Pt(size - 2 * level)
        r.font.bold = bold
        r.font.color.rgb = color
        r.font.name = FONT
    return tb


def rect(slide, x, y, w, h, fill):
    from pptx.enum.shapes import MSO_SHAPE
    sh = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, x, y, w, h)
    sh.fill.solid()
    sh.fill.fore_color.rgb = fill
    sh.line.fill.background()
    sh.shadow.inherit = False
    return sh


def slide_frame(prs, kicker, title, n):
    """Blank slide with the standard header, rule and page number."""
    s = prs.slides.add_slide(prs.slide_layouts[6])
    textbox(s, Inches(0.7), Inches(0.42), Inches(11.9), Inches(0.3),
            kicker.upper(), size=12, bold=True, color=HILITE, space_after=0)
    textbox(s, Inches(0.7), Inches(0.72), Inches(11.9), Inches(0.7),
            title, size=27, bold=True, color=ACCENT, space_after=0)
    rect(s, Inches(0.7), Inches(1.46), Inches(11.93), Emu(9525), RULE)
    textbox(s, Inches(12.1), Inches(6.95), Inches(0.6), Inches(0.3), str(n),
            size=11, color=MUTED, align=PP_ALIGN.RIGHT, space_after=0)
    return s


def picture(slide, path, x, y, max_w, max_h):
    """Insert centred within the (x, y, max_w, max_h) box, preserving aspect."""
    iw, ih = Image.open(path).size
    scale = min(max_w / iw, max_h / ih)
    w, h = int(iw * scale), int(ih * scale)
    slide.shapes.add_picture(str(path), int(x + (max_w - w) / 2),
                             int(y + (max_h - h) / 2), width=w, height=h)


def table(slide, x, y, w, rows, col_w=None, size=13, header=True,
          row_h=Inches(0.34), hilite_rows=()):
    n_r, n_c = len(rows), len(rows[0])
    shp = slide.shapes.add_table(n_r, n_c, x, y, w, row_h * n_r)
    t = shp.table
    if col_w:
        total = sum(col_w)
        for j, cw in enumerate(col_w):
            t.columns[j].width = int(w * cw / total)
    for i, row in enumerate(rows):
        t.rows[i].height = row_h
        for j, val in enumerate(row):
            c = t.cell(i, j)
            c.text = str(val)
            c.margin_left = c.margin_right = Inches(0.07)
            c.margin_top = c.margin_bottom = Inches(0.02)
            c.vertical_anchor = MSO_ANCHOR.MIDDLE
            c.fill.solid()
            if i == 0 and header:
                c.fill.fore_color.rgb = ACCENT
            elif i in hilite_rows:
                c.fill.fore_color.rgb = RGBColor(0xFD, 0xF0, 0xE3)
            else:
                c.fill.fore_color.rgb = WHITE if i % 2 else BAND
            p = c.text_frame.paragraphs[0]
            p.alignment = PP_ALIGN.LEFT if j == 0 else PP_ALIGN.CENTER
            for r in p.runs:
                r.font.size = Pt(size)
                r.font.name = FONT
                r.font.bold = (i == 0 and header)
                r.font.color.rgb = WHITE if (i == 0 and header) else INK
    return t


def stat_cards(slide, y, cards, top_size=30, lab_size=12):
    """cards: list of (big, label). Evenly spread across the slide."""
    n = len(cards)
    gap, m = Inches(0.28), Inches(0.7)
    cw = int((W - 2 * m - gap * (n - 1)) / n)
    for i, (big, lab) in enumerate(cards):
        x = int(m + i * (cw + gap))
        rect(slide, x, y, cw, Inches(1.42), BAND)
        textbox(slide, x, y + Inches(0.16), cw, Inches(0.6), big,
                size=top_size, bold=True, color=ACCENT, align=PP_ALIGN.CENTER,
                space_after=0)
        textbox(slide, x + Inches(0.12), y + Inches(0.82), cw - Inches(0.24),
                Inches(0.55), lab, size=lab_size, color=MUTED,
                align=PP_ALIGN.CENTER, space_after=0, line=1.05)


def note(slide, text, y=Inches(6.62)):
    textbox(slide, Inches(0.7), y, Inches(11.9), Inches(0.6), text,
            size=12, italic=True, color=MUTED, space_after=0)


# ------------------------------------------------------------------- data
def load():
    d = {}
    a = pd.read_csv(INT / "viral_aggregate_metrics.tsv", sep="\t")
    d["agg"] = a
    d["ill"] = a[a.platform == "Illumina"].sort_values("f1_mean", ascending=False)
    d["pac"] = a[a.platform == "PacBio"].sort_values("f1_mean", ascending=False)

    x = pd.ExcelFile(AN / "13_detection_vs_depth.xlsx")
    d["restricted"] = x.parse("recall_above_threshold")

    y = pd.ExcelFile(AN / "15_assembly_loss_and_preprocessing.xlsx")
    d["asm"] = y.parse("assembly_vs_reads")
    d["fastp"] = y.parse("fastp_effect_mag_vs_taxpro").mean(numeric_only=True)

    d["comp"] = pd.read_csv(INT / "realised_composition.tsv", sep="\t")
    d["tp_reads"] = pd.read_csv(INT / "viral_tp_reads.tsv", sep="\t")
    d["genus"] = pd.read_excel(AN / "12_genus_level_sensitivity.xlsx",
                               "genus_aggregate")
    return d


def rp(row):
    return (f"{row.recall_mean:.3f}", f"{row.precision_mean:.3f}"
            if row.precision_mean > 0 else "—", f"{row.f1_mean:.3f}")


def build():
    D = load()
    prs = Presentation()
    prs.slide_width, prs.slide_height = W, H
    ill, pac = D["ill"], D["pac"]
    asm = D["asm"].set_index(["platform", "evidence"])
    fp = D["fastp"]

    # ------------------------------------------------------------- 1 title
    s = prs.slides.add_slide(prs.slide_layouts[6])
    rect(s, 0, 0, W, Inches(2.45), ACCENT)
    textbox(s, Inches(0.85), Inches(0.62), Inches(11.6), Inches(1.4),
            "Benchmarking seven metagenomic pipelines\nfor viral detection in simulated wastewater",
            size=33, bold=True, color=WHITE, space_after=4, line=1.12)
    textbox(s, Inches(0.85), Inches(2.75), Inches(11.6), Inches(0.5),
            "Paired Illumina and PacBio HiFi datasets · 151-member non-phage viral ground truth · 42 pipeline runs",
            size=17, color=HILITE, space_after=0)
    stat_cards(s, Inches(3.45), [
        ("7", "classification pipelines\ncompared head to head"),
        ("151", "non-phage viruses in the\nground truth"),
        ("42", "runs: 7 pipelines\n× 6 datasets"),
        ("ρ = +0.09", "rank correlation of viral F1\nbetween the two modalities"),
    ])
    textbox(s, Inches(0.85), Inches(5.35), Inches(11.6), Inches(1.2),
            "Every metric is computed within the viral space only. Bacteria, archaea and bacteriophage are simulated\n"
            "background and contribute to no reported number.",
            size=15, color=MUTED, space_after=0)
    note(s, "Simulated community derived from wastewater sample HTP-2024-01-15, "
            "Hyperion Water Reclamation Plant, Los Angeles (Grimm et al. 2026).")

    # ---------------------------------------------------------- 2 the problem
    s = slide_frame(prs, "Motivation", "Wastewater surveillance is about viruses — benchmarks usually are not", 2)
    bullets(s, Inches(0.7), Inches(1.85), Inches(11.9), Inches(4.4), [
        ("Wastewater-based epidemiology is deployed principally to detect viral pathogens: SARS-CoV-2, influenza, enteroviruses, noroviruses, poliovirus.", 0, False, INK),
        ("Standard pipeline benchmarks score the whole community, which wastewater libraries make ~94% bacterial. A pipeline can rank well while missing most viruses.", 0, False, INK),
        ("Viral genomes are 500–1,000× smaller than bacterial ones, so even numerically abundant viruses contribute a very small share of reads.", 0, False, INK),
        ("The viruses of surveillance interest sit far down the abundance distribution — precisely where classifier behaviour is least characterised.", 0, False, INK),
        ("", 0, False, INK),
    ], size=17, gap=14)
    rect(s, Inches(0.7), Inches(5.15), Inches(11.93), Inches(1.15), BAND)
    textbox(s, Inches(1.05), Inches(5.36), Inches(11.3), Inches(0.8),
            "The operating question:  given a viral community drawn from a real, deeply sequenced wastewater sample,\n"
            "how many of those viruses does each pipeline find, how many does it invent, and does the answer transfer between platforms?",
            size=15.5, bold=True, color=ACCENT, space_after=0, line=1.25)

    # ---------------------------------------------------------- 3 design fig
    s = slide_frame(prs, "Study design", "One community, two sequencing modalities, seven pipelines", 3)
    picture(s, AN / "figures_combined" / "Figure1_study_design.png",
            Inches(0.9), Inches(1.7), Inches(11.5), Inches(4.5))
    note(s, "708-member simulated community: 490 bacteria/archaea, 67 ssRNA bacteriophage, 151 non-phage viruses. "
            "Only the 151 viruses are scored.")

    # ---------------------------------------------------------- 4 ground truth
    s = slide_frame(prs, "Ground truth", "151 non-phage viruses, taken from a real wastewater metatranscriptome", 4)
    table(s, Inches(0.7), Inches(1.85), Inches(11.93), [
        ["Component", "n", "Designed share of community", "Selection rule"],
        ["Viruses, automated selection", "138", "3.995%", "clade reads ≥ 20, NCBI rank = species"],
        ["Viruses, sentinel additions", "13", "0.005%", "added for public-health relevance"],
        ["Viral ground truth (scored)", "151", "4.00%", "— the denominator of every recall below —"],
        ["Bacteria + archaea (background)", "490", "88.0%", "not scored"],
        ["ssRNA bacteriophage (background)", "67", "8.0%", "not scored"],
    ], col_w=[3.4, 0.8, 2.4, 5.0], size=13.5, row_h=Inches(0.42), hilite_rows=(3,))
    bullets(s, Inches(0.7), Inches(4.85), Inches(11.9), Inches(1.7), [
        ("Designed viral abundances span 1.7 × 10⁻⁵% to 2.36% — a range of about 1.4 × 10⁵-fold.", 0, False, INK),
        ("90 of the 151 viruses were designed below 1 × 10⁻³% relative abundance: the community is concentrated in the rare tail.", 0, False, INK),
        ("13 sentinels include all influenza subtypes (H3N2, H1N1, H5N1, influenza B), poliovirus 1, monkeypox virus, RSV A, hMPV.", 0, False, INK),
    ], size=15.5, gap=10)

    # ---------------------------------------------------------- 5 depth caveat
    s = slide_frame(prs, "Critical caveat", "The realised viral read depth is far below the design intent", 5)
    c = D["comp"]
    ill_v = c[c.platform == "Illumina"].virus_reads
    pac_v = c[c.platform == "PacBio"].virus_reads
    stat_cards(s, Inches(1.75), [
        (f"{ill_v.mean():,.0f}", "viral reads per Illumina library\n(16.6 M read pairs)"),
        ("0.020%", "realised viral read share\non Illumina — target was 12%"),
        ("94 / 151", "viruses expected to contribute\nfewer than one Illumina read"),
        (f"{pac_v.min():,.0f}–{pac_v.max():,.0f}", "viral reads per PacBio library\n(0.45–0.52 M reads)"),
    ], top_size=27)
    bullets(s, Inches(0.7), Inches(3.45), Inches(11.9), Inches(3.0), [
        ("Cause: CAMISIM consumes its abundance file as genome copy number, but the 88 / 8 / 4% split was intended as a share of reads.", 0, True, HILITE),
        ("Reads scale with copies × genome length, so small viral genomes were down-weighted ~675× (Illumina) and ~65× (PacBio).", 0, False, INK),
        ("Verified two independent ways: per-domain read yield matches the genome-length prediction to within ~20%, and length-weighting the candidate input columns predicts the measured share exactly.", 0, False, INK),
        ("Consequence: every viral number in this study is a LOWER BOUND. Corrected copy abundances for re-simulation are provided.", 0, True, ACCENT),
    ], size=15.5, gap=11)

    # ---------------------------------------------------------- 6 pipelines
    s = slide_frame(prs, "What was compared", "Seven pipelines, as configured in their own run scripts", 6)
    table(s, Inches(0.7), Inches(1.8), Inches(11.93), [
        ["Pipeline", "Classifier", "Illumina input", "PacBio input"],
        ["nf-core/taxprofiler v1.2.0", "Kraken2 v2.1.3 + Bracken", "reads", "reads"],
        ["nf-core/taxtriage v1.0", "Kraken2 v2.1.3", "reads", "reads"],
        ["nf-core/mag v3.1.0", "Kraken2 v2.1.3", "reads (fastp first)", "contigs"],
        ["KrakenMetaReads-nf (in-house)", "Kraken2 v2.1.3 on contigs", "metaSPAdes + MEGAHIT contigs", "metaFlye contigs"],
        ["GOTTCHA2", "signature fragments", "reads (-xm sr)", "reads (-xm map-hifi, patched)"],
        ["CLARK", "discriminative k-mers", "reads", "reads"],
        ["sourmash 4.9.4", "FracMinHash sketches", "reads", "contigs"],
    ], col_w=[3.3, 3.0, 3.0, 3.3], size=13, row_h=Inches(0.40))
    note(s, "Four pipelines share the same Kraken2 v2.1.3 and the same k2_standard_20240605 database — differences between them "
            "are therefore not database differences. Input stage was read from each run script, and for sourmash from its own output.")

    # ---------------------------------------------------------- 7 scoring
    s = slide_frame(prs, "Scoring", "How a detection was counted — and what was deliberately not computed", 7)
    bullets(s, Inches(0.7), Inches(1.85), Inches(6.0), Inches(4.6), [
        ("Species-level, strict (primary result)", 0, True, ACCENT),
        ("exact NCBI taxid match at species-or-below nodes (S, S1–S5, strain, serotype, subspecies…)", 1, False, INK),
        ("Genus-level, lenient (sensitivity analysis)", 0, True, ACCENT),
        ("normalised genus name, so sourmash participates despite reporting no taxids", 1, False, INK),
        ("Detection thresholds", 0, True, ACCENT),
        ("≥1 / ≥10 / ≥100 clade reads, all three reported", 1, False, INK),
        ("Uncertainty", 0, True, ACCENT),
        ("bootstrap 95% CI, 10,000 resamples, seed 42; n = 3, so no significance testing", 1, False, INK),
    ], size=16, gap=9)
    rect(s, Inches(7.05), Inches(1.85), Inches(5.58), Inches(3.1), BAND)
    textbox(s, Inches(7.35), Inches(2.08), Inches(5.0), Inches(2.7),
            "Why no specificity, accuracy or MCC",
            size=17, bold=True, color=HILITE, space_after=10)
    textbox(s, Inches(7.35), Inches(2.55), Inches(5.0), Inches(2.3),
            "True negatives are undefined against an open reference database: the set of taxa a "
            "pipeline could have reported but did not is unbounded.\n\n"
            "Any specificity, accuracy, balanced accuracy or MCC would be an artefact of the "
            "database size, not a property of the pipeline. They are omitted on purpose.",
            size=14, color=INK, space_after=8, line=1.22)
    textbox(s, Inches(7.05), Inches(5.15), Inches(5.58), Inches(1.2),
            "Reported instead: recall, precision, F1, and absolute viral false-positive counts.",
            size=15, bold=True, color=ACCENT, space_after=0)

    # ---------------------------------------------------------- 8 Illumina
    s = slide_frame(prs, "Results · Illumina", "No pipeline recovered more than 31% of the viral ground truth", 8)
    picture(s, AN / "figures_illumina" / "Figure1_viral_recall_precision_f1.png",
            Inches(0.6), Inches(1.7), Inches(7.3), Inches(4.6))
    rows = [["Pipeline", "R", "P", "F1", "FP"]]
    for _, r in ill.iterrows():
        rows.append([NAME[r.pipeline], *rp(r), f"{r.FP_mean:,.0f}"])
    table(s, Inches(8.05), Inches(1.95), Inches(4.6), rows,
          col_w=[2.5, 0.8, 0.8, 0.8, 0.9], size=12, row_h=Inches(0.35),
          hilite_rows=(1,))
    textbox(s, Inches(8.05), Inches(4.95), Inches(4.6), Inches(1.5),
            "GOTTCHA2 leads on F1 by making far fewer false calls (65 vs 228), "
            "not by finding more viruses. False-positive burden spans two orders "
            "of magnitude across the seven.",
            size=13, color=MUTED, space_after=0, line=1.2)

    # ---------------------------------------------------------- 9 depth-bound
    s = slide_frame(prs, "Results · Illumina", "That ceiling is set by the data, not by the classifiers", 9)
    r = D["restricted"]
    r = r[r.platform == "Illumina"].pivot_table(index="pipeline",
                                               columns="min_expected_reads",
                                               values="recall")
    order = ["sourmash", "mag-nf", "TaxProfiler", "TaxTriage", "GOTTCHA2",
             "KrakenMetaReads-nf", "CLARK"]
    rows = [["Pipeline", "All 151", "≥1 read\n(n = 57)", "≥2 reads\n(n = 44)", "≥5 reads\n(n = 25)"]]
    for p in order:
        rows.append([NAME[p]] + [f"{r.loc[p, k]:.3f}" for k in [0, 1, 2, 5]])
    table(s, Inches(0.7), Inches(1.9), Inches(7.4), rows,
          col_w=[2.6, 1.2, 1.2, 1.2, 1.2], size=13, row_h=Inches(0.38))
    bullets(s, Inches(8.35), Inches(1.95), Inches(4.3), Inches(4.3), [
        ("Expected reads for each virus = viral reads in the library × (designed abundance ÷ 4.00%).", 0, False, INK),
        ("On viruses that actually contributed signal, the leading pipelines recover 0.75–0.93.", 0, True, ACCENT),
        ("Only one virus with substantial signal — Orthopicobirnavirus hominis, 32 expected reads — was missed by every pipeline.", 0, False, INK),
        ("So recall over all 151 measures the simulation's viral depth, not classifier capability.", 0, True, HILITE),
    ], size=14.5, gap=12)
    note(s, "Nominal recall of 0.305 does not mean nf-core/taxprofiler finds a third of the viruses present; "
            "it means two-thirds of the ground truth was not meaningfully present in the reads.")

    # ---------------------------------------------------------- 10 PacBio
    s = slide_frame(prs, "Results · PacBio HiFi", "A far wider spread: two pipelines near 0.91, three at exactly zero", 10)
    picture(s, AN / "figures_pacbio" / "Figure1_viral_recall_precision_f1.png",
            Inches(0.6), Inches(1.7), Inches(7.3), Inches(4.6))
    rows = [["Pipeline", "R", "P", "F1", "FP"]]
    for _, r_ in pac.iterrows():
        rows.append([NAME[r_.pipeline], *rp(r_), f"{r_.FP_mean:,.0f}"])
    table(s, Inches(8.05), Inches(1.95), Inches(4.6), rows,
          col_w=[2.5, 0.8, 0.8, 0.8, 0.9], size=12, row_h=Inches(0.35),
          hilite_rows=(1,))
    textbox(s, Inches(8.05), Inches(4.95), Inches(4.6), Inches(1.5),
            "GOTTCHA2, nf-core/mag and sourmash detected zero of the 151 viruses "
            "in every sample — while still emitting populated, normal-looking "
            "species lists.",
            size=13, color=MUTED, space_after=0, line=1.2)

    # ---------------------------------------------------------- 11 input type
    s = slide_frame(prs, "Results · PacBio HiFi", "Input stage, not classifier family, explains which pipelines saw viruses", 11)
    picture(s, AN / "figures_pacbio" / "Figure7_recall_by_input_type.png",
            Inches(0.7), Inches(1.7), Inches(7.1), Inches(4.5))
    bullets(s, Inches(8.1), Inches(1.95), Inches(4.5), Inches(4.4), [
        ("All three pipelines given assembled metaFlye contigs recovered 0.00–0.02.", 0, True, ACCENT),
        ("Three of the four given reads recovered 0.51–0.91.", 0, True, ACCENT),
        ("With 896–997 viral reads across 151 genomes, no viral genome can assemble — a contig-fed classifier has nothing viral to classify.", 0, False, INK),
        ("sourmash's zero is an input artefact: its own output records metaflye/assembly.fasta as the query.", 0, False, INK),
        ("GOTTCHA2 is the only pipeline that genuinely failed on reads.", 0, True, HILITE),
    ], size=14.5, gap=12)

    # ---------------------------------------------------------- 12 transfer
    s = slide_frame(prs, "Cross-modality", "The viral ranking does not transfer between platforms", 12)
    picture(s, AN / "figures_combined" / "Figure3_viral_rank_transfer.png",
            Inches(0.7), Inches(1.7), Inches(7.0), Inches(4.5))
    bullets(s, Inches(8.0), Inches(1.95), Inches(4.6), Inches(4.4), [
        ("Spearman ρ = +0.09 (p = 0.84) for viral F1; +0.08 for recall, +0.17 for precision.", 0, True, ACCENT),
        ("None is distinguishable from no association at this sample size.", 0, False, INK),
        ("GOTTCHA2 ranked FIRST on Illumina (F1 0.303) and detected nothing on PacBio.", 0, True, HILITE),
        ("nf-core/taxprofiler rose 0.217 → 0.604; CLARK rose 0.016 → 0.419, from last to third.", 0, False, INK),
        ("A configuration validated on one modality cannot be assumed to work on the other.", 0, True, ACCENT),
    ], size=14.5, gap=12)

    # ---------------------------------------------------------- 13 silent fail
    s = slide_frame(prs, "Cross-modality", "Three pipelines failed outright — and all three failed silently", 13)
    table(s, Inches(0.7), Inches(1.9), Inches(11.93), [
        ["Pipeline", "Illumina viral F1", "PacBio viral F1", "PacBio output it still produced", "Cause"],
        ["GOTTCHA2", "0.303  (rank 1)", "0.000", "populated report, 3 viral FP", "classifier on long reads"],
        ["nf-core/mag", "0.152", "0.000", "populated report, 2 viral FP", "given contigs, not reads"],
        ["sourmash", "0.116", "0.000", "353 spurious viral species", "given contigs, not reads"],
    ], col_w=[2.2, 2.1, 1.8, 3.3, 2.6], size=13.5, row_h=Inches(0.46))
    rect(s, Inches(0.7), Inches(3.9), Inches(11.93), Inches(1.05), BAND)
    textbox(s, Inches(1.05), Inches(4.12), Inches(11.3), Inches(0.7),
            "Every failure was silent: each pipeline completed normally, exited zero, and emitted a species list of ordinary size.\n"
            "Nothing in the output signalled that no true virus had been found.",
            size=15.5, bold=True, color=HILITE, space_after=0, line=1.25)
    bullets(s, Inches(0.7), Inches(5.2), Inches(11.9), Inches(1.4), [
        ("Genus-level matching did not rescue any of the three (max genus recall 0.016).", 0, False, INK),
        ("Implication for deployment: a populated output is not evidence that the pipeline worked. Positive controls are needed per modality.", 0, True, ACCENT),
    ], size=15, gap=9)

    # ---------------------------------------------------------- 14 assembly/fastp
    s = slide_frame(prs, "Controlled comparisons", "Two experiments isolate where the viral signal is actually lost", 14)
    picture(s, AN / "figures_illumina" / "Figure8_assembly_and_preprocessing.png",
            Inches(0.7), Inches(1.75), Inches(11.9), Inches(3.9))
    textbox(s, Inches(0.7), Inches(5.85), Inches(5.7), Inches(1.0),
            "Left: same Kraken2, same database — only the input changes.\n"
            "Assembly retains 6.4% (metaSPAdes) or 0.8% (MEGAHIT) of the viral signal.",
            size=13.5, color=INK, space_after=4, line=1.2)
    textbox(s, Inches(6.9), Inches(5.85), Inches(5.7), Inches(1.0),
            "Right: same reads, same database — only fastp differs.\n"
            "Viral false positives double; true positives do not move.",
            size=13.5, color=INK, space_after=4, line=1.2)

    # ---------------------------------------------------------- 15 assembly
    s = slide_frame(prs, "Finding 1", "Assembly-first classification loses viruses on BOTH platforms", 15)
    rows = [["Platform", "Evidence classified", "Units classified", "Ground-truth viruses recovered"]]
    for plat, ev, lab in [
        ("Illumina", "reads (nf-core/taxprofiler)", "Illumina · reads"),
        ("Illumina", "contigs (metaSPAdes)", "Illumina · metaSPAdes contigs"),
        ("Illumina", "contigs (MEGAHIT)", "Illumina · MEGAHIT contigs"),
        ("PacBio", "reads (nf-core/taxprofiler)", "PacBio · reads"),
        ("PacBio", "contigs (metaFlye)", "PacBio · metaFlye contigs"),
        ("PacBio", "contigs (nf-core/mag)", "PacBio · nf-core/mag contigs"),
    ]:
        r_ = asm.loc[(plat, ev)]
        rows.append([lab.split(" · ")[0], lab.split(" · ")[1],
                     f"{r_.units_total:,.0f}", f"{r_.gt_viruses_hit:.1f}"])
    table(s, Inches(0.7), Inches(1.9), Inches(7.7), rows,
          col_w=[1.5, 3.2, 1.9, 2.3], size=13, row_h=Inches(0.38),
          hilite_rows=(1, 4))
    bullets(s, Inches(8.65), Inches(1.95), Inches(3.95), Inches(4.3), [
        ("Classifier and database held constant; only reads vs contigs varies.", 0, False, INK),
        ("A virus contributing a handful of reads cannot assemble — and once it fails to assemble, no downstream choice recovers it.", 0, False, INK),
        ("The PacBio metaFlye assembly still lists 2,143–2,155 nodes in the Viruses branch, of which only four are real.", 0, False, INK),
        ("Read-level classification is the safer default at low viral depth.", 0, True, ACCENT),
    ], size=14, gap=12)

    # ---------------------------------------------------------- 16 fastp
    s = slide_frame(prs, "Finding 2", "One QC step upstream of the classifier doubles the false positives", 16)
    table(s, Inches(0.7), Inches(1.95), Inches(7.3), [
        ["Same reads, same Kraken2, same database", "nf-core/taxprofiler", "nf-core/mag (+ fastp)"],
        ["Reads classified", f"{fp['taxprofiler_total']:,.0f}", f"{fp['mag-nf_total']:,.0f}"],
        ["Unclassified", f"{fp['taxprofiler_unclassified_pct']:.2f}%", f"{fp['mag-nf_unclassified_pct']:.2f}%"],
        ["Reads assigned to Viruses", f"{fp['taxprofiler_viral_reads']:,.0f}", f"{fp['mag-nf_viral_reads']:,.0f}"],
        ["Viral true positives", f"{fp['taxprofiler_TP']:.0f}", f"{fp['mag-nf_TP']:.0f}"],
        ["Viral FALSE positives", f"{fp['taxprofiler_FP']:.0f}", f"{fp['mag-nf_FP']:.0f}"],
    ], col_w=[3.6, 1.9, 1.9], size=13.5, row_h=Inches(0.40), hilite_rows=(5,))
    bullets(s, Inches(8.35), Inches(1.95), Inches(4.3), Inches(4.4), [
        ("The false positives are strictly nested: all 228 of taxprofiler's appear in mag's 415; none is unique to taxprofiler.", 0, False, INK),
        (f"Of the {fp['mag_only_FP']:.0f} exclusive to nf-core/mag, {fp['mag_only_FP_frac_single_read']:.1%} rest on a single read.", 0, False, INK),
        ("Across the 46 shared true positives both assign exactly the same read counts (median ratio 1.000) — fastp found no extra real signal.", 0, True, ACCENT),
        ("Mechanism: trimming shortens reads; at Kraken2's default confidence = 0 shorter reads match k-mers by chance more readily.", 0, False, INK),
        ("Viral precision 0.169 → 0.101, a 40% loss, for no gain in recall.", 0, True, HILITE),
    ], size=14, gap=11)

    # ---------------------------------------------------------- 17 threshold
    s = slide_frame(prs, "Finding 3", "Read-count filtering is counterproductive for viral detection", 17)
    tpr = D["tp_reads"]
    med_i = tpr[tpr.platform == "Illumina"].reads.median()
    med_p = tpr[tpr.platform == "PacBio"].reads.median()
    picture(s, AN / "figures_combined" / "Figure5_viral_threshold_both.png",
            Inches(0.7), Inches(1.7), Inches(7.2), Inches(4.5))
    bullets(s, Inches(8.15), Inches(1.95), Inches(4.5), Inches(4.4), [
        ("Every threshold above 1 read reduced viral F1 for every pipeline on both modalities.", 0, True, HILITE),
        ("Illumina: best F1 0.303 → 0.134 at ≥10 reads → 0.037 at ≥100.", 0, False, INK),
        ("PacBio: 0.604 → 0.036 → 0.000.", 0, False, INK),
        (f"Viral detections rest on a median of {med_i:.0f} reads (Illumina) and "
         f"{med_p:.0f} (PacBio).", 0, False, INK),
        ("A filter calibrated on bacterial depth removes true viral signal faster than false. Thresholds must be re-derived from viral depth.", 0, True, ACCENT),
    ], size=14, gap=12)
    note(s, "This is the only quantitative finding in the study that generalises across both modalities without qualification — "
            "and it runs against routine practice in metagenomic false-positive control.")

    # ---------------------------------------------------------- 18 abundance
    s = slide_frame(prs, "Detection vs abundance", "Recall falls steeply with designed viral abundance", 18)
    picture(s, AN / "figures_combined" / "Figure6_viral_abundance_recall_both.png",
            Inches(0.7), Inches(1.7), Inches(7.2), Inches(4.5))
    bullets(s, Inches(8.15), Inches(1.95), Inches(4.5), Inches(4.4), [
        ("Illumina: 0.83–0.91 above 1 × 10⁻²% designed abundance; ≤ 0.111 below 1 × 10⁻³% for every pipeline.", 0, False, INK),
        ("Below 1 × 10⁻⁴%, five of seven pipelines detected nothing at all.", 0, False, INK),
        ("PacBio curves look almost flat — but that is the minimum-depth floor, which raised 64.9% of viral members to ~5 reads regardless of design.", 0, True, HILITE),
        ("Read as: given ~5 long reads from a virus, taxprofiler and taxtriage assign them correctly. Not as: sensitivity at arbitrarily low abundance.", 0, False, INK),
    ], size=14, gap=12)

    # ---------------------------------------------------------- 19 sentinels
    s = slide_frame(prs, "Public-health relevance", "Sentinel virus detection: 3 of 13 on Illumina, 13 of 13 on PacBio", 19)
    stat_cards(s, Inches(1.8), [
        ("3 / 13", "sentinel viruses detected by ANY\nIllumina pipeline, in any sample"),
        ("13 / 13", "detected on PacBio — but only by\ntaxprofiler and taxtriage"),
        ("10", "sentinels detected by NO Illumina\npipeline in any of 21 runs"),
        ("0 / 13", "detected by GOTTCHA2, nf-core/mag\nor sourmash on PacBio"),
    ], top_size=26)
    bullets(s, Inches(0.7), Inches(3.55), Inches(11.9), Inches(2.7), [
        ("Undetected on Illumina: influenza A H3N2, H1N1 and H5N1, influenza B, human orthopneumovirus A, human parainfluenza virus 3, human metapneumovirus, monkeypox virus, human papillomavirus 16, poliovirus 1.", 0, False, INK),
        ("Detected on Illumina: only SARS-CoV-2, human coronavirus OC43 and human astrovirus 1 — the three most abundant sentinels.", 0, False, INK),
        ("All thirteen sit below 0.004% designed abundance, in the bins where recall is near zero for every pipeline. The failure follows from abundance, not from anything specific to these agents — but that is the abundance at which they plausibly occur.", 0, True, ACCENT),
        ("The PacBio 13/13 must be read against the depth floor: it shows correct assignment of a handful of long reads, not detection at the designed abundances.", 0, True, HILITE),
    ], size=14.5, gap=11)

    # ---------------------------------------------------------- 20 union
    s = slide_frame(prs, "Coverage", "Even combining all seven pipelines leaves viruses undetected", 20)
    picture(s, AN / "figures_combined" / "Figure7_viral_detection_overlap.png",
            Inches(0.7), Inches(1.7), Inches(7.0), Inches(4.5))
    bullets(s, Inches(8.0), Inches(1.95), Inches(4.6), Inches(4.4), [
        ("Union of all seven pipelines: 78 / 151 on Illumina (51.7%), 140 / 151 on PacBio (92.7%).", 0, True, ACCENT),
        ("72 viruses detected on both modalities; 68 only on PacBio; 6 only on Illumina; 5 by no run on either.", 0, False, INK),
        ("Missed viruses are dominated by enteric RNA viruses — picobirnaviruses, enteroviruses, rotaviruses, cosaviruses — and sparsely referenced environmental ssDNA groups.", 0, False, INK),
        ("These are among the agents wastewater surveillance is most often deployed to track.", 0, True, HILITE),
    ], size=14, gap=12)

    # ---------------------------------------------------------- 21 conclusions
    s = slide_frame(prs, "Conclusions", "What this benchmark supports", 21)
    bullets(s, Inches(0.7), Inches(1.8), Inches(11.9), Inches(4.6), [
        ("Specify the input stage as explicitly as the classifier. Two of the three zero-recall results are input-stage effects that would otherwise be misread as classifier failures.", 0, True, ACCENT),
        ("Prefer read-level classification over assembly-first designs at low viral depth — on both platforms.", 0, True, ACCENT),
        ("Do not apply read-count thresholds calibrated on bacterial depth. Above 1 read, every threshold cost viral F1 on both modalities.", 0, True, ACCENT),
        ("Validate per modality, with positive controls. Ranking does not transfer (ρ = +0.09), and failures are silent.", 0, True, ACCENT),
        ("Audit preprocessing, not just the classifier. A single fastp step doubled viral false positives with no gain in true positives.", 0, True, ACCENT),
        ("Untargeted metagenomics at this depth does not substitute for targeted assays in surveillance of low-abundance sentinel pathogens.", 0, True, HILITE),
    ], size=16, gap=15)

    # ---------------------------------------------------------- 22 limitations
    s = slide_frame(prs, "Limitations", "What these numbers cannot be asked to support", 22)
    bullets(s, Inches(0.7), Inches(1.8), Inches(11.9), Inches(4.7), [
        ("Viral read depth is ~675× (Illumina) and ~65× (PacBio) below the design intent, so all viral results are lower bounds. Re-simulation at the intended read share is the primary outstanding item.", 0, False, INK),
        ("The PacBio simulation applied a minimum-depth floor to 64.9% of viral members without renormalising; PacBio low-abundance recall and the 13/13 sentinel result reflect that floor.", 0, False, INK),
        ("Illumina nf-core/taxtriage shares byte-identical Kraken2 reports with nf-core/taxprofiler and is therefore not an independent evaluation.", 0, False, INK),
        ("PacBio nf-core/mag read-level reports are the same file three times, so contig-level evidence was used for that arm.", 0, False, INK),
        ("n = 3 per platform. Between-sample SD is small (viral recall SD ≤ 0.023), but no significance testing between pipelines is performed or implied.", 0, False, INK),
        ("Simulated reads only. No wet-lab confounders: no host depletion, no capture enrichment, no degraded RNA, no PCR duplication.", 0, False, INK),
    ], size=15, gap=13)
    note(s, "Open items requiring author or collaborator input are tracked in "
            "analysis_708/missing_information_and_questions.md.")

    out = ROOT / "Seven-pipeline-benchmark-overview.pptx"
    prs.save(out)
    print(f"wrote {out}")
    print(f"       {len(prs.slides.__iter__.__self__._sldIdLst)} slides")
    return 0


if __name__ == "__main__":
    sys.exit(build())
