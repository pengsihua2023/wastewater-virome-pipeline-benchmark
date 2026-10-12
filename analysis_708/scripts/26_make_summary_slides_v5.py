#!/usr/bin/env python3
"""Build a concise 10-slide overview deck of the V5 manuscript.

Reuses the drawing primitives and house style (colours, fonts, table/stat-card
helpers) from 23_make_slides.py rather than re-implementing them; imported via
importlib since 23's own filename is not a valid Python module name. Numbers
are pulled from the same verified interim tables as the manuscript wherever a
machine-readable table exists; the two sentinel-virus counts (Table 9 of the
manuscript) are not backed by a saved table and are transcribed verbatim from
the manuscript text, which was itself computed from viral_detection_detail.tsv.

Writes Seven-pipeline-benchmark-summary-V5.pptx at the repository root.
"""
import importlib.util
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
AN = ROOT / "analysis_708"
INT = AN / "interim"

_spec = importlib.util.spec_from_file_location("m23", AN / "scripts" / "23_make_slides.py")
m23 = importlib.util.module_from_spec(_spec)
sys.modules["m23"] = m23
_spec.loader.exec_module(m23)

from pptx import Presentation
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

W, H = m23.W, m23.H
INK, MUTED, ACCENT, HILITE, RULE, BAND, WHITE = (
    m23.INK, m23.MUTED, m23.ACCENT, m23.HILITE, m23.RULE, m23.BAND, m23.WHITE)
textbox, bullets, rect, slide_frame, picture, table, stat_cards, note = (
    m23.textbox, m23.bullets, m23.rect, m23.slide_frame, m23.picture,
    m23.table, m23.stat_cards, m23.note)


def load():
    d = {}
    agg = pd.read_csv(INT / "viral_aggregate_metrics.tsv", sep="\t")
    d["ill"] = agg[agg.platform == "Illumina"].sort_values("f1_mean", ascending=False)
    d["pac"] = agg[agg.platform == "PacBio"].sort_values("f1_mean", ascending=False)
    d["asm"] = pd.read_csv(INT / "assembly_vs_reads.tsv", sep="\t").set_index(["platform", "evidence"])
    d["comp"] = pd.read_csv(INT / "realised_composition.tsv", sep="\t")
    corr = pd.read_csv(INT / "rank_correlations_v4.tsv", sep="\t")
    d["corr"] = corr.set_index(["analysis_set", "metric"])
    return d


def fmt(row):
    r = f"{row.recall_mean:.3f}"
    p = f"{row.precision_mean:.3f}" if row.precision_mean > 0 else "—"
    f1 = f"{row.f1_mean:.3f}"
    return r, p, f1


def build():
    D = load()
    ill, pac = D["ill"], D["pac"]
    asm = D["asm"]
    comp = D["comp"]
    corr = D["corr"]

    prs = Presentation()
    prs.slide_width, prs.slide_height = W, H

    # --------------------------------------------------------------- 1 title
    s = prs.slides.add_slide(prs.slide_layouts[6])
    rect(s, 0, 0, W, Inches(2.75), ACCENT)
    textbox(s, Inches(0.85), Inches(0.55), Inches(11.6), Inches(1.5),
            "Assembly, Not Modality, Determines Whether\nViral Detection Transfers",
            size=32, bold=True, color=WHITE, space_after=4, line=1.12)
    textbox(s, Inches(0.85), Inches(2.15), Inches(11.6), Inches(0.5),
            "Benchmarking seven metagenomic pipelines on paired Illumina / PacBio HiFi wastewater simulations",
            size=16, color=HILITE, space_after=0)
    stat_cards(s, Inches(3.15), [
        ("7", "pipelines\ncompared head to head"),
        ("151", "non-phage viruses\nin the ground truth"),
        ("42", "scored runs\n(7 pipelines × 6 datasets)"),
        ("5 / 7", "pipelines whose viral F1\ntransfers across modalities"),
    ])
    textbox(s, Inches(0.85), Inches(5.55), Inches(11.6), Inches(1.0),
            "Manuscript V5 — Seven-pipeline-comparison-combined_V5.md/.docx. Metrics are computed within the viral "
            "space only; the simulated bacterial/archaeal/phage background contributes to no reported number.",
            size=13.5, color=MUTED, space_after=0, line=1.25)
    note(s, "Simulated community derived from wastewater sample HTP-2024-01-15, Los Angeles (Grimm et al. 2026).", y=Inches(6.85))

    # ----------------------------------------------------------- 2 motivation
    s = slide_frame(prs, "Motivation", "A pipeline validated on one platform is not automatically valid on another", 2)
    bullets(s, Inches(0.7), Inches(1.85), Inches(11.9), Inches(3.4), [
        ("Wastewater surveillance exists to track viruses — SARS-CoV-2, influenza, poliovirus — but almost every classification pipeline in routine use was built and validated on Illumina short reads.", 0, False, INK),
        ("Long-read platforms are now being evaluated for the same purpose. Assemblers, k-mer classifiers and sketching methods all behave differently on long reads, for different mechanistic reasons.", 0, False, INK),
        ("A programme switching platforms needs to know: does the pipeline's viral detection performance carry over, or does it silently fail?", 0, False, INK),
    ], size=17, gap=16)
    rect(s, Inches(0.7), Inches(5.35), Inches(11.93), Inches(1.15), BAND)
    textbox(s, Inches(1.05), Inches(5.56), Inches(11.3), Inches(0.8),
            "This study isolates one variable: is the pipeline given raw reads, or contigs assembled from those reads?\n"
            "That distinction — not the sequencing platform — turns out to predict almost everything below.",
            size=15.5, bold=True, color=ACCENT, space_after=0, line=1.25)

    # -------------------------------------------------------- 3 study design
    s = slide_frame(prs, "Study design", "One 151-virus community, two sequencing modalities, seven pipelines", 3)
    picture(s, AN / "figures_combined" / "Figure1_study_design.png",
            Inches(1.1), Inches(1.75), Inches(11.1), Inches(4.35))
    note(s, "708-member simulated community: 490 bacteria/archaea, 67 ssRNA bacteriophage, 151 non-phage viruses "
            "(the scored ground truth). Both modalities calibrated to the same designed 12% viral+phage read share.")

    # -------------------------------------------------------- 4 data fix
    s = slide_frame(prs, "A data-integrity correction", "The simulated data was regenerated before this analysis was run", 4)
    ill_v = comp[comp.platform == "Illumina"].pct_virus
    pac_v = comp[comp.platform == "PacBio"].pct_virus
    stat_cards(s, Inches(1.85), [
        ("0.009%", "viral reads in the ORIGINAL\nIllumina data (target: 4%)"),
        ("97 / 151", "viruses with ZERO reads\nin the original data"),
        (f"{ill_v.min():.1f}–{ill_v.max():.1f}%", "viral reads in the\nCORRECTED Illumina data"),
        ("151 / 151", "viruses with ≥1 read in\nevery corrected dataset"),
    ], top_size=26)
    bullets(s, Inches(0.7), Inches(3.55), Inches(11.9), Inches(2.9), [
        ("Root cause: the abundance file encodes a read-share target, but CAMISIM (Illumina) consumed it as genome copy number, and the PacBio driver applied an independent, uncorrected depth formula.", 0, False, INK),
        ("Both simulators were corrected and all six datasets regenerated; all seven pipelines were re-run on all six (42 runs). Full derivation in guide-line-regenerat-data.md and Verification-simulated-data-vs-design.md.", 0, False, INK),
        ("Nearly every quantitative result changed as a consequence — the original dataset was mostly measuring data absence, not pipeline behaviour.", 0, True, HILITE),
    ], size=16, gap=14)

    # -------------------------------------------------------- 5 headline table
    s = slide_frame(prs, "Headline result", "GOTTCHA2 leads viral F1 on both platforms; five of seven pipelines transfer", 5)
    rows = [["Pipeline", "Illumina R", "Illumina P", "Illumina F1", "PacBio R", "PacBio P", "PacBio F1"]]
    hilite = []
    order = ["GOTTCHA2", "TaxProfiler", "TaxTriage", "mag-nf", "CLARK", "sourmash", "KrakenMetaReads-nf"]
    for i, pipe in enumerate(order, start=1):
        ri = ill[ill.pipeline == pipe].iloc[0]
        rp = pac[pac.pipeline == pipe].iloc[0]
        r1, p1, f1a = fmt(ri)
        r2, p2, f1b = fmt(rp)
        rows.append([m23.NAME.get(pipe, pipe), r1, p1, f1a, r2, p2, f1b])
        if pipe in ("sourmash", "KrakenMetaReads-nf"):
            hilite.append(i)
    table(s, Inches(0.7), Inches(1.85), Inches(11.93), rows,
          col_w=[3.0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5], size=14, row_h=Inches(0.42),
          hilite_rows=tuple(hilite))
    note(s, "Highlighted rows: the two pipelines given assembled contigs (not reads) on PacBio — the only two whose "
            "viral recall collapses below 0.06.", y=Inches(5.75))

    # -------------------------------------------------------- 6 rank transfer
    s = slide_frame(prs, "Cross-modality transfer", "The ranking transfers well once the input-stage confound is removed", 6)
    picture(s, AN / "figures_combined" / "Figure3_viral_rank_transfer.png",
            Inches(0.7), Inches(1.7), Inches(7.6), Inches(4.6))
    c1 = corr.loc[("All seven pipelines", "viral F1")]
    c2 = corr.loc[("Reads both sides, de-duplicated", "viral F1")]
    stat_cards2_x = Inches(8.55)
    for i, (lab, rho, pval, n) in enumerate([
        ("All 7 pipelines\n(underpowered)", c1.rho, c1.p, int(c1.n)),
        ("Reads on both sides,\nde-duplicated", c2.rho, c2.p, int(c2.n)),
    ]):
        y = Inches(1.9 + i * 2.1)
        rect(s, stat_cards2_x, y, Inches(3.9), Inches(1.85), BAND)
        textbox(s, stat_cards2_x + Inches(0.2), y + Inches(0.15), Inches(3.5), Inches(0.9),
                f"ρ = {rho:+.2f}", size=32, bold=True, color=ACCENT, space_after=0)
        textbox(s, stat_cards2_x + Inches(0.2), y + Inches(0.95), Inches(3.5), Inches(0.85),
                f"{lab}  (n = {n}, p = {pval:.2f})", size=13, color=MUTED, space_after=0, line=1.2)
    note(s, "Restricting to the pipelines given reads on both modalities more than doubles the correlation "
            "(Table 5 of the manuscript).")

    # -------------------------------------------------------- 7 assembly loss
    s = slide_frame(prs, "Why the two collapses happen", "Assembly, not the sequencing platform, removes viral signal", 7)
    rows = [["Modality", "Evidence classified", "Viruses recovered (of 151)", "% of read-level"]]
    ill_base = asm.loc[("Illumina", "reads (nf-core/mag)")].gt_viruses_hit
    pac_base = asm.loc[("PacBio", "reads (nf-core/taxprofiler)")].gt_viruses_hit
    keys = [("Illumina", "reads (nf-core/mag)", "reads (nf-core/mag)", ill_base),
            ("Illumina", "contigs (metaSPAdes)", "metaSPAdes contigs", ill_base),
            ("Illumina", "contigs (MEGAHIT)", "MEGAHIT contigs", ill_base),
            ("PacBio", "reads (nf-core/taxprofiler)", "reads (nf-core/taxprofiler, nf-core/mag)ᵃ", pac_base),
            ("PacBio", "contigs (metaFlye)", "metaFlye contigs", pac_base),
            ("PacBio", "contigs (nf-core/mag)", "nf-core/mag's own contig pipeline", pac_base)]
    for plat, key, label, base in keys:
        row = asm.loc[(plat, key)]
        rows.append([plat, label, f"{row.gt_viruses_hit:.1f}", f"{row.gt_viruses_hit / base * 100:.0f}%"])
    table(s, Inches(0.7), Inches(1.85), Inches(11.93), rows,
          col_w=[1.6, 4.4, 3.0, 2.0], size=13.5, row_h=Inches(0.4), hilite_rows=(5, 6))
    textbox(s, Inches(0.7), Inches(4.75), Inches(11.9), Inches(0.4),
            "ᵃ Identical on PacBio: same classifier, database and reads, no preprocessing difference recorded.",
            size=11.5, italic=True, color=MUTED, space_after=0)
    bullets(s, Inches(0.7), Inches(5.3), Inches(11.9), Inches(1.4), [
        ("This locates both PacBio collapses in the input stage: no long-read classifier was ever asked to classify long reads for either failure.", 0, True, HILITE),
    ], size=16, gap=0)

    # -------------------------------------------------------- 8 F1 vs sentinel
    s = slide_frame(prs, "A counter-intuitive finding", "Aggregate viral F1 does not predict sentinel-pathogen coverage", 8)
    stat_cards(s, Inches(1.85), [
        ("#1", "GOTTCHA2's viral F1 rank\non BOTH platforms"),
        ("1 / 13", "sentinel viruses GOTTCHA2\ndetects on Illumina"),
        ("0 / 13", "sentinel viruses GOTTCHA2\ndetects on PacBio"),
        ("13 / 13", "sentinels detected by\ntaxprofiler / taxtriage / mag"),
    ], top_size=28)
    bullets(s, Inches(0.7), Inches(3.6), Inches(11.9), Inches(2.7), [
        ("The 13 sentinels are public-health-priority viruses: SARS-CoV-2, influenza A/B (incl. H5N1), RSV, hMPV, monkeypox, poliovirus, HPV16 and others.", 0, False, INK),
        ("GOTTCHA2's signature database evidently does not cover most of these specific taxa, even though it classifies the 151-member community overall more accurately than any other pipeline.", 0, False, INK),
        ("A surveillance programme selecting a pipeline purely on aggregate F1 could systematically miss the exact pathogens the programme exists to track.", 0, True, HILITE),
    ], size=16, gap=14)

    # -------------------------------------------------------- 9 recommendations
    s = slide_frame(prs, "Recommendations", "What this means for a programme choosing or migrating a pipeline", 9)
    bullets(s, Inches(0.7), Inches(1.8), Inches(11.9), Inches(4.95), [
        ("Specify and validate the input stage — reads vs. assembled contigs — as part of the pipeline spec, not only the classifier and database.", 0, True, ACCENT),
        ("Prefer read-level classification over assembly-first classification for viral surveillance, especially on long-read platforms.", 0, True, ACCENT),
        ("Verify detection of the specific sentinel pathogens a programme tracks, independently of a pipeline's aggregate F1.", 0, True, ACCENT),
        ("Run viral positive controls in production — both PacBio failures here were completely silent.", 0, True, ACCENT),
        ("Re-measure any read-count filtering threshold empirically; its effect reversed sign for 4 of 7 pipelines once the data were corrected.", 0, True, ACCENT),
        ("Re-benchmark whenever migrating a validated pipeline to a new sequencing modality — even pipelines that transfer well should be re-verified.", 0, True, ACCENT),
    ], size=16, gap=13)

    # -------------------------------------------------------- 10 closing
    s = prs.slides.add_slide(prs.slide_layouts[6])
    rect(s, 0, 0, W, H, ACCENT)
    textbox(s, Inches(0.85), Inches(2.3), Inches(11.6), Inches(1.7),
            "Viral detection is portable across sequencing\nmodalities — provided the pipeline is given\nreads, not contigs.",
            size=30, bold=True, color=WHITE, space_after=6, line=1.25)
    textbox(s, Inches(0.85), Inches(4.35), Inches(11.6), Inches(0.6),
            "The input stage, not the classifier or the platform, is what should be specified and re-validated.",
            size=17, color=HILITE, space_after=0)
    textbox(s, Inches(0.85), Inches(6.55), Inches(11.6), Inches(0.5),
            "Full manuscript: Seven-pipeline-comparison-combined_V5.md / .docx (EN + CH)",
            size=13, color=BAND, space_after=0)

    out = ROOT / "Seven-pipeline-benchmark-summary-V5.pptx"
    prs.save(out)
    print("wrote", out, "—", len(prs.slides.__iter__.__self__._sldIdLst), "slides")
    return out


if __name__ == "__main__":
    build()
    sys.exit(0)
