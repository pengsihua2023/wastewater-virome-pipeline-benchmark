#!/usr/bin/env python3
"""25-slide group-meeting deck summarising manuscript V5.9.

Style helpers come from 23_make_slides.py (via importlib). Every table on the
slides is read from Seven-pipeline-comparison-combined_V5.9.docx so that the
deck cannot drift from the manuscript; coverage numbers come from
analysis_708/interim/coverage/. Writes Seven-pipeline-benchmark-group-meeting-V5.9.pptx.
"""
import importlib.util
import sys
from pathlib import Path

import pandas as pd
from docx import Document

ROOT = Path(__file__).resolve().parents[2]
AN = ROOT / "analysis_708"
FIG = AN / "figures_combined"
COV = AN / "interim" / "coverage"

_spec = importlib.util.spec_from_file_location("m23", AN / "scripts" / "23_make_slides.py")
m23 = importlib.util.module_from_spec(_spec)
sys.modules["m23"] = m23
_spec.loader.exec_module(m23)

from pptx import Presentation
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches

W, H = m23.W, m23.H
INK, MUTED, ACCENT, HILITE, RULE, BAND, WHITE = (
    m23.INK, m23.MUTED, m23.ACCENT, m23.HILITE, m23.RULE, m23.BAND, m23.WHITE)
textbox, bullets, rect, slide_frame, picture, table, stat_cards, note = (
    m23.textbox, m23.bullets, m23.rect, m23.slide_frame, m23.picture,
    m23.table, m23.stat_cards, m23.note)

MS = Document(ROOT / "Seven-pipeline-comparison-combined_V5.9.docx")


def tbl(i):
    return [[c.text.replace("\n", " ") for c in r.cells] for r in MS.tables[i].rows]


T2, T3, T4, T5, T6, T8, T9 = tbl(1), tbl(2), tbl(3), tbl(4), tbl(5), tbl(7), tbl(8)
CLS = pd.read_csv(COV / "class_summary.tsv", sep="\t")
SENT = pd.read_csv(COV / "sentinel_coverage.tsv", sep="\t")
FL = pd.read_csv(COV / "floor_summary.tsv", sep="\t")


def notes(s, text):
    s.notes_slide.notes_text_frame.text = text


def build():
    prs = Presentation()
    prs.slide_width, prs.slide_height = W, H
    n = 0

    def frame(kicker, title):
        nonlocal n
        n += 1
        return slide_frame(prs, kicker, title, n)

    # 1 ------------------------------------------------------------ title
    n += 1
    s = prs.slides.add_slide(prs.slide_layouts[6])
    rect(s, 0, 0, W, Inches(2.9), ACCENT)
    textbox(s, Inches(0.85), Inches(0.5), Inches(11.6), Inches(1.7),
            "Benchmarking Seven Metagenomic Pipelines for Viral Detection Using Simulated Short- and Long-Read Wastewater Data",
            size=29, bold=True, color=WHITE, space_after=4, line=1.12)
    textbox(s, Inches(0.85), Inches(2.3), Inches(11.6), Inches(0.5),
            "Does viral detection transfer from Illumina to PacBio HiFi? The input stage decides.",
            size=16, color=HILITE, space_after=0)
    stat_cards(s, Inches(3.35), [
        ("7", "pipelines compared\nhead to head"),
        ("151", "non-phage viruses\nin the ground truth"),
        ("42", "scored runs\n(7 pipelines × 6 datasets)"),
        ("5 / 7", "pipelines whose viral\nperformance transfers"),
    ])
    textbox(s, Inches(0.85), Inches(5.6), Inches(11.6), Inches(0.9),
            "Sihua Peng, ..., Justin Bahl\nDepartment of Infectious Diseases, College of Veterinary Medicine, University of Georgia\n"
            "Group meeting — based on manuscript V5.9 (target journal: Water Research)",
            size=14, color=MUTED, space_after=2, line=1.2)
    notes(s, "Introduce the question: when a wastewater programme moves from Illumina to long reads, does its validated pipeline still find viruses?")

    # 2 ------------------------------------------------------------ outline
    s = frame("Outline", "What I will cover")
    items = [
        ("1.  Why this question matters: wastewater surveillance and pipeline choice", 0, False, INK),
        ("2.  Study design: a 151-virus simulated community on two modalities", 0, False, INK),
        ("3.  Seven pipelines and how they were scored", 0, False, INK),
        ("4.  Results: performance, transfer, the two silent collapses, assembly loss", 0, False, INK),
        ("5.  Abundance, thresholds and public-health sentinel viruses", 0, False, INK),
        ("6.  Limitations, recommendations and take-home messages", 0, False, INK),
    ]
    bullets(s, Inches(0.9), Inches(1.95), Inches(11.4), Inches(4.4), items, size=21, gap=20)
    notes(s, "Roughly 25 minutes: 5 on background and design, 12 on results, 5 on limitations and discussion.")

    # 3 ------------------------------------------------------------ motivation
    s = frame("Background", "Surveillance pipelines were validated on short reads")
    bullets(s, Inches(0.7), Inches(1.85), Inches(11.9), Inches(3.3), [
        ("Wastewater-based epidemiology (WBE) tracks viruses at population scale without depending on people seeking care — enteroviruses, influenza (human and avian), SARS-CoV-2.", 0, False, INK),
        ("Untargeted metagenomics removes the fixed target panel and can describe the whole virome, but the classification pipeline decides which viruses are reported and at what false-positive cost.", 0, False, INK),
        ("Both errors are costly: a missed pathogen defeats surveillance; a spurious high-consequence detection triggers investigation and public alarm.", 0, False, INK),
    ], size=17, gap=15)
    rect(s, Inches(0.7), Inches(5.3), Inches(11.93), Inches(1.15), BAND)
    textbox(s, Inches(1.05), Inches(5.5), Inches(11.3), Inches(0.8),
            "Almost every pipeline in routine WBE use was developed and validated on Illumina short reads.\n"
            "Long-read sequencing is now being considered — is validation portable?",
            size=16, bold=True, color=ACCENT, space_after=0, line=1.25)
    notes(s, "Key references in the manuscript: Bubba 2023, Xiu & Yin 2024, Tisza 2023, Worp 2025 for WBE and metagenomics.")

    # 4 ------------------------------------------------------------ why doubt
    s = frame("The transferability question", "Mechanistic reasons to doubt pipeline transfer")
    cards = [
        ("Assemblers", "Behave differently on long reads and may fail to produce classifiable contigs at low viral depth (e.g. metaFlye)."),
        ("Signature classifiers", "Built around fixed-length signature fragments; may not accept long reads without modification (e.g. GOTTCHA2)."),
        ("Sketching methods", "Subsample k-mers; may keep too few hashes from a small viral read pool (e.g. sourmash)."),
    ]
    for i, (h, t) in enumerate(cards):
        x = Inches(0.7 + i * 4.05)
        rect(s, x, Inches(1.95), Inches(3.85), Inches(3.2), BAND)
        textbox(s, x + Inches(0.25), Inches(2.15), Inches(3.4), Inches(0.5), h, size=20, bold=True, color=ACCENT, space_after=0)
        textbox(s, x + Inches(0.25), Inches(2.85), Inches(3.4), Inches(2.2), t, size=15.5, color=INK, space_after=0, line=1.25)
    textbox(s, Inches(0.7), Inches(5.5), Inches(11.9), Inches(1.0),
            "Each mechanism acts on a different subset of pipelines, so the viral ranking could in principle be preserved, compressed or scrambled across modalities.",
            size=16, bold=True, color=HILITE, space_after=0, line=1.25)
    notes(s, "Existing long-read benchmarks (Portik 2022, Govender & Eyre 2022) compare classifiers; to our knowledge none separates the input stage in a wastewater viral background.")

    # 5 ------------------------------------------------------------ objectives
    s = frame("Objectives", "Five questions, one controlled comparison")
    bullets(s, Inches(0.8), Inches(1.9), Inches(11.7), Inches(4.6), [
        ("Quantify species-level viral recall, precision, F1 and false-positive counts for each pipeline on each modality.", 0, False, INK),
        ("Measure how well the viral ranking on one modality predicts the ranking on the other.", 0, False, INK),
        ("Identify which pipelines fail on which modality, and in what mode.", 0, False, INK),
        ("Determine which findings are modality-specific and which generalise.", 0, False, INK),
        ("Characterise viral detection blind spots shared across pipelines and modalities.", 0, False, INK),
    ], size=19, gap=18)
    note(s, "Terminology: modality = type of sequencing data (short-read vs long-read); platform = the instrument and its chemistry.")
    notes(s, "Simulation gives an exactly known ground truth, which is what lets us attribute failures to a specific workflow stage.")

    # 6 ------------------------------------------------------------ design
    s = frame("Study design", "One 151-virus community, two modalities, seven pipelines")
    picture(s, FIG / "Figure1_study_design.png", Inches(1.1), Inches(1.7), Inches(11.1), Inches(4.5))
    note(s, "Ground truth derived from the published Los Angeles wastewater profile HTP-2024-01-15 (Grimm et al. 2026); "
            "Illumina by CAMISIM, PacBio HiFi by PBSIM3 + ccs.")
    notes(s, "Walk left to right: source sample, ground truth, simulation, seven pipelines, viral-only scoring.")

    # 7 ------------------------------------------------------------ community
    s = frame("Ground truth", "A realistic wastewater background around 151 target viruses")
    stat_cards(s, Inches(1.85), [
        ("151", "non-phage viral species\n(scored ground truth)"),
        ("13", "public-health\nsentinel viruses"),
        ("490", "bacterial and archaeal\nbackground species"),
        ("67", "ssRNA bacteriophage\nbackground taxa"),
    ])
    bullets(s, Inches(0.7), Inches(3.6), Inches(11.9), Inches(3.0), [
        ("Viruses were selected from the Kraken2 clade-read profile (species rank, ≥ 20 reads); the 13 sentinels were added manually.", 0, False, INK),
        ("Sentinels: SARS-CoV-2, OC43, influenza A (H3N2, H1N1, H5N1) and B, RSV A, HPIV3, hMPV, monkeypox, astrovirus 1, HPV16, poliovirus 1.", 0, False, INK),
        ("Designed abundances span 1.7 × 10⁻⁵% to 2.36%; 90 of 151 viruses are below 1 × 10⁻³%. Background taxa contribute to no metric.", 0, False, INK),
    ], size=16, gap=13)
    notes(s, "The bacterial and phage background is there so classifiers face a realistic non-viral read pool.")

    # 8 ------------------------------------------------------------ simulation
    s = frame("Read simulation", "Both modalities driven by the same abundance files")
    rows = [["", "Illumina (short-read)", "PacBio HiFi (long-read)"],
            ["Simulator", "CAMISIM 2.0", "PBSIM3 (qshmm, RS II model, 10 passes) + ccs"],
            ["Reads per sample", "16.5 M read pairs (≈ 5.0 Gbp)", "≈ 479 k reads (≈ 6.1 Gbp)"],
            ["Read length", "2 × 150 bp", "median ≈ 14 kb; mean predicted accuracy 99.83%"],
            ["Viral (+ phage) reads", "1.87–1.90 M (11.3–11.5%)", "55.8–56.4 k (11.7–11.8%)"],
            ["Design target", "12% viral + phage read share", "12% viral + phage read share"]]
    table(s, Inches(0.7), Inches(1.85), Inches(11.93), rows, col_w=[2.6, 4.2, 5.1], size=14, row_h=Inches(0.5))
    bullets(s, Inches(0.7), Inches(5.1), Inches(11.9), Inches(1.5), [
        ("Abundance-file column = read-share target (verified by exact identity, r = 1.000000 against the source sample).", 0, False, INK),
        ("Illumina: copy number = target share ÷ genome length; PacBio: per-contig depth with globally unique read names.", 0, False, INK),
    ], size=15, gap=8)
    notes(s, "ccs parameters: --min-passes 3 --min-rq 0.99; ccs version still to be filled in the manuscript.")

    # 9 ------------------------------------------------------------ composition / depth
    s = frame("What the pipelines actually received", "Similar read shares, but different bases and depth")
    m = CLS.groupby(["platform", "klass"])[["read_share_pct", "base_share_pct"]].mean().round(2)
    rows = [["Class", "Illumina reads %", "Illumina bases %", "PacBio reads %", "PacBio bases %"]]
    for k, lab in [("non-phage viruses (151 targets)", "151 target viruses"), ("bacteriophages (67)", "67 bacteriophages"),
                   ("bacteria and archaea (490)", "490 bacteria / archaea")]:
        a = m.loc[(m.index.get_level_values(0)[0], k)]
        b = m.loc[("PacBio HiFi", k)]
        rows.append([lab, f"{a.read_share_pct:.2f}", f"{a.base_share_pct:.2f}", f"{b.read_share_pct:.2f}", f"{b.base_share_pct:.2f}"])
    table(s, Inches(0.7), Inches(1.85), Inches(11.93), rows, col_w=[3.4, 2, 2, 2, 2], size=14, row_h=Inches(0.45), hilite_rows=(1,))
    stat_cards(s, Inches(4.05), [
        ("94 vs 3", "median reads per virus\nIllumina pairs vs PacBio reads"),
        ("34 vs 25", "viruses with mean depth < 1×\n(Illumina vs PacBio)"),
        ("59–61", "PacBio viruses whose target was\nraised to a single read"),
        ("2.1% vs 4.0%", "target-virus share of bases\n(PacBio vs Illumina)"),
    ], top_size=24)
    note(s, "Every virus is guaranteed ≥ 1 read by design, so recall is conditional on presence — not a limit of detection.", y=Inches(5.75))
    notes(s, "Be upfront: read-count share is matched but PacBio delivers about half the target-virus bases; 50–54 viruses per dataset have a proportional expectation below one read.")

    # 10 ----------------------------------------------------------- pipelines
    s = frame("Seven pipelines", "Different principles, and different inputs on PacBio")
    rows = [[c for c in (r[0], r[2], r[3], r[4], r[5])] for r in T2]
    rows[0] = ["Pipeline", "Classifier", "Database", "Illumina evidence", "PacBio evidence"]
    for r in rows[1:]:
        for j in range(len(r)):
            r[j] = r[j].replace("\u1d43", "").replace("\u1d47", "").replace("\u1d9c", "").replace("  ", " ").replace(" ,", ",").replace(" )", ")").strip()
        r[2] = r[2].replace("k2_standard_20240605", "Kraken2 standard").replace("GOTTCHA2 species signature DB", "GOTTCHA2 signatures")
        r[2] = r[2].replace("GTDB reps RS226 + NCBI viruses 2025.01", "GTDB + NCBI viruses").replace("CLARK_DB, species targets", "CLARK DB")
    table(s, Inches(0.5), Inches(1.75), Inches(12.33), rows, col_w=[2.2, 2.6, 2.1, 2.6, 2.8], size=12, row_h=Inches(0.5),
          hilite_rows=(6, 7))
    note(s, "Highlighted: the two pipelines classify assembled contigs, not reads, on PacBio (KrakenMetaReads-nf by design; sourmash's PacBio arm). "
            "The three nf-core pipelines share Kraken2 and its database.", y=Inches(5.95))
    notes(s, "Table 2 of the manuscript. Five of seven pipelines see reads on both modalities; two see contigs on at least the PacBio side.")

    # 11 ----------------------------------------------------------- scoring
    s = frame("Scoring", "Viral-only, exact-taxid scoring at three read-count thresholds")
    bullets(s, Inches(0.7), Inches(1.85), Inches(11.9), Inches(4.7), [
        ("True positive: a ground-truth virus whose NCBI taxid appears at a species-or-below node with at least the threshold number of clade reads (1, 10 or 100).", 0, False, INK),
        ("False positive: a reported species-or-below node inside the Viruses branch that is not among the 151 ground-truth members.", 0, False, INK),
        ("A parent or child taxid of a truth virus counts as a miss and a false positive; a genus-level analysis (124 of 151 resolvable) checks this strictness.", 0, False, INK),
        ("Recall = TP/151, precision = TP/(TP+FP), F1 = 2PR/(P+R). Background taxa are ignored; true negatives are undefined, so no specificity or MCC.", 0, False, INK),
        ("GOTTCHA2 is scored on its unfiltered .full.tsv, as for the other six pipelines; sourmash is matched by name (recall is a lower bound).", 0, False, INK),
    ], size=16, gap=14)
    notes(s, "Three replicates per cell, so no significance testing between pipelines; means and SDs reported.")

    # 12 ----------------------------------------------------------- headline table
    s = frame("Headline result", "GOTTCHA2 leads viral F1 on both modalities")
    rows = [["Pipeline", "Illumina R", "Illumina P", "Illumina F1", "PacBio R", "PacBio P", "PacBio F1"]]
    for r in T3[1:]:
        rows.append([r[0]] + [c.split(" ±")[0] for c in (r[1], r[2], r[3], r[5], r[6], r[7])])
    hl = tuple(i for i, r in enumerate(rows) if r[0] in ("sourmash", "KrakenMetaReads-nf"))
    table(s, Inches(0.7), Inches(1.85), Inches(11.93), rows, col_w=[3.0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5], size=14,
          row_h=Inches(0.44), hilite_rows=hl)
    note(s, "Mean of three samples. Highlighted rows: contigs on PacBio — the only two with recall below 0.06 there. "
            "GOTTCHA2's PacBio result uses a locally patched binary.", y=Inches(5.75))
    notes(s, "Precision, not recall, is the discriminating axis on Illumina: five pipelines exceed 0.90 recall; precision ranges 0.026 to 0.628.")

    # 13 ----------------------------------------------------------- figure 2
    s = frame("Performance on both modalities", "Recall, precision and F1 for all seven pipelines")
    picture(s, FIG / "Figure2_viral_performance_both.png", Inches(0.9), Inches(1.65), Inches(11.5), Inches(4.85))
    notes(s, "Bars are means of three samples; orange points are individual samples. Between-sample SD of recall is at most 0.008.")

    # 14 ----------------------------------------------------------- transfer table
    s = frame("Cross-modality transfer", "Five pipelines transfer; two contig-fed ones collapse")
    rows = [["Pipeline", "Evidence (Ill → Pac)", "F1 rank Ill → Pac", "Δ recall", "Δ F1", "Direction"]]
    for r in T4[1:]:
        rows.append([r[0], r[1], f"{r[2]} → {r[3]}", r[4], r[5], r[6]])
    hl = tuple(i for i, r in enumerate(rows) if "collapsed" in r[5])
    table(s, Inches(0.7), Inches(1.85), Inches(11.93), rows, col_w=[2.8, 2.4, 2.0, 1.4, 1.4, 2.6], size=13.5,
          row_h=Inches(0.46), hilite_rows=hl)
    bullets(s, Inches(0.7), Inches(5.9), Inches(11.9), Inches(0.8), [
        ("Read-fed pipelines: recall changes by at most 0.08; F1 falls by at most 0.05 (GOTTCHA2) or rises with PacBio's higher precision.", 0, False, INK),
    ], size=15, gap=0)
    notes(s, "The variable that predicts collapse is the PacBio input stage, not the Illumina baseline.")

    # 15 ----------------------------------------------------------- rank correlation
    s = frame("Rank transfer", "Ranking transfers once the input stage is matched")
    picture(s, FIG / "Figure3_viral_rank_transfer.png", Inches(0.7), Inches(1.65), Inches(6.6), Inches(4.85))
    rows = [["Analysis set", "n", "ρ (recall)", "p (exact)"]]
    for r in T5[1:]:
        rows.append([r[0].replace("Illumina duplicate removed", "Illumina duplicate removed"), r[1],
                     r[2].split(" (")[0], r[2].split("(")[1].rstrip(")")])
    table(s, Inches(7.55), Inches(1.9), Inches(5.1), rows, col_w=[3.3, 0.6, 1.1, 1.1], size=12, row_h=Inches(0.5),
          hilite_rows=(3,))
    textbox(s, Inches(7.55), Inches(4.75), Inches(5.1), Inches(1.7),
            "ρ = 0.34–0.65 (all seven)  →  0.76–0.95 (reads on both sides). With n = 5 none reaches conventional significance "
            "(exact permutation p = 0.067) — read as a consistent trend. The three nf-core pipelines share a classifier and database.",
            size=13, color=MUTED, space_after=0, line=1.2)
    notes(s, "Earlier p = 0.03 was an asymptotic approximation; exact permutation test gives 0.067 for the five read-fed pipelines.")

    # 16 ----------------------------------------------------------- silent failures
    s = frame("The two collapses", "Both failures were silent — no error, plausible reports")
    picture(s, FIG / "Figure4_viral_TP_FN_FP_both.png", Inches(0.7), Inches(1.65), Inches(7.6), Inches(4.7))
    bullets(s, Inches(8.5), Inches(1.9), Inches(4.2), Inches(4.6), [
        ("sourmash F1 0.307 → 0.018; KrakenMetaReads-nf 0.581 → 0.015 (recall below 0.06 on PacBio).", 0, False, INK),
        ("Both classify metaFlye contigs on PacBio; every other pipeline sees reads.", 0, False, INK),
        ("sourmash's PacBio report lists ~330 viral species, almost all wrong (precision 0.013).", 0, False, INK),
        ("Genus-level matching does not rescue them (0.081 and 0.086).", 0, False, INK),
    ], size=14, gap=10)
    notes(s, "No error, no warning: each pipeline finished and produced a populated report — hence the recommendation to run positive controls.")

    # 17 ----------------------------------------------------------- assembly loss
    s = frame("Why the collapses happen", "Assembly, not the platform, removes viral signal")
    rows = [["Modality", "Evidence classified (Kraken2, same DB)", "Viruses recovered", "% of read-level"]]
    rows += [[r[0], r[1], r[2], r[3]] for r in T6[1:]]
    table(s, Inches(0.7), Inches(1.85), Inches(11.93), rows, col_w=[1.6, 5.4, 2.6, 2.3], size=14, row_h=Inches(0.44),
          hilite_rows=(6, 7))
    bullets(s, Inches(0.7), Inches(5.65), Inches(11.9), Inches(1.2), [
        ("Assembly removes 94–95% of viral signal on PacBio but 4–40% on Illumina, with classifier and database held constant.", 0, True, HILITE),
        ("Direct for the Kraken2 path; for sourmash (never run on PacBio reads) this is an inference.", 0, False, INK),
    ], size=15.5, gap=8)
    notes(s, "Table 6. Same three samples, same Kraken2 version and database in every row: the only variable is the input.")

    # 18 ----------------------------------------------------------- abundance
    s = frame("Abundance", "Recall is flat across designed abundance")
    picture(s, FIG / "Figure5_viral_abundance_recall_both.png", Inches(0.7), Inches(1.65), Inches(11.9), Inches(3.9))
    bullets(s, Inches(0.7), Inches(5.6), Inches(11.9), Inches(1.1), [
        ("Viruses with a proportional expectation below one PacBio read were detected at least as often as the rest (83% vs 77%): the one-read minimum does not explain the flat profile.", 0, False, INK),
    ], size=14.5, gap=0)
    notes(s, "But a >=1-read criterion is insensitive to depth by construction. Contig-fed PacBio arms sit near zero in every bin.")

    # 19 ----------------------------------------------------------- thresholds
    s = frame("Read-count filtering", "Whether a threshold helps depends on pipeline and modality")
    picture(s, FIG / "Figure6_viral_threshold_both.png", Inches(0.7), Inches(1.65), Inches(11.9), Inches(3.8))
    bullets(s, Inches(0.7), Inches(5.55), Inches(11.9), Inches(1.2), [
        ("Illumina: four of seven pipelines gain F1 at ≥ 10 reads (a ≥ 10 filter removes far more false positives than true positives); on PacBio every pipeline loses F1 monotonically.", 0, False, INK),
        ("Measure the effect on your own library and pipeline — do not assume a direction.", 0, True, HILITE),
    ], size=14.5, gap=6)
    notes(s, "Table 8 in the manuscript lists F1 at >=1, >=10 and >=100.")

    # 20 ----------------------------------------------------------- sentinel headline
    s = frame("Public-health sentinels", "Aggregate F1 does not predict sentinel-pathogen coverage")
    stat_cards(s, Inches(1.85), [
        ("#1", "GOTTCHA2 viral F1 rank\non both modalities"),
        ("1 / 13", "sentinels GOTTCHA2\ndetects on Illumina"),
        ("0 / 13", "sentinels GOTTCHA2\ndetects on PacBio"),
        ("13 / 13", "detected by taxprofiler,\ntaxtriage and mag"),
    ], top_size=28)
    bullets(s, Inches(0.7), Inches(3.6), Inches(11.9), Inches(3.0), [
        ("All 13 sentinels are detected on both modalities by at least one pipeline — but coverage is uneven: three pipelines detect all 13; the rest detect at most two (KrakenMetaReads-nf: 12 on Illumina).", 0, False, INK),
        ("GOTTCHA2's signature database may not cover these specific taxa (not assessed), even though it classifies the 151-member community best overall.", 0, False, INK),
        ("The three pipelines that detect all 13 share Kraken2 and its database — agreement, not independent cross-validation.", 0, False, INK),
        ("A programme selecting on aggregate F1 alone could miss the pathogens it exists to track.", 0, True, HILITE),
    ], size=15.5, gap=11)
    notes(s, "Table 9. Illumina: KrakenMetaReads-nf detects 12 of 13; CLARK, GOTTCHA2 and sourmash detect influenza B or monkeypox only.")

    # 21 ----------------------------------------------------------- sentinel depth
    s = frame("An honest caveat", "Sentinels were present at only a few PacBio reads")
    rows = [["Sentinel virus", "Illumina read pairs", "Illumina depth", "PacBio reads", "PacBio depth"]]
    for _, r in SENT.iterrows():
        nm = r["name"].replace("Influenza A virus H5N1 (A/goose/Guangdong/1/1996)", "Influenza A H5N1").replace("Influenza A virus", "Influenza A")
        rows.append([nm, f"{r.illumina_reads_mean:.0f}", f"{r.illumina_depth_x:.2f}×", f"{r.pacbio_reads_mean:.1f}", f"{r.pacbio_depth_x:.2f}×"])
    table(s, Inches(0.7), Inches(1.65), Inches(7.5), rows, col_w=[3.4, 1.5, 1.4, 1.3, 1.3], size=11, row_h=Inches(0.33))
    bullets(s, Inches(8.45), Inches(1.75), Inches(4.2), Inches(4.8), [
        ("Nine of 13 sentinels have a proportional expectation below one PacBio read; they are represented by 1–8 reads (influenza: one per segment).", 0, False, INK),
        ("PacBio detection therefore shows a pipeline can classify a single long read — not that it would find a real sentinel at that abundance.", 0, True, HILITE),
    ], size=14, gap=12)
    notes(s, "Supplementary Table S3. Mean of three samples; depth = bases / genome length.")

    # 22 ----------------------------------------------------------- overlap
    s = frame("Detection overlap", "141 of 151 viruses detected on both modalities")
    picture(s, FIG / "Figure7_viral_detection_overlap.png", Inches(0.7), Inches(1.65), Inches(7.2), Inches(4.85))
    bullets(s, Inches(8.2), Inches(1.9), Inches(4.4), Inches(4.5), [
        ("Union of all seven pipelines: 151 / 151 on Illumina (100%), 141 / 151 on PacBio (93.4%).", 0, False, INK),
        ("141 detected on both, 10 only on Illumina, none only on PacBio, none missed by every run.", 0, False, INK),
        ("The ten Illumina-only viruses were not further characterised — a follow-up target.", 0, False, INK),
    ], size=15, gap=12)
    notes(s, "Figure 7.")

    # 23 ----------------------------------------------------------- limitations
    s = frame("Limitations", "What these results do — and do not — support")
    bullets(s, Inches(0.7), Inches(1.75), Inches(11.9), Inches(4.9), [
        ("GOTTCHA2 ran a locally patched, unvalidated binary on PacBio: its PacBio lead is not yet robust cross-modality evidence.", 0, False, INK),
        ("Simulated data: no host nucleic acid, inhibitors, chimeras or quasispecies; RS II error model predates current HiFi chemistry.", 0, False, INK),
        ("Every virus guaranteed ≥ 1 read; PacBio delivers about half the target-virus bases of Illumina — recall is conditional on presence.", 0, False, INK),
        ("Precision compared across different reference databases; database coverage of the 151 targets and 13 sentinels not assessed.", 0, False, INK),
        ("Three replicates, seven pipelines; nf-core pipelines share a classifier and database — rank correlations are directional.", 0, False, INK),
        ("sourmash's collapse attributed to assembly by inference; it was not run on PacBio reads.", 0, False, INK),
    ], size=16, gap=12)
    notes(s, "Also: no true negatives, so no specificity/MCC; genus resolution incomplete (124 of 151).")

    # 24 ----------------------------------------------------------- recommendations
    s = frame("Recommendations", "For a programme choosing or migrating a pipeline")
    bullets(s, Inches(0.7), Inches(1.75), Inches(11.9), Inches(4.9), [
        ("Specify and validate the input stage — reads vs assembled contigs — not only the classifier and database.", 0, True, ACCENT),
        ("Prefer read-level over assembly-first classification, especially with long-read sequencing at low viral depth.", 0, True, ACCENT),
        ("Verify detection of the specific sentinel pathogens you track, independently of aggregate F1.", 0, True, ACCENT),
        ("Run viral positive controls in production — both failures were silent.", 0, True, ACCENT),
        ("Re-measure any read-count threshold on your own library and pipeline.", 0, True, ACCENT),
        ("Re-benchmark when migrating a validated pipeline to a new modality, even if it transfers well.", 0, True, ACCENT),
    ], size=17, gap=14)
    notes(s, "Six recommendations from the Discussion.")

    # 25 ----------------------------------------------------------- closing
    n += 1
    s = prs.slides.add_slide(prs.slide_layouts[6])
    rect(s, 0, 0, W, H, ACCENT)
    textbox(s, Inches(0.85), Inches(1.9), Inches(11.6), Inches(2.0),
            "Viral detection is portable across sequencing modalities — provided the pipeline is given reads, not contigs.",
            size=30, bold=True, color=WHITE, space_after=6, line=1.25)
    textbox(s, Inches(0.85), Inches(4.2), Inches(11.6), Inches(1.4),
            "The input stage, not the classifier or the platform, should be specified and re-validated.\n"
            "Aggregate F1 cannot replace verifying specific sentinel targets.",
            size=17, color=HILITE, space_after=6, line=1.25)
    textbox(s, Inches(0.85), Inches(6.5), Inches(11.6), Inches(0.5),
            "Questions?   Full manuscript: Seven-pipeline-comparison-combined_V5.9.docx (EN + CH); supplementary tables S1–S6.",
            size=13, color=BAND, space_after=0)

    out = ROOT / "Seven-pipeline-benchmark-group-meeting-V5.9.pptx"
    prs.save(out)
    print("wrote", out, "slides:", len(prs.slides))


if __name__ == "__main__":
    build()
