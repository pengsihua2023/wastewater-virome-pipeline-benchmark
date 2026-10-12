#!/usr/bin/env python3
"""Recall of floor-affected viruses vs the rest, and assembly of the coverage supplementary workbook.

Floor-affected = viruses whose proportional expected PacBio read count (designed share x total reads) is < 1
(mean of three samples), i.e. whose PacBio target was raised to one read.
Detection status is taken from interim/viral_detection_detail.tsv (numeric-taxid pipelines only; sourmash is
name-matched and excluded here).
"""
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
COV = ROOT / "analysis_708" / "interim" / "coverage"
INTDIR = ROOT / "analysis_708" / "interim"
d = pd.read_csv(ROOT / "analysis_708" / "interim" / "viral_detection_detail.tsv", sep="\t")
v = pd.read_csv(COV / "per_virus_coverage.tsv", sep="\t")
d = d[d.status.isin(["TP", "FN"])].copy()
d["key"] = pd.to_numeric(d["key"], errors="coerce")
d = d.dropna(subset=["key"])
d["key"] = d["key"].astype(int)
floor = set(v.loc[v.pacbio_floor_applied, "taxid"])
d["group"] = d["key"].map(lambda k: "floor-affected (n=53)" if k in floor else "other (n=98)")
per = (d.groupby(["platform", "pipeline", "group"]).status.apply(lambda s: (s == "TP").mean())
       .rename("recall_at_ge1_read").reset_index())
read_fed = ["TaxProfiler", "TaxTriage", "mag-nf", "GOTTCHA2", "CLARK"]
pooled = (d[d.pipeline.isin(read_fed)].groupby(["platform", "group"]).status
          .apply(lambda s: (s == "TP").mean()).rename("recall_at_ge1_read").reset_index())
pooled["pipeline"] = "five read-fed pipelines (pooled)"
res = pd.concat([per, pooled], ignore_index=True)
res.to_csv(COV / "floor_vs_other_recall.tsv", sep="\t", index=False)

out = ROOT / "Supplementary_Tables_coverage.xlsx"
with pd.ExcelWriter(out) as xw:
    readme = [
        ("SHEETS", "", ""),
        ("S1_class_composition", "Reads and bases per class (151 target viruses, 67 bacteriophages, 490 bacteria/archaea) in each dataset", ""),
        ("S2_per_virus_coverage", "Designed share, expected and observed reads and mean depth of each of the 151 viruses", ""),
        ("S3_sentinel_coverage", "The 13 sentinel viruses (subset of S2)", ""),
        ("S4_pacbio_read_qc", "Length and predicted accuracy of the simulated PacBio HiFi reads", ""),
        ("S5_minimum_read_rule", "Number of PacBio targets raised to a single read", ""),
        ("S6_floor_vs_other_recall", "Detection (>=1 read) of the 53 viruses with a proportional expectation below one PacBio read versus the other 98", ""),
        ("S7_metrics_by_threshold", "Recall, precision, F1, FP and TP at >=1, >=10 and >=100 clade reads (mean and SD of three samples)", ""),
        ("S8_per_virus_detection", "For each of the 151 viruses, the number of the three samples (0-3) in which it is detected, by modality, pipeline and threshold", ""),
        ("S9_full_community_708", "All 708 simulated community members with taxid, rank, genome size, source-sample reads, design reads, designed abundance and simulated reads", ""),
        ("S10_bootstrap_CI", "95% bootstrap confidence intervals over the 151 viruses for recall and F1 (>=1 and >=10 reads), and paired PacBio - Illumina differences; false positives held at their observed means", ""),
        ("S11_floor_excluded", "Metrics recomputed on the 98 viruses not affected by the one-read minimum (false positives unchanged)", ""),
        ("S12_genus_level", "Species-level versus genus-level recall, precision, F1 and false positives (lenient genus matching, >=1 read)", ""),
        ("S13_false_positive_profile", "Read support and genus relatedness of strict false positives; read-support columns are blank for GOTTCHA2 and sourmash, whose count fields are not read counts", ""),
        ("S14_kraken2_misses", "Viruses missed by nf-core/taxprofiler on Illumina: simulated read pairs, whether the taxid exists in the Kraken2 report, clade reads, and detection by the other pipelines", ""),
        ("", "", ""),
        ("COLUMN DEFINITIONS (S9)", "", ""),
        ("source_reads_direct", "Reads of the real wastewater sample HTP-2024-01-15 assigned by Kraken2 directly to this taxonomic node (not to its descendants)", ""),
        ("source_reads_clade", "Reads of HTP-2024-01-15 assigned to this node or any descendant node (cumulative, so clade >= direct)", ""),
        ("design_reads", "Number of reads designed for this member in the community design, in units of the source sample's Kraken2 reads; equals the source-sample reads unless design_basis says it was assigned manually", ""),
        ("design_basis", "Whether design_reads equals the source-sample clade or direct reads, or was assigned manually (virus below the 20-read filter, or not detected in the source sample)", ""),
        ("designed_share_pct", "Designed share of the whole community = domain fraction x design_reads / sum of design_reads within the domain; domain fractions are fixed at 4% (non-phage viruses), 8% (bacteriophages) and 88% (bacteria/archaea)", ""),
        ("illumina_pairs_s0-2, pacbio_reads_s0-2", "Read pairs (Illumina) and reads (PacBio HiFi) of this member in each of the three simulated datasets, counted from the FASTQ read headers", ""),
        ("", "", ""),
        ("THRESHOLDS", "", ""),
        (">=1, >=10, >=100", "A virus counts as detected only if the pipeline reports at least this many clade reads for it; >=1 is the primary threshold used in the main tables and figures", ""),
    ]
    pd.DataFrame(readme, columns=["item", "description", ""]).iloc[:, :2].to_excel(xw, sheet_name="README", index=False)
    pd.read_csv(COV / "class_summary.tsv", sep="\t").to_excel(xw, sheet_name="S1_class_composition", index=False)
    v.to_excel(xw, sheet_name="S2_per_virus_coverage", index=False)
    pd.read_csv(COV / "sentinel_coverage.tsv", sep="\t").to_excel(xw, sheet_name="S3_sentinel_coverage", index=False)
    pd.read_csv(COV / "pacbio_read_qc.tsv", sep="\t").to_excel(xw, sheet_name="S4_pacbio_read_qc", index=False)
    pd.read_csv(COV / "floor_summary.tsv", sep="\t").to_excel(xw, sheet_name="S5_minimum_read_rule", index=False)
    res.to_excel(xw, sheet_name="S6_floor_vs_other_recall", index=False)
    # S7: metrics at the three read-count thresholds (151 non-phage viruses; mean and SD over three samples)
    vm = pd.read_csv(ROOT / "analysis_708" / "interim" / "viral_sample_metrics.tsv", sep="\t")
    vm = vm[vm.ground_truth == "non-phage virus (n=151)"]
    cols = ["recall", "precision", "f1", "FP", "TP"]
    agg = vm.groupby(["platform", "pipeline", "min_reads"])[cols].agg(["mean", "std"])
    agg.columns = [f"{c}_{s}" for c, s in agg.columns]
    agg = agg.reset_index().rename(columns={"min_reads": "threshold_min_clade_reads"})
    pooled = []
    read_fed = ["TaxProfiler", "TaxTriage", "mag-nf", "GOTTCHA2", "CLARK"]
    means = vm.groupby(["platform", "pipeline", "min_reads"])[cols].mean().reset_index()
    for label, sub in (("mean of all seven pipelines", means), ("mean of five read-fed pipelines", means[means.pipeline.isin(read_fed)])):
        t = sub.groupby(["platform", "min_reads"])[cols].mean().reset_index()
        t.insert(1, "pipeline", label)
        pooled.append(t.rename(columns={"min_reads": "threshold_min_clade_reads", **{c: f"{c}_mean" for c in cols}}))
    pd.concat([agg] + pooled, ignore_index=True).to_excel(xw, sheet_name="S7_metrics_by_threshold", index=False)
    # S8: per-virus detection at the three thresholds (from 25_per_virus_threshold_detection.py)
    pd.read_csv(COV / "per_virus_detection_by_threshold.tsv", sep="\t").to_excel(xw, sheet_name="S8_per_virus_detection", index=False)
    # S9: complete simulated community of 708 members (from 26_full_community_table.py)
    pd.read_csv(COV / "full_community_708.tsv", sep="\t").to_excel(xw, sheet_name="S9_full_community_708", index=False)
    # S10-S14: analyses added in response to the pre-submission review (27_revision_analyses.py)
    pd.read_csv(COV / "bootstrap_ci.tsv", sep="\t").to_excel(xw, sheet_name="S10_bootstrap_CI", index=False)
    pd.read_csv(COV / "floor_excluded_metrics.tsv", sep="\t").to_excel(xw, sheet_name="S11_floor_excluded", index=False)
    pd.read_csv(INTDIR / "species_vs_genus_comparison.tsv", sep="\t").to_excel(xw, sheet_name="S12_genus_level", index=False)
    pd.read_csv(COV / "false_positive_profile.tsv", sep="\t").to_excel(xw, sheet_name="S13_false_positive_profile", index=False)
    pd.read_csv(COV / "kraken2_misses_illumina.tsv", sep="\t").to_excel(xw, sheet_name="S14_kraken2_misses", index=False)

    # light formatting: frozen header, filter buttons, and the 13 sentinel viruses highlighted
    from xlsxwriter.utility import xl_col_to_name
    hl = xw.book.add_format({"bg_color": "#FDE9D9"})
    xw.sheets["README"].set_column(0, 0, 38)
    xw.sheets["README"].set_column(1, 1, 140)
    sheet_frames = {
        "S2_per_virus_coverage": pd.read_csv(COV / "per_virus_coverage.tsv", sep="\t"),
        "S3_sentinel_coverage": pd.read_csv(COV / "sentinel_coverage.tsv", sep="\t"),
        "S8_per_virus_detection": pd.read_csv(COV / "per_virus_detection_by_threshold.tsv", sep="\t"),
        "S9_full_community_708": pd.read_csv(COV / "full_community_708.tsv", sep="\t"),
    }
    for name, df in sheet_frames.items():
        ws = xw.sheets[name]
        col = "sentinel_virus" if "sentinel_virus" in df.columns else "sentinel"
        letter = xl_col_to_name(df.columns.get_loc(col))
        ws.conditional_format(1, 0, len(df), len(df.columns) - 1,
                              {"type": "formula", "criteria": f"=${letter}2=TRUE", "format": hl})
    for name, ws in xw.sheets.items():
        if name != "README":
            ws.freeze_panes(1, 0)
print("wrote", out)
