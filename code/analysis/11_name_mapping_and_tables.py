#!/usr/bin/env python3
"""Name-standardisation / manual-review workbook and the manuscript tables."""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
AN = ROOT / "analysis_708"
INT = AN / "interim"
TAB = AN / "_superseded" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

ORDER = ["TaxProfiler", "TaxTriage", "mag-nf", "KrakenMetaReads-nf",
         "GOTTCHA2", "CLARK", "sourmash"]


def norm_key(n):
    s = re.sub(r"^[a-z]__", "", str(n).strip())
    s = re.sub(r"[^a-z0-9]+", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


def main():
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    det = pd.read_csv(INT / "detections_long.tsv", sep="\t",
                      dtype={"taxid": str}, low_memory=False)
    det["name_key"] = det["name"].map(norm_key)
    det["taxid_num"] = pd.to_numeric(det["taxid"], errors="coerce").astype("Int64")

    gt_keys = set(gt["match_key"])
    gt_tax = set(gt["taxid_num"].dropna().astype(int)) if "taxid_num" in gt else set(
        pd.to_numeric(gt["taxid"], errors="coerce").dropna().astype(int))

    # 1. how each ground-truth member is written by each taxid-carrying pipeline
    obs_names = (det[det.taxid_num.notna()]
                 .assign(t=lambda d: d.taxid_num.astype(int))
                 .query("t in @gt_tax")
                 .groupby("t")["name"].agg(lambda s: sorted(set(s))))
    gtx = gt.copy()
    gtx["taxid_int"] = pd.to_numeric(gtx["taxid"], errors="coerce").astype("Int64")
    gtx["names_reported_by_pipelines"] = gtx["taxid_int"].map(
        lambda t: "; ".join(obs_names.get(int(t), [])) if pd.notna(t) else "")
    gtx["name_differs_from_document"] = [
        "yes" if v and norm_key(r) not in {norm_key(x) for x in v.split("; ")} else
        ("no" if v else "not reported by any taxid-based pipeline")
        for v, r in zip(gtx["names_reported_by_pipelines"], gtx["std_name"])]

    # 2. sourmash: which ground-truth members can be matched by name at all
    sm = det[det.parser.isin(["sourmash_prok", "sourmash_viral"])]
    sm_keys = set(sm["name_key"])
    sm_tab = gt[["category", "raw_name", "std_name", "match_key", "taxid"]].copy()
    sm_tab["matched_by_sourmash_name"] = sm_tab["match_key"].isin(sm_keys)
    sm_tab["note"] = np.where(
        sm_tab["matched_by_sourmash_name"], "",
        "sourmash reports GTDB (prokaryote) or NCBI-viral lineage strings; a "
        "missing match may be a genuine non-detection OR a naming-system mismatch "
        "- requires manual confirmation")

    # 3. explicit manual-review queue
    review = []
    for _, r in gt.iterrows():
        reasons = []
        if r["rank_class"] == "strain (S1)":
            reasons.append("ground-truth node is a sub-species (S1) strain node")
        if r["rank_class"] in ("no rank", "serotype", "unspecified"):
            reasons.append(f"ground-truth node rank is '{r['ncbi_rank']}' - not a species node")
        if r["category"].startswith("Virus (manual"):
            reasons.append("manually added to the community after automated selection")
        if reasons:
            review.append(dict(raw_name=r["raw_name"], std_name=r["std_name"],
                               taxid=r["taxid"], ncbi_rank=r["ncbi_rank"],
                               category=r["category"],
                               reason_for_review="; ".join(reasons),
                               resolution="scored by exact NCBI taxid at any "
                                          "species-or-below node; no name-based merging applied"))
    review = pd.DataFrame(review)

    with pd.ExcelWriter(AN / "04_name_mapping_and_manual_review.xlsx") as xw:
        gtx[["category", "raw_name", "std_name", "match_key", "taxid", "ncbi_rank",
             "names_reported_by_pipelines", "name_differs_from_document"]].to_excel(
            xw, sheet_name="taxid_to_reported_names", index=False)
        sm_tab.to_excel(xw, sheet_name="sourmash_name_matching", index=False)
        review.to_excel(xw, sheet_name="manual_review_queue", index=False)
        pd.DataFrame([
            dict(rule="primary matching", detail="exact NCBI taxid equality"),
            dict(rule="rank accepted", detail="S, S1-S5, species, strain, serotype, "
                                              "no rank, subspecies, biotype, isolate, genotype"),
            dict(rule="detection threshold", detail="clade reads > 0 (main analysis); "
                                                   ">=10 and >=100 reported as sensitivity"),
            dict(rule="case/punctuation", detail="lower-cased, non-alphanumerics collapsed "
                                                 "to single spaces (name matching only)"),
            dict(rule="GTDB prefixes", detail="leading d__/p__/c__/o__/f__/g__/s__ stripped"),
            dict(rule="synonym merging", detail="NOT performed - no curated synonym source "
                                                "was available in the project files"),
            dict(rule="genus/family scoring", detail="NOT performed - the ground-truth "
                                                    "document carries no genus or family column"),
        ]).to_excel(xw, sheet_name="matching_rules", index=False)

    # ------------------------------------------------------------------ tables
    agg = pd.read_csv(INT / "aggregate_metrics.tsv", sep="\t")
    met = pd.read_csv(INT / "sample_level_metrics.tsv", sep="\t")
    prim = met[met.min_reads == 1]
    comp = pd.read_csv(INT / "realised_composition.tsv", sep="\t")
    dom = pd.read_csv(INT / "domain_recall.tsv", sep="\t")
    plat = pd.read_csv(INT / "platform_contrast.tsv", sep="\t")
    miss = pd.read_csv(INT / "missed_summary.tsv", sep="\t")

    # Table 1: dataset overview
    t1 = comp[["platform", "sample", "total_reads", "pct_unclassified", "pct_bacteria",
               "pct_archaea", "pct_virus"]].copy()
    t1["designed_pct_bacteria_archaea"] = 88.0
    t1["designed_pct_viral_incl_phage"] = 12.0
    t1.to_csv(TAB / "Table1_dataset_overview.csv", index=False)

    # Table 2: pipelines, versions, databases  (hand-curated from run scripts)
    t2 = pd.DataFrame([
        dict(pipeline="TaxProfiler", implementation="nf-core/taxprofiler v1.2.0",
             classifier="Kraken2 v2.1.3 (+ Bracken, Illumina only)",
             database="k2_standard_20240605 (Kraken2 Standard, 2024-06-05)",
             key_parameters="Bracken -r 150 (Illumina); defaults otherwise",
             platforms="Illumina, PacBio",
             evidence_used="sample-level Kraken2 report"),
        dict(pipeline="TaxTriage", implementation="nf-core/taxtriage v1.0 (local checkout)",
             classifier="Kraken2 v2.1.3",
             database="k2_standard_20240605",
             key_parameters="pipeline defaults", platforms="Illumina, PacBio",
             evidence_used="sample-level Kraken2 report"),
        dict(pipeline="mag-nf", implementation="nf-core/mag v3.1.0",
             classifier="Kraken2 v2.1.3",
             database="k2_standard_20240605",
             key_parameters="--skip_binning --skip_binqc --skip_gtdbtk --skip_prokka "
                            "--skip_quast", platforms="Illumina, PacBio",
             evidence_used="Illumina: read-level Kraken2 report; PacBio: contig-level "
                           "Kraken2 report (read-level file unusable, see audit)"),
        dict(pipeline="KrakenMetaReads-nf", implementation="in-house Nextflow DSL2 workflow",
             classifier="Kraken2 v2.1.3 applied to assembled contigs",
             database="k2_standard_20240605",
             key_parameters="Illumina: metaSPAdes 3.15.5 + MEGAHIT 1.2.9, merged; "
                            "PacBio: metaFlye 2.9.2 --meta; fastp 0.23.4; minimap2 2.26; "
                            "samtools 1.17", platforms="Illumina, PacBio",
             evidence_used="Illumina: merged two-assembler report; PacBio: metaFlye report"),
        dict(pipeline="GOTTCHA2", implementation="gottcha2 (conda env gottcha2_env)",
             classifier="unique signature-fragment profiling",
             database="GOTTCHA2 species-level signature DB (gottcha_db.species.fna)",
             key_parameters="-xm sr (Illumina); -xm map-hifi (PacBio, applied via a local "
                            "patch to gottcha2.py); -t 32", platforms="Illumina, PacBio",
             evidence_used="species-level profile TSV"),
        dict(pipeline="CLARK", implementation="CLARK (local build)",
             classifier="discriminative k-mer classification",
             database="CLARK_DB, targets set at species level",
             key_parameters="k = 20, mode 0 (full); gamma >= 0.2 post-filter reported as "
                            "a sensitivity analysis", platforms="Illumina, PacBio",
             evidence_used="estimate_abundance output (unfiltered)"),
        dict(pipeline="sourmash", implementation="sourmash 4.9.4 (Apptainer container)",
             classifier="FracMinHash gather + tax metagenome",
             database="GTDB reps RS226 k31 (prokaryote) and NCBI viruses 2025.01 k31",
             key_parameters="k = 31; scaled = 50 (viral), scaled = 1000 (GTDB)",
             platforms="Illumina, PacBio",
             evidence_used="species rows of the two summarized.csv files, combined"),
    ])
    t2.to_csv(TAB / "Table2_pipelines_versions_databases.csv", index=False)

    # Table 3/4: main performance
    t3 = (agg[["platform", "pipeline", "TP_mean", "TP_sd", "FP_mean", "FP_sd",
               "FN_mean", "FN_sd", "recall_mean", "recall_sd", "recall_ci_lo",
               "recall_ci_hi", "precision_mean", "precision_sd", "f1_mean", "f1_sd"]]
          .copy())
    t3["pipeline"] = pd.Categorical(t3.pipeline, ORDER, ordered=True)
    t3 = t3.sort_values(["platform", "f1_mean"], ascending=[True, False])
    t3.to_csv(TAB / "Table3_main_performance.csv", index=False)

    plat.to_csv(TAB / "Table4_platform_contrast.csv", index=False)
    dom.to_csv(TAB / "Table5_domain_recall.csv", index=False)
    miss.head(100).to_csv(TAB / "Table6_most_frequently_missed.csv", index=False)
    prim.to_csv(TAB / "TableS1_sample_level_metrics.csv", index=False)
    met.to_csv(TAB / "TableS2_read_threshold_sweep.csv", index=False)

    print("wrote name-mapping workbook and", len(list(TAB.glob("*.csv"))), "tables")
    return 0


if __name__ == "__main__":
    sys.exit(main())
