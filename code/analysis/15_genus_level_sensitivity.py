#!/usr/bin/env python3
"""Genus-level lenient matching as a sensitivity analysis to the strict
species-level taxid scoring.

Motivation: spot-checking the strict false positives showed that many are
taxonomic neighbours of ground-truth members rather than unrelated taxa -
"Mamastrovirus hominis" and "Human astrovirus" against ground-truth
"Mamastrovirus 1", "Norovirus GI" against "Norwalk virus",
"Alphainfluenzavirus influenzae" against the influenza A serotype nodes. Viral
nomenclature is being revised on a large scale (ICTV binomial names), so strict
taxid equality systematically overstates the false-positive count.

The ground-truth document carries no genus column, so the taxid -> genus map is
built from the Kraken2 report with the widest coverage (the nf-core/mag Illumina
report, which was produced with --report-zero-counts and therefore enumerates
the whole database), by indentation-tracked descent within the Viruses branch.

A genus-level TRUE POSITIVE is a ground-truth virus for which the pipeline
reported any species-or-below node in the same genus. A genus-level FALSE
POSITIVE is a reported viral genus containing no ground-truth member.

Main results in the manuscripts remain the strict species-level ones; this
script quantifies how much of the strict false-positive burden is nomenclatural.

Writes 12_genus_level_sensitivity.xlsx.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
AN = ROOT / "analysis_708"
INT = AN / "interim"

MAP_SOURCE = ROOT / "Illumina/mag-nf/sample_0/Sample_708_illumina_0.kraken2_report.txt"
SPECIES_RANKS = {"S", "S1", "S2", "S3", "S4", "S5", "species", "strain", "serotype",
                 "no rank", "subspecies", "biotype", "isolate", "genotype"}
GENUS_RANKS = {"G", "G1", "G2"}
FAMILY_RANKS = {"F", "F1"}


def build_maps(path):
    """taxid -> (genus_taxid, genus_name, family_taxid, family_name) inside Viruses."""
    tax2genus, tax2family, name2genus = {}, {}, {}
    viral_indent = None
    cur_g = cur_gn = cur_f = cur_fn = None
    g_indent = f_indent = None
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        rank, taxid, name = p[3].strip(), p[4].strip(), p[5]
        indent = len(name) - len(name.lstrip())
        nm = name.strip()
        if taxid == "10239":
            viral_indent = indent
            continue
        if viral_indent is None or indent <= viral_indent:
            if viral_indent is not None and indent <= viral_indent:
                viral_indent = None
            continue
        # inside the Viruses branch
        if rank in FAMILY_RANKS:
            cur_f, cur_fn, f_indent = taxid, nm, indent
            cur_g = cur_gn = None
            g_indent = None
        elif rank in GENUS_RANKS:
            if g_indent is None or indent <= g_indent:
                cur_g, cur_gn, g_indent = taxid, nm, indent
        elif rank in SPECIES_RANKS:
            if f_indent is not None and indent <= f_indent:
                cur_f = cur_fn = None
            if g_indent is not None and indent <= g_indent:
                cur_g = cur_gn = None
            tax2genus[taxid] = (cur_g, cur_gn)
            tax2family[taxid] = (cur_f, cur_fn)
            if cur_gn:
                name2genus.setdefault(nm.lower(), cur_gn)
    return tax2genus, tax2family, name2genus


def main():
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    nonphage = gt[gt.category.str.startswith("Virus (")].copy()
    phage = gt[gt.category == "Bacteriophage"].copy()
    assert len(nonphage) == 151

    tax2genus, tax2family, name2genus = build_maps(MAP_SOURCE)
    print(f"genus map built from {MAP_SOURCE.name}: "
          f"{len(tax2genus)} viral species-or-below nodes, "
          f"{len({g for g, _ in tax2genus.values() if g})} distinct genera")

    nonphage["gid"] = nonphage["taxid"].astype(str).map(
        lambda t: (tax2genus.get(t) or (None, None))[0])
    nonphage["gname"] = nonphage["taxid"].astype(str).map(
        lambda t: (tax2genus.get(t) or (None, None))[1])
    resolved = nonphage.gid.notna().sum()
    print(f"ground-truth viruses with a resolvable genus: {resolved}/151")
    unresolved = nonphage[nonphage.gid.isna()][["std_name", "taxid"]]

    def gkey(x):
        return str(x).strip().lower() if x and str(x) != "nan" else None

    nonphage["gkey"] = nonphage["gname"].map(gkey)
    gt_genera = set(nonphage.gkey.dropna())
    phage_tax = set(phage.taxid.astype(str))

    long = pd.read_csv(INT / "viral_detections_long.tsv", sep="\t", low_memory=False)
    long = long[long.is_viral & long["rank"].isin(SPECIES_RANKS) & (long.reads > 0)]
    long["taxid_s"] = long["taxid"].astype(str).str.replace(r"\.0$", "", regex=True)
    long["gkey"] = long["taxid_s"].map(
        lambda t: gkey((tax2genus.get(t) or (None, None))[1]))
    # sourmash carries no taxid but its lineage names the genus directly
    miss = long.gkey.isna()
    long.loc[miss, "gkey"] = long.loc[miss, "genus_from_lineage"].map(gkey)

    rows = []
    for (plat, pipe, samp), g in long.groupby(["platform", "pipeline", "sample"]):
        g = g[~g["taxid_s"].isin(phage_tax)]          # phage out of scope
        obs_gen = set(g.gkey.dropna())
        obs_gen = {x for x in obs_gen if x and not x.startswith("unclassified")}
        tp = obs_gen & gt_genera
        fp = obs_gen - gt_genera
        fn = gt_genera - obs_gen
        TP, FP, FN = len(tp), len(fp), len(fn)
        rec = TP / (TP + FN) if TP + FN else np.nan
        pre = TP / (TP + FP) if TP + FP else np.nan
        rows.append(dict(platform=plat, pipeline=pipe, sample=samp,
                         n_truth_genera=len(gt_genera),
                         TP=TP, FP=FP, FN=FN, recall=rec, precision=pre,
                         f1=(2 * pre * rec / (pre + rec)) if (pre and rec) else 0.0,
                         n_nodes_without_genus=int(g.gkey.isna().sum())))
    gen = pd.DataFrame(rows)
    agg = (gen.groupby(["platform", "pipeline"])
           .agg(TP=("TP", "mean"), FP=("FP", "mean"), FN=("FN", "mean"),
                recall=("recall", "mean"), precision=("precision", "mean"),
                f1=("f1", "mean")).reset_index())

    # ---- compare with the strict species-level result
    strict = pd.read_csv(INT / "viral_aggregate_metrics.tsv", sep="\t")
    cmp = strict[["platform", "pipeline", "recall_mean", "precision_mean",
                  "f1_mean", "FP_mean"]].merge(
        agg, on=["platform", "pipeline"], suffixes=("_strict_species", "_genus"))
    cmp = cmp.rename(columns={"recall_mean": "recall_species",
                              "precision_mean": "precision_species",
                              "f1_mean": "f1_species", "FP_mean": "FP_species",
                              "recall": "recall_genus", "precision": "precision_genus",
                              "f1": "f1_genus", "FP": "FP_genus"})
    cmp["FP_reduction_x"] = cmp["FP_species"] / cmp["FP_genus"].replace(0, np.nan)
    cmp["precision_gain"] = cmp["precision_genus"] - cmp["precision_species"]

    with pd.ExcelWriter(AN / "12_genus_level_sensitivity.xlsx") as xw:
        agg.to_excel(xw, sheet_name="genus_aggregate", index=False)
        gen.to_excel(xw, sheet_name="genus_per_sample", index=False)
        cmp.to_excel(xw, sheet_name="species_vs_genus", index=False)
        nonphage[["std_name", "taxid", "gid", "gname", "gkey"]].to_excel(
            xw, sheet_name="ground_truth_genus_map", index=False)
        unresolved.to_excel(xw, sheet_name="genus_unresolved", index=False)

    agg.to_csv(INT / "genus_aggregate_metrics.tsv", sep="\t", index=False)
    cmp.to_csv(INT / "species_vs_genus_comparison.tsv", sep="\t", index=False)

    print(f"\nground-truth genera: {len(gt_genera)}")
    print("\n=== GENUS-level (lenient) vs SPECIES-level (strict), >=1 read")
    for plat in ["Illumina", "PacBio"]:
        print(f"\n--- {plat}")
        s = cmp[cmp.platform == plat].sort_values("f1_genus", ascending=False)
        for _, r in s.iterrows():
            print(f"  {r.pipeline:20s} R {r.recall_species:.3f}->{r.recall_genus:.3f}  "
                  f"P {r.precision_species:.3f}->{r.precision_genus:.3f}  "
                  f"F1 {r.f1_species:.3f}->{r.f1_genus:.3f}  "
                  f"FP {r.FP_species:7.1f}->{r.FP_genus:6.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
