#!/usr/bin/env python3
"""Two controlled comparisons that isolate causes of viral detection loss.

(1) Assembly discards viral signal, on BOTH platforms.
    KrakenMetaReads-nf is the only Illumina pipeline scored on assembled contigs.
    Comparing the Kraken2 report over its contigs with the Kraken2 report over
    the same sample's reads (nf-core/taxprofiler, identical database) measures
    how much viral signal assembly costs, independently of any classifier
    difference. The same comparison is made on PacBio.

(2) A single QC step doubles the viral false-positive count.
    nf-core/mag and nf-core/taxprofiler classify the SAME reads (identical total
    counts) with the SAME Kraken2 version against the SAME database, and neither
    sets --confidence. The one difference recorded in their run scripts is that
    nf-core/mag runs fastp first. This script quantifies the consequence.

Writes 15_assembly_loss_and_preprocessing.xlsx.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
AN = ROOT / "analysis_708"
INT = AN / "interim"

SPECIES = ["S", "S1", "S2", "S3"]

REPORTS = {
    # (platform, label, input_type) -> path template with {i}
    ("Illumina", "reads (nf-core/taxprofiler)", "reads"):
        "Illumina/TaxProfiler/sample_{i}/sample_708_{i}_null_Viral_ref.kraken2.kraken2.report.txt",
    ("Illumina", "reads (nf-core/mag)", "reads"):
        "Illumina/mag-nf/sample_{i}/Sample_708_illumina_{i}.kraken2_report.txt",
    ("Illumina", "contigs (metaSPAdes)", "contigs"):
        "Illumina/KrakenMetaReads-nf/sample_{i}/sample_{i}_708_Illumina_spades_report.txt",
    ("Illumina", "contigs (MEGAHIT)", "contigs"):
        "Illumina/KrakenMetaReads-nf/sample_{i}/sample_{i}_708_Illumina_megahit_report.txt",
    ("PacBio", "reads (nf-core/taxprofiler)", "reads"):
        "Pacbio/TaxProfiler/sample_{i}/sample_708_{i}_PacBio_null_Viral_ref.kraken2.kraken2.report.txt",
    ("PacBio", "contigs (metaFlye)", "contigs"):
        "Pacbio/KrakenMetaReads-nf/sample_{i}/pacbio_hifi_reads_sample_{i}_flye_report.txt",
    ("PacBio", "contigs (nf-core/mag)", "contigs"):
        "Pacbio/mag-nf/sample_{i}/sample_708_pacbio_{i}_contigs_kraken2_report.txt",
}


def parse(path):
    """Kraken2 report -> rows with a viral-branch flag (indentation-tracked)."""
    rows, vi = [], None
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        _pct, clade, _direct, rank, taxid, name = p[:6]
        ind = len(name) - len(name.lstrip())
        if taxid.strip() == "10239":
            vi, viral = ind, True
        elif vi is not None and ind > vi:
            viral = True
        else:
            if vi is not None and ind <= vi:
                vi = None
            viral = False
        try:
            reads = float(clade)
        except ValueError:
            continue
        rows.append(dict(rank=rank.strip(), taxid=taxid.strip(),
                         reads=reads, viral=viral))
    return pd.DataFrame(rows)


def main():
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    vt = set(pd.to_numeric(gt[gt.category.str.startswith("Virus (")]["taxid"],
                           errors="coerce").dropna().astype(int))
    ph = set(pd.to_numeric(gt[gt.category == "Bacteriophage"]["taxid"],
                           errors="coerce").dropna().astype(int))

    # ------------------------------------------------ (1) assembly signal loss
    rows = []
    for (plat, label, itype), tmpl in REPORTS.items():
        for i in range(3):
            f = ROOT / tmpl.format(i=i)
            if not f.exists():
                continue
            d = parse(f)
            total = d[d["rank"] == "U"].reads.sum() + d[d["rank"] == "R"].reads.sum()
            vnode = d[d.viral & (d.reads > 0)]
            vroot = d[(d.taxid == "10239")].reads.max() if (d.taxid == "10239").any() else 0
            sp = vnode[vnode["rank"].isin(SPECIES)].copy()
            sp["t"] = pd.to_numeric(sp.taxid, errors="coerce")
            obs = set(sp.t.dropna().astype(int))
            rows.append(dict(platform=plat, evidence=label, input_type=itype,
                             sample=f"sample_{i}", units_total=total,
                             viral_units=vroot,
                             viral_species_nodes=len(sp),
                             gt_viruses_hit=len(obs & vt),
                             viral_fp=len(obs - vt - ph)))
    asm = pd.DataFrame(rows)
    asm_agg = (asm.groupby(["platform", "evidence", "input_type"])
               [["units_total", "viral_units", "viral_species_nodes",
                 "gt_viruses_hit", "viral_fp"]].mean().reset_index())
    for plat in ["Illumina", "PacBio"]:
        ref = asm_agg[(asm_agg.platform == plat) &
                      (asm_agg.evidence == "reads (nf-core/taxprofiler)")]
        if len(ref):
            base = ref.viral_units.iloc[0]
            asm_agg.loc[asm_agg.platform == plat, "viral_units_vs_reads"] = (
                asm_agg.loc[asm_agg.platform == plat, "viral_units"] / base)

    # ------------------------------- (2) fastp effect: mag-nf vs taxprofiler
    pre = []
    for i in range(3):
        a = parse(ROOT / f"Illumina/TaxProfiler/sample_{i}/sample_708_{i}_null_Viral_ref.kraken2.kraken2.report.txt")
        b = parse(ROOT / f"Illumina/mag-nf/sample_{i}/Sample_708_illumina_{i}.kraken2_report.txt")
        rec = {"sample": f"sample_{i}"}
        for lab, d in [("taxprofiler", a), ("mag-nf", b)]:
            tot = d[d["rank"] == "U"].reads.sum() + d[d["rank"] == "R"].reads.sum()
            unc = d[d["rank"] == "U"].reads.sum()
            sp = d[d.viral & (d.reads > 0) & d["rank"].isin(SPECIES)].copy()
            sp["t"] = pd.to_numeric(sp.taxid, errors="coerce")
            obs = set(sp.t.dropna().astype(int))
            rec[f"{lab}_total"] = tot
            rec[f"{lab}_unclassified_pct"] = 100 * unc / tot
            rec[f"{lab}_viral_reads"] = d[d.taxid == "10239"].reads.max()
            rec[f"{lab}_TP"] = len(obs & vt)
            rec[f"{lab}_FP"] = len(obs - vt - ph)
            rec[f"{lab}_fp_ids"] = obs - vt - ph
        exclusive = rec["mag-nf_fp_ids"] - rec["taxprofiler_fp_ids"]
        rec["mag_only_FP"] = len(exclusive)
        rec["taxprofiler_only_FP"] = len(rec["taxprofiler_fp_ids"] - rec["mag-nf_fp_ids"])
        spb = b[b.viral & (b.reads > 0) & b["rank"].isin(SPECIES)].copy()
        spb["t"] = pd.to_numeric(spb.taxid, errors="coerce")
        ex = spb[spb.t.isin(exclusive)]
        rec["mag_only_FP_median_reads"] = ex.reads.median()
        rec["mag_only_FP_frac_single_read"] = float((ex.reads == 1).mean())
        # do the shared true positives get more reads under fastp?
        spa = a[a.viral & (a.reads > 0) & a["rank"].isin(SPECIES)].copy()
        spa["t"] = pd.to_numeric(spa.taxid, errors="coerce")
        ta = spa[spa.t.isin(vt)].groupby("t").reads.max()
        tb = spb[spb.t.isin(vt)].groupby("t").reads.max()
        com = ta.index.intersection(tb.index)
        rec["shared_TP"] = len(com)
        rec["TP_read_ratio_median"] = float((tb[com] / ta[com]).median())
        for k in list(rec):
            if k.endswith("_fp_ids"):
                del rec[k]
        pre.append(rec)
    prep = pd.DataFrame(pre)

    with pd.ExcelWriter(AN / "15_assembly_loss_and_preprocessing.xlsx") as xw:
        asm_agg.to_excel(xw, sheet_name="assembly_vs_reads", index=False)
        asm.to_excel(xw, sheet_name="assembly_per_sample", index=False)
        prep.to_excel(xw, sheet_name="fastp_effect_mag_vs_taxpro", index=False)

    asm_agg.to_csv(INT / "assembly_vs_reads.tsv", sep="\t", index=False)
    prep.to_csv(INT / "fastp_effect.tsv", sep="\t", index=False)

    print("=== (1) Viral signal retained by assembly, relative to reads")
    print(asm_agg[["platform", "evidence", "input_type", "units_total",
                   "viral_units", "viral_units_vs_reads", "gt_viruses_hit"]]
          .round(3).to_string(index=False))
    print("\n=== (2) nf-core/mag vs nf-core/taxprofiler: same reads, same DB, "
          "fastp the only recorded difference")
    m = prep.mean(numeric_only=True)
    print(f"  total reads classified   taxprofiler {m['taxprofiler_total']:,.0f}   "
          f"mag-nf {m['mag-nf_total']:,.0f}")
    print(f"  unclassified             taxprofiler {m['taxprofiler_unclassified_pct']:.2f}%   "
          f"mag-nf {m['mag-nf_unclassified_pct']:.2f}%")
    print(f"  viral reads              taxprofiler {m['taxprofiler_viral_reads']:,.0f}   "
          f"mag-nf {m['mag-nf_viral_reads']:,.0f}")
    print(f"  viral TP                 taxprofiler {m['taxprofiler_TP']:.0f}   "
          f"mag-nf {m['mag-nf_TP']:.0f}")
    print(f"  viral FP                 taxprofiler {m['taxprofiler_FP']:.0f}   "
          f"mag-nf {m['mag-nf_FP']:.0f}")
    print(f"  FP unique to mag-nf      {m['mag_only_FP']:.0f}   "
          f"unique to taxprofiler {m['taxprofiler_only_FP']:.0f}")
    print(f"  those FP: median reads   {m['mag_only_FP_median_reads']:.0f}   "
          f"single-read fraction {m['mag_only_FP_frac_single_read']:.1%}")
    print(f"  shared TP read ratio     {m['TP_read_ratio_median']:.3f} "
          f"(mag-nf / taxprofiler, {m['shared_TP']:.0f} viruses)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
