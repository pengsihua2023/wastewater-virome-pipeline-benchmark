#!/usr/bin/env python3
"""Why is the realised viral READ fraction ~0.02% when the design targets 12%?

The design's 88 / 8 / 4 split is intended as a share of READS (confirmed by the
study authors, and consistent with the source column being Kraken2 read counts
from a real metatranscriptome). CAMISIM, however, consumes its abundance file as
GENOME COPIES, and reads scale with copies x genome length. Because the viral
members here have genomes 500-1000x smaller than the bacterial ones, a 12% copy
share delivers only ~0.02% of reads.

This script (a) measures the per-domain read yield per unit designed abundance,
(b) shows it matches the genome-length prediction, and (c) computes the copy
abundances that would be required to hit the intended READ shares.

Hypothesis: the design's "Relative abundance in simulation (%)" column behaves
as a GENOME-COPY abundance (what CAMISIM consumes), whereas the measured 0.02%
is a READ fraction. Reads per genome scale with genome length, and the viral
members of this community have genomes 2-3 orders of magnitude smaller than the
bacterial ones, so a 12% copy-number share becomes a ~0.02% read share.

Test: per domain, compute (total observed clade reads) / (total designed
relative abundance) and normalise to bacteria. Under the "abundance == read
fraction" reading this ratio is 1 for every domain. Under the length-weighted
reading it equals (mean genome length of the domain) / (mean bacterial genome
length).

Writes interim/abundance_scale_diagnosis.tsv.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
INT = ROOT / "analysis_708" / "interim"

# Approximate genome sizes for the ORDER-OF-MAGNITUDE prediction only. These are
# textbook ranges for the dominant members of each block, not measured values,
# and are used solely to show that the observed deficit is of the size that
# genome-length weighting predicts.
APPROX_GENOME_BP = {
    "Bacteria/Archaea": 4.0e6,   # typical wastewater bacterial genome
    "Bacteriophage": 4.0e3,      # ssRNA (Leviviricetes) phage
    "Virus": 7.0e3,              # community dominated by small RNA viruses
}

STREAMS = [("Illumina", "taxprofiler_kraken2"), ("PacBio", "taxprofiler_kraken2")]


def main():
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    gt["dom"] = np.where(gt.category == "Bacteria/Archaea", "Bacteria/Archaea",
                np.where(gt.category == "Bacteriophage", "Bacteriophage", "Virus"))
    gt["t"] = pd.to_numeric(gt["taxid"], errors="coerce")
    abund = gt.groupby("dom").relative_abundance_num.sum()

    det = pd.read_csv(INT / "detections_long.tsv", sep="\t",
                      dtype={"taxid": str}, low_memory=False)
    det["t"] = pd.to_numeric(det["taxid"], errors="coerce")

    rows = []
    for plat, parser in STREAMS:
        for samp in ["sample_0", "sample_1", "sample_2"]:
            s = det[(det.platform == plat) & (det.parser == parser) & (det["sample"] == samp)]

            def clade(tid):
                r = s[s.t == tid]
                return float(r.reads.max()) if len(r) else 0.0

            bact = clade(2) + clade(2157)
            vir_clade = clade(10239)
            gtv = gt[gt.dom != "Bacteria/Archaea"]
            sv = s[s.t.isin(gtv.t)].merge(gtv[["t", "dom"]], on="t")
            byd = sv.groupby("dom").reads.sum()

            base = bact / abund["Bacteria/Archaea"]
            for dom in ["Bacteria/Archaea", "Bacteriophage", "Virus"]:
                reads = bact if dom == "Bacteria/Archaea" else float(byd.get(dom, 0.0))
                per_ab = reads / abund[dom]
                rows.append(dict(
                    platform=plat, sample=samp, domain=dom,
                    designed_abundance_pct=abund[dom], observed_reads=reads,
                    reads_per_abundance=per_ab,
                    ratio_vs_bacteria=per_ab / base,
                    predicted_ratio_from_genome_length=
                        APPROX_GENOME_BP[dom] / APPROX_GENOME_BP["Bacteria/Archaea"],
                    total_viral_clade_reads=vir_clade,
                    viral_read_share_pct=100 * vir_clade / (bact + vir_clade)))

    df = pd.DataFrame(rows)
    df.to_csv(INT / "abundance_scale_diagnosis.tsv", sep="\t", index=False)

    summ = (df.groupby(["platform", "domain"])
            [["ratio_vs_bacteria", "predicted_ratio_from_genome_length"]].mean())
    summ["deficit_factor_observed"] = 1 / summ["ratio_vs_bacteria"]
    summ["deficit_factor_predicted"] = 1 / summ["predicted_ratio_from_genome_length"]
    print(summ.round(5).to_string())

    # closed-form prediction of the viral read share from the observed per-domain
    # yields, as a self-consistency check against the measured share
    print()
    for plat in ["Illumina", "PacBio"]:
        d = df[df.platform == plat].groupby("domain").ratio_vs_bacteria.mean()
        num = sum(abund[k] * d[k] for k in ["Bacteriophage", "Virus"])
        den = num + abund["Bacteria/Archaea"] * d["Bacteria/Archaea"]
        obs = df[df.platform == plat].viral_read_share_pct.mean()
        print(f"{plat}: predicted viral read share {100*num/den:.4f}%  |  measured {obs:.4f}%")

    # ---- copy abundances required to achieve the INTENDED read shares --------
    # reads_share_i  proportional to  copies_i * yield_i   ->   copies_i  proportional to  target_i / yield_i
    target = {"Bacteria/Archaea": 88.0, "Bacteriophage": 8.0, "Virus": 4.0}
    corr = []
    for plat in ["Illumina", "PacBio"]:
        y = df[df.platform == plat].groupby("domain").ratio_vs_bacteria.mean()
        raw = {k: target[k] / y[k] for k in target}
        tot = sum(raw.values())
        for k in target:
            corr.append(dict(platform=plat, domain=k,
                             current_design_abundance_pct=target[k],
                             observed_read_yield_vs_bacteria=y[k],
                             required_copy_abundance_pct=100 * raw[k] / tot,
                             scale_factor=(100 * raw[k] / tot) / target[k]))
    cdf = pd.DataFrame(corr)
    cdf.to_csv(INT / "abundance_correction_required.tsv", sep="\t", index=False)
    print("\nCopy abundances required to hit the intended READ shares:")
    print(cdf.round(4).to_string(index=False))
    print("\nNote: the required abundance file differs between platforms, so a "
          "single abundance file cannot give both platforms a 12% viral read "
          "share. Specifying per-member copy abundance as "
          "(target read share) / (genome length) is the general fix.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
