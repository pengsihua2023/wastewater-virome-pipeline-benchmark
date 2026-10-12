#!/usr/bin/env python3
"""Realised super-kingdom composition of each simulated dataset.

The design document specifies 88% Bacteria/Archaea, 8% bacteriophage and 4%
virus (all three relative-abundance blocks sum to 100%). This script measures
what each dataset actually contains, using the read-level Kraken2 report that
each platform's TaxProfiler run produced against the same database.
"""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
INT = ROOT / "analysis_708" / "interim"

FILES = {
    ("Illumina", "sample_0"): "Illumina/TaxProfiler/sample_0/sample_708_0_null_Viral_ref.kraken2.kraken2.report.txt",
    ("Illumina", "sample_1"): "Illumina/TaxProfiler/sample_1/sample_708_1_null_Viral_ref.kraken2.kraken2.report.txt",
    ("Illumina", "sample_2"): "Illumina/TaxProfiler/sample_2/sample_708_2_null_Viral_ref.kraken2.kraken2.report.txt",
    ("PacBio", "sample_0"): "Pacbio/TaxProfiler/sample_0/sample_708_0_PacBio_null_Viral_ref.kraken2.kraken2.report.txt",
    ("PacBio", "sample_1"): "Pacbio/TaxProfiler/sample_1/sample_708_1_PacBio_null_Viral_ref.kraken2.kraken2.report.txt",
    ("PacBio", "sample_2"): "Pacbio/TaxProfiler/sample_2/sample_708_2_PacBio_null_Viral_ref.kraken2.kraken2.report.txt",
}


def read_report(path):
    rows = []
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 6:
            rows.append(dict(pct=float(p[0]), clade=int(p[1]), direct=int(p[2]),
                             rank=p[3].strip(), taxid=p[4].strip(), name=p[5].strip()))
    return pd.DataFrame(rows)


def main():
    out = []
    for (plat, samp), rel in FILES.items():
        df = read_report(ROOT / rel)
        tot = df.loc[df["rank"] == "U", "clade"].sum() + df.loc[df["rank"] == "R", "clade"].sum()

        def node(taxid):
            s = df[df.taxid == taxid]
            return int(s.clade.iloc[0]) if len(s) else 0

        unclass = int(df.loc[df["rank"] == "U", "clade"].sum())
        bact, arch, vir = node("2"), node("2157"), node("10239")
        out.append(dict(platform=plat, sample=samp, source_file=rel,
                        total_reads=tot, unclassified_reads=unclass,
                        bacteria_reads=bact, archaea_reads=arch, virus_reads=vir,
                        pct_unclassified=100 * unclass / tot,
                        pct_bacteria=100 * bact / tot, pct_archaea=100 * arch / tot,
                        pct_virus=100 * vir / tot))
    df = pd.DataFrame(out)
    df["design_pct_virus_incl_phage"] = 12.0     # 8% phage + 4% virus
    df["design_pct_bacteria_archaea"] = 88.0
    df["realised_over_design_virus"] = df["pct_virus"] / df["design_pct_virus_incl_phage"]
    df.to_csv(INT / "realised_composition.tsv", sep="\t", index=False)
    print(df.drop(columns=["source_file"]).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
