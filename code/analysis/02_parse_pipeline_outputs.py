#!/usr/bin/env python3
"""Parse every pipeline result file into one long-format detection table.

Output: interim/detections_long.tsv with columns
    platform, pipeline, sample, source_file, parser, taxid, name, rank_reported, reads, extra

Only *reads/contig-level taxonomic assignments* are parsed here; scoring
(TP/FP/FN) happens in 03_score.py. Nothing is filtered out at this stage
beyond restricting to the rank the pipeline itself reports.
"""
import csv
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
OUT = ROOT / "analysis_708" / "interim"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLES = ["sample_0", "sample_1", "sample_2"]
PLATFORMS = {"Illumina": "Illumina", "Pacbio": "PacBio"}

rows = []
notes = []


def add(platform, pipeline, sample, src, parser, taxid, name, rank, reads, extra=""):
    rows.append(dict(platform=platform, pipeline=pipeline, sample=sample,
                     source_file=str(src.relative_to(ROOT)), parser=parser,
                     taxid=taxid, name=name, rank_reported=rank, reads=reads, extra=extra))


def note(msg):
    notes.append(msg)
    print("NOTE:", msg)


# ---------------------------------------------------------------- kraken report
def parse_kraken_report(path, platform, pipeline, sample, parser, ranks=None):
    """Standard Kraken2 report: pct, clade_reads, direct_reads, rank, taxid, name."""
    n = 0
    with open(path, newline="") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            pct, clade, direct, rank, taxid, name = parts[:6]
            rank = rank.strip()
            if rank in ("U", "R") or (ranks and rank not in ranks):
                continue
            try:
                reads = float(clade)
            except ValueError:
                continue
            add(platform, pipeline, sample, path, parser, taxid.strip(),
                name.strip(), rank, reads, extra=f"direct={direct.strip()}")
            n += 1
    return n


# ---------------------------------------------------------------- per pipeline
def do_clark(platform, pdir, sample):
    sdir = pdir / sample
    # primary abundance file = the one WITHOUT the .g0.2 confidence filter
    cands = sorted(p for p in sdir.glob("*abundance*.csv") if ".g0.2." not in p.name)
    filt = sorted(sdir.glob("*.g0.2.csv"))
    for path, tag in [(c, "clark_full") for c in cands] + [(f, "clark_g0.2") for f in filt]:
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh):
                if not r.get("Name"):
                    continue
                try:
                    reads = float(r.get("Count") or 0)
                except ValueError:
                    reads = 0.0
                add(platform, "CLARK", sample, path, tag, (r.get("TaxID") or "").strip(),
                    r["Name"].strip(), "species", reads, extra=r.get("Lineage", ""))


def do_gottcha2(platform, pdir, sample):
    sdir = pdir / sample
    for path in sorted(sdir.glob("*.tsv")):
        if path.name.endswith(".full.tsv") or path.name.endswith(".lineage.tsv"):
            continue
        df = pd.read_csv(path, sep="\t", dtype=str)
        if "LEVEL" not in df.columns:
            note(f"unexpected GOTTCHA2 columns in {path}")
            continue
        for _, r in df.iterrows():
            add(platform, "GOTTCHA2", sample, path, "gottcha2_tsv", str(r["TAXID"]).strip(),
                str(r["NAME"]).strip(), str(r["LEVEL"]).strip(),
                float(r["READ_COUNT"]), extra=f"rel_abund={r.get('REL_ABUNDANCE','')}")


def do_kmr(platform, pdir, sample):
    """KrakenMetaReads-nf: assembly-based Kraken2 reports."""
    sdir = pdir / sample
    if platform == "Illumina":
        merged = list(sdir.glob("*merged_report.csv"))
        if merged:
            df = pd.read_csv(merged[0], dtype=str)
            for _, r in df.iterrows():
                tot = pd.to_numeric(r.get("total_reads"), errors="coerce")
                add(platform, "KrakenMetaReads-nf", sample, merged[0], "kmr_merged",
                    str(r["tax_id"]).strip(), str(r["name"]).strip(), str(r["rank"]).strip(),
                    float(tot) if pd.notna(tot) else 0.0,
                    extra=f"spades={r.get('spades_reads')};megahit={r.get('megahit_reads')}")
        for path in sorted(sdir.glob("*_report.txt")):
            if "merged" in path.name:
                continue
            asm = "spades" if "spades" in path.name else "megahit"
            parse_kraken_report(path, platform, "KrakenMetaReads-nf", sample,
                                f"kmr_{asm}")
    else:
        for path in sorted(sdir.glob("*_flye_report.txt")):
            parse_kraken_report(path, platform, "KrakenMetaReads-nf", sample,
                                "kmr_flye", ranks=None)


def do_mag(platform, pdir, sample):
    sdir = pdir / sample
    if platform == "Illumina":
        for path in sorted(sdir.glob("*kraken2_report.txt")):
            parse_kraken_report(path, platform, "mag-nf", sample, "mag_reads", ranks=None)
    else:
        for path in sorted(sdir.glob("*contigs_kraken2_report.txt")):
            parse_kraken_report(path, platform, "mag-nf", sample, "mag_contigs", ranks=None)
        for path in sorted(sdir.glob("*reads_kraken2_report.txt")):
            parse_kraken_report(path, platform, "mag-nf", sample, "mag_reads_SUSPECT",
                                ranks=None)


def do_taxprofiler(platform, pdir, sample):
    sdir = pdir / sample
    for path in sorted(sdir.glob("*abundance.tsv")):
        df = pd.read_csv(path, sep="\t", dtype=str)
        for _, r in df.iterrows():
            add(platform, "TaxProfiler", sample, path, "taxprofiler_abundance",
                str(r["Taxonomy_ID"]).strip(), str(r["Species"]).strip(), "species",
                float(pd.to_numeric(r["Assigned_Reads"], errors="coerce") or 0),
                extra=f"fraction={r.get('Fraction','')}")
    # NOTE: newer nf-core/taxprofiler runs also emit "bracken_Viral_ref.tsv" and
    # "kraken2_Viral_ref-bracken.tsv" alongside the real per-sample bracken file.
    # Those two lack the taxonomy_id/name/taxonomy_lvl/new_est_reads columns this
    # parser expects, so the glob must match the exact "_null_Viral_ref.bracken.tsv"
    # suffix rather than the loose "*bracken.tsv".
    for path in sorted(sdir.glob("*_null_Viral_ref.bracken.tsv")):
        df = pd.read_csv(path, sep="\t", dtype=str)
        for _, r in df.iterrows():
            add(platform, "TaxProfiler", sample, path, "taxprofiler_bracken",
                str(r["taxonomy_id"]).strip(), str(r["name"]).strip(),
                str(r["taxonomy_lvl"]).strip(),
                float(pd.to_numeric(r["new_est_reads"], errors="coerce") or 0))
    for path in sorted(sdir.glob("*kraken2.report.txt")):
        parse_kraken_report(path, platform, "TaxProfiler", sample,
                            "taxprofiler_kraken2", ranks=None)


def do_taxtriage(platform, pdir, sample):
    sdir = pdir / sample
    for path in sorted(sdir.glob("*.kraken2.report.txt")):
        parse_kraken_report(path, platform, "TaxTriage", sample, "taxtriage_kraken2",
                            ranks=None)
    for path in sorted(sdir.glob("krakenreport.merged_mqc.tsv")):
        df = pd.read_csv(path, sep="\t", dtype=str)
        for _, r in df.iterrows():
            add(platform, "TaxTriage", sample, path, "taxtriage_mqc", "",
                str(r["name"]).strip(), str(r["rank"]).strip(),
                float(pd.to_numeric(r["clade_fragments_covered"], errors="coerce") or 0),
                extra=f"abundance={r.get('abundance','')}")


def do_sourmash(platform, pdir, sample):
    sdir = pdir / sample
    for path in sorted(sdir.glob("*.summarized.csv")):
        db = "viral" if ".viral." in path.name else "prok"
        df = pd.read_csv(path, dtype=str)
        sub = df[df["rank"] == "species"]
        for _, r in sub.iterrows():
            lin = str(r["lineage"])
            name = lin.split(";")[-1].strip()
            add(platform, "sourmash", sample, path, f"sourmash_{db}", "", name, "species",
                float(pd.to_numeric(r["f_weighted_at_rank"], errors="coerce") or 0) * 1e6,
                extra=f"fraction={r['fraction']};lineage={lin}")


HANDLERS = {
    "CLARK": do_clark,
    "GOTTCHA2": do_gottcha2,
    "KrakenMetaReads-nf": do_kmr,
    "mag-nf": do_mag,
    "TaxProfiler": do_taxprofiler,
    "TaxTriage": do_taxtriage,
    "sourmash": do_sourmash,
}


def main():
    for pdir_name, platform in PLATFORMS.items():
        base = ROOT / pdir_name
        for pipeline, fn in HANDLERS.items():
            pdir = base / pipeline
            if not pdir.is_dir():
                note(f"missing pipeline dir {pdir}")
                continue
            for sample in SAMPLES:
                if not (pdir / sample).is_dir():
                    note(f"missing {platform}/{pipeline}/{sample}")
                    continue
                fn(platform, pdir, sample)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "detections_long.tsv", sep="\t", index=False)
    print("\nrows:", len(df))
    print(df.groupby(["platform", "pipeline", "parser"]).size().to_string())
    Path(OUT / "parse_notes.txt").write_text("\n".join(notes))
    return 0


if __name__ == "__main__":
    sys.exit(main())
