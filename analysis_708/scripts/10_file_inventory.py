#!/usr/bin/env python3
"""Full recursive inventory of every input file, with read status and analysis use."""
import hashlib
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
AN = ROOT / "analysis_708"

USED_PARSERS = {
    "clark_full": "PRIMARY (CLARK)", "clark_g0.2": "sensitivity (CLARK, g>=0.2 filter)",
    "gottcha2_tsv": "PRIMARY (GOTTCHA2)",
    "kmr_merged": "PRIMARY (KrakenMetaReads-nf, Illumina)",
    "kmr_flye": "PRIMARY (KrakenMetaReads-nf, PacBio)",
    "kmr_spades": "sensitivity (single assembler)",
    "kmr_megahit": "sensitivity (single assembler)",
    "mag_reads": "PRIMARY (mag-nf, Illumina)",
    "mag_contigs": "PRIMARY (mag-nf, PacBio)",
    "mag_reads_SUSPECT": "EXCLUDED - identical file re-used across samples",
    "taxprofiler_kraken2": "PRIMARY (TaxProfiler / TaxTriage)",
    "taxprofiler_bracken": "sensitivity (Bracken re-estimation, Illumina only)",
    "taxprofiler_abundance": "sensitivity (pipeline-derived abundance table)",
    "taxtriage_kraken2": "PRIMARY (TaxTriage)",
    "taxtriage_mqc": "context only (no taxid column)",
    "sourmash_prok": "PRIMARY (sourmash, GTDB prokaryote DB)",
    "sourmash_viral": "PRIMARY (sourmash, viral DB)",
}

CONTENT_HINTS = [
    (".out", "SLURM stdout log", "context only - run provenance, no taxonomic table"),
    (".err", "SLURM stderr log", "context only - run provenance"),
    (".sh", "job submission / driver script", "context only - parameters and versions"),
    (".nf", "Nextflow workflow definition", "context only - pipeline structure"),
    (".config", "Nextflow configuration", "context only - resources and DB paths"),
    (".py", "helper script shipped with the run", "context only"),
    ("samplesheet", "pipeline sample sheet", "context only - input FASTQ paths"),
    ("databases.csv", "nf-core/taxprofiler database sheet", "context only - DB identity"),
]


def md5(p, limit=None):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    det = pd.read_csv(AN / "interim" / "detections_long.tsv", sep="\t",
                      dtype={"taxid": str}, low_memory=False)
    parsed = det.groupby("source_file").agg(rows_parsed=("name", "size"),
                                            parser=("parser", "first")).to_dict("index")

    rows = []
    targets = [p for p in ROOT.rglob("*") if p.is_file()
               and "analysis_708" not in p.parts and not p.name.startswith("~$")]
    for p in sorted(targets):
        rel = p.relative_to(ROOT)
        parts = rel.parts
        platform = {"Illumina": "Illumina", "Pacbio": "PacBio"}.get(parts[0], "n/a (project root)")
        pipeline = parts[1] if len(parts) > 2 else "n/a"
        sample = next((x for x in parts if x.startswith("sample_")), "n/a (pipeline level)")
        info = parsed.get(str(rel))
        if info:
            ftype = "taxonomic result table"
            content = f"{info['rows_parsed']} taxon rows parsed"
            read_ok, included = "yes", USED_PARSERS.get(info["parser"], "yes")
            reason = ""
        else:
            ftype, content, included = "other", "", "no"
            read_ok = "yes"
            reason = "not a taxon-per-row result table"
            for key, t, use in CONTENT_HINTS:
                if p.name.endswith(key) or key in p.name:
                    ftype, content, reason = t, t, use
                    included = "context only"
                    break
            else:
                if p.suffix in (".docx", ".pdf"):
                    ftype = "source document"
                    included = "yes"
                    reason = "ground truth / background"
                elif "full.tsv" in p.name:
                    ftype = "GOTTCHA2 unfiltered signature table"
                    reason = "superseded by the filtered .tsv used as PRIMARY"
                    included = "checked, not scored"
                elif "lineage.tsv" in p.name:
                    ftype = "GOTTCHA2 lineage table"
                    reason = "redundant with the scored .tsv"
                    included = "checked, not scored"
                elif "krona" in p.name:
                    ftype = "Krona input table"
                    reason = "visualisation input, redundant with the Kraken2 report"
                    included = "checked, not scored"
                elif "hierarchy" in p.name:
                    ftype = "TaxTriage merged hierarchy CSV"
                    reason = "byte-identical to the sibling .single.csv"
                    included = "checked, not scored"
                elif "abundance_summary" in p.name or "consensus" in p.name:
                    ftype = "pipeline summary table"
                    reason = "derived summary; primary table scored instead"
                    included = "checked, not scored"
                elif "combined_reports" in p.name:
                    ftype = "multi-sample combined report"
                    reason = "single-sample report scored instead"
                    included = "checked, not scored"
                elif p.name.endswith("_abundance.txt"):
                    ftype = "per-contig abundance table"
                    reason = "contig-level, carries no taxon assignment"
                    included = "checked, not scored"
                elif "single.csv" in p.name:
                    ftype = "TaxTriage per-sample taxon CSV"
                    reason = "no read counts; Kraken2 report scored instead"
                    included = "checked, not scored"
                elif p.name.endswith(".csv") and "assignments" in p.name:
                    ftype = "CLARK per-read assignment file"
                    reason = "read-level, aggregated form scored instead"
                    included = "checked, not scored"
                elif "merged_report.txt" in p.name:
                    ftype = "KrakenMetaReads-nf merged summary (text)"
                    reason = "3-line summary; merged_report.csv scored instead"
                    included = "checked, not scored"
                elif "report_bracken.txt" in p.name:
                    ftype = "Bracken-adjusted Kraken2 report"
                    reason = "sensitivity stream taken from the .bracken.tsv instead"
                    included = "checked, not scored"
                elif "top_viruses_summary" in p.name:
                    ftype = "TaxProfiler top-taxa summary"
                    reason = "truncated view of the abundance table"
                    included = "checked, not scored"
                elif "flye_classification" in p.name:
                    ftype = "per-contig Kraken2 assignment"
                    reason = "contig-level; the flye Kraken2 report is scored instead"
                    included = "checked, not scored"
                elif p.name == "prompt.txt":
                    ftype = "analysis brief supplied by the author"
                    reason = "task specification, not data"
                    included = "context only"
                elif "human.txt" in p.name:
                    ftype = "sourmash human-readable summary"
                    reason = "redundant with summarized.csv"
                    included = "checked, not scored"
                elif "kreport.txt" in p.name:
                    ftype = "sourmash Kraken-style report"
                    reason = "no NCBI taxids; summarized.csv scored instead"
                    included = "checked, not scored"
                elif p.name.endswith(".csv") and "virus_results" in p.name:
                    ftype = "CLARK raw per-read result"
                    reason = "read-level, aggregated abundance file scored instead"
                    included = "checked, not scored"
        rows.append(dict(path=str(rel), platform=platform, pipeline=pipeline,
                         sample=sample, size_bytes=p.stat().st_size,
                         file_type=ftype, main_content=content or ftype,
                         read_successfully=read_ok, included_in_analysis=included,
                         reason_if_excluded=reason,
                         md5=md5(p) if p.stat().st_size < 200_000_000 else "not hashed (>200 MB)"))

    df = pd.DataFrame(rows)
    dup = df[df.md5.duplicated(keep=False) & (df.md5 != "not hashed (>200 MB)")]
    dup = dup.sort_values("md5")

    with pd.ExcelWriter(AN / "01_file_inventory.xlsx") as xw:
        df.to_excel(xw, sheet_name="all_files", index=False)
        dup.to_excel(xw, sheet_name="byte_identical_files", index=False)
        (df.groupby(["platform", "pipeline"]).size().rename("n_files").reset_index()
         .to_excel(xw, sheet_name="counts_by_pipeline", index=False))
    print("files inventoried:", len(df))
    print(df.included_in_analysis.value_counts().to_string())
    print("byte-identical file rows:", len(dup))
    return 0


if __name__ == "__main__":
    sys.exit(main())
