#!/usr/bin/env python3
"""Build the traceability workbook and convert the manuscript Markdown to .docx."""
import re
import sys
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

ROOT = Path(__file__).resolve().parents[2]
AN = ROOT / "analysis_708"
SUP = AN / "_superseded"          # whole-community scoring, superseded
SUP.mkdir(parents=True, exist_ok=True)

TRACE = [
    # (claim, value, source_file, field, script, formula, intermediate, where_used)
    ("Ground-truth community size", "708 members (490 bacteria/archaea, 67 phage, 151 non-phage virus)",
     "all_species_list_for_simulation_data.docx tables 1-6", "Name / NCBI taxid",
     "01_extract_ground_truth.py", "row count per section table",
     "interim/ground_truth_raw.tsv", "Abstract; M&M 2.2; Table 1"),
    ("Ground-truth uniqueness", "708 unique names, 708 unique taxids, 708 unique match keys",
     "same", "raw_name, taxid, match_key", "01_extract_ground_truth.py",
     "nunique()", "03_ground_truth_708_viruses.xlsx sheet 'audit'", "M&M 2.2"),
    ("Designed composition", "88.00% bacteria/archaea, 8.00% phage, 4.00% virus",
     "same", "Relative abundance in simulation (%)", "01_extract_ground_truth.py",
     "sum by section", "interim/ground_truth_raw.tsv", "M&M 2.2; Fig 10"),
    ("Illumina library size", "16,642,601-16,643,343 reads per sample",
     "Illumina/TaxProfiler/sample_N/*.kraken2.report.txt", "clade reads at U + R",
     "05_realised_composition.py", "unclassified + root", "interim/realised_composition.tsv",
     "M&M 2.3; Table 1"),
    ("PacBio library size", "449,134-518,608 reads per sample",
     "Pacbio/TaxProfiler/sample_N/*.kraken2.report.txt", "clade reads at U + R",
     "05_realised_composition.py", "unclassified + root", "interim/realised_composition.tsv",
     "M&M 2.3; Table 1"),
    ("Abundance scale diagnosis", "viral reads-per-abundance = 1/645 of bacteria (Illumina), predicted 1/571 from genome length; predicted viral read share 0.018% vs 0.021% measured",
     "interim/detections_long.tsv; interim/ground_truth_raw.tsv", "clade reads at taxid 2/2157/10239; relative_abundance_num",
     "06_abundance_scale_diagnosis.py", "(reads_domain / abundance_domain) normalised to bacteria",
     "interim/abundance_scale_diagnosis.tsv", "Abstract; M&M 2.4; Results 3.5; Discussion 4.4, 4.5; Fig 10"),
    ("Realised viral fraction, Illumina", "0.0198-0.0201% of reads",
     "Illumina/TaxProfiler/sample_N/*.kraken2.report.txt", "clade reads at taxid 10239",
     "05_realised_composition.py", "100 * virus_reads / total_reads",
     "interim/realised_composition.tsv", "Abstract; M&M 2.4; Results 3.5; Discussion 4.5; Fig 10"),
    ("Realised viral fraction, PacBio", "0.183-0.222% of reads",
     "Pacbio/TaxProfiler/sample_N/*.kraken2.report.txt", "clade reads at taxid 10239",
     "05_realised_composition.py", "100 * virus_reads / total_reads",
     "interim/realised_composition.tsv", "Abstract; M&M 2.4; Results 3.5; Discussion 4.5; Fig 10"),
    ("Illumina recall, mag-nf", "0.853 +/- 0.002",
     "Illumina/mag-nf/sample_N/*.kraken2_report.txt", "clade reads, rank S/S1-S5",
     "03_score.py -> 04_aggregate_errors_abundance.py", "TP / (TP + FN)",
     "interim/sample_level_metrics.tsv", "Abstract; Results 3.1; Table 3; Fig 2"),
    ("Illumina recall, taxprofiler / taxtriage", "0.852 +/- 0.001 (identical source file)",
     "Illumina/TaxProfiler/sample_N/*.kraken2.report.txt", "clade reads, rank S/S1-S5",
     "03_score.py", "TP / (TP + FN)", "interim/sample_level_metrics.tsv",
     "Results 3.1; M&M 2.6; Table 3"),
    ("Illumina precision, GOTTCHA2", "0.721 +/- 0.006; FP = 119.7 +/- 3.5",
     "Illumina/GOTTCHA2/sample_N/*_01.tsv", "READ_COUNT where LEVEL = species",
     "03_score.py", "TP / (TP + FP)", "interim/sample_level_metrics.tsv",
     "Abstract; Results 3.1; Table 3; Fig 2, Fig 3"),
    ("PacBio recall, taxprofiler", "0.974 +/- 0.006 (689.7 +/- 4.5 of 708)",
     "Pacbio/TaxProfiler/sample_N/*.kraken2.report.txt", "clade reads, rank S/S1-S5",
     "03_score.py", "TP / (TP + FN)", "interim/sample_level_metrics.tsv",
     "Abstract; Results 3.2; Table 3; Fig 2"),
    ("PacBio precision, GOTTCHA2", "0.902 +/- 0.010; FP = 13.0 +/- 1.7",
     "Pacbio/GOTTCHA2/sample_N/*.tsv", "READ_COUNT where LEVEL = species",
     "03_score.py", "TP / (TP + FP)", "interim/sample_level_metrics.tsv",
     "Abstract; Results 3.2; Table 3"),
    ("PacBio F1, mag-nf (highest of study)", "0.596 +/- 0.004",
     "Pacbio/mag-nf/sample_N/*contigs_kraken2_report.txt", "clade reads, rank S/S1-S5",
     "03_score.py", "2PR/(P+R)", "interim/sample_level_metrics.tsv",
     "Abstract; Results 3.2; Table 3"),
    ("KrakenMetaReads-nf PacBio recall", "0.007 +/- 0.001 (5 of 708); FP = 1278 +/- 13",
     "Pacbio/KrakenMetaReads-nf/sample_N/*_flye_report.txt", "clade reads, rank S/S1-S5",
     "03_score.py", "TP / (TP + FN)", "interim/sample_level_metrics.tsv",
     "Abstract; Results 3.2; Discussion 4.3"),
    ("Platform contrast, taxprofiler recall", "+0.122 (bootstrap 95% CI +0.116 to +0.127)",
     "interim/sample_level_metrics.tsv", "recall", "04_aggregate_errors_abundance.py",
     "mean(PacBio_i - Illumina_i), 10000-resample bootstrap, seed 42",
     "interim/platform_contrast.tsv", "Results 3.3; Table 4; Fig 5"),
    ("Domain recall, Illumina non-phage virus", "0.305 (taxprofiler/taxtriage), 0.309 (mag-nf)",
     "interim/detection_detail.tsv", "status, category", "04_aggregate_errors_abundance.py",
     "TP_domain / n_truth_domain", "interim/domain_recall.tsv",
     "Abstract; Results 3.4; Table 5; Fig 8"),
    ("Domain recall, PacBio non-phage virus", "0.909 (taxprofiler), 0.903 (taxtriage); 0.000 for GOTTCHA2, mag-nf, sourmash",
     "interim/detection_detail.tsv", "status, category", "04_aggregate_errors_abundance.py",
     "TP_domain / n_truth_domain", "interim/domain_recall.tsv",
     "Abstract; Results 3.4; Table 5; Fig 8"),
    ("Abundance-stratified recall", "PacBio taxprofiler 0.916-0.985 in every decade; Illumina 0.000 below 1e-4%",
     "interim/detection_detail.tsv", "rel_abundance_pct, status",
     "04_aggregate_errors_abundance.py", "TP_bin / n_truth_bin",
     "interim/abundance_recall.tsv", "Abstract; Results 3.5; Fig 6"),
    ("Abundance correlation", "Spearman rho 0.22-0.77 by pipeline and platform",
     "interim/detection_detail.tsv", "rel_abundance_pct vs observed_reads",
     "04_aggregate_errors_abundance.py", "scipy.stats.spearmanr over true positives",
     "interim/abundance_correlation.tsv", "Results 3.6; Fig 7; 08_statistical_analysis.md"),
    ("Members missed by 36 of 42 runs", "8 non-phage viruses (Orthopicobirnavirus hominis, Enterovirus C, Rotavirus A, Rotavirus C, Avastrovirus 2, Cardiovirus D, Cosavirus D, Rabovirus A)",
     "interim/detection_detail.tsv", "status = FN", "04_aggregate_errors_abundance.py",
     "count of runs per ground-truth member", "interim/missed_summary.tsv; 07_error_analysis.xlsx",
     "Results 3.8; Discussion 4.2; Table 6"),
    ("No member missed by all 42 runs", "0",
     "interim/detection_detail.tsv", "status = FN", "04_aggregate_errors_abundance.py",
     "max(n_runs_missed) = 36 < 42", "07_error_analysis.xlsx sheet 'missed_by_all_runs'",
     "Results 3.8"),
    ("Platform-exclusive viral detections", "70 PacBio-only, 52 Illumina-only, 70 both",
     "interim/detection_detail.tsv", "status = TP, platform, category",
     "04_aggregate_errors_abundance.py", "set membership per ground-truth member",
     "07_error_analysis.xlsx sheet 'platform_exclusive'", "Results 3.8"),
    ("Sentinel virus recovery", "11 of 13 recovered only by read-level Kraken2-Standard pipelines; 9 only on PacBio",
     "interim/detection_detail.tsv", "status = TP, category = Virus (manual, *)",
     "04_aggregate_errors_abundance.py", "platform x pipeline combinations per taxid",
     "07_error_analysis.xlsx", "Abstract; Results 3.8"),
    ("Maximum between-sample recall SD", "0.009 across all 14 platform x pipeline cells",
     "interim/sample_level_metrics.tsv", "recall", "04_aggregate_errors_abundance.py",
     "std(ddof=1) within cell", "interim/aggregate_metrics.tsv",
     "Abstract; Results 3.2; 08_statistical_analysis.md"),
    ("PacBio mag-nf read reports excluded", "3 byte-identical files (MD5 70de11486f26...)",
     "Pacbio/mag-nf/sample_N/*reads_kraken2_report.txt", "whole file",
     "10_file_inventory.py", "md5", "01_file_inventory.xlsx sheet 'byte_identical_files'",
     "M&M 2.6; 02_data_audit_report.md 3.1"),
    ("Illumina taxtriage = taxprofiler", "3 byte-identical Kraken2 report pairs",
     "Illumina/{TaxTriage,TaxProfiler}/sample_N/*.report.txt", "whole file",
     "10_file_inventory.py", "md5", "01_file_inventory.xlsx sheet 'byte_identical_files'",
     "M&M 2.6; Results 3.1; Discussion 4.5"),
    ("Files inventoried", "302, all opened successfully",
     "Illumina/, Pacbio/, project root", "path, size, md5", "10_file_inventory.py",
     "recursive walk", "01_file_inventory.xlsx", "M&M 2.6; 02_data_audit_report.md"),
]


def build_trace():
    df = pd.DataFrame(TRACE, columns=[
        "claim_in_manuscript", "value", "source_file", "field_used", "script",
        "formula", "intermediate_table", "where_used"])
    with pd.ExcelWriter(SUP / "manuscript_traceability.xlsx") as xw:
        df.to_excel(xw, sheet_name="claims_to_sources", index=False)
        pd.DataFrame([
            dict(script="01_extract_ground_truth.py", produces="03_ground_truth_708_viruses.xlsx, interim/ground_truth_raw.tsv"),
            dict(script="02_parse_pipeline_outputs.py", produces="interim/detections_long.tsv (587,824 rows)"),
            dict(script="03_score.py", produces="interim/sample_level_metrics.tsv, detection_detail.tsv, sensitivity_streams.tsv"),
            dict(script="04_aggregate_errors_abundance.py", produces="05/06/07 workbooks and interim aggregates"),
            dict(script="05_realised_composition.py", produces="interim/realised_composition.tsv"),
            dict(script="09_whole_community_figures.py", produces="figures/Figure1-10 (PNG 600 dpi + PDF)"),
            dict(script="10_file_inventory.py", produces="01_file_inventory.xlsx"),
            dict(script="11_name_mapping_and_tables.py", produces="04_name_mapping_and_manual_review.xlsx, tables/"),
            dict(script="12_traceability_and_docx.py", produces="manuscript_traceability.xlsx, manuscript-708-viruses-draft.docx"),
            dict(script="13_consistency_check.py", produces="verification of 46 headline numbers against interim tables"),
            dict(script="06_abundance_scale_diagnosis.py", produces="interim/abundance_scale_diagnosis.tsv"),
        ]).to_excel(xw, sheet_name="script_provenance", index=False)
    print("wrote manuscript_traceability.xlsx")


# ----------------------------------------------------------------- md -> docx
def add_runs(par, text):
    """Render **bold**, *italic* and `code` inline spans."""
    for tok in re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)", text):
        if not tok:
            continue
        if tok.startswith("**") and tok.endswith("**"):
            par.add_run(tok[2:-2]).bold = True
        elif tok.startswith("*") and tok.endswith("*"):
            par.add_run(tok[1:-1]).italic = True
        elif tok.startswith("`") and tok.endswith("`"):
            r = par.add_run(tok[1:-1])
            r.font.name = "Consolas"
            r.font.size = Pt(9)
        else:
            par.add_run(tok)


def md_to_docx(md_path, out_path, figures_dir):
    lines = md_path.read_text().split("\n")
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(11)
    inserted_figures = set()

    i = 0
    while i < len(lines):
        line = lines[i]
        s = line.rstrip()
        if not s.strip():
            i += 1
            continue
        if s.strip() == "---":
            i += 1
            continue
        # table block
        if s.lstrip().startswith("|"):
            block = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                block.append(lines[i].strip())
                i += 1
            rows = [[c.strip() for c in r.strip("|").split("|")] for r in block]
            rows = [r for r in rows if not all(set(c) <= set("-: ") for c in r)]
            if rows:
                t = doc.add_table(rows=len(rows), cols=max(len(r) for r in rows))
                t.style = "Table Grid"
                for ri, r in enumerate(rows):
                    for ci, c in enumerate(r):
                        cell = t.cell(ri, ci)
                        cell.text = ""
                        add_runs(cell.paragraphs[0], c)
                        for run in cell.paragraphs[0].runs:
                            run.font.size = Pt(8)
                            if ri == 0:
                                run.bold = True
            doc.add_paragraph()
            continue
        m = re.match(r"^(#{1,4})\s+(.*)$", s)
        if m:
            doc.add_heading(m.group(2), level=min(len(m.group(1)), 4))
            i += 1
            continue
        m = re.match(r"^[-*]\s+(.*)$", s.strip())
        if m:
            p = doc.add_paragraph(style="List Bullet")
            add_runs(p, m.group(1))
            i += 1
            continue
        m = re.match(r"^(\d+)\.\s+(.*)$", s.strip())
        if m:
            p = doc.add_paragraph(style="List Number")
            add_runs(p, m.group(2))
            i += 1
            continue
        # A "**Figure N. ...**" caption line: insert the actual image right
        # here, in the body, immediately above its own caption -- rather than
        # only listing captions in a "Figure Legends" section and dumping all
        # images separately at the very end, which left the two disconnected.
        fm = re.match(r"^\*\*Figure\s+(\d+)\.", s.strip())
        if fm:
            n = fm.group(1)
            hits = sorted(figures_dir.glob(f"Figure{n}_*.png"))
            if hits:
                doc.add_picture(str(hits[0]), width=Inches(6.3))
                inserted_figures.add(n)
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        add_runs(p, s.strip())
        i += 1

    # Any figure not referenced by a "**Figure N." caption line in the body
    # (e.g. a markdown file that still uses the old end-of-document Figure
    # Legends convention) is appended here so nothing is silently dropped.
    remaining = [n for n in range(1, 11)
                 if str(n) not in inserted_figures and sorted(figures_dir.glob(f"Figure{n}_*.png"))]
    if remaining:
        doc.add_page_break()
        doc.add_heading("Figures", level=1)
        for n in remaining:
            hits = sorted(figures_dir.glob(f"Figure{n}_*.png"))
            doc.add_heading(f"Figure {n}", level=2)
            doc.add_picture(str(hits[0]), width=Inches(6.3))
    doc.save(out_path)
    print("wrote", out_path)


def main():
    build_trace()
    # The first-pass combined draft has been superseded by the three virus-only
    # manuscripts (see _superseded/README.md). Convert it only if it is present.
    src = SUP / "manuscript-708-viruses-draft.md"
    if src.exists():
        md_to_docx(src, SUP / "manuscript-708-viruses-draft.docx", SUP / "figures")
    else:
        print("superseded draft not present; skipping its .docx conversion")
    return 0


if __name__ == "__main__":
    sys.exit(main())
