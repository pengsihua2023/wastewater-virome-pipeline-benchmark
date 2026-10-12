#!/usr/bin/env python3
"""Extract the CAMISIM ground-truth species list from all_species_list_for_simulation_data.docx.

Writes:
  interim/ground_truth_raw.tsv     -- every table row, verbatim
  03_ground_truth_708_viruses.xlsx -- standardised ground truth workbook
Does not modify any input file.
"""
import re
import sys
from pathlib import Path

import pandas as pd
from docx import Document

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis_708"
DOCX = ROOT / "all_species_list_for_simulation_data.docx"

# Table index -> (section label, category)
SECTIONS = {
    1: ("Section 1: Bacteria/Archaea", "Bacteria/Archaea"),
    2: ("Section 2: Bacteriophage", "Bacteriophage"),
    3: ("Section 3: Viruses - automated selection", "Virus (automated)"),
    4: ("Section 4.1: Viruses - manually added, excluded by rank filter", "Virus (manual, no-rank)"),
    5: ("Section 4.2: Viruses - manually added, excluded by read-count filter", "Virus (manual, low-read)"),
    6: ("Section 4.3: Viruses - manually added, not detected in sample", "Virus (manual, not detected)"),
}


def cell_rows(table):
    rows = []
    for r in table.rows:
        rows.append([c.text.strip() for c in r.cells])
    return rows


def to_num(x):
    if x is None:
        return None
    s = str(x).replace("%", "").replace(",", "").strip()
    m = re.match(r"^-?\d+(\.\d+)?", s)
    return float(m.group(0)) if m else None


def main():
    doc = Document(DOCX)
    records = []
    for ti, (section, category) in SECTIONS.items():
        t = doc.tables[ti]
        rows = cell_rows(t)
        header = rows[0]
        for row in rows[1:]:
            d = dict(zip(header, row))
            name = d.get("Name") or d.get("Virus")
            if not name:
                continue
            # design reads column varies slightly in wording
            design = next((v for k, v in d.items() if "design" in k.lower()), None)
            rel = next((v for k, v in d.items() if "relative abundance" in k.lower()), None)
            # actual reads: prefer clade, else direct, else generic
            actual_clade = next((v for k, v in d.items()
                                 if "actual reads" in k.lower() and "clade" in k.lower()), None)
            actual_direct = next((v for k, v in d.items()
                                  if "actual reads" in k.lower() and "direct" in k.lower()), None)
            actual_any = next((v for k, v in d.items() if "actual reads" in k.lower()), None)
            records.append({
                "section": section,
                "category": category,
                "raw_name": name,
                "taxid": d.get("NCBI taxid", "").strip(),
                "ncbi_rank": d.get("NCBI rank", "").strip(),
                "actual_reads_clade": actual_clade,
                "actual_reads_direct": actual_direct,
                "actual_reads_any": actual_any,
                "design_reads": design,
                "relative_abundance_pct": rel,
                "parent_species_node": d.get("Parent species node", "").strip(),
                "rationale": d.get("Rationale", "").strip(),
            })

    df = pd.DataFrame(records)
    df["design_reads_num"] = df["design_reads"].map(to_num)
    df["relative_abundance_num"] = df["relative_abundance_pct"].map(to_num)
    df["taxid_num"] = pd.to_numeric(df["taxid"], errors="coerce").astype("Int64")

    # --- standardisation -------------------------------------------------
    def norm(n):
        """Drop a trailing parenthetical *acronym* (e.g. "... (hMPV)") but keep
        parentheticals that carry taxonomic information such as influenza
        subtypes (H3N2) or strain designations."""
        s = n.strip()
        m = re.search(r"\s*\(([^()]{1,40})\)\s*$", s)
        if m:
            inner = m.group(1).strip()
            keep = bool(re.fullmatch(r"H\d+N\d+", inner, flags=re.I)) or "/" in inner
            if not keep:
                s = s[: m.start()].strip()
        return s

    def key(n):
        s = n.lower()
        s = re.sub(r"[^a-z0-9]+", " ", s)
        return re.sub(r"\s+", " ", s).strip()

    df["std_name"] = df["raw_name"].map(norm)
    df["match_key"] = df["std_name"].map(key)

    # is this entry viral?
    df["is_viral"] = df["category"] != "Bacteria/Archaea"

    # rank class actually used by the source table
    def rank_class(r):
        r = (r or "").strip().lower()
        if r in ("s", "species"):
            return "species"
        if r == "s1":
            return "strain (S1)"
        if r == "serotype":
            return "serotype"
        if r == "no rank":
            return "no rank"
        return r or "unspecified"

    df["rank_class"] = df["ncbi_rank"].map(rank_class)

    OUT.mkdir(exist_ok=True)
    (OUT / "interim").mkdir(exist_ok=True)
    df.to_csv(OUT / "interim" / "ground_truth_raw.tsv", sep="\t", index=False)

    # --- audit -----------------------------------------------------------
    audit = []
    audit.append(("total rows extracted", len(df)))
    audit.append(("unique raw names", df["raw_name"].nunique()))
    audit.append(("unique taxids", df["taxid_num"].nunique()))
    audit.append(("unique match keys", df["match_key"].nunique()))
    for cat, n in df["category"].value_counts().items():
        audit.append((f"category: {cat}", n))
    audit.append(("viral entries (all non-bacterial)", int(df["is_viral"].sum())))
    audit.append(("bacteria/archaea entries", int((~df["is_viral"]).sum())))
    for rc, n in df["rank_class"].value_counts().items():
        audit.append((f"rank class: {rc}", n))
    dup_name = df[df.duplicated("raw_name", keep=False)]
    dup_tax = df[df["taxid_num"].notna() & df.duplicated("taxid_num", keep=False)]
    audit.append(("duplicate raw names (rows)", len(dup_name)))
    audit.append(("duplicate taxids (rows)", len(dup_tax)))
    audit.append(("total design reads", int(df["design_reads_num"].fillna(0).sum())))
    audit.append(("viral design reads", int(df.loc[df.is_viral, "design_reads_num"].fillna(0).sum())))
    audit.append(("sum of relative abundance (%)", round(df["relative_abundance_num"].fillna(0).sum(), 6)))
    audit_df = pd.DataFrame(audit, columns=["metric", "value"])

    with pd.ExcelWriter(OUT / "03_ground_truth_708_viruses.xlsx", engine="openpyxl") as xw:
        cols = ["section", "category", "raw_name", "std_name", "match_key", "taxid", "taxid_num",
                "ncbi_rank", "rank_class", "is_viral", "actual_reads_clade", "actual_reads_direct",
                "design_reads", "design_reads_num", "relative_abundance_pct",
                "relative_abundance_num", "parent_species_node", "rationale"]
        df[cols].to_excel(xw, sheet_name="ground_truth_all_708", index=False)
        df.loc[df.is_viral, cols].to_excel(xw, sheet_name="viral_subset", index=False)
        df.loc[~df.is_viral, cols].to_excel(xw, sheet_name="bacteria_archaea_subset", index=False)
        audit_df.to_excel(xw, sheet_name="audit", index=False)
        if len(dup_name):
            dup_name[cols].to_excel(xw, sheet_name="duplicate_names", index=False)
        if len(dup_tax):
            dup_tax[cols].to_excel(xw, sheet_name="duplicate_taxids", index=False)

    print(audit_df.to_string(index=False))
    print("\nwrote", OUT / "03_ground_truth_708_viruses.xlsx")
    return 0


if __name__ == "__main__":
    sys.exit(main())
