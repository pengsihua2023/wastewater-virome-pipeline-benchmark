# Analysis summary

## Scope

Seven metagenomic classification pipelines were benchmarked for **viral
detection** on six simulated wastewater datasets (three Illumina, three PacBio
HiFi), against a ground truth of **151 non-phage viral species** drawn from the
metatranscriptomic profile of Los Angeles wastewater sample HTP-2024-01-15.

All reported metrics are computed within the viral space only. The 490
bacterial/archaeal species and 67 ssRNA phage taxa in the simulated community
are background: they are excluded from the ground truth, and their false
positives are excluded from precision. An earlier pass that scored all 708
members together is retained in `_superseded/` and must not be quoted — it
overstates viral capability by 1.5–7.1× on Illumina and ranks a pipeline first
on PacBio that detects no viruses at all.

## Deliverables

| File | Contents |
|---|---|
| `../Seven-pipleline-comparison-on-Illumina.md` / `.docx` | Illumina manuscript, 10 tables, 5 figures |
| `../Seven-pipleline-comparison-on-Pacbio.md` / `.docx` | PacBio manuscript, 10 tables, 5 figures |
| `../Seven-pipeline-comparison-combined.md` / `.docx` | cross-modality manuscript, 7 tables, 7 figures |
| `11_viral_only_metrics.xlsx` | viral metrics, threshold sweep, per-virus error analysis |
| `12_genus_level_sensitivity.xlsx` | strict species vs lenient genus matching |
| `13_detection_vs_depth.xlsx` | Illumina: detection vs reads each virus actually contributed |
| `01_file_inventory.xlsx` | all 302 input files, MD5, role, exclusion reason |
| `02_data_audit_report.md` | data-structure audit, anomalies, analysis limits |
| `03_ground_truth_708_viruses.xlsx` | standardised community and viral ground truth |
| `04_name_mapping_and_manual_review.xlsx` | matching rules, name mapping, review queue |
| `09_designed_vs_observed_reads.xlsx` | designed vs delivered reads per member |
| `10_simulation_config_audit.xlsx` | per-genome simulation audit (depth floor, config) |
| `missing_information_and_questions.md` | what the authors still need to confirm |
| `figures_illumina/`, `figures_pacbio/`, `figures_combined/` | current figures, 600 dpi PNG + vector PDF |
| `scripts/` | 20 Python scripts, run 01→20 |
| `_superseded/` | first-pass whole-community analysis, provenance only |

## Headline results

**Illumina.** Viral recall 0.113–0.309 over the full 151-member ground truth —
but that ceiling is set by the data, not the classifiers: each library carries
only ~3,300 viral reads, so **94 of 151 viruses are expected to contribute fewer
than one read**. Restricted to the 25 viruses with ≥5 expected reads, recall is
**0.51–0.93** (sourmash 0.93; mag-nf, taxprofiler, taxtriage 0.84; GOTTCHA2
0.75; KrakenMetaReads-nf 0.61; CLARK 0.51). Only *Orthopicobirnavirus hominis*
(32 expected reads) was missed by every pipeline despite adequate signal. Viral precision 0.008–0.203, so even the most
conservative pipeline reported four spurious viruses for every real one. Best
viral F1 0.217 (nf-core/taxprofiler, nf-core/taxtriage). The union of all seven
pipelines recovered 78/151 (51.7%). **Ten of the thirteen public-health sentinel
viruses were detected by no pipeline in any sample**, including every influenza
subtype, poliovirus 1 and monkeypox virus.

**PacBio HiFi.** Two groups with nothing in between. nf-core/taxprofiler and
nf-core/taxtriage reached viral recall 0.909 and 0.903 (F1 0.604, 0.601) and
recovered all 13 sentinel viruses; CLARK reached 0.506; and **GOTTCHA2,
nf-core/mag and sourmash detected zero of the 151 viruses in every sample** while
producing normal-looking output. KrakenMetaReads-nf recovered four while
reporting 1,278 spurious viral species.

**Cross-modality.** Viral F1 rank correlation between the two modalities is
Spearman ρ = +0.43 (p = 0.34). Three pipelines went from detecting viruses to
detecting none; CLARK went from last to third. **Viral detection performance on
one modality does not predict performance on the other.**

**Generalises across both modalities.** (i) No pipeline combined viral recall
with viral precision. (ii) **Every read-count threshold above 1 read reduced
viral F1** — at ≥10 reads the best F1 fell from 0.217 to 0.134 (Illumina) and
0.604 to 0.036 (PacBio) — the opposite of standard metagenomic practice, because
viral detections rest on a median of 3–8 reads.

**Genus-level sensitivity analysis.** Lenient genus matching at least halved the
false-positive count for every pipeline, confirming that much of the strict
false-positive burden is ICTV renaming rather than misclassification. But genus
recall capped at 0.508, and it did not rescue the three PacBio zero-recall
pipelines, so the detection gap is real.

## Constraints on the conclusions

1. **Viral read depth is low** — ~3,300 reads per Illumina library and ~950 per
   PacBio library across 151 viruses. All viral results are lower bounds.
2. **The two modalities are not the same viral community in read space.**
   CAMISIM consumed abundance as genome copies while the PacBio driver treated it
   as a share of bases, and only the PacBio run applied a minimum-depth floor.
   No claim is made about which platform detects more viruses.
3. **The PacBio depth floor** raised 64.9% of viral members to ~5 reads
   regardless of designed abundance, so PacBio recall at the low-abundance end,
   and the recovery of all 13 sentinel viruses, are not measures of sensitivity
   at the designed abundances.
4. **No true negatives** are definable, so specificity, accuracy and MCC are not
   reported.
5. **n = 3** simulation replicates; no significance testing.
6. nf-core/taxtriage is not independently evaluated on Illumina (byte-identical
   Kraken2 report to nf-core/taxprofiler); nf-core/mag was scored on contigs on
   PacBio; GOTTCHA2 ran a locally patched binary on PacBio.

## Recommended next steps

1. **Re-simulate at the intended viral read share.** Copy abundances must be set
   to (target read share) ÷ (genome length); block-level targets are in
   `interim/abundance_correction_required.tsv`. Remove the PacBio `--min_reads`
   floor, or renormalise after applying it.
2. **Confirm or regenerate** the PacBio nf-core/mag read-level reports and the
   Illumina TaxTriage run, so both are independently evaluated.
3. **Report viral recall against viral read depth**, not against the platform
   label, so the modality comparison becomes interpretable.
4. **Record the missing versions** (CAMISIM, GOTTCHA2, CLARK, database build
   dates) — see `missing_information_and_questions.md` §B3.
5. **Do not adopt a read-count threshold above 1 read** for viral surveillance at
   these depths without re-deriving it from viral depth.
