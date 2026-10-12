# Statistical analysis — methods, effect sizes and what was deliberately not computed

Covers the **virus-only** analysis that the three manuscripts report. All
numbers are produced by `scripts/14_viral_only_scoring.py`,
`15_genus_level_sensitivity.py` and `16_viral_stratified_tables.py` from
`interim/viral_detections_long.tsv`.

---

## 1. Unit of analysis and design

- **Unit**: one (platform, pipeline, sample) run. 2 modalities × 7 pipelines ×
  3 samples = **42 runs**.
- **Replicates**: n = 3 per modality × pipeline cell. These are three
  independently simulated samples from the same design — *simulation*
  replicates, not biological replicates.
- **Design**: fully crossed; every pipeline was run on every sample of both
  modalities.
- **Ground truth**: 151 non-phage viral species.

## 2. Metric definitions

Every reported node is first classified as viral or non-viral (indentation-
tracked descent from the `Viruses` node for Kraken2 reports; the `Lineage`
column for CLARK; the lineage string for sourmash; the sibling `.lineage.tsv`
for GOTTCHA2). Validation confirmed that no ground-truth virus was misclassified
as non-viral.

Within the viral space, per run:

- **TP** — a ground-truth virus whose NCBI taxid appears at a species-or-below
  node (S, S1–S5, species, strain, serotype, no rank, subspecies, biotype,
  isolate, genotype) with clade reads above the threshold.
- **FP** — a reported species-or-below node **inside the Viruses branch** that
  is not in the 151-member ground truth.
- **FN** — a ground-truth virus with no such node.
- Recall = TP/(TP+FN) with TP+FN = 151; Precision = TP/(TP+FP); F1 = 2PR/(P+R).

Two exclusions. Detections of the 490 bacterial/archaeal background taxa are
ignored, as are their false positives: they bear on no question these
manuscripts ask. The 67 ssRNA phage members of the background are removed from
the reported node set, so they count neither as TP nor as FP; metrics against
the alternative 218-member viral ground truth including phage are in
`11_viral_only_metrics.xlsx`, sheet `all_thresholds_both_GT`.

### 2.1 Metrics deliberately not computed

**Specificity, accuracy, balanced accuracy and MCC are not reported.** All
require true negatives — a defined set of viruses that could have been reported
and correctly were not. That set is the viral portion of each tool's reference
database: the Kraken2 Standard 2024-06-05 build, the CLARK species target set,
the GOTTCHA2 signature database and NCBI-viruses 2025.01 all differ in size, and
none was supplied with the project. Any TN count would be an arbitrary choice
that would dominate the resulting statistic.

### 2.2 sourmash is scored on a different basis

sourmash emits GTDB and NCBI-viral lineage strings without NCBI taxids. It is
matched on normalised species name after stripping GTDB rank prefixes, and is
flagged `match_mode = name` throughout. Its recall is a **lower bound** and is
not strictly comparable with the six taxid-matched pipelines.

## 3. Detection-threshold sweep

Every metric was recomputed at **≥1, ≥10 and ≥100 clade reads**. The
manuscripts report ≥1 read as the primary operating point and give the full
sweep, because on viruses the threshold does not behave as it does on whole
communities: at these viral depths every threshold above 1 read *reduced* viral
F1 on both modalities. With a median of 2 viral reads per
detected virus on Illumina and 3 on PacBio, a read-count filter removes true
positives faster than false ones.

## 4. Uncertainty

With n = 3 per cell, **no significance testing is performed**. Reported instead:

- the three individual sample values (never only the mean);
- the mean and sample standard deviation;
- where useful, a non-parametric bootstrap 95% interval of the mean
  (10,000 resamples, `numpy.random.default_rng(42)`).

A bootstrap interval from three observations is itself unstable; it describes
the observed spread, not an inferential claim. Between-sample variability was
small in every cell (viral recall SD ≤ 0.023 on Illumina, ≤ 0.020 on PacBio), so
the *ranking* is stable across the three simulation replicates even though no
formal test supports it.

## 5. Cross-modality comparison

Illumina and PacBio `sample_N` directories correspond by index and derive from
the same abundance files, but the two simulators interpret those files
differently and only the PacBio run applied a minimum-depth floor
(`02_data_audit_report.md`, `10_simulation_config_audit.xlsx`). The two
modalities therefore do **not** carry the same viral community in read space.

Consequently:

- No paired significance test is claimed, and no statement is made about which
  platform detects more viruses.
- Cross-modality statements concern the behaviour of a given pipeline when its
  input type changes, which is well-posed regardless of dataset equivalence.
- Rank agreement is summarised by Spearman correlation over the seven pipelines:
  **viral F1 ρ = +0.09 (p = 0.84)**, viral recall ρ = −0.12 (p = 0.81), viral
  precision ρ = +0.17 (p = 0.72), all n = 7. None is distinguishable from no
  association at this sample size, and the substantive argument rests on the
  magnitude of individual movements (three pipelines to zero) rather than on the
  coefficient.

## 6. Genus-level sensitivity analysis

Strict taxid equality counts *Mamastrovirus hominis* as a false positive against
ground-truth *Mamastrovirus 1*, and *Alphainfluenzavirus influenzae* against the
influenza A serotype nodes. Because viral nomenclature is under large-scale ICTV
revision, strict matching systematically overstates false positives.

A taxid → genus map was built by indentation-tracked descent within the Viruses
branch of the Kraken2 report with the widest coverage (23,059 viral nodes,
2,640 genera). A genus-level TP is a ground-truth virus for which any
species-or-below node in the same genus was reported; a genus-level FP is a
reported viral genus containing no ground-truth member.

**125 of 151 ground-truth viruses had a resolvable genus**, defining 63
ground-truth genera; the analysis therefore covers 83% of the ground truth, and
the 26 unresolved viruses are listed in `12_genus_level_sensitivity.xlsx`.
Main results remain the strict species-level ones.

## 7. Abundance stratification

Viral recall was stratified by designed relative abundance into four decades
(< 1e-4%, 1e-4–1e-3%, 1e-3–1e-2%, ≥ 1e-2%; n = 15, 75, 43, 18). **Bins contain
only viral ground-truth members**, so denominators differ from the
whole-community binning used in the superseded first pass.

Cells in which a pipeline detected nothing are included as zeros before
averaging over the three samples. Omitting them — as an earlier draft of this
analysis did — inflates recall in sparse bins; the correction changed, for
example, CLARK's Illumina lowest-abundance recall from 0.133 to 0.044.

## 8. What is not claimed

1. No pipeline is described as "best" without naming the metric, modality and
   detection threshold.
2. No difference is called significant; none was tested.
3. Viral precision differences between tools with different viral database
   scope are not interpreted as differences in classifier stringency alone.
4. No statement is made about which sequencing platform detects more viruses.
5. PacBio recall in the low-abundance bins, and the recovery of all 13 sentinel
   viruses on PacBio, are not presented as sensitivity at the designed
   abundances: 64.9% of viral members were raised to a ~5-read floor by the
   simulation.
6. All viral results are lower bounds, measured at ~3,300 (Illumina) and ~950
   (PacBio) viral reads per library.
