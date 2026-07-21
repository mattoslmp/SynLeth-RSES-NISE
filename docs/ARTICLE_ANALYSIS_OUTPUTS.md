# Empirical outputs for the RSES-Onco manuscript

This document maps reproducible outputs to manuscript results, figures and
supplementary tables. It does not prescribe conclusions before the real-data
workflow is complete.

## Main empirical tables

### `Table_main_top_rses_onco_by_cancer.tsv`

Use for the main cancer-specific prioritization table. Report observed RSES-Onco,
evidence coverage and coverage-adjusted score together. Do not report an adjusted
score without its coverage.

### `Table_main_candidates_with_empirical_evidence.tsv`

Separates candidates with at least one direct DepMap or TCGA component from
literature-only or unevaluable hypotheses.

### `Table_main_significant_synthetic_lethal_dependencies.tsv`

Contains DepMap contrasts with more negative target gene effect in the loss group
and Benjamini-Hochberg adjusted `q < 0.05`. This table is empty when no contrast
passes the predefined criterion; an empty result must be reported rather than
replaced by nominal P values.

### `Table_main_significant_expression_compensation.tsv`

Contains target-expression induction in loss models with adjusted `q < 0.05`.
Expression induction is a compensatory signature, not causal proof.

## Supplementary tables

- `Table_S_DepMap_dependency_contrasts.tsv`: every eligible CRISPR contrast.
- `Table_S_DepMap_expression_contrasts.tsv`: every eligible expression contrast.
- `Table_S_TCGA_homozygous_deletion_events.tsv`: event counts and frequencies.
- `Table_S_DepMap_vs_TCGA_DepMap_score_comparison.tsv`: change after adding TCGA.
- `Table_S_complex_biomarkers_not_in_simple_CN_analysis.tsv`: intentionally skipped
  MSI/MMR, HRD, dual-loss and other composite biomarkers.
- `gdc_matrix_qc.tsv`: matrix dimensions, allowed values and duplicate checks.
- `tcga_gene_event_summary.tsv`: gene-level event availability and frequencies.

## Empirical figures

### Figure empirical 1 — cancer-specific ranking

Horizontal ranking of coverage-adjusted RSES-Onco scores for colorectal, gastric
and lung cancer. The legend should state that the score integrates empirical
cohort components with literature-based relation and validation priors.

### Figure empirical 2 — DepMap dependency contrasts

Displays `delta_effect = loss - intact`. Negative values indicate greater target
dependency in the loss cohort. The caption must report the one-sided Mann-Whitney
test and Benjamini-Hochberg correction.

### Figure empirical 3 — expression compensation

Displays target-expression differences between loss and intact cohorts. Positive
values indicate increased target expression in the loss cohort.

### Figure empirical 4 — tumor-event prevalence and selectivity

Integrates TCGA homozygous-deletion frequency with the DepMap selectivity component.
This is a prioritization visualization and not a causal model.

### Figure empirical 5 — cross-cancer heatmap

Compares coverage-adjusted scores across colorectal, gastric and lung contexts.
Missing cells must remain visually distinguishable from measured zeros.

## Required manuscript updates after execution

1. Replace pilot-only counts with release-specific sample and model counts.
2. State DepMap release, exact file names and acquisition date.
3. State the GDC query, ASCAT3 workflow and manifest checksum.
4. Describe default-model selection for duplicate DepMap profile rows.
5. Report the 0.30 linear copy-number threshold as a strong loss/LoF-like grouping,
   not a universal homozygous-deletion definition.
6. Describe TCGA event coding as total integer copy number equal to zero.
7. Report all tests, alternatives, minimum group size and FDR procedure.
8. Distinguish empirical components from literature priors in RSES-Onco.
9. Report unavailable domains, including the SOD2 limitation, as missing.
10. Include negative and non-significant findings; do not retain only favorable pairs.
11. Add a limitations paragraph on cell-line versus primary-tumor differences.
12. Add experimental validation plans for the highest supported candidates.

## Recommended results-section structure

1. **Release-specific data coverage and quality control**
2. **Cancer-specific frequency of candidate loss events in TCGA**
3. **Loss-associated target dependencies in DepMap**
4. **Expression compensation in loss-defined models**
5. **Integrated RSES-Onco prioritization**
6. **Robustness, missing evidence and complex biomarkers**
7. **Experimental validation roadmap**

## Evidence language

Use:

- “compatible with synthetic-lethal selectivity”;
- “prioritized for experimental validation”;
- “loss-associated dependency”;
- “expression pattern compatible with compensation”.

Avoid without direct validation:

- “proved synthetic lethal”;
- “clinically effective”;
- “tumor-specific target”;
- “causal backup mechanism”.
