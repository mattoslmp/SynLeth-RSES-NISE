# RSES-Onco v0.10.9: WGCNA correlation-policy correction

## Correction

The primary WGCNA correlation remains biweight midcorrelation (`bicor`). The implementation now applies one explicit policy consistently to soft-threshold selection, signed adjacency and signed module membership (`kME`):

```text
primary correlation: bicor
maxPOutliers: 0.10
pearsonFallback: individual
```

Pearson is used only for an individual gene or module eigengene with zero or non-finite median absolute deviation. It is not applied globally. Every fallback-eligible entity is exported with its MAD and reason.

The correction also replaces unsafe named-vector indexing with `match`, so candidate genes removed during WGCNA quality control are recorded as `gene_missing_after_wgcna_qc` instead of terminating the run. Per-cancer WGCNA outputs are resumable only when all required outputs, including the fallback audit, are complete.

## New audit outputs

- `data/processed/regulatory/wgcna/wgcna_correlation_fallback_all_cancers.tsv`
- `data/processed/regulatory/wgcna/wgcna_run_diagnostics_all_cancers.tsv`
- per-cancer `wgcna_correlation_fallback.tsv`

## Required rerun

Because the correlation policy and kME settings changed, remove the generated WGCNA outputs for `colon`, `stomach` and `lung`, then rerun the regulatory workflow and all publication assets. Cached Ensembl promoters, JASPAR motifs, FIMO results, DepMap inputs, functional-source caches, pharmacology and structures may be reused.

The recalculated ranking must contain:

```text
expression_regulatory_semantics_version=eligibility-aware-wgcna-regulatory-v3
score_version=RSES-Onco-expanded-v0.10.9
```
