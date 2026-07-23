# WGCNA, promoter and transcription-factor layer for RSES-Onco v0.10.9

## Scientific design

This implementation keeps the original functional-microniche domain weights. It does not add WGCNA, pairwise coexpression, DoRothEA and promoter motifs as four independent full-weight domains.

### Expression-context domain

The existing expression-context domain is internally decomposed into:

- cancer-specific pairwise expression divergence: 0.50;
- cancer-specific signed WGCNA context: 0.50.

WGCNA uses DepMap log2(TPM+1), cancer-compatible models, all mapped candidate genes plus cancer-specific highly variable genes, and biweight midcorrelation (`bicor`) as the primary correlation. The same explicit settings are used for soft-threshold selection, signed adjacency and kME: `maxPOutliers=0.10` and `pearsonFallback="individual"`. Pearson is used only for an individual gene or module eigengene with zero or non-finite MAD; every affected entity is exported in a correlation-fallback audit table. Signed TOM, dynamic tree cut, merged modules, kME and intramodular connectivity are then calculated. The pair subscore combines TOM divergence (0.40), module divergence (0.30) and kME divergence (0.30). Each cancer-specific RSES-Onco row receives only the WGCNA network derived from the corresponding cancer models. A separate median consensus table is retained for pair-level source compatibility and descriptive reporting, not as an additional evidence unit.

### Regulatory-network domain

The existing regulatory-network domain is internally decomposed into:

- DoRothEA regulator-set divergence: 0.40;
- cancer-specific TF-target expression-profile divergence: 0.35;
- JASPAR promoter motif divergence: 0.25.

Canonical promoter windows are based on Ensembl canonical-transcript TSS coordinates, 2 kb upstream and 500 bp downstream. JASPAR 2026 CORE vertebrate non-redundant motifs are scanned with FIMO. Each cancer score receives the DoRothEA/TF-expression regulatory context calculated for that cancer; promoter motif annotations are sequence based and therefore shared across cancers without being multiplied into separate evidence units.

A JASPAR motif is predicted cis-regulatory support. It is not direct TF binding, ChIP evidence, promoter occupancy, causal regulation or experimental validation. The pipeline sets `direct_promoter_binding_claim = false` and rejects outputs that violate this rule.

### Missingness, overlap and ablation

Missing subcomponents are not zero-imputed. They lower internal coverage. The total parent-domain weights remain unchanged, preventing expression-derived and regulatory-derived evidence from being counted multiple times.

Ablations remove the selected subcomponent from the eligible internal model and renormalize the retained subcomponents. This differs from treating the ablated source as missing. The workflow evaluates rankings without WGCNA, pairwise expression, DoRothEA regulator sets, TF-expression consistency, promoter motifs and the complete regulatory domain.

## New source tables

- `data/processed/regulatory/wgcna/wgcna_pair_metrics_all_cancers.tsv`
- `data/processed/regulatory/wgcna/wgcna_input_preparation.tsv`
- `data/processed/regulatory/wgcna/wgcna_correlation_fallback_all_cancers.tsv`
- `data/processed/regulatory/wgcna/wgcna_run_diagnostics_all_cancers.tsv`
- `data/processed/regulatory/promoter_tf_regulatory_pair_metrics.tsv`
- `data/processed/regulatory/expanded_pair_functional_evidence_by_cancer.tsv`
- `data/raw/regulatory/ensembl_promoters.tsv`
- `data/raw/regulatory/ensembl_promoters.fa`
- `data/processed/regulatory/jaspar_promoter_motif_hits.tsv`
- `data/processed/regulatory/jaspar_promoter_tf_summary.tsv`
- `article_outputs/tables/robustness/wgcna_regulatory_ablation_scores.tsv`
- `article_outputs/tables/robustness/wgcna_regulatory_ablation_summary.tsv`

## Environment update

```bash
conda activate rses-onco
conda install -c conda-forge -c bioconda \
  'r-base>=4.3' \
  r-wgcna \
  r-dynamictreecut \
  r-fastcluster \
  meme
```

Validate the runtime:

```bash
bash scripts/resume_wgcna_regulatory_pipeline.sh check-runtime
```

## Resume from cached functional evidence

The STRING, DoRothEA, HPA and UniProt caches do not need to be downloaded again.

```bash
bash scripts/resume_wgcna_regulatory_pipeline.sh resume-regulatory
```

The command preserves the pre-WGCNA functional evidence, downloads and caches promoter and JASPAR data, runs WGCNA and FIMO, recalculates the DepMap-only and integrated TCGA/DepMap scores, regenerates pharmacology and publication assets, performs sublayer ablations, validates the package, runs the complete test suite and rebuilds checksums. Cached pharmacology and structure outputs are reused through `assets-only` when present; the full publication acquisition path is selected automatically only when a mandatory cache is absent.

## Required validation

The final ranking retains `scoring_semantics_version=eligibility-aware-v1` for compatibility with the eligibility audit and adds:

```text
expression_regulatory_semantics_version=eligibility-aware-wgcna-regulatory-v3
score_version=RSES-Onco-expanded-v0.10.9
```

The WGCNA/regulatory validator requires all values to be bounded, all promoter evidence to be labelled as motif prediction, and all direct promoter-binding claims to remain false unless a separate traceable direct-binding source is explicitly added in a future version.
