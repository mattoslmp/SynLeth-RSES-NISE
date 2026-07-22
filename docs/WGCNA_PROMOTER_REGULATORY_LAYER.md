# WGCNA, promoter and transcription-factor layer for RSES-Onco v0.10.8

## Scientific design

This implementation keeps the original functional-microniche domain weights. It does not add WGCNA, pairwise coexpression, DoRothEA and promoter motifs as four independent full-weight domains.

### Expression-context domain

The existing expression-context domain is internally decomposed into:

- cancer-specific pairwise expression divergence: 0.50;
- consensus signed WGCNA context: 0.50.

WGCNA uses DepMap log2(TPM+1), cancer-compatible models, all mapped candidate genes plus cancer-specific highly variable genes, biweight midcorrelation, signed adjacency, signed TOM, dynamic tree cut, merged modules, kME and intramodular connectivity. The pair subscore combines TOM divergence (0.40), module divergence (0.30) and kME divergence (0.30). Cancer-specific results remain available, but the functional prior uses their median so cancers are not counted as independent evidence units.

### Regulatory-network domain

The existing regulatory-network domain is internally decomposed into:

- DoRothEA regulator-set divergence: 0.40;
- cancer-specific TF-target expression-profile divergence: 0.35;
- JASPAR promoter motif divergence: 0.25.

Canonical promoter windows are based on Ensembl canonical transcript TSS coordinates, 2 kb upstream and 500 bp downstream. JASPAR 2026 CORE vertebrate non-redundant motifs are scanned with FIMO.

A JASPAR motif is predicted cis-regulatory support. It is not direct TF binding, ChIP evidence, promoter occupancy, causal regulation or experimental validation. The pipeline sets `direct_promoter_binding_claim = false` and rejects outputs that violate this rule.

### Missingness and overlap

Missing subcomponents are not zero-imputed. They lower internal coverage. The total parent-domain weights remain unchanged, preventing expression-derived and regulatory-derived evidence from being counted multiple times.

## New source tables

- `data/processed/regulatory/wgcna/wgcna_pair_metrics_all_cancers.tsv`
- `data/processed/regulatory/wgcna/wgcna_input_preparation.tsv`
- `data/processed/regulatory/promoter_tf_regulatory_pair_metrics.tsv`
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
  'r-wgcna>=1.74' \
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

The command preserves the pre-WGCNA functional evidence, downloads/caches promoter and JASPAR data, runs WGCNA and FIMO, recalculates the DepMap-only and integrated TCGA/DepMap scores, regenerates pharmacology and publication assets, performs sublayer ablation analyses, validates the package, runs the complete test suite and rebuilds checksums.

## Required validation

The final ranking retains `scoring_semantics_version=eligibility-aware-v1` for compatibility with the eligibility audit and adds:

```text
expression_regulatory_semantics_version=eligibility-aware-wgcna-regulatory-v2
score_version=RSES-Onco-expanded-v0.10.8
```

The WGCNA/regulatory validator requires all values to be bounded, all promoter evidence to be labelled as motif prediction, and all direct promoter-binding claims to remain false unless a separate traceable direct-binding source is explicitly added in a future version.
