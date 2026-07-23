# RSES-Onco v0.11.1: promoter methylation data and score integration

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

## 1. Scientific role

Promoter methylation is incorporated as epigenetic context within the existing
regulatory-network domain of the functional microniche. It does not receive a new
independent top-level RSES-Onco weight. This design limits double counting because
methylation, promoter motifs, transcription-factor associations and TF-expression
consistency all describe related aspects of gene regulation.

The methylation layer addresses two distinct questions:

1. do the two genes have divergent promoter-methylation profiles within the same
   cancer lineage;
2. is loss of the origin gene associated with lower promoter methylation of the
   proposed target gene, consistent with possible epigenetic derepression?

The second question is an association test. It is not proof that methylation causes
silencing, derepression, compensation or synthetic lethality.

## 2. Supported source

The primary supported source is the DepMap custom-download dataset labelled
`Methylation (1kb upstream TSS)`. A traceable historical CCLE RRBS matrix is also
supported for reproducibility:

```text
CCLE_RRBS_TSS1kb_20181022.txt.gz
```

The historical CCLE measurement represents promoter/TSS methylation ratios for CpG
sites within 1 kb upstream of transcription start sites. Values are expected to be
beta-like methylation ratios in the interval [0,1].

TCGA/GDC methylation beta-value data may be used in a future external-validation
layer. They are not silently substituted for cell-line methylation in the current
score because platform, sample type and cohort semantics differ.

## 3. Accepted local filenames

The resume pipeline searches for the following names under `DEPMAP_DIR`:

```text
Methylation_(1kb_upstream_TSS)_subsetted_NAsdropped.csv
Methylation_1kb_upstream_TSS.csv
CCLE_RRBS_TSS1kb_20181022.txt.gz
CCLE_RRBS_TSS1kb_20181022.txt
```

A different path can be supplied explicitly:

```bash
export METHYLATION="/absolute/path/to/promoter_methylation.csv"
```

## 4. Supported table layouts

`scripts/integrate_methylation_regulatory_layer.py` and
`src/rses_onco/methylation.py` support:

- ModelID-by-promoter-feature matrices;
- long ModelID/gene/methylation tables;
- historical promoter-feature-by-cell-line matrices.

Cell-line names are mapped to `ModelID` through `Model.csv`. Promoter-feature labels
are mapped to gene symbols. If more than one 1-kb-upstream-TSS feature maps to the
same gene, the per-model median is used and the number of original promoter
features is retained.

## 5. Pair-profile divergence

For lost gene `L`, target gene `T` and cancer context `c`, complete methylation
profiles are selected among compatible DepMap models. Let `rho` be the Spearman
correlation between promoter methylation of `L` and `T`, and let `D_beta` be the
median absolute beta-value difference.

```text
correlation_divergence = (1 - rho) / 2
methylation_separation = clip(D_beta / 0.25, 0, 1)
pair_profile_divergence = mean(available terms)
```

The saturation value 0.25 is explicit and configurable through
`--difference-saturation` or `METHYLATION_DIFFERENCE_SATURATION`.

## 6. Conditional target hypomethylation

Models are stratified using the same copy-number loss threshold used in the
conditional dependency and expression analyses:

```text
loss group:   copy number of L < 0.30
intact group: copy number of L >= 0.30
```

For target promoter methylation:

```text
delta_beta = median_beta_loss - median_beta_intact
hypomethylation_support = clip((-delta_beta) / 0.25, 0, 1)
```

A positive support value therefore requires lower target-promoter methylation in
the loss group. The loss and intact groups must each meet the declared minimum
sample size. The Mann–Whitney U test is used for the group comparison, and
Benjamini–Hochberg values are exported globally and within each cancer context.

## 7. Methylation subscore

The methylation context is coverage-aware:

```text
0.50 × pair promoter-methylation profile divergence
0.50 × conditional target-promoter hypomethylation support
```

Missing subcomponents are omitted from the observed numerator and denominator and
lower methylation subcoverage. A measured value of zero remains a real observed
zero.

## 8. Regulatory-network subscore

The existing regulatory-network domain weight remains unchanged. Its internal
composition becomes:

```text
0.32 × DoRothEA TF-association divergence
0.28 × TF-expression-profile divergence
0.20 × JASPAR/FIMO promoter-motif divergence
0.20 × promoter-methylation context
```

The previous 0.40/0.35/0.25 regulatory proportions are rescaled to 80% of the
internal domain and methylation receives the remaining 20%. The entire regulatory
network continues to occupy 0.15 of the functional microniche, which itself
occupies 0.16 of the top-level RSES-Onco score.

When every top-level and microniche domain is eligible and observed, the maximum
direct contribution represented by the methylation subcomponent is:

```text
0.20 × 0.15 × 0.16 = 0.0048
```

That is 0.48% of the fully observed top-level weighted score before the effects of
nested coverage adjustment. The intentionally small contribution prevents an
epigenetic association from dominating dependency or selectivity evidence.

## 9. Missing-data states

If the source is absent, the pipeline still creates an auditable methylation table
with:

```text
component_promoter_methylation_context = NA
methylation_coverage = 0
methylation_absence_reason = methylation_source_not_provided_or_file_absent
```

This is eligible missing evidence, not negative biological evidence. Other absence
reasons include unmapped genes, unmapped cell lines, insufficient complete profiles
and insufficient loss/intact group sizes.

## 10. Outputs

```text
data/processed/regulatory/promoter_methylation_pair_metrics.tsv
data/processed/regulatory/promoter_methylation_status.json
data/processed/regulatory/expanded_pair_functional_evidence_by_cancer.tsv
article_outputs/tables/supporting_evidence/expression_regulatory/promoter_methylation_pair_metrics.tsv
```

The final rankings include methylation raw values, coverage, group sizes, effect,
P value, FDR, absence reason, formula and semantics version.

## 11. Execution

After the current pre-v0.11.1 run finishes, update the repository and provide the
methylation source before recalculating:

```bash
export DEPMAP_DIR="$PWD/data/raw/depmap"
export METHYLATION="$DEPMAP_DIR/Methylation_(1kb_upstream_TSS)_subsetted_NAsdropped.csv"

bash scripts/resume_wgcna_regulatory_pipeline.sh resume-regulatory
```

The corrected methylation-aware ranking must report:

```text
score_version = RSES-Onco-expanded-v0.10.10
methylation_semantics_version = promoter-methylation-context-v1
regulatory_layer_version = wgcna-promoter-methylation-regulatory-v3
```

## 12. Interpretation boundary

Promoter methylation is epigenetic association evidence. It must not be described
as direct causal silencing, proof of compensatory transcription, proof of promoter
occupancy, therapeutic validation or clinical efficacy without separate evidence.
