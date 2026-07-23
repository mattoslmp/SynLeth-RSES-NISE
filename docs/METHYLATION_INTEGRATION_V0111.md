# RSES-Onco v0.11.1: GDC promoter methylation integration and complete rerun protocol

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

## 1. Purpose

This document describes the methylation extension introduced in RSES-Onco v0.11.1. It defines the data source, acquisition process, promoter aggregation, pair-level evidence, score formula, missing-data rules, publication outputs and exact commands required to reproduce the analysis.

The implementation does not create an independent global methylation domain. Promoter methylation is integrated inside the existing expression-compensation domain so that epigenetic and expression evidence do not receive duplicated full-domain weights.

## 2. Current RSES-Onco domains

The global score uses the following weights:

| Domain | Internal name | Weight |
|---|---|---:|
| Tumor event | `tumor_event` | 0.16 |
| Conditional dependency | `dependency` | 0.22 |
| Selectivity | `selectivity` | 0.14 |
| Expression compensation | `expression_compensation` | 0.08 |
| Functional relation | `functional_relation` | 0.06 |
| Functional microniche | `functional_microniche` | 0.16 |
| Validation and tractability | `validation_tractability` | 0.18 |

The functional-microniche component contains:

| Microniche subdomain | Internal name | Weight within the microniche |
|---|---|---:|
| Expression context | `expression_context` | 0.20 |
| Localization | `localization` | 0.15 |
| Biochemical/structural evidence | `biochemical_structural` | 0.15 |
| Genetic phenotype | `genetic_phenotype` | 0.20 |
| Interaction network | `interaction_network` | 0.15 |
| Regulatory network | `regulatory_network` | 0.15 |

Within expression context, pairwise expression divergence and WGCNA share the existing expression-context weight equally. Within regulatory network, DoRothEA regulator divergence, TF-expression-profile divergence and JASPAR/FIMO promoter-motif divergence share the existing regulatory-network weight.

## 3. Coverage-aware score

For candidate $i$, cancer context $c$ and eligible domain $d$, define:

- $x_{icd}$: normalized observed domain value;
- $w_d$: domain weight;
- $m_{icd}$: availability indicator;
- $E_{ic}$: eligible domains;
- $O_{ic}$: observed eligible domains.

The observed-domain score is:

$$
S_{ic}=\frac{\sum_{d\in O_{ic}}w_dx_{icd}}{\sum_{d\in O_{ic}}w_d}
$$

Coverage is:

$$
C_{ic}=\frac{\sum_{d\in O_{ic}}w_d}{\sum_{d\in E_{ic}}w_d}
$$

The coverage-adjusted score is:

$$
A_{ic}=S_{ic}C_{ic}=\frac{\sum_{d\in O_{ic}}w_dx_{icd}}{\sum_{d\in E_{ic}}w_d}
$$

Missing eligible evidence is omitted from the observed-score denominator but lowers coverage. Non-eligible evidence is excluded from both denominators. An observed numeric zero remains a real observed zero.

## 4. Methylation source

The methylation layer uses open primary-tumor **Methylation Beta Value** files from the NCI Genomic Data Commons for:

- TCGA-COAD;
- TCGA-READ;
- TCGA-STAD;
- TCGA-LUAD;
- TCGA-LUSC.

The GDC harmonized methylation workflow uses SeSAMe and provides beta-values for HM27, HM450 and EPIC methylation arrays. The pipeline also downloads the official GDC GENCODE-v36 array annotation manifests.

Official annotation resources used by the downloader:

| Platform | GDC UUID | File |
|---|---|---|
| EPIC | `5ce8ae8f-3386-4d12-9035-152742aa07e0` | `EPIC.hg38.manifest.gencode.v36.tsv.gz` |
| HM27 | `e5182c42-bdc6-433e-9b4a-7b7c6696ce89` | `HM27.hg38.manifest.gencode.v36.tsv.gz` |
| HM450 | `021a2330-951d-474f-af24-1acd77e7664f` | `HM450.hg38.manifest.gencode.v36.tsv.gz` |

## 5. Promoter definition and aggregation

The aggregation script retains probes that meet at least one promoter-proximity criterion:

1. annotated distance from the transcription start site between -2000 and +500 bp; or
2. an annotation group indicating TSS1500, TSS200, 5' UTR, first exon or promoter.

Probes marked by the official general mask are excluded when the annotation contains `MASK_general`.

For each sample and gene, the promoter beta-value is the median of all retained promoter probes mapped to that gene:

$$
\beta_{g,s}=\operatorname{median}_{p\in P_g}(\beta_{p,s})
$$

where $P_g$ is the set of retained promoter probes for gene $g$.

## 6. Directional pair-level methylation evidence

For a directional hypothesis `lost_gene -> target_gene`, the cancer-specific summary includes:

- paired sample count;
- lost-gene median promoter beta;
- target-gene median promoter beta;
- median paired delta beta;
- Spearman correlation between promoter beta profiles;
- one-sided paired Wilcoxon test of lost-gene beta greater than target-gene beta;
- global and within-cancer Benjamini-Hochberg FDR;
- explicit missingness or non-eligibility reason.

The directional methylation context score is:

$$
M_{ic}=\widetilde{\beta}_{lost,ic}\left(1-\widetilde{\beta}_{target,ic}\right)
$$

where the medians are clipped to [0,1]. The score is high only when the hypothesized lost gene is more promoter-methylated and the target/paralog remains less promoter-methylated.

This score is not proof that the lost gene is silenced. Methylation is treated as epigenetic context and must be interpreted together with expression and event data.

## 7. Integration into expression compensation

The global expression-compensation weight remains 0.08. Its internal composition becomes:

| Subcomponent | Internal weight |
|---|---:|
| Event-stratified expression compensation | 0.70 |
| Promoter methylation context | 0.30 |

The observed internal score is:

$$
E^{obs}_{ic}=\frac{0.70X_{ic}m_X+0.30M_{ic}m_M}{0.70m_X+0.30m_M}
$$

Internal coverage is:

$$
C^{expr+meth}_{ic}=\frac{0.70m_X+0.30m_M}{1.00}
$$

The expression-compensation value supplied to the global RSES-Onco score is:

$$
E^{adj}_{ic}=E^{obs}_{ic}C^{expr+meth}_{ic}
$$

Thus:

- when both expression and methylation are observed, both contribute with 70/30 weights;
- when only expression is observed, the raw internal value remains the expression value, but internal coverage is 0.70;
- missing methylation is never converted to zero;
- methylation does not receive an additional global weight;
- TCGA methylation is not added to the DepMap-only ranking.

## 8. Pipeline modes

`METHYLATION_MODE` accepts:

- `auto`: use an existing promoter-methylation matrix; otherwise continue without methylation;
- `download`: query, download, validate, aggregate and integrate GDC methylation;
- `require`: require a usable methylation layer and fail otherwise;
- `off`: disable methylation explicitly.

## 9. Exact acquisition commands

### Create the GDC methylation manifest

```bash
python -u scripts/download_gdc_methylation.py \
  --stage manifest \
  --output-dir data/raw/methylation
```

### Download and validate all open methylation files

```bash
python -u scripts/download_gdc_methylation.py \
  --stage all \
  --output-dir data/raw/methylation \
  --workers 4 \
  --retries 4
```

### Aggregate promoter beta-values

```bash
python -u scripts/aggregate_gdc_methylation.py \
  --manifest data/raw/methylation/gdc_methylation_manifest.tsv \
  --annotation-dir data/raw/methylation/annotations \
  --output data/processed/methylation/gdc_promoter_methylation_gene_sample.tsv \
  --gene-summary data/processed/methylation/gdc_promoter_methylation_gene_summary.tsv \
  --status-output data/processed/methylation/gdc_promoter_methylation_aggregation_status.tsv
```

### Build candidate-pair methylation evidence

```bash
python -u scripts/build_methylation_pair_evidence.py \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --gene-sample data/processed/methylation/gdc_promoter_methylation_gene_sample.tsv \
  --output data/processed/methylation/pair_promoter_methylation_evidence.tsv \
  --min-samples 10
```

## 10. Complete integrated rerun

Run this after any previous runner has finished and after updating the local repository to v0.11.1:

```bash
NEW="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010"
cd "$NEW" || exit 1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco
python -m pip install -e .

METHYLATION_MODE=download \
MPLBACKEND=Agg \
STRICT_LAYOUT=1 \
PYTHONUNBUFFERED=1 \
bash scripts/resume_wgcna_regulatory_pipeline.sh resume-regulatory \
  2>&1 | tee logs/run_rses_v0111_methylation.log

status=${PIPESTATUS[0]}
printf '%s\n' "$status" > logs/run_rses_v0111_methylation.exitcode
test "$status" -eq 0
```

For an existing aggregated matrix, use `METHYLATION_MODE=auto`. To run only methylation acquisition/aggregation/evidence construction, use:

```bash
METHYLATION_MODE=download \
bash scripts/resume_wgcna_regulatory_pipeline.sh methylation-only
```

## 11. Validation

```bash
python -u scripts/validate_methylation_evidence.py \
  --evidence data/processed/methylation/pair_promoter_methylation_evidence.tsv \
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv
```

The integrated ranking must contain:

```text
score_version=RSES-Onco-expanded-v0.11.1
expression_regulatory_semantics_version=eligibility-aware-wgcna-regulatory-methylation-v4
```

The DepMap-only ranking must retain its non-methylation version and must not contain methylation-derived score columns.

## 12. Publication outputs

The v0.11.1 publication package retains 77 figures and adds methylation to Supplementary Figure S69. It contains 51 registered tables, including:

- Table S45: pair-level promoter methylation evidence;
- Table S46: gene/cancer promoter methylation summary;
- Table S47: methylation source and aggregation status.

When methylation was not acquired, these tables and Figure S69 report explicit unavailability instead of biological zero.

## 13. References

1. National Cancer Institute. Genomic Data Commons: DNA methylation array harmonization and Methylation Beta Value data documentation.
2. National Cancer Institute. GDC Reference Files: GENCODE-v36 HM27, HM450 and EPIC methylation probe annotation manifests.
3. Zhou W, Triche TJ Jr, Laird PW, Shen H. SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. *Nucleic Acids Research*. 2018;46:e123.
4. Bibikova M, Barnes B, Tsan C, et al. High density DNA methylation array with single CpG site resolution. *Genomics*. 2011;98:288-295.
5. Moran S, Arribas C, Esteller M. Validation of a DNA methylation microarray for 850,000 CpG sites of the human genome enriched in enhancer sequences. *Epigenomics*. 2016;8:389-399.
6. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*. 1995;57:289-300.
