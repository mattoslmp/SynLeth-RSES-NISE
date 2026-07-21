# SynLeth-RSES-Onco / RSES-Onco

**RSES-Onco** is a reproducible framework for prioritizing cancer-selective
synthetic-lethal dependencies created by hidden functional backups, analogous
enzymes, homologous paralogs, pathway redundancy and collateral deletion. The
initial disease scope is colorectal, gastric and lung cancer.

The framework adapts the coverage-aware Role Specialization Evidence Score
(RSES): unavailable evidence is omitted from the observed-domain denominator and
reduces explicit coverage rather than being treated as biological similarity or as
a negative result. Functional/evolutionary relation confidence remains separate
from empirical synthetic-lethality evidence.

## Included

- Curated table of **70 human proteins in 15 bona fide intragenomic analogous-enzyme activities**.
- **101 cross-cluster human analogue pairs**.
- Literature benchmark of **25 validated vulnerabilities and explicitly labelled hypotheses**.
- DepMap Public 26Q1 readers for CRISPR gene effect, WGS copy number, model metadata and expression.
- Current OncoTree cohort definitions for colorectal (`Bowel`/colorectal adenocarcinoma), gastric (`Esophagus/Stomach` restricted to stomach/STAD) and lung.
- Selection of unique default rows for duplicate DepMap profiles.
- Open GDC ASCAT3 primary-tumor gene-level copy-number acquisition with resume, `.part` files, retries, size validation and MD5 validation.
- Aggregation of TCGA-COAD/READ, TCGA-STAD and TCGA-LUAD/LUSC into deletion-only matrices.
- Cancer-specific RSES-Onco integration, Benjamini-Hochberg adjusted contrasts, article tables, empirical figures and supplementary workbook generation.
- Tests, GitHub Actions and a documented end-to-end shell orchestrator.

## Scientific boundary

The repository separates three classes of output:

1. literature-anchored priors;
2. synthetic software-verification fixtures;
3. release-specific empirical TCGA/DepMap results generated locally.

Only the third class may be described as a real cohort analysis. Computational
candidates still require biomarker-confirmed isogenic, pharmacological, rescue and
mechanistic validation before therapeutic or clinical interpretation.

## Installation

```bash
conda env create -f environment.yml
conda activate rses-onco
python -m pip install -e .
python -m pytest -q -p no:cacheprovider
```

## Reproduce the bundled literature pilot

```bash
rses-onco score-literature \
  --input data/curated/synthetic_lethality_reference_pairs.tsv \
  --output results/literature_anchored_candidates.tsv
```

## Real-data workflow

The complete documented workflow is in:

- [`docs/REAL_DATA_WORKFLOW.md`](docs/REAL_DATA_WORKFLOW.md)
- [`docs/ARTICLE_ANALYSIS_OUTPUTS.md`](docs/ARTICLE_ANALYSIS_OUTPUTS.md)

The main orchestrator is:

```bash
bash scripts/run_real_data_pipeline.sh all
```

When the GDC download has already started or finished, continue without restarting
it using:

```bash
bash scripts/run_real_data_pipeline.sh after-download
```

Available stages:

```bash
bash scripts/run_real_data_pipeline.sh validate-depmap
bash scripts/run_real_data_pipeline.sh run-depmap
bash scripts/run_real_data_pipeline.sh manifest-gdc
bash scripts/run_real_data_pipeline.sh download-gdc
bash scripts/run_real_data_pipeline.sh validate-gdc
bash scripts/run_real_data_pipeline.sh aggregate-gdc
bash scripts/run_real_data_pipeline.sh validate-matrices
bash scripts/run_real_data_pipeline.sh run-full
bash scripts/run_real_data_pipeline.sh summarize
bash scripts/run_real_data_pipeline.sh figures
bash scripts/run_real_data_pipeline.sh workbook
bash scripts/run_real_data_pipeline.sh finalize
```

## Expected DepMap files

Place these under `data/raw/depmap/`:

```text
CRISPRGeneEffect.csv
OmicsCNGeneWGS.csv
Model.csv
OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
OmicsGlobalSignatures.csv                         # recommended
OmicsInferredMolecularSubtypes.csv                # recommended
README_26Q1.txt                                   # recommended
```

Validate them with:

```bash
python -u scripts/download_depmap.py \
  --input-dir data/raw/depmap \
  --write-checksums

python -u scripts/validate_real_inputs.py \
  --gene-effect data/raw/depmap/CRISPRGeneEffect.csv \
  --copy-number data/raw/depmap/OmicsCNGeneWGS.csv \
  --models data/raw/depmap/Model.csv \
  --expression data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
```

## GDC acquisition

Create and review the manifest:

```bash
python -u scripts/download_gdc.py \
  --manifest-only \
  --workflow ASCAT3 \
  --output-dir data/raw/gdc
```

Download or resume from the reviewed manifest:

```bash
python -u scripts/download_gdc.py \
  --use-existing-manifest \
  --manifest data/raw/gdc/gdc_gene_level_copy_number_manifest.json \
  --output-dir data/raw/gdc \
  --retries 3
```

Validate every file:

```bash
python -u scripts/download_gdc.py \
  --validate-only \
  --manifest data/raw/gdc/gdc_gene_level_copy_number_manifest.json \
  --output-dir data/raw/gdc
```

## Empirical outputs

The integrated workflow creates:

```text
results/empirical_26Q1/depmap_only/
results/empirical_26Q1/full/
results/empirical_26Q1/full/article_tables/
figures/empirical_26Q1/                    # PDF, PNG and SVG
supplementary/RSES_Onco_Empirical_26Q1.xlsx
logs/empirical_26Q1/
```

Raw third-party matrices are intentionally excluded from version control.

## RSES-Onco score

For available domains `d`, weights `w_d`, values `D_d` and availability indicators
`m_d`:

```text
RSES-Onco = sum(w_d m_d D_d) / sum(w_d m_d)
Coverage = sum(w_d m_d) / sum(w_d)
Adjusted score = RSES-Onco × Coverage
```

Default domains are tumor event, target dependency, loss selectivity, expression
compensation, functional-relation confidence and validation/tractability.

## Statistical conventions

- DepMap strong-loss cohort: linear WGS relative copy number `< 0.30`.
- Minimum group size: three loss and three intact models.
- Dependency test: one-sided Mann-Whitney, testing more negative target effect in loss models.
- Expression test: one-sided Mann-Whitney, testing increased target expression in loss models.
- Multiple testing: Benjamini-Hochberg across eligible contrasts.
- TCGA homozygous deletion: ASCAT3 integer total copy number equal to zero.

The DepMap threshold is a strong loss/LoF-like grouping and must not be described
as a universal pure homozygous-deletion definition. The TCGA output is a
deletion-only event matrix, not a GISTIC call set.

## Current known limitation

In the evaluated DepMap Public 26Q1 files, `SOD2` is present in CRISPR gene effect
but absent from WGS copy number, protein-coding expression and inferred LoF
annotations. Therefore `SOD2 -> SOD1` remains an exploratory hypothesis with
missing empirical loss-cohort domains. Missingness is preserved and does not imply
a negative biological result.

## Main references

- Piergiorge RM et al. *Genome Biology and Evolution*. 2017. DOI: 10.1093/gbe/evx119.
- de Oliveira FC et al. *BMC Research Notes*. 2026. DOI: 10.1186/s13104-026-07742-5.
- Kryukov GV et al. *Science*. 2016. DOI: 10.1126/science.aad5214.
- Chan EM et al. *Nature*. 2019. DOI: 10.1038/s41586-019-1102-x.

## License

MIT for code. Third-party data remain under their original terms and are acquired
by scripts rather than redistributed.
