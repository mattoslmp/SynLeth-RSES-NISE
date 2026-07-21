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
- **101 cross-cluster human analogue pairs**, evaluated in both directions as **202 loss-to-target hypotheses**.
- Literature benchmark of **25 validated vulnerabilities and explicitly labelled hypotheses**.
- Ensembl Compara expansion of homologous paralogs for every analyzed seed gene.
- DepMap Public 26Q1 readers for CRISPR gene effect, WGS copy number, model metadata and expression.
- Human CRISPR phenotype-profile divergence analogous to the mutant-phenotype domain used in EcoNISE.
- Cancer-specific expression-context divergence analogous to the PRECISE expression domain.
- STRING functional-interaction neighborhoods.
- OmniPath/DoRothEA transcriptional-regulatory neighborhoods.
- Human Protein Atlas subcellular localization.
- UniProt catalytic, cofactor, localization and PDB traceability.
- Open GDC ASCAT3 primary-tumor gene-level copy-number acquisition with resume, `.part` files, retries, size validation and MD5 validation.
- Aggregation of TCGA-COAD/READ, TCGA-STAD and TCGA-LUAD/LUSC into deletion-only matrices.
- Cancer-specific RSES-Onco integration, Benjamini-Hochberg adjusted contrasts, article tables, empirical figures and supplementary workbook generation.
- Tests, GitHub Actions and documented benchmark and expanded shell workflows.

## Scientific boundary

The repository separates three classes of output:

1. literature-anchored priors;
2. synthetic software-verification fixtures;
3. release-specific empirical TCGA/DepMap and human-network results generated locally.

Only the third class may be described as a real cohort analysis. Computational
candidates still require biomarker-confirmed isogenic, pharmacological, rescue and
mechanistic validation before therapeutic or clinical interpretation.

“All NISEs” is exhaustive for the bundled curated human NISE catalogue. “All
paralogs” is source-bounded to the Ensembl release and the analyzed seed genes.
There is no canonical resource that exhaustively defines every possible pathway
backup, downstream dependency or collateral vulnerability in the human proteome;
those classes must retain explicit source and release provenance.

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

## Real-data documentation

- [`docs/REAL_DATA_WORKFLOW.md`](docs/REAL_DATA_WORKFLOW.md)
- [`docs/ARTICLE_ANALYSIS_OUTPUTS.md`](docs/ARTICLE_ANALYSIS_OUTPUTS.md)
- [`docs/EXPANDED_HUMAN_EVIDENCE_WORKFLOW.md`](docs/EXPANDED_HUMAN_EVIDENCE_WORKFLOW.md)

## Benchmark workflow

```bash
bash scripts/run_real_data_pipeline.sh all
```

When the GDC download already exists:

```bash
bash scripts/run_real_data_pipeline.sh after-download
```

## Expanded all-NISE workflow

Run the complete workflow, including all NISE directions, Ensembl paralogs,
STRING, DoRothEA, HPA, UniProt, DepMap and TCGA:

```bash
bash scripts/run_real_data_pipeline.sh expanded-all
```

When the GDC download is already running or complete, do not restart it. Continue
with:

```bash
bash scripts/run_real_data_pipeline.sh expanded-after-download
```

The expanded setup before TCGA can be run independently:

```bash
bash scripts/run_real_data_pipeline.sh expanded-setup
```

Expanded stages:

```bash
bash scripts/run_real_data_pipeline.sh build-universe
bash scripts/run_real_data_pipeline.sh expand-paralogs
bash scripts/run_real_data_pipeline.sh functional-evidence
bash scripts/run_real_data_pipeline.sh run-expanded-depmap
bash scripts/run_real_data_pipeline.sh run-expanded-full
bash scripts/run_real_data_pipeline.sh summarize-expanded
bash scripts/run_real_data_pipeline.sh figures-expanded
bash scripts/run_real_data_pipeline.sh workbook-expanded
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

## Build the expanded candidate universe

```bash
python -u scripts/build_expanded_candidate_universe.py \
  --output data/processed/expanded_candidate_universe.tsv \
  --members-output data/processed/expanded_class_member_inventory.tsv
```

This includes every direction of every curated cross-cluster NISE pair plus all
benchmark classes. Add Ensembl paralogs:

```bash
python -u scripts/download_ensembl_paralogs.py \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --output data/raw/ensembl/human_seed_paralogs.tsv

python -u scripts/build_expanded_candidate_universe.py \
  --additional data/raw/ensembl/human_seed_paralogs.tsv \
  --output data/processed/expanded_candidate_universe.tsv \
  --members-output data/processed/expanded_class_member_inventory.tsv
```

Additional source-bounded classes can be included with repeated `--additional`
arguments.

## Acquire human functional evidence

```bash
python -u scripts/download_human_functional_evidence.py \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --raw-dir data/raw/human_functional_evidence \
  --output data/processed/expanded_pair_functional_evidence.tsv
```

The acquisition stage records:

```text
STRING functional interaction partners
OmniPath/DoRothEA TF-target interactions
Human Protein Atlas localization
UniProt reviewed annotations and PDB cross-references
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

## Two-level scoring

### Human functional-microniche RSES

```text
Expression/context       0.20
Localization             0.15
Biochemical/structural   0.15
Genetic/phenotype        0.20
STRING network           0.15
Regulatory network       0.15
```

### Expanded RSES-Onco

```text
Tumor event              0.16
Conditional dependency   0.22
Loss selectivity         0.14
Expression compensation  0.08
Functional relation      0.06
Functional microniche    0.16
Validation/tractability  0.18
```

For available domains `d`, weights `w_d`, values `D_d` and availability indicators
`m_d`:

```text
Score = sum(w_d m_d D_d) / sum(w_d m_d)
Coverage = sum(w_d m_d) / sum(w_d)
Adjusted score = Score × Coverage
```

## Statistical conventions

- DepMap strong-loss cohort: linear WGS relative copy number `< 0.30`.
- Minimum group size: three loss and three intact models.
- Dependency test: one-sided Mann-Whitney, testing more negative target effect in loss models.
- Expression-compensation test: one-sided Mann-Whitney, testing increased target expression in loss models.
- Expression-context domain: cancer-specific Spearman decorrelation and median absolute abundance separation.
- Genetic/phenotype domain: cancer-specific CRISPR profile decorrelation, score separation and dependency-set non-overlap.
- STRING and regulatory domains: pairwise neighborhood/regulator-set divergence.
- Multiple testing: Benjamini-Hochberg across eligible contrasts.
- TCGA homozygous deletion: ASCAT3 integer total copy number equal to zero.

The DepMap threshold is a strong loss/LoF-like grouping and must not be described
as a universal pure homozygous-deletion definition. The TCGA output is a
deletion-only event matrix, not a GISTIC call set.

## Expanded outputs

```text
results/expanded_26Q1/depmap_only/
results/expanded_26Q1/full/
results/expanded_26Q1/full/article_tables/
figures/expanded_26Q1/                    # PDF, PNG and SVG
supplementary/RSES_Onco_Expanded_All_NISE_26Q1.xlsx
data/processed/expanded_candidate_universe.tsv
data/processed/expanded_class_member_inventory.tsv
data/processed/expanded_pair_functional_evidence.tsv
```

Raw third-party matrices are intentionally excluded from version control.

## Current known limitation

In the evaluated DepMap Public 26Q1 files, `SOD2` is present in CRISPR gene effect
but absent from WGS copy number, protein-coding expression and inferred LoF
annotations. Therefore `SOD2 -> SOD1` remains an exploratory hypothesis with
missing empirical loss-cohort domains. Missingness is preserved and does not imply
a negative biological result.

## Main data resources

- Curated human NISE catalogue and cross-cluster pairs.
- DepMap Public 26Q1.
- TCGA/GDC ASCAT3 gene-level copy number.
- STRING functional associations.
- OmniPath/DoRothEA transcriptional regulation.
- Human Protein Atlas subcellular localization.
- UniProtKB reviewed protein annotations.
- Ensembl Compara human paralogs.

## License

MIT for code. Third-party data remain under their original terms and are acquired
by scripts rather than redistributed.
