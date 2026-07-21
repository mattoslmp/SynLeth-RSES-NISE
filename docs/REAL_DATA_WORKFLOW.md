# Complete real-data workflow: DepMap Public 26Q1 and TCGA/GDC ASCAT3

This document records the complete command-line workflow used to acquire, validate,
process and integrate real DepMap and TCGA/GDC data for RSES-Onco. It is designed
for Linux or WSL and assumes the repository root is the current directory.

## Scientific scope

The initial cancer cohorts are:

- colorectal: DepMap colorectal models and TCGA-COAD plus TCGA-READ;
- gastric: DepMap stomach/gastric models and TCGA-STAD;
- lung: DepMap lung models and TCGA-LUAD plus TCGA-LUSC.

The empirical workflow analyzes simple single-gene loss biomarkers. Composite
biomarkers such as MSI/MMR deficiency, HRD, replication-gap phenotype, high
alkylation state and dual-gene loss require explicit annotations and are not
silently approximated by a single copy-number column.

## 1. Installation

```bash
cd /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE
conda env create -f environment.yml
conda activate rses-onco
python -m pip install -e .
python -m pytest -q -p no:cacheprovider
```

If the environment already exists:

```bash
conda activate rses-onco
python -m pip install -e .
```

## 2. Required DepMap Public 26Q1 files

Place the following files in `data/raw/depmap/`:

```text
CRISPRGeneEffect.csv
OmicsCNGeneWGS.csv
Model.csv
OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
OmicsGlobalSignatures.csv                         # recommended
OmicsInferredMolecularSubtypes.csv                # recommended
README_26Q1.txt                                   # recommended
```

The main analysis uses the non-stranded protein-coding expression matrix. Do not
mix stranded and non-stranded expression matrices within the same execution.

## 3. Validate DepMap files and checksums

```bash
mkdir -p logs/empirical_26Q1

PYTHONUNBUFFERED=1 \
python -u scripts/download_depmap.py \
  --input-dir data/raw/depmap \
  --write-checksums \
  2>&1 | tee logs/empirical_26Q1/01_validate_depmap_files.log
```

Check the generated hashes:

```bash
(
  cd data/raw/depmap
  sha256sum -c SHA256SUMS.txt
)
```

The validator recognizes a literal `ModelID` column and unnamed first columns
containing `ACH-*` identifiers. In `OmicsCNGeneWGS.csv`, duplicate ModelIDs are
accepted only when exactly one row is marked as the default model entry.

## 4. Validate model crosswalk and cancer cohorts

```bash
PYTHONUNBUFFERED=1 \
python -u scripts/validate_real_inputs.py \
  --gene-effect data/raw/depmap/CRISPRGeneEffect.csv \
  --copy-number data/raw/depmap/OmicsCNGeneWGS.csv \
  --models data/raw/depmap/Model.csv \
  --expression data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
  2>&1 | tee logs/empirical_26Q1/02_validate_depmap_crosswalk.log
```

A validated 26Q1 run used:

```text
Gene effect models: 1,208
Copy-number models: 1,118 unique ModelIDs after default-entry selection
Metadata models: 2,154
Common ModelIDs: 858
Colorectal: 41 common models
Gastric: 26 common models
Lung: 100 common models
```

These exact counts are release-specific. Future releases may differ.

### SOD2 limitation

In the analyzed 26Q1 files, SOD2 was present in `CRISPRGeneEffect.csv` but absent
from WGS copy number, both protein-coding expression matrices and inferred LoF
annotations. Therefore `SOD2 loss/low expression -> SOD1` is retained as an
exploratory hypothesis with missing empirical loss-cohort components. Absence is
not converted to zero or interpreted as a negative result.

## 5. Run the DepMap-only empirical analysis

```bash
mkdir -p results/empirical_26Q1/depmap_only

set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/run_empirical_rses_onco.py \
  --gene-effect data/raw/depmap/CRISPRGeneEffect.csv \
  --copy-number data/raw/depmap/OmicsCNGeneWGS.csv \
  --models data/raw/depmap/Model.csv \
  --expression data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
  --loss-threshold 0.30 \
  --min-group-size 3 \
  --output results/empirical_26Q1/depmap_only/empirical_rses_onco_by_cancer.tsv \
  2>&1 | tee logs/empirical_26Q1/03_run_depmap_only.log

status=${PIPESTATUS[0]}
echo "DepMap-only exit code: $status"
test "$status" -eq 0
```

The threshold `< 0.30` is applied to the current linear WGS relative copy-number
matrix. It defines a strong loss/LoF-like cohort; it is not described as a pure
homozygous-deletion call.

Expected outputs:

```text
results/empirical_26Q1/depmap_only/
├── empirical_rses_onco_by_cancer.tsv
├── dependency_contrasts.tsv
├── expression_contrasts.tsv
└── skipped_complex_biomarkers.tsv
```

## 6. Query the GDC manifest

```bash
mkdir -p data/raw/gdc data/processed

set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/download_gdc.py \
  --manifest-only \
  --workflow ASCAT3 \
  --output-dir data/raw/gdc \
  2>&1 | tee logs/empirical_26Q1/04_gdc_manifest.log

status=${PIPESTATUS[0]}
echo "GDC manifest exit code: $status"
test "$status" -eq 0
```

The reviewed 21 July 2026 manifest contained:

```text
TCGA-COAD: 422 files
TCGA-READ: 153 files
TCGA-STAD: 429 files
TCGA-LUAD: 503 files
TCGA-LUSC: 490 files
Total: 1,997 files; approximately 6.410 GiB
```

Preserve a dated copy:

```bash
cp -av \
  data/raw/gdc/gdc_gene_level_copy_number_manifest.json \
  data/raw/gdc/gdc_gene_level_copy_number_manifest_$(date +%Y%m%d).json
```

## 7. Download or resume GDC files

Use the reviewed existing manifest so that a long download does not silently
change its cohort during a later API query:

```bash
set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/download_gdc.py \
  --use-existing-manifest \
  --manifest data/raw/gdc/gdc_gene_level_copy_number_manifest.json \
  --workflow ASCAT3 \
  --output-dir data/raw/gdc \
  --retries 3 \
  2>&1 | tee logs/empirical_26Q1/05_gdc_download.log

status=${PIPESTATUS[0]}
echo "GDC download exit code: $status"
test "$status" -eq 0
```

The downloader:

- prints progress without buffering;
- writes temporary `.part` files;
- retries failed requests;
- validates size and MD5 before accepting a file;
- skips existing valid files, allowing safe restart.

Monitor from a second terminal:

```bash
watch -n 15 '
echo -n "Complete files: "
find data/raw/gdc/TCGA-* -type f ! -name "*.part" 2>/dev/null | wc -l
echo -n "Partial files:  "
find data/raw/gdc/TCGA-* -type f -name "*.part" 2>/dev/null | wc -l
echo -n "Current volume:  "
du -sh data/raw/gdc 2>/dev/null
'
```

Do not start a second downloader while one is already running.

## 8. Validate all downloaded GDC files

```bash
set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/download_gdc.py \
  --validate-only \
  --manifest data/raw/gdc/gdc_gene_level_copy_number_manifest.json \
  --output-dir data/raw/gdc \
  2>&1 | tee logs/empirical_26Q1/06_gdc_validate.log

status=${PIPESTATUS[0]}
echo "GDC validation exit code: $status"
test "$status" -eq 0
```

Verify that no partial file remains:

```bash
find data/raw/gdc -type f -name '*.part'
```

The command should print nothing. The final validation message must state that all
selected files passed size and MD5 validation.

## 9. Aggregate ASCAT3 files

```bash
set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/aggregate_gdc_gene_cna.py \
  --raw-dir data/raw/gdc \
  --manifest data/raw/gdc/gdc_gene_level_copy_number_manifest.json \
  --output-dir data/processed \
  2>&1 | tee logs/empirical_26Q1/07_gdc_aggregate.log

status=${PIPESTATUS[0]}
echo "GDC aggregation exit code: $status"
test "$status" -eq 0
```

The aggregation produces project-specific and combined cancer matrices. Integer
total copy number is converted as follows:

```text
copy_number == 0  -> -2
copy_number > 0   ->  0
missing           -> NA
```

These files are deletion-only matrices compatible with the RSES-Onco event
reader. They are not GISTIC call sets.

Expected combined files:

```text
data/processed/TCGA_COLON_homdel_discrete.tsv
data/processed/TCGA_STOMACH_homdel_discrete.tsv
data/processed/TCGA_LUNG_homdel_discrete.tsv
```

## 10. Validate aggregated matrices

```bash
mkdir -p results/empirical_26Q1/full

set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/validate_gdc_matrices.py \
  --output results/empirical_26Q1/full/gdc_matrix_qc.tsv \
  --event-output results/empirical_26Q1/full/tcga_gene_event_summary.tsv \
  2>&1 | tee logs/empirical_26Q1/08_validate_gdc_matrices.log

status=${PIPESTATUS[0]}
echo "GDC matrix validation exit code: $status"
test "$status" -eq 0
```

The validator checks:

- allowed values `-2`, `0` and missing;
- duplicated genes and samples;
- expected sample counts from the reviewed manifest;
- candidate-gene availability;
- homozygous-deletion counts and frequencies.

## 11. Run integrated TCGA plus DepMap scoring

```bash
set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/run_empirical_rses_onco.py \
  --gene-effect data/raw/depmap/CRISPRGeneEffect.csv \
  --copy-number data/raw/depmap/OmicsCNGeneWGS.csv \
  --models data/raw/depmap/Model.csv \
  --expression data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
  --loss-threshold 0.30 \
  --min-group-size 3 \
  --tcga colon=data/processed/TCGA_COLON_homdel_discrete.tsv \
  --tcga stomach=data/processed/TCGA_STOMACH_homdel_discrete.tsv \
  --tcga lung=data/processed/TCGA_LUNG_homdel_discrete.tsv \
  --output results/empirical_26Q1/full/empirical_rses_onco_by_cancer.tsv \
  2>&1 | tee logs/empirical_26Q1/09_run_full_tcga_depmap.log

status=${PIPESTATUS[0]}
echo "Full RSES-Onco exit code: $status"
test "$status" -eq 0
```

The integrated ranking preserves missing empirical domains and retains only the
functional-relation and validation/tractability literature priors when empirical
data are unavailable.

## 12. Build article tables

```bash
PYTHONUNBUFFERED=1 \
python -u scripts/summarize_empirical_results.py \
  --input results/empirical_26Q1/full/empirical_rses_onco_by_cancer.tsv \
  --depmap-only results/empirical_26Q1/depmap_only/empirical_rses_onco_by_cancer.tsv \
  --output-dir results/empirical_26Q1/full/article_tables \
  2>&1 | tee logs/empirical_26Q1/10_summarize_empirical_results.log
```

This creates:

- top RSES-Onco rankings by cancer;
- candidates with at least one empirical evidence domain;
- all DepMap dependency and expression contrasts;
- FDR-significant supportive contrasts;
- TCGA homozygous-deletion event tables;
- DepMap-only versus integrated-score comparison;
- skipped complex-biomarker table.

## 13. Generate publication figures

```bash
PYTHONUNBUFFERED=1 \
python -u scripts/make_empirical_figures.py \
  --input results/empirical_26Q1/full/empirical_rses_onco_by_cancer.tsv \
  --output-dir figures/empirical_26Q1 \
  2>&1 | tee logs/empirical_26Q1/11_make_empirical_figures.log
```

Every empirical figure is written in PDF, PNG and SVG.

## 14. Build the empirical supplementary workbook

```bash
PYTHONUNBUFFERED=1 \
python -u scripts/build_empirical_workbook.py \
  --input results/empirical_26Q1/full/empirical_rses_onco_by_cancer.tsv \
  --output supplementary/RSES_Onco_Empirical_26Q1.xlsx \
  2>&1 | tee logs/empirical_26Q1/12_build_empirical_workbook.log
```

## 15. Final validation and checksums

```bash
python -m pytest -q -p no:cacheprovider \
  2>&1 | tee logs/empirical_26Q1/13_pytest.log

find \
  results/empirical_26Q1/depmap_only \
  results/empirical_26Q1/full \
  figures/empirical_26Q1 \
  supplementary \
  -type f ! -name SHA256SUMS.txt -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > results/empirical_26Q1/full/SHA256SUMS.txt
```

Validate later with:

```bash
sha256sum -c results/empirical_26Q1/full/SHA256SUMS.txt
```

## 16. One-command continuation after the current GDC download

Once the 1,997-file download has completed successfully, all subsequent stages can
be run with:

```bash
bash scripts/run_real_data_pipeline.sh after-download
```

The complete pipeline from DepMap validation through final assets is:

```bash
bash scripts/run_real_data_pipeline.sh all
```

Do not use `all` while a separate GDC downloader is already running.

## 17. Statistical interpretation

For DepMap dependency contrasts:

```text
delta_effect = median CRISPR effect in loss models
               - median CRISPR effect in intact models
```

A negative value supports stronger dependency in the loss cohort. The strongest
screen-level evidence combines:

```text
n_loss >= 3
n_intact >= 3
delta_effect < 0
q_value_bh < 0.05
```

For expression compensation:

```text
delta_expression > 0
```

indicates higher target expression in the loss cohort. This is supportive evidence,
not proof of synthetic lethality.

For TCGA, `tcga_homdel_frequency` is the proportion of evaluable primary tumors
with total integer copy number equal to zero for the lost gene.

## 18. Claims that must not be made automatically

The pipeline does not establish:

- clinical efficacy;
- direct pairwise rescue;
- causality from association alone;
- MSI, HRD or other complex phenotypes from a single copy-number variable;
- a negative biological result when a required empirical domain is missing.

All high-priority candidates require biomarker-confirmed isogenic, pharmacological,
rescue and mechanistic validation.

## 19. Data-management rules

Do not commit large third-party raw data to GitHub. The repository tracks:

- acquisition and validation scripts;
- release identifiers and reviewed manifests;
- checksums and provenance records;
- analysis code;
- small result summaries and publication assets when redistribution is permitted.

The raw DepMap and GDC files remain subject to their original terms.
