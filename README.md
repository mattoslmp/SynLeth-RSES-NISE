# SynLeth-RSES-Onco / RSES-Onco

**RSES-Onco v0.9** is a coverage-aware framework for discovering and prioritizing
cancer-selective dependencies created by non-homologous isofunctional enzymes
(NISEs), homologous paralogs, pathway backups, collateral deletions and downstream
vulnerabilities. The initial disease scope is colorectal, gastric and lung cancer.

The repository now includes a complete scripted publication workflow and a
pharmacology layer for generating target-compound hypotheses for experimental
validation. It does **not** make clinical recommendations or claims of treatment
efficacy or cure.

## Core scope

- 70 proteins in 15 curated human NISE activities;
- 101 cross-cluster NISE pairs;
- 202 directional NISE hypotheses;
- Ensembl Compara paralog expansion;
- all-target DepMap conditional-dependency discovery;
- DepMap CRISPR, expression and WGS copy number;
- TCGA/GDC ASCAT3 gene-level copy number;
- STRING functional networks;
- OmniPath/DoRothEA regulatory networks;
- Human Protein Atlas localization;
- UniProt/PDB biochemical and structural traceability;
- Open Targets, ChEMBL, DGIdb, MyChem, Pharos/TCRD and CIViC pharmacology evidence;
- optional PRISM, GDSC and CTRP drug-response integration;
- 7 main and 14 supplementary figures generated only by scripts;
- 4 main and 15 supplementary tables;
- source data, figure legends, Excel workbook, provenance and SHA-256 manifests.

## Scientific boundary

Unavailable evidence is not converted to zero. It is excluded from the observed-
domain denominator and lowers explicit coverage. NISE status, paralogy, network
connectivity, druggability or a statistical dependency does not by itself prove
synthetic lethality or therapeutic efficacy. High-priority candidates require
biomarker-matched isogenic perturbation, rescue, orthogonal pharmacology,
mechanistic assays and in vivo validation.

“All NISEs” is exhaustive relative to the bundled curated human NISE catalogue.
Other classes are exhaustive only relative to an explicit source and release.
There is no canonical database containing every possible biological backup or
downstream dependency in the human proteome.

## Installation

```bash
conda env create -f environment.yml
conda activate rses-onco
python -m pip install -e .
python -m pytest -q -p no:cacheprovider
```

## Documentation

- [`docs/REAL_DATA_WORKFLOW.md`](docs/REAL_DATA_WORKFLOW.md)
- [`docs/EXPANDED_HUMAN_EVIDENCE_WORKFLOW.md`](docs/EXPANDED_HUMAN_EVIDENCE_WORKFLOW.md)
- [`docs/ALL_CLASS_AND_ALL_TARGET_DISCOVERY.md`](docs/ALL_CLASS_AND_ALL_TARGET_DISCOVERY.md)
- [`docs/PUBLICATION_PHARMACOLOGY_WORKFLOW.md`](docs/PUBLICATION_PHARMACOLOGY_WORKFLOW.md)
- [`docs/ARTICLE_ANALYSIS_OUTPUTS.md`](docs/ARTICLE_ANALYSIS_OUTPUTS.md)

## Complete workflow after an existing GDC download

Do not restart a GDC download that is already active. After it finishes and
returns exit code zero:

```bash
cd /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE
conda activate rses-onco

git remote set-url origin \
  https://github.com/mattoslmp/SynLeth-RSES-Onco.git

git fetch origin
git checkout main
git pull --ff-only origin main

python -m pip install -e .
python -m pytest -q -p no:cacheprovider

set -o pipefail
PYTHONUNBUFFERED=1 \
bash scripts/run_expanded_pipeline.sh after-download \
  2>&1 | tee logs/run_expanded_after_download_v09.log

status=${PIPESTATUS[0]}
echo "Exit code: $status"
test "$status" -eq 0
```

This command runs all-NISE construction, paralog expansion, all-target discovery,
human network evidence, DepMap, GDC validation and aggregation, TCGA integration,
pharmacology, all article tables, all figures, workbook, manifests and tests.

## Full workflow including GDC acquisition

```bash
bash scripts/run_expanded_pipeline.sh all
```

## Rebuild pharmacology and publication assets only

When the integrated score already exists:

```bash
bash scripts/run_publication_pipeline.sh all
```

Generate only all figures:

```bash
bash scripts/run_publication_pipeline.sh figures
```

The figure orchestrator is:

```bash
python -u scripts/make_all_article_figures.py \
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --discovery results/expanded_26Q1/discovery/all_target_dependency_screen.tsv \
  --pharmacology results/expanded_26Q1/pharmacology/pharmacology_ranked_hypotheses.tsv \
  --output-root article_outputs \
  --strict-layout
```

## Figure generation and quality control

Every registered figure is written as:

```text
PNG: 600 dpi
PDF: vector-compatible
SVG: editable text
```

Each figure also receives a `*.layout_audit.json` file. Strict mode fails on:

- panel/axes overlap;
- tick-label collision;
- clipped text;
- legends outside the figure;
- missing PNG, PDF or SVG;
- missing figure source-data table.

Final visual inspection at 100% zoom remains mandatory before submission.

## Main figures

1. expanded RSES-Onco framework;
2. candidate universe by mechanistic class;
3. colorectal, gastric and lung rankings;
4. TCGA event versus DepMap selectivity;
5. human functional-microniche map;
6. class-specific and all-target discoveries;
7. pharmacological actionability and target-compound hypotheses.

## Supplementary figures

Figures S1–S14 cover evidence availability, all NISE activities, the complete
all-target screen, dependency and expression heatmaps, CRISPR phenotype profiles,
expression context, STRING, regulation, localization, TCGA events, pharmacology
source coverage, PRISM/GDSC/CTRP selectivity and layout/reproducibility quality
control.

## Pharmacology sources

### Live/cached evidence

```text
Open Targets GraphQL
ChEMBL REST
DGIdb GraphQL
MyChem.info REST
Pharos/TCRD GraphQL
CIViC gene record resolver
```

Cached API responses are stored under:

```text
data/raw/pharmacology/api_cache/
```

### Optional drug-response releases

Place release files under:

```text
data/raw/pharmacology/prism/
data/raw/pharmacology/gdsc/
data/raw/pharmacology/ctrp/
```

Column mappings are configured in:

```text
config/drug_sensitivity_sources.yaml
```

If these files are absent, the pipeline records the missing sources and continues
with reduced pharmacology coverage.

## Publication output structure

```text
article_outputs/
├── figures/main/
├── figures/supplementary/
├── tables/main/
├── tables/supplementary/
├── source_data/
├── manuscript_assets/
├── workbooks/
└── manifests/
```

The package includes exactly:

```text
7 main figures × PNG/PDF/SVG
14 supplementary figures × PNG/PDF/SVG
4 main tables
15 supplementary tables
```

## Validate the final package

```bash
python -u scripts/validate_publication_outputs.py \
  --article-root article_outputs

(
  cd article_outputs
  sha256sum -c manifests/SHA256SUMS.txt
)
```

## Main scoring layers

### Human functional-microniche RSES

```text
Expression/context       0.20
Localization             0.15
Biochemical/structural   0.15
Genetic/CRISPR phenotype 0.20
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

### Pharmacology actionability

```text
Target tractability             0.18
Direct target-drug interaction  0.18
Compound potency                0.18
Clinical maturity               0.14
Cancer relevance                0.12
Biomarker-selective response    0.20
```

The final therapeutic-hypothesis score combines coverage-adjusted vulnerability
and coverage-adjusted pharmacology by geometric concordance.

## License

MIT for repository code. Third-party data remain under their original terms and
are acquired or supplied locally rather than redistributed.
