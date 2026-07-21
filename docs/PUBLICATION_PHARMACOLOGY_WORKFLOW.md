# Publication and pharmacology workflow (RSES-Onco v0.9)

## Scope

This workflow generates every main figure, supplementary figure, main table,
supplementary table, source-data file, workbook and reproducibility manifest from
version-controlled scripts. It also adds a coverage-aware pharmacology layer for
prioritizing experimental target-compound hypotheses.

The pharmacology output is not a clinical recommendation, efficacy claim or
promise of cure. It is a transparent, source-linked prioritization layer for
biomarker-matched experimental validation.

## New scripts

### Pharmacology

| Script | Purpose |
|---|---|
| `scripts/acquire_pharmacology_evidence.py` | Queries and caches Open Targets, ChEMBL, DGIdb, MyChem, Pharos/TCRD and CIViC evidence for prioritized targets. |
| `scripts/standardize_drug_sensitivity.py` | Converts locally downloaded PRISM, GDSC and CTRP releases into one ModelID-linked long table. |
| `scripts/analyze_drug_response_selectivity.py` | Tests whether candidate compounds show stronger response in biomarker-loss than intact cancer models. |
| `scripts/prioritize_pharmacology.py` | Calculates coverage-aware pharmacology and combined therapeutic-hypothesis scores. |

### Figures and publication assets

| Script | Purpose |
|---|---|
| `scripts/make_main_figures.py` | Generates Figures 1–7 from source tables only. |
| `scripts/make_supplementary_figures.py` | Generates Figures S1–S14 from source tables only. |
| `scripts/make_all_article_figures.py` | Runs all figure scripts and verifies all 21 figures and all 63 PNG/PDF/SVG files. |
| `scripts/export_article_tables.py` | Exports 4 main and 15 supplementary tables. |
| `scripts/build_article_workbook.py` | Builds a formatted Excel workbook containing all article tables and manifests. |
| `scripts/build_publication_manifest.py` | Writes inventory, provenance, input fingerprints and SHA-256 checksums. |
| `scripts/validate_publication_outputs.py` | Fails when a figure, format, source-data file, table, layout audit or workbook is missing. |
| `scripts/run_publication_pipeline.sh` | Orchestrates pharmacology, all figures, all tables, workbook, manifests and validation. |
| `scripts/run_expanded_pipeline.sh` | Runs the complete all-NISE analysis and invokes the publication pipeline after integrated TCGA/DepMap scoring. |

### Shared modules and configuration

| File | Purpose |
|---|---|
| `src/rses_onco/publication.py` | Publication style, label wrapping, dynamic sizing, triplet export and automated overlap/clipping audit. |
| `src/rses_onco/pharmacology.py` | Pharmacology-domain normalization, coverage-aware score and vulnerability-actionability concordance. |
| `config/article_assets.yaml` | Registry of all main/supplementary figures and tables. |
| `config/drug_sensitivity_sources.yaml` | Release-agnostic PRISM/GDSC/CTRP file and column mappings. |

## Figure-quality requirements enforced by code

Every figure is generated in:

```text
PNG: 600 dpi
PDF: vector-compatible embedded TrueType text
SVG: editable text retained
```

The style module sets journal-safe font sizes and dynamic dimensions. Every
figure receives a `*.layout_audit.json` file. In strict mode, generation fails if
the rendered figure contains:

- overlapping axes/panels;
- tick-label collisions;
- text outside the rendered canvas;
- legends outside the canvas;
- missing figure formats.

Automated layout auditing is a safeguard, not a replacement for final visual
inspection at 100% zoom. The final package must still be reviewed manually before
submission.

## Final directory structure

```text
article_outputs/
├── README.md
├── figures/
│   ├── main/
│   │   ├── Figure_1_framework.{png,pdf,svg}
│   │   ├── Figure_2_candidate_universe.{png,pdf,svg}
│   │   ├── Figure_3_cancer_specific_ranking.{png,pdf,svg}
│   │   ├── Figure_4_tcga_depmap_integration.{png,pdf,svg}
│   │   ├── Figure_5_functional_microniches.{png,pdf,svg}
│   │   ├── Figure_6_class_discoveries.{png,pdf,svg}
│   │   └── Figure_7_pharmacology.{png,pdf,svg}
│   └── supplementary/
│       ├── Figure_S1_data_coverage.{png,pdf,svg}
│       ├── ...
│       └── Figure_S14_layout_and_reproducibility_qc.{png,pdf,svg}
├── tables/
│   ├── main/                   # 4 tables
│   └── supplementary/          # 15 tables
├── source_data/
│   ├── figures/
│   │   ├── main/
│   │   └── supplementary/
│   └── tables/
├── manuscript_assets/
│   ├── main_figure_legends.md
│   ├── supplementary_figure_legends.md
│   └── all_figure_legends.md
├── workbooks/
│   └── RSES_Onco_Article_Tables_and_Evidence.xlsx
└── manifests/
    ├── figure_manifest.tsv
    ├── main_figure_manifest.tsv
    ├── supplementary_figure_manifest.tsv
    ├── table_manifest.tsv
    ├── analysis_input_fingerprints.tsv
    ├── publication_file_inventory.tsv
    ├── publication_provenance.json
    └── SHA256SUMS.txt
```

## Main figures

1. **Figure 1 — Expanded RSES-Onco framework.** Candidate universe, human
   microniches, cancer evidence, pharmacology and validation path.
2. **Figure 2 — Candidate universe across classes.** All NISE directions,
   paralogs, curated classes and discovered dependencies.
3. **Figure 3 — Cancer-specific rankings.** Top colorectal, gastric and lung
   vulnerability scores.
4. **Figure 4 — TCGA + DepMap integration.** Tumor event frequency versus
   conditional target selectivity.
5. **Figure 5 — Functional microniches.** Expression, localization,
   biochemical/structural, CRISPR phenotype, STRING and regulation.
6. **Figure 6 — Class discoveries.** Best evidence by class and FDR-supported
   all-target discoveries.
7. **Figure 7 — Pharmacology.** Vulnerability-actionability concordance and top
   experimental target-compound hypotheses.

## Supplementary figures

1. Figure S1 — evidence coverage by domain;
2. Figure S2 — complete NISE activity catalogue;
3. Figure S3 — all-target discovery volcano panels;
4. Figure S4 — complete dependency heatmap;
5. Figure S5 — expression compensation;
6. Figure S6 — CRISPR mutant-phenotype profiles;
7. Figure S7 — expression-context profiles;
8. Figure S8 — STRING network evidence;
9. Figure S9 — DoRothEA regulatory evidence;
10. Figure S10 — HPA localization divergence;
11. Figure S11 — TCGA deletion-event landscape;
12. Figure S12 — pharmacology source coverage;
13. Figure S13 — PRISM/GDSC/CTRP biomarker-selective response;
14. Figure S14 — layout and reproducibility quality control.

## Main and supplementary tables

The registry in `config/article_assets.yaml` defines exactly:

```text
4 main tables
15 supplementary tables
```

These include the complete candidate universe, all directed NISEs, member
inventory, complete score components, dependency and expression contrasts,
CRISPR phenotype and expression profiles, STRING/regulation/localization,
TCGA events, all-target discoveries, raw pharmacology evidence, pharmacology
ranking, drug-response selectivity and source-status/coverage reports.

## Pharmacology evidence sources

### API/cached sources

- Open Targets Platform GraphQL: target tractability and known drug context;
- ChEMBL REST: target mechanisms and quantitative bioactivity;
- DGIdb GraphQL: aggregated drug-gene interactions;
- MyChem.info REST: compound identifier and annotation enrichment;
- Pharos/TCRD GraphQL: target development level and drug ligands;
- CIViC: precision-oncology gene record resolver.

Every request is cached under:

```text
data/raw/pharmacology/api_cache/
```

Source failures are written to:

```text
data/processed/pharmacology/pharmacology_source_status.tsv
```

A failed source remains missing. It is never assigned score zero.

### PRISM, GDSC and CTRP

These resources are handled as local release files because the analysis needs
release-specific response matrices and cell-line metadata. Place files in:

```text
data/raw/pharmacology/prism/
data/raw/pharmacology/gdsc/
data/raw/pharmacology/ctrp/
```

The standardizer reads `config/drug_sensitivity_sources.yaml`. It attempts common
column names and records every accepted, skipped and failed file. For a release
with different column names, update the YAML rather than editing raw data.

The pipeline continues when these optional matrices are absent. The resulting
pharmacology coverage is lower and Figure S13 contains a scripted no-data panel.

## Pharmacology score

The pharmacology score uses available domains only:

```text
target tractability             0.18
direct target-drug interaction  0.18
compound potency                0.18
clinical maturity               0.14
cancer relevance                0.12
biomarker-selective response    0.20
```

The coverage-adjusted pharmacology score is combined with the coverage-adjusted
RSES-Onco vulnerability through geometric concordance. This prevents a highly
druggable but biologically weak target, or a strong vulnerability without any
actionable compound evidence, from being presented as a mature hypothesis.

## Complete execution order after the GDC download

Do not stop the current GDC download. After it finishes and reports exit code
zero, synchronize the code and run the integrated orchestrator:

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
echo "Expanded publication pipeline exit code: $status"
test "$status" -eq 0
```

This performs, in order:

1. DepMap validation;
2. all-NISE candidate construction;
3. Ensembl paralog expansion;
4. all-target DepMap discovery;
5. STRING/DoRothEA/HPA/UniProt acquisition;
6. expanded DepMap scoring;
7. validation of all downloaded GDC files;
8. GDC aggregation;
9. integrated TCGA/DepMap scoring;
10. pharmacology API acquisition;
11. optional PRISM/GDSC/CTRP standardization;
12. biomarker-matched drug-response testing;
13. pharmacology prioritization;
14. all main/supplementary tables;
15. all main/supplementary figures;
16. workbook;
17. manifests, checksums and tests.

## Rebuild only publication assets

After the expanded score exists, pharmacology and article assets can be rebuilt
without repeating TCGA/DepMap analysis:

```bash
bash scripts/run_publication_pipeline.sh all
```

To regenerate only all figures from existing outputs:

```bash
bash scripts/run_publication_pipeline.sh figures
```

or directly:

```bash
python -u scripts/make_all_article_figures.py \
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --discovery results/expanded_26Q1/discovery/all_target_dependency_screen.tsv \
  --pharmacology results/expanded_26Q1/pharmacology/pharmacology_ranked_hypotheses.tsv \
  --output-root article_outputs \
  --strict-layout
```

## Final checks

```bash
python -u scripts/validate_publication_outputs.py \
  --article-root article_outputs

(
  cd article_outputs
  sha256sum -c manifests/SHA256SUMS.txt
)
```

Inspect every final PDF or PNG at 100% zoom before submission. The automated
audit converts common overlap and clipping defects into failures, but human
editorial review remains necessary for density, scientific clarity and journal
formatting.
