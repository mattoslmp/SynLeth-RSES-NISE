# RSES-Onco v0.11.0: canonical end-to-end data, analysis and article protocol

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

## Purpose and scientific boundary

This is the canonical command-by-command protocol for obtaining or validating the required public data, constructing the candidate universe, calculating the coverage-aware RSES-Onco score, generating all supporting analyses, producing every publication table and figure, building the editable documents, rendering the PDFs, performing manual inspection and generating the final package.

RSES-Onco prioritizes experimental hypotheses. It does not establish clinical efficacy, treatment suitability, safety, patient benefit or cure. Missing evidence is not converted to biological zero, and a non-eligible domain does not enter the eligible score denominator.

## Current version contract

```text
Repository and publication framework: 0.11.0
Scientific score: RSES-Onco-expanded-v0.10.9
WGCNA/regulatory semantics: eligibility-aware-wgcna-regulatory-v3
8 main figures
69 supplementary figures
77 registered figures
231 PNG/PDF/SVG exports
4 main tables
44 supplementary tables
48 registered tables
```

Supplementary Figures S68 and S69 must be rendered on different pages.

## 1. Repository preparation

Do not update the working tree while a pipeline from that same directory is active.

```bash
NEW="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010"
OLD="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE"

cd "$NEW" || exit 1

git status --short --branch
git fetch origin
git switch main
git pull --ff-only origin main
git log -1 --decorate --oneline

git rev-parse HEAD | tee logs/analysis_git_commit.txt
```

For a fresh clone:

```bash
git clone https://github.com/mattoslmp/SynLeth-RSES-Onco.git \
  /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010
```

## 2. Environment preparation

```bash
cd "$NEW" || exit 1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco
python -m pip install -e .
```

For a new environment:

```bash
conda env create -f environment.yml
conda activate rses-onco
python -m pip install -e .
```

Required external executables include R/Rscript, WGCNA, FIMO/MEME, PyMOL, LibreOffice, Graphviz and Poppler utilities.

```bash
Rscript -e 'cat(R.version.string, "\n"); cat("WGCNA ", as.character(packageVersion("WGCNA")), "\n", sep="")'
fimo --version
pymol -cq <<'PYMOL'
reinitialize
print "PyMOL headless execution: OK"
quit
PYMOL
command -v libreoffice
command -v pdftotext
command -v pdftoppm
command -v dot
```

## 3. Syntax and test validation

```bash
python -m compileall -q src scripts tests

for script in \
  scripts/run_real_data_pipeline.sh \
  scripts/run_expanded_pipeline.sh \
  scripts/run_expanded_pipeline_portable.sh \
  scripts/resume_wgcna_regulatory_pipeline.sh \
  scripts/run_publication_pipeline.sh \
  scripts/publication_pipeline_complete.sh \
  scripts/publication_pipeline_steps.sh \
  scripts/run_structural_pipeline.sh \
  scripts/verify_complete_article_run.sh
  do
  bash -n "$script" || exit 1
  done

Rscript -e 'parse(file="scripts/run_wgcna_expression_network.R"); cat("R syntax OK\n")'

PYTHONDONTWRITEBYTECODE=1 \
MPLBACKEND=Agg \
python -m pytest -q -p no:cacheprovider
```

## 4. Data-acquisition protocol

The source-specific commands, expected files, provenance rules and recovery procedures are documented in:

- [`DATA_ACQUISITION_AND_REPRODUCTION_V0110.md`](DATA_ACQUISITION_AND_REPRODUCTION_V0110.md)
- [`STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md`](STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md)
- [`DOROTHEA_RECOVERY_WORKFLOW.md`](DOROTHEA_RECOVERY_WORKFLOW.md)
- [`STRUCTURAL_ATLAS_WORKFLOW.md`](STRUCTURAL_ATLAS_WORKFLOW.md)
- [`PUBLICATION_PHARMACOLOGY_WORKFLOW.md`](PUBLICATION_PHARMACOLOGY_WORKFLOW.md)

### 4.1 Required DepMap matrices

Place or symlink these files under `data/raw/depmap/` or set `DEPMAP_DIR`:

```text
CRISPRGeneEffect.csv
OmicsCNGeneWGS.csv
Model.csv
OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
```

Validate them:

```bash
export DEPMAP_DIR="$NEW/data/raw/depmap"

python -u scripts/download_depmap.py \
  --input-dir "$DEPMAP_DIR" \
  --write-checksums

python -u scripts/validate_real_inputs.py \
  --gene-effect "$DEPMAP_DIR/CRISPRGeneEffect.csv" \
  --copy-number "$DEPMAP_DIR/OmicsCNGeneWGS.csv" \
  --models "$DEPMAP_DIR/Model.csv" \
  --expression "$DEPMAP_DIR/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
```

### 4.2 GDC/TCGA acquisition

```bash
export GDC_DIR="$OLD/data/raw/gdc"
mkdir -p "$GDC_DIR" logs/expanded_26Q1

python -u scripts/download_gdc.py \
  --manifest-only \
  --workflow ASCAT3 \
  --output-dir "$GDC_DIR" \
  2>&1 | tee logs/expanded_26Q1/03_gdc_manifest.log

python -u scripts/download_gdc.py \
  --use-existing-manifest \
  --manifest "$GDC_DIR/gdc_gene_level_copy_number_manifest.json" \
  --workflow ASCAT3 \
  --output-dir "$GDC_DIR" \
  --retries 3 \
  2>&1 | tee logs/expanded_26Q1/04_gdc_download.log

python -u scripts/download_gdc.py \
  --validate-only \
  --manifest "$GDC_DIR/gdc_gene_level_copy_number_manifest.json" \
  --output-dir "$GDC_DIR"

find "$GDC_DIR" -type f -name '*.part'
```

No `.part` file may remain before downstream aggregation.

### 4.3 Candidate and evidence acquisition

The full acquisition sequence includes:

```bash
python -u scripts/download_human_nise.py
python -u scripts/build_expanded_candidate_universe.py
python -u scripts/download_ensembl_paralogs.py
python -u scripts/discover_conditional_dependencies.py
python -u scripts/download_human_functional_evidence_resilient.py
python -u scripts/download_hpa_subcellular_current.py
python -u scripts/download_ensembl_promoters.py
python -u scripts/download_jaspar_core_vertebrates.py
python -u scripts/scan_promoter_motifs.py
```

Use the source-specific documents for required arguments and recovery modes. Source failures must be recorded as missing or technical failure, never as biological zero.

## 5. Required input check before the final run

```bash
REQUIRED_INPUTS=(
  "data/raw/depmap/CRISPRGeneEffect.csv"
  "data/raw/depmap/OmicsCNGeneWGS.csv"
  "data/raw/depmap/Model.csv"
  "data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
  "data/processed/expanded_candidate_universe.tsv"
  "data/processed/expanded_class_member_inventory.tsv"
  "data/processed/expanded_pair_functional_evidence_pre_wgcna.tsv"
  "data/processed/TCGA_COLON_homdel_discrete.tsv"
  "data/processed/TCGA_STOMACH_homdel_discrete.tsv"
  "data/processed/TCGA_LUNG_homdel_discrete.tsv"
  "data/raw/human_functional_evidence/omnipath_dorothea.tsv"
  "data/raw/regulatory/ensembl_promoters.tsv"
  "data/raw/regulatory/ensembl_promoters.fa"
  "results/expanded_26Q1/discovery/all_target_dependency_screen.tsv"
)

for path in "${REQUIRED_INPUTS[@]}"; do
  [[ -s "$path" ]] || { echo "Missing or empty: $path" >&2; exit 1; }
  ls -lh "$path"
done
```

Release-specific sample counts may change and must be reported as observed rather than forced to match an older release.

## 6. Safe backup and removal of derived outputs

```bash
STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="$NEW/backups/pre_v0110_${STAMP}"
mkdir -p "$BACKUP/results" "$BACKUP/regulatory" "$BACKUP/article_outputs"

[[ -d results/expanded_26Q1 ]] && rsync -a results/expanded_26Q1/ "$BACKUP/results/"
[[ -d data/processed/regulatory ]] && rsync -a data/processed/regulatory/ "$BACKUP/regulatory/"
[[ -d article_outputs ]] && rsync -a --exclude 'structure_atlas/' article_outputs/ "$BACKUP/article_outputs/"

cp -f \
  data/processed/expanded_pair_functional_evidence_pre_wgcna.tsv \
  data/processed/expanded_pair_functional_evidence.tsv

rm -rf data/processed/regulatory/wgcna
rm -f \
  data/processed/regulatory/expanded_pair_functional_evidence_by_cancer.tsv \
  data/processed/regulatory/promoter_tf_regulatory_pair_metrics.tsv \
  data/processed/regulatory/wgcna_regulatory_layer_status.json
rm -f \
  results/expanded_26Q1/depmap_only/expanded_rses_onco.tsv \
  results/expanded_26Q1/full/expanded_rses_onco.tsv \
  results/expanded_26Q1/full/SHA256SUMS.txt

for directory in figures tables source_data manifests workbooks manuscript_assets review_records documents; do
  rm -rf "article_outputs/$directory"
done
```

The cached promoter, JASPAR/FIMO, pharmacology and structural resources may be preserved when valid.

## 7. Corrected WGCNA/regulatory policy

The primary WGCNA correlation is `bicor` with:

```text
maxPOutliers=0.10
pearsonFallback=individual
networkType=signed
TOMType=signed
```

Pearson is used only for an individual gene or module eigengene where bicor is undefined because the MAD is zero or non-finite. Every fallback is exported in an audit table.

The functional-microniche weights are:

```text
expression_context        0.20
localization              0.15
biochemical_structural    0.15
genetic_phenotype         0.20
interaction_network       0.15
regulatory_network        0.15
```

The final RSES-Onco weights are:

```text
tumor_event               0.16
dependency                0.22
selectivity               0.14
expression_compensation   0.08
functional_relation       0.06
functional_microniche     0.16
validation_tractability   0.18
```

## 8. Execute the complete workflow in tmux

```bash
STAMP="$(date +%Y%m%d_%H%M%S)"
SESSION="rses_v0110_${STAMP}"
RUNNER="$NEW/logs/run_rses_v0110_${STAMP}.sh"
RUN_LOG="$NEW/logs/run_rses_v0110_${STAMP}.log"
EXITCODE_FILE="$NEW/logs/run_rses_v0110_${STAMP}.exitcode"

cat > "$RUNNER" <<EOF
#!/usr/bin/env bash
set -Eeuo pipefail
set -o pipefail
cd "$NEW" || exit 1
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco
export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg
export STRICT_LAYOUT=1
export PUBLICATION_STAGE=assets-only
export DEPMAP_DIR="$NEW/data/raw/depmap"
export GENE_EFFECT="\$DEPMAP_DIR/CRISPRGeneEffect.csv"
export COPY_NUMBER="\$DEPMAP_DIR/OmicsCNGeneWGS.csv"
export MODELS="\$DEPMAP_DIR/Model.csv"
export EXPRESSION="\$DEPMAP_DIR/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
set +e
bash scripts/resume_wgcna_regulatory_pipeline.sh resume-regulatory 2>&1 | tee "$RUN_LOG"
status=\${PIPESTATUS[0]}
set -e
printf '%s\n' "\$status" > "$EXITCODE_FILE"
exit "\$status"
EOF

chmod +x "$RUNNER"
printf '%s\n' "$SESSION" > logs/last_rses_v0110_session.txt
printf '%s\n' "$RUN_LOG" > logs/last_rses_v0110_log.txt
printf '%s\n' "$EXITCODE_FILE" > logs/last_rses_v0110_exitcode_path.txt

tmux new-session -d -s "$SESSION" "bash '$RUNNER'"
```

Attach with:

```bash
tmux attach -t "$(cat logs/last_rses_v0110_session.txt)"
```

Detach without stopping the run with `Ctrl+B`, then `D`.

Monitor with:

```bash
RUN_LOG="$(cat logs/last_rses_v0110_log.txt)"
tail -f "$RUN_LOG"
```

## 9. Completion and WGCNA checks

```bash
EXITCODE_FILE="$(cat logs/last_rses_v0110_exitcode_path.txt)"
RUN_LOG="$(cat logs/last_rses_v0110_log.txt)"
cat "$EXITCODE_FILE"
tail -n 300 "$RUN_LOG"
```

The exit code must be `0`.

For each cancer (`colon`, `stomach`, `lung`), the following must exist:

```text
wgcna_gene_modules.tsv
wgcna_module_eigengenes.tsv
wgcna_pair_metrics.tsv
wgcna_soft_threshold_diagnostics.tsv
wgcna_run_diagnostics.tsv
wgcna_correlation_fallback.tsv
```

## 10. Publication assets

The regulatory resume workflow calls the publication workflow. To regenerate the publication assets independently from existing rankings and caches:

```bash
MPLBACKEND=Agg \
STRICT_LAYOUT=1 \
PYTHONUNBUFFERED=1 \
bash scripts/run_publication_pipeline.sh assets-only
```

The `assets-only` command must generate all figures and tables, exact figure-source TSVs, score-component tables, coverage/missingness audits, overlap controls, robustness analyses, networks, expression and phenotype tables, manifests, workbooks and validation reports.

## 11. Automated validation

```bash
python -u scripts/validate_wgcna_regulatory_evidence.py \
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv \
  --functional-evidence data/processed/expanded_pair_functional_evidence.tsv \
  --article-root article_outputs

python -u scripts/validate_extended_supporting_evidence.py \
  --article-root article_outputs

python -u scripts/validate_publication_scientific_integrity.py \
  --article-root article_outputs \
  --run-marker logs/publication_26Q1/assets_only_run.marker

python -u scripts/validate_publication_outputs.py \
  --article-root article_outputs \
  --run-marker logs/publication_26Q1/assets_only_run.marker

PYTHONDONTWRITEBYTECODE=1 \
MPLBACKEND=Agg \
python -m pytest -q -p no:cacheprovider
```

## 12. Build and render the article documents

```bash
bash scripts/run_publication_pipeline.sh documents
```

Expected outputs:

```text
article_outputs/documents/RSES_Onco_manuscript.docx
article_outputs/documents/RSES_Onco_manuscript.pdf
article_outputs/documents/RSES_Onco_supplementary_material.docx
article_outputs/documents/RSES_Onco_supplementary_material.pdf
article_outputs/documents/RSES_Onco_manuscript_pages/
article_outputs/documents/RSES_Onco_supplementary_material_pages/
```

Validate them:

```bash
python -u scripts/validate_publication_documents.py \
  --article-root article_outputs \
  --document-dir article_outputs/documents \
  --require-page-renders
```

Confirm S68 and S69 are on different pages using `article_outputs/documents/document_figure_page_map.tsv`.

## 13. Manual inspection

Inspect every rendered figure and document page at 100% zoom. Complete:

```text
article_outputs/review_records/MANUAL_VISUAL_INSPECTION_CHECKLIST.tsv
```

Then run:

```bash
python -u scripts/validate_manual_inspection_completion.py \
  --article-root article_outputs
```

Do not auto-fill the checklist.

## 14. Final verification

```bash
PIPELINE_EXITCODE_FILE="$(cat logs/last_rses_v0110_exitcode_path.txt)" \
GDC_DIR="$GDC_DIR" \
MPLBACKEND=Agg \
bash scripts/verify_complete_article_run.sh \
  2>&1 | tee logs/verify_complete_article_run_v0110.log
```

## 15. Final package

```bash
PACKAGE="RSES_Onco_v0110_complete_submission_package.zip"
rm -f "$PACKAGE" "$PACKAGE.sha256"

zip -r -9 "$PACKAGE" \
  article_outputs \
  results/expanded_26Q1 \
  data/processed \
  data/curated \
  data/raw/regulatory \
  config docs manuscript supplementary scripts src tests \
  README.md CITATION.cff pyproject.toml environment.yml LICENSE \
  -x '*/__pycache__/*' '*.pyc' '.git/*' 'backups/*'

sha256sum "$PACKAGE" > "$PACKAGE.sha256"
sha256sum -c "$PACKAGE.sha256"
```

## 16. Companion documentation and manuscript sources

- [`DATA_ACQUISITION_AND_REPRODUCTION_V0110.md`](DATA_ACQUISITION_AND_REPRODUCTION_V0110.md)
- [`../supplementary/Supplementary_Methods_RSES_Onco_v0110.md`](../supplementary/Supplementary_Methods_RSES_Onco_v0110.md)
- [`../manuscript/RSES_Onco_intro_methods_draft_v0110.md`](../manuscript/RSES_Onco_intro_methods_draft_v0110.md)
- [`figures/RSES_Onco_workflow_and_applications.svg`](figures/RSES_Onco_workflow_and_applications.svg)
- `figures/RSES_Onco_workflow_and_applications.png`
