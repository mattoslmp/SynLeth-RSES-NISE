#!/usr/bin/env bash
# Reproducible real-data workflow for DepMap Public 26Q1 and open TCGA/GDC ASCAT3 data.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DEPMAP_DIR="${DEPMAP_DIR:-data/raw/depmap}"
GDC_DIR="${GDC_DIR:-data/raw/gdc}"
PROCESSED_DIR="${PROCESSED_DIR:-data/processed}"
DEPMAP_RESULTS="${DEPMAP_RESULTS:-results/empirical_26Q1/depmap_only}"
FULL_RESULTS="${FULL_RESULTS:-results/empirical_26Q1/full}"
EXPANDED_ROOT="${EXPANDED_ROOT:-results/expanded_26Q1}"
EXPANDED_DEPMAP_RESULTS="${EXPANDED_DEPMAP_RESULTS:-$EXPANDED_ROOT/depmap_only}"
EXPANDED_FULL_RESULTS="${EXPANDED_FULL_RESULTS:-$EXPANDED_ROOT/full}"
FIGURE_DIR="${FIGURE_DIR:-figures/empirical_26Q1}"
EXPANDED_FIGURE_DIR="${EXPANDED_FIGURE_DIR:-figures/expanded_26Q1}"
LOG_DIR="${LOG_DIR:-logs/empirical_26Q1}"
LOSS_THRESHOLD="${LOSS_THRESHOLD:-0.30}"
MIN_GROUP_SIZE="${MIN_GROUP_SIZE:-3}"
GDC_WORKFLOW="${GDC_WORKFLOW:-ASCAT3}"

GENE_EFFECT="${GENE_EFFECT:-$DEPMAP_DIR/CRISPRGeneEffect.csv}"
COPY_NUMBER="${COPY_NUMBER:-$DEPMAP_DIR/OmicsCNGeneWGS.csv}"
MODELS="${MODELS:-$DEPMAP_DIR/Model.csv}"
EXPRESSION="${EXPRESSION:-$DEPMAP_DIR/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv}"
GDC_MANIFEST="${GDC_MANIFEST:-$GDC_DIR/gdc_gene_level_copy_number_manifest.json}"
EXPANDED_CANDIDATES="${EXPANDED_CANDIDATES:-$PROCESSED_DIR/expanded_candidate_universe.tsv}"
EXPANDED_EVIDENCE="${EXPANDED_EVIDENCE:-$PROCESSED_DIR/expanded_pair_functional_evidence.tsv}"
ENSEMBL_PARALOGS="${ENSEMBL_PARALOGS:-data/raw/ensembl/human_seed_paralogs.tsv}"

mkdir -p "$DEPMAP_DIR" "$GDC_DIR" "$PROCESSED_DIR" \
  "$DEPMAP_RESULTS" "$FULL_RESULTS" "$EXPANDED_DEPMAP_RESULTS" \
  "$EXPANDED_FULL_RESULTS" "$FIGURE_DIR" "$EXPANDED_FIGURE_DIR" \
  "$LOG_DIR" supplementary

log_stage() {
  printf '\n[%s] %s\n' "$(date -Iseconds)" "$*"
}

run_logged() {
  local log_file="$1"
  shift
  set -o pipefail
  PYTHONUNBUFFERED=1 "$@" 2>&1 | tee "$log_file"
  local status=${PIPESTATUS[0]}
  if [[ $status -ne 0 ]]; then
    echo "Command failed with status $status; see $log_file" >&2
    return "$status"
  fi
}

validate_depmap() {
  log_stage "Validate DepMap files and checksums"
  run_logged "$LOG_DIR/01_validate_depmap_files.log" \
    python -u scripts/download_depmap.py \
      --input-dir "$DEPMAP_DIR" \
      --write-checksums

  log_stage "Validate DepMap ModelID and cancer-cohort crosswalk"
  run_logged "$LOG_DIR/02_validate_depmap_crosswalk.log" \
    python -u scripts/validate_real_inputs.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --expression "$EXPRESSION"
}

run_depmap() {
  log_stage "Run benchmark DepMap-only empirical RSES-Onco"
  run_logged "$LOG_DIR/03_run_depmap_only.log" \
    python -u scripts/run_empirical_rses_onco.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --expression "$EXPRESSION" \
      --loss-threshold "$LOSS_THRESHOLD" \
      --min-group-size "$MIN_GROUP_SIZE" \
      --output "$DEPMAP_RESULTS/empirical_rses_onco_by_cancer.tsv"
}

build_universe() {
  log_stage "Build every directed human NISE pair plus curated benchmark classes"
  local command=(
    python -u scripts/build_expanded_candidate_universe.py
    --output "$EXPANDED_CANDIDATES"
    --members-output "$PROCESSED_DIR/expanded_class_member_inventory.tsv"
  )
  if [[ -s "$ENSEMBL_PARALOGS" ]]; then
    command+=(--additional "$ENSEMBL_PARALOGS")
  fi
  run_logged "$LOG_DIR/03a_build_expanded_candidate_universe.log" "${command[@]}"
}

expand_paralogs() {
  log_stage "Identify Ensembl Compara paralogs for every seed gene"
  run_logged "$LOG_DIR/03b_download_ensembl_paralogs.log" \
    python -u scripts/download_ensembl_paralogs.py \
      --candidates "$EXPANDED_CANDIDATES" \
      --output "$ENSEMBL_PARALOGS"

  log_stage "Rebuild candidate universe with the Ensembl paralog class"
  run_logged "$LOG_DIR/03c_rebuild_expanded_candidate_universe.log" \
    python -u scripts/build_expanded_candidate_universe.py \
      --additional "$ENSEMBL_PARALOGS" \
      --output "$EXPANDED_CANDIDATES" \
      --members-output "$PROCESSED_DIR/expanded_class_member_inventory.tsv"
}

functional_evidence() {
  log_stage "Acquire STRING, DoRothEA, HPA and UniProt evidence for all candidate genes"
  run_logged "$LOG_DIR/03d_download_human_functional_evidence.log" \
    python -u scripts/download_human_functional_evidence.py \
      --candidates "$EXPANDED_CANDIDATES" \
      --output "$EXPANDED_EVIDENCE"
}

run_expanded_depmap() {
  log_stage "Run all-NISE and multi-class expanded DepMap analysis"
  run_logged "$LOG_DIR/03e_run_expanded_depmap.log" \
    python -u scripts/run_expanded_rses_onco.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --expression "$EXPRESSION" \
      --candidates "$EXPANDED_CANDIDATES" \
      --functional-evidence "$EXPANDED_EVIDENCE" \
      --loss-threshold "$LOSS_THRESHOLD" \
      --min-group-size "$MIN_GROUP_SIZE" \
      --output "$EXPANDED_DEPMAP_RESULTS/expanded_rses_onco.tsv"
}

manifest_gdc() {
  log_stage "Query open GDC ASCAT3 primary-tumor gene-level copy-number files"
  run_logged "$LOG_DIR/04_gdc_manifest.log" \
    python -u scripts/download_gdc.py \
      --manifest-only \
      --workflow "$GDC_WORKFLOW" \
      --output-dir "$GDC_DIR"

  cp -f "$GDC_MANIFEST" \
    "$GDC_DIR/gdc_gene_level_copy_number_manifest_$(date +%Y%m%d).json"
}

download_gdc() {
  log_stage "Download or resume GDC files using the existing reviewed manifest"
  run_logged "$LOG_DIR/05_gdc_download.log" \
    python -u scripts/download_gdc.py \
      --use-existing-manifest \
      --manifest "$GDC_MANIFEST" \
      --workflow "$GDC_WORKFLOW" \
      --output-dir "$GDC_DIR" \
      --retries 3
}

validate_gdc() {
  log_stage "Validate every GDC file by size and MD5"
  run_logged "$LOG_DIR/06_gdc_validate.log" \
    python -u scripts/download_gdc.py \
      --validate-only \
      --manifest "$GDC_MANIFEST" \
      --output-dir "$GDC_DIR"
}

aggregate_gdc() {
  log_stage "Aggregate GDC files into deletion-only matrices"
  run_logged "$LOG_DIR/07_gdc_aggregate.log" \
    python -u scripts/aggregate_gdc_gene_cna.py \
      --raw-dir "$GDC_DIR" \
      --manifest "$GDC_MANIFEST" \
      --output-dir "$PROCESSED_DIR"
}

validate_matrices() {
  log_stage "Validate aggregated TCGA matrices and candidate-gene coverage"
  run_logged "$LOG_DIR/08_validate_gdc_matrices.log" \
    python -u scripts/validate_gdc_matrices.py \
      --output "$FULL_RESULTS/gdc_matrix_qc.tsv" \
      --event-output "$FULL_RESULTS/tcga_gene_event_summary.tsv"
}

run_full() {
  log_stage "Run benchmark integrated TCGA plus DepMap RSES-Onco"
  run_logged "$LOG_DIR/09_run_full_tcga_depmap.log" \
    python -u scripts/run_empirical_rses_onco.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --expression "$EXPRESSION" \
      --loss-threshold "$LOSS_THRESHOLD" \
      --min-group-size "$MIN_GROUP_SIZE" \
      --tcga "colon=$PROCESSED_DIR/TCGA_COLON_homdel_discrete.tsv" \
      --tcga "stomach=$PROCESSED_DIR/TCGA_STOMACH_homdel_discrete.tsv" \
      --tcga "lung=$PROCESSED_DIR/TCGA_LUNG_homdel_discrete.tsv" \
      --output "$FULL_RESULTS/empirical_rses_onco_by_cancer.tsv"
}

run_expanded_full() {
  log_stage "Run all-NISE, paralog and benchmark integrated TCGA plus DepMap RSES-Onco"
  run_logged "$LOG_DIR/09a_run_expanded_tcga_depmap.log" \
    python -u scripts/run_expanded_rses_onco.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --expression "$EXPRESSION" \
      --candidates "$EXPANDED_CANDIDATES" \
      --functional-evidence "$EXPANDED_EVIDENCE" \
      --loss-threshold "$LOSS_THRESHOLD" \
      --min-group-size "$MIN_GROUP_SIZE" \
      --tcga "colon=$PROCESSED_DIR/TCGA_COLON_homdel_discrete.tsv" \
      --tcga "stomach=$PROCESSED_DIR/TCGA_STOMACH_homdel_discrete.tsv" \
      --tcga "lung=$PROCESSED_DIR/TCGA_LUNG_homdel_discrete.tsv" \
      --output "$EXPANDED_FULL_RESULTS/expanded_rses_onco.tsv"
}

summarize() {
  log_stage "Create manuscript-ready benchmark result tables"
  run_logged "$LOG_DIR/10_summarize_empirical_results.log" \
    python -u scripts/summarize_empirical_results.py \
      --input "$FULL_RESULTS/empirical_rses_onco_by_cancer.tsv" \
      --depmap-only "$DEPMAP_RESULTS/empirical_rses_onco_by_cancer.tsv" \
      --output-dir "$FULL_RESULTS/article_tables"
}

summarize_expanded() {
  log_stage "Create manuscript-ready all-NISE and multi-class result tables"
  run_logged "$LOG_DIR/10a_summarize_expanded_results.log" \
    python -u scripts/summarize_empirical_results.py \
      --input "$EXPANDED_FULL_RESULTS/expanded_rses_onco.tsv" \
      --depmap-only "$EXPANDED_DEPMAP_RESULTS/expanded_rses_onco.tsv" \
      --output-dir "$EXPANDED_FULL_RESULTS/article_tables"
}

figures() {
  log_stage "Generate benchmark empirical figures in PDF, PNG and SVG"
  run_logged "$LOG_DIR/11_make_empirical_figures.log" \
    python -u scripts/make_empirical_figures.py \
      --input "$FULL_RESULTS/empirical_rses_onco_by_cancer.tsv" \
      --output-dir "$FIGURE_DIR"
}

figures_expanded() {
  log_stage "Generate all-NISE and multi-class figures in PDF, PNG and SVG"
  run_logged "$LOG_DIR/11a_make_expanded_figures.log" \
    python -u scripts/make_empirical_figures.py \
      --input "$EXPANDED_FULL_RESULTS/expanded_rses_onco.tsv" \
      --output-dir "$EXPANDED_FIGURE_DIR"
}

workbook() {
  log_stage "Build benchmark empirical supplementary workbook"
  run_logged "$LOG_DIR/12_build_empirical_workbook.log" \
    python -u scripts/build_empirical_workbook.py \
      --input "$FULL_RESULTS/empirical_rses_onco_by_cancer.tsv" \
      --output supplementary/RSES_Onco_Empirical_26Q1.xlsx
}

workbook_expanded() {
  log_stage "Build all-NISE and multi-class supplementary workbook"
  run_logged "$LOG_DIR/12a_build_expanded_workbook.log" \
    python -u scripts/build_empirical_workbook.py \
      --input "$EXPANDED_FULL_RESULTS/expanded_rses_onco.tsv" \
      --output supplementary/RSES_Onco_Expanded_All_NISE_26Q1.xlsx
}

finalize() {
  log_stage "Run tests and generate reproducibility manifests"
  run_logged "$LOG_DIR/13_pytest.log" \
    python -m pytest -q -p no:cacheprovider

  {
    printf 'path\tsize_bytes\n'
    find "$DEPMAP_RESULTS" "$FULL_RESULTS" "$EXPANDED_ROOT" \
      "$FIGURE_DIR" "$EXPANDED_FIGURE_DIR" supplementary \
      -type f -printf '%p\t%s\n' | sort
  } > "$EXPANDED_FULL_RESULTS/empirical_output_manifest.tsv"

  find "$DEPMAP_RESULTS" "$FULL_RESULTS" "$EXPANDED_ROOT" \
    "$FIGURE_DIR" "$EXPANDED_FIGURE_DIR" supplementary \
    -type f ! -name 'SHA256SUMS.txt' -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "$EXPANDED_FULL_RESULTS/SHA256SUMS.txt"

  echo "Wrote $EXPANDED_FULL_RESULTS/empirical_output_manifest.tsv"
  echo "Wrote $EXPANDED_FULL_RESULTS/SHA256SUMS.txt"
}

usage() {
  cat <<'EOF'
Usage: bash scripts/run_real_data_pipeline.sh STAGE

Benchmark stages:
  validate-depmap     Validate DepMap files, schemas, IDs and cancer cohorts
  run-depmap          Run benchmark DepMap-only analysis
  manifest-gdc        Query and save the GDC manifest
  download-gdc        Download/resume files from the existing GDC manifest
  validate-gdc        Validate all downloaded files by size and MD5
  aggregate-gdc       Build project and cancer-level homdel matrices
  validate-matrices   Validate values, dimensions, duplicates and gene coverage
  run-full            Run benchmark integrated TCGA plus DepMap scoring
  summarize           Create benchmark manuscript-ready tables
  figures             Generate benchmark PDF/PNG/SVG figures
  workbook            Build benchmark supplementary workbook

Expanded v0.8 stages:
  build-universe      Add all 202 directed NISE hypotheses plus benchmark classes
  expand-paralogs     Add all Ensembl Compara paralogs of analyzed seed genes
  functional-evidence Acquire STRING, DoRothEA, HPA and UniProt evidence
  run-expanded-depmap Run expanded DepMap expression/phenotype profiles and dependencies
  run-expanded-full   Run expanded integrated TCGA plus DepMap scoring
  summarize-expanded Create all-NISE and multi-class manuscript tables
  figures-expanded   Generate expanded PDF/PNG/SVG figures
  workbook-expanded  Build expanded supplementary workbook
  expanded-setup     Build universe, paralogs, functional evidence and DepMap profiles
  expanded-after-download Run all stages after an existing GDC download
  expanded-all       Run the entire expanded workflow including GDC download

Shared:
  finalize            Run tests and write manifests/checksums
  after-download      Run the legacy benchmark workflow after GDC download
  all                 Run the legacy benchmark workflow from the beginning
EOF
}

stage="${1:-}"
case "$stage" in
  validate-depmap) validate_depmap ;;
  run-depmap) run_depmap ;;
  build-universe) build_universe ;;
  expand-paralogs) build_universe; expand_paralogs ;;
  functional-evidence) functional_evidence ;;
  run-expanded-depmap) run_expanded_depmap ;;
  manifest-gdc) manifest_gdc ;;
  download-gdc) download_gdc ;;
  validate-gdc) validate_gdc ;;
  aggregate-gdc) aggregate_gdc ;;
  validate-matrices) validate_matrices ;;
  run-full) run_full ;;
  run-expanded-full) run_expanded_full ;;
  summarize) summarize ;;
  summarize-expanded) summarize_expanded ;;
  figures) figures ;;
  figures-expanded) figures_expanded ;;
  workbook) workbook ;;
  workbook-expanded) workbook_expanded ;;
  finalize) finalize ;;
  expanded-setup)
    build_universe
    expand_paralogs
    functional_evidence
    run_expanded_depmap
    ;;
  after-download)
    validate_gdc
    aggregate_gdc
    validate_matrices
    run_full
    summarize
    figures
    workbook
    finalize
    ;;
  expanded-after-download)
    build_universe
    expand_paralogs
    functional_evidence
    run_expanded_depmap
    validate_gdc
    aggregate_gdc
    validate_matrices
    run_expanded_full
    summarize_expanded
    figures_expanded
    workbook_expanded
    finalize
    ;;
  all)
    validate_depmap
    run_depmap
    manifest_gdc
    download_gdc
    validate_gdc
    aggregate_gdc
    validate_matrices
    run_full
    summarize
    figures
    workbook
    finalize
    ;;
  expanded-all)
    validate_depmap
    build_universe
    expand_paralogs
    functional_evidence
    run_expanded_depmap
    manifest_gdc
    download_gdc
    validate_gdc
    aggregate_gdc
    validate_matrices
    run_expanded_full
    summarize_expanded
    figures_expanded
    workbook_expanded
    finalize
    ;;
  *) usage; exit 2 ;;
esac
