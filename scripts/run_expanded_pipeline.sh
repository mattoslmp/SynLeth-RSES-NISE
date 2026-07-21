#!/usr/bin/env bash
# End-to-end expanded RSES-Onco v0.8 workflow.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DEPMAP_DIR="${DEPMAP_DIR:-data/raw/depmap}"
GDC_DIR="${GDC_DIR:-data/raw/gdc}"
PROCESSED_DIR="${PROCESSED_DIR:-data/processed}"
RESULT_ROOT="${RESULT_ROOT:-results/expanded_26Q1}"
DEPMAP_RESULTS="${DEPMAP_RESULTS:-$RESULT_ROOT/depmap_only}"
FULL_RESULTS="${FULL_RESULTS:-$RESULT_ROOT/full}"
DISCOVERY_RESULTS="${DISCOVERY_RESULTS:-$RESULT_ROOT/discovery}"
FIGURE_DIR="${FIGURE_DIR:-figures/expanded_26Q1}"
LOG_DIR="${LOG_DIR:-logs/expanded_26Q1}"

GENE_EFFECT="${GENE_EFFECT:-$DEPMAP_DIR/CRISPRGeneEffect.csv}"
COPY_NUMBER="${COPY_NUMBER:-$DEPMAP_DIR/OmicsCNGeneWGS.csv}"
MODELS="${MODELS:-$DEPMAP_DIR/Model.csv}"
EXPRESSION="${EXPRESSION:-$DEPMAP_DIR/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv}"
GDC_MANIFEST="${GDC_MANIFEST:-$GDC_DIR/gdc_gene_level_copy_number_manifest.json}"
CANDIDATES="${CANDIDATES:-$PROCESSED_DIR/expanded_candidate_universe.tsv}"
MEMBERS="${MEMBERS:-$PROCESSED_DIR/expanded_class_member_inventory.tsv}"
PARALOGS="${PARALOGS:-data/raw/ensembl/human_seed_paralogs.tsv}"
DISCOVERED="${DISCOVERED:-data/raw/discovery/depmap_discovered_candidate_pairs.tsv}"
FUNCTIONAL_EVIDENCE="${FUNCTIONAL_EVIDENCE:-$PROCESSED_DIR/expanded_pair_functional_evidence.tsv}"
LOSS_THRESHOLD="${LOSS_THRESHOLD:-0.30}"
MIN_GROUP_SIZE="${MIN_GROUP_SIZE:-3}"
DISCOVERY_FDR="${DISCOVERY_FDR:-0.10}"
DISCOVERY_MIN_DELTA="${DISCOVERY_MIN_DELTA:-0.15}"

mkdir -p "$PROCESSED_DIR" "$DEPMAP_RESULTS" "$FULL_RESULTS" \
  "$DISCOVERY_RESULTS" "$FIGURE_DIR" "$LOG_DIR" supplementary

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
  log_stage "Validate DepMap 26Q1 files, checksums and cancer cohorts"
  run_logged "$LOG_DIR/01_validate_depmap_files.log" \
    python -u scripts/download_depmap.py \
      --input-dir "$DEPMAP_DIR" \
      --write-checksums
  run_logged "$LOG_DIR/02_validate_depmap_crosswalk.log" \
    python -u scripts/validate_real_inputs.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --expression "$EXPRESSION"
}

build_base_universe() {
  log_stage "Build all 202 directed NISE hypotheses and curated benchmark classes"
  run_logged "$LOG_DIR/03_build_all_nise_universe.log" \
    python -u scripts/build_expanded_candidate_universe.py \
      --output "$CANDIDATES" \
      --members-output "$MEMBERS"
}

expand_paralogs() {
  log_stage "Retrieve all Ensembl Compara paralogs for current seed genes"
  run_logged "$LOG_DIR/04_ensembl_paralogs.log" \
    python -u scripts/download_ensembl_paralogs.py \
      --candidates "$CANDIDATES" \
      --output "$PARALOGS"

  log_stage "Rebuild candidate universe with Ensembl paralogs"
  run_logged "$LOG_DIR/05_rebuild_with_paralogs.log" \
    python -u scripts/build_expanded_candidate_universe.py \
      --additional "$PARALOGS" \
      --output "$CANDIDATES" \
      --members-output "$MEMBERS"
}

discover_all_targets() {
  log_stage "Screen every DepMap CRISPR target for every analyzable candidate loss gene"
  run_logged "$LOG_DIR/06_all_target_dependency_discovery.log" \
    python -u scripts/discover_conditional_dependencies.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --candidates "$CANDIDATES" \
      --loss-universe candidates \
      --loss-threshold "$LOSS_THRESHOLD" \
      --min-group-size "$MIN_GROUP_SIZE" \
      --minimum-delta "$DISCOVERY_MIN_DELTA" \
      --fdr "$DISCOVERY_FDR" \
      --output "$DISCOVERY_RESULTS/all_target_dependency_screen.tsv" \
      --candidate-output "$DISCOVERED"

  log_stage "Rebuild candidate universe with supported all-target discoveries"
  local command=(
    python -u scripts/build_expanded_candidate_universe.py
    --additional "$PARALOGS"
    --output "$CANDIDATES"
    --members-output "$MEMBERS"
  )
  if [[ -s "$DISCOVERED" ]]; then
    command+=(--additional "$DISCOVERED")
  fi
  run_logged "$LOG_DIR/07_rebuild_with_discoveries.log" "${command[@]}"
}

acquire_functional_evidence() {
  log_stage "Acquire STRING, DoRothEA, HPA and UniProt evidence for every candidate gene"
  run_logged "$LOG_DIR/08_human_functional_evidence.log" \
    python -u scripts/download_human_functional_evidence.py \
      --candidates "$CANDIDATES" \
      --raw-dir data/raw/human_functional_evidence \
      --output "$FUNCTIONAL_EVIDENCE"
}

run_depmap_expanded() {
  log_stage "Compute expression, CRISPR phenotype, network and expanded DepMap scores"
  run_logged "$LOG_DIR/09_expanded_depmap.log" \
    python -u scripts/run_expanded_rses_onco.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --expression "$EXPRESSION" \
      --candidates "$CANDIDATES" \
      --functional-evidence "$FUNCTIONAL_EVIDENCE" \
      --loss-threshold "$LOSS_THRESHOLD" \
      --min-group-size "$MIN_GROUP_SIZE" \
      --output "$DEPMAP_RESULTS/expanded_rses_onco.tsv"
}

validate_and_aggregate_gdc() {
  log_stage "Validate all reviewed-manifest GDC files"
  run_logged "$LOG_DIR/10_validate_gdc.log" \
    python -u scripts/download_gdc.py \
      --validate-only \
      --manifest "$GDC_MANIFEST" \
      --output-dir "$GDC_DIR"

  log_stage "Aggregate ASCAT3 files into project and cancer deletion-only matrices"
  run_logged "$LOG_DIR/11_aggregate_gdc.log" \
    python -u scripts/aggregate_gdc_gene_cna.py \
      --raw-dir "$GDC_DIR" \
      --manifest "$GDC_MANIFEST" \
      --output-dir "$PROCESSED_DIR"

  log_stage "Validate TCGA matrix dimensions, values, duplicates and gene coverage"
  run_logged "$LOG_DIR/12_validate_gdc_matrices.log" \
    python -u scripts/validate_gdc_matrices.py \
      --output "$FULL_RESULTS/gdc_matrix_qc.tsv" \
      --event-output "$FULL_RESULTS/tcga_gene_event_summary.tsv"
}

run_full_expanded() {
  log_stage "Run all-NISE/multi-class TCGA plus DepMap scoring"
  run_logged "$LOG_DIR/13_expanded_tcga_depmap.log" \
    python -u scripts/run_expanded_rses_onco.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --expression "$EXPRESSION" \
      --candidates "$CANDIDATES" \
      --functional-evidence "$FUNCTIONAL_EVIDENCE" \
      --loss-threshold "$LOSS_THRESHOLD" \
      --min-group-size "$MIN_GROUP_SIZE" \
      --tcga "colon=$PROCESSED_DIR/TCGA_COLON_homdel_discrete.tsv" \
      --tcga "stomach=$PROCESSED_DIR/TCGA_STOMACH_homdel_discrete.tsv" \
      --tcga "lung=$PROCESSED_DIR/TCGA_LUNG_homdel_discrete.tsv" \
      --output "$FULL_RESULTS/expanded_rses_onco.tsv"
}

build_article_assets() {
  log_stage "Build expanded manuscript-ready tables"
  run_logged "$LOG_DIR/14_summarize_expanded.log" \
    python -u scripts/summarize_empirical_results.py \
      --input "$FULL_RESULTS/expanded_rses_onco.tsv" \
      --depmap-only "$DEPMAP_RESULTS/expanded_rses_onco.tsv" \
      --output-dir "$FULL_RESULTS/article_tables"

  log_stage "Build expanded PDF, PNG and SVG figures"
  run_logged "$LOG_DIR/15_figures_expanded.log" \
    python -u scripts/make_empirical_figures.py \
      --input "$FULL_RESULTS/expanded_rses_onco.tsv" \
      --output-dir "$FIGURE_DIR"

  log_stage "Build expanded supplementary workbook"
  run_logged "$LOG_DIR/16_workbook_expanded.log" \
    python -u scripts/build_empirical_workbook.py \
      --input "$FULL_RESULTS/expanded_rses_onco.tsv" \
      --output supplementary/RSES_Onco_Expanded_All_NISE_26Q1.xlsx
}

finalize() {
  log_stage "Run tests, output manifest and checksums"
  run_logged "$LOG_DIR/17_pytest.log" \
    python -m pytest -q -p no:cacheprovider

  {
    printf 'path\tsize_bytes\n'
    find "$RESULT_ROOT" "$FIGURE_DIR" supplementary \
      -type f -printf '%p\t%s\n' | sort
  } > "$FULL_RESULTS/expanded_output_manifest.tsv"

  find "$RESULT_ROOT" "$FIGURE_DIR" supplementary \
    -type f ! -name SHA256SUMS.txt -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "$FULL_RESULTS/SHA256SUMS.txt"

  echo "Wrote $FULL_RESULTS/expanded_output_manifest.tsv"
  echo "Wrote $FULL_RESULTS/SHA256SUMS.txt"
}

setup() {
  validate_depmap
  build_base_universe
  expand_paralogs
  discover_all_targets
  acquire_functional_evidence
  run_depmap_expanded
}

after_download() {
  setup
  validate_and_aggregate_gdc
  run_full_expanded
  build_article_assets
  finalize
}

all() {
  setup
  log_stage "Create reviewed GDC manifest"
  run_logged "$LOG_DIR/10a_manifest_gdc.log" \
    python -u scripts/download_gdc.py \
      --manifest-only \
      --workflow ASCAT3 \
      --output-dir "$GDC_DIR"
  log_stage "Download/resume reviewed GDC manifest"
  run_logged "$LOG_DIR/10b_download_gdc.log" \
    python -u scripts/download_gdc.py \
      --use-existing-manifest \
      --manifest "$GDC_MANIFEST" \
      --workflow ASCAT3 \
      --output-dir "$GDC_DIR" \
      --retries 3
  validate_and_aggregate_gdc
  run_full_expanded
  build_article_assets
  finalize
}

usage() {
  cat <<'EOF'
Usage: bash scripts/run_expanded_pipeline.sh STAGE

Stages:
  setup           Validate DepMap, build all NISEs/paralogs, discover all-target
                  conditional dependencies, acquire networks, run expanded DepMap
  after-download  Run setup and then continue from an already completed GDC download
  all             Run the complete expanded workflow including GDC acquisition
  discover-all-cn Optional broad screen using every CN gene with sufficient losses
EOF
}

stage="${1:-}"
case "$stage" in
  setup) setup ;;
  after-download) after_download ;;
  all) all ;;
  discover-all-cn)
    build_base_universe
    run_logged "$LOG_DIR/06b_all_cn_all_target_discovery.log" \
      python -u scripts/discover_conditional_dependencies.py \
        --gene-effect "$GENE_EFFECT" \
        --copy-number "$COPY_NUMBER" \
        --models "$MODELS" \
        --candidates "$CANDIDATES" \
        --loss-universe all-cn \
        --loss-threshold "$LOSS_THRESHOLD" \
        --min-group-size "$MIN_GROUP_SIZE" \
        --minimum-delta "$DISCOVERY_MIN_DELTA" \
        --fdr "$DISCOVERY_FDR" \
        --output "$DISCOVERY_RESULTS/all_cn_all_target_dependency_screen.tsv" \
        --candidate-output data/raw/discovery/depmap_all_cn_discovered_pairs.tsv
    ;;
  *) usage; exit 2 ;;
esac
