#!/usr/bin/env bash
# Build, score, visualize and validate the RSES-Onco v0.12.0 extended multi-omics layer.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

EXTENDED_DATA_DIR="${EXTENDED_DATA_DIR:-dmap_data}"
EXTENDED_CONFIG="${EXTENDED_CONFIG:-config/extended_multiomics_sources.yaml}"
EXTENDED_PROCESSED_DIR="${EXTENDED_PROCESSED_DIR:-data/processed/extended_multiomics}"
DEPMAP_DIR="${DEPMAP_DIR:-data/raw/depmap}"
MODELS="${MODELS:-$DEPMAP_DIR/Model.csv}"
COPY_NUMBER="${COPY_NUMBER:-$DEPMAP_DIR/OmicsCNGeneWGS.csv}"
FULL_RANKING="${FULL_RANKING:-results/expanded_26Q1/full/expanded_rses_onco.tsv}"
DEPMAP_RANKING="${DEPMAP_RANKING:-results/expanded_26Q1/depmap_only/expanded_rses_onco.tsv}"
ARTICLE_ROOT="${ARTICLE_ROOT:-article_outputs}"
LOG_DIR="${EXTENDED_LOG_DIR:-logs/extended_multiomics_v0120}"
MIN_GROUP_SIZE="${MIN_GROUP_SIZE:-3}"
LOSS_THRESHOLD="${LOSS_THRESHOLD:-0.30}"
STRICT_EXTENDED_SOURCES="${STRICT_EXTENDED_SOURCES:-1}"
STRICT_LAYOUT="${STRICT_LAYOUT:-1}"

mkdir -p "$LOG_DIR" "$EXTENDED_PROCESSED_DIR"

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

require_file() {
  [[ -s "$1" ]] || {
    echo "Mandatory input is absent or empty: $1" >&2
    exit 1
  }
}

preflight() {
  require_file "$EXTENDED_CONFIG"
  require_file "$MODELS"
  require_file "$COPY_NUMBER"
  require_file "$FULL_RANKING"
  [[ -d "$EXTENDED_DATA_DIR" ]] || {
    echo "Extended DepMap directory is absent: $EXTENDED_DATA_DIR" >&2
    exit 1
  }
  bash -n scripts/run_extended_multiomics_pipeline.sh
  python -m compileall -q src scripts tests
  echo "Extended data directory: $EXTENDED_DATA_DIR"
  echo "Full ranking: $FULL_RANKING"
}

analysis() {
  preflight
  local strict_flag=()
  if [[ "$STRICT_EXTENDED_SOURCES" == "1" ]]; then
    strict_flag=(--strict)
  fi

  log_stage "Validate all local extended-source schemas and duplicate content"
  run_logged "$LOG_DIR/00_validate_extended_input_schemas.log" \
    python -u scripts/validate_extended_input_schemas.py \
      --config "$EXTENDED_CONFIG" \
      --input-dir "$EXTENDED_DATA_DIR" \
      --output "$EXTENDED_PROCESSED_DIR/extended_input_schema_validation.tsv" \
      "${strict_flag[@]}"

  log_stage "Build standardized extended multi-omics evidence"
  run_logged "$LOG_DIR/01_build_extended_multiomics_evidence.log" \
    python -u scripts/build_extended_multiomics_evidence.py \
      --config "$EXTENDED_CONFIG" \
      --input-dir "$EXTENDED_DATA_DIR" \
      --models "$MODELS" \
      --copy-number "$COPY_NUMBER" \
      --ranking "$FULL_RANKING" \
      --output-dir "$EXTENDED_PROCESSED_DIR" \
      --min-group-size "$MIN_GROUP_SIZE" \
      --loss-threshold "$LOSS_THRESHOLD" \
      "${strict_flag[@]}"

  log_stage "Annotate gene- and platform-specific score eligibility"
  run_logged "$LOG_DIR/02_annotate_extended_layer_eligibility.log" \
    python -u scripts/annotate_extended_layer_eligibility.py \
      --config "$EXTENDED_CONFIG" \
      --input-dir "$EXTENDED_DATA_DIR" \
      --models "$MODELS" \
      --evidence "$EXTENDED_PROCESSED_DIR/extended_pair_evidence_by_cancer.tsv"

  log_stage "Recompute integrated TCGA-DepMap ranking as RSES-Onco v0.12.0"
  run_logged "$LOG_DIR/03_recompute_full_extended_rses.log" \
    python -u scripts/recompute_rses_with_extended_multiomics.py \
      --ranking "$FULL_RANKING" \
      --evidence "$EXTENDED_PROCESSED_DIR/extended_pair_evidence_by_cancer.tsv" \
      --output "$FULL_RANKING"

  if [[ -s "$DEPMAP_RANKING" ]]; then
    log_stage "Recompute DepMap-only ranking with the same causal multi-omics semantics"
    run_logged "$LOG_DIR/04_recompute_depmap_extended_rses.log" \
      python -u scripts/recompute_rses_with_extended_multiomics.py \
        --ranking "$DEPMAP_RANKING" \
        --evidence "$EXTENDED_PROCESSED_DIR/extended_pair_evidence_by_cancer.tsv" \
        --output "$DEPMAP_RANKING"
  fi

  log_stage "Validate extended multi-omics analysis integrity"
  local source_flag=()
  if [[ "$STRICT_EXTENDED_SOURCES" == "1" ]]; then
    source_flag=(--strict-sources)
  fi
  run_logged "$LOG_DIR/05_validate_extended_multiomics_analysis.log" \
    python -u scripts/validate_extended_multiomics_integrity.py \
      --processed-dir "$EXTENDED_PROCESSED_DIR" \
      --ranking "$FULL_RANKING" \
      --article-root "$ARTICLE_ROOT" \
      "${source_flag[@]}"
}

publication() {
  preflight
  require_file "$EXTENDED_PROCESSED_DIR/extended_pair_evidence_by_cancer.tsv"
  log_stage "Generate Supplementary Figures S71-S78"
  local layout_flag="--strict-layout"
  [[ "$STRICT_LAYOUT" == "1" ]] || layout_flag="--no-strict-layout"
  run_logged "$LOG_DIR/06_make_extended_multiomics_figures.log" \
    env MPLBACKEND=Agg python -u scripts/make_extended_multiomics_figures.py \
      --config config/extended_multiomics_asset.yaml \
      --ranking "$FULL_RANKING" \
      --processed-dir "$EXTENDED_PROCESSED_DIR" \
      --output-root "$ARTICLE_ROOT" \
      "$layout_flag"

  log_stage "Register Figures S71-S78 and Supplementary Tables S53-S64"
  run_logged "$LOG_DIR/07_register_extended_multiomics_assets.log" \
    python -u scripts/register_extended_multiomics_assets.py \
      --article-root "$ARTICLE_ROOT"

  log_stage "Recatalog exact source data for every registered figure"
  run_logged "$LOG_DIR/08_recatalog_figure_source_data.log" \
    python -u scripts/catalog_figure_source_data.py \
      --article-root "$ARTICLE_ROOT"

  log_stage "Validate extended multi-omics publication integration"
  local source_flag=()
  if [[ "$STRICT_EXTENDED_SOURCES" == "1" ]]; then
    source_flag=(--strict-sources)
  fi
  run_logged "$LOG_DIR/09_validate_extended_multiomics_publication.log" \
    python -u scripts/validate_extended_multiomics_integrity.py \
      --processed-dir "$EXTENDED_PROCESSED_DIR" \
      --ranking "$FULL_RANKING" \
      --article-root "$ARTICLE_ROOT" \
      --require-publication-assets \
      "${source_flag[@]}"
}

stage="${1:-all}"
case "$stage" in
  check-runtime)
    preflight
    ;;
  analysis)
    analysis
    ;;
  publication)
    publication
    ;;
  all)
    analysis
    publication
    ;;
  *)
    echo "Usage: bash scripts/run_extended_multiomics_pipeline.sh [check-runtime|analysis|publication|all]" >&2
    exit 2
    ;;
esac
