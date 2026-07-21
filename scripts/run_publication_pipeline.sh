#!/usr/bin/env bash
# Generate pharmacology priorities and every article asset after expanded scoring.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

RANKING="${RANKING:-results/expanded_26Q1/full/expanded_rses_onco.tsv}"
CANDIDATES="${CANDIDATES:-data/processed/expanded_candidate_universe.tsv}"
MEMBERS="${MEMBERS:-data/processed/expanded_class_member_inventory.tsv}"
DISCOVERY="${DISCOVERY:-results/expanded_26Q1/discovery/all_target_dependency_screen.tsv}"
FUNCTIONAL_EVIDENCE="${FUNCTIONAL_EVIDENCE:-data/processed/expanded_pair_functional_evidence.tsv}"
PHARMACOLOGY_DATA="${PHARMACOLOGY_DATA:-data/processed/pharmacology}"
PHARMACOLOGY_RESULTS="${PHARMACOLOGY_RESULTS:-results/expanded_26Q1/pharmacology}"
ARTICLE_ROOT="${ARTICLE_ROOT:-article_outputs}"
LOG_DIR="${LOG_DIR:-logs/publication_26Q1}"
DEPMAP_DIR="${DEPMAP_DIR:-data/raw/depmap}"
GENE_EFFECT="${GENE_EFFECT:-$DEPMAP_DIR/CRISPRGeneEffect.csv}"
COPY_NUMBER="${COPY_NUMBER:-$DEPMAP_DIR/OmicsCNGeneWGS.csv}"
MODELS="${MODELS:-$DEPMAP_DIR/Model.csv}"
LOSS_THRESHOLD="${LOSS_THRESHOLD:-0.30}"
MIN_GROUP_SIZE="${MIN_GROUP_SIZE:-3}"
PHARMACOLOGY_MINIMUM_SCORE="${PHARMACOLOGY_MINIMUM_SCORE:-0.15}"
PHARMACOLOGY_MAX_TARGETS="${PHARMACOLOGY_MAX_TARGETS:-500}"
STRICT_LAYOUT="${STRICT_LAYOUT:-1}"

mkdir -p "$PHARMACOLOGY_DATA" "$PHARMACOLOGY_RESULTS" "$ARTICLE_ROOT" "$LOG_DIR"

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

require_ranking() {
  if [[ ! -s "$RANKING" ]]; then
    echo "Expanded ranking is absent: $RANKING" >&2
    echo "Run the expanded TCGA plus DepMap stage first." >&2
    exit 1
  fi
}

acquire_pharmacology() {
  require_ranking
  log_stage "Acquire Open Targets, ChEMBL, DGIdb, MyChem, Pharos and CIViC evidence"
  run_logged "$LOG_DIR/01_acquire_pharmacology.log" \
    python -u scripts/acquire_pharmacology_evidence.py \
      --ranking "$RANKING" \
      --output-dir "$PHARMACOLOGY_DATA" \
      --cache-dir data/raw/pharmacology/api_cache \
      --minimum-vulnerability-score "$PHARMACOLOGY_MINIMUM_SCORE" \
      --max-targets "$PHARMACOLOGY_MAX_TARGETS"
}

standardize_sensitivity() {
  log_stage "Standardize any locally available PRISM, GDSC and CTRP releases"
  run_logged "$LOG_DIR/02_standardize_drug_sensitivity.log" \
    python -u scripts/standardize_drug_sensitivity.py \
      --config config/drug_sensitivity_sources.yaml \
      --models "$MODELS" \
      --output "$PHARMACOLOGY_DATA/drug_sensitivity_long.tsv" \
      --status-output "$PHARMACOLOGY_DATA/drug_sensitivity_source_status.tsv"
}

analyze_sensitivity() {
  log_stage "Test biomarker-matched drug-response selectivity"
  run_logged "$LOG_DIR/03_analyze_drug_response_selectivity.log" \
    python -u scripts/analyze_drug_response_selectivity.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --ranking "$RANKING" \
      --pharmacology-evidence "$PHARMACOLOGY_DATA/pharmacology_evidence_long.tsv" \
      --sensitivity "$PHARMACOLOGY_DATA/drug_sensitivity_long.tsv" \
      --loss-threshold "$LOSS_THRESHOLD" \
      --min-group-size "$MIN_GROUP_SIZE" \
      --output "$PHARMACOLOGY_DATA/drug_response_selectivity.tsv" \
      --evidence-output "$PHARMACOLOGY_DATA/pharmacology_sensitivity_evidence.tsv"
}

prioritize_pharmacology() {
  require_ranking
  log_stage "Build coverage-aware experimental target-drug priorities"
  run_logged "$LOG_DIR/04_prioritize_pharmacology.log" \
    python -u scripts/prioritize_pharmacology.py \
      --ranking "$RANKING" \
      --evidence "$PHARMACOLOGY_DATA/pharmacology_evidence_long.tsv" \
      --sensitivity "$PHARMACOLOGY_DATA/drug_response_selectivity.tsv" \
      --output-dir "$PHARMACOLOGY_RESULTS"
}

export_tables() {
  require_ranking
  log_stage "Export all main and supplementary article tables"
  run_logged "$LOG_DIR/05_export_article_tables.log" \
    python -u scripts/export_article_tables.py \
      --ranking "$RANKING" \
      --candidates "$CANDIDATES" \
      --members "$MEMBERS" \
      --functional-evidence "$FUNCTIONAL_EVIDENCE" \
      --discovery "$DISCOVERY" \
      --pharmacology-evidence "$PHARMACOLOGY_DATA/pharmacology_evidence_long.tsv" \
      --pharmacology-ranking "$PHARMACOLOGY_RESULTS/pharmacology_ranked_hypotheses.tsv" \
      --drug-sensitivity "$PHARMACOLOGY_DATA/drug_response_selectivity.tsv" \
      --pharmacology-source-status "$PHARMACOLOGY_DATA/pharmacology_source_status.tsv" \
      --pharmacology-source-coverage "$PHARMACOLOGY_RESULTS/pharmacology_source_coverage.tsv" \
      --output-root "$ARTICLE_ROOT"
}

generate_figures() {
  require_ranking
  log_stage "Generate every main and supplementary figure from scripts"
  local layout_flag="--strict-layout"
  if [[ "$STRICT_LAYOUT" != "1" ]]; then
    layout_flag="--no-strict-layout"
  fi
  run_logged "$LOG_DIR/06_make_all_article_figures.log" \
    python -u scripts/make_all_article_figures.py \
      --ranking "$RANKING" \
      --candidates "$CANDIDATES" \
      --discovery "$DISCOVERY" \
      --pharmacology "$PHARMACOLOGY_RESULTS/pharmacology_ranked_hypotheses.tsv" \
      --output-root "$ARTICLE_ROOT" \
      "$layout_flag"
}

build_workbook() {
  log_stage "Build organized article and supplementary workbook"
  run_logged "$LOG_DIR/07_build_article_workbook.log" \
    python -u scripts/build_article_workbook.py \
      --article-root "$ARTICLE_ROOT" \
      --output "$ARTICLE_ROOT/workbooks/RSES_Onco_Article_Tables_and_Evidence.xlsx"
}

build_manifests() {
  log_stage "Build publication inventory, provenance and checksums"
  run_logged "$LOG_DIR/08_build_publication_manifest.log" \
    python -u scripts/build_publication_manifest.py \
      --article-root "$ARTICLE_ROOT" \
      --input "$RANKING" \
      --input "$CANDIDATES" \
      --input "$DISCOVERY" \
      --input "$FUNCTIONAL_EVIDENCE" \
      --input "$PHARMACOLOGY_DATA/pharmacology_evidence_long.tsv" \
      --input "$PHARMACOLOGY_RESULTS/pharmacology_ranked_hypotheses.tsv"
}

validate_outputs() {
  log_stage "Validate publication package and software tests"
  run_logged "$LOG_DIR/09_validate_publication_outputs.log" \
    python -u scripts/validate_publication_outputs.py \
      --article-root "$ARTICLE_ROOT"
  run_logged "$LOG_DIR/10_pytest.log" \
    python -m pytest -q -p no:cacheprovider
}

all() {
  acquire_pharmacology
  standardize_sensitivity
  analyze_sensitivity
  prioritize_pharmacology
  export_tables
  generate_figures
  build_workbook
  build_manifests
  validate_outputs
}

assets_only() {
  prioritize_pharmacology
  export_tables
  generate_figures
  build_workbook
  build_manifests
  validate_outputs
}

usage() {
  cat <<'EOF'
Usage: bash scripts/run_publication_pipeline.sh STAGE

Stages:
  acquire-pharmacology    Query/cached pharmacology evidence APIs
  standardize-sensitivity Standardize local PRISM, GDSC and CTRP releases
  analyze-sensitivity     Test biomarker-matched drug response
  prioritize              Build target-drug actionability priorities
  tables                  Export main and supplementary tables
  figures                 Generate all 7 main and 14 supplementary figures
  workbook                Build organized Excel workbook
  manifests               Build file inventory, provenance and SHA-256 checksums
  validate                Validate figure triplets, layout audits, tables and tests
  assets-only             Rebuild priorities and all publication assets from cached data
  all                     Run the complete pharmacology and publication workflow

Optional sensitivity files are discovered under:
  data/raw/pharmacology/prism/
  data/raw/pharmacology/gdsc/
  data/raw/pharmacology/ctrp/

If no local drug-sensitivity releases are present, the pipeline records them as
not found and continues with API-based target/compound evidence. Missing evidence
remains missing and lowers coverage.
EOF
}

stage="${1:-}"
case "$stage" in
  acquire-pharmacology) acquire_pharmacology ;;
  standardize-sensitivity) standardize_sensitivity ;;
  analyze-sensitivity) analyze_sensitivity ;;
  prioritize) prioritize_pharmacology ;;
  tables) export_tables ;;
  figures) generate_figures ;;
  workbook) build_workbook ;;
  manifests) build_manifests ;;
  validate) validate_outputs ;;
  assets-only) assets_only ;;
  all) all ;;
  *) usage; exit 2 ;;
esac
