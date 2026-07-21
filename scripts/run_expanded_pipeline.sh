#!/usr/bin/env bash
# End-to-end expanded RSES-Onco analysis, pharmacology and publication workflow.
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
  "$DISCOVERY_RESULTS" "$LOG_DIR" supplementary

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
  local path="$1"
  local description="$2"
  if [[ ! -s "$path" ]]; then
    echo "Required resume input is absent or empty: $description: $path" >&2
    exit 1
  fi
}

validate_depmap() {
  log_stage "Validate DepMap files, checksums and cancer cohorts"
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
  log_stage "Retrieve Ensembl Compara paralogs"
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
  log_stage "Screen every measured DepMap CRISPR target"
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

  log_stage "Rebuild universe with supported discovery candidates"
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
  log_stage "Acquire STRING, DoRothEA, HPA and UniProt evidence"
  run_logged "$LOG_DIR/08_human_functional_evidence.log" \
    python -u scripts/download_human_functional_evidence.py \
      --candidates "$CANDIDATES" \
      --raw-dir data/raw/human_functional_evidence \
      --output "$FUNCTIONAL_EVIDENCE" \
      --strict-string-requests
}

run_depmap_expanded() {
  log_stage "Run expanded DepMap-only RSES-Onco"
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
  log_stage "Validate all GDC files from the reviewed manifest"
  run_logged "$LOG_DIR/10_validate_gdc.log" \
    python -u scripts/download_gdc.py \
      --validate-only \
      --manifest "$GDC_MANIFEST" \
      --output-dir "$GDC_DIR"
  log_stage "Aggregate ASCAT3 files"
  run_logged "$LOG_DIR/11_aggregate_gdc.log" \
    python -u scripts/aggregate_gdc_gene_cna.py \
      --raw-dir "$GDC_DIR" \
      --manifest "$GDC_MANIFEST" \
      --output-dir "$PROCESSED_DIR"
  log_stage "Validate combined cancer matrices"
  run_logged "$LOG_DIR/12_validate_gdc_matrices.log" \
    python -u scripts/validate_gdc_matrices.py \
      --output "$FULL_RESULTS/gdc_matrix_qc.tsv" \
      --event-output "$FULL_RESULTS/tcga_gene_event_summary.tsv"
}

run_full_expanded() {
  log_stage "Run integrated TCGA plus DepMap expanded RSES-Onco"
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

run_publication() {
  log_stage "Run pharmacology and complete publication-asset workflow"
  RANKING="$FULL_RESULTS/expanded_rses_onco.tsv" \
  CANDIDATES="$CANDIDATES" \
  MEMBERS="$MEMBERS" \
  DISCOVERY="$DISCOVERY_RESULTS/all_target_dependency_screen.tsv" \
  FUNCTIONAL_EVIDENCE="$FUNCTIONAL_EVIDENCE" \
  GENE_EFFECT="$GENE_EFFECT" \
  COPY_NUMBER="$COPY_NUMBER" \
  MODELS="$MODELS" \
  LOSS_THRESHOLD="$LOSS_THRESHOLD" \
  MIN_GROUP_SIZE="$MIN_GROUP_SIZE" \
  bash scripts/run_publication_pipeline.sh all \
    2>&1 | tee "$LOG_DIR/14_publication_pipeline.log"
  local status=${PIPESTATUS[0]}
  if [[ $status -ne 0 ]]; then
    return "$status"
  fi
}

finalize() {
  log_stage "Run final tests and expanded result checksums"
  run_logged "$LOG_DIR/15_pytest.log" \
    python -m pytest -q -p no:cacheprovider
  mkdir -p "$FULL_RESULTS"
  find "$RESULT_ROOT" article_outputs \
    -type f ! -name SHA256SUMS.txt -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "$FULL_RESULTS/SHA256SUMS.txt"
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

resume_functional() {
  require_file "$CANDIDATES" "expanded candidate universe"
  require_file "$MEMBERS" "class-member inventory"
  require_file "$PARALOGS" "Ensembl paralog catalogue"
  require_file \
    "$DISCOVERY_RESULTS/all_target_dependency_screen.tsv" \
    "all-target dependency screen"
  acquire_functional_evidence
  run_depmap_expanded
  validate_and_aggregate_gdc
  run_full_expanded
  run_publication
  finalize
}

after_download() {
  setup
  validate_and_aggregate_gdc
  run_full_expanded
  run_publication
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
  log_stage "Download or resume the reviewed GDC manifest"
  run_logged "$LOG_DIR/10b_download_gdc.log" \
    python -u scripts/download_gdc.py \
      --use-existing-manifest \
      --manifest "$GDC_MANIFEST" \
      --workflow ASCAT3 \
      --output-dir "$GDC_DIR" \
      --retries 3
  validate_and_aggregate_gdc
  run_full_expanded
  run_publication
  finalize
}

usage() {
  cat <<'EOF'
Usage: bash scripts/run_expanded_pipeline.sh STAGE

Stages:
  setup             Build all NISEs/paralogs, all-target discoveries and DepMap scores
  resume-functional Resume after a completed all-target discovery at functional evidence
  after-download    Continue after the current GDC download and generate every article asset
  all               Complete workflow including GDC acquisition
  publication       Rebuild pharmacology, main/supplementary figures and tables from existing scores
  discover-all-cn   Optional broad screen with every eligible copy-number loss gene
EOF
}

stage="${1:-}"
case "$stage" in
  setup) setup ;;
  resume-functional) resume_functional ;;
  after-download) after_download ;;
  all) all ;;
  publication) run_publication ;;
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
