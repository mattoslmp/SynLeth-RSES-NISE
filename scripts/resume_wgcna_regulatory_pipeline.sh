#!/usr/bin/env bash
# Resume from cached STRING/DoRothEA/HPA/UniProt evidence and rebuild RSES-Onco
# with signed WGCNA, TF-expression consistency and promoter motif support.
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
LOG_DIR="${LOG_DIR:-logs/expanded_26Q1_wgcna}"

GENE_EFFECT="${GENE_EFFECT:-$DEPMAP_DIR/CRISPRGeneEffect.csv}"
COPY_NUMBER="${COPY_NUMBER:-$DEPMAP_DIR/OmicsCNGeneWGS.csv}"
MODELS="${MODELS:-$DEPMAP_DIR/Model.csv}"
EXPRESSION="${EXPRESSION:-$DEPMAP_DIR/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv}"
CANDIDATES="${CANDIDATES:-$PROCESSED_DIR/expanded_candidate_universe.tsv}"
MEMBERS="${MEMBERS:-$PROCESSED_DIR/expanded_class_member_inventory.tsv}"
FUNCTIONAL_EVIDENCE="${FUNCTIONAL_EVIDENCE:-$PROCESSED_DIR/expanded_pair_functional_evidence.tsv}"
BASE_FUNCTIONAL_EVIDENCE="${BASE_FUNCTIONAL_EVIDENCE:-$PROCESSED_DIR/expanded_pair_functional_evidence_pre_wgcna.tsv}"
CANCER_SPECIFIC_EVIDENCE="${CANCER_SPECIFIC_EVIDENCE:-$PROCESSED_DIR/regulatory/expanded_pair_functional_evidence_by_cancer.tsv}"
DOROTHEA="${DOROTHEA:-data/raw/human_functional_evidence/omnipath_dorothea.tsv}"
PROMOTERS="${PROMOTERS:-data/raw/regulatory/ensembl_promoters.tsv}"
PROMOTER_FASTA="${PROMOTER_FASTA:-data/raw/regulatory/ensembl_promoters.fa}"
JASPAR_MEME="${JASPAR_MEME:-data/raw/regulatory/JASPAR2026_CORE_vertebrates_non-redundant.meme}"
PROMOTER_MOTIFS="${PROMOTER_MOTIFS:-data/processed/regulatory/jaspar_promoter_tf_summary.tsv}"
LOSS_THRESHOLD="${LOSS_THRESHOLD:-0.30}"
MIN_GROUP_SIZE="${MIN_GROUP_SIZE:-3}"

mkdir -p "$LOG_DIR" "$DEPMAP_RESULTS" "$FULL_RESULTS" \
  "$PROCESSED_DIR/regulatory"

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
  if [[ ! -s "$path" ]]; then
    echo "Mandatory input is absent or empty: $path" >&2
    exit 1
  fi
}

check_runtime() {
  command -v Rscript >/dev/null 2>&1 || {
    echo "Rscript is missing from the rses-onco environment." >&2
    return 1
  }
  Rscript -e "stopifnot(requireNamespace('WGCNA', quietly=TRUE)); cat(as.character(packageVersion('WGCNA')), '\n')"
  command -v fimo >/dev/null 2>&1 || {
    echo "FIMO from the MEME suite is missing." >&2
    return 1
  }
  bash -n scripts/resume_wgcna_regulatory_pipeline.sh
  python -m compileall -q src scripts tests
  echo "WGCNA/promoter regulatory runtime is ready."
}

build_regulatory_layer() {
  for path in "$GENE_EFFECT" "$COPY_NUMBER" "$MODELS" "$EXPRESSION" \
    "$CANDIDATES" "$MEMBERS" "$FUNCTIONAL_EVIDENCE" "$DOROTHEA"; do
    require_file "$path"
  done

  if [[ ! -s "$BASE_FUNCTIONAL_EVIDENCE" ]]; then
    cp -f "$FUNCTIONAL_EVIDENCE" "$BASE_FUNCTIONAL_EVIDENCE"
    echo "Preserved pre-WGCNA functional evidence: $BASE_FUNCTIONAL_EVIDENCE"
  fi

  log_stage "Acquire canonical Ensembl promoters and promoter sequences"
  run_logged "$LOG_DIR/01_ensembl_promoters.log" \
    python -u scripts/download_ensembl_promoters.py \
      --candidates "$CANDIDATES" \
      --output "$PROMOTERS" \
      --fasta "$PROMOTER_FASTA"

  log_stage "Acquire official JASPAR 2026 CORE vertebrate motifs"
  run_logged "$LOG_DIR/02_jaspar_core.log" \
    python -u scripts/download_jaspar_core_vertebrates.py \
      --output "$JASPAR_MEME"

  log_stage "Scan canonical promoters with JASPAR motifs using FIMO"
  run_logged "$LOG_DIR/03_fimo_promoter_scan.log" \
    python -u scripts/scan_promoter_motifs.py \
      --motifs "$JASPAR_MEME" \
      --promoters "$PROMOTER_FASTA" \
      --summary-output "$PROMOTER_MOTIFS"

  log_stage "Build cancer-specific signed WGCNA and promoter-aware TF evidence"
  run_logged "$LOG_DIR/04_wgcna_regulatory_layer.log" \
    python -u scripts/build_wgcna_regulatory_layer.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --expression "$EXPRESSION" \
      --candidates "$CANDIDATES" \
      --functional-evidence "$BASE_FUNCTIONAL_EVIDENCE" \
      --dorothea "$DOROTHEA" \
      --promoter-motifs "$PROMOTER_MOTIFS" \
      --output "$CANCER_SPECIFIC_EVIDENCE"

  log_stage "Build consensus functional prior without counting cancers as independent evidence"
  run_logged "$LOG_DIR/05_consensus_regulatory_layer.log" \
    python -u scripts/aggregate_wgcna_regulatory_layer.py \
      --base "$BASE_FUNCTIONAL_EVIDENCE" \
      --cancer-specific "$CANCER_SPECIFIC_EVIDENCE" \
      --output "$FUNCTIONAL_EVIDENCE"
}

score_depmap() {
  log_stage "Recalculate DepMap-only RSES-Onco and inject WGCNA/regulatory sublayers"
  run_logged "$LOG_DIR/06_depmap_score_base.log" \
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
  run_logged "$LOG_DIR/07_depmap_score_wgcna.log" \
    python -u scripts/recompute_rses_with_wgcna_regulatory.py \
      --ranking "$DEPMAP_RESULTS/expanded_rses_onco.tsv" \
      --functional-evidence "$FUNCTIONAL_EVIDENCE" \
      --output "$DEPMAP_RESULTS/expanded_rses_onco.tsv"
}

score_full() {
  for path in \
    "$PROCESSED_DIR/TCGA_COLON_homdel_discrete.tsv" \
    "$PROCESSED_DIR/TCGA_STOMACH_homdel_discrete.tsv" \
    "$PROCESSED_DIR/TCGA_LUNG_homdel_discrete.tsv"; do
    require_file "$path"
  done
  log_stage "Recalculate integrated TCGA-DepMap RSES-Onco"
  run_logged "$LOG_DIR/08_full_score_base.log" \
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
  run_logged "$LOG_DIR/09_full_score_wgcna.log" \
    python -u scripts/recompute_rses_with_wgcna_regulatory.py \
      --ranking "$FULL_RESULTS/expanded_rses_onco.tsv" \
      --functional-evidence "$FUNCTIONAL_EVIDENCE" \
      --output "$FULL_RESULTS/expanded_rses_onco.tsv"
}

publication() {
  require_file "$DISCOVERY_RESULTS/all_target_dependency_screen.tsv"
  log_stage "Regenerate pharmacology, evidence audits, figures, tables and manifests"
  RANKING="$FULL_RESULTS/expanded_rses_onco.tsv" \
  CANDIDATES="$CANDIDATES" \
  MEMBERS="$MEMBERS" \
  DISCOVERY="$DISCOVERY_RESULTS/all_target_dependency_screen.tsv" \
  FUNCTIONAL_EVIDENCE="$FUNCTIONAL_EVIDENCE" \
  GENE_EFFECT="$GENE_EFFECT" \
  COPY_NUMBER="$COPY_NUMBER" \
  MODELS="$MODELS" \
  EXPRESSION="$EXPRESSION" \
  LOSS_THRESHOLD="$LOSS_THRESHOLD" \
  MIN_GROUP_SIZE="$MIN_GROUP_SIZE" \
  bash scripts/run_publication_pipeline.sh all \
    2>&1 | tee "$LOG_DIR/10_publication_pipeline.log"
  local status=${PIPESTATUS[0]}
  [[ $status -eq 0 ]]
}

finalize() {
  log_stage "Run complete tests and rebuild expanded SHA-256 manifest"
  run_logged "$LOG_DIR/11_pytest.log" \
    python -m pytest -q -p no:cacheprovider
  find "$RESULT_ROOT" article_outputs \
    -type f ! -name SHA256SUMS.txt -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "$FULL_RESULTS/SHA256SUMS.txt"
  echo "Wrote $FULL_RESULTS/SHA256SUMS.txt"
}

stage="${1:-all}"
case "$stage" in
  --check-runtime|check-runtime)
    check_runtime
    ;;
  regulatory-only)
    check_runtime
    build_regulatory_layer
    ;;
  all|resume-regulatory)
    check_runtime
    build_regulatory_layer
    score_depmap
    score_full
    publication
    finalize
    ;;
  *)
    echo "Usage: bash scripts/resume_wgcna_regulatory_pipeline.sh [check-runtime|regulatory-only|resume-regulatory]" >&2
    exit 2
    ;;
esac
