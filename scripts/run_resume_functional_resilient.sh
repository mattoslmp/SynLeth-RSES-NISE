#!/usr/bin/env bash
# Resume after completed candidate/Ensembl/all-target stages using resilient
# STRING and DoRothEA acquisition, then generate every publication asset.
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
LOG_DIR="${LOG_DIR:-logs/expanded_26Q1_resilient}"

GENE_EFFECT="${GENE_EFFECT:-$DEPMAP_DIR/CRISPRGeneEffect.csv}"
COPY_NUMBER="${COPY_NUMBER:-$DEPMAP_DIR/OmicsCNGeneWGS.csv}"
MODELS="${MODELS:-$DEPMAP_DIR/Model.csv}"
EXPRESSION="${EXPRESSION:-$DEPMAP_DIR/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv}"
GDC_MANIFEST="${GDC_MANIFEST:-$GDC_DIR/gdc_gene_level_copy_number_manifest.json}"
CANDIDATES="${CANDIDATES:-$PROCESSED_DIR/expanded_candidate_universe.tsv}"
MEMBERS="${MEMBERS:-$PROCESSED_DIR/expanded_class_member_inventory.tsv}"
PARALOGS="${PARALOGS:-data/raw/ensembl/human_seed_paralogs.tsv}"
DISCOVERY="${DISCOVERY:-$DISCOVERY_RESULTS/all_target_dependency_screen.tsv}"
FUNCTIONAL_EVIDENCE="${FUNCTIONAL_EVIDENCE:-$PROCESSED_DIR/expanded_pair_functional_evidence.tsv}"
LOSS_THRESHOLD="${LOSS_THRESHOLD:-0.30}"
MIN_GROUP_SIZE="${MIN_GROUP_SIZE:-3}"

mkdir -p "$PROCESSED_DIR" "$DEPMAP_RESULTS" "$FULL_RESULTS" \
  "$DISCOVERY_RESULTS" "$LOG_DIR" article_outputs

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

require_nonempty() {
  local path="$1"
  if [[ ! -s "$path" ]]; then
    echo "Required completed-stage input is absent or empty: $path" >&2
    exit 1
  fi
}

log_stage "Validate prerequisites from completed stages"
for required in \
  "$GENE_EFFECT" \
  "$COPY_NUMBER" \
  "$MODELS" \
  "$EXPRESSION" \
  "$GDC_MANIFEST" \
  "$CANDIDATES" \
  "$MEMBERS" \
  "$PARALOGS" \
  "$DISCOVERY"; do
  require_nonempty "$required"
done

echo "Candidate universe: $CANDIDATES"
echo "All-target screen: $DISCOVERY"
echo "DepMap directory: $DEPMAP_DIR"
echo "GDC directory: $GDC_DIR"

log_stage "Acquire STRING, DoRothEA, HPA and UniProt evidence resiliently"
run_logged "$LOG_DIR/08_human_functional_evidence_resilient.log" \
  python -u scripts/download_human_functional_evidence_resilient.py \
    --candidates "$CANDIDATES" \
    --raw-dir data/raw/human_functional_evidence \
    --output "$FUNCTIONAL_EVIDENCE" \
    --string-required-score 700 \
    --string-limit 100 \
    --string-sleep 1.0 \
    --string-retries 7 \
    --string-map-chunk-size 200 \
    --strict-string-requests

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

log_stage "Run pharmacology, structures and complete publication assets"
set -o pipefail
RANKING="$FULL_RESULTS/expanded_rses_onco.tsv" \
CANDIDATES="$CANDIDATES" \
MEMBERS="$MEMBERS" \
DISCOVERY="$DISCOVERY" \
FUNCTIONAL_EVIDENCE="$FUNCTIONAL_EVIDENCE" \
GENE_EFFECT="$GENE_EFFECT" \
COPY_NUMBER="$COPY_NUMBER" \
MODELS="$MODELS" \
DEPMAP_DIR="$DEPMAP_DIR" \
LOSS_THRESHOLD="$LOSS_THRESHOLD" \
MIN_GROUP_SIZE="$MIN_GROUP_SIZE" \
bash scripts/run_publication_pipeline.sh all \
  2>&1 | tee "$LOG_DIR/14_publication_pipeline.log"
publication_status=${PIPESTATUS[0]}
if [[ $publication_status -ne 0 ]]; then
  echo "Publication pipeline failed with status $publication_status" >&2
  exit "$publication_status"
fi

log_stage "Run final tests and expanded-result checksums"
run_logged "$LOG_DIR/15_pytest.log" \
  python -m pytest -q -p no:cacheprovider

find "$RESULT_ROOT" article_outputs \
  -type f ! -name SHA256SUMS.txt -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "$FULL_RESULTS/SHA256SUMS.txt"

echo "Wrote $FULL_RESULTS/SHA256SUMS.txt"
log_stage "Resilient resume completed successfully"
