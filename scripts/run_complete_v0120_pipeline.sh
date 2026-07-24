#!/usr/bin/env bash
# End-to-end RSES-Onco v0.12.0: WGCNA/regulation/methylation, extended multi-omics, figures, tables and documents.
set -Eeuo pipefail
set -o pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

EXTENDED_DATA_DIR="${EXTENDED_DATA_DIR:-dmap_data}"
DEPMAP_DIR="${DEPMAP_DIR:-data/raw/depmap}"
RESULT_ROOT="${RESULT_ROOT:-results/expanded_26Q1}"
ARTICLE_ROOT="${ARTICLE_ROOT:-article_outputs}"
LOG_DIR="${COMPLETE_V0120_LOG_DIR:-logs/complete_v0120}"
STRICT_LAYOUT="${STRICT_LAYOUT:-1}"
STRICT_EXTENDED_SOURCES="${STRICT_EXTENDED_SOURCES:-1}"

mkdir -p "$LOG_DIR"
RUN_STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_LOG="$LOG_DIR/complete_v0120_${RUN_STAMP}.log"
EXITCODE_FILE="$LOG_DIR/complete_v0120_${RUN_STAMP}.exitcode"

finish() {
  local status=$?
  trap - EXIT
  printf '%s\n' "$status" > "$EXITCODE_FILE"
  echo "Complete v0.12.0 exit code: $status"
  echo "Log: $RUN_LOG"
  echo "Exitcode: $EXITCODE_FILE"
  exit "$status"
}
trap finish EXIT
exec > >(tee "$RUN_LOG") 2>&1

export EXTENDED_DATA_DIR DEPMAP_DIR RESULT_ROOT ARTICLE_ROOT
export STRICT_LAYOUT STRICT_EXTENDED_SOURCES
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1

METHYLATION="${METHYLATION:-}"
if [[ -z "$METHYLATION" ]]; then
  for candidate in \
    "$EXTENDED_DATA_DIR/Methylation_1kb_upstream_TSS_subsetted_NAsdropped.csv" \
    "$EXTENDED_DATA_DIR/Methylation_(1kb_upstream_TSS)_subsetted_NAsdropped.csv" \
    "$DEPMAP_DIR/Methylation_1kb_upstream_TSS.csv" \
    "$DEPMAP_DIR/CCLE_RRBS_TSS1kb_20181022.txt.gz"; do
    if [[ -s "$candidate" ]]; then
      METHYLATION="$candidate"
      break
    fi
  done
fi
if [[ -z "$METHYLATION" || ! -s "$METHYLATION" ]]; then
  echo "Promoter methylation input was not found." >&2
  exit 1
fi
export METHYLATION

echo "RSES-Onco complete v0.12.0 run"
echo "Commit: $(git rev-parse HEAD)"
echo "Extended data: $EXTENDED_DATA_DIR"
echo "Methylation: $METHYLATION"

bash scripts/resume_wgcna_regulatory_pipeline.sh check-runtime

PUBLICATION_STAGE=all \
METHYLATION="$METHYLATION" \
bash scripts/resume_wgcna_regulatory_pipeline.sh resume-regulatory

bash scripts/run_extended_multiomics_pipeline.sh all

# Rebuild workbook and editable documents after registering S71-S78 and S53-S64.
bash scripts/run_publication_pipeline.sh workbook
bash scripts/run_publication_pipeline.sh documents
bash scripts/run_publication_pipeline.sh manifests

# Rebuild result checksums after the score was replaced by v0.12.0.
find "$RESULT_ROOT" "$ARTICLE_ROOT" \
  -type f ! -name SHA256SUMS.txt -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "$RESULT_ROOT/full/SHA256SUMS.txt"

bash scripts/verify_complete_article_run_v0120.sh

echo "RSES-Onco v0.12.0 completed successfully."
