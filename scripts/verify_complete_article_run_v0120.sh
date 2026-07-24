#!/usr/bin/env bash
# Verify the complete RSES-Onco v0.12.0 extended multi-omics submission package.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

ARTICLE_ROOT="${ARTICLE_ROOT:-article_outputs}"
RESULT_ROOT="${RESULT_ROOT:-results/expanded_26Q1}"
GDC_DIR="${GDC_DIR:-data/raw/gdc}"
STRICT_EXTENDED_SOURCES="${STRICT_EXTENDED_SOURCES:-1}"

log_stage() {
  printf '\n[%s] %s\n' "$(date -Iseconds)" "$*"
}

log_stage "Verify no incomplete GDC download remains"
partial_count="$(
  find "$GDC_DIR" -type f -name '*.part' 2>/dev/null \
    | wc -l \
    | tr -d '[:space:]'
)"
echo "GDC .part files: $partial_count"
test "$partial_count" -eq 0

log_stage "Validate v0.12.0 multi-omics scientific integrity"
source_flag=()
if [[ "$STRICT_EXTENDED_SOURCES" == "1" ]]; then
  source_flag=(--strict-sources)
fi
python -u scripts/validate_extended_multiomics_integrity.py \
  --article-root "$ARTICLE_ROOT" \
  "${source_flag[@]}"

log_stage "Validate complete 86-figure and 68-table package"
MPLBACKEND=Agg \
python -u scripts/validate_publication_outputs_v0120.py \
  --article-root "$ARTICLE_ROOT" \
  "${source_flag[@]}"

log_stage "Validate rendered manuscript and supplementary documents"
python -u scripts/validate_publication_documents.py \
  --article-root "$ARTICLE_ROOT" \
  --document-dir "$ARTICLE_ROOT/documents" \
  --require-page-renders

log_stage "Verify publication-package SHA-256 checksums"
(
  cd "$ARTICLE_ROOT"
  sha256sum -c manifests/SHA256SUMS.txt
)

if [[ -s "$RESULT_ROOT/full/SHA256SUMS.txt" ]]; then
  log_stage "Verify expanded-result SHA-256 checksums"
  sha256sum -c "$RESULT_ROOT/full/SHA256SUMS.txt"
fi

log_stage "Run complete software test suite"
PYTHONDONTWRITEBYTECODE=1 \
MPLBACKEND=Agg \
python -m pytest -q -p no:cacheprovider

log_stage "Complete RSES-Onco v0.12.0 validation passed"
echo "Main figures: 8"
echo "Supplementary figures: 78"
echo "Figure exports: 258"
echo "Main tables: 4"
echo "Supplementary tables: 64"
echo "Manual inspection of every rendered figure and document page at 100% zoom remains mandatory."
