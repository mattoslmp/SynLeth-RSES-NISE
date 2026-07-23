#!/usr/bin/env bash
# Verify the complete real-data, Circos, structural and publication workflow.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

ARTICLE_ROOT="${ARTICLE_ROOT:-article_outputs}"
RESULT_ROOT="${RESULT_ROOT:-results/expanded_26Q1}"
GDC_DIR="${GDC_DIR:-data/raw/gdc}"
ENSEMBL_METADATA="${ENSEMBL_METADATA:-data/raw/ensembl/ensembl_acquisition_metadata.json}"
PIPELINE_EXITCODE_FILE="${PIPELINE_EXITCODE_FILE:-}"

log_stage() {
  printf '\n[%s] %s\n' "$(date -Iseconds)" "$*"
}

if [[ -n "$PIPELINE_EXITCODE_FILE" ]]; then
  log_stage "Check recorded pipeline exit code"
  [[ -f "$PIPELINE_EXITCODE_FILE" ]] || {
    echo "Missing pipeline exit-code file: $PIPELINE_EXITCODE_FILE" >&2
    exit 1
  }
  pipeline_status="$(tr -d '[:space:]' < "$PIPELINE_EXITCODE_FILE")"
  [[ "$pipeline_status" == "0" ]] || {
    echo "Pipeline exit code is $pipeline_status, expected 0" >&2
    exit 1
  }
  echo "Pipeline exit code: 0"
fi

log_stage "Verify Ensembl acquisition completeness"
python - "$ENSEMBL_METADATA" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
if not path.exists():
  raise SystemExit(f"Missing Ensembl metadata: {path}")
data = json.loads(path.read_text(encoding="utf-8"))
for key, expected in {
  "failed_seed_homology_queries": 0,
  "unresolved_target_identifier_count": 0,
  "complete": True,
}.items():
  observed = data.get(key)
  print(f"{key}: {observed}")
  if observed != expected:
    raise SystemExit(
      f"Ensembl completeness failed: {key}={observed!r}, "
      f"expected {expected!r}"
    )
PY

log_stage "Verify that no GDC partial download remains"
partial_count="$(
  find "$GDC_DIR" -type f -name '*.part' 2>/dev/null \
    | wc -l \
    | tr -d '[:space:]'
)"
echo "GDC .part files: $partial_count"
test "$partial_count" -eq 0

log_stage "Verify the base scientific-integrity report from the 77-figure core"
test -s "$ARTICLE_ROOT/manifests/scientific_integrity_validation.json"

log_stage "Validate Figure S70 and Supplementary Tables S45-S52"
MPLBACKEND=Agg \
python -u scripts/validate_genomic_circos_integrity.py \
  --article-root "$ARTICLE_ROOT"

log_stage "Validate the complete 78-figure publication package"
MPLBACKEND=Agg \
python -u scripts/validate_publication_outputs.py \
  --article-root "$ARTICLE_ROOT"

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

log_stage "Verify expanded-result SHA-256 checksums"
sha256sum -c "$RESULT_ROOT/full/SHA256SUMS.txt"

log_stage "Run the complete software test suite"
PYTHONDONTWRITEBYTECODE=1 \
MPLBACKEND=Agg \
python -m pytest -q -p no:cacheprovider

log_stage "Summarize final publication assets"
python - "$ARTICLE_ROOT" <<'PY'
import sys
from pathlib import Path
import pandas as pd

root = Path(sys.argv[1])
figures = pd.read_csv(
  root / "manifests/figure_manifest.tsv",
  sep="\t",
  low_memory=False,
)
tables = pd.read_csv(
  root / "manifests/table_manifest.tsv",
  sep="\t",
  low_memory=False,
)
summary = {
  "main_figures": int(
    figures["category"].astype(str).eq("main").sum()
  ),
  "supplementary_figures": int(
    figures["category"].astype(str).eq("supplementary").sum()
  ),
  "exported_figure_files": sum(
    1
    for path in (root / "figures").rglob("*")
    if path.is_file()
    and path.suffix.lower() in {".png", ".pdf", ".svg"}
  ),
  "main_tables": int(
    tables["category"].astype(str).eq("main").sum()
  ),
  "supplementary_tables": int(
    tables["category"].astype(str).eq("supplementary").sum()
  ),
  "individual_structure_renders": sum(
    1
    for path in (root / "structure_atlas/individual").rglob("*.png")
    if path.is_file()
  ),
}
for key, value in summary.items():
  print(f"{key}: {value}")
expected = {
  "main_figures": 8,
  "supplementary_figures": 70,
  "exported_figure_files": 234,
  "main_tables": 4,
  "supplementary_tables": 52,
}
for key, value in expected.items():
  if summary[key] != value:
    raise SystemExit(
      f"Final asset count failed: {key}={summary[key]}, "
      f"expected {value}"
    )
if summary["individual_structure_renders"] < 140:
  raise SystemExit("Expected at least 140 individual structure renders")
required = [
  root / "tables/supplementary/Table_S45_genomic_circos_gene_coordinates.tsv",
  root / "tables/supplementary/Table_S46_genomic_circos_pair_links.tsv",
  root / "tables/supplementary/Table_S47_genomic_circos_ring_values.tsv",
  root / "tables/supplementary/Table_S48_genomic_circos_track_definitions.tsv",
  root / "tables/supplementary/Table_S49_genomic_circos_expression_summary.tsv",
  root / "tables/supplementary/Table_S50_genomic_circos_expression_model_values.tsv",
  root / "tables/supplementary/Table_S51_pipeline_script_catalog.tsv",
  root / "tables/supplementary/Table_S52_genomic_circos_source_provenance.tsv",
  root / "tables/figure_data/supplementary/Figure_S70_source_data.tsv",
  root / "manifests/genomic_circos_integrity_validation.json",
  Path("docs/SCRIPT_CATALOG.md"),
  Path("docs/script_manifest.tsv"),
]
for path in required:
  if not path.exists() or path.stat().st_size == 0:
    raise SystemExit(f"Missing or empty required asset: {path}")
PY

log_stage "Complete run validation passed"
echo "Automated validation is complete. Manual inspection of every rendered page and figure at 100% zoom remains mandatory."
