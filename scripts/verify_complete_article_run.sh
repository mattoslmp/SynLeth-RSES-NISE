#!/usr/bin/env bash
# Verify that the complete real-data, structural and publication workflow finished correctly.
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
  if [[ ! -f "$PIPELINE_EXITCODE_FILE" ]]; then
    echo "Missing pipeline exit-code file: $PIPELINE_EXITCODE_FILE" >&2
    exit 1
  fi
  pipeline_status="$(tr -d '[:space:]' < "$PIPELINE_EXITCODE_FILE")"
  if [[ "$pipeline_status" != "0" ]]; then
    echo "Pipeline exit code is $pipeline_status, expected 0" >&2
    exit 1
  fi
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
required = {
  "failed_seed_homology_queries": 0,
  "unresolved_target_identifier_count": 0,
  "complete": True,
}
for key, expected in required.items():
  observed = data.get(key)
  print(f"{key}: {observed}")
  if observed != expected:
    raise SystemExit(f"Ensembl completeness check failed: {key}={observed!r}, expected {expected!r}")
print(f"seed_gene_count: {data.get('seed_gene_count')}")
print(f"directed_paralog_count: {data.get('directed_paralog_count')}")
PY

log_stage "Verify that no GDC partial download remains"
partial_count="$(find "$GDC_DIR" -type f -name '*.part' 2>/dev/null | wc -l | tr -d '[:space:]')"
echo "GDC .part files: $partial_count"
test "$partial_count" -eq 0

log_stage "Validate the complete publication package"
MPLBACKEND=Agg \
python -u scripts/validate_publication_outputs.py \
  --article-root "$ARTICLE_ROOT"

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
import json
import sys
from pathlib import Path

import pandas as pd

root = Path(sys.argv[1])
figures = pd.read_csv(root / "manifests/figure_manifest.tsv", sep="\t")
tables = pd.read_csv(root / "manifests/table_manifest.tsv", sep="\t")

main_figures = int(figures["figure_id"].astype(str).str.match(r"^Figure_[1-8]$").sum())
supplementary_figures = int(figures["figure_id"].astype(str).str.match(r"^Figure_S(?:[1-9]|[12][0-9]|3[0-2])$").sum())
image_files = sum(
  1
  for path in (root / "figures").rglob("*")
  if path.is_file() and path.suffix.lower() in {".png", ".pdf", ".svg"}
)
structure_renders = sum(
  1 for path in (root / "structure_atlas/individual").rglob("*.png") if path.is_file()
)
main_tables = int(tables["category"].astype(str).eq("main").sum())
supplementary_tables = int(tables["category"].astype(str).eq("supplementary").sum())

summary = {
  "main_figures": main_figures,
  "supplementary_figures": supplementary_figures,
  "exported_figure_files": image_files,
  "main_tables": main_tables,
  "supplementary_tables": supplementary_tables,
  "individual_structure_renders": structure_renders,
}
for key, value in summary.items():
  print(f"{key}: {value}")

expected = {
  "main_figures": 8,
  "supplementary_figures": 32,
  "exported_figure_files": 120,
  "main_tables": 4,
  "supplementary_tables": 18,
}
for key, value in expected.items():
  if summary[key] != value:
    raise SystemExit(f"Final asset count failed: {key}={summary[key]}, expected {value}")
if structure_renders < 140:
  raise SystemExit(
    f"Final asset count failed: individual_structure_renders={structure_renders}, expected at least 140"
  )

provenance_path = root / "manifests/publication_provenance.json"
if provenance_path.exists():
  provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
  print(f"git_commit: {provenance.get('git_commit')}")
  print(f"git_status_porcelain: {provenance.get('git_status_porcelain')!r}")
PY

log_stage "Complete run validation passed"
echo "Automated validation is complete. Manual inspection of every figure at 100% zoom remains mandatory."
