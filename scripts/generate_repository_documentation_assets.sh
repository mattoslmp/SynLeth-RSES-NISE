#!/usr/bin/env bash

set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

require_command() {
  local command_name="$1"
  command -v "$command_name" >/dev/null 2>&1 || {
    echo "Required command is unavailable: $command_name" >&2
    exit 1
  }
}

require_command python
require_command pandoc
require_command libreoffice
require_command dot

python -u scripts/sync_methylation_documentation_v0111.py
python -u scripts/sync_circos_documentation_v0112.py

mkdir -p \
  supplementary \
  manuscript \
  docs/figures \
  logs/documentation_assets

SUPPLEMENT_MD="supplementary/Supplementary_Methods_RSES_Onco_v0110.md"
SUPPLEMENT_DOCX="supplementary/Supplementary_Methods_RSES_Onco_v0110.docx"
SUPPLEMENT_PDF="supplementary/Supplementary_Methods_RSES_Onco_v0110.pdf"

MANUSCRIPT_MD="manuscript/RSES_Onco_intro_methods_draft_v0110.md"
MANUSCRIPT_DOCX="manuscript/RSES_Onco_intro_methods_draft_v0110.docx"
MANUSCRIPT_PDF="manuscript/RSES_Onco_intro_methods_draft_v0110.pdf"

WORKFLOW_DOT="docs/figures/RSES_Onco_workflow_and_applications.dot"
WORKFLOW_SVG="docs/figures/RSES_Onco_workflow_and_applications.svg"
WORKFLOW_PNG="docs/figures/RSES_Onco_workflow_and_applications.png"

for source in "$SUPPLEMENT_MD" "$MANUSCRIPT_MD" "$WORKFLOW_DOT"; do
  [[ -s "$source" ]] || {
    echo "Required documentation source is missing or empty: $source" >&2
    exit 1
  }
done

pandoc \
  "$SUPPLEMENT_MD" \
  --from markdown+tex_math_dollars \
  --standalone \
  --metadata title="Supplementary Methods: RSES-Onco v0.11.2" \
  --output "$SUPPLEMENT_DOCX"

pandoc \
  "$MANUSCRIPT_MD" \
  --from markdown+tex_math_dollars \
  --standalone \
  --metadata title="RSES-Onco Introduction and Materials and Methods v0.11.2" \
  --output "$MANUSCRIPT_DOCX"

LO_PROFILE="file:///tmp/rses-onco-libreoffice-${UID}-$$"

libreoffice \
  -env:UserInstallation="$LO_PROFILE" \
  --headless \
  --convert-to pdf \
  --outdir supplementary \
  "$SUPPLEMENT_DOCX" \
  > logs/documentation_assets/supplementary_pdf.log \
  2>&1

libreoffice \
  -env:UserInstallation="$LO_PROFILE" \
  --headless \
  --convert-to pdf \
  --outdir manuscript \
  "$MANUSCRIPT_DOCX" \
  > logs/documentation_assets/manuscript_pdf.log \
  2>&1

dot -Tsvg "$WORKFLOW_DOT" -o "$WORKFLOW_SVG"
dot -Tpng "$WORKFLOW_DOT" -o "$WORKFLOW_PNG"

required_outputs=(
  "$SUPPLEMENT_DOCX"
  "$SUPPLEMENT_PDF"
  "$MANUSCRIPT_DOCX"
  "$MANUSCRIPT_PDF"
  "$WORKFLOW_SVG"
  "$WORKFLOW_PNG"
)

for output in "${required_outputs[@]}"; do
  [[ -s "$output" ]] || {
    echo "Documentation output is missing or empty: $output" >&2
    exit 1
  }
  ls -lh "$output"
done

python - <<'PY'
from pathlib import Path

required = [
  Path("docs/END_TO_END_ARTICLE_PROTOCOL.md"),
  Path("docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md"),
  Path("docs/METHYLATION_DATA_AND_SCORING_V0111.md"),
  Path("docs/GENOMIC_CIRCOS_WORKFLOW_V0112.md"),
  Path("supplementary/Supplementary_Methods_RSES_Onco_v0110.md"),
  Path("supplementary/Supplementary_Methods_RSES_Onco_v0110.docx"),
  Path("supplementary/Supplementary_Methods_RSES_Onco_v0110.pdf"),
  Path("manuscript/RSES_Onco_intro_methods_draft_v0110.md"),
  Path("manuscript/RSES_Onco_intro_methods_draft_v0110.docx"),
  Path("manuscript/RSES_Onco_intro_methods_draft_v0110.pdf"),
  Path("docs/figures/RSES_Onco_workflow_and_applications.svg"),
  Path("docs/figures/RSES_Onco_workflow_and_applications.png"),
]
for path in required:
  if not path.exists() or path.stat().st_size == 0:
    raise SystemExit(f"Missing generated documentation asset: {path}")

protocol = Path("docs/END_TO_END_ARTICLE_PROTOCOL.md").read_text(
  encoding="utf-8"
)
for token in (
  "RSES-Onco-expanded-v0.10.10",
  "promoter-methylation-context-v1",
  "Methylation (1kb upstream TSS)",
  "Figure S70",
  "Supplementary Tables S45-S52",
  "70 supplementary figures",
  "52 supplementary tables",
):
  if token not in protocol:
    raise SystemExit(f"Canonical protocol lacks required token: {token}")

print("Circos-aware documentation asset validation passed.")
PY
