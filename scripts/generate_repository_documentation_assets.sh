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

require_command pandoc
require_command libreoffice
require_command dot

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
  --metadata title="Supplementary Methods: RSES-Onco v0.11.0" \
  --output "$SUPPLEMENT_DOCX"

pandoc \
  "$MANUSCRIPT_MD" \
  --from markdown+tex_math_dollars \
  --standalone \
  --metadata title="RSES-Onco Introduction and Materials and Methods" \
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

print("Documentation asset validation passed.")
PY
