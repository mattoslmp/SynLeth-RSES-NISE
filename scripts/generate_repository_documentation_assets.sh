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
SUPPLEMENT_METHYLATION_MD="supplementary/Supplementary_Methylation_Methods_RSES_Onco_v0111.md"
SUPPLEMENT_DOCX="supplementary/Supplementary_Methods_RSES_Onco_v0110.docx"
SUPPLEMENT_PDF="supplementary/Supplementary_Methods_RSES_Onco_v0110.pdf"
SUPPLEMENT_METHYLATION_DOCX="supplementary/Supplementary_Methylation_Methods_RSES_Onco_v0111.docx"
SUPPLEMENT_METHYLATION_PDF="supplementary/Supplementary_Methylation_Methods_RSES_Onco_v0111.pdf"

MANUSCRIPT_MD="manuscript/RSES_Onco_intro_methods_draft_v0110.md"
MANUSCRIPT_METHYLATION_MD="manuscript/RSES_Onco_methylation_methods_addendum_v0111.md"
MANUSCRIPT_DOCX="manuscript/RSES_Onco_intro_methods_draft_v0110.docx"
MANUSCRIPT_PDF="manuscript/RSES_Onco_intro_methods_draft_v0110.pdf"
MANUSCRIPT_METHYLATION_DOCX="manuscript/RSES_Onco_methylation_methods_addendum_v0111.docx"
MANUSCRIPT_METHYLATION_PDF="manuscript/RSES_Onco_methylation_methods_addendum_v0111.pdf"

WORKFLOW_DOT="docs/figures/RSES_Onco_workflow_and_applications.dot"
WORKFLOW_SVG="docs/figures/RSES_Onco_workflow_and_applications.svg"
WORKFLOW_PNG="docs/figures/RSES_Onco_workflow_and_applications.png"

for source in \
  "$SUPPLEMENT_MD" \
  "$SUPPLEMENT_METHYLATION_MD" \
  "$MANUSCRIPT_MD" \
  "$MANUSCRIPT_METHYLATION_MD" \
  "$WORKFLOW_DOT"; do
  [[ -s "$source" ]] || {
    echo "Required documentation source is missing or empty: $source" >&2
    exit 1
  }
done

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
COMBINED_SUPPLEMENT="$TMP_DIR/Supplementary_Methods_RSES_Onco_v0111_combined.md"
COMBINED_MANUSCRIPT="$TMP_DIR/RSES_Onco_intro_methods_v0111_combined.md"

{
  cat "$SUPPLEMENT_MD"
  printf '\n\n---\n\n'
  cat "$SUPPLEMENT_METHYLATION_MD"
} > "$COMBINED_SUPPLEMENT"

{
  cat "$MANUSCRIPT_MD"
  printf '\n\n---\n\n'
  cat "$MANUSCRIPT_METHYLATION_MD"
} > "$COMBINED_MANUSCRIPT"

pandoc \
  "$COMBINED_SUPPLEMENT" \
  --from markdown+tex_math_dollars \
  --standalone \
  --metadata title="Supplementary Methods: RSES-Onco v0.11.1" \
  --output "$SUPPLEMENT_DOCX"

pandoc \
  "$SUPPLEMENT_METHYLATION_MD" \
  --from markdown+tex_math_dollars \
  --standalone \
  --metadata title="Supplementary Methylation Methods: RSES-Onco v0.11.1" \
  --output "$SUPPLEMENT_METHYLATION_DOCX"

pandoc \
  "$COMBINED_MANUSCRIPT" \
  --from markdown+tex_math_dollars \
  --standalone \
  --metadata title="RSES-Onco Introduction and Materials and Methods v0.11.1" \
  --output "$MANUSCRIPT_DOCX"

pandoc \
  "$MANUSCRIPT_METHYLATION_MD" \
  --from markdown+tex_math_dollars \
  --standalone \
  --metadata title="RSES-Onco methylation methods addendum v0.11.1" \
  --output "$MANUSCRIPT_METHYLATION_DOCX"

LO_PROFILE="file:///tmp/rses-onco-libreoffice-${UID}-$$"

libreoffice \
  -env:UserInstallation="$LO_PROFILE" \
  --headless \
  --convert-to pdf \
  --outdir supplementary \
  "$SUPPLEMENT_DOCX" \
  "$SUPPLEMENT_METHYLATION_DOCX" \
  > logs/documentation_assets/supplementary_pdf.log \
  2>&1

libreoffice \
  -env:UserInstallation="$LO_PROFILE" \
  --headless \
  --convert-to pdf \
  --outdir manuscript \
  "$MANUSCRIPT_DOCX" \
  "$MANUSCRIPT_METHYLATION_DOCX" \
  > logs/documentation_assets/manuscript_pdf.log \
  2>&1

dot -Tsvg "$WORKFLOW_DOT" -o "$WORKFLOW_SVG"
dot -Tpng "$WORKFLOW_DOT" -o "$WORKFLOW_PNG"

required_outputs=(
  "$SUPPLEMENT_DOCX"
  "$SUPPLEMENT_PDF"
  "$SUPPLEMENT_METHYLATION_DOCX"
  "$SUPPLEMENT_METHYLATION_PDF"
  "$MANUSCRIPT_DOCX"
  "$MANUSCRIPT_PDF"
  "$MANUSCRIPT_METHYLATION_DOCX"
  "$MANUSCRIPT_METHYLATION_PDF"
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
  Path("supplementary/Supplementary_Methylation_Methods_RSES_Onco_v0111.md"),
  Path("supplementary/Supplementary_Methods_RSES_Onco_v0110.docx"),
  Path("supplementary/Supplementary_Methods_RSES_Onco_v0110.pdf"),
  Path("supplementary/Supplementary_Methylation_Methods_RSES_Onco_v0111.docx"),
  Path("supplementary/Supplementary_Methylation_Methods_RSES_Onco_v0111.pdf"),
  Path("manuscript/RSES_Onco_intro_methods_draft_v0110.md"),
  Path("manuscript/RSES_Onco_methylation_methods_addendum_v0111.md"),
  Path("manuscript/RSES_Onco_intro_methods_draft_v0110.docx"),
  Path("manuscript/RSES_Onco_intro_methods_draft_v0110.pdf"),
  Path("manuscript/RSES_Onco_methylation_methods_addendum_v0111.docx"),
  Path("manuscript/RSES_Onco_methylation_methods_addendum_v0111.pdf"),
  Path("docs/figures/RSES_Onco_workflow_and_applications.svg"),
  Path("docs/figures/RSES_Onco_workflow_and_applications.png"),
]

for path in required:
  if not path.exists() or path.stat().st_size == 0:
    raise SystemExit(f"Missing generated documentation asset: {path}")

print("Documentation asset validation passed.")
PY
