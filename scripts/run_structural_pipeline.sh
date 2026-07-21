#!/usr/bin/env bash
# Download AlphaFold DB structures, collect known residues, render PyMOL images,
# and generate all structural article figures.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PROTEINS="${PROTEINS:-data/curated/human_nise_bonafide_2017.tsv}"
RANKING="${RANKING:-results/expanded_26Q1/full/expanded_rses_onco.tsv}"
STRUCTURE_MANIFEST="${STRUCTURE_MANIFEST:-data/processed/structures/alphafold_structure_manifest.tsv}"
ANNOTATIONS="${ANNOTATIONS:-data/processed/structures/nise_structural_residue_annotations.tsv}"
COVERAGE="${COVERAGE:-data/processed/structures/nise_structural_annotation_coverage.tsv}"
RENDER_MANIFEST="${RENDER_MANIFEST:-data/processed/structures/nise_structure_render_manifest.tsv}"
ARTICLE_ROOT="${ARTICLE_ROOT:-article_outputs}"
LOG_DIR="${LOG_DIR:-logs/structures_26Q1}"
PYMOL_EXECUTABLE="${PYMOL_EXECUTABLE:-pymol}"
STRICT_LAYOUT="${STRICT_LAYOUT:-1}"
CURATED_STRUCTURAL_RESIDUES="${CURATED_STRUCTURAL_RESIDUES:-data/curated/nise_exact_drug_contact_residues.tsv}"

mkdir -p data/raw/structures data/processed/structures \
  "$ARTICLE_ROOT/structure_atlas/individual" "$LOG_DIR"

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

download_structures() {
  run_logged "$LOG_DIR/01_download_alphafold.log" \
    python -u scripts/download_alphafold_nise_structures.py \
      --proteins "$PROTEINS" \
      --output-dir data/raw/structures/alphafold \
      --manifest "$STRUCTURE_MANIFEST"
}

collect_annotations() {
  local command=(
    python -u scripts/collect_nise_structural_annotations.py
    --proteins "$PROTEINS"
    --output "$ANNOTATIONS"
    --coverage-output "$COVERAGE"
    --cache-dir data/raw/structures/annotation_cache
  )
  if [[ -s "$CURATED_STRUCTURAL_RESIDUES" ]]; then
    command+=(--curated-residues "$CURATED_STRUCTURAL_RESIDUES")
  fi
  run_logged "$LOG_DIR/02_collect_structural_annotations.log" "${command[@]}"
}

render_structures() {
  run_logged "$LOG_DIR/03_render_structures.log" \
    python -u scripts/render_nise_structures.py \
      --structure-manifest "$STRUCTURE_MANIFEST" \
      --annotations "$ANNOTATIONS" \
      --output-dir "$ARTICLE_ROOT/structure_atlas/individual" \
      --render-manifest "$RENDER_MANIFEST" \
      --pymol "$PYMOL_EXECUTABLE"
}

generate_figures() {
  local layout_flag="--strict-layout"
  if [[ "$STRICT_LAYOUT" != "1" ]]; then
    layout_flag="--no-strict-layout"
  fi
  run_logged "$LOG_DIR/04_make_structural_figures.log" \
    python -u scripts/make_nise_structure_figures.py \
      --ranking "$RANKING" \
      --proteins "$PROTEINS" \
      --structure-manifest "$STRUCTURE_MANIFEST" \
      --render-manifest "$RENDER_MANIFEST" \
      --annotations "$ANNOTATIONS" \
      --coverage "$COVERAGE" \
      --output-root "$ARTICLE_ROOT" \
      "$layout_flag"
}

all() {
  download_structures
  collect_annotations
  render_structures
  generate_figures
}

usage() {
  cat <<'EOF'
Usage: bash scripts/run_structural_pipeline.sh STAGE

Stages:
  download      Download all curated human NISE AlphaFold DB models
  annotations   Collect exact UniProt-numbered M-CSA/UniProt/PDBe residues
  render        Produce 600-dpi PyMOL whole-structure and site images
  figures       Generate Figure 8 and Figures S15-S32
  all           Run all stages

Install the rendering engine with:
  conda install -c conda-forge pymol-open-source pillow

Optional exact drug-contact mappings:
  cp data/curated/nise_exact_drug_contact_residues.template.tsv \
     data/curated/nise_exact_drug_contact_residues.tsv
  # Fill only exact UniProt-numbered residues with source and mapping status.
EOF
}

case "${1:-}" in
  download) download_structures ;;
  annotations) collect_annotations ;;
  render) render_structures ;;
  figures) generate_figures ;;
  all) all ;;
  *) usage; exit 2 ;;
esac
