#!/usr/bin/env bash
# Generate pharmacology priorities, evidence audits, robustness analyses, structures and every publication asset.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

RANKING="${RANKING:-results/expanded_26Q1/full/expanded_rses_onco.tsv}"
CANDIDATES="${CANDIDATES:-data/processed/expanded_candidate_universe.tsv}"
MEMBERS="${MEMBERS:-data/processed/expanded_class_member_inventory.tsv}"
DISCOVERY="${DISCOVERY:-results/expanded_26Q1/discovery/all_target_dependency_screen.tsv}"
FUNCTIONAL_EVIDENCE="${FUNCTIONAL_EVIDENCE:-data/processed/expanded_pair_functional_evidence.tsv}"
PHARMACOLOGY_DATA="${PHARMACOLOGY_DATA:-data/processed/pharmacology}"
PHARMACOLOGY_RESULTS="${PHARMACOLOGY_RESULTS:-results/expanded_26Q1/pharmacology}"
ARTICLE_ROOT="${ARTICLE_ROOT:-article_outputs}"
LOG_DIR="${LOG_DIR:-logs/publication_26Q1}"
DEPMAP_DIR="${DEPMAP_DIR:-data/raw/depmap}"
GENE_EFFECT="${GENE_EFFECT:-$DEPMAP_DIR/CRISPRGeneEffect.csv}"
COPY_NUMBER="${COPY_NUMBER:-$DEPMAP_DIR/OmicsCNGeneWGS.csv}"
MODELS="${MODELS:-$DEPMAP_DIR/Model.csv}"
LOSS_THRESHOLD="${LOSS_THRESHOLD:-0.30}"
MIN_GROUP_SIZE="${MIN_GROUP_SIZE:-3}"
PHARMACOLOGY_MINIMUM_SCORE="${PHARMACOLOGY_MINIMUM_SCORE:-0.15}"
PHARMACOLOGY_MAX_TARGETS="${PHARMACOLOGY_MAX_TARGETS:-500}"
STRICT_LAYOUT="${STRICT_LAYOUT:-1}"
PYMOL_EXECUTABLE="${PYMOL_EXECUTABLE:-pymol}"

PROTEINS="${PROTEINS:-data/curated/human_nise_bonafide_2017.tsv}"
STRUCTURE_MANIFEST="${STRUCTURE_MANIFEST:-data/processed/structures/alphafold_structure_manifest.tsv}"
STRUCTURAL_ANNOTATIONS="${STRUCTURAL_ANNOTATIONS:-data/processed/structures/nise_structural_residue_annotations.tsv}"
STRUCTURAL_COVERAGE="${STRUCTURAL_COVERAGE:-data/processed/structures/nise_structural_annotation_coverage.tsv}"
STRUCTURE_RENDER_MANIFEST="${STRUCTURE_RENDER_MANIFEST:-data/processed/structures/nise_structure_render_manifest.tsv}"
STRUCTURE_LOG_DIR="${STRUCTURE_LOG_DIR:-logs/structures_26Q1}"
RUN_MARKER="${RUN_MARKER:-$LOG_DIR/assets_only_run.marker}"

mkdir -p "$PHARMACOLOGY_DATA" "$PHARMACOLOGY_RESULTS" "$ARTICLE_ROOT" \
  "$LOG_DIR" "$STRUCTURE_LOG_DIR" data/processed/structures data/raw/structures

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

preflight_assets() {
  log_stage "Validate mandatory publication inputs"
  require_file "$RANKING"
  require_file "$CANDIDATES"
  require_file "$DISCOVERY"
  require_file "$FUNCTIONAL_EVIDENCE"
  require_file "$PHARMACOLOGY_DATA/pharmacology_evidence_long.tsv"
  require_file "$STRUCTURE_MANIFEST"
  require_file "$STRUCTURAL_ANNOTATIONS"
  require_file "$STRUCTURE_RENDER_MANIFEST"
  python -m compileall -q src scripts tests
  bash -n scripts/run_publication_pipeline.sh
  printf '%s\n' "$(date -Iseconds)" > "$RUN_MARKER"
  echo "Assets-only run marker: $RUN_MARKER"
}

require_ranking() {
  require_file "$RANKING"
}

acquire_pharmacology() {
  require_ranking
  log_stage "Acquire Open Targets, ChEMBL, DGIdb, MyChem, Pharos and CIViC evidence"
  run_logged "$LOG_DIR/01_acquire_pharmacology.log" \
    python -u scripts/acquire_pharmacology_evidence.py \
      --ranking "$RANKING" \
      --output-dir "$PHARMACOLOGY_DATA" \
      --cache-dir data/raw/pharmacology/api_cache \
      --minimum-vulnerability-score "$PHARMACOLOGY_MINIMUM_SCORE" \
      --max-targets "$PHARMACOLOGY_MAX_TARGETS"
}

standardize_sensitivity() {
  log_stage "Standardize any locally available PRISM, GDSC and CTRP releases"
  run_logged "$LOG_DIR/02_standardize_drug_sensitivity.log" \
    python -u scripts/standardize_drug_sensitivity.py \
      --config config/drug_sensitivity_sources.yaml \
      --models "$MODELS" \
      --output "$PHARMACOLOGY_DATA/drug_sensitivity_long.tsv" \
      --status-output "$PHARMACOLOGY_DATA/drug_sensitivity_source_status.tsv"
}

analyze_sensitivity() {
  log_stage "Test biomarker-matched drug-response selectivity"
  run_logged "$LOG_DIR/03_analyze_drug_response_selectivity.log" \
    python -u scripts/analyze_drug_response_selectivity.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --ranking "$RANKING" \
      --pharmacology-evidence "$PHARMACOLOGY_DATA/pharmacology_evidence_long.tsv" \
      --sensitivity "$PHARMACOLOGY_DATA/drug_sensitivity_long.tsv" \
      --loss-threshold "$LOSS_THRESHOLD" \
      --min-group-size "$MIN_GROUP_SIZE" \
      --output "$PHARMACOLOGY_DATA/drug_response_selectivity.tsv" \
      --evidence-output "$PHARMACOLOGY_DATA/pharmacology_sensitivity_evidence.tsv"
}

prioritize_pharmacology() {
  require_ranking
  log_stage "Build coverage-aware experimental target-drug priorities"
  run_logged "$LOG_DIR/04_prioritize_pharmacology.log" \
    python -u scripts/prioritize_pharmacology.py \
      --ranking "$RANKING" \
      --evidence "$PHARMACOLOGY_DATA/pharmacology_evidence_long.tsv" \
      --sensitivity "$PHARMACOLOGY_DATA/drug_response_selectivity.tsv" \
      --output-dir "$PHARMACOLOGY_RESULTS"
}

build_evidence_audit() {
  log_stage "Audit candidate-domain eligibility, missingness, coverage and evidence overlap"
  run_logged "$LOG_DIR/04b_build_publication_evidence_audit.log" \
    python -u scripts/build_publication_evidence_audit.py \
      --ranking "$RANKING" \
      --functional-evidence "$FUNCTIONAL_EVIDENCE" \
      --pharmacology-evidence "$PHARMACOLOGY_DATA/pharmacology_evidence_long.tsv" \
      --source-metadata-root data/raw \
      --output-root "$ARTICLE_ROOT"
}

run_robustness() {
  log_stage "Run raw-versus-adjusted, leave-one-domain-out and weight-sensitivity analyses"
  run_logged "$LOG_DIR/04c_run_rses_robustness.log" \
    python -u scripts/run_rses_robustness_analyses.py \
      --ranking "$RANKING" \
      --output-root "$ARTICLE_ROOT" \
      --top-k 20 \
      --weight-delta 0.20
}

export_supporting_evidence() {
  log_stage "Export expression, network, phenotype, genomic, structural and pharmacology support tables"
  run_logged "$LOG_DIR/04d_export_supporting_evidence.log" \
    python -u scripts/export_supporting_evidence_tables.py \
      --ranking "$RANKING" \
      --functional-evidence "$FUNCTIONAL_EVIDENCE" \
      --output-root "$ARTICLE_ROOT"
}

prepare_structures() {
  require_ranking
  log_stage "Download every curated human NISE AlphaFold DB structure"
  run_logged "$STRUCTURE_LOG_DIR/01_download_alphafold.log" \
    python -u scripts/download_alphafold_nise_structures.py \
      --proteins "$PROTEINS" --output-dir data/raw/structures/alphafold --manifest "$STRUCTURE_MANIFEST"
  log_stage "Collect exact-numbered M-CSA, UniProt and PDBe residue annotations"
  run_logged "$STRUCTURE_LOG_DIR/02_collect_structural_annotations.log" \
    python -u scripts/collect_nise_structural_annotations.py \
      --proteins "$PROTEINS" --output "$STRUCTURAL_ANNOTATIONS" \
      --coverage-output "$STRUCTURAL_COVERAGE" --cache-dir data/raw/structures/annotation_cache
  log_stage "Render whole structures and enlarged functional-site views with PyMOL"
  run_logged "$STRUCTURE_LOG_DIR/03_render_structures.log" \
    python -u scripts/render_nise_structures.py \
      --structure-manifest "$STRUCTURE_MANIFEST" --annotations "$STRUCTURAL_ANNOTATIONS" \
      --output-dir "$ARTICLE_ROOT/structure_atlas/individual" \
      --render-manifest "$STRUCTURE_RENDER_MANIFEST" --pymol "$PYMOL_EXECUTABLE"
}

export_tables() {
  require_ranking
  log_stage "Export all main, supplementary, audit and robustness article tables"
  run_logged "$LOG_DIR/05_export_article_tables.log" \
    python -u scripts/export_article_tables.py \
      --ranking "$RANKING" --candidates "$CANDIDATES" --members "$MEMBERS" \
      --functional-evidence "$FUNCTIONAL_EVIDENCE" --discovery "$DISCOVERY" \
      --pharmacology-evidence "$PHARMACOLOGY_DATA/pharmacology_evidence_long.tsv" \
      --pharmacology-ranking "$PHARMACOLOGY_RESULTS/pharmacology_ranked_hypotheses.tsv" \
      --drug-sensitivity "$PHARMACOLOGY_DATA/drug_response_selectivity.tsv" \
      --pharmacology-source-status "$PHARMACOLOGY_DATA/pharmacology_source_status.tsv" \
      --pharmacology-source-coverage "$PHARMACOLOGY_RESULTS/pharmacology_source_coverage.tsv" \
      --structure-manifest "$STRUCTURE_MANIFEST" \
      --structural-annotations "$STRUCTURAL_ANNOTATIONS" \
      --structure-render-manifest "$STRUCTURE_RENDER_MANIFEST" \
      --output-root "$ARTICLE_ROOT"
}

generate_figures() {
  require_ranking
  log_stage "Generate every main, supplementary, structural and audit figure from scripts"
  local layout_flag="--strict-layout"
  if [[ "$STRICT_LAYOUT" != "1" ]]; then
    layout_flag="--no-strict-layout"
  fi
  run_logged "$LOG_DIR/06_make_all_article_figures.log" \
    python -u scripts/make_all_article_figures.py \
      --ranking "$RANKING" --candidates "$CANDIDATES" --discovery "$DISCOVERY" \
      --pharmacology "$PHARMACOLOGY_RESULTS/pharmacology_ranked_hypotheses.tsv" \
      --proteins "$PROTEINS" --structure-manifest "$STRUCTURE_MANIFEST" \
      --render-manifest "$STRUCTURE_RENDER_MANIFEST" \
      --structural-annotations "$STRUCTURAL_ANNOTATIONS" \
      --structural-coverage "$STRUCTURAL_COVERAGE" \
      --output-root "$ARTICLE_ROOT" "$layout_flag"
}

catalog_figure_data() {
  log_stage "Catalog exact source tables and reproduction commands for every figure"
  run_logged "$LOG_DIR/06b_catalog_figure_source_data.log" \
    python -u scripts/catalog_figure_source_data.py --article-root "$ARTICLE_ROOT"
}

validate_scientific_integrity() {
  log_stage "Validate missingness, score formula, overlap control, FDR and figure semantics"
  run_logged "$LOG_DIR/06c_validate_scientific_integrity.log" \
    python -u scripts/validate_publication_scientific_integrity.py \
      --article-root "$ARTICLE_ROOT" --run-marker "$RUN_MARKER"
}

build_workbook() {
  log_stage "Build organized article and supplementary workbook"
  run_logged "$LOG_DIR/07_build_article_workbook.log" \
    python -u scripts/build_article_workbook.py \
      --article-root "$ARTICLE_ROOT" \
      --output "$ARTICLE_ROOT/workbooks/RSES_Onco_Article_Tables_and_Evidence.xlsx"
}

build_manifests() {
  log_stage "Build publication inventory, provenance and checksums"
  run_logged "$LOG_DIR/08_build_publication_manifest.log" \
    python -u scripts/build_publication_manifest.py \
      --article-root "$ARTICLE_ROOT" \
      --input "$RANKING" --input "$CANDIDATES" --input "$DISCOVERY" \
      --input "$FUNCTIONAL_EVIDENCE" \
      --input "$PHARMACOLOGY_DATA/pharmacology_evidence_long.tsv" \
      --input "$PHARMACOLOGY_RESULTS/pharmacology_ranked_hypotheses.tsv" \
      --input "$STRUCTURE_MANIFEST" --input "$STRUCTURAL_ANNOTATIONS" \
      --input "$STRUCTURE_RENDER_MANIFEST"
}

validate_outputs() {
  log_stage "Validate publication package and software tests"
  run_logged "$LOG_DIR/09_validate_publication_outputs.log" \
    python -u scripts/validate_publication_outputs.py --article-root "$ARTICLE_ROOT"
  run_logged "$LOG_DIR/10_pytest.log" \
    python -m pytest -q -p no:cacheprovider
}

all() {
  preflight_assets
  acquire_pharmacology
  standardize_sensitivity
  analyze_sensitivity
  prioritize_pharmacology
  build_evidence_audit
  run_robustness
  export_supporting_evidence
  prepare_structures
  export_tables
  generate_figures
  catalog_figure_data
  validate_scientific_integrity
  build_workbook
  build_manifests
  validate_outputs
}

assets_only() {
  preflight_assets
  prioritize_pharmacology
  build_evidence_audit
  run_robustness
  export_supporting_evidence
  export_tables
  generate_figures
  catalog_figure_data
  validate_scientific_integrity
  build_workbook
  build_manifests
  validate_outputs
}

usage() {
  cat <<'EOF'
Usage: bash scripts/run_publication_pipeline.sh STAGE

Stages:
  acquire-pharmacology     Query/cached pharmacology evidence APIs
  standardize-sensitivity Standardize local PRISM, GDSC and CTRP releases
  analyze-sensitivity     Test biomarker-matched drug response
  prioritize              Build target-drug actionability priorities
  evidence-audit          Build candidate-domain coverage, missingness and overlap audits
  robustness              Run leave-one-domain-out and controlled weight sensitivity
  supporting-evidence     Export expression, network, phenotype, structural and pharmacology tables
  structures              Download, annotate and render all human NISE structures
  tables                  Export all main and supplementary tables
  figures                 Generate all main, supplementary, structural and audit figures
  figure-data             Catalog the exact table used by every figure
  scientific-validate     Validate score, missingness, overlap, FDR and figure semantics
  workbook                Build organized Excel workbook
  manifests               Build file inventory, provenance and SHA-256 checksums
  validate                Validate figure triplets, audits, tables and tests
  assets-only             Rebuild all publication assets from cached evidence and structures
  all                     Run complete pharmacology, structural and publication workflow
EOF
}

stage="${1:-}"
case "$stage" in
  acquire-pharmacology) acquire_pharmacology ;;
  standardize-sensitivity) standardize_sensitivity ;;
  analyze-sensitivity) analyze_sensitivity ;;
  prioritize) prioritize_pharmacology ;;
  evidence-audit) build_evidence_audit ;;
  robustness) run_robustness ;;
  supporting-evidence) export_supporting_evidence ;;
  structures) prepare_structures ;;
  tables) export_tables ;;
  figures) generate_figures ;;
  figure-data) catalog_figure_data ;;
  scientific-validate) validate_scientific_integrity ;;
  workbook) build_workbook ;;
  manifests) build_manifests ;;
  validate) validate_outputs ;;
  assets-only) assets_only ;;
  all) all ;;
  *) usage; exit 2 ;;
esac
