#!/usr/bin/env bash
# Complete publication workflow wrapper with genomic Circos and script documentation.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CORE="$ROOT/scripts/publication_pipeline_steps.sh"
RANKING="${RANKING:-results/expanded_26Q1/full/expanded_rses_onco.tsv}"
CANDIDATES="${CANDIDATES:-data/processed/expanded_candidate_universe.tsv}"
FUNCTIONAL_EVIDENCE="${FUNCTIONAL_EVIDENCE:-data/processed/expanded_pair_functional_evidence.tsv}"
ARTICLE_ROOT="${ARTICLE_ROOT:-article_outputs}"
LOG_DIR="${LOG_DIR:-logs/publication_26Q1}"
DEPMAP_DIR="${DEPMAP_DIR:-data/raw/depmap}"
GENE_EFFECT="${GENE_EFFECT:-$DEPMAP_DIR/CRISPRGeneEffect.csv}"
COPY_NUMBER="${COPY_NUMBER:-$DEPMAP_DIR/OmicsCNGeneWGS.csv}"
MODELS="${MODELS:-$DEPMAP_DIR/Model.csv}"
EXPRESSION="${EXPRESSION:-$DEPMAP_DIR/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv}"
PROMOTERS="${PROMOTERS:-data/raw/regulatory/ensembl_promoters.tsv}"
WGCNA_PAIR_METRICS="${WGCNA_PAIR_METRICS:-data/processed/regulatory/wgcna/wgcna_pair_metrics_all_cancers.tsv}"
FUNCTIONAL_RAW_DIR="${FUNCTIONAL_RAW_DIR:-data/raw/human_functional_evidence}"
LOSS_THRESHOLD="${LOSS_THRESHOLD:-0.30}"
STRICT_LAYOUT="${STRICT_LAYOUT:-1}"

mkdir -p "$LOG_DIR" "$ARTICLE_ROOT"

log_stage() {
  printf '\n[%s] %s\n' "$(date -Iseconds)" "$*"
}

run_logged() {
  local log_file="$1"
  shift
  set -o pipefail
  PYTHONUNBUFFERED=1 "$@" 2>&1 | tee "$log_file"
  local status=${PIPESTATUS[0]}
  [[ $status -eq 0 ]] || {
    echo "Command failed with status $status; see $log_file" >&2
    return "$status"
  }
}

require_file() {
  [[ -s "$1" ]] || {
    echo "Mandatory publication input is absent or empty: $1" >&2
    exit 1
  }
}

build_extended_supporting_evidence() {
  for path in "$RANKING" "$FUNCTIONAL_EVIDENCE" "$GENE_EFFECT" \
    "$COPY_NUMBER" "$MODELS" "$EXPRESSION"; do
    require_file "$path"
  done
  log_stage "Build model-level expression, copy-number, coexpression and CRISPR evidence"
  run_logged "$LOG_DIR/04e_build_model_level_supporting_evidence.log" \
    python -u scripts/build_model_level_supporting_evidence.py \
      --gene-effect "$GENE_EFFECT" --copy-number "$COPY_NUMBER" \
      --models "$MODELS" --expression "$EXPRESSION" \
      --ranking "$RANKING" --output-root "$ARTICLE_ROOT" \
      --loss-threshold "$LOSS_THRESHOLD"
  log_stage "Export raw STRING, DoRothEA, HPA and UniProt candidate evidence"
  run_logged "$LOG_DIR/04f_export_raw_functional_network_evidence.log" \
    python -u scripts/export_raw_functional_network_evidence.py \
      --ranking "$RANKING" --raw-dir "$FUNCTIONAL_RAW_DIR" \
      --output-root "$ARTICLE_ROOT"
  log_stage "Export WGCNA, promoter and TF regulatory evidence"
  run_logged "$LOG_DIR/04g_export_wgcna_regulatory_evidence.log" \
    python -u scripts/export_wgcna_regulatory_supporting_evidence.py \
      --article-root "$ARTICLE_ROOT"
  log_stage "Run WGCNA, pairwise-expression, TF and promoter ablations"
  run_logged "$LOG_DIR/04h_wgcna_regulatory_ablation.log" \
    python -u scripts/run_wgcna_regulatory_ablation.py \
      --ranking "$RANKING" --article-root "$ARTICLE_ROOT" --top-k 20
}

build_repository_documentation_and_circos_data() {
  for path in "$RANKING" "$CANDIDATES" "$PROMOTERS" \
    "$EXPRESSION" "$MODELS" "$WGCNA_PAIR_METRICS"; do
    require_file "$path"
  done
  log_stage "Document every Python, Bash and R pipeline script/module"
  run_logged "$LOG_DIR/04k_build_script_documentation.log" \
    python -u scripts/build_script_documentation.py \
      --output-md docs/SCRIPT_CATALOG.md \
      --output-tsv docs/script_manifest.tsv \
      --processed-output data/processed/documentation/pipeline_script_catalog.tsv
  log_stage "Build GRCh38 coordinates, NISE/paralog links, score rings and expression tables"
  run_logged "$LOG_DIR/04l_build_genomic_circos_inputs.log" \
    python -u scripts/build_genomic_circos_inputs.py \
      --ranking "$RANKING" --candidates "$CANDIDATES" \
      --promoters "$PROMOTERS" --expression "$EXPRESSION" \
      --models "$MODELS" --output-dir data/processed/circos
  log_stage "Add WGCNA, methylation and validation internal layers to Circos rings"
  run_logged "$LOG_DIR/04l1_enrich_genomic_circos_internal_layers.log" \
    python -u scripts/enrich_genomic_circos_internal_layers.py \
      --ranking "$RANKING" --candidates "$CANDIDATES" \
      --wgcna "$WGCNA_PAIR_METRICS" \
      --coordinates data/processed/circos/genomic_circos_gene_coordinates.tsv \
      --rings data/processed/circos/genomic_circos_ring_values.tsv \
      --tracks data/processed/circos/genomic_circos_track_definitions.tsv
  log_stage "Complete expression availability for every Circos gene and cancer"
  run_logged "$LOG_DIR/04la_complete_genomic_circos_expression.log" \
    python -u scripts/complete_genomic_circos_expression_summary.py \
      --coordinates data/processed/circos/genomic_circos_gene_coordinates.tsv \
      --model-values data/processed/circos/genomic_circos_expression_model_values.tsv \
      --summary data/processed/circos/genomic_circos_expression_summary.tsv
  log_stage "Guarantee one chord for every simple NISE/paralog candidate pair"
  run_logged "$LOG_DIR/04m_complete_genomic_circos_links.log" \
    python -u scripts/complete_genomic_circos_links.py \
      --candidates "$CANDIDATES" \
      --coordinates data/processed/circos/genomic_circos_gene_coordinates.tsv \
      --ranking "$RANKING" \
      --links data/processed/circos/genomic_circos_pair_links.tsv
  log_stage "Finalize Circos counts and full SHA-256 source provenance"
  run_logged "$LOG_DIR/04n_finalize_genomic_circos_status.log" \
    python -u scripts/finalize_genomic_circos_status.py \
      --output-dir data/processed/circos \
      --ranking "$RANKING" --candidates "$CANDIDATES" \
      --promoters "$PROMOTERS" --expression "$EXPRESSION" \
      --models "$MODELS" --wgcna "$WGCNA_PAIR_METRICS"
  log_stage "Validate finalized Circos counts and complete provenance"
  run_logged "$LOG_DIR/04o_validate_genomic_circos_final_status.log" \
    python -u scripts/validate_genomic_circos_final_status.py \
      --circos-dir data/processed/circos
}

run_core_without_terminal_validation() {
  local requested_stage="$1"
  if [[ "$requested_stage" == "all" ]]; then
    for stage in acquire-pharmacology standardize-sensitivity \
      analyze-sensitivity prioritize evidence-audit robustness \
      supporting-evidence extended-tables structures; do
      bash "$CORE" "$stage"
    done
  else
    for stage in prioritize evidence-audit robustness supporting-evidence \
      extended-tables; do
      bash "$CORE" "$stage"
    done
  fi
  for stage in tables figures figure-data scientific-validate workbook manifests; do
    bash "$CORE" "$stage"
  done
}

integrate_genomic_circos_assets() {
  local layout_flag="--strict-layout"
  [[ "$STRICT_LAYOUT" == "1" ]] || layout_flag="--no-strict-layout"
  log_stage "Generate genomic Circos Figure S70 from exact source tables"
  run_logged "$LOG_DIR/06i_make_genomic_circos_figure.log" \
    python -u scripts/make_genomic_circos_figure_resilient.py \
      --config config/genomic_circos_asset.yaml \
      --coordinates data/processed/circos/genomic_circos_gene_coordinates.tsv \
      --links data/processed/circos/genomic_circos_pair_links.tsv \
      --ring-values data/processed/circos/genomic_circos_ring_values.tsv \
      --tracks data/processed/circos/genomic_circos_track_definitions.tsv \
      --output-root "$ARTICLE_ROOT" "$layout_flag"
  log_stage "Register Figure S70 and Supplementary Tables S45-S52"
  run_logged "$LOG_DIR/06j_register_genomic_circos_assets.log" \
    python -u scripts/register_genomic_circos_assets.py \
      --article-root "$ARTICLE_ROOT"
  log_stage "Generate genomic Circos supplementary methods"
  run_logged "$LOG_DIR/06k_build_genomic_circos_methods.log" \
    python -u scripts/build_genomic_circos_methods.py \
      --tracks data/processed/circos/genomic_circos_track_definitions.tsv \
      --output "$ARTICLE_ROOT/manuscript_assets/supplementary_methods/GENOMIC_CIRCOS_METHODS.md"
  log_stage "Recatalog exact source data for all 78 figures"
  run_logged "$LOG_DIR/06l_recatalog_figure_source_data.log" \
    python -u scripts/catalog_figure_source_data.py \
      --article-root "$ARTICLE_ROOT"
  log_stage "Validate genomic Circos coordinates, links, rings, expression and script catalogue"
  run_logged "$LOG_DIR/06m_validate_genomic_circos_integrity.log" \
    python -u scripts/validate_genomic_circos_integrity.py \
      --article-root "$ARTICLE_ROOT" --candidates "$CANDIDATES"
}

finalize_extended_publication_assets() {
  log_stage "Validate model-level and raw-source supporting evidence"
  run_logged "$LOG_DIR/06d_validate_extended_supporting_evidence.log" \
    python -u scripts/validate_extended_supporting_evidence.py \
      --article-root "$ARTICLE_ROOT"
  log_stage "Validate WGCNA, promoter and methylation-aware regulatory integration"
  run_logged "$LOG_DIR/06e_validate_wgcna_regulatory_evidence.log" \
    python -u scripts/validate_wgcna_regulatory_evidence.py \
      --ranking "$RANKING" --functional-evidence "$FUNCTIONAL_EVIDENCE" \
      --article-root "$ARTICLE_ROOT"
  log_stage "Generate score-formula, evidence-category and reproduction documentation"
  run_logged "$LOG_DIR/06f_build_publication_methods_documentation.log" \
    python -u scripts/build_publication_methods_documentation.py \
      --article-root "$ARTICLE_ROOT"
  log_stage "Generate WGCNA and promoter-aware regulatory supplementary methods"
  run_logged "$LOG_DIR/06g_build_wgcna_regulatory_methods.log" \
    python -u scripts/build_wgcna_regulatory_methods.py \
      --article-root "$ARTICLE_ROOT"
  log_stage "Create mandatory manual 100-percent-zoom inspection checklist"
  run_logged "$LOG_DIR/06h_create_manual_visual_inspection_checklist.log" \
    python -u scripts/create_manual_visual_inspection_checklist.py \
      --article-root "$ARTICLE_ROOT"
  bash "$CORE" workbook
  bash "$CORE" manifests
  bash "$CORE" validate
}

stage="${1:-}"
case "$stage" in
  all|assets-only)
    build_extended_supporting_evidence
    build_repository_documentation_and_circos_data
    run_core_without_terminal_validation "$stage"
    integrate_genomic_circos_assets
    finalize_extended_publication_assets
    ;;
  *) exec bash "$CORE" "$@" ;;
esac
