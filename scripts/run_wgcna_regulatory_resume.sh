#!/usr/bin/env bash
# Resume after completed STRING/DoRothEA/HPA/UniProt acquisition and rebuild
# WGCNA/promoter-aware RSES-Onco scores plus all publication assets.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DEPMAP_DIR="${DEPMAP_DIR:-data/raw/depmap}"
GDC_DIR="${GDC_DIR:-data/raw/gdc}"
PROCESSED_DIR="${PROCESSED_DIR:-data/processed}"
RESULT_ROOT="${RESULT_ROOT:-results/expanded_26Q1}"
DEPMAP_RESULTS="${DEPMAP_RESULTS:-$RESULT_ROOT/depmap_only}"
FULL_RESULTS="${FULL_RESULTS:-$RESULT_ROOT/full}"
DISCOVERY_RESULTS="${DISCOVERY_RESULTS:-$RESULT_ROOT/discovery}"
LOG_DIR="${LOG_DIR:-logs/wgcna_regulatory_26Q1}"
ARTICLE_ROOT="${ARTICLE_ROOT:-article_outputs}"

GENE_EFFECT="${GENE_EFFECT:-$DEPMAP_DIR/CRISPRGeneEffect.csv}"
COPY_NUMBER="${COPY_NUMBER:-$DEPMAP_DIR/OmicsCNGeneWGS.csv}"
MODELS="${MODELS:-$DEPMAP_DIR/Model.csv}"
EXPRESSION="${EXPRESSION:-$DEPMAP_DIR/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv}"
CANDIDATES="${CANDIDATES:-$PROCESSED_DIR/expanded_candidate_universe.tsv}"
MEMBERS="${MEMBERS:-$PROCESSED_DIR/expanded_class_member_inventory.tsv}"
DISCOVERY="${DISCOVERY:-$DISCOVERY_RESULTS/all_target_dependency_screen.tsv}"
DOROTHEA="${DOROTHEA:-data/raw/human_functional_evidence/omnipath_dorothea.tsv}"

FUNCTIONAL_EVIDENCE="${FUNCTIONAL_EVIDENCE:-$PROCESSED_DIR/expanded_pair_functional_evidence.tsv}"
FUNCTIONAL_EVIDENCE_BASE="${FUNCTIONAL_EVIDENCE_BASE:-$PROCESSED_DIR/expanded_pair_functional_evidence_base.tsv}"
FUNCTIONAL_EVIDENCE_CANCER="${FUNCTIONAL_EVIDENCE_CANCER:-$PROCESSED_DIR/regulatory/expanded_pair_functional_evidence_cancer_specific.tsv}"

PROMOTER_TABLE="${PROMOTER_TABLE:-data/raw/regulatory/ensembl_promoters.tsv}"
PROMOTER_FASTA="${PROMOTER_FASTA:-data/raw/regulatory/ensembl_promoters.fa}"
JASPAR_MEME="${JASPAR_MEME:-data/raw/regulatory/JASPAR2026_CORE_vertebrates_non-redundant.meme}"
PROMOTER_MOTIF_HITS="${PROMOTER_MOTIF_HITS:-data/processed/regulatory/jaspar_promoter_motif_hits.tsv}"
PROMOTER_TF_SUMMARY="${PROMOTER_TF_SUMMARY:-data/processed/regulatory/jaspar_promoter_tf_summary.tsv}"

LOSS_THRESHOLD="${LOSS_THRESHOLD:-0.30}"
MIN_GROUP_SIZE="${MIN_GROUP_SIZE:-3}"
STRICT_LAYOUT="${STRICT_LAYOUT:-1}"

mkdir -p "$LOG_DIR" "$DEPMAP_RESULTS" "$FULL_RESULTS" \
  "$PROCESSED_DIR/regulatory" "$ARTICLE_ROOT"

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
  local description="$2"
  if [[ ! -s "$path" ]]; then
    echo "Required input is absent or empty: $description: $path" >&2
    exit 1
  fi
}

preflight() {
  log_stage "Validate WGCNA/promoter regulatory prerequisites"
  for item in \
    "$GENE_EFFECT" \
    "$COPY_NUMBER" \
    "$MODELS" \
    "$EXPRESSION" \
    "$CANDIDATES" \
    "$MEMBERS" \
    "$DISCOVERY" \
    "$FUNCTIONAL_EVIDENCE" \
    "$DOROTHEA"
  do
    require_file "$item" "mandatory real-data input"
  done

  command -v Rscript >/dev/null 2>&1 || {
    echo "Rscript is missing. Update the rses-onco Conda environment." >&2
    exit 1
  }
  command -v fimo >/dev/null 2>&1 || {
    echo "FIMO is missing. Install the MEME suite in the rses-onco environment." >&2
    exit 1
  }
  Rscript -e "stopifnot(requireNamespace('WGCNA', quietly=TRUE)); stopifnot(requireNamespace('dynamicTreeCut', quietly=TRUE)); cat('WGCNA R dependencies: OK\n')"
  python -m compileall -q src scripts tests
  bash -n scripts/run_wgcna_regulatory_resume.sh
  bash -n scripts/run_publication_pipeline.sh
}

preserve_base_evidence() {
  log_stage "Preserve immutable pre-WGCNA functional evidence"
  if [[ -s "$FUNCTIONAL_EVIDENCE_BASE" ]]; then
    echo "Reusing $FUNCTIONAL_EVIDENCE_BASE"
    return
  fi
  python - "$FUNCTIONAL_EVIDENCE" "$FUNCTIONAL_EVIDENCE_BASE" <<'PY'
import sys
from pathlib import Path
import pandas as pd

source = Path(sys.argv[1])
output = Path(sys.argv[2])
frame = pd.read_csv(source, sep="\t", low_memory=False)
new_columns = {
  "component_wgcna_expression_network",
  "component_regulatory_network_composite",
  "regulatory_layer_version",
}
if new_columns & set(frame.columns):
  raise SystemExit(
    "The current functional-evidence table is already enriched, but the immutable "
    "base table is absent. Restore or reacquire the pre-WGCNA evidence table."
  )
output.parent.mkdir(parents=True, exist_ok=True)
temporary = output.with_suffix(output.suffix + ".tmp")
frame.to_csv(temporary, sep="\t", index=False)
temporary.replace(output)
print(f"Preserved {len(frame):,} rows: {output}")
PY
}

build_promoter_evidence() {
  log_stage "Acquire canonical Ensembl promoter windows and sequences"
  run_logged "$LOG_DIR/01_download_ensembl_promoters.log" \
    python -u scripts/download_ensembl_promoters.py \
      --candidates "$CANDIDATES" \
      --output "$PROMOTER_TABLE" \
      --fasta "$PROMOTER_FASTA"

  log_stage "Acquire official JASPAR 2026 CORE vertebrate motifs"
  run_logged "$LOG_DIR/02_download_jaspar_core.log" \
    python -u scripts/download_jaspar_core_vertebrates.py \
      --output "$JASPAR_MEME"

  log_stage "Scan promoters with JASPAR motifs using FIMO"
  run_logged "$LOG_DIR/03_scan_promoter_motifs.log" \
    python -u scripts/scan_promoter_motifs.py \
      --motifs "$JASPAR_MEME" \
      --promoters "$PROMOTER_FASTA" \
      --output "$PROMOTER_MOTIF_HITS" \
      --summary-output "$PROMOTER_TF_SUMMARY"
}

build_network_evidence() {
  log_stage "Build cancer-specific signed WGCNA and promoter-aware TF evidence"
  run_logged "$LOG_DIR/04_build_wgcna_regulatory_layer.log" \
    python -u scripts/build_wgcna_regulatory_layer.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --expression "$EXPRESSION" \
      --candidates "$CANDIDATES" \
      --functional-evidence "$FUNCTIONAL_EVIDENCE_BASE" \
      --dorothea "$DOROTHEA" \
      --promoter-motifs "$PROMOTER_TF_SUMMARY" \
      --output "$FUNCTIONAL_EVIDENCE_CANCER"

  log_stage "Build consensus pair evidence for non-cancer-specific downstream assets"
  run_logged "$LOG_DIR/05_aggregate_wgcna_regulatory_layer.log" \
    python -u scripts/aggregate_wgcna_regulatory_layer.py \
      --base "$FUNCTIONAL_EVIDENCE_BASE" \
      --cancer-specific "$FUNCTIONAL_EVIDENCE_CANCER" \
      --output "$FUNCTIONAL_EVIDENCE"
}

recompute_one_ranking() {
  local output="$1"
  shift
  local pre_wgcna="${output%.tsv}_pre_wgcna_regulatory.tsv"
  local stem
  stem="$(basename "${output%.tsv}")"

  run_logged "$LOG_DIR/${stem}_base_scoring.log" \
    python -u scripts/run_expanded_rses_onco.py \
      --gene-effect "$GENE_EFFECT" \
      --copy-number "$COPY_NUMBER" \
      --models "$MODELS" \
      --expression "$EXPRESSION" \
      --candidates "$CANDIDATES" \
      --functional-evidence "$FUNCTIONAL_EVIDENCE" \
      --loss-threshold "$LOSS_THRESHOLD" \
      --min-group-size "$MIN_GROUP_SIZE" \
      "$@" \
      --output "$output"

  cp -f "$output" "$pre_wgcna"

  run_logged "$LOG_DIR/${stem}_wgcna_regulatory_recompute.log" \
    python -u scripts/recompute_rses_with_wgcna_regulatory.py \
      --ranking "$pre_wgcna" \
      --functional-evidence "$FUNCTIONAL_EVIDENCE_CANCER" \
      --output "$output"

  run_logged "$LOG_DIR/${stem}_ranking_contract.log" \
    python -u scripts/stamp_wgcna_ranking_contract.py \
      --ranking "$output"
}

ensure_tcga_matrices() {
  local missing=0
  for path in \
    "$PROCESSED_DIR/TCGA_COLON_homdel_discrete.tsv" \
    "$PROCESSED_DIR/TCGA_STOMACH_homdel_discrete.tsv" \
    "$PROCESSED_DIR/TCGA_LUNG_homdel_discrete.tsv"
  do
    [[ -s "$path" ]] || missing=1
  done
  if [[ $missing -eq 0 ]]; then
    echo "Reusing existing TCGA/GDC cancer matrices."
    return
  fi

  log_stage "Rebuild missing TCGA/GDC cancer matrices"
  run_logged "$LOG_DIR/06_validate_gdc.log" \
    python -u scripts/download_gdc.py \
      --validate-only \
      --manifest "$GDC_DIR/gdc_gene_level_copy_number_manifest.json" \
      --output-dir "$GDC_DIR"
  run_logged "$LOG_DIR/07_aggregate_gdc.log" \
    python -u scripts/aggregate_gdc_gene_cna.py \
      --raw-dir "$GDC_DIR" \
      --manifest "$GDC_DIR/gdc_gene_level_copy_number_manifest.json" \
      --output-dir "$PROCESSED_DIR"
}

recompute_rankings() {
  log_stage "Recompute DepMap-only WGCNA/promoter-aware RSES-Onco"
  recompute_one_ranking "$DEPMAP_RESULTS/expanded_rses_onco.tsv"

  ensure_tcga_matrices

  log_stage "Recompute integrated TCGA plus DepMap RSES-Onco"
  recompute_one_ranking \
    "$FULL_RESULTS/expanded_rses_onco.tsv" \
    --tcga "colon=$PROCESSED_DIR/TCGA_COLON_homdel_discrete.tsv" \
    --tcga "stomach=$PROCESSED_DIR/TCGA_STOMACH_homdel_discrete.tsv" \
    --tcga "lung=$PROCESSED_DIR/TCGA_LUNG_homdel_discrete.tsv"
}

rebuild_publication() {
  log_stage "Rebuild publication assets from cached pharmacology and structures"
  require_file "data/processed/pharmacology/pharmacology_evidence_long.tsv" \
    "cached pharmacology evidence"
  require_file "data/processed/structures/alphafold_structure_manifest.tsv" \
    "cached AlphaFold manifest"
  require_file "data/processed/structures/nise_structural_residue_annotations.tsv" \
    "cached structural annotations"
  require_file "data/processed/structures/nise_structure_render_manifest.tsv" \
    "cached structure-render manifest"

  set -o pipefail
  RANKING="$FULL_RESULTS/expanded_rses_onco.tsv" \
  CANDIDATES="$CANDIDATES" \
  MEMBERS="$MEMBERS" \
  DISCOVERY="$DISCOVERY" \
  FUNCTIONAL_EVIDENCE="$FUNCTIONAL_EVIDENCE" \
  DEPMAP_DIR="$DEPMAP_DIR" \
  GENE_EFFECT="$GENE_EFFECT" \
  COPY_NUMBER="$COPY_NUMBER" \
  MODELS="$MODELS" \
  EXPRESSION="$EXPRESSION" \
  ARTICLE_ROOT="$ARTICLE_ROOT" \
  STRICT_LAYOUT="$STRICT_LAYOUT" \
  bash scripts/run_publication_pipeline.sh assets-only \
    2>&1 | tee "$LOG_DIR/08_publication_assets_only.log"
  local status=${PIPESTATUS[0]}
  if [[ $status -ne 0 ]]; then
    return "$status"
  fi

  run_logged "$LOG_DIR/09_export_wgcna_regulatory_support.log" \
    python -u scripts/export_wgcna_regulatory_supporting_evidence.py \
      --article-root "$ARTICLE_ROOT"

  ARTICLE_ROOT="$ARTICLE_ROOT" bash scripts/publication_pipeline_steps.sh workbook
  ARTICLE_ROOT="$ARTICLE_ROOT" bash scripts/publication_pipeline_steps.sh manifests
  ARTICLE_ROOT="$ARTICLE_ROOT" bash scripts/publication_pipeline_steps.sh validate
}

finalize() {
  log_stage "Validate ranking extension and write final checksums"
  python - "$FULL_RESULTS/expanded_rses_onco.tsv" <<'PY'
import sys
import pandas as pd

frame = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)
required = {
  "scoring_extension_version",
  "component_wgcna_expression_network",
  "regulatory_tf_association_divergence",
  "regulatory_tf_expression_profile_divergence",
  "regulatory_promoter_motif_divergence",
  "direct_promoter_binding_claim",
}
missing = sorted(required - set(frame.columns))
if missing:
  raise SystemExit("WGCNA/regulatory ranking fields missing: " + ", ".join(missing))
assert set(frame["scoring_extension_version"].dropna().astype(str)) == {
  "wgcna-promoter-regulatory-v1"
}
assert not frame["direct_promoter_binding_claim"].fillna(False).astype(bool).any()
print(f"WGCNA/promoter-aware ranking validated: {len(frame):,} rows")
PY

  python -m pytest -q -p no:cacheprovider
  find "$RESULT_ROOT" "$ARTICLE_ROOT" \
    -type f ! -name SHA256SUMS.txt -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "$FULL_RESULTS/SHA256SUMS.txt"
  echo "Wrote $FULL_RESULTS/SHA256SUMS.txt"
}

main() {
  preflight
  preserve_base_evidence
  build_promoter_evidence
  build_network_evidence
  recompute_rankings
  rebuild_publication
  finalize
  log_stage "WGCNA/promoter-aware RSES-Onco resume completed"
}

main "$@"
