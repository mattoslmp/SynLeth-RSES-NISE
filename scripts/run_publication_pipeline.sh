#!/usr/bin/env bash
# Canonical publication-pipeline entry point.
#
# The implementation is kept in publication_pipeline_steps.sh so each stage remains
# directly testable while this stable command continues to support:
#
#   bash scripts/run_publication_pipeline.sh assets-only
#
# The assets-only contract executed by publication_pipeline_steps.sh includes:
#
#   build_model_level_supporting_evidence.py
#   export_raw_functional_network_evidence.py
#   build_publication_evidence_audit_complete.py
#   run_rses_robustness_analyses.py
#   export_supporting_evidence_tables.py
#   export_article_tables.py
#   make_all_article_figures.py
#   catalog_figure_source_data.py
#   validate_publication_scientific_integrity.py
#   build_article_workbook.py
#   build_publication_manifest.py
#   validate_publication_outputs.py
#
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec bash "$ROOT/scripts/publication_pipeline_steps.sh" "$@"
