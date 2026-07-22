#!/usr/bin/env bash
# Canonical publication-pipeline entry point.
#
# The implementation is kept in publication_pipeline_steps.sh so each stage remains
# directly testable while this stable command continues to support:
#
#   bash scripts/run_publication_pipeline.sh assets-only
#
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec bash "$ROOT/scripts/publication_pipeline_steps.sh" "$@"
