#!/usr/bin/env bash
# Canonical publication-pipeline entry point.
#
# `publication_pipeline_complete.sh` integrates model-level evidence, raw functional
# source exports, the established publication stages, extended validation,
# reproduction documentation and the mandatory manual-inspection checklist.
#
#   bash scripts/run_publication_pipeline.sh assets-only
#   bash scripts/run_publication_pipeline.sh all
#
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec bash "$ROOT/scripts/publication_pipeline_complete.sh" "$@"
