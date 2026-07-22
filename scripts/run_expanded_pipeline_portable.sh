#!/usr/bin/env bash
# Portable launcher for the expanded RSES-Onco pipeline.
#
# Direct execution of Python files under scripts/ makes Python place the scripts
# directory, rather than the repository root, at sys.path[0]. Some recovery
# entry points import sibling modules through the `scripts.*` namespace. Exporting
# the repository root through PYTHONPATH makes those imports deterministic in
# Linux, WSL, tmux, cron and non-interactive shells.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"

if [[ "${1:-}" == "--check-python-imports" ]]; then
  cd "$ROOT"
  python - <<'PY'
import scripts.download_human_functional_evidence
import scripts.download_human_functional_evidence_resilient
print("RSES-Onco script namespace imports: OK")
PY
  exit 0
fi

exec bash "$ROOT/scripts/run_expanded_pipeline.sh" "$@"
