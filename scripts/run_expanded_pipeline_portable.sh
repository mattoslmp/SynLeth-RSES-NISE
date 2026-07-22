#!/usr/bin/env bash
# Portable launcher for the expanded RSES-Onco pipeline.
#
# Direct execution of Python files under scripts/ makes Python place the scripts
# directory, rather than the repository root, at sys.path[0]. Some recovery
# entry points import sibling modules through the `scripts.*` namespace. Exporting
# the repository root through PYTHONPATH makes those imports deterministic in
# Linux, WSL, tmux, cron and non-interactive shells.
#
# The launcher also prepares the official OmniPath static DoRothEA backup before
# stages that acquire functional evidence. This avoids repeated long waits when
# the primary OmniPath REST service returns HTTP 5xx or connection timeouts.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"

if [[ "${1:-}" == "--check-python-imports" ]]; then
  cd "$ROOT"
  python - <<'PY'
import scripts.download_dorothea_official_static
import scripts.download_human_functional_evidence
import scripts.download_human_functional_evidence_resilient
print("RSES-Onco script namespace imports: OK")
PY
  exit 0
fi

stage="${1:-}"
case "$stage" in
  setup|resume-functional|after-download|all)
    if [[ -z "${DOROTHEA_FILE:-}" && "${DOROTHEA_DISABLE_STATIC:-0}" != "1" ]]; then
      dorothea_static="$ROOT/data/raw/human_functional_evidence/dorothea_official_static_9606.tsv"
      dorothea_status="$ROOT/data/raw/human_functional_evidence/dorothea_official_static_9606_status.json"
      command=(
        python -u "$ROOT/scripts/download_dorothea_official_static.py"
        --output "$dorothea_static"
        --status-output "$dorothea_status"
        --levels "${DOROTHEA_LEVELS:-A,B,C}"
        --retries "${DOROTHEA_STATIC_RETRIES:-4}"
        --timeout "${DOROTHEA_STATIC_TIMEOUT:-120}"
      )
      if [[ "${DOROTHEA_REFRESH_STATIC:-0}" == "1" ]]; then
        command+=(--refresh)
      fi

      echo "Preparing official OmniPath static DoRothEA fallback..."
      if "${command[@]}"; then
        export DOROTHEA_FILE="$dorothea_static"
        echo "DoRothEA source selected: $DOROTHEA_FILE"
      elif [[ "${DOROTHEA_STRICT:-0}" == "1" ]]; then
        echo "Official static DoRothEA acquisition failed in strict mode." >&2
        exit 1
      else
        echo "WARNING: official static DoRothEA acquisition failed; falling back to the existing non-strict REST recovery logic." >&2
      fi
    fi
    ;;
esac

exec bash "$ROOT/scripts/run_expanded_pipeline.sh" "$@"
