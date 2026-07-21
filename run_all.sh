#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

python scripts/run_literature_pilot.py
python scripts/build_relation_forest.py
python scripts/make_figures.py
PYTHONDONTWRITEBYTECODE=1 pytest -q -p no:cacheprovider

cat <<'MSG'
Bundled literature pilot, relation forest, figures and tests completed.
The current empirical TCGA/DepMap analysis is intentionally not run unless the official release files are placed under data/raw/depmap and processed GDC tables are supplied.
MSG
