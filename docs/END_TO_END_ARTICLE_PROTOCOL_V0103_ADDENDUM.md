# RSES-Onco v0.10.3 addendum: DoRothEA/OmniPath outage recovery

This addendum supplements `docs/END_TO_END_ARTICLE_PROTOCOL.md` and
`docs/STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md`.

Use it when STRING completed but the same functional-evidence stage stopped because
OmniPath returned persistent HTTP 5xx responses while acquiring DoRothEA.

The complete scientific rationale, source provenance, cache files, recovery command
and post-run validation are documented in:

```text
docs/DOROTHEA_RECOVERY_WORKFLOW.md
```

The one-command recovery is:

```bash
OLD="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE"
NEW="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010"

cd "$NEW" || exit 1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco

export DEPMAP_DIR="$OLD/data/raw/depmap"
export GDC_DIR="$OLD/data/raw/gdc"

mkdir -p logs
set -o pipefail

MPLBACKEND=Agg \
PYTHONUNBUFFERED=1 \
DEPMAP_DIR="$DEPMAP_DIR" \
GDC_DIR="$GDC_DIR" \
bash scripts/run_resume_functional_resilient.sh \
  2>&1 | tee logs/run_resume_functional_resilient_v0103.log

status=${PIPESTATUS[0]}
echo "$status" > logs/run_resume_functional_resilient_v0103.exitcode
echo "Resilient resume exit code: $status"
test "$status" -eq 0
```

The workflow reuses completed STRING caches, attempts OmniPath, and falls back only
when necessary to the official pinned DoRothEA human regulon from
`saezlab/dorothea-py` commit
`833165d3c790ced3a3e3852899e93412c63f0f44`. It retains confidence levels A, B and C
and records the exact fallback commit and observed SHA-256 hash.

After code 0, run:

```bash
export PIPELINE_EXITCODE_FILE="logs/run_resume_functional_resilient_v0103.exitcode"

MPLBACKEND=Agg \
GDC_DIR="$GDC_DIR" \
PIPELINE_EXITCODE_FILE="$PIPELINE_EXITCODE_FILE" \
bash scripts/verify_complete_article_run.sh \
  2>&1 | tee logs/verify_complete_article_run_v0103.log
```

The publication target remains 8 main figures, 32 supplementary figures, 120
PNG/PDF/SVG exports, 4 main tables, 18 supplementary tables and at least 140
individual structure renders. Manual inspection of all 40 figures at 100% zoom is
still required.
