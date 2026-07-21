# RSES-Onco v0.10.2 addendum to the end-to-end article protocol

This addendum supplements `docs/END_TO_END_ARTICLE_PROTOCOL.md` for executions that
reach the functional-evidence stage after completing all-target DepMap discovery.

RSES-Onco v0.10.2 replaces direct symbol-level STRING partner queries with:

1. stable-version discovery;
2. batched mapping to exact STRING identifiers;
3. per-gene interaction-partner queries;
4. per-gene resumable caches;
5. explicit mapping and acquisition-status tables;
6. strict failure only for persistent request errors, not for genuine unmapped
   identifiers.

The complete implementation and recovery commands are documented in:

```text
docs/STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md
```

After a previous run stopped in `download_human_functional_evidence.py`, update to
v0.10.2 and resume without repeating the completed all-target CRISPR screen:

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
bash scripts/run_expanded_pipeline.sh resume-functional \
  2>&1 | tee logs/run_expanded_resume_functional_v0102.log

status=${PIPESTATUS[0]}
echo "$status" > logs/run_expanded_resume_functional_v0102.exitcode
echo "Resume pipeline exit code: $status"
test "$status" -eq 0
```

After a successful resume, use the same post-run verifier documented in the
canonical protocol, pointing it to the resume exit-code file:

```bash
export PIPELINE_EXITCODE_FILE="logs/run_expanded_resume_functional_v0102.exitcode"

MPLBACKEND=Agg \
GDC_DIR="$GDC_DIR" \
PIPELINE_EXITCODE_FILE="$PIPELINE_EXITCODE_FILE" \
bash scripts/verify_complete_article_run.sh \
  2>&1 | tee logs/verify_complete_article_run.log
```

The publication target remains 8 main figures, 32 supplementary figures, 120
PNG/PDF/SVG exports, 4 main tables, 18 supplementary tables and at least 140
individual structural renders. Manual inspection of all 40 figures at 100% zoom
remains mandatory.
