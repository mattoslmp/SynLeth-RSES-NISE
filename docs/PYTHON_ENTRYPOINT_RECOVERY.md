# Recovery from `ModuleNotFoundError: No module named 'scripts'`

## Symptom

A resumed run may stop immediately at the resilient functional-evidence entry point:

```text
ModuleNotFoundError: No module named 'scripts'
```

This is an execution-path problem, not a scientific-data failure. When Python runs a
file such as `python scripts/download_human_functional_evidence_resilient.py`, it
places the `scripts/` directory at `sys.path[0]`. The repository root is therefore
not guaranteed to be available for imports using the `scripts.*` namespace.

## Canonical portable launcher

Use the repository launcher, which exports the repository root through `PYTHONPATH`
before delegating to the complete pipeline:

```bash
bash scripts/run_expanded_pipeline_portable.sh resume-functional
```

Validate imports before starting a long run:

```bash
bash scripts/run_expanded_pipeline_portable.sh --check-python-imports
```

Expected output:

```text
RSES-Onco script namespace imports: OK
```

## Validated WSL recovery command

```bash
OLD="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE"
NEW="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010"

cd "$NEW" || exit 1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco

export DEPMAP_DIR="$OLD/data/raw/depmap"
export GDC_DIR="$OLD/data/raw/gdc"
unset DOROTHEA_STRICT
unset DOROTHEA_FILE

RUN_LOG="logs/run_expanded_resume_functional_v0103_portable.log"
EXITCODE_FILE="logs/run_expanded_resume_functional_v0103_portable.exitcode"
mkdir -p logs
set -o pipefail

MPLBACKEND=Agg \
PYTHONUNBUFFERED=1 \
DEPMAP_DIR="$DEPMAP_DIR" \
GDC_DIR="$GDC_DIR" \
bash scripts/run_expanded_pipeline_portable.sh resume-functional \
  2>&1 | tee "$RUN_LOG"

status=${PIPESTATUS[0]}
echo "$status" > "$EXITCODE_FILE"
echo "Portable resume exit code: $status"
test "$status" -eq 0
```

The complete STRING acquisition is reused when its aggregate and status files cover
all current candidate genes without persistent request failures. DoRothEA follows
the resilient cache/local-file/multiple-endpoint behavior documented in
`docs/DOROTHEA_RECOVERY_WORKFLOW.md`.
