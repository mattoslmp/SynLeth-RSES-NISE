# DoRothEA/OmniPath acquisition and recovery

This document records the RSES-Onco v0.10.3 recovery procedure when the official
OmniPath interaction service returns persistent HTTP 5xx errors while downloading
DoRothEA transcription-factor target interactions.

## Why the execution stopped

The v0.10.2 execution completed all 513 STRING queries and then requested the
DoRothEA dataset from the OmniPath `interactions` web service. Both accepted
organism parameter spellings returned HTTP 502 after all retries, so the previous
entry point stopped before HPA, UniProt and downstream scoring.

An HTTP 502 is a service/gateway failure. It is not scientific evidence that a
TF-target relation is absent.

## v0.10.3 behavior

`scripts/download_human_functional_evidence_resilient.py` now:

1. reuses a complete STRING aggregate/status table and does not repeat the 513
   partner queries;
2. reuses a valid cached DoRothEA TSV when present;
3. accepts a local standardized DoRothEA TSV through `--dorothea-file`;
4. tries both official OmniPath domains, with and without a trailing slash;
5. tries both `organisms=9606` and `organism=9606`;
6. validates that the returned table contains recognized TF and target columns;
7. writes `omnipath_dorothea_status.json` with every attempted endpoint and error;
8. preserves a persistent outage as explicit missing regulatory coverage by default;
9. never converts unavailable DoRothEA evidence to a score of zero;
10. supports strict mode when complete regulatory evidence is mandatory.

## Resume the interrupted complete pipeline

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
  2>&1 | tee logs/run_expanded_resume_functional_v0103.log

status=${PIPESTATUS[0]}
echo "$status" > logs/run_expanded_resume_functional_v0103.exitcode
echo "Resume exit code: $status"
test "$status" -eq 0
```

The pipeline continues through HPA, UniProt, DepMap-only scoring, GDC validation and
aggregation, integrated TCGA plus DepMap scoring, pharmacology, AlphaFold/PyMOL,
all 8 main and 32 supplementary figures, all tables, workbook, manifests,
checksums and tests.

## Require complete DoRothEA evidence

Set strict mode before the same command:

```bash
export DOROTHEA_STRICT=1
```

In strict mode a persistent outage remains fatal, but the status JSON is written
before exit.

## Use a local DoRothEA TSV

The local TSV must be non-empty and contain at least one recognized source column:

```text
source_genesymbol
source
tf
```

and one recognized target column:

```text
target_genesymbol
target
gene
```

Run:

```bash
export DOROTHEA_FILE=/absolute/path/to/dorothea_human_ABC.tsv
export DOROTHEA_STRICT=1
bash scripts/run_expanded_pipeline.sh resume-functional
```

The original target file path, row count and acquisition status are recorded in:

```text
data/raw/human_functional_evidence/omnipath_dorothea_status.json
```

## Inspect the regulatory source status

```bash
python - <<'PY'
import json
from pathlib import Path

path = Path(
  "data/raw/human_functional_evidence/omnipath_dorothea_status.json"
)
data = json.loads(path.read_text(encoding="utf-8"))
print(json.dumps(data, indent=2))
PY
```

Interpretation:

- `downloaded`: official OmniPath service returned a validated DoRothEA table;
- `cache`: a prior validated table was reused;
- `local_file`: a user-supplied validated table was used;
- `unavailable`: every official endpoint failed; regulatory evidence is missing,
  not negative.

The pair-level evidence table records:

```text
regulatory_source_available
regulatory_source_status
component_regulatory_network
```

When DoRothEA is unavailable, `component_regulatory_network` remains missing and the
coverage-aware score is reduced accordingly.

## Refresh after OmniPath recovers

Run the evidence entry point directly:

```bash
python -u scripts/download_human_functional_evidence_resilient.py \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --raw-dir data/raw/human_functional_evidence \
  --output data/processed/expanded_pair_functional_evidence.tsv \
  --refresh-dorothea \
  --strict-dorothea
```

Then rerun the downstream complete stage:

```bash
bash scripts/run_expanded_pipeline.sh resume-functional
```

The STRING aggregate is reused as long as its acquisition status covers the current
gene universe without persistent request failures.
