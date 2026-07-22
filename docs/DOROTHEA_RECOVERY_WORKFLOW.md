# DoRothEA/OmniPath outage recovery and complete article execution

This document records the RSES-Onco v0.10.3 recovery procedure for an execution
that completed STRING acquisition and then stopped because the OmniPath
`/interactions?datasets=dorothea` service returned persistent HTTP 5xx responses.

## Scientific source policy

The preferred regulatory source remains the OmniPath webservice DoRothEA dataset
for human TF-target interactions at confidence levels A, B and C.

When the OmniPath service is unavailable after bounded retries, RSES-Onco uses the
pinned official `saezlab/dorothea-py` human regulon distributed by the same Saez
Laboratory organization:

```text
Repository: saezlab/dorothea-py
Commit: 833165d3c790ced3a3e3852899e93412c63f0f44
File: dorothea/data/dorothea_hs.pkl
```

The fallback is not a substitution with a different regulatory database. It is an
official DoRothEA human regulon snapshot. Only confidence levels A, B and C are
retained. Its exact URL, commit and observed SHA-256 hash are written to acquisition
metadata.

The pinned pickle is trusted only because it comes from the official repository and
an immutable commit. Arbitrary pickle files must never be supplied to this workflow.

## Why the previous execution failed

The completed STRING stage wrote all per-gene caches. The next request to:

```text
https://omnipathdb.org/interactions
```

returned HTTP 502 after repeated attempts for both historical organism parameter
spellings. The previous implementation treated the unavailable regulatory source as
fatal and stopped before HPA, UniProt, DepMap scoring, TCGA integration and article
asset generation.

## RSES-Onco v0.10.3 design

The resilient workflow:

1. reuses the 513 completed STRING per-gene caches;
2. reuses a normalized DoRothEA TSV from a previous successful run when present;
3. otherwise attempts OmniPath using bounded retries;
4. immediately switches to the pinned official DoRothEA snapshot after persistent
   OmniPath 5xx failure;
5. normalizes TF and target symbols and retains A/B/C confidence levels;
6. writes a normalized TSV and acquisition metadata atomically;
7. continues with HPA, UniProt and all downstream analyses;
8. does not repeat the completed all-target CRISPR discovery.

## Files written

```text
data/raw/human_functional_evidence/
├── string_interaction_partners.tsv
├── string_id_mapping.tsv
├── string_acquisition_status.tsv
├── string_cache/
├── omnipath_dorothea.tsv
├── dorothea_hs_official.pkl
├── dorothea_acquisition_metadata.json
├── hpa_subcellular_location.tsv
├── uniprot_reviewed_annotations.tsv
└── source_metadata.json
```

`omnipath_dorothea.tsv` is the normalized regulatory table used by the analysis. The
filename is retained for backward compatibility; `dorothea_acquisition_metadata.json`
records whether its source was OmniPath, the normalized cache or the pinned official
fallback.

## Update the repository

Do this only after the previous pipeline has terminated:

```bash
cd /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco

git status --short
git fetch origin
git checkout main
git pull --ff-only origin main

python -m pip install -e .
```

Confirm:

```bash
python - <<'PY'
from importlib.metadata import version
print(version("rses-onco"))
PY
```

Expected version:

```text
0.10.3
```

Run tests:

```bash
PYTHONDONTWRITEBYTECODE=1 \
MPLBACKEND=Agg \
python -m pytest -q -p no:cacheprovider
```

## One-command resilient recovery

Use the existing tmux session or create one only when not already inside tmux.

```bash
OLD="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE"
NEW="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010"

cd "$NEW" || exit 1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco

export DEPMAP_DIR="$OLD/data/raw/depmap"
export GDC_DIR="$OLD/data/raw/gdc"

RUN_LOG="logs/run_resume_functional_resilient_v0103.log"
EXITCODE_FILE="logs/run_resume_functional_resilient_v0103.exitcode"

mkdir -p logs
set -o pipefail

MPLBACKEND=Agg \
PYTHONUNBUFFERED=1 \
DEPMAP_DIR="$DEPMAP_DIR" \
GDC_DIR="$GDC_DIR" \
bash scripts/run_resume_functional_resilient.sh \
  2>&1 | tee "$RUN_LOG"

status=${PIPESTATUS[0]}
echo "$status" > "$EXITCODE_FILE"
echo "Resilient resume exit code: $status"
test "$status" -eq 0
```

This command continues through:

```text
resilient STRING/DoRothEA/HPA/UniProt acquisition
expanded DepMap-only scoring
GDC validation and aggregation
integrated TCGA plus DepMap scoring
pharmacology evidence and prioritization
AlphaFold structures and exact residue annotations
PyMOL whole/site renders
8 main figures
32 supplementary figures
4 main tables
18 supplementary tables
Excel workbook
source-data files and legends
figure/table manifests
provenance and SHA-256 checksums
full software test suite
```

## Expected recovery messages

STRING should primarily report cached genes, because the previous run completed all
513 queries.

DoRothEA should report one of:

```text
DoRothEA: cache (... A/B/C interactions)
DoRothEA: OmniPath (... A/B/C interactions)
DoRothEA: official pinned fallback (... A/B/C interactions)
```

All three are valid when their metadata file is present. The pinned fallback message
is expected while OmniPath remains unavailable.

## Inspect DoRothEA provenance

```bash
cat \
  data/raw/human_functional_evidence/dorothea_acquisition_metadata.json
```

For the fallback route, verify:

```text
status: fallback
fallback_used: true
fallback_repository: saezlab/dorothea-py
fallback_commit: 833165d3c790ced3a3e3852899e93412c63f0f44
row_count: greater than zero
fallback_sha256: non-empty
```

Programmatic check:

```bash
python - <<'PY'
import json
from pathlib import Path

path = Path(
  "data/raw/human_functional_evidence/"
  "dorothea_acquisition_metadata.json"
)
data = json.loads(path.read_text(encoding="utf-8"))
print(json.dumps(data, indent=2, sort_keys=True))
assert data["status"] in {"ok", "cache", "fallback"}
assert int(data["row_count"]) > 0
if data["status"] == "fallback":
  assert data["fallback_used"] is True
  assert data["fallback_commit"] == (
    "833165d3c790ced3a3e3852899e93412c63f0f44"
  )
  assert data["fallback_sha256"]
print("DoRothEA acquisition provenance: OK")
PY
```

## Monitor the run

```bash
tail -f \
  /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010/logs/run_resume_functional_resilient_v0103.log
```

## Final validation

After the resilient resume returns code 0:

```bash
OLD="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE"
NEW="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010"

cd "$NEW" || exit 1
conda activate rses-onco

export GDC_DIR="$OLD/data/raw/gdc"
export PIPELINE_EXITCODE_FILE="logs/run_resume_functional_resilient_v0103.exitcode"

set -o pipefail
MPLBACKEND=Agg \
GDC_DIR="$GDC_DIR" \
PIPELINE_EXITCODE_FILE="$PIPELINE_EXITCODE_FILE" \
bash scripts/verify_complete_article_run.sh \
  2>&1 | tee logs/verify_complete_article_run_v0103.log

status=${PIPESTATUS[0]}
echo "Final verification exit code: $status"
test "$status" -eq 0
```

The final validator requires:

```text
8 main figures
32 supplementary figures
120 PNG/PDF/SVG exports
4 main tables
18 supplementary tables
at least 140 individual structural renders
all registered layout audits passed
publication checksums passed
expanded-result checksums passed
all tests passed
```

Manual inspection of all 40 figures at 100% zoom remains mandatory.

## Run only the resilient functional-evidence stage

For diagnosis or cache preparation:

```bash
set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/download_human_functional_evidence_resilient.py \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --raw-dir data/raw/human_functional_evidence \
  --output data/processed/expanded_pair_functional_evidence.tsv \
  --string-required-score 700 \
  --string-limit 100 \
  --string-sleep 1.0 \
  --string-retries 7 \
  --string-map-chunk-size 200 \
  --strict-string-requests \
  2>&1 | tee logs/expanded_26Q1/08_human_functional_evidence_resilient.log

status=${PIPESTATUS[0]}
echo "Functional-evidence exit code: $status"
test "$status" -eq 0
```

Environment controls:

```text
RSES_DOROTHEA_RETRIES=4              default bounded OmniPath/fallback retries
RSES_REFRESH_DOROTHEA=1              ignore normalized and pickle caches
RSES_DISABLE_DOROTHEA_FALLBACK=1     require live OmniPath and fail otherwise
```

Do not set `RSES_REFRESH_DOROTHEA=1` during normal recovery.
