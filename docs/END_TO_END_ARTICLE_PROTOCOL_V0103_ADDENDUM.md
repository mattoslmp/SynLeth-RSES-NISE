# RSES-Onco v0.10.3 protocol addendum: OmniPath/DoRothEA outage

This addendum supplements `docs/END_TO_END_ARTICLE_PROTOCOL.md`.

When a completed STRING acquisition is followed by persistent HTTP 502 responses
from the OmniPath DoRothEA service, update to v0.10.3 and resume from functional
evidence. Do not repeat Ensembl or the all-target CRISPR screen.

```bash
cd /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco

git fetch origin
git checkout main
git pull --ff-only origin main
python -m pip install -e .

PYTHONDONTWRITEBYTECODE=1 \
MPLBACKEND=Agg \
python -m pytest -q -p no:cacheprovider
```

Resume:

```bash
OLD="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE"
NEW="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010"

cd "$NEW" || exit 1
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

The v0.10.3 functional stage reuses the completed STRING acquisition. It tries a
validated DoRothEA cache, optional local TSV and multiple official OmniPath service
addresses. If all services remain unavailable, default mode records the regulatory
source as unavailable and continues with missing regulatory coverage. It does not
turn missing evidence into zero.

For a submission run that requires complete DoRothEA coverage:

```bash
export DOROTHEA_STRICT=1
bash scripts/run_expanded_pipeline.sh resume-functional
```

For a validated local table:

```bash
export DOROTHEA_FILE=/absolute/path/to/dorothea_human_ABC.tsv
export DOROTHEA_STRICT=1
bash scripts/run_expanded_pipeline.sh resume-functional
```

After successful completion, validate using the new exit-code file:

```bash
export PIPELINE_EXITCODE_FILE=logs/run_expanded_resume_functional_v0103.exitcode
export GDC_DIR="$OLD/data/raw/gdc"

MPLBACKEND=Agg \
GDC_DIR="$GDC_DIR" \
PIPELINE_EXITCODE_FILE="$PIPELINE_EXITCODE_FILE" \
bash scripts/verify_complete_article_run.sh \
  2>&1 | tee logs/verify_complete_article_run.log
```

Inspect the source status:

```bash
cat data/raw/human_functional_evidence/omnipath_dorothea_status.json
```

A status of `unavailable` must be disclosed as missing regulatory-network coverage
in the methods, results limitations and data provenance. A later successful refresh
requires rerunning the downstream scores, tables and all figures.
