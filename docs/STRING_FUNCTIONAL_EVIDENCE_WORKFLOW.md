# STRING functional-evidence acquisition and recovery

This document records the RSES-Onco v0.10.2 workflow for acquiring STRING
functional-interaction evidence and resuming the complete real-data pipeline after
an interruption at the functional-evidence stage.

## Why the previous execution failed

The earlier downloader sent one gene symbol directly to:

```text
https://string-db.org/api/tsv/interaction_partners
```

STRING documents that a request returns HTTP 404 when none of the submitted
identifiers can be resolved. A 404 therefore does not necessarily mean that the API
endpoint is absent. The previous code retried every HTTP error and then aborted the
entire analysis when one symbol-level request returned 404.

STRING also recommends:

1. mapping input identifiers to exact STRING identifiers first;
2. using mapped STRING identifiers for subsequent queries;
3. using a version-specific stable STRING address in production;
4. waiting approximately one second between calls.

RSES-Onco v0.10.2 implements all four requirements.

## Current acquisition design

`scripts/download_human_functional_evidence.py` now:

- discovers the current stable STRING version through `/api/json/version`;
- falls back to the pinned STRING 12.0 address when version discovery fails;
- maps candidate genes in bounded POST batches through `/api/tsv/get_string_ids`;
- queries `/api/tsv/interaction_partners` using exact mapped STRING IDs;
- retries only network failures, HTTP 408/425/429 and HTTP 5xx responses;
- does not repeatedly retry non-retryable HTTP 404 responses;
- writes a per-gene TSV/JSON cache after every successful query;
- reuses completed gene caches on a later run;
- records unmapped genes as missing evidence rather than zero evidence;
- writes mapping and acquisition-status tables;
- fails strict execution only when persistent request failures remain.

## Output and cache files

```text
data/raw/human_functional_evidence/
├── string_interaction_partners.tsv
├── string_id_mapping.tsv
├── string_acquisition_status.tsv
├── string_cache/
│   ├── AADAC.tsv
│   ├── AADAC.json
│   ├── ...
│   └── WRN.json
├── omnipath_dorothea.tsv
├── hpa_subcellular_location.tsv
├── uniprot_reviewed_annotations.tsv
└── source_metadata.json
```

The per-gene JSON cache records the exact stable API root, mapped STRING ID,
required score, partner limit, network type, status, number of rows and any error.
A cache is reused only when these acquisition parameters match.

## Run only the functional-evidence acquisition

From the repository root:

```bash
conda activate rses-onco

mkdir -p logs/expanded_26Q1
set -o pipefail

PYTHONUNBUFFERED=1 \
python -u scripts/download_human_functional_evidence.py \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --raw-dir data/raw/human_functional_evidence \
  --output data/processed/expanded_pair_functional_evidence.tsv \
  --string-required-score 700 \
  --string-limit 100 \
  --string-sleep 1.0 \
  --string-retries 7 \
  --string-map-chunk-size 200 \
  --strict-string-requests \
  2>&1 | tee logs/expanded_26Q1/08_human_functional_evidence.log

status=${PIPESTATUS[0]}
echo "Functional-evidence exit code: $status"
test "$status" -eq 0
```

Do not use `--refresh-string` during a normal recovery. Without that option, every
completed per-gene cache is reused and only absent or failed queries are repeated.

## Inspect STRING mapping and acquisition status

```bash
python - <<'PY'
import pandas as pd

mapping = pd.read_csv(
  "data/raw/human_functional_evidence/string_id_mapping.tsv",
  sep="\t",
)
status = pd.read_csv(
  "data/raw/human_functional_evidence/string_acquisition_status.tsv",
  sep="\t",
)

print("Mapping status:")
print(mapping["mapping_status"].value_counts(dropna=False).to_string())
print("\nAcquisition status:")
print(status["status"].value_counts(dropna=False).to_string())
print("\nPersistent request failures:")
failed = status.loc[
  status["status"].isin(["mapping_request_failed", "request_failed"])
]
print(failed.to_string(index=False) if not failed.empty else "none")
PY
```

A strict successful execution requires no `mapping_request_failed` or
`request_failed` rows. An `unmapped` row is retained as explicit missing source
coverage and is not converted to a network score of zero.

## Resume the complete pipeline after all-target discovery

Use this stage when the previous run already completed:

- candidate-universe construction;
- Ensembl paralog acquisition;
- all-target DepMap discovery;
- rebuilding the candidate universe with supported discoveries;

and then stopped while acquiring STRING/DoRothEA/HPA/UniProt evidence.

For the validated WSL layout:

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

The resume stage verifies that the expanded candidate universe, class-member
inventory, Ensembl paralog catalogue and all-target dependency screen exist. It
then executes:

```text
functional evidence
expanded DepMap-only scoring
GDC validation and aggregation
integrated TCGA plus DepMap scoring
pharmacology
AlphaFold acquisition and structural annotations
PyMOL renders
all main and supplementary tables
all 8 main and 32 supplementary figures
workbook, manifests and checksums
final tests
```

It does not repeat the completed all-target CRISPR discovery.

## Monitor a resumed run

```bash
tail -f \
  /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010/logs/run_expanded_resume_functional_v0102.log
```

The STRING stage prints one of these statuses for each gene:

```text
ok
no_interactions
unmapped
cache ok
cache no_interactions
request_failed
mapping_request_failed
```

A rerun reuses `cache ok` and `cache no_interactions` genes.

## Post-run validation

After the resumed pipeline returns exit code 0:

```bash
export GDC_DIR="$OLD/data/raw/gdc"
export PIPELINE_EXITCODE_FILE="logs/run_expanded_resume_functional_v0102.exitcode"

MPLBACKEND=Agg \
GDC_DIR="$GDC_DIR" \
PIPELINE_EXITCODE_FILE="$PIPELINE_EXITCODE_FILE" \
bash scripts/verify_complete_article_run.sh \
  2>&1 | tee logs/verify_complete_article_run.log
```

The final publication target remains:

```text
8 main figures
32 supplementary figures
120 PNG/PDF/SVG exports
4 main tables
18 supplementary tables
at least 140 individual structural renders
```

Manual inspection of all 40 figures at 100% zoom remains mandatory.
