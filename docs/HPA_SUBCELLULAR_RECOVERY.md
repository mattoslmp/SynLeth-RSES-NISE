# Human Protein Atlas subcellular-data recovery

The Human Protein Atlas current downloadable-data page uses:

```text
https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip
```

The previous RSES-Onco URL omitted the `/tsv/` path segment and now returns HTTP
404. The current official archive contains the per-gene subcellular localization
columns required by the RSES-Onco localization-divergence component.

## Recovery after an HPA HTTP 404

From the repository root:

```bash
conda activate rses-onco

python -u scripts/download_hpa_subcellular_current.py \
  --output data/raw/human_functional_evidence/hpa_subcellular_location.tsv \
  --status-output data/raw/human_functional_evidence/hpa_subcellular_status.json
```

Expected output:

```text
HPA: downloaded ... rows from https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip
HPA: wrote .../data/raw/human_functional_evidence/hpa_subcellular_location.tsv
```

Validate the cached table:

```bash
python - <<'PY'
import pandas as pd

path = "data/raw/human_functional_evidence/hpa_subcellular_location.tsv"
frame = pd.read_csv(path, sep="\t")
required = {
  "Gene",
  "Gene name",
  "Reliability",
  "Enhanced",
  "Supported",
  "Approved",
  "Uncertain",
}
missing = sorted(required - set(frame.columns))
print("Rows:", len(frame))
print("Missing required columns:", missing)
assert not frame.empty
assert not missing
print("HPA cache: OK")
PY
```

Then resume through the portable launcher:

```bash
bash scripts/run_expanded_pipeline_portable.sh resume-functional
```

The resilient functional-evidence script detects the cached HPA TSV and does not
request the obsolete URL. STRING remains cached and is not repeated.
