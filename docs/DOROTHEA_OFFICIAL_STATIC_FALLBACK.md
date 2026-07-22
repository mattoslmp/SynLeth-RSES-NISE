# Official OmniPath static DoRothEA fallback

## Why the REST URL was not replaced

The official live DoRothEA REST query remains:

```text
https://omnipathdb.org/interactions?datasets=dorothea
```

The RSES-Onco parameters (`genesymbols`, `dorothea_levels`, `fields` and human
organism) are valid. HTTP 502 responses indicate a temporary failure in the
OmniPath service or its upstream gateway, not an invalid DoRothEA endpoint.
Trying slash variants, the deprecated `omnipath.org` host, or both `organism`
and `organisms` does not repair a server-side outage and can delay a resume run
for many minutes.

## Official backup selected by RSES-Onco

OmnipathR provides static TSV tables specifically as a backup when the primary
server or a client computer has connectivity problems. For human DoRothEA, the
official static table is:

```text
http://no-tls.static.omnipathdb.org/resources/interactions_dorothea_9606.tsv.gz
```

The portable launcher now downloads and validates this table before functional
evidence acquisition. It writes:

```text
data/raw/human_functional_evidence/dorothea_official_static_9606.tsv
data/raw/human_functional_evidence/dorothea_official_static_9606_status.json
```

The launcher exports the validated TSV through `DOROTHEA_FILE`, so the existing
resilient functional-evidence stage uses it as a local source and does not wait
for the failing REST hosts.

## Manual acquisition

```bash
python -u scripts/download_dorothea_official_static.py \
  --output data/raw/human_functional_evidence/dorothea_official_static_9606.tsv \
  --status-output data/raw/human_functional_evidence/dorothea_official_static_9606_status.json \
  --levels A,B,C
```

## Resume the complete workflow

```bash
bash scripts/run_expanded_pipeline_portable.sh resume-functional
```

Expected beginning:

```text
Preparing official OmniPath static DoRothEA fallback...
DoRothEA official static: downloaded ... rows
DoRothEA source selected: .../dorothea_official_static_9606.tsv
Acquire STRING, DoRothEA, HPA and UniProt evidence
STRING: reused complete acquisition (... edge rows)
DoRothEA: local file (... rows)
```

## Controls

Force a new static download:

```bash
DOROTHEA_REFRESH_STATIC=1 \
  bash scripts/run_expanded_pipeline_portable.sh resume-functional
```

Provide another validated local table:

```bash
DOROTHEA_FILE=/absolute/path/dorothea.tsv \
  bash scripts/run_expanded_pipeline_portable.sh resume-functional
```

Disable automatic static acquisition and restore the older REST recovery path:

```bash
DOROTHEA_DISABLE_STATIC=1 \
  bash scripts/run_expanded_pipeline_portable.sh resume-functional
```

Strict mode fails if neither the selected local/static source nor the fallback
source is available:

```bash
DOROTHEA_STRICT=1 \
  bash scripts/run_expanded_pipeline_portable.sh resume-functional
```

## Scientific handling

The downloader validates source, target and `dorothea_level` columns, retains
confidence levels A, B and C by default, removes exact duplicate records, writes
the cache atomically, and records the archive SHA-256 plus retrieval time. The
source remains explicitly identified as the official OmniPath static DoRothEA
backup in the status JSON.
