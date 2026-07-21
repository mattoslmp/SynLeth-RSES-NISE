# RSES-Onco v0.10.1: end-to-end real-data and article-generation protocol

This is the canonical execution protocol for acquiring the required public data,
validating every input, running the complete RSES-Onco analysis, generating every
main and supplementary figure, building all tables and workbooks, and validating
the final publication package.

The commands are designed for Linux or WSL. Run them from the repository root.
The example paths below match the validated WSL execution used during development:

```text
Repository: /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010
Existing large-data directory: /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE
Conda environment: rses-onco
```

A fresh installation may keep all data inside the repository. The `DEPMAP_DIR` and
`GDC_DIR` environment variables permit large source data to remain outside the Git
working tree.

## 1. What the complete protocol produces

A successful run performs:

1. validation and fingerprinting of the required DepMap matrices;
2. construction of all 202 directed human NISE hypotheses;
3. preservation of curated composite biomarkers without inventing single-gene loss calls;
4. expansion of explicit multi-target fields such as `CDK4/CDK6` into traceable atomic targets;
5. resumable Ensembl Compara human-paralog acquisition;
6. all-target DepMap conditional-dependency discovery;
7. STRING, DoRothEA, Human Protein Atlas and UniProt functional evidence acquisition;
8. DepMap-only scoring;
9. validation and aggregation of TCGA/GDC ASCAT3 gene-level copy-number files;
10. integrated TCGA plus DepMap scoring for colorectal, gastric and lung cancer;
11. Open Targets, ChEMBL, DGIdb, MyChem, Pharos/TCRD and CIViC evidence acquisition;
12. optional PRISM, GDSC and CTRP drug-response standardization and analysis;
13. AlphaFold DB acquisition for all 70 curated human NISE proteins;
14. exact-numbered M-CSA, UniProt and PDBe residue annotation;
15. PyMOL generation of whole-structure and functional-site renders;
16. generation of all article tables, figures, legends, source-data files and layout audits;
17. generation of the article workbook, provenance inventory and SHA-256 manifests;
18. final publication-package validation and the complete software test suite.

The validated publication target is:

```text
8 main figures
32 supplementary figures
40 registered figures total
120 figure exports: PNG, PDF and SVG
4 main tables
18 supplementary tables
at least 140 individual AlphaFold/PyMOL PNG renders
```

## 2. Clone or update the repository

### Fresh clone

```bash
git clone https://github.com/mattoslmp/SynLeth-RSES-Onco.git \
  /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010

cd /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010
```

### Existing clone

Do not update files while a pipeline from the same working tree is running. After
that process finishes:

```bash
cd /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010

git status --short
git fetch origin
git checkout main
git pull --ff-only origin main

git log -1 --oneline
```

Record the exact commit used for a release-specific analysis:

```bash
git rev-parse HEAD | tee logs/analysis_git_commit.txt
```

## 3. Create and validate the software environment

### New environment

```bash
cd /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010

conda env create -f environment.yml
conda activate rses-onco
python -m pip install -e .
```

### Existing environment

```bash
cd /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco
python -m pip install -e .
```

Confirm the installed version and PyMOL availability:

```bash
python - <<'PY'
from importlib.metadata import version
print("RSES-Onco:", version("rses-onco"))
PY

pymol -cq <<'PYMOL'
reinitialize
print "PyMOL headless execution: OK"
quit
PYMOL
```

Run the complete test suite before using real data:

```bash
PYTHONDONTWRITEBYTECODE=1 \
MPLBACKEND=Agg \
python -m pytest -q -p no:cacheprovider
```

The command must return exit code `0`.

## 4. Prepare directories and environment variables

### Validated WSL layout with large data in the previous directory

```bash
OLD="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE"
NEW="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010"

cd "$NEW" || exit 1

export DEPMAP_DIR="$OLD/data/raw/depmap"
export GDC_DIR="$OLD/data/raw/gdc"

mkdir -p \
  "$DEPMAP_DIR" \
  "$GDC_DIR" \
  data/raw/ensembl \
  data/raw/pharmacology \
  data/raw/structures \
  data/processed \
  results \
  article_outputs \
  logs
```

### Self-contained layout

```bash
cd /path/to/SynLeth-RSES-Onco

export DEPMAP_DIR="$PWD/data/raw/depmap"
export GDC_DIR="$PWD/data/raw/gdc"

mkdir -p "$DEPMAP_DIR" "$GDC_DIR" logs
```

Keep these variables exported in the same shell used to launch the pipeline.

## 5. Acquire the required DepMap release files

The repository validates DepMap files but does not bypass the official DepMap
download process. Download the current release through the official DepMap data
portal and place these four required files in `$DEPMAP_DIR`:

```text
CRISPRGeneEffect.csv
OmicsCNGeneWGS.csv
Model.csv
OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
```

The validator also accepts the documented legacy or compressed filenames listed in
`scripts/download_depmap.py`.

Recommended release metadata files may be archived beside the matrices, but they do
not replace the four required inputs.

Validate file presence, schemas and SHA-256 fingerprints:

```bash
mkdir -p logs/expanded_26Q1

set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/download_depmap.py \
  --input-dir "$DEPMAP_DIR" \
  --write-checksums \
  2>&1 | tee logs/expanded_26Q1/01_validate_depmap_files.log

status=${PIPESTATUS[0]}
echo "DepMap validation exit code: $status"
test "$status" -eq 0
```

Verify the generated checksums:

```bash
(
  cd "$DEPMAP_DIR"
  sha256sum -c SHA256SUMS.txt
)
```

Validate model identifiers, matrix crosswalks and cancer cohorts:

```bash
set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/validate_real_inputs.py \
  --gene-effect "$DEPMAP_DIR/CRISPRGeneEffect.csv" \
  --copy-number "$DEPMAP_DIR/OmicsCNGeneWGS.csv" \
  --models "$DEPMAP_DIR/Model.csv" \
  --expression "$DEPMAP_DIR/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv" \
  2>&1 | tee logs/expanded_26Q1/02_validate_depmap_crosswalk.log

status=${PIPESTATUS[0]}
echo "DepMap crosswalk exit code: $status"
test "$status" -eq 0
```

Release-specific model counts may change. They must be recorded, not forced to match
an older release.

## 6. Optional local drug-sensitivity releases

Open Targets, ChEMBL, DGIdb, MyChem, Pharos/TCRD and CIViC are queried by the
publication pipeline. PRISM, GDSC and CTRP are optional local releases.

Place optional files under:

```text
data/raw/pharmacology/prism/
data/raw/pharmacology/gdsc/
data/raw/pharmacology/ctrp/
```

Accepted filename patterns and column aliases are defined in:

```text
config/drug_sensitivity_sources.yaml
```

Absence of these optional local releases must remain explicit. The pipeline writes a
source-status table and does not invent drug-response measurements.

## 7. Create and review the GDC manifest

The target projects are:

```text
TCGA-COAD
TCGA-READ
TCGA-STAD
TCGA-LUAD
TCGA-LUSC
```

The workflow requests open, primary-tumor, gene-level copy-number files produced by
ASCAT3.

Create the manifest without downloading files:

```bash
mkdir -p "$GDC_DIR" logs/expanded_26Q1

set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/download_gdc.py \
  --manifest-only \
  --workflow ASCAT3 \
  --output-dir "$GDC_DIR" \
  2>&1 | tee logs/expanded_26Q1/03_gdc_manifest.log

status=${PIPESTATUS[0]}
echo "GDC manifest exit code: $status"
test "$status" -eq 0
```

Preserve a dated copy before downloading:

```bash
cp -av \
  "$GDC_DIR/gdc_gene_level_copy_number_manifest.json" \
  "$GDC_DIR/gdc_gene_level_copy_number_manifest_$(date +%Y%m%d).json"
```

Summarize the reviewed manifest:

```bash
python - "$GDC_DIR/gdc_gene_level_copy_number_manifest.json" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
data = json.loads(path.read_text(encoding="utf-8"))
total_files = 0
total_bytes = 0
for project, records in data.items():
  count = len(records)
  size = sum(int(record.get("file_size") or 0) for record in records)
  total_files += count
  total_bytes += size
  print(f"{project}: {count} files; {size / 1024**3:.3f} GiB")
print(f"Total: {total_files} files; {total_bytes / 1024**3:.3f} GiB")
PY
```

For the manifest reviewed on 21 July 2026, the counts were 422 COAD, 153 READ,
429 STAD, 503 LUAD and 490 LUSC, for 1,997 files. Later GDC releases may differ and
must be reported using their actual reviewed manifest.

## 8. Download or resume the reviewed GDC manifest

Use the existing reviewed manifest rather than silently querying a different cohort
mid-download:

```bash
set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/download_gdc.py \
  --use-existing-manifest \
  --manifest "$GDC_DIR/gdc_gene_level_copy_number_manifest.json" \
  --workflow ASCAT3 \
  --output-dir "$GDC_DIR" \
  --retries 3 \
  2>&1 | tee logs/expanded_26Q1/04_gdc_download.log

status=${PIPESTATUS[0]}
echo "GDC download exit code: $status"
test "$status" -eq 0
```

The downloader skips valid files, retries failed requests, writes temporary `.part`
files and validates size plus MD5 before accepting each file.

Do not start a second downloader while this process is active.

Monitor from a second terminal:

```bash
watch -n 15 "
  echo -n 'Complete files: '
  find '$GDC_DIR'/TCGA-* -type f ! -name '*.part' 2>/dev/null | wc -l
  echo -n 'Partial files:  '
  find '$GDC_DIR'/TCGA-* -type f -name '*.part' 2>/dev/null | wc -l
  echo -n 'Current volume:  '
  du -sh '$GDC_DIR' 2>/dev/null
"
```

Validate the completed download:

```bash
set -o pipefail
PYTHONUNBUFFERED=1 \
python -u scripts/download_gdc.py \
  --validate-only \
  --manifest "$GDC_DIR/gdc_gene_level_copy_number_manifest.json" \
  --output-dir "$GDC_DIR" \
  2>&1 | tee logs/expanded_26Q1/05_gdc_validate.log

status=${PIPESTATUS[0]}
echo "GDC validation exit code: $status"
test "$status" -eq 0
```

Confirm that no partial file remains:

```bash
find "$GDC_DIR" -type f -name '*.part'
```

The command must print nothing.

## 9. Start the complete analysis after the reviewed GDC download

Use `tmux` for the long execution. Do not nest a new session when `$TMUX` is already
set.

Outside tmux:

```bash
tmux new-session -A -s rses_v010_full
```

Already inside tmux:

```bash
tmux display-message -p 'Session: #S | Window: #I:#W | Pane: #P'
```

Before launching, verify that no other copy of the pipeline is active:

```bash
pgrep -af \
  'run_expanded_pipeline|discover_conditional_dependencies|run_expanded_rses_onco|run_publication_pipeline|render_nise_structures' \
  || true
```

Launch the canonical complete workflow:

```bash
OLD="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE"
NEW="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010"

cd "$NEW" || exit 1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco

export DEPMAP_DIR="$OLD/data/raw/depmap"
export GDC_DIR="$OLD/data/raw/gdc"

mkdir -p logs

RUN_LOG="logs/run_expanded_after_download_v0101.log"
EXITCODE_FILE="logs/run_expanded_after_download_v0101.exitcode"

printf 'Project: %s\n' "$PWD"
printf 'Conda environment: %s\n' "$CONDA_DEFAULT_ENV"
printf 'DepMap: %s\n' "$DEPMAP_DIR"
printf 'GDC: %s\n' "$GDC_DIR"
printf 'Start: %s\n' "$(date -Iseconds)"

set -o pipefail

MPLBACKEND=Agg \
PYTHONUNBUFFERED=1 \
DEPMAP_DIR="$DEPMAP_DIR" \
GDC_DIR="$GDC_DIR" \
bash scripts/run_expanded_pipeline.sh after-download \
  2>&1 | tee "$RUN_LOG"

status=${PIPESTATUS[0]}
printf '%s\n' "$status" > "$EXITCODE_FILE"
printf 'Finish: %s\n' "$(date -Iseconds)"
printf 'Pipeline exit code: %s\n' "$status"

test "$status" -eq 0
```

To detach without stopping the run, press `Ctrl+B`, release, then press `D`.

Monitor from another terminal:

```bash
tail -f \
  /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010/logs/run_expanded_after_download_v0101.log
```

Return to the session with:

```bash
tmux attach -t rses_v010_full
```

## 10. Complete pipeline stage order

`scripts/run_expanded_pipeline.sh after-download` runs these stages:

```text
01  Validate DepMap files and write checksums
02  Validate DepMap model crosswalks and cancer cohorts
03  Build the curated directed NISE and benchmark universe
04  Acquire/cache Ensembl Compara paralogs
05  Rebuild the universe with Ensembl paralogs
06  Screen every measured DepMap CRISPR target
07  Rebuild the universe with supported discoveries
08  Acquire STRING, DoRothEA, HPA and UniProt evidence
09  Run expanded DepMap-only scoring
10  Validate all GDC files against the reviewed manifest
11  Aggregate ASCAT3 gene-level copy-number files
12  Validate combined colorectal, gastric and lung matrices
13  Run integrated TCGA plus DepMap scoring
14  Run pharmacology, structures and publication assets
15  Re-run tests and write expanded-result checksums
```

The publication subpipeline then runs:

```text
01  Acquire pharmacology evidence
02  Standardize optional PRISM/GDSC/CTRP releases
03  Analyze biomarker-matched drug response
04  Build coverage-aware target-drug priorities
05  Download all 70 AlphaFold structures
06  Collect exact-numbered structural annotations
07  Render whole and functional-site PyMOL views
08  Export all main and supplementary tables
09  Generate all main, supplementary and structural figures
10  Build the organized Excel workbook
11  Build inventories, provenance and SHA-256 manifests
12  Validate all publication outputs and run tests
```

## 11. Ensembl completeness requirement

The Ensembl downloader caches each seed-gene response and resolves target identifiers
in batches. Explicit target lists are expanded into atomic symbols before acquisition.

After the Ensembl stage, verify:

```bash
python - <<'PY'
import json
from pathlib import Path

path = Path("data/raw/ensembl/ensembl_acquisition_metadata.json")
data = json.loads(path.read_text(encoding="utf-8"))
for key in (
  "seed_gene_count",
  "successful_seed_homology_queries",
  "failed_seed_homology_queries",
  "target_identifier_count",
  "resolved_target_identifier_count",
  "unresolved_target_identifier_count",
  "directed_paralog_count",
  "complete",
):
  print(f"{key}: {data.get(key)}")
PY
```

A complete run requires:

```text
failed_seed_homology_queries: 0
unresolved_target_identifier_count: 0
complete: True
```

The validated v0.10.1 execution contained 95 seed genes, 471 resolved target
identifiers and 719 directed Ensembl paralog candidates. These values are
release-specific and may change when the curated inputs or Ensembl release changes.

## 12. Exact commands to execute after the complete pipeline finishes

Do not assume that a returned prompt means success. Execute every command below.

### 12.1 Read the recorded exit code and final log

```bash
cd /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010
conda activate rses-onco

cat logs/run_expanded_after_download_v0101.exitcode

tail -n 120 logs/run_expanded_after_download_v0101.log
```

The recorded exit code must be exactly:

```text
0
```

Search for failures:

```bash
grep -nEi \
  'traceback|error|exception|failed|failure|fatal|command failed|no such file' \
  logs/run_expanded_after_download_v0101.log \
  | tail -n 100 \
  || true
```

Warnings from optional unavailable sources must be interpreted using their source
status tables; they must not be silently treated as observed evidence.

### 12.2 Run the canonical final verifier

For the validated WSL layout:

```bash
OLD="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE"
NEW="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010"

cd "$NEW" || exit 1
conda activate rses-onco

export GDC_DIR="$OLD/data/raw/gdc"
export PIPELINE_EXITCODE_FILE="logs/run_expanded_after_download_v0101.exitcode"

mkdir -p logs
set -o pipefail

MPLBACKEND=Agg \
GDC_DIR="$GDC_DIR" \
PIPELINE_EXITCODE_FILE="$PIPELINE_EXITCODE_FILE" \
bash scripts/verify_complete_article_run.sh \
  2>&1 | tee logs/verify_complete_article_run.log

status=${PIPESTATUS[0]}
echo "Final verification exit code: $status"
test "$status" -eq 0
```

This command verifies:

- the recorded pipeline exit code;
- Ensembl completeness;
- absence of GDC `.part` files;
- all figure/table/manuscript assets;
- publication and expanded-result SHA-256 checksums;
- the full test suite;
- 8 main and 32 supplementary figures;
- 120 PNG/PDF/SVG figure exports;
- 4 main and 18 supplementary tables;
- at least 140 individual structural renders.

### 12.3 Independently validate publication outputs

The helper already performs this command, but it is retained here for transparent
manual reproduction:

```bash
python -u scripts/validate_publication_outputs.py \
  --article-root article_outputs
```

Expected final messages:

```text
Publication package validation passed.
Main figures: 8; supplementary figures: 32; exported image files: 120
Main tables: 4; supplementary tables: 18
All 70 curated NISE proteins have whole and site structural renders.
All registered figures passed automated layout audits.
```

### 12.4 Independently verify checksums

Publication assets:

```bash
(
  cd article_outputs
  sha256sum -c manifests/SHA256SUMS.txt
)
```

Integrated analysis outputs:

```bash
sha256sum -c \
  results/expanded_26Q1/full/SHA256SUMS.txt
```

Every listed file must report `OK`.

### 12.5 Inspect the output inventory

```bash
column -t -s $'\t' \
  article_outputs/manifests/figure_manifest.tsv \
  | less -S
```

```bash
column -t -s $'\t' \
  article_outputs/manifests/table_manifest.tsv \
  | less -S
```

```bash
python - <<'PY'
import pandas as pd

figures = pd.read_csv(
  "article_outputs/manifests/figure_manifest.tsv",
  sep="\t",
)
tables = pd.read_csv(
  "article_outputs/manifests/table_manifest.tsv",
  sep="\t",
)

print("Registered figures:", len(figures))
print("Layout pass:", int(figures["layout_status"].eq("pass").sum()))
print("Registered tables:", len(tables))
print("Main tables:", int(tables["category"].eq("main").sum()))
print(
  "Supplementary tables:",
  int(tables["category"].eq("supplementary").sum()),
)
PY
```

### 12.6 Record final provenance

```bash
cat article_outputs/manifests/publication_provenance.json
```

```bash
git rev-parse HEAD

git status --short
```

The generated provenance records the commit, branch, Python/platform information,
key package versions and working-tree status. A dirty working tree must be explained
before archiving or publication.

## 13. Mandatory manual inspection at 100% zoom

Automated layout validation is necessary but does not replace scientific and visual
review. Inspect every main and supplementary PDF at 100% zoom.

Open the directories:

```bash
explorer.exe "$(wslpath -w "$PWD/article_outputs/figures/main")"
explorer.exe "$(wslpath -w "$PWD/article_outputs/figures/supplementary")"
```

For every figure verify:

- all panel letters are present and correctly ordered;
- axis labels, tick labels, sample names and gene names are readable at 100% zoom;
- legends and color bars are readable and are not clipped;
- heatmap values, dendrogram labels and statistical annotations are legible;
- no text overlaps another text element, axis, panel or structure;
- no scientific label has been truncated;
- PNG, PDF and SVG show the same scientific content;
- the displayed values agree with the associated source-data TSV;
- the figure legend accurately describes every panel and symbol;
- missing evidence is displayed as missing and not as zero;
- AlphaFold figures do not imply experimental ligand poses that were not observed.

Record the inspection in a dated text file:

```bash
cat > article_outputs/manifests/MANUAL_VISUAL_INSPECTION.txt <<EOF
RSES-Onco manual visual inspection
Date: $(date -Iseconds)
Git commit: $(git rev-parse HEAD)
Zoom: 100%
Main figures inspected: Figure 1 through Figure 8
Supplementary figures inspected: Figure S1 through Figure S32
PNG/PDF/SVG correspondence checked: yes
Source-data correspondence checked: yes
Clipping/overlap detected: no
Reviewer: Leandro de Mattos Pereira
EOF
```

Edit the statements before saving when any problem is detected. Do not mark the
package as submission-ready until all problems are corrected by the generating
scripts and the complete publication workflow is rerun.

Because this manual record is created after the generated checksum manifest, rebuild
the inventory and checksums after adding it:

```bash
RANKING="results/expanded_26Q1/full/expanded_rses_onco.tsv" \
CANDIDATES="data/processed/expanded_candidate_universe.tsv" \
DISCOVERY="results/expanded_26Q1/discovery/all_target_dependency_screen.tsv" \
FUNCTIONAL_EVIDENCE="data/processed/expanded_pair_functional_evidence.tsv" \
PHARMACOLOGY_DATA="data/processed/pharmacology" \
PHARMACOLOGY_RESULTS="results/expanded_26Q1/pharmacology" \
STRUCTURE_MANIFEST="data/processed/structures/alphafold_structure_manifest.tsv" \
STRUCTURAL_ANNOTATIONS="data/processed/structures/nise_structural_residue_annotations.tsv" \
STRUCTURE_RENDER_MANIFEST="data/processed/structures/nise_structure_render_manifest.tsv" \
bash scripts/run_publication_pipeline.sh manifests
```

Then re-run:

```bash
MPLBACKEND=Agg \
GDC_DIR="$GDC_DIR" \
bash scripts/verify_complete_article_run.sh
```

## 14. Regenerate only publication assets from completed scores and caches

When the integrated ranking, pharmacology caches and structural renders already
exist, rebuild tables, figures, workbook and manifests with:

```bash
MPLBACKEND=Agg \
PYTHONUNBUFFERED=1 \
bash scripts/run_publication_pipeline.sh assets-only \
  2>&1 | tee logs/run_publication_assets_only.log
```

Regenerate only the 40 registered figures:

```bash
MPLBACKEND=Agg \
PYTHONUNBUFFERED=1 \
bash scripts/run_publication_pipeline.sh figures \
  2>&1 | tee logs/regenerate_all_article_figures.log
```

Run the publication validator afterward:

```bash
python -u scripts/validate_publication_outputs.py \
  --article-root article_outputs
```

## 15. Re-run selected stages after interruption

### Ensembl

The per-seed homology cache is stored in:

```text
data/raw/ensembl/homology_cache/
```

Re-run without `--refresh` to reuse successful responses:

```bash
python -u scripts/download_ensembl_paralogs.py \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --output data/raw/ensembl/human_seed_paralogs.tsv \
  --retries 10 \
  --batch-size 500 \
  --strict-completeness
```

### GDC

Re-run the reviewed manifest command from Section 8. Existing valid files are
skipped; invalid or partial files are replaced.

### Pharmacology

API responses are cached under:

```text
data/raw/pharmacology/api_cache/
```

Re-run:

```bash
bash scripts/run_publication_pipeline.sh acquire-pharmacology
bash scripts/run_publication_pipeline.sh prioritize
```

### Structures

AlphaFold structures and annotation caches are stored under:

```text
data/raw/structures/
data/processed/structures/
```

Re-run selected structural stages:

```bash
bash scripts/run_structural_pipeline.sh download
bash scripts/run_structural_pipeline.sh annotations
bash scripts/run_structural_pipeline.sh render
bash scripts/run_structural_pipeline.sh figures
```

## 16. Output locations

```text
results/expanded_26Q1/
├── depmap_only/
├── discovery/
├── full/
└── pharmacology/

article_outputs/
├── figures/main/
├── figures/supplementary/
├── structure_atlas/individual/
├── tables/main/
├── tables/supplementary/
├── source_data/
├── manuscript_assets/
├── workbooks/
└── manifests/
```

Important final files include:

```text
results/expanded_26Q1/full/expanded_rses_onco.tsv
results/expanded_26Q1/discovery/all_target_dependency_screen.tsv
results/expanded_26Q1/pharmacology/pharmacology_ranked_hypotheses.tsv
article_outputs/workbooks/RSES_Onco_Article_Tables_and_Evidence.xlsx
article_outputs/manuscript_assets/all_figure_legends.md
article_outputs/manifests/figure_manifest.tsv
article_outputs/manifests/table_manifest.tsv
article_outputs/manifests/publication_file_inventory.tsv
article_outputs/manifests/publication_provenance.json
article_outputs/manifests/SHA256SUMS.txt
```

## 17. Minimal command summary

After the four required DepMap files and the reviewed GDC download are present:

```bash
export DEPMAP_DIR=/absolute/path/to/data/raw/depmap
export GDC_DIR=/absolute/path/to/data/raw/gdc

MPLBACKEND=Agg \
PYTHONUNBUFFERED=1 \
DEPMAP_DIR="$DEPMAP_DIR" \
GDC_DIR="$GDC_DIR" \
bash scripts/run_expanded_pipeline.sh after-download
```

After it finishes:

```bash
export PIPELINE_EXITCODE_FILE=logs/run_expanded_after_download_v0101.exitcode

MPLBACKEND=Agg \
GDC_DIR="$GDC_DIR" \
PIPELINE_EXITCODE_FILE="$PIPELINE_EXITCODE_FILE" \
bash scripts/verify_complete_article_run.sh
```

The second command does not replace manual inspection of all 40 figures at 100%
zoom.
