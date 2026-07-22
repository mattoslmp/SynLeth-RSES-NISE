# SynLeth-RSES-Onco / RSES-Onco

**RSES-Onco v0.10.3** is a coverage-aware framework for discovering and
prioritizing cancer-selective dependencies created by non-homologous
isofunctional enzymes (NISEs), homologous paralogs, pathway backups, collateral
deletions and downstream vulnerabilities. The initial disease scope is
colorectal, gastric and lung cancer.

The repository combines:

- all 70 proteins in 15 curated human NISE activities;
- all 101 cross-cluster NISE pairs in both directions (202 hypotheses);
- Ensembl Compara paralog expansion;
- DepMap CRISPR, expression and WGS copy-number evidence;
- TCGA/GDC ASCAT3 homozygous-deletion events;
- STRING functional and DoRothEA regulatory networks;
- Human Protein Atlas localization;
- UniProt/PDB biochemical and structural traceability;
- all-target conditional-dependency discovery;
- Open Targets, ChEMBL, DGIdb, MyChem, Pharos/TCRD and CIViC evidence;
- optional PRISM, GDSC and CTRP drug-response analysis;
- a complete AlphaFold DB human NISE structural atlas;
- exact-numbered M-CSA, UniProt and PDBe functional-residue highlighting;
- a fully scripted publication package.

## Scientific boundary

The software prioritizes experimental hypotheses. It does not establish clinical
efficacy, patient benefit, safety, treatment suitability or cure.

AlphaFold DB models do not contain experimental substrates, cofactors or drug
poses. RSES-Onco highlights known residues on AlphaFold models but does not infer
that a ligand binds in a displayed pose. Only exact UniProt-numbered residues are
projected by default. Missing evidence is not converted to zero and reduces
explicit coverage.

## Canonical end-to-end protocol

The complete, command-by-command protocol from DepMap/GDC acquisition through all
analyses, all 40 figures, all tables, checksums and post-run validation is:

- [`docs/END_TO_END_ARTICLE_PROTOCOL.md`](docs/END_TO_END_ARTICLE_PROTOCOL.md)

STRING mapping, per-gene caching and recovery after a functional-evidence failure
are documented in:

- [`docs/STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md`](docs/STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md)
- [`docs/DOROTHEA_RECOVERY_WORKFLOW.md`](docs/DOROTHEA_RECOVERY_WORKFLOW.md)

The canonical post-run verification command is:

```bash
MPLBACKEND=Agg \
GDC_DIR=/absolute/path/to/data/raw/gdc \
PIPELINE_EXITCODE_FILE=logs/run_expanded_after_download_v0103.exitcode \
bash scripts/verify_complete_article_run.sh
```

Automated validation does not replace manual inspection of every figure at 100%
zoom.

## Installation

```bash
conda env create -f environment.yml
conda activate rses-onco
python -m pip install -e .
python -m pytest -q -p no:cacheprovider
```

For an existing environment, install structural rendering support:

```bash
conda activate rses-onco
conda install -c conda-forge pymol-open-source pillow
python -m pip install -e .
```

## Complete workflow after an existing GDC download

Do not restart an active GDC download. After all reviewed files finish and validate:

```bash
OLD="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-NISE"
NEW="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010"

cd "$NEW" || exit 1
conda activate rses-onco

git fetch origin
git checkout main
git pull --ff-only origin main

python -m pip install -e .
python -m pytest -q -p no:cacheprovider

export DEPMAP_DIR="$OLD/data/raw/depmap"
export GDC_DIR="$OLD/data/raw/gdc"

mkdir -p logs
set -o pipefail

MPLBACKEND=Agg \
PYTHONUNBUFFERED=1 \
DEPMAP_DIR="$DEPMAP_DIR" \
GDC_DIR="$GDC_DIR" \
bash scripts/run_expanded_pipeline.sh after-download \
  2>&1 | tee logs/run_expanded_after_download_v0103.log

status=${PIPESTATUS[0]}
echo "$status" > logs/run_expanded_after_download_v0103.exitcode
echo "Exit code: $status"
test "$status" -eq 0
```

This command performs all-NISE construction, composite-target normalization,
paralog expansion, all-target DepMap discovery, human network evidence, TCGA
integration, pharmacology, AlphaFold structure acquisition, functional-residue
annotation, PyMOL rendering, all tables, all figures, workbook creation, manifests,
checksums and tests.

After it finishes, execute:

```bash
MPLBACKEND=Agg \
GDC_DIR="$GDC_DIR" \
PIPELINE_EXITCODE_FILE=logs/run_expanded_after_download_v0103.exitcode \
bash scripts/verify_complete_article_run.sh \
  2>&1 | tee logs/verify_complete_article_run.log
```

## Resume after a functional-evidence interruption

Use this only when candidate construction, Ensembl expansion and all-target DepMap
discovery already completed, and the previous execution stopped at STRING,
DoRothEA, HPA or UniProt acquisition:

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

The STRING downloader maps symbols to exact STRING IDs, uses the version-specific
stable API, waits one second between network calls, writes per-gene caches and
reuses completed queries after interruption. It records unmapped identifiers as
missing coverage rather than score zero. Persistent STRING request failures remain
fatal after status and resume files have been written.

The resilient DoRothEA stage first reuses a valid local cache, then tries multiple
official OmniPath service addresses and both supported organism parameter spellings.
A persistent OmniPath outage is recorded in
`data/raw/human_functional_evidence/omnipath_dorothea_status.json`; regulatory-network
evidence remains missing and is not converted to zero. Set `DOROTHEA_STRICT=1` to
make the outage fatal, or set `DOROTHEA_FILE=/absolute/path/to/dorothea.tsv` to use a
validated local table.

## Structural atlas only

```bash
set -o pipefail
MPLBACKEND=Agg \
PYTHONUNBUFFERED=1 \
bash scripts/run_structural_pipeline.sh all \
  2>&1 | tee logs/run_structural_atlas.log

status=${PIPESTATUS[0]}
test "$status" -eq 0
```

Manual stages:

```bash
bash scripts/run_structural_pipeline.sh download
bash scripts/run_structural_pipeline.sh annotations
bash scripts/run_structural_pipeline.sh render
bash scripts/run_structural_pipeline.sh figures
```

## Publication workflow only

When the integrated ranking already exists:

```bash
bash scripts/run_publication_pipeline.sh all
```

Regenerate all figures from existing evidence and structural renders:

```bash
MPLBACKEND=Agg \
bash scripts/run_publication_pipeline.sh figures
```

The figure orchestrator is:

```bash
python -u scripts/make_all_article_figures.py \
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --discovery results/expanded_26Q1/discovery/all_target_dependency_screen.tsv \
  --pharmacology results/expanded_26Q1/pharmacology/pharmacology_ranked_hypotheses.tsv \
  --structure-manifest data/processed/structures/alphafold_structure_manifest.tsv \
  --render-manifest data/processed/structures/nise_structure_render_manifest.tsv \
  --structural-annotations data/processed/structures/nise_structural_residue_annotations.tsv \
  --structural-coverage data/processed/structures/nise_structural_annotation_coverage.tsv \
  --output-root article_outputs \
  --strict-layout
```

## Publication assets

The final package contains:

```text
8 main figures
32 supplementary figures
120 PNG/PDF/SVG files
4 main tables
18 supplementary tables
```

Structural organization:

- **Figure 8** — representative high-priority human NISE structures;
- **Figures S15-S29** — one legible AlphaFold structural atlas per NISE activity;
- **Figure S30** — AlphaFold confidence and annotation coverage;
- **Figure S31** — catalytic/binding/ligand-residue provenance;
- **Figure S32** — pairwise structural evidence.

Every composite figure is exported as:

```text
PNG: 600 dpi
PDF: vector-compatible layout
SVG: editable text and layout
```

Every registered figure receives:

```text
source-data TSV
layout-audit JSON
manifest entry
script/input/output provenance
```

Strict mode fails on panel overlap, tick collisions, clipped text, legends outside
the canvas, missing formats or missing source data. Manual inspection at 100% zoom
remains mandatory before submission.

## Final output structure

```text
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

## Documentation

- [`docs/END_TO_END_ARTICLE_PROTOCOL.md`](docs/END_TO_END_ARTICLE_PROTOCOL.md)
- [`docs/STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md`](docs/STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md)
- [`docs/DOROTHEA_RECOVERY_WORKFLOW.md`](docs/DOROTHEA_RECOVERY_WORKFLOW.md)
- [`docs/REAL_DATA_WORKFLOW.md`](docs/REAL_DATA_WORKFLOW.md)
- [`docs/EXPANDED_HUMAN_EVIDENCE_WORKFLOW.md`](docs/EXPANDED_HUMAN_EVIDENCE_WORKFLOW.md)
- [`docs/ALL_CLASS_AND_ALL_TARGET_DISCOVERY.md`](docs/ALL_CLASS_AND_ALL_TARGET_DISCOVERY.md)
- [`docs/PUBLICATION_PHARMACOLOGY_WORKFLOW.md`](docs/PUBLICATION_PHARMACOLOGY_WORKFLOW.md)
- [`docs/STRUCTURAL_ATLAS_WORKFLOW.md`](docs/STRUCTURAL_ATLAS_WORKFLOW.md)

## Key structural data resources

- AlphaFold Protein Structure Database;
- M-CSA Mechanism and Catalytic Site Atlas;
- UniProtKB reviewed residue features;
- PDBe/SIFTS/Arpeggio experimental ligand-binding evidence;
- optional user-curated exact UniProt residue mappings.

## License

MIT for code. Third-party data and structures retain their original licenses and
terms. The repository downloads them through source APIs rather than
redistributing the raw datasets.
