# Human NISE AlphaFold structural atlas workflow

## Scope

RSES-Onco v0.10 generates a scripted structural atlas for every protein in the
curated human NISE catalogue:

- 70 human proteins;
- 15 activity groups;
- AlphaFold DB PDB, mmCIF, metadata and PAE downloads;
- exact UniProt-numbered functional residues from UniProtKB reviewed features;
- exact-accession M-CSA catalytic-residue queries;
- PDBe/Arpeggio ligand-contact residues only when UniProt numbering is explicitly
  available;
- optional user-curated exact residue mappings;
- whole-structure and enlarged functional-site PyMOL renders;
- Figure 8 and Figures S15-S32 in PNG, PDF and SVG;
- three structural supplementary tables and full SHA-256 provenance.

## Scientific boundary

AlphaFold DB models contain predicted protein coordinates and confidence values,
not substrates, cofactors or drug poses. The pipeline therefore highlights known
residues on the AlphaFold coordinate model but does not claim that a ligand or drug
was predicted to bind in a specific pose.

Residue projection rules are conservative:

1. UniProtKB reviewed features are used in canonical UniProt numbering.
2. M-CSA is queried by exact UniProt accession.
3. PDBe binding residues are used only when the API explicitly returns UniProt
   residue numbers.
4. PDB-only numbering and homology-transfer residues are not projected without an
   explicit mapping.
5. User-curated residues must state accession, residue number, source and mapping
   status.

M-CSA homologue files are not used by default because the database explicitly
warns that homologous residues may not perform the same catalytic function.

## Installation

For an existing environment:

```bash
conda activate rses-onco

conda install -c conda-forge \
  pymol-open-source \
  pillow

python -m pip install -e .
```

Verify:

```bash
pymol -cq <<'EOF'
print('PyMOL headless mode available')
quit
EOF

bash -n scripts/run_structural_pipeline.sh
python -m pytest -q -p no:cacheprovider tests/test_structural.py
```

## One-command structural workflow

After the expanded TCGA plus DepMap ranking exists:

```bash
set -o pipefail

MPLBACKEND=Agg \
PYTHONUNBUFFERED=1 \
bash scripts/run_structural_pipeline.sh all \
  2>&1 | tee logs/run_structural_atlas.log

status=${PIPESTATUS[0]}
echo "Structural pipeline exit code: $status"
test "$status" -eq 0
```

The complete publication orchestrator also runs this stage automatically:

```bash
bash scripts/run_publication_pipeline.sh all
```

The complete post-GDC workflow is:

```bash
bash scripts/run_expanded_pipeline.sh after-download
```

## Manual execution order

### 1. Download AlphaFold DB models

```bash
python -u scripts/download_alphafold_nise_structures.py \
  --proteins data/curated/human_nise_bonafide_2017.tsv \
  --output-dir data/raw/structures/alphafold \
  --manifest data/processed/structures/alphafold_structure_manifest.tsv
```

The downloader:

- queries the AlphaFold DB prediction API per UniProt accession;
- preserves fragments for long proteins;
- downloads PDB, mmCIF and PAE JSON when available;
- writes via `.part` files;
- validates non-empty files;
- records model version, fragment range, mean pLDDT and SHA-256.

### 2. Collect known functional residues

```bash
python -u scripts/collect_nise_structural_annotations.py \
  --proteins data/curated/human_nise_bonafide_2017.tsv \
  --output data/processed/structures/nise_structural_residue_annotations.tsv \
  --coverage-output data/processed/structures/nise_structural_annotation_coverage.tsv \
  --cache-dir data/raw/structures/annotation_cache
```

The output distinguishes:

```text
mcsa_catalytic
uniprot_active_site
uniprot_binding
uniprot_metal
uniprot_site
pdbe_ligand
drug_binding
curated_user
```

### 3. Add optional exact drug-contact residues

A user-curated TSV may be supplied:

```text
gene_symbol
uniprot_accession
residue_number
residue_name
annotation_type
source
description
ligand_code
ligand_name
pdb_id
evidence_level
mapping_status
```

Run:

```bash
python -u scripts/collect_nise_structural_annotations.py \
  --curated-residues data/curated/nise_exact_drug_contact_residues.tsv \
  --output data/processed/structures/nise_structural_residue_annotations.tsv \
  --coverage-output data/processed/structures/nise_structural_annotation_coverage.tsv
```

Only exact UniProt-numbered residues should be entered.

### 4. Render all proteins with PyMOL

```bash
python -u scripts/render_nise_structures.py \
  --structure-manifest data/processed/structures/alphafold_structure_manifest.tsv \
  --annotations data/processed/structures/nise_structural_residue_annotations.tsv \
  --output-dir article_outputs/structure_atlas/individual \
  --render-manifest data/processed/structures/nise_structure_render_manifest.tsv \
  --pymol pymol
```

For every protein, the script creates:

```text
<gene>_<accession>_whole.png
<gene>_<accession>_site.png
<gene>_<accession>_whole.pml
<gene>_<accession>_site.pml
PyMOL logs
```

Raw molecular renders are 2800 × 2200 pixels at 600 dpi. Known residues are
shown as sticks and CA spheres. Text labels are excluded from the molecular image
and placed in dedicated figure regions later, preventing text-structure overlap.

### 5. Generate the article structural figures

```bash
MPLBACKEND=Agg \
python -u scripts/make_nise_structure_figures.py \
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv \
  --proteins data/curated/human_nise_bonafide_2017.tsv \
  --structure-manifest data/processed/structures/alphafold_structure_manifest.tsv \
  --render-manifest data/processed/structures/nise_structure_render_manifest.tsv \
  --annotations data/processed/structures/nise_structural_residue_annotations.tsv \
  --coverage data/processed/structures/nise_structural_annotation_coverage.tsv \
  --output-root article_outputs \
  --strict-layout
```

## Figure organization

### Main figure

- **Figure 8** — representative top-ranked human NISE pairs, showing whole
  AlphaFold models and enlarged known functional-residue views for both members.

### Activity-specific supplementary structural atlas

- **Figures S15-S29** — one high-resolution figure per NISE activity group,
  covering all 70 proteins without placing all models in one unreadably dense
  panel.

### Structural QC and evidence figures

- **Figure S30** — mean AlphaFold pLDDT and annotation coverage for all proteins.
- **Figure S31** — residue-evidence provenance by source and annotation class.
- **Figure S32** — pairwise structural evidence across directional NISE
  vulnerabilities.

## Color policy

```text
AlphaFold very high confidence     marine
AlphaFold confident                cyan
AlphaFold low confidence           orange
AlphaFold very low confidence      salmon
M-CSA / active-site residues       red
UniProt binding residues           orange
Metal-binding residues             violet
Other important sites              yellow
PDBe experimental ligand contacts  magenta
Curated exact drug contacts        cyan
```

## Publication quality

All composite structural figures use:

- dynamic figure height based on the number of proteins;
- no text inside molecular render panels;
- dedicated title and residue-description areas;
- large fonts for 100% zoom inspection;
- `constrained_layout`;
- 600-dpi PNG;
- PDF and SVG exports;
- strict layout audits;
- source-data TSVs;
- figure manifests and legends.

The final automated package contains:

```text
8 main figures
32 supplementary figures
120 PNG/PDF/SVG files
4 main tables
18 supplementary tables
```

## Structural outputs

```text
data/raw/structures/
├── alphafold/
│   ├── pdb/
│   ├── cif/
│   ├── pae/
│   └── metadata/
└── annotation_cache/

data/processed/structures/
├── alphafold_structure_manifest.tsv
├── nise_structural_residue_annotations.tsv
├── nise_structural_annotation_coverage.tsv
├── nise_structural_annotation_source_status.tsv
└── nise_structure_render_manifest.tsv

article_outputs/
├── structure_atlas/individual/
├── figures/main/Figure_8_human_nise_structures.*
├── figures/supplementary/Figure_S15_... through Figure_S32_...
├── source_data/figures/structures/
└── tables/supplementary/Table_S16...Table_S18...
```

## Required manual review

Automated layout checks detect common clipping and overlap problems, but all
figures must still be opened at 100% zoom before submission. Review every PDF and
PNG for:

- residue visibility;
- adequate contrast;
- appropriate molecular orientation;
- sufficiently enlarged active-site views;
- no panel, legend, title or colorbar overlap;
- no misleading interpretation of low-confidence regions;
- exact correspondence between residue annotations, source tables and captions.
