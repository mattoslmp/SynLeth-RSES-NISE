# RSES-Onco v0.11.0: data acquisition, provenance and reproduction guide

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

## 1. Purpose

This document describes the origin, acquisition or recovery route, expected local path, validation procedure and downstream use of the principal data resources used by RSES-Onco. It complements `END_TO_END_ARTICLE_PROTOCOL.md`, which contains the complete execution order.

The guide follows four rules:

1. raw public data are acquired from their official provider or from a traceable validated local copy;
2. release, access date, accession or source-status metadata are recorded whenever available;
3. missing or technically unavailable evidence remains missing and is not replaced by zero;
4. the same original evidence unit is not counted independently merely because it appears in more than one aggregator.

## 2. Recommended directory layout

```text
data/
├── curated/
├── raw/
│   ├── depmap/
│   ├── gdc/
│   ├── ensembl/
│   ├── human_functional_evidence/
│   ├── pharmacology/
│   ├── regulatory/
│   └── structures/
└── processed/
    ├── regulatory/
    ├── pharmacology/
    └── structures/
```

Large DepMap or GDC releases may remain outside the Git working tree through `DEPMAP_DIR` and `GDC_DIR`, but the exact paths and release identifiers must be recorded.

## 3. Human NISE catalogue

### Scientific role

The curated NISE catalogue defines human proteins that catalyze equivalent or closely related functions despite belonging to different structural or evolutionary groups. It is used to construct directed candidate hypotheses.

### Scripts and paths

```text
scripts/download_human_nise.py
data/curated/
data/processed/expanded_candidate_universe.tsv
data/processed/expanded_class_member_inventory.tsv
```

### Validation requirements

- gene symbols and UniProt identifiers must remain traceable;
- activity/EC grouping must be preserved;
- both directions of a pair must receive distinct directional identifiers;
- absence of validation must not be rewritten as evidence against the hypothesis.

## 4. Ensembl Compara paralogs

### Scientific role

Ensembl paralogs expand the candidate universe with homologous backup relationships that can be compared with NISE-based hypotheses.

### Script

```bash
python -u scripts/download_ensembl_paralogs.py
```

### Provenance

Record Ensembl release, species, query date, stable gene identifiers, relationship type and any unmapped symbol. The downstream candidate table must preserve the Ensembl source fields.

## 5. DepMap dependency, copy number, models and expression

### Required files

```text
data/raw/depmap/CRISPRGeneEffect.csv
data/raw/depmap/OmicsCNGeneWGS.csv
data/raw/depmap/Model.csv
data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
```

### Acquisition

The files are obtained from the official DepMap data portal for a consistent release. `scripts/download_depmap.py` validates and fingerprints local files; it does not bypass official access requirements.

```bash
export DEPMAP_DIR="$PWD/data/raw/depmap"

python -u scripts/download_depmap.py \
  --input-dir "$DEPMAP_DIR" \
  --write-checksums

python -u scripts/validate_real_inputs.py \
  --gene-effect "$DEPMAP_DIR/CRISPRGeneEffect.csv" \
  --copy-number "$DEPMAP_DIR/OmicsCNGeneWGS.csv" \
  --models "$DEPMAP_DIR/Model.csv" \
  --expression "$DEPMAP_DIR/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
```

### Downstream use

- general and context-specific dependency;
- event-stratified selectivity;
- gene expression and compensation;
- pairwise expression and phenotype divergence;
- cancer-specific WGCNA networks;
- TF-expression consistency.

Release-specific model counts may change. Record the observed values rather than forcing older counts.

## 6. TCGA/GDC tumor-event data

### Projects

```text
TCGA-COAD
TCGA-READ
TCGA-STAD
TCGA-LUAD
TCGA-LUSC
```

### Manifest and download

```bash
export GDC_DIR="$PWD/data/raw/gdc"
mkdir -p "$GDC_DIR" logs/expanded_26Q1

python -u scripts/download_gdc.py \
  --manifest-only \
  --workflow ASCAT3 \
  --output-dir "$GDC_DIR"

python -u scripts/download_gdc.py \
  --use-existing-manifest \
  --manifest "$GDC_DIR/gdc_gene_level_copy_number_manifest.json" \
  --workflow ASCAT3 \
  --output-dir "$GDC_DIR" \
  --retries 3

python -u scripts/download_gdc.py \
  --validate-only \
  --manifest "$GDC_DIR/gdc_gene_level_copy_number_manifest.json" \
  --output-dir "$GDC_DIR"
```

### Aggregation

```bash
python -u scripts/aggregate_gdc_gene_cna.py
python -u scripts/validate_gdc_matrices.py
```

Expected processed outputs include:

```text
data/processed/TCGA_COLON_homdel_discrete.tsv
data/processed/TCGA_STOMACH_homdel_discrete.tsv
data/processed/TCGA_LUNG_homdel_discrete.tsv
```

No `.part` file may remain before aggregation. Frequencies used in figures must remain within the biologically valid interval and graphical jitter, when used, must never overwrite the source value.

## 7. STRING functional-network evidence

### Scripts

```text
scripts/download_human_functional_evidence.py
scripts/download_human_functional_evidence_resilient.py
scripts/export_raw_functional_network_evidence.py
scripts/export_supporting_evidence_tables.py
```

### Interpretation

The STRING combined score is not equivalent to direct experimental binding. Preserve and export the available channels, including experimental, database, coexpression, neighborhood, cooccurrence, fusion and text-mining evidence. Per-gene caches and mapping status must be retained for recovery after interruption.

See `STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md` for the recovery procedure.

## 8. DoRothEA and OmniPath regulatory evidence

### Scripts

```text
scripts/download_human_functional_evidence_resilient.py
scripts/download_dorothea_official_static.py
scripts/export_wgcna_regulatory_supporting_evidence.py
```

### Expected source table

```text
data/raw/human_functional_evidence/omnipath_dorothea.tsv
```

A validated official static DoRothEA table may be used when OmniPath is unavailable. The source-status JSON and acquisition route must be preserved. Regulatory association must not be presented as direct promoter binding unless the underlying source provides such evidence.

## 9. Human Protein Atlas localization

### Script

```bash
python -u scripts/download_hpa_subcellular_current.py
```

### Downstream use

Localization evidence records each protein's compartments, source reliability, shared or distinct compartments and compatibility with the proposed specialization or compensatory relationship. The existence of any annotation does not automatically imply maximum localization support.

## 10. UniProt, biochemical and residue evidence

UniProt identifiers, molecular functions, reactions, cofactors, catalytic residues and reviewed feature annotations are retained as traceable supporting evidence. Absence of a reviewed feature is missing evidence, not evidence that a feature does not exist.

Relevant scripts include:

```text
scripts/download_human_functional_evidence.py
scripts/collect_nise_structural_annotations.py
scripts/export_supporting_evidence_tables.py
```

## 11. Ensembl promoters, JASPAR and FIMO

### Promoter acquisition

```bash
python -u scripts/download_ensembl_promoters.py
```

Expected outputs:

```text
data/raw/regulatory/ensembl_promoters.tsv
data/raw/regulatory/ensembl_promoters.fa
```

The promoter table must retain canonical transcript, TSS, genomic coordinates, strand, window definition, source release and retrieval status.

### Motif acquisition and scanning

```bash
python -u scripts/download_jaspar_core_vertebrates.py
python -u scripts/scan_promoter_motifs.py
```

Expected resource and summaries:

```text
data/raw/regulatory/JASPAR2026_CORE_vertebrates_non-redundant.meme
data/processed/regulatory/jaspar_promoter_motif_hits.tsv
data/processed/regulatory/jaspar_promoter_tf_summary.tsv
```

A motif occurrence is sequence-based cis-regulatory support. It is not ChIP evidence, occupancy, causal regulation or direct binding proof.

## 12. WGCNA expression-network evidence

### Scripts

```text
scripts/build_wgcna_regulatory_layer.py
scripts/run_wgcna_expression_network.R
scripts/aggregate_wgcna_regulatory_layer.py
```

### Policy

```text
primary correlation: bicor
maxPOutliers: 0.10
pearsonFallback: individual
network type: signed
TOM type: signed
```

Pearson is used only for individual zero/non-finite-MAD entities and all fallbacks are exported. Required all-cancer outputs include:

```text
data/processed/regulatory/wgcna/wgcna_pair_metrics_all_cancers.tsv
data/processed/regulatory/wgcna/wgcna_input_preparation.tsv
data/processed/regulatory/wgcna/wgcna_correlation_fallback_all_cancers.tsv
data/processed/regulatory/wgcna/wgcna_run_diagnostics_all_cancers.tsv
```

## 13. Pharmacology and tractability

### Sources and scripts

Open Targets, ChEMBL, DGIdb, MyChem, Pharos/TCRD and CIViC are acquired or queried through:

```text
scripts/acquire_pharmacology_evidence.py
scripts/prioritize_pharmacology.py
```

Optional PRISM, GDSC and CTRP releases are standardized and analyzed through:

```text
scripts/standardize_drug_sensitivity.py
scripts/analyze_drug_response_selectivity.py
```

Tractability is not clinical efficacy. Compound mechanism, assay type, model, development status, selectivity and source must remain separate columns.

## 14. Structural resources

### Scripts

```text
scripts/download_alphafold_nise_structures.py
scripts/collect_nise_structural_annotations.py
scripts/render_nise_structures.py
scripts/make_nise_structure_figures.py
scripts/run_structural_pipeline.sh
```

### Sources

The workflow uses AlphaFold DB models and traceable annotations from resources such as UniProt, M-CSA and PDBe when available. AlphaFold confidence is not experimental validation. Missing experimental structure is not evidence of structural absence.

## 15. Candidate construction and all-target discovery

```bash
python -u scripts/build_expanded_candidate_universe.py
python -u scripts/discover_conditional_dependencies.py
```

Expected outputs include the directed candidate universe, class inventory and all-target dependency screen. Directionality, composite events, source class, cancer context and mapping status must be preserved.

## 16. Scoring and source overlap

The pipeline distinguishes evidence used for:

- direct score calculation;
- prioritization;
- independent validation;
- interpretation only.

Representations sharing a publication, original dataset or evidence identifier are grouped into one evidence unit. Multiple aggregators may corroborate an interpretation but may not receive multiple full independent weights for the same underlying observation.

## 17. Publication asset generation

```bash
MPLBACKEND=Agg \
STRICT_LAYOUT=1 \
bash scripts/run_publication_pipeline.sh assets-only
```

Every registered figure must receive:

- exact source-data TSV;
- generating script and input paths;
- reproduction command;
- PNG, PDF and SVG;
- layout-audit JSON;
- manifest entry.

The current contract is 77 figures, 231 image exports and 48 registered tables.

## 18. Document generation and packaging

```bash
bash scripts/run_publication_pipeline.sh documents
python -u scripts/validate_publication_documents.py \
  --article-root article_outputs \
  --document-dir article_outputs/documents \
  --require-page-renders
```

Manual inspection remains mandatory before packaging. Complete execution and packaging commands are in `END_TO_END_ARTICLE_PROTOCOL.md`.

## 19. Minimum provenance fields

Each acquired or generated dataset should record, when applicable:

```text
resource_name
source_provider
release_or_version
access_date
original_accession
license_or_terms
acquisition_script
raw_path
processed_path
mapping_status
rows_or_records
checksum
score_role
figures_using_source
tables_using_source
known_limitations
```

## 20. Source-specific companion documents

- `END_TO_END_ARTICLE_PROTOCOL.md`
- `STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md`
- `DOROTHEA_RECOVERY_WORKFLOW.md`
- `HPA_SUBCELLULAR_RECOVERY.md`
- `PUBLICATION_PHARMACOLOGY_WORKFLOW.md`
- `STRUCTURAL_ATLAS_WORKFLOW.md`
- `PUBLICATION_EVIDENCE_AUDIT_AND_REPRODUCTION.md`
