# Expanded human functional-microniche workflow (RSES-Onco v0.8)

## Purpose

The expanded workflow no longer limits empirical analysis to the 25 literature
benchmarks. It evaluates:

1. every one of the 101 curated cross-cluster human NISE pairs;
2. both directional hypotheses for each NISE pair (202 directions);
3. every protein in the 70-protein/15-activity curated human NISE catalogue;
4. the curated synthetic-lethality benchmark classes;
5. all Ensembl Compara paralogs returned for the analyzed seed genes;
6. optional additional source-bounded class catalogues supplied as TSV/CSV.

A directional pair `A loss -> B dependency` is different from `B loss -> A
dependency`; both are retained and tested independently.

## Scientific scope and the meaning of “all”

“All NISEs” means exhaustive coverage of the bundled curated human NISE table.
“All paralogs” means all human paralogs returned by the selected Ensembl release
for the analyzed seed genes. There is no canonical database that exhaustively
defines every possible pathway backup, downstream dependency or collateral
metabolic vulnerability in the human proteome. Those classes are therefore
source-bounded and must retain their provenance. Additional catalogues can be
added with repeated `--additional` arguments to
`scripts/build_expanded_candidate_universe.py`.

The software must never label an untested pair as synthetic lethal merely because
it is a NISE, paralog, pathway neighbor or interaction partner.

## Human equivalents of the EcoNISE evidence domains

| EcoNISE evidence concept | Human RSES-Onco implementation |
|---|---|
| Expression/context | Cancer-specific DepMap expression correlation and abundance separation; conditional target induction in loss versus intact models |
| Localization | Human Protein Atlas subcellular-localization annotations |
| Biochemical/structural | Curated NISE EC activity and structural-cluster difference; UniProt catalytic, cofactor, domain, localization and PDB traceability |
| Genetic/mutant phenotype | DepMap Chronos CRISPR gene-effect profiles, essential-set overlap and cancer-specific loss-versus-intact contrasts |
| Functional interaction network | STRING taxon 9606 functional associations and pairwise neighborhood divergence |
| Regulatory network | OmniPath/DoRothEA TF-target interactions and pairwise regulator-set divergence |
| Tumor event | TCGA/GDC ASCAT3 zero-copy homozygous-deletion frequency |

Human cell-line CRISPR profiles are the closest scalable analogue of the
condition-resolved mutant phenotypes used in the E. coli analysis. They do not
represent whole-organism knockout phenotypes and must be described accordingly.

## Two-level score

### Functional-microniche RSES

The inner score integrates six domains:

```text
expression context        0.20
localization              0.15
biochemical/structural    0.15
genetic/phenotype         0.20
STRING interaction        0.15
regulatory network        0.15
```

### Expanded RSES-Onco

The cancer-priority score integrates:

```text
tumor event               0.16
conditional dependency    0.22
loss selectivity          0.14
expression compensation   0.08
functional relation       0.06
functional-microniche     0.16
validation/tractability   0.18
```

For both levels, unavailable evidence is omitted from the observed-domain
denominator and reduces explicit coverage. Missing values are never converted to
zero.

## Installation or synchronization

```bash
conda activate rses-onco
python -m pip install -e .
python -m pytest -q -p no:cacheprovider
```

## Build all directed NISE candidates and benchmark classes

```bash
python -u scripts/build_expanded_candidate_universe.py \
  --output data/processed/expanded_candidate_universe.tsv \
  --members-output data/processed/expanded_class_member_inventory.tsv
```

Expected minimum before paralog expansion:

```text
202 directed NISE hypotheses
70 unique NISE proteins
25 literature benchmark rows before overlap handling
```

The exact combined direction count can be smaller than the arithmetic sum when a
benchmark duplicates a direction already present in the systematic NISE table.

## Expand homologous paralogs with Ensembl Compara

```bash
python -u scripts/download_ensembl_paralogs.py \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --output data/raw/ensembl/human_seed_paralogs.tsv

python -u scripts/build_expanded_candidate_universe.py \
  --additional data/raw/ensembl/human_seed_paralogs.tsv \
  --output data/processed/expanded_candidate_universe.tsv \
  --members-output data/processed/expanded_class_member_inventory.tsv
```

The Ensembl output is a candidate class, not proof of functional redundancy or
synthetic lethality.

## Acquire STRING, regulation, localization and biochemical traceability

```bash
PYTHONUNBUFFERED=1 \
python -u scripts/download_human_functional_evidence.py \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --raw-dir data/raw/human_functional_evidence \
  --output data/processed/expanded_pair_functional_evidence.tsv \
  2>&1 | tee logs/empirical_26Q1/03d_download_human_functional_evidence.log
```

Raw source tables are retained under:

```text
data/raw/human_functional_evidence/
├── string_interaction_partners.tsv
├── omnipath_dorothea.tsv
├── hpa_subcellular_location.tsv
└── uniprot_reviewed_annotations.tsv
```

The pair table includes STRING neighborhood Jaccard/divergence, direct STRING
score, shared regulators, regulator-set Jaccard/divergence, HPA localization
divergence and structural-cluster evidence.

To rebuild pair metrics from already downloaded raw evidence without new network
requests:

```bash
python -u scripts/download_human_functional_evidence.py \
  --skip-download \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --raw-dir data/raw/human_functional_evidence \
  --output data/processed/expanded_pair_functional_evidence.tsv
```

## Run expanded DepMap analysis

```bash
python -u scripts/run_expanded_rses_onco.py \
  --gene-effect data/raw/depmap/CRISPRGeneEffect.csv \
  --copy-number data/raw/depmap/OmicsCNGeneWGS.csv \
  --models data/raw/depmap/Model.csv \
  --expression data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --functional-evidence data/processed/expanded_pair_functional_evidence.tsv \
  --loss-threshold 0.30 \
  --min-group-size 3 \
  --output results/expanded_26Q1/depmap_only/expanded_rses_onco.tsv
```

Outputs include:

```text
expanded_rses_onco.tsv
expanded_dependency_contrasts.tsv
expanded_expression_compensation.tsv
expanded_expression_context_profiles.tsv
expanded_crispr_phenotype_profiles.tsv
```

Every candidate remains in the ranking even when its loss cohort cannot be
defined. In that case empirical components remain missing and coverage is lower.

## Run expanded TCGA plus DepMap analysis

After GDC download, validation and aggregation:

```bash
python -u scripts/run_expanded_rses_onco.py \
  --gene-effect data/raw/depmap/CRISPRGeneEffect.csv \
  --copy-number data/raw/depmap/OmicsCNGeneWGS.csv \
  --models data/raw/depmap/Model.csv \
  --expression data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --functional-evidence data/processed/expanded_pair_functional_evidence.tsv \
  --loss-threshold 0.30 \
  --min-group-size 3 \
  --tcga colon=data/processed/TCGA_COLON_homdel_discrete.tsv \
  --tcga stomach=data/processed/TCGA_STOMACH_homdel_discrete.tsv \
  --tcga lung=data/processed/TCGA_LUNG_homdel_discrete.tsv \
  --output results/expanded_26Q1/full/expanded_rses_onco.tsv
```

## One-command expanded continuation after the current GDC download

```bash
bash scripts/run_real_data_pipeline.sh expanded-after-download
```

This runs:

1. all-NISE candidate construction;
2. Ensembl paralog expansion;
3. STRING, DoRothEA, HPA and UniProt acquisition;
4. expanded DepMap expression and CRISPR phenotype profiles;
5. GDC validation and aggregation;
6. expanded TCGA plus DepMap scoring;
7. manuscript-ready tables;
8. figures in PDF, PNG and SVG;
9. expanded supplementary workbook;
10. tests, manifests and SHA-256 checksums.

## Optional additional classes

A standardized additional table may contain:

```text
pair_id
lost_feature
lost_gene
target_gene
source_class
relation_type
mechanism
colon
stomach
lung
relation_confidence
genetic_screen
isogenic_validation
in_vivo
clinical_tractability
lineage_relevance
evidence_stage
primary_doi
supporting_doi
status
```

Add one or more catalogues:

```bash
python -u scripts/build_expanded_candidate_universe.py \
  --additional data/curated/pathway_backups.tsv \
  --additional data/curated/collateral_deletions.tsv \
  --additional data/curated/downstream_dependencies.tsv \
  --output data/processed/expanded_candidate_universe.tsv
```

Each catalogue must state its source, release, inclusion rule and evidence stage.

## Interpretation boundary

Network divergence, common regulation, shared localization, paralogy or reaction
equivalence can prioritize a hypothesis but cannot establish synthetic lethality.
The strongest candidates require convergent tumor-event, dependency, selectivity,
functional-microniche and experimental-validation evidence. Final confirmation
requires biomarker-matched isogenic perturbation, rescue, pharmacology and
mechanistic assays.
