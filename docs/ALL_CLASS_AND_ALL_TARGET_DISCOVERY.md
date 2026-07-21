# All-class and all-target discovery workflow

This workflow complements the curated RSES-Onco benchmark by systematically
expanding the tested human functional-backup universe.

## What is exhaustive

The following layers are exhaustive relative to an explicit source and release:

- all 101 cross-cluster pairs in the bundled curated human NISE table;
- both loss-to-target directions for each pair, producing 202 NISE hypotheses;
- all 70 proteins and 15 NISE activity groups represented in that table;
- all human paralogs returned by the selected Ensembl Compara release for the
  analyzed seed genes;
- every CRISPR target measured by the selected DepMap release for each analyzable
  loss gene and cancer cohort.

No canonical database exhaustively defines every biological pathway backup,
collateral deletion or downstream dependency in the human proteome. Those classes
are therefore identified through a combination of source-bounded catalogues and
the all-target DepMap discovery screen. A statistical dependency is not assigned
a mechanism unless independent network, paralogy, reaction or literature evidence
supports that classification.

## Preferred continuation after an existing GDC download

Do not restart the active GDC transfer. After all files have finished and passed
size/MD5 validation, synchronize the v0.8 code and run:

```bash
bash scripts/run_expanded_pipeline.sh after-download
```

This single command performs:

1. DepMap schema and cohort validation;
2. complete directed NISE universe construction;
3. Ensembl Compara paralog expansion;
4. all-target DepMap conditional-dependency discovery;
5. universe reconstruction with supported discoveries;
6. STRING, DoRothEA, HPA and UniProt acquisition;
7. cancer-specific expression-context profiles;
8. cancer-specific CRISPR mutant-phenotype profiles;
9. conditional dependency and expression-compensation contrasts;
10. GDC file validation and ASCAT3 aggregation;
11. integrated TCGA plus DepMap scoring;
12. article tables, PDF/PNG/SVG figures and supplementary workbook;
13. tests, output manifest and SHA-256 checksums.

## Run only the expanded setup before TCGA

```bash
bash scripts/run_expanded_pipeline.sh setup
```

This stage does not require completed TCGA/GDC files. It creates the complete
NISE/paralog universe, performs all-target DepMap discovery and calculates the
human functional-microniche evidence available from DepMap and public networks.

## All-target discovery method

The default discovery screen uses every analyzable loss gene already present in
the NISE, paralog and benchmark universe. For each loss gene and cancer cohort:

- models are divided using linear WGS relative copy number `< 0.30`;
- at least three loss and three intact models are required;
- every measured DepMap CRISPR target is tested;
- the alternative hypothesis is a more negative target gene effect in loss
  models;
- one-sided Mann-Whitney tests are vectorized across all targets;
- Benjamini-Hochberg correction is performed across the complete target family
  for that loss-gene/cancer analysis;
- effect-size requirements are applied only after P and Q values are calculated;
- supported novel pairs are written as unvalidated discovery candidates.

Run manually:

```bash
python -u scripts/discover_conditional_dependencies.py \
  --gene-effect data/raw/depmap/CRISPRGeneEffect.csv \
  --copy-number data/raw/depmap/OmicsCNGeneWGS.csv \
  --models data/raw/depmap/Model.csv \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --loss-universe candidates \
  --loss-threshold 0.30 \
  --min-group-size 3 \
  --minimum-delta 0.15 \
  --maximum-median-loss-effect -0.25 \
  --fdr 0.10 \
  --output results/expanded_26Q1/discovery/all_target_dependency_screen.tsv \
  --candidate-output data/raw/discovery/depmap_discovered_candidate_pairs.tsv
```

Rebuild the candidate universe with discoveries:

```bash
python -u scripts/build_expanded_candidate_universe.py \
  --additional data/raw/ensembl/human_seed_paralogs.tsv \
  --additional data/raw/discovery/depmap_discovered_candidate_pairs.tsv \
  --output data/processed/expanded_candidate_universe.tsv \
  --members-output data/processed/expanded_class_member_inventory.tsv
```

## Optional all-copy-number-gene screen

A much broader and more computationally expensive analysis can use every copy-
number gene with enough loss models as the loss-gene universe:

```bash
bash scripts/run_expanded_pipeline.sh discover-all-cn
```

This screen may involve many thousands of loss-gene/target families and should be
run only with sufficient CPU time, memory and disk space. An optional safety cap
can be supplied directly to `discover_conditional_dependencies.py` with
`--max-loss-genes`.

## Human functional-microniche evidence

For every pair, the inner RSES integrates:

- expression/context: Spearman decorrelation and abundance separation within each
  cancer lineage;
- localization: HPA subcellular-location set divergence;
- biochemical/structural: curated NISE structural clusters and UniProt/PDB
  traceability;
- genetic/mutant phenotype: DepMap CRISPR profile decorrelation, score separation
  and dependency-set non-overlap;
- interaction network: STRING functional-neighborhood divergence and direct
  association score;
- regulatory network: DoRothEA regulator-set divergence.

Human CRISPR profiles are the scalable cancer-cell analogue of the condition-
resolved mutant phenotypes used in EcoNISE. They are not whole-organism knockout
phenotypes.

## Expanded outputs

```text
data/processed/expanded_candidate_universe.tsv
data/processed/expanded_class_member_inventory.tsv
data/processed/expanded_pair_functional_evidence.tsv

data/raw/ensembl/human_seed_paralogs.tsv
data/raw/discovery/depmap_discovered_candidate_pairs.tsv
data/raw/human_functional_evidence/

results/expanded_26Q1/discovery/
results/expanded_26Q1/depmap_only/
results/expanded_26Q1/full/
results/expanded_26Q1/full/article_tables/
figures/expanded_26Q1/
supplementary/RSES_Onco_Expanded_All_NISE_26Q1.xlsx
```

## Interpretation rules

- Missing evidence is omitted and lowers coverage; it is never assigned zero.
- NISE status, paralogy, STRING connectivity or shared regulation does not prove
  synthetic lethality.
- A novel all-target dependency remains an empirical discovery candidate until
  mechanism and reproducibility are independently established.
- The all-target Q value controls the tested target family for one loss gene and
  one cancer. It must be reported with the family size, effect size, group sizes
  and analysis release.
- Final candidates require biomarker-matched isogenic perturbation, rescue,
  orthogonal pharmacology and mechanistic validation.
