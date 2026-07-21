# Data dictionary

## Human NISE catalogue

`data/curated/human_nise_bonafide_2017.tsv` contains the 70 protein sequences in 15 enzymatic activities classified as the bona fide human intragenomic analogous-enzyme set by Piergiorge et al. (2017). `structural_cluster` is the cluster label reported in Table 1, not a new phylogeny.

## Synthetic-lethality reference pairs

`data/curated/synthetic_lethality_reference_pairs.tsv` contains literature-anchored positive controls and RSES-Onco discovery hypotheses. Numeric fields are transparent evidence priors on [0,1], not measured TCGA or DepMap values. The empirical pipeline replaces the event, dependency and selectivity components when compatible matrices are supplied.

## Expected DepMap inputs

- `CRISPRGeneEffect.csv`: ModelID rows, gene columns; more negative values indicate stronger dependency.
- `OmicsCNGene.csv`: ModelID rows, gene-level copy-number values.
- `Model.csv`: model metadata, including lineage.
- `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`: optional expression matrix.

## Expected TCGA/GISTIC input

A gene-by-sample matrix with a gene-symbol first column and discrete values -2, -1, 0, 1, 2. A value of -2 is treated as homozygous deletion. Do not substitute continuous segment means without changing the threshold and recording that decision.
