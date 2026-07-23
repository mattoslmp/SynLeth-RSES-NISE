# Genomic Circos implementation contract — RSES-Onco v0.11.2

This repository contract is enforced by the canonical publication pipeline and its
test suite.

## Figure contract

- Supplementary Figure S70 is generated in PNG, PDF and SVG.
- A strict layout audit must return `status=pass` and no warnings.
- GRCh38 chromosomes 1-22, X, Y and MT are supported.
- Every coordinate-complete simple NISE/paralog gene is represented as a genomic
  tick.
- Every simple NISE/paralog candidate pair is represented by exactly one chord,
  including pairs with missing score values.
- NISE chords are red (`#C62828`).
- homologous-paralog chords are black (`#111111`).
- pairs with missing scores remain in the figure with `link_status=score_missing`.

## Ring contract

Panel A contains the coverage-adjusted score, evidence coverage and all seven
RSES-Onco domains. Panel B contains all six functional-microniche domains,
pairwise expression, WGCNA, DoRothEA, TF-expression consistency, JASPAR/FIMO
promoter motifs, promoter methylation and nested coverage. Missing/non-eligible
values remain `NA` and are rendered as hollow markers rather than zero.

## Supplementary-data contract

- S45: genomic coordinates and gene class;
- S46: all pair links and rendering fields;
- S47: all gene-by-ring values and missingness status;
- S48: ring definitions and source columns;
- S49: expression summary for every gene × cancer context;
- S50: every observed model-level expression value plus explicit `NA` sentinel
  rows for unavailable gene/context combinations;
- S51: complete generated catalogue of every Python, Bash and R pipeline source;
- S52: source paths, sizes and SHA-256 provenance.

The exact combined source table used to render S70 is copied byte-for-byte to:

```text
article_outputs/tables/figure_data/supplementary/Figure_S70_source_data.tsv
```

## Validation contract

The following commands must pass:

```bash
python -u scripts/validate_genomic_circos_integrity.py \
  --article-root article_outputs

python -u scripts/validate_publication_outputs.py \
  --article-root article_outputs

python -m pytest -q -p no:cacheprovider
```

The final publication package contains 8 main figures, 70 supplementary figures,
234 PNG/PDF/SVG files, 4 main tables and 52 supplementary tables.
