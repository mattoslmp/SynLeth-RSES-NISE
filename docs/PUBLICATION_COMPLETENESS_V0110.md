# Publication completeness contract — v0.11.0

Version 0.11.0 extends the executable publication contract without changing the RSES-Onco v0.10.9 score semantics.

## Registered assets

- 8 main figures;
- 69 supplementary figures;
- 231 figure files across PNG, PDF and SVG;
- 4 main tables;
- 44 supplementary tables;
- exact source TSV and layout audit for every figure;
- DOCX and rendered PDF manuscript/supplement;
- page-render directories for manual inspection;
- explicit page separation for Figures S68 and S69.

## New source-backed supplementary series

Figures S39–S69 cover score decomposition, weight and missing-data sensitivity, expression and coexpression conditioning, WGCNA modules and correlation fallback, regulatory subcomponents, STRING evidence channels, conditional CRISPR phenotypes, localization, biochemical/structural coverage, TCGA event context, pharmacology, NISE/paralog comparisons, controls, eligibility/coverage and integrated network context.

No figure is generated solely to increase the count. When a registered source has no eligible observations, the exact source table records that state and the figure displays an explicit evidence-unavailable boundary rather than fabricated values.

## Canonical commands

```bash
MPLBACKEND=Agg bash scripts/run_publication_pipeline.sh assets-only
bash scripts/run_publication_pipeline.sh documents
bash scripts/verify_complete_article_run.sh
```

Automated validation does not replace manual page inspection.
