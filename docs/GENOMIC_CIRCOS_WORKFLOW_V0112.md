# RSES-Onco v0.11.2 genomic Circos workflow

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

## Purpose

Supplementary Figure S70 is a source-backed genomic Circos representation of every
coordinate-complete simple-gene hypothesis classified as a non-homologous
isofunctional enzyme (NISE), a homologous paralog, or both. It combines genomic
position, every pair relationship, all top-level RSES-Onco domains, all
functional-microniche domains, WGCNA internal terms, regulatory terms, promoter
methylation, validation terms and nested evidence coverage.

## Genomic coordinates

Canonical coordinates are read from:

```text
data/raw/regulatory/ensembl_promoters.tsv
```

The canonical-transcript TSS is used as the genomic position. Chromosomes 1-22, X,
Y and MT are supported on GRCh38. A gene without a supported coordinate is not
assigned an invented location; the input stage terminates and reports the gene.

## Included hypotheses and links

The universe is `data/processed/expanded_candidate_universe.tsv`. Every simple pair
identified by `source_class`, `relation_type` or `ensembl_homology_type` as a NISE
or paralog/homolog is retained exactly once.

- NISE chords: red (`#C62828`);
- homologous-paralog chords: black (`#111111`);
- genes represented in both classes: purple genomic ticks;
- chord width/transparency: maximum cancer-specific
  `coverage_adjusted_rses`;
- pairs without a score row remain as low-opacity chords with
  `link_status=score_missing`;
- composite features that cannot be mapped to two genes are not assigned fictitious
  positions.

## Thirty-five ring tracks

### Panel A — observed score, top-level domains and validation terms

1. observed-domain RSES-Onco;
2. coverage-adjusted RSES-Onco;
3. evidence coverage;
4. tumor event;
5. conditional dependency;
6. selectivity;
7. expression compensation;
8. functional relation;
9. functional microniche;
10. validation and tractability;
11. genetic-screen validation;
12. isogenic validation;
13. in vivo validation;
14. clinical tractability.

### Panel B — functional microniche and all score-internal layers

1. expression context;
2. localization;
3. biochemical/structural evidence;
4. genetic phenotype;
5. interaction network;
6. regulatory network;
7. pairwise expression context;
8. WGCNA expression-network composite;
9. DoRothEA TF-association divergence;
10. TF-expression-profile divergence;
11. JASPAR/FIMO promoter-motif divergence;
12. promoter-methylation context;
13. functional-microniche coverage;
14. expression-context subcoverage;
15. regulatory-network subcoverage;
16. methylation coverage;
17. WGCNA TOM divergence;
18. WGCNA module divergence;
19. WGCNA kME divergence;
20. promoter-methylation profile divergence;
21. conditional target-promoter hypomethylation support.

For each gene and ring, the plotted value is the maximum observed value across all
associated pair-by-cancer rows. Supplementary Table S47 also retains the median,
minimum, observed-row count, eligible-row count and evidence status. Missing or
non-eligible values remain `NA` and are drawn as hollow markers, not zero.

## Expression data

Only required gene columns are read from:

```text
data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
```

- Table S49 contains every Circos gene × colorectal/gastric/lung context, including
  explicit unavailable rows;
- Table S50 contains every observed model-level `log2(TPM+1)` value;
- when a gene/context has no expression value, S50 contains one sentinel row with
  `expression_log2_tpm_plus_1=NA`, `is_measurement=false` and an explicit reason;
- no unavailable expression is represented as numeric zero.

## Exact commands

Build base coordinates, links, rings and expression tables:

```bash
python -u scripts/build_genomic_circos_inputs.py \
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --promoters data/raw/regulatory/ensembl_promoters.tsv \
  --expression data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
  --models data/raw/depmap/Model.csv \
  --output-dir data/processed/circos
```

Add WGCNA, methylation and validation internal layers:

```bash
python -u scripts/enrich_genomic_circos_internal_layers.py \
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --wgcna data/processed/regulatory/wgcna/wgcna_pair_metrics_all_cancers.tsv \
  --coordinates data/processed/circos/genomic_circos_gene_coordinates.tsv \
  --rings data/processed/circos/genomic_circos_ring_values.tsv \
  --tracks data/processed/circos/genomic_circos_track_definitions.tsv
```

Complete expression availability and pair links:

```bash
python -u scripts/complete_genomic_circos_expression_summary.py
python -u scripts/complete_genomic_circos_links.py
```

Generate Figure S70 with all 35 tracks under strict layout:

```bash
MPLBACKEND=Agg \
python -u scripts/make_genomic_circos_figure_resilient.py \
  --config config/genomic_circos_asset.yaml \
  --coordinates data/processed/circos/genomic_circos_gene_coordinates.tsv \
  --links data/processed/circos/genomic_circos_pair_links.tsv \
  --ring-values data/processed/circos/genomic_circos_ring_values.tsv \
  --tracks data/processed/circos/genomic_circos_track_definitions.tsv \
  --output-root article_outputs \
  --strict-layout
```

Register and validate:

```bash
python -u scripts/register_genomic_circos_assets.py \
  --article-root article_outputs

python -u scripts/catalog_figure_source_data.py \
  --article-root article_outputs

python -u scripts/validate_genomic_circos_integrity.py \
  --article-root article_outputs
```

The canonical entrypoint performs every step automatically:

```bash
MPLBACKEND=Agg \
STRICT_LAYOUT=1 \
bash scripts/run_publication_pipeline.sh assets-only
```

## Supplementary tables

| Table | Contents |
|---|---|
| S45 | gene class, Ensembl identifiers and genomic coordinates |
| S46 | every NISE/paralog chord, including score-missing pairs |
| S47 | every gene-by-ring value, status and coverage count |
| S48 | all 35 ring definitions and source columns |
| S49 | complete gene × cancer expression summary |
| S50 | observed model-level expression plus explicit NA sentinels |
| S51 | complete generated catalogue of every pipeline script/module |
| S52 | paths, sizes and SHA-256 provenance of every Circos source |

The exact combined source table used to render S70 is copied byte-for-byte to:

```text
article_outputs/tables/figure_data/supplementary/Figure_S70_source_data.tsv
```

## Complete script documentation

```bash
python -u scripts/build_script_documentation.py
```

This scans every Python, Bash and R file in `scripts/` and `src/rses_onco/` and
writes:

```text
docs/SCRIPT_CATALOG.md
docs/script_manifest.tsv
data/processed/documentation/pipeline_script_catalog.tsv
```

Omission of any pipeline source causes a regression-test and integrity-validation
failure.
