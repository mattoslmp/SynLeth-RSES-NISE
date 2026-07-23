# RSES-Onco v0.11.2 genomic Circos workflow

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

## Purpose

The pipeline generates Supplementary Figure S70 as a source-backed genomic Circos
representation of every coordinate-complete simple-gene hypothesis classified as a
non-homologous isofunctional enzyme (NISE), a homologous paralog, or both. The
figure combines genomic position, pair relationships, all RSES-Onco score domains,
all functional-microniche subdomains, expression-network layers, regulatory layers,
promoter methylation and nested evidence coverage.

## Genomic coordinates

Canonical coordinates are read from the Ensembl promoter acquisition table:

```text
data/raw/regulatory/ensembl_promoters.tsv
```

The plotted point uses the canonical-transcript TSS or the gene midpoint when the
latter is present. Chromosome labels are normalized to GRCh38 chromosomes 1-22, X
and Y. A gene without a coordinate is not silently placed at an arbitrary position;
the Circos input stage fails and lists the unresolved genes.

## Included hypotheses

The input universe is:

```text
data/processed/expanded_candidate_universe.tsv
```

The Circos stage includes every simple pair whose `source_class`, `relation_type`
or `ensembl_homology_type` identifies a NISE or paralog/homolog relationship.
Composite features that cannot be represented as two genomic genes are not placed
artificially on the chromosome ideogram.

## Link representation

- NISE links: red (`#C62828`);
- homologous-paralog links: black (`#111111`);
- genes participating in both relationship classes: purple genomic ticks;
- chord width and transparency: proportional to the maximum cancer-specific
  `coverage_adjusted_rses` observed for the pair;
- every coordinate-complete NISE/paralog pair is retained, including pairs with
  missing score components.

## Ring panels

### Panel A — top-level RSES-Onco

1. coverage-adjusted RSES-Onco;
2. evidence coverage;
3. tumor event;
4. conditional dependency;
5. selectivity;
6. expression compensation;
7. functional relation;
8. functional microniche;
9. validation and tractability.

### Panel B — functional microniche and internal layers

1. expression context;
2. localization;
3. biochemical/structural evidence;
4. genetic phenotype;
5. interaction network;
6. regulatory network;
7. pairwise expression context;
8. WGCNA expression network;
9. DoRothEA TF-association divergence;
10. TF-expression-profile divergence;
11. JASPAR/FIMO promoter-motif divergence;
12. promoter-methylation context;
13. functional-microniche coverage;
14. expression-context subcoverage;
15. regulatory-network subcoverage;
16. methylation coverage.

For each gene and ring, the plotted value is the maximum observed value across all
associated pair-by-cancer rows. Supplementary Table S47 additionally retains the
median, minimum, number of observed rows, number of eligible rows and evidence
status. Missing or non-eligible evidence remains `NA` and is drawn as a hollow
marker rather than zero.

## Expression data

The Circos stage reads only the required gene columns from:

```text
data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
```

The following outputs ensure that every expression value used is available for
review and reproduction:

- Supplementary Table S49: cancer-by-gene expression summary;
- Supplementary Table S50: every model-level `log2(TPM+1)` value for every Circos
  gene in colorectal, gastric and lung models.

## Commands

Build all Circos source tables:

```bash
python -u scripts/build_genomic_circos_inputs.py \
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --promoters data/raw/regulatory/ensembl_promoters.tsv \
  --expression data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
  --models data/raw/depmap/Model.csv \
  --output-dir data/processed/circos
```

Generate Figure S70:

```bash
MPLBACKEND=Agg \
python -u scripts/make_genomic_circos_figure.py \
  --config config/genomic_circos_asset.yaml \
  --coordinates data/processed/circos/genomic_circos_gene_coordinates.tsv \
  --links data/processed/circos/genomic_circos_pair_links.tsv \
  --ring-values data/processed/circos/genomic_circos_ring_values.tsv \
  --tracks data/processed/circos/genomic_circos_track_definitions.tsv \
  --output-root article_outputs \
  --strict-layout
```

Register S70 and Tables S45-S52:

```bash
python -u scripts/register_genomic_circos_assets.py \
  --article-root article_outputs
```

The canonical entrypoint performs these steps automatically:

```bash
MPLBACKEND=Agg \
STRICT_LAYOUT=1 \
bash scripts/run_publication_pipeline.sh assets-only
```

## Supplementary tables

| Table | Contents |
|---|---|
| S45 | gene class, Ensembl identifiers and genomic coordinates |
| S46 | every NISE/paralog chord and link rendering value |
| S47 | every gene-by-ring value, status and coverage count |
| S48 | ring ID, label, source column, parent domain and aggregation rule |
| S49 | cancer-by-gene expression summary |
| S50 | complete model-level expression values used for Circos genes |
| S51 | complete generated catalogue of every pipeline script/module |
| S52 | paths, sizes and SHA-256 provenance of every Circos source |

The exact combined source table used to render Figure S70 is also copied to:

```text
article_outputs/tables/figure_data/supplementary/Figure_S70_source_data.tsv
```

## Complete script documentation

Run:

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

The catalogue is tested against the repository tree; omission of any script/module
causes a test failure.
