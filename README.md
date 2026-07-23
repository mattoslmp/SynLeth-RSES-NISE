# SynLeth-RSES-Onco / RSES-Onco

**RSES-Onco v0.11.2** is a coverage-aware framework for discovering and
prioritizing cancer-selective dependencies created by non-homologous
isofunctional enzymes (NISEs), homologous paralogs, pathway backups, collateral
deletions and downstream vulnerabilities. The initial disease scope is
colorectal, gastric and lung cancer.

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

The repository combines curated NISE and paralog hypotheses, DepMap/TCGA evidence,
functional and regulatory networks, promoter methylation, pharmacology, structural
context, explicit missingness and overlap control, a genomic Circos representation,
and a fully scripted publication package.

## Scientific boundary

The software prioritizes experimental hypotheses. It does not establish clinical
efficacy, patient benefit, safety, treatment suitability or cure. Missing evidence
is not converted to zero and non-eligible domains do not enter the eligible score
denominator.

## What enters the RSES-Onco calculation

| Top-level domain | Weight |
|---|---:|
| Tumor event | 0.16 |
| Conditional dependency | 0.22 |
| Selectivity | 0.14 |
| Expression compensation | 0.08 |
| Functional relation | 0.06 |
| Functional microniche | 0.16 |
| Validation and tractability | 0.18 |

The functional-microniche domain contains expression context (0.20), localization
(0.15), biochemical/structural evidence (0.15), genetic phenotype (0.20),
interaction network (0.15) and regulatory network (0.15).

Expression context is divided between pairwise expression divergence (0.50) and
cancer-specific signed WGCNA context (0.50). The regulatory-network subscore is
divided into DoRothEA TF-association divergence (0.32), TF-expression-profile
divergence (0.28), JASPAR/FIMO promoter-motif divergence (0.20) and
promoter-methylation context (0.20). Methylation is not counted as a new independent
full RSES-Onco domain.

## Genomic Circos

The canonical publication pipeline generates **Supplementary Figure S70**, a
GRCh38 Circos representation containing every coordinate-complete simple NISE and
homologous-paralog hypothesis. Every simple candidate pair receives exactly one
chord, including pairs whose score is unavailable; those pairs retain
`link_status=score_missing` instead of disappearing.

- red chords: NISE relationships;
- black chords: homologous-paralog relationships;
- purple genomic ticks: genes represented in both classes;
- **35 rings** in total;
- Panel A: observed and coverage-adjusted RSES-Onco, evidence coverage, all seven
  top-level domains and four individual validation/tractability terms;
- Panel B: all six microniche domains, pairwise expression, WGCNA composite, TOM,
  module and kME divergence, DoRothEA, TF-expression consistency, JASPAR/FIMO
  motifs, promoter-methylation composite, methylation-profile divergence,
  conditional target hypomethylation and nested coverage values;
- hollow ring markers: missing or non-eligible evidence, never numeric zero.

The Circos stage exports:

```text
data/processed/circos/genomic_circos_gene_coordinates.tsv
data/processed/circos/genomic_circos_pair_links.tsv
data/processed/circos/genomic_circos_ring_values.tsv
data/processed/circos/genomic_circos_track_definitions.tsv
data/processed/circos/genomic_circos_expression_summary.tsv
data/processed/circos/genomic_circos_expression_model_values.tsv
data/processed/circos/genomic_circos_source_provenance.tsv
data/processed/circos/genomic_circos_status.json
```

All observed model-level DepMap expression values used for Circos genes are retained
in Supplementary Table S50. A gene/cancer context without an observed expression
value receives an explicit sentinel row with expression `NA`,
`is_measurement=false` and a documented absence reason; missing expression is never
represented as zero. The exact combined source TSV used for Figure S70 is copied to
`article_outputs/tables/figure_data/supplementary/Figure_S70_source_data.tsv`.

The finalized status and source-provenance contract records 35 rings, all pair
chords, missing-score chords, observed expression measurements, unavailable
expression sentinels, and SHA-256 values for the final ranking, candidate universe,
Ensembl coordinates, DepMap expression/models and WGCNA pair metrics.

See [`docs/GENOMIC_CIRCOS_WORKFLOW_V0112.md`](docs/GENOMIC_CIRCOS_WORKFLOW_V0112.md).

## Complete script documentation

Every Python, Bash and R source under `scripts/` and `src/rses_onco/` is catalogued
automatically from the repository source code. The generated outputs are:

```text
docs/SCRIPT_CATALOG.md
docs/script_manifest.tsv
data/processed/documentation/pipeline_script_catalog.tsv
```

The catalogue records purpose, language, pipeline stage, command, CLI options,
declared paths, line count and SHA-256 for every script/module. A test and the
Circos integrity validator fail if any pipeline source is omitted.

## Methylation input

Use the official DepMap custom-download dataset **Methylation (1kb upstream TSS)**
or the traceable historical CCLE RRBS file. Set the path explicitly when necessary:

```bash
export METHYLATION="$DEPMAP_DIR/Methylation_(1kb_upstream_TSS)_subsetted_NAsdropped.csv"
```

See [`docs/METHYLATION_DATA_AND_SCORING_V0111.md`](docs/METHYLATION_DATA_AND_SCORING_V0111.md).

## Installation

```bash
conda env create -f environment.yml
conda activate rses-onco
python -m pip install -e .
python -m pytest -q -p no:cacheprovider
```

## Complete execution and data acquisition

- [`docs/END_TO_END_ARTICLE_PROTOCOL.md`](docs/END_TO_END_ARTICLE_PROTOCOL.md) — canonical command-by-command pipeline tutorial;
- [`docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md`](docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md) — source acquisition, provenance, validation and recovery;
- [`docs/METHYLATION_DATA_AND_SCORING_V0111.md`](docs/METHYLATION_DATA_AND_SCORING_V0111.md) — methylation source, formulas and missingness;
- [`docs/GENOMIC_CIRCOS_WORKFLOW_V0112.md`](docs/GENOMIC_CIRCOS_WORKFLOW_V0112.md) — genomic coordinates, all pair links, 35 rings, expression tables and exact commands;
- [`docs/GENOMIC_CIRCOS_IMPLEMENTATION_CONTRACT_V0112.md`](docs/GENOMIC_CIRCOS_IMPLEMENTATION_CONTRACT_V0112.md) — executable completeness contract;
- [`docs/SCRIPT_CATALOG.md`](docs/SCRIPT_CATALOG.md) — generated complete code catalogue;
- [`supplementary/Supplementary_Methods_RSES_Onco_v0110.md`](supplementary/Supplementary_Methods_RSES_Onco_v0110.md) — scientific methods, formulas, evidence rules and references;
- [`manuscript/RSES_Onco_intro_methods_draft_v0110.md`](manuscript/RSES_Onco_intro_methods_draft_v0110.md) — editable Introduction and Materials and Methods draft.

## Publication workflow

```bash
MPLBACKEND=Agg \
STRICT_LAYOUT=1 \
bash scripts/run_publication_pipeline.sh assets-only
```

Build and render the editable documents:

```bash
bash scripts/run_publication_pipeline.sh documents
```

Validate the completed package:

```bash
bash scripts/verify_complete_article_run.sh
```

## Publication asset contract

```text
8 main figures
70 supplementary figures
78 registered figures
234 PNG/PDF/SVG files
4 main tables
52 supplementary tables
56 registered tables
```

Supplementary Tables S45-S52 contain the Circos coordinates, every pair link, all
35 ring values, track definitions, complete expression summary, observed
model-level expression plus explicit NA sentinels, complete script catalogue and
SHA-256 source provenance. Every registered figure has an exact source TSV,
generator script, input list, reproduction command, layout audit and PNG/PDF/SVG
export.

The document pipeline creates editable DOCX files, rendered PDFs and page PNGs.
Every supplementary figure starts on a separate page.

## Output structure

```text
article_outputs/
├── figures/main/
├── figures/supplementary/
├── structure_atlas/individual/
├── tables/main/
├── tables/supplementary/
├── tables/qc/
├── tables/score_components/
├── tables/robustness/
├── tables/figure_data/
├── tables/supporting_evidence/
├── source_data/
├── manuscript_assets/
├── workbooks/
├── documents/
├── review_records/
└── manifests/
```

## Additional documentation

- [`docs/PUBLICATION_EVIDENCE_AUDIT_AND_REPRODUCTION.md`](docs/PUBLICATION_EVIDENCE_AUDIT_AND_REPRODUCTION.md)
- [`docs/WGCNA_CORRELATION_POLICY_V0109.md`](docs/WGCNA_CORRELATION_POLICY_V0109.md)
- [`docs/STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md`](docs/STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md)
- [`docs/DOROTHEA_RECOVERY_WORKFLOW.md`](docs/DOROTHEA_RECOVERY_WORKFLOW.md)
- [`docs/STRUCTURAL_ATLAS_WORKFLOW.md`](docs/STRUCTURAL_ATLAS_WORKFLOW.md)

## License

MIT for code. Third-party data retain their original licenses and terms.
