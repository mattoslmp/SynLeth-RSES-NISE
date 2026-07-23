# SynLeth-RSES-Onco / RSES-Onco

**RSES-Onco v0.11.1** is a coverage-aware framework for discovering and
prioritizing cancer-selective dependencies created by non-homologous
isofunctional enzymes (NISEs), homologous paralogs, pathway backups, collateral
deletions and downstream vulnerabilities. The initial disease scope is
colorectal, gastric and lung cancer.

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

The repository combines curated NISE and paralog hypotheses, DepMap/TCGA evidence,
functional and regulatory networks, optional GDC/TCGA promoter methylation,
pharmacology, structural context, explicit missingness and overlap control, and a
scripted publication package.

## Scientific boundary

The software prioritizes experimental hypotheses. It does not establish clinical
efficacy, patient benefit, safety, treatment suitability or cure. Missing evidence
is not converted to zero and non-eligible domains do not enter the eligible score
denominator. Promoter methylation beta values provide epigenetic context and are not
direct proof of transcriptional silencing.

## Current score domains

The global RSES-Onco weights remain unchanged:

```text
tumor event                 0.16
dependency                  0.22
selectivity                 0.14
expression compensation     0.08
functional relation         0.06
functional microniche       0.16
validation/tractability     0.18
```

The functional-microniche domain contains:

```text
expression context          0.20
localization                0.15
biochemical/structural      0.15
genetic phenotype           0.20
interaction network         0.15
regulatory network          0.15
```

When GDC promoter methylation is available, it is integrated **inside** the existing
`expression_compensation` domain rather than receiving a new independent global
weight:

```text
expression compensation                  0.70
promoter methylation context             0.30
```

The directional promoter-methylation context subscore is:

```text
lost-gene median promoter beta × (1 − target-gene median promoter beta)
```

This construction supports a context in which the hypothesized lost gene is more
promoter-methylated while the target/paralog remains less methylated. The result is
coverage-adjusted and must be interpreted together with expression evidence.

## Installation

```bash
conda env create -f environment.yml
conda activate rses-onco
python -m pip install -e .
python -m pytest -q -p no:cacheprovider
```

## Complete execution and data acquisition

The repository contains a complete protocol from source acquisition to final article
packaging:

- [`docs/END_TO_END_ARTICLE_PROTOCOL.md`](docs/END_TO_END_ARTICLE_PROTOCOL.md) — canonical command-by-command pipeline tutorial;
- [`docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md`](docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md) — source acquisition, provenance, validation and recovery;
- [`docs/METHYLATION_INTEGRATION_V0111.md`](docs/METHYLATION_INTEGRATION_V0111.md) — GDC methylation acquisition, promoter aggregation, formulas, scoring semantics and rerun commands;
- [`supplementary/Supplementary_Methods_RSES_Onco_v0110.md`](supplementary/Supplementary_Methods_RSES_Onco_v0110.md) — scientific methods, formulas, evidence rules and references;
- [`supplementary/Supplementary_Methylation_Methods_RSES_Onco_v0111.md`](supplementary/Supplementary_Methylation_Methods_RSES_Onco_v0111.md) — methylation-specific supplementary methods and references;
- [`manuscript/RSES_Onco_intro_methods_draft_v0110.md`](manuscript/RSES_Onco_intro_methods_draft_v0110.md) — editable Introduction and Materials and Methods draft;
- [`docs/figures/RSES_Onco_workflow_and_applications.svg`](docs/figures/RSES_Onco_workflow_and_applications.svg) — vector workflow and practical-application figure.

The repository versions the corresponding DOCX, PDF and PNG derivatives. They are
generated reproducibly with:

```bash
bash scripts/generate_repository_documentation_assets.sh
```

## Methylation modes

The regulatory resume pipeline supports four modes:

```text
auto       use an existing promoter-methylation matrix; otherwise continue without it
download   query/download/validate GDC methylation and integrate it
require    require a usable methylation layer and fail if it cannot be produced
off        explicitly disable methylation
```

Run the complete methylation-aware workflow:

```bash
METHYLATION_MODE=download \
MPLBACKEND=Agg \
STRICT_LAYOUT=1 \
bash scripts/resume_wgcna_regulatory_pipeline.sh resume-regulatory
```

Build only the methylation evidence layer:

```bash
METHYLATION_MODE=download \
bash scripts/resume_wgcna_regulatory_pipeline.sh methylation-only
```

The DepMap-only ranking remains free of TCGA/GDC methylation. Methylation is added
only to the integrated TCGA–DepMap ranking.

## Publication workflow

Rebuild all publication assets from cached evidence and structure renders:

```bash
MPLBACKEND=Agg \
STRICT_LAYOUT=1 \
bash scripts/run_publication_pipeline.sh assets-only
```

Build and render the editable Word documents and inspection PDFs:

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
69 supplementary figures
77 registered figures
231 PNG/PDF/SVG files
4 main tables
47 supplementary tables
51 registered tables
```

Every registered figure has an exact source TSV, generator script, input list,
reproduction command, layout audit and PNG/PDF/SVG export. Figures S39-S69 provide
source-backed score decomposition, robustness, expression, coexpression, WGCNA,
regulatory, protein-network, CRISPR, localization, biochemical/structural, TCGA,
pharmacology, NISE/paralog, control, stability and promoter-methylation support.
Unavailable optional evidence is displayed explicitly and is never invented.

Supplementary Figure S69 reports cancer-specific promoter methylation context.
Supplementary Tables S45-S47 retain pair-level methylation evidence, gene-level
promoter summaries and source/aggregation status.

The document pipeline creates editable DOCX files, rendered PDFs and page PNGs.
Supplementary Figures S68 and S69 are forced onto separate pages and verified from
the rendered PDF.

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
- [`docs/PUBLICATION_COMPLETENESS_V0110.md`](docs/PUBLICATION_COMPLETENESS_V0110.md)
- [`docs/WGCNA_CORRELATION_POLICY_V0109.md`](docs/WGCNA_CORRELATION_POLICY_V0109.md)
- [`docs/STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md`](docs/STRING_FUNCTIONAL_EVIDENCE_WORKFLOW.md)
- [`docs/DOROTHEA_RECOVERY_WORKFLOW.md`](docs/DOROTHEA_RECOVERY_WORKFLOW.md)
- [`docs/STRUCTURAL_ATLAS_WORKFLOW.md`](docs/STRUCTURAL_ATLAS_WORKFLOW.md)

## License

MIT for code. Third-party data retain their original licenses and terms.
