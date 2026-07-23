# SynLeth-RSES-Onco / RSES-Onco

**RSES-Onco v0.11.0** is a coverage-aware framework for discovering and
prioritizing cancer-selective dependencies created by non-homologous
isofunctional enzymes (NISEs), homologous paralogs, pathway backups, collateral
deletions and downstream vulnerabilities. The initial disease scope is
colorectal, gastric and lung cancer.

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

The repository combines curated NISE and paralog hypotheses, DepMap/TCGA evidence,
functional and regulatory networks, pharmacology, structural context, explicit
missingness and overlap control, and a completely scripted publication package.

## Scientific boundary

The software prioritizes experimental hypotheses. It does not establish clinical
efficacy, patient benefit, safety, treatment suitability or cure. Missing evidence
is not converted to zero and non-eligible domains do not enter the eligible score
denominator.

## Installation

```bash
conda env create -f environment.yml
conda activate rses-onco
python -m pip install -e .
python -m pytest -q -p no:cacheprovider
```

## Complete execution and data acquisition

The repository contains a complete protocol from source acquisition to final article packaging:

- [`docs/END_TO_END_ARTICLE_PROTOCOL.md`](docs/END_TO_END_ARTICLE_PROTOCOL.md) — canonical command-by-command pipeline tutorial;
- [`docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md`](docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md) — source acquisition, provenance, validation and recovery;
- [`supplementary/Supplementary_Methods_RSES_Onco_v0110.md`](supplementary/Supplementary_Methods_RSES_Onco_v0110.md) — scientific methods, formulas, evidence rules and references;
- [`manuscript/RSES_Onco_intro_methods_draft_v0110.md`](manuscript/RSES_Onco_intro_methods_draft_v0110.md) — editable Introduction and Materials and Methods draft;
- [`docs/figures/RSES_Onco_workflow_and_applications.svg`](docs/figures/RSES_Onco_workflow_and_applications.svg) — vector workflow and practical-application figure.

The repository versions the corresponding DOCX, PDF and PNG derivatives. They are generated reproducibly with:

```bash
bash scripts/generate_repository_documentation_assets.sh
```

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
44 supplementary tables
48 registered tables
```

Every registered figure has an exact source TSV, generator script, input list,
reproduction command, layout audit and PNG/PDF/SVG export. Figures S39-S69 provide
source-backed score decomposition, robustness, expression, coexpression, WGCNA,
regulatory, protein-network, CRISPR, localization, biochemical/structural, TCGA,
pharmacology, NISE/paralog, control and stability support. Unavailable optional
evidence is displayed explicitly and is never invented.

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
