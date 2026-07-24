# SynLeth-RSES-Onco / RSES-Onco

**RSES-Onco v0.12.0** is a coverage-aware multi-omics framework for discovering and prioritizing cancer-selective dependencies created by non-homologous isofunctional enzymes (NISEs), homologous paralogs, pathway backups, collateral deletions and downstream vulnerabilities. The initial disease scope is colorectal, gastric and lung cancer.

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

The repository combines curated NISE and paralog hypotheses, DepMap and TCGA evidence, integrated functional-loss states, CRISPR and RNAi perturbation evidence, multi-platform proteomics, expression and WGCNA, functional and regulatory networks, promoter methylation, pharmacology, structural context, explicit missingness and overlap control, genomic Circos visualization and a scripted publication package.

## Scientific boundary

The software prioritizes experimental hypotheses. It does not establish clinical efficacy, patient benefit, safety, treatment suitability or cure. Missing evidence is not converted to zero. Non-eligible domains do not enter the eligible score denominator. Exploratory context and pharmacology are kept separate from causal discovery evidence when a mechanistic mapping is unavailable.

## RSES-Onco domains

The seven top-level domains are preserved in v0.12.0:

| Top-level domain | Weight |
|---|---:|
| Tumor event | 0.16 |
| Conditional dependency | 0.22 |
| Selectivity | 0.14 |
| Expression compensation | 0.08 |
| Functional relation | 0.06 |
| Functional microniche | 0.16 |
| Validation and tractability | 0.18 |

The functional-microniche domain contains expression context (0.20), localization (0.15), biochemical/structural evidence (0.15), genetic phenotype (0.20), interaction network (0.15) and regulatory network (0.15).

Expression context is divided between pairwise expression divergence and cancer-specific signed WGCNA. Regulatory evidence includes DoRothEA TF associations, TF-expression consistency, JASPAR/FIMO promoter motifs and promoter methylation. Methylation remains inside the regulatory domain and is not counted as a new independent top-level domain.

## Extended multi-omics scoring in v0.12.0

New datasets do not become independent top-level domains merely because they are available. Correlated sources are normalized separately and combined inside existing biological domains.

### Layers entering the primary score

1. **Integrated functional loss**, inside `tumor_event`:
   - relative copy-number loss;
   - absolute copy number;
   - loss of heterozygosity;
   - damaging mutations;
   - clear variant-level loss-of-function observations.

   Overlapping events are consolidated into one model–gene loss state instead of being summed as independent evidence.

2. **CRISPR Gene Dependency**, inside `dependency`:
   - Chronos Gene Effect remains the magnitude channel;
   - Gene Dependency functions as a confidence channel.

3. **Protein compensation**, inside `expression_compensation`:
   - Gygi mass spectrometry;
   - Sanger mass spectrometry;
   - Olink;
   - CCLE RPPA;
   - MCLP RPPA500.

   Platforms are normalized and contrasted separately. Raw values from different technologies are never directly pooled.

4. **DEMETER2 RNAi**, inside the genetic-phenotype microniche:
   - provides orthogonal perturbation support;
   - does not receive a second full dependency-domain weight.

### Context and translation layers not directly scored

The following sources are standardized, exported, summarized and visualized but do not automatically increase the primary score without a causal mapping:

- metabolomics without a gene–reaction–metabolite map;
- miRNA expression without a curated miRNA–target relation;
- global chromatin profiling;
- ssGSEA pathway states;
- molecular subtypes and genomic signatures;
- MetMap phenotypes;
- hotspot mutations without loss-of-function annotation;
- fusions without breakpoint-level evidence of gene disruption;
- single-drug and combination-drug response.

This separation prevents circular prioritization and preserves the distinction between discovery evidence, biological interpretation and translational validation.

## Extended-layer robustness

The complete ranking retains the pre-extension baseline and calculates leave-one-extended-layer-out ablations for:

```text
integrated functional loss
CRISPR dependency probability
protein compensation
RNAi orthogonal support
```

The final ranking includes baseline scores, extended scores, rank changes and cancer-specific ranks under each ablation. Supplementary Figure S78 reports Spearman rank stability relative to the complete v0.12.0 method.

## Genomic Circos

Supplementary Figure S70 is a GRCh38 Circos representation containing every coordinate-complete simple NISE and homologous-paralog hypothesis. Every simple pair receives one chord, including candidates whose score is unavailable.

- red chords: NISE relationships;
- black chords: homologous-paralog relationships;
- purple ticks: genes represented in both classes;
- 35 score and evidence rings;
- hollow markers: missing or non-eligible evidence, never numeric zero.

The Circos stage exports coordinates, pair links, ring values, track definitions, complete model-level expression, explicit missing-expression sentinels, source provenance and final status under `data/processed/circos/`.

See [`docs/GENOMIC_CIRCOS_WORKFLOW_V0112.md`](docs/GENOMIC_CIRCOS_WORKFLOW_V0112.md).

## Input directories

Core DepMap matrices remain under:

```text
data/raw/depmap/
```

The additional DepMap custom-download files are placed under:

```text
dmap_data/
```

Canonical filenames and scientific roles are declared in:

```text
config/extended_multiomics_sources.yaml
```

Raw files are never modified. The pipeline records size, SHA-256, standardization status, source role and whether the source is eligible for direct scoring.

## Installation

```bash
conda env create -f environment.yml
conda activate rses-onco
python -m pip install -e .
python -m pytest -q -p no:cacheprovider
```

For an existing environment:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco
python -m pip install -e .
```

## Runtime checks

```bash
bash scripts/resume_wgcna_regulatory_pipeline.sh check-runtime
bash scripts/run_extended_multiomics_pipeline.sh check-runtime
```

The regulatory runtime requires R/WGCNA and FIMO from the MEME Suite. Document rendering uses LibreOffice, Poppler, Pandoc and Graphviz.

## Complete v0.12.0 execution

The canonical command executes WGCNA, regulatory networks, promoter motifs, methylation, baseline scoring, multi-omics integration, pharmacology, structural analyses, all figures, all tables, workbooks, editable documents, manifests, tests and checksums:

```bash
STRICT_LAYOUT=1 \
STRICT_EXTENDED_SOURCES=1 \
EXTENDED_DATA_DIR=dmap_data \
bash scripts/run_complete_v0120_pipeline.sh
```

Run it in a persistent `tmux` session:

```bash
SESSION="rses_v0120_$(date +%Y%m%d_%H%M%S)"

tmux new-session -d -s "$SESSION" \
  "bash -lc 'cd /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010 && source \"\$(conda info --base)/etc/profile.d/conda.sh\" && conda activate rses-onco && STRICT_LAYOUT=1 STRICT_EXTENDED_SOURCES=1 EXTENDED_DATA_DIR=dmap_data bash scripts/run_complete_v0120_pipeline.sh; status=\$?; echo RUNNER_EXIT_CODE=\$status; exec bash'"

tmux attach -t "$SESSION"
```

The shell remains open after success or failure and displays `RUNNER_EXIT_CODE`.

## Staged extended execution

Build and score the extended layer:

```bash
bash scripts/run_extended_multiomics_pipeline.sh analysis
```

Generate and register only Figures S71–S78 and Tables S53–S64:

```bash
bash scripts/run_extended_multiomics_pipeline.sh publication
```

Run both stages:

```bash
bash scripts/run_extended_multiomics_pipeline.sh all
```

## Main extended outputs

```text
data/processed/extended_multiomics/
├── extended_multiomics_source_inventory.tsv
├── extended_multiomics_source_status.tsv
├── functional_loss_states.tsv
├── extended_pair_evidence_by_cancer.tsv
├── proteomics_pair_evidence_by_source.tsv
├── extended_covariate_context.tsv
├── custom_drug_sensitivity_long.tsv
├── gdsc_combination_evidence_long.tsv
├── extended_multiomics_source_provenance.tsv
├── extended_multiomics_status.json
└── extended_multiomics_integrity_validation.json
```

The rankings are updated atomically and retain their baseline columns:

```text
results/expanded_26Q1/full/expanded_rses_onco.tsv
results/expanded_26Q1/depmap_only/expanded_rses_onco.tsv
```

## New supplementary assets

```text
Figure S71  Source availability and scientific role
Figure S72  Baseline versus extended ranking
Figure S73  Candidate-level scored-layer contributions
Figure S74  Cross-platform proteomic compensation
Figure S75  Integrated functional-loss states
Figure S76  CRISPR and RNAi orthogonal support
Figure S77  Non-causal context-layer shifts
Figure S78  Leave-one-extended-layer-out rank stability

Tables S53–S64
```

Every generated figure has PNG, PDF and SVG exports, an exact source TSV, generator script, input list, reproduction command and layout audit.

## Publication asset contract

```text
8 main figures
78 supplementary figures
86 registered figures
258 PNG/PDF/SVG files
4 main tables
64 supplementary tables
68 registered tables
```

## Document and package validation

```bash
STRICT_EXTENDED_SOURCES=1 \
bash scripts/verify_complete_article_run_v0120.sh
```

This checks scientific integrity, complete source standardization, score version, ablations, figure/table counts, layout status, rendered documents, SHA-256 manifests, incomplete GDC downloads and the complete software test suite.

Automated validation does not replace manual visual inspection. All figures and rendered document pages must be reviewed at 100% zoom using the generated checklist.

## Complete documentation

- [`docs/EXTENDED_MULTIOMICS_V0120.md`](docs/EXTENDED_MULTIOMICS_V0120.md) — v0.12.0 scientific integration, outputs and commands;
- [`docs/END_TO_END_ARTICLE_PROTOCOL.md`](docs/END_TO_END_ARTICLE_PROTOCOL.md) — core command-by-command workflow;
- [`docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md`](docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md) — source acquisition, provenance and recovery;
- [`docs/METHYLATION_DATA_AND_SCORING_V0111.md`](docs/METHYLATION_DATA_AND_SCORING_V0111.md) — methylation formulas and missingness;
- [`docs/GENOMIC_CIRCOS_WORKFLOW_V0112.md`](docs/GENOMIC_CIRCOS_WORKFLOW_V0112.md) — genomic Circos workflow;
- [`docs/GENOMIC_CIRCOS_IMPLEMENTATION_CONTRACT_V0112.md`](docs/GENOMIC_CIRCOS_IMPLEMENTATION_CONTRACT_V0112.md) — Circos completeness contract;
- [`docs/SCRIPT_CATALOG.md`](docs/SCRIPT_CATALOG.md) — generated code catalogue;
- [`supplementary/Supplementary_Methods_RSES_Onco_v0110.md`](supplementary/Supplementary_Methods_RSES_Onco_v0110.md) — scientific methods and formulas;
- [`manuscript/RSES_Onco_intro_methods_draft_v0110.md`](manuscript/RSES_Onco_intro_methods_draft_v0110.md) — editable Introduction and Materials and Methods draft.

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

## License

MIT for code. Third-party data retain their original licenses and terms.
