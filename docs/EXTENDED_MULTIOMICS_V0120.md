# RSES-Onco v0.12.0 extended multi-omics workflow

## Scientific design

RSES-Onco v0.12.0 preserves the seven top-level domains and prevents the number of available datasets from becoming an artificial score advantage. New sources are integrated as internal subcomponents or retained as non-causal validation/context layers.

### Layers entering the primary discovery score

1. **Integrated functional loss** within `tumor_event`:
   - relative copy-number loss;
   - absolute copy number;
   - LOH;
   - damaging mutations;
   - variant-level loss-of-function observations when available.

   Overlapping events are consolidated into one model–gene loss state and are not summed as independent evidence.

2. **CRISPR Gene Dependency** within `dependency`:
   - Chronos Gene Effect remains the magnitude channel;
   - Gene Dependency is a confidence channel inside the same domain.

3. **Protein compensation** within `expression_compensation`:
   - Gygi MS;
   - Sanger MS;
   - Olink;
   - CCLE RPPA;
   - MCLP RPPA500.

   Every platform is normalized and contrasted separately. The score uses a coverage-aware median consensus; raw platform values are never directly pooled.

4. **DEMETER2 RNAi** within the genetic-phenotype microniche:
   - used as orthogonal perturbation evidence;
   - does not receive a second full dependency-domain weight.

### Layers retained outside the primary score

The following sources are exported, summarized and visualized but do not automatically increase the discovery score without an explicit mechanistic mapping:

- metabolomics without a curated gene–reaction–metabolite map;
- miRNA expression without a curated miRNA–target relation;
- global chromatin profiling;
- ssGSEA pathway states;
- molecular subtypes;
- genomic signatures;
- MetMap phenotypes;
- hotspot mutations without loss-of-function classification;
- fusions without breakpoint-level evidence of gene disruption;
- single-drug and combination-drug response.

This separation prevents circular prioritization and distinguishes **discovery evidence** from **translation and biological-context evidence**.

## Input directory

Place the DepMap custom-download files in:

```text
dmap_data/
```

The canonical file names are declared in:

```text
config/extended_multiomics_sources.yaml
```

Raw files are never modified. Every available input receives a SHA-256 value and a standardization status.

## Main commands

Activate the environment and update the editable package:

```bash
cd /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco
python -m pip install -e .
```

Check the complete runtime:

```bash
bash scripts/resume_wgcna_regulatory_pipeline.sh check-runtime
bash scripts/run_extended_multiomics_pipeline.sh check-runtime
```

Run the complete v0.12.0 workflow:

```bash
STRICT_LAYOUT=1 \
STRICT_EXTENDED_SOURCES=1 \
EXTENDED_DATA_DIR=dmap_data \
bash scripts/run_complete_v0120_pipeline.sh
```

Run inside tmux while keeping the shell open after completion or failure:

```bash
SESSION="rses_v0120_$(date +%Y%m%d_%H%M%S)"

tmux new-session -d -s "$SESSION" \
  "bash -lc 'cd /mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010 && source \"\$(conda info --base)/etc/profile.d/conda.sh\" && conda activate rses-onco && STRICT_LAYOUT=1 STRICT_EXTENDED_SOURCES=1 EXTENDED_DATA_DIR=dmap_data bash scripts/run_complete_v0120_pipeline.sh; status=\$?; echo RUNNER_EXIT_CODE=\$status; exec bash'"

tmux attach -t "$SESSION"
```

## Staged execution

Build and score only the extended layer:

```bash
bash scripts/run_extended_multiomics_pipeline.sh analysis
```

Generate and register only the extended publication assets:

```bash
bash scripts/run_extended_multiomics_pipeline.sh publication
```

Execute both:

```bash
bash scripts/run_extended_multiomics_pipeline.sh all
```

## Principal processed outputs

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

## Ranking outputs

The existing ranking paths are updated atomically and retain baseline columns:

```text
results/expanded_26Q1/full/expanded_rses_onco.tsv
results/expanded_26Q1/depmap_only/expanded_rses_onco.tsv
```

Important columns include:

```text
baseline_coverage_adjusted_rses
coverage_adjusted_rses
extended_score_delta
baseline_rank_within_cancer
extended_rank_within_cancer
extended_rank_change
extended_scored_layer_count
extended_scored_layer_coverage
score_version
```

The score version is:

```text
RSES-Onco-expanded-v0.12.0
```

## New publication assets

### Supplementary figures

```text
Figure S71  Extended source availability and scientific role
Figure S72  Baseline versus extended ranking
Figure S73  Candidate-level scored layer contributions
Figure S74  Cross-platform proteomic compensation
Figure S75  Integrated functional-loss states
Figure S76  CRISPR and RNAi orthogonal support
Figure S77  Non-causal multi-omics context shifts
Figure S78  Rank stability after extension
```

### Supplementary tables

```text
Tables S53-S64
```

These include source inventory/status, functional-loss states, pair-level evidence, platform-resolved proteomics, covariate context, standardized single-drug and combination-drug data, provenance, complete extended ranking, baseline-versus-extended comparison and scientific layer definitions.

Every figure is exported as PNG, PDF and SVG, and every figure source table is retained.

## Expected final package

```text
8 main figures
78 supplementary figures
86 registered figures
258 PNG/PDF/SVG exports
4 main tables
64 supplementary tables
68 registered tables
```

## Final validation

```bash
STRICT_EXTENDED_SOURCES=1 \
bash scripts/verify_complete_article_run_v0120.sh
```

Automated validation does not replace visual inspection. Every figure and rendered document page must still be inspected at 100% zoom using the generated manual checklist.
