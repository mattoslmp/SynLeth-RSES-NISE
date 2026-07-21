# SynLeth-RSES / RSES-Onco

**RSES-Onco** is a reproducible computational module for identifying cancer-selective synthetic-lethal dependencies created by hidden functional backups, non-homologous isofunctional enzymes, homologous paralogs, pathway redundancy and collateral deletion. The initial scope is colorectal, gastric and lung cancer.

The design adapts the coverage-aware Role Specialization Evidence Score (RSES): missing evidence is omitted from the observed-domain denominator and reduces explicit coverage rather than being treated as biological similarity. Reaction/evolutionary validity remains separate from synthetic-lethality evidence.

## Included

- A curated table of **70 human proteins in 15 bona fide intragenomic analogous-enzyme activities** from Piergiorge et al. (2017).
- Acquisition script for the 2026 Swiss-Prot NISE dataset and Homo sapiens extraction.
- DepMap readers for CRISPR gene effect, copy number, model metadata and expression.
- GDC public copy-number acquisition script for TCGA-COAD, TCGA-READ, TCGA-STAD, TCGA-LUAD and TCGA-LUSC.
- Coverage-aware RSES-Onco scoring.
- A transparent benchmark of **25 literature-anchored vulnerabilities and discovery hypotheses**, including MTAP-PRMT5/MAT2A, MSI-WRN, HRD-PARP1/POLQ/USP1, SMARCA4-SMARCA2, dual-SMARCA4/2-MCL1, ENO1-ENO2, ME2-ME3, VPS4B-VPS4A, STAG2-STAG1 and high-alkylation-NTHL1/CHEK1.
- Exploratory DNA-repair NISE hypotheses involving NTHL1, OGG1, NEIL1/2 and APEX1/2.
- Manuscript, supplementary material, an interactive supplementary workbook, publication figures in PDF/PNG/SVG, tests and GitHub Actions CI.

## Scientific boundary

The bundled ranked table is a **literature-anchored pilot**, not a claim that the complete current TCGA and DepMap matrices were processed inside the distribution. Run the empirical command with the official matrices to produce cohort- and release-specific results. Computational candidates require experimental validation before pharmacological or clinical interpretation.

## Installation

```bash
conda env create -f environment.yml
conda activate rses-onco
pip install -e .
pytest -q
```

Alternatively:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e . pytest
pytest -q
```

## Reproduce the bundled pilot

```bash
rses-onco score-literature   --input data/curated/synthetic_lethality_reference_pairs.tsv   --output results/literature_anchored_candidates.tsv
```

## Obtain current human NISE records

```bash
python scripts/download_human_nise.py
```

This downloads `SwissProt_NISE.tsv.gz` from Zenodo record `18008936` and filters Homo sapiens. The bundled 2017 table is retained because it is a specifically curated, bona fide human set with 70 proteins and 15 EC activities.

## Prepare DepMap

Download the current public release from the official DepMap downloads page. The code is configured for the 26Q1 naming convention but accepts documented aliases.

```bash
python scripts/download_depmap.py --input-dir data/raw/depmap
```

Expected files:

- `CRISPRGeneEffect.csv`
- `OmicsCNGene.csv`
- `Model.csv`
- `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` (optional but recommended)

## Obtain public GDC copy-number files

```bash
python scripts/download_gdc.py --manifest-only
python scripts/download_gdc.py
```

The first command writes a complete query manifest without downloading. Review it before bulk transfer.

## Run empirical analysis

```bash
python scripts/run_empirical_rses_onco.py   --gene-effect data/raw/depmap/CRISPRGeneEffect.csv   --copy-number data/raw/depmap/OmicsCNGene.csv   --models data/raw/depmap/Model.csv   --expression data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv   --tcga colon=data/processed/TCGA_COADREAD_discrete_cna.tsv   --tcga stomach=data/processed/TCGA_STAD_discrete_cna.tsv   --tcga lung=data/processed/TCGA_LUNG_discrete_cna.tsv   --output results/empirical_rses_onco.tsv
```

Complex biomarkers such as MSI, MMR deficiency and HRD should be supplied as explicit model annotations in a future extension; they must not be inferred from a single copy-number column.

## Build a sequence tree

```bash
python scripts/build_sequence_tree.py --threads 8
```

This retrieves UniProt sequences, aligns with MAFFT and builds a FastTree Newick tree. A whole-set tree is exploratory because proteins from different EC activities are not expected to be globally alignable; activity-specific trees are recommended for interpretation.

## RSES-Onco score

For available domains `d` with weights `w_d`, values `D_d` and availability indicators `m_d`:

`RSES-Onco = sum(w_d m_d D_d) / sum(w_d m_d)`

`Coverage = sum(w_d m_d) / sum(w_d)`

`Adjusted score = RSES-Onco x Coverage`

The default domains are tumor event, target dependency, loss-selectivity, expression compensation, functional-relation confidence, and validation/tractability.

## Priority validation panels

- **Colorectal:** MSI/MMR loss-WRN; MTAP-PRMT5/MAT2A; HRD-PARP1/POLQ/USP1; NTHL1-centered base-excision repair hypotheses.
- **Gastric:** MSI-WRN; ARID1A-ATR; MTAP-PRMT5/MAT2A; HRD-related targets.
- **Lung:** SMARCA4-SMARCA2; SMARCA4-CDK4/6; dual-SMARCA4/2-MCL1; MTAP-PRMT5/MAT2A; ATM-ATR; HRD-USP1/POLQ.

## Main references

- Piergiorge RM et al. Genome Biology and Evolution. 2017. DOI: 10.1093/gbe/evx119.
- de Oliveira FC et al. BMC Research Notes. 2026. DOI: 10.1186/s13104-026-07742-5; Zenodo: 10.5281/zenodo.18008936.
- Kryukov GV et al. Science. 2016. DOI: 10.1126/science.aad5214.
- Chan EM et al. Nature. 2019. DOI: 10.1038/s41586-019-1102-x.
- Bryant HE et al. Nature. 2005. DOI: 10.1038/nature03443.
- Ceccaldi R et al. Nature. 2015. DOI: 10.1038/nature14184.
- Muller FL et al. Nature. 2012. DOI: 10.1038/nature11331.
- Dey P et al. Nature. 2017. DOI: 10.1038/nature21052.

## License

MIT for code. Third-party datasets remain under their original terms and are acquired by scripts rather than silently redistributed.
