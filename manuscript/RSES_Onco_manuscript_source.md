# RSES-Onco: a coverage-aware framework for prioritizing cancer synthetic lethality from hidden functional enzyme backups

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratory of Bioinformatics and Data Science, WBPereira, Itaperuna, Rio de Janeiro, Brazil

This Markdown file is the editable source companion to `RSES_Onco_manuscript.docx`. The full formatted manuscript includes five main figures, tables, equations and references.

## Core formula

$RSES_{Onco} = \frac{\sum_d w_d m_d D_d}{\sum_d w_d m_d}$

$Coverage = \frac{\sum_d w_d m_d}{\sum_d w_d}$

$Adjusted\ score = RSES_{Onco} \times Coverage$

Missing domains are excluded from the observed-domain denominator and reduce coverage.

## Resource summary

- 70 human proteins in 15 curated bona fide activity groups.
- 101 cross-structural-cluster analogue pairs.
- 25 literature-anchored or discovery-hypothesis pairs.
- Initial disease scope: colorectal, gastric and lung cancer.
- Empirical data inputs: GDC/TCGA copy number or explicit LOF biomarkers, DepMap CRISPR gene effect/copy number/expression, and GTEx normal-tissue expression.

## Scientific boundary

The bundled rankings are literature priors and synthetic software verification. They are not presented as a completed current-release TCGA/DepMap discovery analysis.
