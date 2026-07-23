# RSES-Onco: a coverage-aware framework for prioritizing context-dependent cancer vulnerabilities from hidden functional enzyme backups

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

> This file is the editable manuscript source for the Introduction and Materials and Methods. Numerical statements about results must be synchronized with the final validated real-data outputs before submission.

## Introduction

Synthetic lethality and related context-dependent genetic interactions offer a route to identify vulnerabilities that are selectively exposed by defined molecular states. In cancer, these relationships can reveal dependencies created by gene loss, altered pathway activity, metabolic adaptation or network rewiring [1,2]. Genome-scale CRISPR screening and resources such as DepMap have greatly expanded the systematic search for such dependencies, while also highlighting the need to distinguish general essentiality from context-specific effects [3-5].

Most computational and experimental strategies emphasize homologous paralogs, direct pathway partners or canonical genetic interactions. Functional redundancy, however, is not restricted to homologous proteins. **Non-homologous isofunctional enzymes (NISEs)** can catalyze equivalent or closely related biochemical reactions despite limited sequence or structural homology. Such proteins may represent hidden alternatives that are missed by homology-centered searches and may become relevant only under particular cellular or tumor conditions.

The existence of the same nominal biochemical activity does not imply universal interchangeability. Enzymes can differ in expression, regulation, localization, cofactor usage, interaction partners, structural organization and response to genetic perturbation. Consequently, the biological relevance of a possible backup relationship may be restricted to a **functional microniche**: the local and context-dependent combination of molecular and phenotypic evidence within which a directional relationship becomes plausible. In the present framework, a functional microniche integrates expression and compensation, coexpression and WGCNA organization, regulatory and protein-network context, dependency profiles, localization and biochemical/structural evidence.

Integrating these layers is complicated by unequal coverage. Some candidate pairs are richly annotated, whereas others lack structures, regulatory evidence, event-positive cell lines or sufficient samples. Treating missing evidence as zero can confound absence of data with negative biological evidence. Conversely, averaging only observed components without a coverage penalty can favor sparsely characterized candidates. A rigorous prioritization system therefore requires explicit eligibility rules, missing-data states, source-overlap controls and a coverage-aware mathematical formulation.

We developed **RSES-Onco**, a reproducible framework for prioritizing cancer-context hypotheses derived from curated human NISEs, homologous paralogs, collateral-loss relationships and all-target conditional-dependency discovery. RSES-Onco integrates DepMap dependency, copy-number and expression data; TCGA/GDC tumor-event matrices; STRING functional networks; DoRothEA/OmniPath regulatory relationships; Ensembl promoters and JASPAR/FIMO motif predictions; Human Protein Atlas localization; biochemical and structural annotations; and pharmacology or tractability evidence [6-18].

The framework separates evidence used directly in scoring from prioritization, external validation and interpretation. It also records source overlap so that multiple representations of the same original observation do not receive independent full weights. The resulting score is not presented as proof of clinical efficacy or synthetic lethality. Instead, it provides an auditable ranking of hypotheses for CRISPR validation, rescue experiments, expression perturbation, biochemical follow-up, structural interpretation and exploratory pharmacology.

## Materials and Methods

### Computational framework and reproducibility

All analyses were implemented in the SynLeth-RSES-Onco repository. The canonical execution protocol is provided in `docs/END_TO_END_ARTICLE_PROTOCOL.md`, and source-specific acquisition and provenance procedures are provided in `docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md`. The software environment is defined by `environment.yml` and `pyproject.toml`. Publication assets are generated through `scripts/run_publication_pipeline.sh`.

The repository and publication framework are versioned as v0.11.0. The scientific ranking produced by the corrected scoring workflow is labelled `RSES-Onco-expanded-v0.10.9`, with expression/regulatory semantics `eligibility-aware-wgcna-regulatory-v3`.

### Candidate-universe construction

Curated human NISE activities and proteins were used to construct directed NISE hypotheses. For each cross-group pair, both directional hypotheses were represented separately because loss of gene A followed by dependence on gene B is not assumed to be equivalent to the reverse direction. Ensembl Compara paralogs were incorporated as a distinct mechanistic class [11]. Additional candidates included curated benchmark relationships, collateral-loss hypotheses and all-target conditional-dependency discovery candidates.

Each candidate retained a stable pair identifier, lost/source feature, target gene, source class, relationship type, tumor-context eligibility and available source-specific identifiers. Composite biomarkers were preserved as composite events unless the source supported traceable expansion into atomic targets.

### Cancer contexts

The primary disease contexts were colorectal, gastric and lung cancer. TCGA projects included COAD and READ for the colorectal context, STAD for the gastric context, and LUAD and LUSC for the lung context. Compatible DepMap cell lines were selected using the model annotation table. Actual release-specific model and tumor counts were recorded from the input files and were not forced to match previous releases.

### DepMap inputs and validation

DepMap evidence included CRISPR gene-effect scores, whole-genome copy-number data, model annotations and protein-coding expression. The expected matrices were `CRISPRGeneEffect.csv`, `OmicsCNGeneWGS.csv`, `Model.csv` and `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`. Input validation checked file availability, schema, model identifiers, gene-column mapping and cohort compatibility. Checksums were generated for release traceability.

### TCGA/GDC tumor-event processing

Gene-level copy-number files were acquired from GDC using a reviewed ASCAT3 manifest for TCGA-COAD, TCGA-READ, TCGA-STAD, TCGA-LUAD and TCGA-LUSC. Downloads were validated by expected file size and checksum, and incomplete `.part` files were excluded. Gene-level data were aggregated into cancer-context event matrices, including context-specific homozygous-deletion tables.

Tumor-event supporting tables retained event definition, cancer, number of samples, frequency and source. Frequencies were restricted to biologically valid values. Any jitter used in plotting was applied only to graphical positions and was not stored as the true event frequency.

### Dependency and selectivity analysis

DepMap gene-effect values were aligned to model annotations by `ModelID`. Candidate target dependency was evaluated within compatible cancer models and, where applicable, after stratification by the source/lost-gene event. General dependency and event-associated selectivity were treated as distinct concepts.

Supporting outputs retained the event-positive and event-intact sample sizes, group summaries, effect direction, effect-size estimate, nominal P value and FDR when computable. Candidates without sufficient compatible models were labelled as insufficient sample rather than assigned a zero component.

### Expression and compensation

Expression values were obtained from the DepMap protein-coding expression matrix using the release-defined transformation. For each directional candidate and cancer context, the workflow retained expression of both genes, paired-model count, event-positive and event-intact summaries and pairwise expression metrics when available.

Spearman correlation and median absolute expression difference were used to characterize expression-profile divergence. Correlation alone was not considered proof of compensation. Compensation required an operational event-stratified pattern consistent with preservation or increased expression of the target after the relevant source-gene event. Non-computable comparisons remained missing.

### Functional-microniche definition

The functional microniche represented the context-specific evidence environment of each candidate pair. It comprised six parent domains:

| Functional-microniche domain | Weight |
|---|---:|
| Expression context | 0.20 |
| Localization | 0.15 |
| Biochemical/structural evidence | 0.15 |
| Genetic phenotype | 0.20 |
| Interaction network | 0.15 |
| Regulatory network | 0.15 |

Missing eligible components reduced coverage. Non-eligible domains were excluded from the eligible denominator.

### Pairwise genetic-phenotype profile

Cancer-specific DepMap gene-effect profiles were compared using Spearman correlation, median absolute profile difference and divergence of dependency-positive model sets. Dependency-positive sets were defined by the threshold implemented in the source code, and set overlap was summarized by Jaccard similarity. The available divergence subcomponents were combined without replacing missing values by zero.

### WGCNA expression-network analysis

Cancer-specific signed WGCNA networks were calculated separately for colorectal, gastric and lung models [7,8]. Candidate genes were supplemented with cancer-specific highly variable genes. Genes had to satisfy observation and non-zero-variance criteria. Remaining missing values were median-imputed only for network estimation, and the imputation count was recorded.

The primary correlation was biweight midcorrelation (`bicor`), with `maxPOutliers=0.10` and `pearsonFallback="individual"`. Pearson was therefore used only for an individual gene or module eigengene with zero or non-finite median absolute deviation. The fallback was not applied globally and every affected entity was exported in an audit table.

The workflow selected a soft-thresholding power from the scale-free-topology diagnostics, calculated signed adjacency and signed topological overlap, performed hierarchical clustering and dynamic tree cutting, merged similar modules, and calculated module eigengenes, module membership (`kME`) and intramodular connectivity.

The cancer-specific WGCNA pair subscore combined:

| WGCNA subcomponent | Weight |
|---|---:|
| TOM divergence | 0.40 |
| Module divergence | 0.30 |
| kME divergence | 0.30 |

### Regulatory-network analysis

Regulatory evidence combined DoRothEA/OmniPath TF-target associations, cancer-specific TF-expression profile divergence and promoter-motif divergence [9,10]. The sublayer weights were:

| Regulatory sublayer | Weight |
|---|---:|
| TF-association divergence | 0.40 |
| TF-expression-profile divergence | 0.35 |
| Promoter-motif divergence | 0.25 |

Canonical promoter windows were obtained from Ensembl. JASPAR vertebrate motifs were scanned with FIMO [11-13]. Motif predictions were labelled as sequence-based cis-regulatory support and were not described as direct binding, occupancy or causal regulation.

### STRING functional-network evidence

STRING functional-network data were acquired with identifier mapping and resumable per-gene caching [10]. Evidence channels were retained where available, including experiments, curated databases, coexpression, neighborhood, cooccurrence, fusion and text mining. The STRING combined score was not interpreted as direct experimental interaction.

### Localization evidence

Subcellular localization was derived from Human Protein Atlas source tables [14]. The localization component considered shared and distinct compartments, source reliability and compatibility with the proposed directional relationship. The presence of any annotation did not automatically produce maximum support.

### Biochemical and structural evidence

Biochemical evidence retained molecular function, reaction, substrate, product, cofactor, domain and residue annotations when available. Structural evidence incorporated experimental structures and AlphaFold models, together with residue-level provenance from sources such as UniProt, M-CSA and PDBe [15-18]. AlphaFold confidence was not treated as experimental validation, and absence of an experimental structure was not treated as evidence of structural absence.

### Pharmacology and tractability

Pharmacology and tractability evidence was assembled from sources including Open Targets, ChEMBL, DGIdb, MyChem, Pharos/TCRD and CIViC, with optional PRISM, GDSC and CTRP analyses [19-21]. Compound identity, mechanism, assay, model, response, development status, experimental tractability and clinical evidence were maintained as separate fields. Tractability was not equated with clinical efficacy.

### Coverage-aware score

For candidate \(i\), context \(c\) and eligible domain set \(E_{ic}\), let \(w_d\) denote the domain weight, \(x_{icd}\in[0,1]\) the normalized observed value and \(m_{icd}\) an indicator of observed usable evidence.

The observed eligible-domain score was:

$$
S_{ic}^{observed}=\frac{\sum_{d\in E_{ic}}w_dm_{icd}x_{icd}}{\sum_{d\in E_{ic}}w_dm_{icd}}
$$

Coverage was:

$$
C_{ic}=\frac{\sum_{d\in E_{ic}}w_dm_{icd}}{\sum_{d\in E_{ic}}w_d}
$$

The final adjusted score was:

$$
S_{ic}^{adjusted}=S_{ic}^{observed}C_{ic}=\frac{\sum_{d\in E_{ic}}w_dm_{icd}x_{icd}}{\sum_{d\in E_{ic}}w_d}
$$

The final parent-domain weights were:

| RSES-Onco domain | Weight |
|---|---:|
| Tumor event | 0.16 |
| Dependency | 0.22 |
| Selectivity | 0.14 |
| Expression compensation | 0.08 |
| Functional relation | 0.06 |
| Functional microniche | 0.16 |
| Validation and tractability | 0.18 |

A real observed zero was retained as zero. Missing components were omitted from the observed denominator and penalized through coverage. Non-eligible components were excluded from both the observed and eligible denominators.

### Source overlap and evidence roles

Evidence records sharing the same original publication, dataset or stable evidence identifier were grouped into one evidence unit. Multiple aggregators could confirm or annotate an observation but could not receive independent full weights for the same underlying event.

Evidence was classified as direct score evidence, prioritization evidence, independent validation evidence, overlapping/confirmatory evidence or interpretative evidence. This distinction was used to avoid double counting and circular validation.

### Multiple testing and evidence categories

When multiple hypotheses formed a statistical test family, false discovery rate was controlled using the Benjamini-Hochberg method [22]. Nominal significance, FDR support and external validation were reported separately. A high integrated score was not interpreted as statistical significance, and statistical significance was not interpreted as clinical efficacy.

Candidates were classified as candidate-universe entries, computational hypotheses, prioritized hypotheses, hypotheses with microniche support, hypotheses with conditional-dependency support, nominally significant results, FDR-supported results, externally validated results, tractable candidates or candidates with clinical evidence.

### Robustness analyses

The supporting workflow included raw versus adjusted score comparison, leave-one-domain-out analysis, individual component ablation, controlled weight perturbation, rank-correlation analysis, top-candidate stability, mechanistic-class and cancer-specific stability, and positive/negative control comparisons when supported by the available data.

### Publication assets and quality control

All figures were generated programmatically. Each registered figure received an exact source-data TSV, generator script, input paths, reproduction command, PNG/PDF/SVG exports, layout-audit JSON and manifest entry. Tables received source and generator provenance through the table manifest and the human-readable asset reproduction registry.

The publication contract comprised 8 main figures, 69 supplementary figures, 231 image exports, 4 main tables and 44 supplementary tables. The documents were generated as editable DOCX files, rendered to PDF and page images, and subjected to automated and manual inspection. Supplementary Figures S68 and S69 were required to appear on separate pages.

## References

1. Kaelin WG Jr. The concept of synthetic lethality in the context of anticancer therapy. *Nature Reviews Cancer*. 2005;5:689-698.
2. O'Neil NJ, Bailey ML, Hieter P. Synthetic lethality and cancer. *Nature Reviews Genetics*. 2017;18:613-623.
3. Meyers RM, Bryan JG, McFarland JM, et al. Computational correction of copy-number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. *Nature Genetics*. 2017;49:1779-1784.
4. Behan FM, Iorio F, Picco G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. *Nature*. 2019;568:511-516.
5. Dempster JM, Rossen J, Kazachkova M, et al. Extracting biological insights from the Project Achilles genome-scale CRISPR screens in cancer cell lines. *bioRxiv*. 2019:720243.
6. Tomczak K, Czerwińska P, Wiznerowicz M. The Cancer Genome Atlas: an immeasurable source of knowledge. *Contemporary Oncology*. 2015;19:A68-A77.
7. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*. 2008;9:559.
8. Langfelder P, Horvath S. Fast R functions for robust correlations and hierarchical clustering. *Journal of Statistical Software*. 2012;46(11):1-17.
9. Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. Benchmark and integration of resources for the estimation of human transcription factor activities. *Genome Research*. 2019;29:1363-1375.
10. Türei D, Valdeolivas A, Gul L, et al. Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. *Molecular Systems Biology*. 2021;17:e9923.
11. Szklarczyk D, Gable AL, Nastou KC, et al. The STRING database in 2021. *Nucleic Acids Research*. 2021;49:D605-D612.
12. Howe KL, Achuthan P, Allen J, et al. Ensembl 2021. *Nucleic Acids Research*. 2021;49:D884-D891.
13. Fornes O, Castro-Mondragon JA, Khan A, et al. JASPAR 2020. *Nucleic Acids Research*. 2020;48:D87-D92.
14. Grant CE, Bailey TL, Noble WS. FIMO: scanning for occurrences of a given motif. *Bioinformatics*. 2011;27:1017-1018.
15. Uhlén M, Fagerberg L, Hallström BM, et al. Tissue-based map of the human proteome. *Science*. 2015;347:1260419.
16. The UniProt Consortium. UniProt: the Universal Protein Knowledgebase. *Nucleic Acids Research*. 2023;51:D523-D531.
17. Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. *Nature*. 2021;596:583-589.
18. Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database. *Nucleic Acids Research*. 2022;50:D439-D444.
19. Ribeiro AJM, Holliday GL, Furnham N, et al. Mechanism and Catalytic Site Atlas. *Nucleic Acids Research*. 2018;46:D618-D623.
20. Ochoa D, Hercules A, Carmona M, et al. Open Targets Platform. *Nucleic Acids Research*. 2021;49:D1302-D1310.
21. Gaulton A, Hersey A, Nowotka M, et al. The ChEMBL database in 2017. *Nucleic Acids Research*. 2017;45:D945-D954.
22. Benjamini Y, Hochberg Y. Controlling the false discovery rate. *Journal of the Royal Statistical Society Series B*. 1995;57:289-300.

## TCGA/GDC DNA methylation layer (v0.11.1)

RSES-Onco now acquires gene-associated CpG methylation beta values from the NCI Genomic Data Commons through the UCSC Xena GDC hub. Repbase is not used because it is a reference library of repetitive DNA sequences rather than a sample-level methylation resource.

```bash
python -u scripts/acquire_tcga_nise_methylation.py \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --output-dir data/processed/epigenetics/methylation
```

Methylation is integrated inside the existing regulatory-network microniche weight. It does not receive a new independent top-level RSES-Onco domain. Regulatory subweights are DoRothEA regulator divergence 0.32, TF-expression-profile divergence 0.28, JASPAR motif divergence 0.20 and TCGA/GDC methylation context 0.20. A technical source failure makes methylation non-eligible and preserves the original three-component regulatory score. Available source data with missing pair-level probes or samples reduce internal regulatory coverage.

New outputs are Figures S70-S72 and Supplementary Tables S45-S47.
