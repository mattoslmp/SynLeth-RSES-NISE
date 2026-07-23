# Supplementary Methods: RSES-Onco v0.11.0

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

## S1. Scientific objective

RSES-Onco is a coverage-aware computational framework for prioritizing experimentally testable, cancer-context hypotheses involving non-homologous isofunctional enzymes (NISEs), homologous paralogs, pathway alternatives, collateral-loss relationships and downstream dependencies. The framework integrates heterogeneous evidence without interpreting unavailable evidence as biological absence.

The system does not provide a clinical diagnosis, treatment recommendation or proof of therapeutic efficacy. Its output is a ranked and auditable set of hypotheses for subsequent experimental validation.

## S2. Functional microniche concept

A **functional microniche** is the local and context-dependent evidence environment surrounding a directional candidate pair in a defined cancer setting. It is not limited to physical location. Instead, it combines the subset of molecular and phenotypic conditions in which one gene or protein may specialize, compensate for, or become selectively required after alteration of another.

The microniche may include:

- expression and expression-profile divergence;
- transcriptional compensation or lack of compensation;
- WGCNA module and topological-overlap context;
- regulatory-network divergence;
- functional protein-network evidence;
- genetic phenotype and dependency-profile divergence;
- subcellular localization;
- biochemical and structural evidence.

A pair may therefore be functionally interchangeable in one condition, divergent in another, or non-evaluable when the relevant evidence is missing.

## S3. Candidate universe

The candidate universe combines several explicitly labelled sources:

1. directed NISE hypotheses derived from curated human isofunctional activity groups;
2. directed paralog hypotheses derived from Ensembl Compara;
3. curated benchmark or literature-anchored candidates;
4. collateral-loss and downstream-vulnerability hypotheses;
5. all-target conditional-dependency discovery candidates.

Directionality is retained because `gene A loss → gene B dependency` is not assumed to be equivalent to `gene B loss → gene A dependency`. Composite biomarkers remain composite unless the source supports an atomic single-gene interpretation.

## S4. Cancer contexts

The primary contexts are:

```text
colon: TCGA-COAD, TCGA-READ and compatible DepMap colorectal models
stomach: TCGA-STAD and compatible DepMap gastric models
lung: TCGA-LUAD, TCGA-LUSC and compatible DepMap lung models
```

Release-specific sample and model counts are recorded from the actual source files and are not forced to match historical values.

## S5. Evidence states and eligibility

For each candidate, cancer context and evidence domain, the pipeline distinguishes:

- `observed_evidence`: a numeric value was calculated from traceable evidence;
- `observed_negative_evidence`: an eligible, observed component equals zero;
- `observed_neutral_evidence`: an eligible, observed component equals the operational midpoint;
- `missing_evidence`: no eligible observation was available;
- `insufficient_sample`: the minimum sample requirement was not satisfied;
- `technical_or_source_failure`: acquisition, parsing or mapping failed;
- `not_eligible`: the candidate/context cannot be evaluated by that domain;
- `analysis_not_executable`: required inputs or contrasts were not definable.

Missing, insufficient, failed and non-eligible states are not silently converted to zero.

## S6. Coverage-aware scoring

Let \(E_{ic}\) be the set of domains eligible for candidate \(i\) in cancer context \(c\). Let \(O_{ic} \subseteq E_{ic}\) be the subset with observed numeric evidence. For domain \(d\):

- \(w_d\) is the predefined weight;
- \(x_{icd} \in [0,1]\) is the normalized observed component;
- \(m_{icd}=1\) when the component is observed and usable, otherwise \(m_{icd}=0\).

The observed eligible-domain score is:

$$
S_{ic}^{observed}=
\frac{\sum_{d \in E_{ic}} w_d m_{icd}x_{icd}}
{\sum_{d \in E_{ic}} w_d m_{icd}}
$$

The weighted coverage is:

$$
C_{ic}=
\frac{\sum_{d \in E_{ic}} w_d m_{icd}}
{\sum_{d \in E_{ic}} w_d}
$$

The coverage-adjusted score is:

$$
S_{ic}^{adjusted}=S_{ic}^{observed}C_{ic}
=\frac{\sum_{d \in E_{ic}} w_dm_{icd}x_{icd}}
{\sum_{d \in E_{ic}} w_d}
$$

This formulation omits missing components from the observed-score denominator while penalizing sparse evidence through coverage. Ineligible domains are excluded from both denominators. A real observed value of zero remains zero.

## S7. Functional-microniche weights

The functional-microniche score uses the following parent-domain weights, which sum to 1.00:

| Domain | Weight |
|---|---:|
| Expression context | 0.20 |
| Localization | 0.15 |
| Biochemical/structural evidence | 0.15 |
| Genetic phenotype | 0.20 |
| Interaction network | 0.15 |
| Regulatory network | 0.15 |

The functional-microniche output retains the observed score, eligible coverage, adjusted score, observed-domain count, eligible-domain count, observed weight and eligible weight.

## S8. Final RSES-Onco weights

The expanded cancer score uses the following weights, which sum to 1.00:

| Domain | Weight |
|---|---:|
| Tumor event | 0.16 |
| Dependency | 0.22 |
| Selectivity | 0.14 |
| Expression compensation | 0.08 |
| Functional relation | 0.06 |
| Functional microniche | 0.16 |
| Validation and tractability | 0.18 |

The scientific score version is `RSES-Onco-expanded-v0.10.9`.

## S9. Candidate normalization and identifier mapping

Gene symbols are canonicalized before integration. Source-specific identifiers are retained in provenance columns whenever possible. Mapping failure is recorded explicitly and does not imply absent biology. Pair identifiers preserve direction, source class and, when available, activity or relationship identifiers.

Joins are audited for row loss, unexpected duplication and one-to-many expansion. Composite targets and biomarkers are expanded only when the original semantics permit traceable atomic representation.

## S10. DepMap dependency analysis

DepMap gene-effect data are aligned to the model table by `ModelID`. Cancer-compatible models are selected using lineage and disease annotations. Gene-effect values are interpreted according to the DepMap release conventions; more negative values generally represent stronger dependency.

For each candidate/context, the workflow retains when available:

- target-gene effect in event-positive models;
- target-gene effect in event-intact models;
- sample sizes for both groups;
- effect-size estimate;
- test statistic, nominal P value and FDR;
- general dependency and event-conditioned dependency status;
- mapping and eligibility state.

Conditional dependency is distinguished from general essentiality.

## S11. Selectivity analysis

Selectivity quantifies differential vulnerability associated with a defined event or state. It is not inferred solely from a strongly negative overall gene-effect value. The event-positive and event-intact distributions are compared using the analysis implemented by the scoring pipeline, and effect direction, sample size, nominal significance and multiple-testing correction are retained.

High selectivity component values do not constitute clinical efficacy.

## S12. Expression and compensation

Expression data are derived from the DepMap protein-coding log-transformed expression matrix and aligned to compatible models. The source matrix is treated according to its documented release transformation rather than transformed a second time without justification.

The supporting tables may include:

- expression of the lost/source gene and target gene;
- cancer-specific model values;
- event-positive and event-intact summaries;
- median or mean group difference, as defined by the script;
- pairwise Spearman correlation;
- sample size;
- compensation status;
- missingness or insufficiency reason.

Simple correlation is not treated as automatic evidence of compensation. Compensation requires an operationally defined event-stratified expression pattern consistent with increased or preserved target expression after the relevant event.

## S13. Pairwise expression-profile divergence

For cancer-compatible models with both genes observed, the implementation calculates:

- Spearman correlation \(
ho\);
- median absolute expression difference;
- correlation divergence \((1-\rho)/2\);
- saturated expression separation, clipped to \([0,1]\);
- mean of the available divergence subcomponents.

At least three paired models are required; otherwise the pairwise metric remains missing.

## S14. Genetic-phenotype divergence

Genetic-phenotype divergence compares cancer-specific DepMap gene-effect profiles. The implementation may combine:

- Spearman-profile divergence;
- median absolute profile separation;
- divergence of dependency-positive model sets using one minus Jaccard similarity.

The dependency-set threshold and saturation parameters are recorded by the generating code and exported supporting tables.

## S15. WGCNA expression-network analysis

Cancer-specific WGCNA networks are built separately for colon, stomach and lung models. Candidate genes are combined with cancer-specific highly variable genes. Genes must satisfy the script-defined observation and non-zero-variance criteria. Remaining missing values may be median-imputed for network estimation only; the number of imputed values is recorded.

The correlation policy is:

```text
primary correlation: bicor
maxPOutliers: 0.10
pearsonFallback: individual
network type: signed
TOM type: signed
```

Pearson is used only for an individual gene or module eigengene with zero or non-finite MAD. It is not a global replacement for bicor. Every fallback-eligible entity is exported.

The soft-threshold power is selected using the scale-free-topology diagnostics implemented in `scripts/run_wgcna_expression_network.R`. The workflow calculates signed adjacency, signed topological overlap, average-linkage clustering, dynamic tree cutting, module merging, module eigengenes, module membership (`kME`) and intramodular connectivity.

### WGCNA pair subscore

The WGCNA pair subscore uses:

| Subcomponent | Weight |
|---|---:|
| TOM divergence | 0.40 |
| Module divergence | 0.30 |
| kME divergence | 0.30 |

Missing subcomponents reduce internal coverage and are not zero-imputed.

## S16. Regulatory-network analysis

The regulatory domain integrates three sublayers:

| Sublayer | Weight |
|---|---:|
| TF-association divergence | 0.40 |
| TF-expression-profile divergence | 0.35 |
| Promoter-motif divergence | 0.25 |

DoRothEA/OmniPath associations provide regulator-target support. Cancer-specific TF and target expression provide expression-consistency context. Ensembl promoter sequences are scanned against JASPAR motifs with FIMO.

A motif hit is sequence-based cis-regulatory support. It is not direct TF binding, ChIP occupancy or causal regulation. The pipeline labels promoter evidence as motif prediction and keeps direct-binding claims false unless a separate direct-binding source is introduced.

## S17. STRING interaction-network evidence

STRING evidence is preserved with channel-level interpretation where available:

- experiments;
- curated databases;
- coexpression;
- neighborhood;
- cooccurrence;
- gene fusion;
- text mining;
- combined score.

The combined score is not described as direct experimental binding. Node and edge tables are exported for network figures.

## S18. Localization evidence

Localization annotations are derived from the Human Protein Atlas or other explicitly recorded source tables. The component considers shared and distinct compartments, evidence reliability and compatibility with the proposed directional relationship. The mere presence of a localization annotation does not automatically produce maximum support.

## S19. Biochemical and structural evidence

Biochemical and structural support may include:

- activity or reaction annotation;
- substrate, product and cofactor information;
- conserved domains;
- reviewed catalytic or binding residues;
- experimental structures;
- AlphaFold models;
- structural coverage and confidence;
- pairwise structural evidence.

AlphaFold confidence is not experimental validation. Absence of an experimental structure is missing evidence, not evidence of absent structural relatedness. Residue provenance must remain traceable to UniProt, M-CSA, PDBe or another declared source.

## S20. Tumor-event evidence

Tumor-event evidence is derived from validated TCGA/GDC matrices and candidate-specific event definitions. Supporting tables retain cancer, event type, sample count, frequency, genomic context and source. Homozygous-deletion frequencies are bounded to the biologically valid interval. Plotting jitter, when used, is graphical only and never written as the true source value.

## S21. Pharmacology and tractability

Pharmacology and tractability evidence may integrate Open Targets, ChEMBL, DGIdb, MyChem, Pharos/TCRD, CIViC and optional PRISM, GDSC or CTRP data. The following are kept distinct:

- known compound or probe;
- mechanism of action;
- assay or model;
- response measurement;
- development status;
- experimental tractability;
- preclinical evidence;
- clinical evidence.

Tractability is not presented as efficacy.

## S22. Source overlap and independence

Evidence representations sharing an original publication, dataset or stable evidence identifier are grouped as one evidence unit. A source aggregator may corroborate or annotate an underlying observation but does not receive a second full independent weight for the same event.

Each source record should be classified as one of:

- independent score evidence;
- overlapping/confirmatory evidence;
- prioritization evidence;
- independent validation evidence;
- interpretative evidence only.

This prevents circular validation and double counting.

## S23. Multiple testing

When multiple hypotheses form a test family, false-discovery-rate control is applied using the Benjamini-Hochberg procedure. The outputs distinguish:

- no statistical support;
- nominal significance;
- FDR-supported result;
- external validation.

An integrated score is not itself a P value, and a high score does not automatically imply statistical significance.

## S24. Evidence-category terminology

The publication uses explicit categories:

1. candidate in the universe;
2. computational hypothesis;
3. prioritized hypothesis;
4. hypothesis with microniche support;
5. hypothesis with conditional-dependency support;
6. nominally significant result;
7. FDR-supported result;
8. result with external validation;
9. tractable candidate;
10. candidate with clinical evidence.

The terms “discovery”, “validated”, “selective”, “synthetic lethal” and “therapeutic target” are restricted to cases satisfying the corresponding evidence criteria.

## S25. Robustness analyses

The publication workflow includes, when supported by the real data:

- raw versus coverage-adjusted score comparison;
- leave-one-domain-out analysis;
- individual component ablation;
- controlled weight perturbation;
- rank-correlation analysis;
- top-candidate stability;
- cancer-specific and mechanistic-class stability;
- sensitivity to missing-data treatment;
- positive and negative controls;
- null or permutation analysis when statistically appropriate.

Robustness analyses are supplementary and do not replace the primary scoring results.

## S26. Figure source data and reproducibility

Every registered figure must have:

- an exact source-data TSV;
- the generator script;
- input paths;
- the exact reproduction command;
- PNG, PDF and SVG exports;
- layout-audit JSON;
- figure-manifest entry.

Every registered table must have a manifest entry with source path, script, row count, column count and status. Table S44 is the human-readable asset reproduction registry.

## S27. Document generation and manual inspection

The document pipeline creates editable DOCX files, rendered PDFs and page PNGs. Automated validation checks file existence, size, page rendering and figure-page mapping. Supplementary Figures S68 and S69 must be on separate pages.

Automated checks do not replace inspection of each figure and page at 100% zoom. The manual checklist must not be auto-filled.

## S28. Practical applications

RSES-Onco supports prioritization for:

- CRISPR knockout or knockdown experiments;
- event-stratified dependency validation;
- rescue experiments;
- expression perturbation;
- mechanistic comparison of NISE and paralog backups;
- collateral-vulnerability studies;
- structural and biochemical follow-up;
- exploratory compound or probe testing.

The framework organizes hypotheses and evidence; it does not replace experimental validation.

## S29. Software and principal scripts

The principal execution and publication entry points are:

```text
scripts/run_expanded_pipeline.sh
scripts/resume_wgcna_regulatory_pipeline.sh
scripts/run_structural_pipeline.sh
scripts/run_publication_pipeline.sh
scripts/verify_complete_article_run.sh
```

The complete command sequence is documented in `docs/END_TO_END_ARTICLE_PROTOCOL.md`, and source acquisition is documented in `docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md`.

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
11. Szklarczyk D, Gable AL, Nastou KC, et al. The STRING database in 2021: customizable protein-protein networks and functional characterization of user-uploaded gene/measurement sets. *Nucleic Acids Research*. 2021;49:D605-D612.
12. Howe KL, Achuthan P, Allen J, et al. Ensembl 2021. *Nucleic Acids Research*. 2021;49:D884-D891.
13. Fornes O, Castro-Mondragon JA, Khan A, et al. JASPAR 2020: update of the open-access database of transcription factor binding profiles. *Nucleic Acids Research*. 2020;48:D87-D92.
14. Grant CE, Bailey TL, Noble WS. FIMO: scanning for occurrences of a given motif. *Bioinformatics*. 2011;27:1017-1018.
15. Uhlén M, Fagerberg L, Hallström BM, et al. Tissue-based map of the human proteome. *Science*. 2015;347:1260419.
16. The UniProt Consortium. UniProt: the Universal Protein Knowledgebase. *Nucleic Acids Research*. 2023;51:D523-D531.
17. Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. *Nature*. 2021;596:583-589.
18. Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. *Nucleic Acids Research*. 2022;50:D439-D444.
19. Ribeiro AJM, Holliday GL, Furnham N, et al. Mechanism and Catalytic Site Atlas (M-CSA): a database of enzyme reaction mechanisms and active sites. *Nucleic Acids Research*. 2018;46:D618-D623.
20. Burley SK, Bhikadiya C, Bi C, et al. RCSB Protein Data Bank: powerful new tools for exploring 3D structures of biological macromolecules. *Nucleic Acids Research*. 2021;49:D437-D451.
21. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*. 1995;57:289-300.
22. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science*. 2020;369:1318-1330.
23. Ochoa D, Hercules A, Carmona M, et al. Open Targets Platform: supporting systematic drug-target identification and prioritisation. *Nucleic Acids Research*. 2021;49:D1302-D1310.
24. Gaulton A, Hersey A, Nowotka M, et al. The ChEMBL database in 2017. *Nucleic Acids Research*. 2017;45:D945-D954.
25. Freshour SL, Kiwala S, Cotto KC, et al. Integration of the Drug-Gene Interaction Database (DGIdb 4.0) with open crowdsource efforts. *Nucleic Acids Research*. 2021;49:D1144-D1151.


<!-- BEGIN V0.11.1 METHYLATION ADDENDUM -->

## Promoter methylation integration in RSES-Onco v0.11.1

# Supplementary Methods Addendum: promoter methylation integration in RSES-Onco v0.11.1

**Author:** Leandro de Mattos Pereira  
**Affiliation:** Databiomics, Laboratório de Bioinformática e Ciências de Dados, WBPereira  
**Address:** Av. Coronel José Bastos, Itaperuna, RJ, Brazil

## Data source and cohort definition

Open DNA methylation beta-value files were queried from the NCI Genomic Data Commons for primary tumors in TCGA-COAD, TCGA-READ, TCGA-STAD, TCGA-LUAD and TCGA-LUSC. COAD and READ were combined as the colorectal context, STAD represented the gastric context, and LUAD/LUSC were combined as the lung context. File identifiers, names, sizes, MD5 checksums, platform, workflow, project, case and sample identifiers were retained in a machine-readable manifest.

The workflow targets GDC `Methylation Beta Value` files produced by the harmonized SeSAMe methylation-array pipeline. Official GDC GENCODE-v36 probe annotation manifests were acquired for the HumanMethylation27, HumanMethylation450 and EPIC platforms.

## Download and validation

Files were downloaded from the GDC data endpoint with temporary `.part` files, retry logic, file-size validation and MD5 verification. A completed file was reused only when its size and MD5 matched the GDC manifest. Source metadata and per-file status were written to `data/raw/methylation/`.

## Probe-to-promoter mapping

Probes were retained as promoter-proximal if the official annotation placed them between 2000 bp upstream and 500 bp downstream of a transcription start site, or if the platform annotation classified the probe as TSS1500, TSS200, 5' UTR, first exon or promoter-associated. Probes marked by the official general mask were removed when `MASK_general` was available.

For each sample $s$ and gene $g$, promoter methylation was summarized as the median beta-value across retained promoter probes:

$$
\beta_{g,s}=\operatorname{median}_{p\in P_g}(\beta_{p,s})
$$

where $P_g$ is the set of retained promoter probes mapped to gene $g$.

## Pair-level methylation metrics

For each directional candidate pair and cancer context, samples with promoter beta-values for both the hypothesized lost gene and the target/paralog were retained. Analyses required at least 10 paired samples by default. The following metrics were exported:

- paired sample count;
- median promoter beta for the lost gene;
- median promoter beta for the target/paralog;
- median paired delta beta;
- Spearman correlation of promoter beta profiles;
- one-sided paired Wilcoxon test for lost-gene beta greater than target-gene beta;
- global Benjamini-Hochberg FDR;
- within-cancer Benjamini-Hochberg FDR;
- availability and absence reason.

The directional promoter-methylation context score was:

$$
M_{ic}=\widetilde{\beta}_{lost,ic}\left(1-\widetilde{\beta}_{target,ic}\right)
$$

where the median beta-values were bounded to [0,1]. This score represents an epigenetic configuration compatible with greater promoter methylation of the hypothesized lost gene and lower promoter methylation of the target. It is not interpreted as direct evidence of transcriptional silencing.

## Integration with expression compensation

The global RSES-Onco domain weights were not changed. Methylation was integrated as an internal subcomponent of the existing expression-compensation domain, which retains a global weight of 0.08.

Internal subweights were:

- event-stratified expression compensation: 0.70;
- promoter methylation context: 0.30.

For availability indicators $m_X$ and $m_M$, expression-compensation evidence $X_{ic}$ and methylation evidence $M_{ic}$, the observed internal score was:

$$
E^{obs}_{ic}=\frac{0.70X_{ic}m_X+0.30M_{ic}m_M}{0.70m_X+0.30m_M}
$$

Internal coverage was:

$$
C^{expr+meth}_{ic}=0.70m_X+0.30m_M
$$

and the value entering the global score was:

$$
E^{adj}_{ic}=E^{obs}_{ic}C^{expr+meth}_{ic}
$$

If methylation was missing for an otherwise eligible candidate, the expression value remained the observed internal score but internal coverage was reduced to 0.70. Missing methylation was not converted to zero. If the entire GDC methylation layer was unavailable in `auto` mode, the prior score semantics were retained rather than retroactively penalizing every candidate.

## Eligibility and missingness

Composite events that could not be represented as a simple lost-gene/target-gene pair were marked non-eligible for pair-level promoter methylation analysis. Other explicit states included missing genes or cancer context, unsupported platform annotation, insufficient paired samples, absent promoter probes, technical file failure and observed evidence.

## Separation of analysis scopes

GDC methylation was applied only to the integrated TCGA-DepMap ranking. The DepMap-only ranking remained free of GDC methylation evidence. This preserves the intended interpretation of the two score products.

## Publication and reproducibility outputs

The methylation extension generated:

- `data/raw/methylation/gdc_methylation_manifest.tsv`;
- `data/raw/methylation/gdc_methylation_download_status.tsv`;
- `data/processed/methylation/gdc_promoter_methylation_gene_sample.tsv`;
- `data/processed/methylation/gdc_promoter_methylation_gene_summary.tsv`;
- `data/processed/methylation/gdc_promoter_methylation_aggregation_status.tsv`;
- `data/processed/methylation/pair_promoter_methylation_evidence.tsv`;
- Supplementary Figure S69;
- Supplementary Tables S45-S47.

## References

1. National Cancer Institute. Genomic Data Commons DNA methylation array harmonization documentation.
2. National Cancer Institute. GDC Reference Files: GENCODE-v36 methylation array probe annotation manifests.
3. Zhou W, Triche TJ Jr, Laird PW, Shen H. SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. *Nucleic Acids Research*. 2018;46:e123.
4. Bibikova M, Barnes B, Tsan C, et al. High density DNA methylation array with single CpG site resolution. *Genomics*. 2011;98:288-295.
5. Moran S, Arribas C, Esteller M. Validation of a DNA methylation microarray for 850,000 CpG sites of the human genome enriched in enhancer sequences. *Epigenomics*. 2016;8:389-399.
6. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*. 1995;57:289-300.

<!-- END V0.11.1 METHYLATION ADDENDUM -->
