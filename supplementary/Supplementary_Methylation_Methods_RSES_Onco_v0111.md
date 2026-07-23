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
