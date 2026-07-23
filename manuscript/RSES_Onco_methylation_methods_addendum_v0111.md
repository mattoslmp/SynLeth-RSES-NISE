## Promoter methylation analysis

To incorporate epigenetic context without creating an independent evidence domain, open primary-tumor DNA methylation beta-value files were queried from the NCI Genomic Data Commons for TCGA-COAD, TCGA-READ, TCGA-STAD, TCGA-LUAD and TCGA-LUSC. The analysis used GDC-harmonized SeSAMe methylation-array outputs and the official GDC GENCODE-v36 probe annotations for HM27, HM450 and EPIC arrays.

Promoter-associated probes were defined using the official annotation as probes located from 2 kb upstream to 500 bp downstream of a transcription start site or classified in promoter/TSS-related groups. Masked probes were removed when the GDC general-mask field was available. Promoter methylation for each gene and sample was represented by the median beta-value across retained promoter probes.

For each directional candidate pair and cancer context, promoter beta-values for the hypothesized lost gene and target/paralog were matched within samples. Analyses required at least 10 paired samples. The workflow calculated median beta-values, median paired delta beta, Spearman correlation, a one-sided paired Wilcoxon test and Benjamini-Hochberg-adjusted P values. The directional promoter-methylation context score was defined as the median promoter beta of the hypothesized lost gene multiplied by one minus the median promoter beta of the target gene.

Promoter methylation was integrated inside the existing expression-compensation domain. The global expression-compensation weight remained 0.08, while its internal subweights were 0.70 for event-stratified expression compensation and 0.30 for promoter methylation context. The internal score was coverage-adjusted, so missing methylation reduced subcoverage but was not interpreted as biological zero. GDC methylation was applied only to the integrated TCGA-DepMap ranking and not to the DepMap-only ranking. Methylation beta-values were treated as epigenetic context and not as direct proof of gene silencing.

### Methylation references

1. National Cancer Institute. Genomic Data Commons DNA methylation array harmonization documentation.
2. Zhou W, Triche TJ Jr, Laird PW, Shen H. SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. *Nucleic Acids Research*. 2018;46:e123.
3. Bibikova M, Barnes B, Tsan C, et al. High density DNA methylation array with single CpG site resolution. *Genomics*. 2011;98:288-295.
4. Moran S, Arribas C, Esteller M. Validation of a DNA methylation microarray for 850,000 CpG sites of the human genome enriched in enhancer sequences. *Epigenomics*. 2016;8:389-399.
