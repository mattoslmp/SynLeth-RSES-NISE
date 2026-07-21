# Prespecified analysis plan

1. **Relation gate.** Keep reaction/evolutionary validity separate from cancer dependency evidence. Strict NISE, homologous paralog, pathway backup and collateral-deletion mechanisms are distinct labels.
2. **Tumor event.** Calculate project-specific loss prevalence in TCGA-COAD, TCGA-READ, TCGA-STAD, TCGA-LUAD and TCGA-LUSC. Homozygous deletion is GISTIC -2; mutation-based loss requires a separate pathogenic loss-of-function annotation.
3. **Cell-line dependency.** Within each DepMap lineage, compare the target gene effect between event-positive and event-negative models. Require at least three models per group for pilot analysis; report medians, delta, one-sided Mann-Whitney P and BH-adjusted Q.
4. **Expression compensation.** Test whether the proposed backup is expressed and whether expression increases in the loss group. Expression is supportive, not proof of flux or dependency.
5. **Tumor-normal safety.** Add GTEx or matched-normal expression only when tissue labels and transformations are harmonized. Missing normal coverage reduces score coverage rather than becoming zero.
6. **RSES-Onco.** Calculate an observed-domain weighted mean, evidence coverage and coverage-adjusted score. Report all three together.
7. **Validation.** Prioritize isogenic knockout/knockdown, rescue, small-molecule perturbation, orthogonal viability assays, and in-vivo validation. Computational ranking is hypothesis generation, not clinical evidence.
