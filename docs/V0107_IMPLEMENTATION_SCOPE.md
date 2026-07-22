# RSES-Onco v0.10.7 implementation scope

This release implements the publication-evidence audit requested for the real-data manuscript workflow.

Implemented in `run_publication_pipeline.sh assets-only`:

- candidate × cancer × domain eligibility and missingness audit;
- explicit evidence states and absence reasons;
- source/version/date/status fields when recorded by acquisition metadata;
- score decomposition and independent formula recomputation;
- raw-versus-adjusted, leave-one-domain-out and controlled weight sensitivity analyses;
- evidence overlap registry with one-unit total-weight cap for repeated representations;
- organized expression, CRISPR, tumor-event, network, structural and pharmacology support tables;
- exact table used by every figure, plus SHA-256 and reproduction command;
- expanded Figure S1 and new data-supported Figures S33–S38;
- corrected Figures 1, 2, 4 and 5;
- strict scientific-integrity validation before final package validation;
- 8 main figures, 38 supplementary figures, 4 main tables and 25 supplementary tables.

Scientific boundaries:

- no missing value is converted to zero;
- technical failure is not negative biological evidence;
- composite events remain composite features and are not converted into single genes;
- STRING combined evidence is not called experimental interaction;
- DoRothEA regulatory association is not called direct promoter binding without a direct source;
- tractability is not clinical efficacy;
- robustness analyses use existing observations and do not create evidence;
- unsupported figures, promoter annotations, phenotypes, compound names or document pages are not fabricated.

The repository does not create Figures S39–S69 merely to satisfy numbering. The requirement that S68 and S69 occupy separate pages applies only after real, registered, data-supported S68/S69 figures and manuscript source documents exist. Final manuscript packaging remains blocked until real-data execution succeeds and every rendered page is manually inspected at 100% zoom.
