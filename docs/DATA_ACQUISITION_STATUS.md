# Data acquisition status

## Included directly

- Curated human bona fide analogue catalogue from Piergiorge et al. (2017): 70 proteins, 15 EC activities.
- Complete cross-cluster pair expansion from that catalogue: 101 human analogue pairs.
- Literature-anchored benchmark: 25 validated vulnerabilities or explicitly labelled discovery hypotheses.
- Synthetic verification fixtures for software testing only.

## Acquired by scripts at the point of use

- Current Swiss-Prot NISE file from Zenodo record 18008936.
- DepMap Public release files (`CRISPRGeneEffect.csv`, `OmicsCNGene.csv`, `Model.csv`, expression matrix).
- Public GDC/TCGA copy-number files for COAD, READ, STAD, LUAD and LUSC.

Large release-specific third-party matrices are not redistributed because their size, versioning and terms require release-aware acquisition. The execution environment used to assemble this package did not provide direct network access to the raw Zenodo/DepMap/GDC payloads. The repository therefore includes official URLs, acquisition code, input validation and an empirical analysis command, while labelling the bundled ranking as a literature prior rather than a completed release-specific TCGA/DepMap discovery analysis.
