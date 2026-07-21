# Third-party input directory

Large or release-specific third-party matrices are intentionally not redistributed in this repository.

- `SwissProt_NISE.tsv.gz`: run `python scripts/download_human_nise.py` to obtain the 2026 Swiss-Prot NISE source and filter Homo sapiens.
- `depmap/`: place the current official public DepMap files here and validate with `python scripts/download_depmap.py`.
- `gdc/`: run `python scripts/download_gdc.py --manifest-only` and then the download command for public GDC files.

The bundled `data/demo/` matrices are synthetic software-verification fixtures and must not be interpreted as patient or cancer-cell-line measurements.
