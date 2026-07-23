#!/usr/bin/env python3
"""Idempotently synchronize canonical tutorials and methods with Circos v0.11.2."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BEGIN = "<!-- BEGIN GENOMIC CIRCOS V0.11.2 -->"
END = "<!-- END GENOMIC CIRCOS V0.11.2 -->"


def read(path: str) -> str:
  return (ROOT / path).read_text(encoding="utf-8")


def write(path: str, text: str) -> None:
  (ROOT / path).write_text(text.rstrip() + "\n", encoding="utf-8")


def replace_block(text: str, heading: str, body: str) -> str:
  block = f"{BEGIN}\n\n## {heading}\n\n{body.strip()}\n\n{END}"
  if BEGIN in text and END in text:
    prefix, remainder = text.split(BEGIN, 1)
    _, suffix = remainder.split(END, 1)
    return prefix.rstrip() + "\n\n" + block + suffix
  return text.rstrip() + "\n\n" + block + "\n"


def patch_protocol() -> None:
  path = "docs/END_TO_END_ARTICLE_PROTOCOL.md"
  text = read(path)
  replacements = {
    "# RSES-Onco v0.11.1:": "# RSES-Onco v0.11.2:",
    "Repository and publication framework: 0.11.1": (
      "Repository and publication framework: 0.11.2"
    ),
    "69 supplementary figures": "70 supplementary figures",
    "77 registered figures": "78 registered figures",
    "231 PNG/PDF/SVG exports": "234 PNG/PDF/SVG exports",
    "44 supplementary tables": "52 supplementary tables",
    "48 registered tables": "56 registered tables",
    "rses_v0111_": "rses_v0112_",
    "run_rses_v0111_": "run_rses_v0112_",
    "last_rses_v0111_": "last_rses_v0112_",
    "verify_complete_article_run_v0111.log": (
      "verify_complete_article_run_v0112.log"
    ),
    "RSES_Onco_v0111_complete_submission_package.zip": (
      "RSES_Onco_v0112_complete_submission_package.zip"
    ),
  }
  for old, new in replacements.items():
    text = text.replace(old, new)
  link = (
    "- [`GENOMIC_CIRCOS_WORKFLOW_V0112.md`]"
    "(GENOMIC_CIRCOS_WORKFLOW_V0112.md)"
  )
  if link not in text:
    anchor = (
      "- [`METHYLATION_DATA_AND_SCORING_V0111.md`]"
      "(METHYLATION_DATA_AND_SCORING_V0111.md)"
    )
    text = text.replace(anchor, anchor + "\n" + link)
  body = r"""
The canonical `all` and `assets-only` workflows generate Supplementary Figure S70 as a GRCh38 genomic Circos containing every coordinate-complete simple NISE and homologous-paralog hypothesis. Red chords denote NISEs and black chords denote homologous paralogs. Panel A contains the complete top-level RSES-Onco score, coverage and seven domains. Panel B contains all functional-microniche, expression-network, WGCNA, regulatory, promoter-motif, methylation and nested-coverage layers.

Before publication assets are generated, the wrapper executes:

```bash
python -u scripts/build_script_documentation.py

python -u scripts/build_genomic_circos_inputs.py \
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --promoters data/raw/regulatory/ensembl_promoters.tsv \
  --expression data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
  --models data/raw/depmap/Model.csv \
  --output-dir data/processed/circos
```

The final publication stage generates and registers S70 and Supplementary Tables S45-S52. Table S50 contains every model-level expression value used for every Circos gene. Table S51 is the complete source-derived catalogue of every Python, Bash and R script/module under `scripts/` and `src/rses_onco/`.

The canonical command remains:

```bash
MPLBACKEND=Agg \
STRICT_LAYOUT=1 \
bash scripts/run_publication_pipeline.sh assets-only
```

Detailed ring definitions, source tables and individual commands are documented in [`GENOMIC_CIRCOS_WORKFLOW_V0112.md`](GENOMIC_CIRCOS_WORKFLOW_V0112.md).
"""
  write(path, replace_block(text, "Genomic Circos and complete code catalogue", body))


def patch_acquisition() -> None:
  path = "docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md"
  text = read(path)
  body = r"""
The genomic Circos requires the completed candidate universe, final cancer-specific ranking, Ensembl canonical-promoter coordinate table, DepMap model metadata and protein-coding expression matrix. The Ensembl table provides chromosome and canonical TSS coordinates; no gene lacking a canonical coordinate is assigned an invented position.

All expression columns for Circos genes are read directly from `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`. The pipeline exports both cancer-by-gene summaries and every model-level `log2(TPM+1)` value. These outputs become Supplementary Tables S49 and S50.

Run the source-data stage with:

```bash
python -u scripts/build_genomic_circos_inputs.py \
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --promoters data/raw/regulatory/ensembl_promoters.tsv \
  --expression data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
  --models data/raw/depmap/Model.csv \
  --output-dir data/processed/circos
```

The generated source-provenance table retains every source path, byte size and SHA-256. Missing coordinates terminate the stage with the unresolved gene names; missing domain values remain NA and are never converted to zero.
"""
  write(path, replace_block(text, "Genomic Circos source acquisition and provenance", body))


def patch_supplement() -> None:
  path = "supplementary/Supplementary_Methods_RSES_Onco_v0110.md"
  text = read(path)
  body = r"""
Supplementary Figure S70 was generated from every coordinate-complete simple-gene candidate classified as a NISE, homologous paralog or both. Canonical positions were derived from the Ensembl canonical-transcript lookup used for promoter acquisition and displayed on a GRCh38 chromosome ideogram. Every gene was represented as a genomic tick; NISE relationships were connected in red and homologous-paralog relationships in black. Chord width and transparency were proportional to the maximum cancer-specific coverage-adjusted RSES-Onco score for the pair.

Panel A included rings for coverage-adjusted RSES-Onco, evidence coverage, tumor event, dependency, selectivity, expression compensation, functional relation, functional microniche and validation/tractability. Panel B included expression context, localization, biochemical/structural evidence, genetic phenotype, interaction network, regulatory network, pairwise expression, WGCNA, DoRothEA, TF-expression profiles, JASPAR/FIMO promoter motifs, promoter methylation and nested coverage values.

For each gene and ring, the plotted value was the maximum observed value across associated pair-by-cancer records. Median, minimum, observed-row count, eligible-row count and missingness status were preserved in Supplementary Table S47. Missing or non-eligible evidence was rendered as a hollow marker and not as numeric zero.

All model-level DepMap `log2(TPM+1)` values used for Circos genes were exported in Supplementary Table S50. Supplementary Tables S45-S52 contain coordinates, links, rings, ring definitions, expression summaries, model-level expression, complete script documentation and source provenance. The exact combined source TSV used to render S70 was registered and copied byte-for-byte to the figure-data directory.
"""
  write(path, replace_block(text, "Genomic Circos methods and source-data transparency", body))


def patch_manuscript() -> None:
  path = "manuscript/RSES_Onco_intro_methods_draft_v0110.md"
  text = read(path)
  body = r"""
A genomic Circos representation was generated to integrate chromosomal position, NISE/paralog relationships and every RSES-Onco evidence layer. Canonical Ensembl/GRCh38 positions were assigned only to coordinate-resolved genes. NISE and homologous-paralog links were displayed separately, while concentric rings summarized top-level score domains, functional-microniche domains, expression-network and regulatory subcomponents, promoter methylation and evidence coverage. Missing evidence was retained as missing and displayed using hollow markers. Complete model-level expression values and exact figure-source tables were exported as supplementary data.
"""
  write(path, replace_block(text, "Genomic Circos visualization", body))


def main() -> None:
  patch_protocol()
  patch_acquisition()
  patch_supplement()
  patch_manuscript()
  print("Canonical documentation synchronized with genomic Circos v0.11.2.")


if __name__ == "__main__":
  main()
