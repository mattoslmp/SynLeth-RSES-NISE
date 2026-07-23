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
The canonical `all` and `assets-only` workflows generate Supplementary Figure S70 as a GRCh38 genomic Circos containing every coordinate-complete simple NISE and homologous-paralog hypothesis. Red chords denote NISEs and black chords denote homologous paralogs. Every simple pair remains represented even when the final score is unavailable.

Figure S70 contains 35 score/domain/internal tracks. Panel A contains observed and coverage-adjusted RSES-Onco, evidence coverage, all seven top-level domains and the four individual validation/tractability terms. Panel B contains all six functional-microniche domains, pairwise expression, WGCNA composite, TOM, module and kME divergence, DoRothEA, TF-expression profiles, JASPAR/FIMO promoter motifs, methylation composite, methylation-profile divergence, conditional target hypomethylation and nested coverage values.

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

python -u scripts/enrich_genomic_circos_internal_layers.py
python -u scripts/complete_genomic_circos_expression_summary.py
python -u scripts/complete_genomic_circos_links.py
```

The final stage generates and registers S70 and Supplementary Tables S45-S52. Table S50 contains every observed model-level expression value plus explicit `NA` sentinel rows for unavailable gene/context combinations. Table S51 is the complete source-derived catalogue of every Python, Bash and R script/module under `scripts/` and `src/rses_onco/`.

```bash
MPLBACKEND=Agg \
STRICT_LAYOUT=1 \
bash scripts/run_publication_pipeline.sh assets-only
```

Detailed ring definitions, source tables and commands are documented in [`GENOMIC_CIRCOS_WORKFLOW_V0112.md`](GENOMIC_CIRCOS_WORKFLOW_V0112.md).
"""
  write(path, replace_block(text, "Genomic Circos and complete code catalogue", body))


def patch_acquisition() -> None:
  path = "docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md"
  text = read(path)
  body = r"""
The genomic Circos requires the completed candidate universe, final cancer-specific ranking, Ensembl canonical-promoter coordinate table, DepMap model metadata, protein-coding expression matrix and WGCNA pair-metrics table. No gene lacking a supported GRCh38 coordinate is assigned an invented position.

All required expression columns are read directly from `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`. The pipeline exports cancer-by-gene summaries and every observed model-level `log2(TPM+1)` value. Gene/context combinations without values receive explicit `NA` sentinel rows with `is_measurement=false`, never numeric zero.

```bash
python -u scripts/build_genomic_circos_inputs.py \
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --promoters data/raw/regulatory/ensembl_promoters.tsv \
  --expression data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
  --models data/raw/depmap/Model.csv \
  --output-dir data/processed/circos
```

The source-provenance table retains every source path, byte size and SHA-256. Missing coordinates terminate the stage with unresolved gene names; missing domain values remain NA. A completion stage verifies that the pair-link table contains exactly one chord for every simple NISE/paralog candidate pair.
"""
  write(path, replace_block(text, "Genomic Circos source acquisition and provenance", body))


def patch_supplement() -> None:
  path = "supplementary/Supplementary_Methods_RSES_Onco_v0110.md"
  text = read(path)
  body = r"""
Supplementary Figure S70 was generated from every coordinate-complete simple-gene candidate classified as a NISE, homologous paralog or both. Canonical positions were derived from the Ensembl canonical-transcript lookup and displayed on a GRCh38 ideogram including chromosomes 1-22, X, Y and MT. Every gene was a genomic tick; NISE relationships were red and homologous-paralog relationships black. Every simple candidate pair was retained exactly once, including score-missing pairs.

The figure contained 35 tracks. Panel A included observed and coverage-adjusted RSES-Onco, evidence coverage, all seven top-level domains and genetic-screen, isogenic, in vivo and clinical-tractability terms. Panel B included all six functional-microniche domains; pairwise expression; WGCNA composite, TOM, module and kME terms; DoRothEA; TF-expression profiles; JASPAR/FIMO motifs; promoter-methylation composite, methylation-profile divergence and conditional target hypomethylation; and nested coverage values.

For each gene and ring, the plotted value was the maximum observed value across associated pair-by-cancer records. Median, minimum, observed-row count, eligible-row count and missingness status were retained in Supplementary Table S47. Missing or non-eligible evidence was a hollow marker, not numeric zero.

All observed model-level DepMap `log2(TPM+1)` values were exported in Table S50. Unavailable gene/context combinations were represented by `NA` sentinel rows with `is_measurement=false`. Tables S45-S52 contain coordinates, all links, all ring values, 35 ring definitions, expression summaries, model-level expression, complete script documentation and source provenance. The exact combined source TSV used for S70 was copied byte-for-byte to the figure-data directory.
"""
  write(path, replace_block(text, "Genomic Circos methods and source-data transparency", body))


def patch_manuscript() -> None:
  path = "manuscript/RSES_Onco_intro_methods_draft_v0110.md"
  text = read(path)
  body = r"""
A genomic Circos representation integrated chromosomal position, every simple NISE/paralog relationship and 35 RSES-Onco score/domain/internal layers. Canonical Ensembl/GRCh38 positions were assigned only to coordinate-resolved genes. NISE and homologous-paralog links were displayed separately, while concentric rings summarized top-level domains, validation terms, functional-microniche domains, pairwise expression, WGCNA internal terms, regulatory subcomponents, promoter-methylation internal terms and evidence coverage. Missing evidence remained missing and was displayed using hollow markers. Complete model-level expression values, explicit unavailable sentinels and exact figure-source tables were exported as supplementary data.
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
