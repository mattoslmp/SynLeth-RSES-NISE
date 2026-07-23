#!/usr/bin/env python3
"""Idempotently update canonical RSES-Onco documentation for v0.11.1."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MARKER = "<!-- BEGIN V0.11.1 METHYLATION ADDENDUM -->"
END_MARKER = "<!-- END V0.11.1 METHYLATION ADDENDUM -->"


def read(path: str) -> str:
  return (ROOT / path).read_text(encoding="utf-8")


def write(path: str, text: str) -> None:
  (ROOT / path).write_text(text.rstrip() + "\n", encoding="utf-8")


def append_once(path: str, heading: str, body: str) -> None:
  text = read(path)
  if MARKER in text:
    return
  block = (
    f"\n\n{MARKER}\n\n## {heading}\n\n"
    + body.strip()
    + f"\n\n{END_MARKER}\n"
  )
  write(path, text + block)


def patch_protocol() -> None:
  path = "docs/END_TO_END_ARTICLE_PROTOCOL.md"
  text = read(path)
  replacements = {
    "# RSES-Onco v0.11.0:": "# RSES-Onco v0.11.1:",
    "Repository and publication framework: 0.11.0": (
      "Repository and publication framework: 0.11.1"
    ),
    "Scientific score: RSES-Onco-expanded-v0.10.9": (
      "Scientific score: RSES-Onco-expanded-v0.11.1 (integrated ranking)"
    ),
    "WGCNA/regulatory semantics: eligibility-aware-wgcna-regulatory-v3": (
      "Expression/regulatory semantics: "
      "eligibility-aware-wgcna-regulatory-methylation-v4"
    ),
    "44 supplementary tables": "47 supplementary tables",
    "48 registered tables": "51 registered tables",
    "pre_v0110_": "pre_v0111_",
    "SESSION=\"rses_v0110_": "SESSION=\"rses_v0111_",
    "run_rses_v0110_": "run_rses_v0111_",
    "last_rses_v0110_": "last_rses_v0111_",
  }
  for old, new in replacements.items():
    text = text.replace(old, new)
  acquisition_link = (
    "- [`METHYLATION_INTEGRATION_V0111.md`]"
    "(METHYLATION_INTEGRATION_V0111.md)"
  )
  if acquisition_link not in text:
    anchor = (
      "- [`DATA_ACQUISITION_AND_REPRODUCTION_V0110.md`]"
      "(DATA_ACQUISITION_AND_REPRODUCTION_V0110.md)"
    )
    text = text.replace(anchor, anchor + "\n" + acquisition_link)
  if "export METHYLATION_MODE=download" not in text:
    text = text.replace(
      "export PUBLICATION_STAGE=assets-only\n",
      "export PUBLICATION_STAGE=assets-only\n"
      "export METHYLATION_MODE=download\n",
    )
  write(path, text)
  append_once(
    path,
    "GDC promoter methylation extension",
    """
The v0.11.1 integrated ranking may include open GDC/TCGA promoter methylation beta-values. Methylation is not a new global score domain. It is integrated inside the existing expression-compensation domain with internal weights of 0.70 for event-stratified expression compensation and 0.30 for directional promoter methylation context. The DepMap-only ranking remains free of TCGA methylation.

Acquire, aggregate and integrate methylation with:

```bash
METHYLATION_MODE=download \\
MPLBACKEND=Agg \\
STRICT_LAYOUT=1 \\
bash scripts/resume_wgcna_regulatory_pipeline.sh resume-regulatory
```

Validate the integrated ranking with:

```bash
python -u scripts/validate_methylation_evidence.py \\
  --evidence data/processed/methylation/pair_promoter_methylation_evidence.tsv \\
  --ranking results/expanded_26Q1/full/expanded_rses_onco.tsv
```

Detailed source UUIDs, promoter rules, formulas, missingness semantics and stage-specific commands are provided in [`METHYLATION_INTEGRATION_V0111.md`](METHYLATION_INTEGRATION_V0111.md).
""",
  )


def patch_acquisition() -> None:
  path = "docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md"
  append_once(
    path,
    "GDC/TCGA promoter methylation acquisition",
    """
RSES-Onco v0.11.1 supports open GDC `Methylation Beta Value` files for primary tumors in TCGA-COAD, TCGA-READ, TCGA-STAD, TCGA-LUAD and TCGA-LUSC. The downloader retains file UUID, file name, size, MD5, platform, project, case and sample provenance and retrieves official GDC GENCODE-v36 HM27, HM450 and EPIC probe annotations.

Create the manifest:

```bash
python -u scripts/download_gdc_methylation.py \\
  --stage manifest \\
  --output-dir data/raw/methylation
```

Download and validate open files:

```bash
python -u scripts/download_gdc_methylation.py \\
  --stage all \\
  --output-dir data/raw/methylation \\
  --workers 4 \\
  --retries 4
```

Aggregate promoter beta-values:

```bash
python -u scripts/aggregate_gdc_methylation.py \\
  --manifest data/raw/methylation/gdc_methylation_manifest.tsv \\
  --annotation-dir data/raw/methylation/annotations \\
  --output data/processed/methylation/gdc_promoter_methylation_gene_sample.tsv \\
  --gene-summary data/processed/methylation/gdc_promoter_methylation_gene_summary.tsv \\
  --status-output data/processed/methylation/gdc_promoter_methylation_aggregation_status.tsv
```

Build candidate-pair evidence:

```bash
python -u scripts/build_methylation_pair_evidence.py \\
  --candidates data/processed/expanded_candidate_universe.tsv \\
  --gene-sample data/processed/methylation/gdc_promoter_methylation_gene_sample.tsv \\
  --output data/processed/methylation/pair_promoter_methylation_evidence.tsv \\
  --min-samples 10
```

Beta-values are epigenetic context and are not direct proof of transcriptional silencing. Missing methylation remains missing and is not converted to zero. The complete methodology is documented in [`METHYLATION_INTEGRATION_V0111.md`](METHYLATION_INTEGRATION_V0111.md).
""",
  )


def patch_supplement() -> None:
  addendum = read(
    "supplementary/Supplementary_Methylation_Methods_RSES_Onco_v0111.md"
  )
  append_once(
    "supplementary/Supplementary_Methods_RSES_Onco_v0110.md",
    "Promoter methylation integration in RSES-Onco v0.11.1",
    addendum,
  )


def patch_manuscript() -> None:
  addendum = read(
    "manuscript/RSES_Onco_methylation_methods_addendum_v0111.md"
  )
  append_once(
    "manuscript/RSES_Onco_intro_methods_draft_v0110.md",
    "Promoter methylation analysis",
    addendum,
  )


def main() -> None:
  patch_protocol()
  patch_acquisition()
  patch_supplement()
  patch_manuscript()
  print("Canonical v0.11.1 documentation sources updated.")


if __name__ == "__main__":
  main()
