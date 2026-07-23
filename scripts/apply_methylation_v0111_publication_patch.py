#!/usr/bin/env python3
"""Apply publication, table, validation and documentation methylation patches."""
from __future__ import annotations

from pathlib import Path
import re

import yaml

ROOT = Path(__file__).resolve().parents[1]


def read(path: str) -> str:
  return (ROOT / path).read_text(encoding="utf-8")


def write(path: str, text: str) -> None:
  (ROOT / path).write_text(text, encoding="utf-8")


def replace_once(path: str, old: str, new: str) -> None:
  text = read(path)
  if new in text:
    return
  if old not in text:
    raise RuntimeError(f"Patch anchor missing in {path}: {old[:100]!r}")
  write(path, text.replace(old, new, 1))


def append_once(path: str, marker: str, block: str) -> None:
  text = read(path)
  if marker not in text:
    write(path, text.rstrip() + "\n\n" + block.strip() + "\n")


def patch_config() -> None:
  path = ROOT / "config/article_assets.yaml"
  config = yaml.safe_load(path.read_text(encoding="utf-8"))
  figures = config["supplementary_figures"]
  existing = {record["id"] for record in figures}
  additions = [
    {
      "id": "Figure_S70",
      "file": "Figure_S70_methylation_gene_coverage",
      "title": "TCGA/GDC methylation coverage for candidate genes",
      "caption": (
        "Gene-associated CpG methylation coverage and primary-tumor "
        "sample support by cancer context."
      ),
    },
    {
      "id": "Figure_S71",
      "file": "Figure_S71_pair_methylation_context",
      "title": "NISE and paralog pair methylation context",
      "caption": (
        "Coverage-adjusted methylation-context divergence for directed "
        "candidate pairs; methylation is not direct proof of silencing."
      ),
    },
    {
      "id": "Figure_S72",
      "file": "Figure_S72_methylation_regulatory_integration",
      "title": "Integration of methylation into regulatory microniches",
      "caption": (
        "Relationship between the methylation subcomponent and the "
        "integrated regulatory-network component."
      ),
    },
  ]
  figures.extend(item for item in additions if item["id"] not in existing)
  tables = config["supplementary_tables"]
  for name in (
    "Table_S45_tcga_nise_methylation_gene_summary.tsv",
    "Table_S46_tcga_nise_methylation_pair_metrics.tsv",
    "Table_S47_tcga_nise_methylation_source_status.tsv",
  ):
    if name not in tables:
      tables.append(name)
  path.write_text(
    yaml.safe_dump(
      config,
      sort_keys=False,
      allow_unicode=True,
      width=120,
    ),
    encoding="utf-8",
  )


def patch_orchestrator() -> None:
  path = "scripts/make_all_article_figures.py"
  replace_once(
    path,
    '  run([\n    sys.executable, "-u", "scripts/make_extended_supporting_figures.py",\n    "--config", args.config,\n    "--output-root", args.output_root,\n    strict_flag,\n  ])',
    '  run([\n    sys.executable, "-u", "scripts/make_extended_supporting_figures.py",\n    "--config", args.config,\n    "--output-root", args.output_root,\n    strict_flag,\n  ])\n  run([\n    sys.executable, "-u", "scripts/make_methylation_figures.py",\n    "--config", args.config,\n    "--output-root", args.output_root,\n    strict_flag,\n  ])',
  )
  replace_once(
    path,
    '    "extended_supplementary_figure_manifest.tsv",\n  ]',
    '    "extended_supplementary_figure_manifest.tsv",\n    "methylation_supplementary_figure_manifest.tsv",\n  ]',
  )
  replace_once(path, "  expected = 77", "  expected = 80")
  replace_once(
    path,
    '    "extended_supplementary_figure_legends.md",\n  ]',
    '    "extended_supplementary_figure_legends.md",\n    "methylation_supplementary_figure_legends.md",\n  ]',
  )


def patch_tables() -> None:
  path = "scripts/materialize_extended_supplementary_tables.py"
  replace_once(
    path,
    '    "Table_S42_pharmacology_source_coverage.tsv": ROOT / "results/expanded_26Q1/pharmacology/pharmacology_source_coverage.tsv",\n  }',
    '    "Table_S42_pharmacology_source_coverage.tsv": ROOT / "results/expanded_26Q1/pharmacology/pharmacology_source_coverage.tsv",\n    "Table_S45_tcga_nise_methylation_gene_summary.tsv": ROOT / "data/processed/epigenetics/methylation/tcga_nise_methylation_gene_summary.tsv",\n    "Table_S46_tcga_nise_methylation_pair_metrics.tsv": ROOT / "data/processed/epigenetics/methylation/tcga_nise_methylation_pair_metrics.tsv",\n    "Table_S47_tcga_nise_methylation_source_status.tsv": ROOT / "data/processed/epigenetics/methylation/tcga_nise_methylation_source_status.tsv",\n  }',
  )


def patch_counts() -> None:
  outputs = "scripts/validate_publication_outputs.py"
  replacements = (
    ("len(figures) != 77", "len(figures) != 80"),
    ("expected_77", "expected_80"),
    ("range(1, 70)", "range(1, 73)"),
    ("len(tables) != 48", "len(tables) != 51"),
    ("expected_48", "expected_51"),
    ("supplementary_count != 44", "supplementary_count != 47"),
    ("expected_44", "expected_47"),
    ("supplementary figures: 69", "supplementary figures: 72"),
    ("exported image files: 231", "exported image files: 240"),
    ("supplementary tables: 44", "supplementary tables: 47"),
  )
  text = read(outputs)
  for old, new in replacements:
    text = text.replace(old, new)
  write(outputs, text)
  replace_once(
    outputs,
    '    ROOT / "data/processed/regulatory/wgcna/wgcna_run_diagnostics_all_cancers.tsv",',
    '    ROOT / "data/processed/regulatory/wgcna/wgcna_run_diagnostics_all_cancers.tsv",\n    ROOT / "data/processed/epigenetics/methylation/tcga_nise_methylation_pair_metrics.tsv",\n    ROOT / "data/processed/epigenetics/methylation/tcga_nise_methylation_source_status.tsv",\n    article_root / "tables/supplementary/Table_S45_tcga_nise_methylation_gene_summary.tsv",\n    article_root / "tables/supplementary/Table_S46_tcga_nise_methylation_pair_metrics.tsv",\n    article_root / "tables/supplementary/Table_S47_tcga_nise_methylation_source_status.tsv",',
  )

  integrity = "scripts/validate_publication_scientific_integrity.py"
  text = read(integrity)
  text = text.replace("len(manifest) != 77", "len(manifest) != 80")
  text = text.replace("Expected 77 registered figures", "Expected 80 registered figures")
  text = text.replace("!= 44", "!= 47")
  text = text.replace("Expected 44 supplementary tables", "Expected 47 supplementary tables")
  write(integrity, text)
  replace_once(
    integrity,
    '    "Figure_S69_integrated_regulatory_context_source_data.tsv",\n  ):',
    '    "Figure_S69_integrated_regulatory_context_source_data.tsv",\n    "Figure_S70_methylation_gene_coverage_source_data.tsv",\n    "Figure_S71_pair_methylation_context_source_data.tsv",\n    "Figure_S72_methylation_regulatory_integration_source_data.tsv",\n  ):',
  )

  verify = "scripts/verify_complete_article_run.sh"
  text = read(verify)
  text = text.replace(
    r"^Figure_S(?:[1-9]|[1-5][0-9]|6[0-9])$",
    r"^Figure_S(?:[1-9]|[1-6][0-9]|7[0-2])$",
  )
  text = text.replace('"supplementary_figures": 69', '"supplementary_figures": 72')
  text = text.replace('"exported_figure_files": 231', '"exported_figure_files": 240')
  text = text.replace('"supplementary_tables": 44', '"supplementary_tables": 47')
  text = text.replace(
    '  root / "source_data/figures/supplementary/Figure_S69_integrated_regulatory_context_source_data.tsv",',
    '  root / "source_data/figures/supplementary/Figure_S69_integrated_regulatory_context_source_data.tsv",\n  root / "source_data/figures/supplementary/Figure_S70_methylation_gene_coverage_source_data.tsv",\n  root / "source_data/figures/supplementary/Figure_S71_pair_methylation_context_source_data.tsv",\n  root / "source_data/figures/supplementary/Figure_S72_methylation_regulatory_integration_source_data.tsv",',
  )
  write(verify, text)

  readme = read("README.md")
  readme = readme.replace("69 supplementary figures", "72 supplementary figures")
  readme = readme.replace("231 PNG/PDF/SVG files", "240 PNG/PDF/SVG files")
  readme = readme.replace("44 supplementary tables", "47 supplementary tables")
  write("README.md", readme)

  for path in (ROOT / "tests").glob("test_*.py"):
    content = path.read_text(encoding="utf-8")
    content = content.replace("range(1, 70)", "range(1, 73)")
    content = re.sub(r"(?<![A-Za-z0-9_])77(?![A-Za-z0-9_])", "80", content)
    content = re.sub(r"(?<![A-Za-z0-9_])231(?![A-Za-z0-9_])", "240", content)
    content = re.sub(r"(?<![A-Za-z0-9_])44(?![A-Za-z0-9_])", "47", content)
    content = re.sub(r"(?<![A-Za-z0-9_])48(?![A-Za-z0-9_])", "51", content)
    content = content.replace("S39-S69", "S39-S72")
    path.write_text(content, encoding="utf-8")


def append_once(path: str, marker: str, block: str) -> None:
  text = read(path)
  if marker not in text:
    write(path, text.rstrip() + "\n\n" + block.strip() + "\n")


def patch_docs() -> None:
  appendix = '''## TCGA/GDC DNA methylation layer (v0.11.1)

RSES-Onco now acquires gene-associated CpG methylation beta values from the NCI Genomic Data Commons through the UCSC Xena GDC hub. Repbase is not used because it is a reference library of repetitive DNA sequences rather than a sample-level methylation resource.

```bash
python -u scripts/acquire_tcga_nise_methylation.py \\
  --candidates data/processed/expanded_candidate_universe.tsv \\
  --output-dir data/processed/epigenetics/methylation
```

Methylation is integrated inside the existing regulatory-network microniche weight. It does not receive a new independent top-level RSES-Onco domain. Regulatory subweights are DoRothEA regulator divergence 0.32, TF-expression-profile divergence 0.28, JASPAR motif divergence 0.20 and TCGA/GDC methylation context 0.20. A technical source failure makes methylation non-eligible and preserves the original three-component regulatory score. Available source data with missing pair-level probes or samples reduce internal regulatory coverage.

New outputs are Figures S70-S72 and Supplementary Tables S45-S47.
'''
  for path in (
    "docs/END_TO_END_ARTICLE_PROTOCOL.md",
    "docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md",
    "supplementary/Supplementary_Methods_RSES_Onco_v0110.md",
    "manuscript/RSES_Onco_intro_methods_draft_v0110.md",
  ):
    append_once(
      path,
      "## TCGA/GDC DNA methylation layer (v0.11.1)",
      appendix,
    )
  append_once(
    "supplementary/Supplementary_Methods_RSES_Onco_v0110.md",
    "## SeSAMe methylation references (v0.11.1)",
    '''## SeSAMe methylation references (v0.11.1)

- Zhou W, Triche TJ Jr, Laird PW, Shen H. SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. *Nucleic Acids Research*. 2018;46:e123.
- Zhou W, Laird PW, Shen H. Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. *Nucleic Acids Research*. 2017;45:e22.
- Goldman MJ, Craft B, Hastie M, et al. Visualizing and interpreting cancer genomics data via the Xena platform. *Nature Biotechnology*. 2020;38:675-678.
''',
  )
  append_once(
    "README.md",
    "## TCGA/GDC methylation integration",
    '''## TCGA/GDC methylation integration

Version 0.11.1 adds gene-associated TCGA/GDC CpG methylation through the UCSC Xena GDC hub as an internal regulatory-network sublayer. Repbase is not used because it is a repetitive-element sequence library. The total top-level and microniche domain weights remain unchanged; the regulatory weight is internally shared with methylation to prevent double counting.
''',
  )


def main() -> None:
  patch_config()
  patch_orchestrator()
  patch_tables()
  patch_counts()
  patch_docs()
  print("Applied publication and documentation methylation patch.")


if __name__ == "__main__":
  main()
