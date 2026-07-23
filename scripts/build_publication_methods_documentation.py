#!/usr/bin/env python3
"""Generate human-readable scientific methods and asset-reproduction documentation."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from rses_onco.expanded import EXPANDED_ONCO_WEIGHTS, FUNCTIONAL_MICRONICHE_WEIGHTS
from rses_onco.evidence_categories import evidence_category_definitions

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()
  article_root = resolve_path(args.article_root)
  methods_dir = article_root / "manuscript_assets" / "supplementary_methods"
  methods_dir.mkdir(parents=True, exist_ok=True)

  formula = r"""# RSES-Onco score definition and evidence rules

## Coverage-aware score

For an eligible component $i$ with weight $w_i$ and observed normalized value $x_i \in [0,1]$:

$S_{observed} = \frac{\sum_{i \in O} w_i x_i}{\sum_{i \in O} w_i}$

$C = \frac{\sum_{i \in O} w_i}{\sum_{i \in E} w_i}$

$S_{adjusted} = S_{observed} C = \frac{\sum_{i \in O} w_i x_i}{\sum_{i \in E} w_i}$

Here, $O$ contains observed eligible components and $E$ contains all components eligible under the score definition. A missing component is not converted to zero. A non-eligible component does not enter the observed denominator. Technical source failure is not biological negative evidence. An observed value of zero is retained as observed negative evidence.

## Evidence states

- **Observed evidence:** a numeric component was calculated from available evidence.
- **Observed negative evidence:** an eligible observed component equals zero.
- **Observed neutral evidence:** an eligible observed component equals the operational midpoint.
- **Missing evidence:** no eligible evidence was available and no technical failure was documented.
- **Insufficient sample:** the analysis could not satisfy the required sample size.
- **Technical/source failure:** acquisition or mapping failed; this is not scored as zero.
- **Not eligible:** the candidate/context cannot be evaluated by that domain, for example a composite event that cannot be represented as a single-gene loss.

## Independence and overlap

Evidence representations sharing a publication, original dataset or traceable evidence identifier are grouped into one evidence unit. Repeated aggregators may confirm interpretation but cannot receive more than one combined evidence weight. Score evidence, prioritization evidence, independent validation and interpretative evidence are explicitly distinguished.

## Terminology boundaries

A candidate in the universe is not automatically a discovery. A prioritized hypothesis is not automatically statistically significant. Nominal significance is distinct from FDR support. Conditional dependency is distinct from general essentiality. Experimental tractability is not clinical efficacy.
"""
  (methods_dir / "RSES_Onco_score_formula_and_missing_data_rules.md").write_text(formula, encoding="utf-8")

  weights = []
  for family, mapping in (
    ("RSES-Onco", EXPANDED_ONCO_WEIGHTS),
    ("Functional microniche", FUNCTIONAL_MICRONICHE_WEIGHTS),
  ):
    for domain, weight in mapping.items():
      weights.append({"score_family": family, "domain": domain, "weight": weight})
  pd.DataFrame(weights).to_csv(methods_dir / "RSES_Onco_domain_weights.tsv", sep="\t", index=False)

  inventory_path = article_root / "tables/figure_data/figure_source_data_inventory.tsv"
  table_manifest_path = article_root / "manifests/table_manifest.tsv"
  rows = []
  if inventory_path.exists() and inventory_path.stat().st_size:
    inventory = pd.read_csv(inventory_path, sep="\t", low_memory=False)
    for record in inventory.to_dict("records"):
      rows.append({
        "asset": record.get("figure_id"),
        "asset_type": "figure",
        "script": record.get("generator_script"),
        "inputs": record.get("input_paths"),
        "intermediate_or_source_table": record.get("source_table"),
        "command": record.get("reproduction_command"),
        "outputs": record.get("figure_base_path"),
        "dependencies": "Python environment defined by environment.yml and pyproject.toml",
        "article_location": record.get("category"),
      })
  if table_manifest_path.exists() and table_manifest_path.stat().st_size:
    tables = pd.read_csv(table_manifest_path, sep="\t", low_memory=False)
    for record in tables.to_dict("records"):
      rows.append({
        "asset": record.get("table_id"),
        "asset_type": "table",
        "script": record.get("script"),
        "inputs": record.get("source_paths"),
        "intermediate_or_source_table": record.get("path"),
        "command": "MPLBACKEND=Agg bash scripts/run_publication_pipeline.sh assets-only",
        "outputs": record.get("path"),
        "dependencies": "Python environment defined by environment.yml and pyproject.toml",
        "article_location": record.get("category"),
      })
  reproduction = pd.DataFrame(rows)
  reproduction.to_csv(methods_dir / "Table_asset_script_input_output_reproduction.tsv", sep="\t", index=False)
  markdown = [
    "# Publication asset reproduction", "",
    "| Asset | Type | Script | Inputs | Exact source/intermediate table | Command | Output |",
    "|---|---|---|---|---|---|---|",
  ]
  for row in rows:
    markdown.append(
      f"| {row['asset']} | {row['asset_type']} | `{row['script']}` | `{row['inputs']}` | "
      f"`{row['intermediate_or_source_table']}` | `{row['command']}` | `{row['outputs']}` |"
    )
  (methods_dir / "PUBLICATION_ASSET_REPRODUCTION.md").write_text("\n".join(markdown) + "\n", encoding="utf-8")

  definitions = evidence_category_definitions()
  definitions.to_csv(
    methods_dir / "evidence_category_definitions.tsv", sep="\t", index=False
  )

  registry_path = (
    article_root
    / "tables/supplementary/Table_S44_asset_reproduction_registry.tsv"
  )
  registry_path.parent.mkdir(parents=True, exist_ok=True)
  reproduction.to_csv(registry_path, sep="\t", index=False)
  if table_manifest_path.exists() and table_manifest_path.stat().st_size:
    table_manifest = pd.read_csv(table_manifest_path, sep="\t", low_memory=False)
    mask = table_manifest["table_id"].astype(str).eq(registry_path.stem)
    if mask.any():
      table_manifest.loc[mask, "rows"] = len(reproduction)
      table_manifest.loc[mask, "columns"] = len(reproduction.columns)
      table_manifest.loc[mask, "source_paths"] = str(inventory_path)
      table_manifest.loc[mask, "script"] = "scripts/build_publication_methods_documentation.py"
      table_manifest.loc[mask, "status"] = "ok" if not reproduction.empty else "empty_no_eligible_records"
      table_manifest.to_csv(table_manifest_path, sep="\t", index=False)
  print(f"Wrote scientific methods and reproduction documentation to {methods_dir}")


if __name__ == "__main__":
  main()
