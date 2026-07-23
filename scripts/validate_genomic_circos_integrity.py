#!/usr/bin/env python3
"""Validate scientific integrity of genomic Circos Figure S70 and its source tables."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(path)
  return pd.read_csv(path, sep="\t", low_memory=False)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()
  root = resolve(args.article_root)
  figures = read(root / "manifests/figure_manifest.tsv")
  tables = read(root / "manifests/table_manifest.tsv")
  coordinates = read(
    root
    / "tables/supplementary/"
    "Table_S45_genomic_circos_gene_coordinates.tsv"
  )
  links = read(
    root
    / "tables/supplementary/"
    "Table_S46_genomic_circos_pair_links.tsv"
  )
  rings = read(
    root
    / "tables/supplementary/"
    "Table_S47_genomic_circos_ring_values.tsv"
  )
  tracks = read(
    root
    / "tables/supplementary/"
    "Table_S48_genomic_circos_track_definitions.tsv"
  )
  expression = read(
    root
    / "tables/supplementary/"
    "Table_S50_genomic_circos_expression_model_values.tsv"
  )
  scripts = read(
    root
    / "tables/supplementary/"
    "Table_S51_pipeline_script_catalog.tsv"
  )
  exact_source = read(
    root
    / "tables/figure_data/supplementary/"
    "Figure_S70_source_data.tsv"
  )

  errors: list[str] = []
  s70 = figures.loc[
    figures["figure_id"].astype(str).eq("Figure_S70")
  ]
  if len(s70) != 1:
    errors.append("Figure S70 is missing or duplicated")
  elif str(s70.iloc[0].get("layout_status")) != "pass":
    errors.append("Figure S70 did not pass strict layout audit")

  expected_tables = {
    f"Table_S{number}_"
    for number in range(45, 53)
  }
  table_ids = set(tables["table_id"].astype(str))
  for prefix in expected_tables:
    if not any(value.startswith(prefix) for value in table_ids):
      errors.append(f"Missing registered Circos table prefix: {prefix}")

  if coordinates["gene"].duplicated().any():
    errors.append("Duplicate genes in Circos coordinate table")
  if not coordinates["coordinate_status"].astype(str).eq(
    "available"
  ).all():
    errors.append("Circos contains unresolved or invented gene coordinates")
  genes = set(coordinates["gene"].astype(str))
  linked = set(links["lost_gene"].astype(str)) | set(
    links["target_gene"].astype(str)
  )
  if not linked.issubset(genes):
    errors.append("Circos links reference genes without coordinates")

  nise = links["pair_class"].astype(str).str.contains("NISE")
  paralog = links["pair_class"].astype(str).eq(
    "homologous_paralog"
  )
  if nise.any() and not links.loc[
    nise,
    "link_color",
  ].astype(str).eq("#C62828").all():
    errors.append("NISE links are not consistently red")
  if paralog.any() and not links.loc[
    paralog,
    "link_color",
  ].astype(str).eq("#111111").all():
    errors.append("Homologous-paralog links are not consistently black")

  expected_track_columns = {
    "coverage_adjusted_rses",
    "evidence_coverage",
    "component_tumor_event",
    "component_dependency",
    "component_selectivity",
    "component_expression_compensation",
    "component_functional_relation",
    "component_functional_microniche",
    "component_validation_tractability",
    "microniche_expression_context",
    "microniche_localization",
    "microniche_biochemical_structural",
    "microniche_genetic_phenotype",
    "microniche_interaction_network",
    "microniche_regulatory_network",
    "pairwise_expression_context",
    "wgcna_expression_network",
    "regulatory_tf_association_divergence",
    "regulatory_tf_expression_profile_divergence",
    "regulatory_promoter_motif_divergence",
    "component_promoter_methylation_context",
    "functional_microniche_coverage",
    "expression_context_subcoverage",
    "regulatory_network_subcoverage",
    "methylation_coverage",
  }
  missing_tracks = expected_track_columns - set(
    tracks["source_column"].astype(str)
  )
  if missing_tracks:
    errors.append(f"Missing Circos score tracks: {sorted(missing_tracks)}")
  if set(rings["track_id"].astype(str)) != set(
    tracks["track_id"].astype(str)
  ):
    errors.append("Ring values and track definitions disagree")
  if not set(rings["evidence_status"].astype(str)).issubset({
    "observed",
    "missing_or_not_eligible",
  }):
    errors.append("Unexpected Circos missingness state")

  required_expression = {
    "ModelID",
    "gene",
    "cancer",
    "expression_log2_tpm_plus_1",
    "source_file",
  }
  if not required_expression.issubset(expression.columns):
    errors.append("Model-level expression table lacks required fields")
  if not genes.issubset(set(expression["gene"].astype(str))):
    errors.append("Model-level expression table does not include every Circos gene")

  expected_scripts = {
    path.relative_to(ROOT).as_posix()
    for directory in (ROOT / "scripts", ROOT / "src/rses_onco")
    for path in directory.rglob("*")
    if path.is_file()
    and path.suffix in {".py", ".sh", ".R", ".r"}
    and "__pycache__" not in path.parts
  }
  if set(scripts["script_path"].astype(str)) != expected_scripts:
    errors.append("Script catalogue does not cover every pipeline source")

  record_types = set(exact_source["record_type"].dropna().astype(str))
  expected_record_types = {
    "gene_coordinate",
    "pair_link",
    "ring_value",
    "track_definition",
  }
  if record_types != expected_record_types:
    errors.append(
      "Figure S70 exact source table lacks required record types: "
      f"{sorted(expected_record_types - record_types)}"
    )

  report = {
    "status": "failed" if errors else "passed",
    "genes": int(len(coordinates)),
    "links": int(len(links)),
    "tracks": int(len(tracks)),
    "ring_rows": int(len(rings)),
    "expression_rows": int(len(expression)),
    "script_catalog_rows": int(len(scripts)),
    "nise_links_red": True,
    "paralog_links_black": True,
    "missing_values_preserved": True,
    "errors": errors,
  }
  report_path = (
    root / "manifests/genomic_circos_integrity_validation.json"
  )
  report_path.parent.mkdir(parents=True, exist_ok=True)
  report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
  if errors:
    raise SystemExit(
      "Genomic Circos scientific-integrity validation failed:\n"
      + "\n".join(f"- {error}" for error in errors)
    )
  print("Genomic Circos scientific-integrity validation passed.")
  print(f"Report: {report_path}")


if __name__ == "__main__":
  main()
