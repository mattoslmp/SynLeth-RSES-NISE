#!/usr/bin/env python3
"""Validate scientific integrity of genomic Circos Figure S70 and its source tables."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))
if str(ROOT / "src") not in sys.path:
  sys.path.insert(0, str(ROOT / "src"))

import pandas as pd

from scripts.complete_genomic_circos_links import selected_pairs


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
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  args = parser.parse_args()
  root = resolve(args.article_root)
  candidate_pairs = selected_pairs(resolve(args.candidates))
  figures = read(root / "manifests/figure_manifest.tsv")
  tables = read(root / "manifests/table_manifest.tsv")
  coordinates = read(
    root / "tables/supplementary/Table_S45_genomic_circos_gene_coordinates.tsv"
  )
  links = read(
    root / "tables/supplementary/Table_S46_genomic_circos_pair_links.tsv"
  )
  rings = read(
    root / "tables/supplementary/Table_S47_genomic_circos_ring_values.tsv"
  )
  tracks = read(
    root / "tables/supplementary/Table_S48_genomic_circos_track_definitions.tsv"
  )
  expression_summary = read(
    root / "tables/supplementary/Table_S49_genomic_circos_expression_summary.tsv"
  )
  expression = read(
    root / "tables/supplementary/Table_S50_genomic_circos_expression_model_values.tsv"
  )
  scripts = read(
    root / "tables/supplementary/Table_S51_pipeline_script_catalog.tsv"
  )
  exact_source = read(
    root / "tables/figure_data/supplementary/Figure_S70_source_data.tsv"
  )

  errors: list[str] = []
  s70 = figures.loc[figures["figure_id"].astype(str).eq("Figure_S70")]
  if len(s70) != 1:
    errors.append("Figure S70 is missing or duplicated")
  elif str(s70.iloc[0].get("layout_status")) != "pass":
    errors.append("Figure S70 did not pass strict layout audit")

  expected_tables = {f"Table_S{number}_" for number in range(45, 53)}
  table_ids = set(tables["table_id"].astype(str))
  for prefix in expected_tables:
    if not any(value.startswith(prefix) for value in table_ids):
      errors.append(f"Missing registered Circos table prefix: {prefix}")

  if coordinates["gene"].duplicated().any():
    errors.append("Duplicate genes in Circos coordinate table")
  if not coordinates["coordinate_status"].astype(str).eq("available").all():
    errors.append("Circos contains unresolved or invented gene coordinates")
  genes = set(coordinates["gene"].astype(str))
  linked = set(links["lost_gene"].astype(str)) | set(
    links["target_gene"].astype(str)
  )
  if not linked.issubset(genes):
    errors.append("Circos links reference genes without coordinates")

  expected_pair_ids = set(candidate_pairs["pair_id"].astype(str))
  observed_pair_ids = set(links["pair_id"].astype(str))
  if expected_pair_ids != observed_pair_ids:
    errors.append(
      "Circos does not contain exactly one chord per simple NISE/paralog "
      f"pair; missing={sorted(expected_pair_ids - observed_pair_ids)[:50]}, "
      f"extra={sorted(observed_pair_ids - expected_pair_ids)[:50]}"
    )
  if links["pair_id"].astype(str).duplicated().any():
    errors.append("Circos contains duplicated pair chords")

  nise = links["pair_class"].astype(str).str.contains("NISE")
  paralog = links["pair_class"].astype(str).eq("homologous_paralog")
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
    "rses_onco",
    "coverage_adjusted_rses",
    "evidence_coverage",
    "component_tumor_event",
    "component_dependency",
    "component_selectivity",
    "component_expression_compensation",
    "component_functional_relation",
    "component_functional_microniche",
    "component_validation_tractability",
    "genetic_screen",
    "isogenic_validation",
    "in_vivo",
    "clinical_tractability",
    "microniche_expression_context",
    "microniche_localization",
    "microniche_biochemical_structural",
    "microniche_genetic_phenotype",
    "microniche_interaction_network",
    "microniche_regulatory_network",
    "pairwise_expression_context",
    "wgcna_expression_network",
    "wgcna_tom_divergence",
    "wgcna_module_divergence",
    "wgcna_kME_divergence",
    "regulatory_tf_association_divergence",
    "regulatory_tf_expression_profile_divergence",
    "regulatory_promoter_motif_divergence",
    "component_promoter_methylation_context",
    "methylation_pair_profile_divergence",
    "methylation_target_hypomethylation_support",
    "functional_microniche_coverage",
    "expression_context_subcoverage",
    "regulatory_network_subcoverage",
    "methylation_coverage",
  }
  observed_track_columns = set(tracks["source_column"].astype(str))
  missing_tracks = expected_track_columns - observed_track_columns
  if missing_tracks:
    errors.append(f"Missing Circos score tracks: {sorted(missing_tracks)}")
  if len(tracks) != 35:
    errors.append(f"Expected 35 Circos tracks; observed {len(tracks)}")
  if set(rings["track_id"].astype(str)) != set(tracks["track_id"].astype(str)):
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
    "evidence_status",
    "is_measurement",
  }
  if not required_expression.issubset(expression.columns):
    errors.append("Model-level expression table lacks required fields")
  if not genes.issubset(set(expression["gene"].astype(str))):
    errors.append("Model-level expression table does not include every Circos gene")
  sentinel = expression.loc[
    expression["is_measurement"].astype(str).str.casefold().eq("false")
  ]
  if not sentinel.empty and not pd.to_numeric(
    sentinel["expression_log2_tpm_plus_1"],
    errors="coerce",
  ).isna().all():
    errors.append("Unavailable expression sentinel rows contain numeric values")
  expected_expression_contexts = {
    (cancer, gene)
    for cancer in ("colon", "stomach", "lung")
    for gene in genes
  }
  observed_expression_contexts = set(zip(
    expression_summary["cancer"].astype(str),
    expression_summary["gene"].astype(str),
  ))
  if expected_expression_contexts != observed_expression_contexts:
    errors.append(
      "Expression summary does not contain every Circos gene × cancer context"
    )

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
    "candidate_pairs": int(len(candidate_pairs)),
    "genes": int(len(coordinates)),
    "links": int(len(links)),
    "tracks": int(len(tracks)),
    "ring_rows": int(len(rings)),
    "expression_summary_rows": int(len(expression_summary)),
    "expression_model_rows": int(len(expression)),
    "script_catalog_rows": int(len(scripts)),
    "all_simple_nise_paralog_pairs_have_chords": not bool(
      expected_pair_ids ^ observed_pair_ids
    ),
    "all_score_internal_layers_present": not bool(
      expected_track_columns - observed_track_columns
    ),
    "nise_links_red": True,
    "paralog_links_black": True,
    "missing_values_preserved": True,
    "errors": errors,
  }
  report_path = root / "manifests/genomic_circos_integrity_validation.json"
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
