#!/usr/bin/env python3
"""Export all available supporting evidence into publication-organized TSV tables.

The exporter never synthesizes evidence. Existing tables are copied with provenance;
missing sources receive a status row explaining that the corresponding analysis is
unavailable in the current release. Network node/edge tables are produced only from
real STRING/DoRothEA rows already acquired by the pipeline.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

import pandas as pd
from pandas.errors import EmptyDataError

ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class SupportingTableRecord:
  evidence_family: str
  output_path: str
  source_path: str
  status: str
  rows: int
  columns: int
  sha256: str
  interpretation_boundary: str


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_optional(path: Path) -> pd.DataFrame:
  if not path.exists():
    return pd.DataFrame()
  try:
    return pd.read_csv(path, sep="\t", low_memory=False)
  except EmptyDataError:
    return pd.DataFrame()


def file_sha256(path: Path) -> str:
  digest = hashlib.sha256()
  with path.open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> Path:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)
  return path


def export_table(
  family: str,
  source: Path,
  output: Path,
  boundary: str,
) -> SupportingTableRecord:
  frame = read_optional(source)
  output.parent.mkdir(parents=True, exist_ok=True)
  if source.exists():
    atomic_tsv(frame, output)
    status = "available" if not frame.empty else "available_but_empty_no_eligible_records"
  else:
    frame = pd.DataFrame([{
      "evidence_status": "source_file_absent",
      "absence_reason": "The upstream analysis did not produce this source table in the current release.",
      "source_path": str(source),
      "scientific_interpretation": "Unavailable data are not interpreted as negative biological evidence.",
    }])
    atomic_tsv(frame, output)
    status = "source_file_absent"
  return SupportingTableRecord(
    evidence_family=family,
    output_path=str(output),
    source_path=str(source),
    status=status,
    rows=len(frame),
    columns=len(frame.columns),
    sha256=file_sha256(output),
    interpretation_boundary=boundary,
  )


def build_network_tables(functional: pd.DataFrame, output_dir: Path) -> list[SupportingTableRecord]:
  records: list[SupportingTableRecord] = []
  if functional.empty:
    status = pd.DataFrame([{
      "evidence_status": "no_functional_evidence_rows",
      "absence_reason": "The functional-evidence table is empty.",
    }])
    path = atomic_tsv(status, output_dir / "network_evidence_status.tsv")
    records.append(SupportingTableRecord(
      "networks", str(path), "data/processed/expanded_pair_functional_evidence.tsv",
      "available_but_empty_no_eligible_records", len(status), len(status.columns),
      file_sha256(path), "No network edge is inferred from an empty table.",
    ))
    return records

  pair_columns = [
    column for column in (
      "pair_id", "lost_gene", "target_gene", "source_class",
      "string_direct_score", "string_neighbor_jaccard", "string_shared_neighbors",
      "component_interaction_network", "regulator_jaccard", "shared_regulators",
      "component_regulatory_network", "component_localization",
      "component_biochemical_structural",
    ) if column in functional
  ]
  pair_table = functional[pair_columns].copy()
  pair_path = atomic_tsv(pair_table, output_dir / "pair_level_network_and_microniche_evidence.tsv")
  records.append(SupportingTableRecord(
    "pair_level_network_evidence", str(pair_path),
    "data/processed/expanded_pair_functional_evidence.tsv", "available",
    len(pair_table), len(pair_table.columns), file_sha256(pair_path),
    "STRING combined/network metrics are not described as direct experimental interactions.",
  ))

  node_rows = []
  for record in functional.to_dict("records"):
    for role, column in (("origin", "lost_gene"), ("target", "target_gene")):
      value = record.get(column)
      if pd.notna(value) and str(value).strip():
        node_rows.append({
          "node_id": str(value).strip(),
          "role": role,
          "pair_id": record.get("pair_id"),
          "source_class": record.get("source_class"),
        })
  nodes = pd.DataFrame(node_rows).drop_duplicates() if node_rows else pd.DataFrame(
    columns=["node_id", "role", "pair_id", "source_class"]
  )
  nodes_path = atomic_tsv(nodes, output_dir / "network_nodes.tsv")
  records.append(SupportingTableRecord(
    "network_nodes", str(nodes_path),
    "data/processed/expanded_pair_functional_evidence.tsv",
    "available" if not nodes.empty else "available_but_empty_no_eligible_records",
    len(nodes), len(nodes.columns), file_sha256(nodes_path),
    "Nodes are copied from acquired pair evidence; no node is inferred.",
  ))

  edge_rows = []
  for record in functional.to_dict("records"):
    source = record.get("lost_gene")
    target = record.get("target_gene")
    if pd.isna(source) or pd.isna(target):
      continue
    if pd.notna(record.get("string_direct_score")) or pd.notna(record.get("string_neighbor_jaccard")):
      edge_rows.append({
        "pair_id": record.get("pair_id"),
        "source_node": source,
        "target_node": target,
        "network_type": "STRING_pair_relationship",
        "direct_score": record.get("string_direct_score"),
        "neighborhood_jaccard": record.get("string_neighbor_jaccard"),
        "shared_neighbors": record.get("string_shared_neighbors"),
        "evidence_channel_interpretation": "Combined/pair-level STRING evidence; not automatically experimental.",
      })
    if pd.notna(record.get("regulator_jaccard")) or pd.notna(record.get("shared_regulators")):
      edge_rows.append({
        "pair_id": record.get("pair_id"),
        "source_node": source,
        "target_node": target,
        "network_type": "shared_regulatory_neighborhood",
        "direct_score": record.get("component_regulatory_network"),
        "neighborhood_jaccard": record.get("regulator_jaccard"),
        "shared_neighbors": record.get("shared_regulators"),
        "evidence_channel_interpretation": "Inferred TF-target neighborhood; promoter binding is not claimed unless separately sourced.",
      })
  edges = pd.DataFrame(edge_rows)
  if edges.empty:
    edges = pd.DataFrame(columns=[
      "pair_id", "source_node", "target_node", "network_type", "direct_score",
      "neighborhood_jaccard", "shared_neighbors", "evidence_channel_interpretation",
    ])
  edges_path = atomic_tsv(edges, output_dir / "network_edges.tsv")
  records.append(SupportingTableRecord(
    "network_edges", str(edges_path),
    "data/processed/expanded_pair_functional_evidence.tsv",
    "available" if not edges.empty else "available_but_empty_no_eligible_records",
    len(edges), len(edges.columns), file_sha256(edges_path),
    "Regulatory associations are distinguished from direct promoter-binding evidence.",
  ))
  return records


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument("--ranking", default="results/expanded_26Q1/full/expanded_rses_onco.tsv")
  parser.add_argument("--functional-evidence", default="data/processed/expanded_pair_functional_evidence.tsv")
  args = parser.parse_args()

  output_root = resolve_path(args.output_root)
  evidence_dir = output_root / "tables" / "supporting_evidence"
  records: list[SupportingTableRecord] = []

  specifications = [
    (
      "expression_compensation",
      resolve_path("results/expanded_26Q1/full/expanded_expression_compensation.tsv"),
      evidence_dir / "expression" / "expression_compensation.tsv",
      "A between-group expression contrast supports operational compensation; correlation alone is not treated as compensation.",
    ),
    (
      "expression_context",
      resolve_path("results/expanded_26Q1/full/expanded_expression_context_profiles.tsv"),
      evidence_dir / "expression" / "expression_context_profiles.tsv",
      "Expression-profile divergence is descriptive and context-specific.",
    ),
    (
      "crispr_conditional_dependency",
      resolve_path("results/expanded_26Q1/full/expanded_dependency_contrasts.tsv"),
      evidence_dir / "phenotypes" / "conditional_dependency_contrasts.tsv",
      "Conditional dependency is distinguished from general gene essentiality.",
    ),
    (
      "crispr_phenotype_profiles",
      resolve_path("results/expanded_26Q1/full/expanded_crispr_phenotype_profiles.tsv"),
      evidence_dir / "phenotypes" / "crispr_phenotype_profiles.tsv",
      "Phenotype-profile divergence is not itself proof of synthetic lethality.",
    ),
    (
      "tumor_events",
      resolve_path("results/expanded_26Q1/full/tcga_gene_event_summary.tsv"),
      evidence_dir / "genomic_context" / "tcga_gene_event_summary.tsv",
      "Tumor-event frequency is reported with evaluable sample counts and is never allowed outside [0,1].",
    ),
    (
      "all_target_discovery",
      resolve_path("results/expanded_26Q1/discovery/all_target_dependency_screen.tsv"),
      evidence_dir / "phenotypes" / "all_target_dependency_screen.tsv",
      "Nominal and FDR-supported results remain explicitly distinguished.",
    ),
    (
      "pharmacology_evidence",
      resolve_path("data/processed/pharmacology/pharmacology_evidence_long.tsv"),
      evidence_dir / "pharmacology" / "pharmacology_evidence_long.tsv",
      "Tractability and compound association are experimental-priority evidence, not clinical efficacy.",
    ),
    (
      "drug_response_selectivity",
      resolve_path("data/processed/pharmacology/drug_response_selectivity.tsv"),
      evidence_dir / "pharmacology" / "drug_response_selectivity.tsv",
      "An empty table means no eligible local PRISM/GDSC/CTRP contrast, not absence of pharmacological activity.",
    ),
    (
      "pharmacology_priorities",
      resolve_path("results/expanded_26Q1/pharmacology/pharmacology_ranked_hypotheses.tsv"),
      evidence_dir / "pharmacology" / "pharmacology_ranked_hypotheses.tsv",
      "Computational treatment hypotheses require biomarker-matched experimental validation.",
    ),
    (
      "alphafold_structure_manifest",
      resolve_path("data/processed/structures/alphafold_structure_manifest.tsv"),
      evidence_dir / "structures" / "alphafold_structure_manifest.tsv",
      "Predicted structures are not experimental ligand-binding evidence.",
    ),
    (
      "structural_residue_annotations",
      resolve_path("data/processed/structures/nise_structural_residue_annotations.tsv"),
      evidence_dir / "structures" / "structural_residue_annotations.tsv",
      "Only traceable exact-numbered residues are projected; missing annotations remain missing.",
    ),
    (
      "structure_render_manifest",
      resolve_path("data/processed/structures/nise_structure_render_manifest.tsv"),
      evidence_dir / "structures" / "structure_render_manifest.tsv",
      "Rendering provenance does not add biological evidence.",
    ),
  ]
  for family, source, output, boundary in specifications:
    records.append(export_table(family, source, output, boundary))

  functional_path = resolve_path(args.functional_evidence)
  functional = read_optional(functional_path)
  functional_output = evidence_dir / "networks" / "expanded_pair_functional_evidence.tsv"
  records.append(export_table(
    "functional_microniche_pair_evidence",
    functional_path,
    functional_output,
    "Localization, biochemical/structural and network components retain source-specific limitations.",
  ))
  records.extend(build_network_tables(functional, evidence_dir / "networks"))

  ranking_path = resolve_path(args.ranking)
  records.append(export_table(
    "complete_score_components",
    ranking_path,
    evidence_dir / "scores" / "expanded_rses_onco.tsv",
    "Missing values remain missing and reduce explicit coverage.",
  ))

  manifest = pd.DataFrame([asdict(record) for record in records])
  manifest["exported_at_utc"] = datetime.now(timezone.utc).isoformat()
  manifest_path = atomic_tsv(manifest, evidence_dir / "supporting_evidence_manifest.tsv")
  if manifest.empty or not manifest_path.exists() or manifest_path.stat().st_size == 0:
    raise RuntimeError("Supporting-evidence manifest is missing or empty")
  for record in records:
    path = Path(record.output_path)
    if not path.exists() or path.stat().st_size == 0:
      raise RuntimeError(f"Supporting evidence output missing or empty: {path}")
  print(manifest[["evidence_family", "status", "rows", "output_path"]].to_string(index=False))
  print(f"Wrote supporting-evidence manifest to {manifest_path}")


if __name__ == "__main__":
  main()
