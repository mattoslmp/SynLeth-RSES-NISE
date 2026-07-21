#!/usr/bin/env python3
"""Build a non-phylogenetic human enzyme relation forest from the curated NISE table.

Leaves are grouped first by reported structural cluster and then by EC activity.
This is a functional/annotation topology for navigation, not an evolutionary tree.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def safe_label(text: object) -> str:
  return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text).strip()).strip("_")


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--catalog", default="data/curated/human_nise_bonafide_2017.tsv")
  parser.add_argument("--output", default="results/human_nise_relation_forest.nwk")
  parser.add_argument("--edges", default="results/human_nise_relation_edges.tsv")
  args = parser.parse_args()
  catalog = Path(args.catalog)
  if not catalog.is_absolute():
    catalog = ROOT / catalog
  frame = pd.read_csv(catalog, sep="\t")
  activity_nodes = []
  edges = []
  for group_id, group in frame.groupby("group_id", sort=True):
    cluster_nodes = []
    activity_label = f"{safe_label(group_id)}_EC_{safe_label(group.ec_number.iloc[0])}"
    for cluster, members in group.groupby("structural_cluster", sort=True):
      leaves = []
      cluster_label = f"{activity_label}_cluster_{safe_label(cluster)}"
      for row in members.itertuples():
        leaf = f"{safe_label(row.gene_symbol)}_{safe_label(row.uniprot_accession)}"
        leaves.append(leaf)
        edges.append({"parent": cluster_label, "child": leaf, "edge_type": "reported_structural_cluster_membership", "group_id": group_id, "ec_number": row.ec_number})
      cluster_nodes.append(f"({','.join(leaves)}){cluster_label}")
      edges.append({"parent": activity_label, "child": cluster_label, "edge_type": "activity_to_structural_cluster", "group_id": group_id, "ec_number": group.ec_number.iloc[0]})
    activity_nodes.append(f"({','.join(cluster_nodes)}){activity_label}")
  forest = f"({','.join(activity_nodes)})Human_NISE_functional_forest;\n"
  output = Path(args.output)
  if not output.is_absolute():
    output = ROOT / output
  output.parent.mkdir(parents=True, exist_ok=True)
  output.write_text(forest, encoding="utf-8")
  edge_path = Path(args.edges)
  if not edge_path.is_absolute():
    edge_path = ROOT / edge_path
  pd.DataFrame(edges).to_csv(edge_path, sep="\t", index=False)
  print(f"Wrote annotation relation forest to {output}")
  print(f"Wrote {len(edges)} relation edges to {edge_path}")
  print("This topology is not a sequence phylogeny.")


if __name__ == "__main__":
  main()
