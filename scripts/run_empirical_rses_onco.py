#!/usr/bin/env python3
"""Run RSES-Onco against user-supplied DepMap and TCGA matrices."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from rses_onco.depmap import dependency_components, dependency_contrast, read_depmap_inputs
from rses_onco.integrate import load_reference_candidates, score_candidate_table
from rses_onco.tcga import event_component, homozygous_deletion_frequency, read_gistic_matrix

LINEAGES = {"colon": "Colorectal", "stomach": "Gastric", "lung": "Lung"}
ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--gene-effect", required=True)
  parser.add_argument("--copy-number", required=True)
  parser.add_argument("--models", required=True)
  parser.add_argument("--expression")
  parser.add_argument("--tcga", action="append", default=[], help="cancer=path/to/discrete_cna.tsv")
  parser.add_argument("--candidates", default="data/curated/synthetic_lethality_reference_pairs.tsv")
  parser.add_argument("--output", default="results/empirical_rses_onco.tsv")
  args = parser.parse_args()
  effect, copy_number, models, expression = read_depmap_inputs(
    resolve_path(args.gene_effect), resolve_path(args.copy_number), resolve_path(args.models), resolve_path(args.expression) if args.expression else None
  )
  tcga = {}
  for item in args.tcga:
    cancer, path = item.split("=", 1)
    tcga[cancer] = read_gistic_matrix(resolve_path(path))
  candidates = load_reference_candidates(resolve_path(args.candidates))
  empirical = {}
  details = []
  for row in candidates.to_dict("records"):
    lost = str(row["lost_feature"])
    target = str(row["target_gene"])
    if any(token in lost for token in ("/", " ", "deficiency", "phenotype")):
      continue
    components = {}
    for cancer, lineage in LINEAGES.items():
      if int(row[cancer]) != 1:
        continue
      contrast = dependency_contrast(effect, copy_number, models, lost, target, lineage)
      if contrast:
        details.append(contrast.__dict__)
        dep = dependency_components(contrast)
        components["dependency"] = max(components.get("dependency", 0), dep["dependency"])
        components["selectivity"] = max(components.get("selectivity", 0), dep["selectivity"])
      if cancer in tcga:
        freq = homozygous_deletion_frequency(tcga[cancer], lost)
        if freq:
          components["tumor_event"] = max(components.get("tumor_event", 0), event_component(freq[2]))
    if components:
      empirical[str(row["pair_id"])] = components
  result = score_candidate_table(candidates, empirical)
  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  result.to_csv(output, sep="\t", index=False)
  pd.DataFrame(details).to_csv(output.with_name("dependency_contrasts.tsv"), sep="\t", index=False)
  print(f"Wrote {len(result)} scored pairs to {output}")


if __name__ == "__main__":
  main()
