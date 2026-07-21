#!/usr/bin/env python3
"""Run cancer-specific empirical RSES-Onco with current DepMap/GDC inputs."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from rses_onco.depmap import (
  dependency_components,
  dependency_contrast,
  expression_component,
  expression_contrast,
  read_depmap_inputs,
)
from rses_onco.integrate import load_reference_candidates, score_candidate_table
from rses_onco.tcga import event_component, homozygous_deletion_frequency, read_gistic_matrix
from rses_onco.utils import bh_adjust

LINEAGES = {"colon": "colon", "stomach": "stomach", "lung": "lung"}
ROOT = Path(__file__).resolve().parents[1]
SIMPLE_GENE = re.compile(r"^[A-Za-z0-9-]+$")
SIMPLE_LOSS_PATTERNS = [
  re.compile(r"^([A-Za-z0-9-]+)$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+loss$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+homozygous deletion$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+loss/low expression$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+loss or low expression$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+loss/low activity$", re.I),
]


def extract_single_lost_gene(feature: str) -> str | None:
  text = str(feature).strip()
  for pattern in SIMPLE_LOSS_PATTERNS:
    match = pattern.fullmatch(text)
    if match:
      return match.group(1).upper()
  return None


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--gene-effect", required=True)
  parser.add_argument("--copy-number", required=True)
  parser.add_argument("--models", required=True)
  parser.add_argument("--expression")
  parser.add_argument("--tcga", action="append", default=[], help="cancer=path/to/homdel_discrete.tsv")
  parser.add_argument("--candidates", default="data/curated/synthetic_lethality_reference_pairs.tsv")
  parser.add_argument("--output", default="results/empirical_rses_onco_by_cancer.tsv")
  parser.add_argument("--loss-threshold", type=float, default=0.30)
  parser.add_argument("--min-group-size", type=int, default=3)
  parser.add_argument("--expression-saturation-log2", type=float, default=2.0)
  args = parser.parse_args()

  effect, copy_number, models, expression = read_depmap_inputs(
    resolve_path(args.gene_effect),
    resolve_path(args.copy_number),
    resolve_path(args.models),
    resolve_path(args.expression) if args.expression else None,
  )
  tcga = {}
  for item in args.tcga:
    cancer, path = item.split("=", 1)
    tcga[cancer] = read_gistic_matrix(resolve_path(path))
  candidates = load_reference_candidates(resolve_path(args.candidates))

  scored_rows = []
  dependency_rows = []
  expression_rows = []
  skipped_rows = []

  for record in candidates.to_dict("records"):
    lost_feature = str(record["lost_feature"])
    lost = extract_single_lost_gene(lost_feature)
    target = str(record["target_gene"]).strip().upper()
    for cancer, lineage in LINEAGES.items():
      if int(record[cancer]) != 1:
        continue
      if lost is None or not SIMPLE_GENE.fullmatch(target):
        skipped_rows.append({
          "pair_id": record["pair_id"], "cancer": cancer, "lost_feature": lost_feature,
          "target_gene": target, "reason": "complex biomarker or multi-gene target requires explicit annotation",
        })
        continue
      components = {
        "tumor_event": None,
        "dependency": None,
        "selectivity": None,
        "expression_compensation": None,
      }
      expr = None
      contrast = dependency_contrast(
        effect, copy_number, models, lost, target, lineage,
        loss_threshold=args.loss_threshold, min_group_size=args.min_group_size,
      )
      if contrast:
        dependency_rows.append({"pair_id": record["pair_id"], "cancer": cancer, **contrast.__dict__})
        components.update(dependency_components(contrast))
      if expression is not None:
        expr = expression_contrast(
          expression, copy_number, models, lost, target, lineage,
          loss_threshold=args.loss_threshold, min_group_size=args.min_group_size,
        )
        if expr:
          expression_rows.append({"pair_id": record["pair_id"], "cancer": cancer, **expr.__dict__})
          components["expression_compensation"] = expression_component(
            expr, saturation_log2=args.expression_saturation_log2
          )
      if cancer in tcga:
        frequency = homozygous_deletion_frequency(tcga[cancer], lost)
        if frequency:
          components["tumor_event"] = event_component(frequency[2])
          components["tcga_homdel_n"] = frequency[0]
          components["tcga_evaluable_n"] = frequency[1]
          components["tcga_homdel_frequency"] = frequency[2]
      one = pd.DataFrame([record])
      empirical = {str(record["pair_id"]): {k: v for k, v in components.items() if not k.startswith("tcga_")}}
      scored = score_candidate_table(one, empirical).iloc[0].to_dict()
      scored.update({
        "cancer": cancer,
        "analysis_lost_gene": lost,
        "has_empirical_dependency": contrast is not None,
        "has_empirical_expression": expr is not None,
        "has_empirical_tcga": cancer in tcga and components.get("tumor_event") is not None,
        "score_basis": "empirical cohort components plus literature relation/validation priors",
        **{k: v for k, v in components.items() if k.startswith("tcga_")},
      })
      scored_rows.append(scored)

  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  result = pd.DataFrame(scored_rows)
  if not result.empty:
    result = result.sort_values(["cancer", "coverage_adjusted_rses"], ascending=[True, False])
  result.to_csv(output, sep="\t", index=False)

  dep = pd.DataFrame(dependency_rows)
  if not dep.empty:
    dep["q_value_bh"] = bh_adjust(dep["p_value"])
  dep.to_csv(output.with_name("dependency_contrasts.tsv"), sep="\t", index=False)

  expr = pd.DataFrame(expression_rows)
  if not expr.empty:
    expr["q_value_bh"] = bh_adjust(expr["p_value"])
  expr.to_csv(output.with_name("expression_contrasts.tsv"), sep="\t", index=False)
  pd.DataFrame(skipped_rows).to_csv(output.with_name("skipped_complex_biomarkers.tsv"), sep="\t", index=False)
  print(f"Wrote {len(result)} cancer-specific scored rows to {output}")


if __name__ == "__main__":
  main()
