#!/usr/bin/env python3
"""Discover cancer-specific conditional dependencies across all DepMap CRISPR targets.

For each analyzable loss gene, every measured CRISPR target is tested between
copy-number-loss and intact models. One-sided Mann-Whitney tests are vectorized
across targets, and Benjamini-Hochberg correction is applied across the complete
target family for each loss-gene/cancer analysis before effect-size filtering.

This discovery layer can identify downstream, pathway-backup and other hidden
conditional dependencies not present in the curated benchmark. It does not assign
a mechanistic class without independent relation evidence.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from rses_onco.depmap import cancer_model_ids, read_depmap_inputs
from rses_onco.utils import bh_adjust, canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
CANCERS = ("colon", "stomach", "lung")
METADATA_COLUMNS = {
  "ModelID", "ProfileID", "PROFILEID", "is_default_entry",
  "IsDefaultEntryForModel", "IsDefaultEntryForMC", "ModelConditionID",
  "SequencingID",
}


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def candidate_loss_genes(candidates: pd.DataFrame) -> list[str]:
  if "lost_gene" not in candidates:
    raise ValueError("Candidate universe lacks lost_gene")
  return sorted({
    canonical_gene_name(value)
    for value in candidates["lost_gene"].dropna()
    if canonical_gene_name(value)
  })


def all_copy_number_loss_genes(
  copy_number: pd.DataFrame,
  model_ids: set[str],
  threshold: float,
  minimum_loss_models: int,
  maximum_genes: int | None,
) -> list[str]:
  genes = [column for column in copy_number.columns if column not in METADATA_COLUMNS]
  subset = copy_number.loc[
    copy_number["ModelID"].astype(str).isin(model_ids),
    genes,
  ].apply(pd.to_numeric, errors="coerce")
  counts = subset.lt(threshold).sum(axis=0)
  eligible = counts.loc[counts >= minimum_loss_models].sort_values(ascending=False)
  if maximum_genes is not None:
    eligible = eligible.head(maximum_genes)
  return [canonical_gene_name(gene) for gene in eligible.index]


def known_class_map(candidates: pd.DataFrame) -> dict[tuple[str, str], str]:
  result: dict[tuple[str, str], set[str]] = {}
  for record in candidates.to_dict("records"):
    lost = canonical_gene_name(record.get("lost_gene"))
    target = canonical_gene_name(record.get("target_gene"))
    raw_class = record.get("source_class")
    if raw_class is None or (isinstance(raw_class, float) and not np.isfinite(raw_class)):
      raw_class = record.get("relation_type") or "curated"
    source_class = str(raw_class)
    if not lost or not target:
      continue
    result.setdefault((lost, target), set()).add(source_class)
  return {key: ";".join(sorted(values)) for key, values in result.items()}


def vectorized_target_tests(
  loss_matrix: pd.DataFrame,
  intact_matrix: pd.DataFrame,
  minimum_group_size: int,
) -> pd.DataFrame:
  """Test every target and adjust across the whole target family.

  The target family is defined by one loss gene in one cancer cohort. Effect-size
  filters are applied only after P values and within-family BH Q values have been
  computed, avoiding a data-dependent P-value prefilter.
  """
  targets = list(loss_matrix.columns)
  loss_values = loss_matrix.to_numpy(dtype=float)
  intact_values = intact_matrix.to_numpy(dtype=float)
  n_loss = np.sum(np.isfinite(loss_values), axis=0)
  n_intact = np.sum(np.isfinite(intact_values), axis=0)
  eligible = (n_loss >= minimum_group_size) & (n_intact >= minimum_group_size)
  median_loss = np.nanmedian(loss_values, axis=0)
  median_intact = np.nanmedian(intact_values, axis=0)
  delta = median_loss - median_intact
  p_values = np.full(len(targets), np.nan, dtype=float)
  if eligible.any():
    tested = mannwhitneyu(
      loss_values[:, eligible],
      intact_values[:, eligible],
      alternative="less",
      axis=0,
      method="auto",
      nan_policy="omit",
    )
    p_values[eligible] = np.asarray(tested.pvalue, dtype=float)
  q_values = bh_adjust(p_values)
  return pd.DataFrame({
    "target_gene": [canonical_gene_name(target) for target in targets],
    "n_loss": n_loss,
    "n_intact": n_intact,
    "median_effect_loss": median_loss,
    "median_effect_intact": median_intact,
    "delta_effect": delta,
    "p_value": p_values,
    "q_value_bh_within_loss_cancer": q_values,
    "target_family_size": int(eligible.sum()),
  })


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--gene-effect", required=True)
  parser.add_argument("--copy-number", required=True)
  parser.add_argument("--models", required=True)
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--loss-universe",
    choices=("candidates", "all-cn"),
    default="candidates",
  )
  parser.add_argument("--loss-threshold", type=float, default=0.30)
  parser.add_argument("--min-group-size", type=int, default=3)
  parser.add_argument("--minimum-delta", type=float, default=0.15)
  parser.add_argument("--maximum-median-loss-effect", type=float, default=-0.25)
  parser.add_argument(
    "--max-loss-genes",
    type=int,
    default=None,
    help="Optional safety cap for --loss-universe all-cn; omit for all eligible CN genes.",
  )
  parser.add_argument("--fdr", type=float, default=0.10)
  parser.add_argument(
    "--output",
    default="results/expanded_26Q1/discovery/all_target_dependency_screen.tsv",
  )
  parser.add_argument(
    "--candidate-output",
    default="data/raw/discovery/depmap_discovered_candidate_pairs.tsv",
  )
  args = parser.parse_args()

  effect, copy_number, models, _ = read_depmap_inputs(
    resolve_path(args.gene_effect),
    resolve_path(args.copy_number),
    resolve_path(args.models),
    None,
  )
  candidates = pd.read_csv(resolve_path(args.candidates), sep="\t")
  class_map = known_class_map(candidates)
  target_genes = [
    column for column in effect.columns
    if column not in METADATA_COLUMNS
  ]
  effect_by_id = effect.set_index("ModelID")
  copy_by_id = copy_number.set_index("ModelID")
  rows: list[pd.DataFrame] = []

  for cancer in CANCERS:
    cohort_ids = set(cancer_model_ids(models, cancer))
    common_ids = sorted(
      cohort_ids
      & set(effect_by_id.index.astype(str))
      & set(copy_by_id.index.astype(str))
    )
    if args.loss_universe == "candidates":
      lost_genes = candidate_loss_genes(candidates)
    else:
      lost_genes = all_copy_number_loss_genes(
        copy_number,
        set(common_ids),
        args.loss_threshold,
        args.min_group_size,
        args.max_loss_genes,
      )

    for index, lost_gene in enumerate(lost_genes, start=1):
      if lost_gene not in copy_by_id.columns:
        continue
      cn = pd.to_numeric(copy_by_id.loc[common_ids, lost_gene], errors="coerce")
      loss_ids = cn.index[cn < args.loss_threshold].astype(str).tolist()
      intact_ids = cn.index[cn >= args.loss_threshold].astype(str).tolist()
      if len(loss_ids) < args.min_group_size or len(intact_ids) < args.min_group_size:
        continue

      loss_matrix = effect_by_id.loc[
        effect_by_id.index.astype(str).isin(loss_ids), target_genes
      ].apply(pd.to_numeric, errors="coerce")
      intact_matrix = effect_by_id.loc[
        effect_by_id.index.astype(str).isin(intact_ids), target_genes
      ].apply(pd.to_numeric, errors="coerce")
      tested = vectorized_target_tests(
        loss_matrix,
        intact_matrix,
        args.min_group_size,
      )
      tested = tested.loc[
        (tested["target_gene"] != lost_gene)
        & (tested["q_value_bh_within_loss_cancer"] < args.fdr)
        & (tested["delta_effect"] <= -abs(args.minimum_delta))
        & (tested["median_effect_loss"] <= args.maximum_median_loss_effect)
      ].copy()
      if not tested.empty:
        tested.insert(0, "lost_gene", lost_gene)
        tested.insert(0, "cancer", cancer)
        tested["known_relation_class"] = [
          class_map.get((lost_gene, target))
          for target in tested["target_gene"]
        ]
        tested["is_preexisting_candidate"] = tested["known_relation_class"].notna()
        tested["loss_universe"] = args.loss_universe
        rows.append(tested)
      print(
        f"[{cancer} {index}/{len(lost_genes)}] {lost_gene}: "
        f"loss={len(loss_ids)}, intact={len(intact_ids)}, "
        f"all_targets={len(target_genes)}, supported={len(tested)}",
        flush=True,
      )

  result = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
  if not result.empty:
    result = result.sort_values(
      ["q_value_bh_within_loss_cancer", "delta_effect"],
      ascending=[True, True],
    )
  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  result.to_csv(output, sep="\t", index=False)

  discovered_rows: list[dict[str, object]] = []
  for record in result.to_dict("records"):
    cancer = str(record["cancer"])
    lost = canonical_gene_name(record["lost_gene"])
    target = canonical_gene_name(record["target_gene"])
    known_class = record.get("known_relation_class")
    source_class = (
      str(known_class)
      if known_class is not None and pd.notna(known_class)
      else "empirical_conditional_dependency"
    )
    discovered_rows.append({
      "pair_id": f"DEPMAP_DISC_{cancer}_{lost}_TO_{target}",
      "lost_feature": f"{lost} loss",
      "lost_gene": lost,
      "target_gene": target,
      "source_class": source_class,
      "relation_type": "depmap_all_target_conditional_dependency",
      "mechanism": (
        f"DepMap {cancer} models with {lost} copy-number loss show stronger "
        f"dependency on {target}; mechanism requires independent network and experimental annotation."
      ),
      "colon": int(cancer == "colon"),
      "stomach": int(cancer == "stomach"),
      "lung": int(cancer == "lung"),
      "relation_confidence": 0.0,
      "genetic_screen": 0.0,
      "isogenic_validation": 0.0,
      "in_vivo": 0.0,
      "clinical_tractability": 0.0,
      "lineage_relevance": 0.0,
      "evidence_stage": "DepMap all-target discovery screen",
      "primary_doi": "",
      "supporting_doi": "",
      "status": "unvalidated empirical discovery",
      "discovery_n_loss": record["n_loss"],
      "discovery_n_intact": record["n_intact"],
      "discovery_delta_effect": record["delta_effect"],
      "discovery_p_value": record["p_value"],
      "discovery_q_value_bh_within_loss_cancer": record[
        "q_value_bh_within_loss_cancer"
      ],
      "discovery_target_family_size": record["target_family_size"],
    })
  discovered = pd.DataFrame(discovered_rows)
  if not discovered.empty:
    discovered = discovered.drop_duplicates(
      ["lost_gene", "target_gene", "source_class", "colon", "stomach", "lung"]
    )

  candidate_output = resolve_path(args.candidate_output)
  candidate_output.parent.mkdir(parents=True, exist_ok=True)
  discovered.to_csv(candidate_output, sep="\t", index=False)
  print(f"FDR- and effect-supported all-target contrasts: {len(result):,}")
  print(f"Standardized discovered candidate rows: {len(discovered):,}")
  print(f"Wrote {output}")
  print(f"Wrote {candidate_output}")


if __name__ == "__main__":
  main()
