#!/usr/bin/env python3
"""Discover cancer-specific conditional dependencies across all DepMap CRISPR targets.

For each analyzable loss gene, the script compares every measured CRISPR target
between copy-number-loss and intact models. Vectorized median effects are used as
a transparent prefilter; one-sided Mann-Whitney tests and global
Benjamini-Hochberg adjustment are then applied to retained targets.

This discovery layer can identify downstream, pathway-backup and other hidden
conditional dependencies not present in the curated benchmark. It does not assign
mechanistic class without independent relation evidence.
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
    source_class = str(record.get("source_class") or record.get("relation_type") or "curated")
    if not lost or not target:
      continue
    result.setdefault((lost, target), set()).add(source_class)
  return {key: ";".join(sorted(values)) for key, values in result.items()}


def test_one_target(
  effect: pd.DataFrame,
  loss_ids: list[str],
  intact_ids: list[str],
  target: str,
) -> tuple[int, int, float, float, float, float] | None:
  series = effect[["ModelID", target]].copy()
  series[target] = pd.to_numeric(series[target], errors="coerce")
  loss = series.loc[series["ModelID"].astype(str).isin(loss_ids), target].dropna().astype(float)
  intact = series.loc[series["ModelID"].astype(str).isin(intact_ids), target].dropna().astype(float)
  if len(loss) < 3 or len(intact) < 3:
    return None
  median_loss = float(np.median(loss))
  median_intact = float(np.median(intact))
  delta = median_loss - median_intact
  p_value = float(mannwhitneyu(loss, intact, alternative="less").pvalue)
  return len(loss), len(intact), median_loss, median_intact, delta, p_value


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
  parser.add_argument("--prefilter-targets", type=int, default=250)
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
  rows: list[dict[str, object]] = []

  for cancer in CANCERS:
    cohort_ids = set(cancer_model_ids(models, cancer))
    common_ids = sorted(cohort_ids & set(effect_by_id.index.astype(str)) & set(copy_by_id.index.astype(str)))
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
      median_loss = loss_matrix.median(axis=0, skipna=True)
      median_intact = intact_matrix.median(axis=0, skipna=True)
      delta = median_loss - median_intact
      prefilter = pd.DataFrame({
        "target_gene": delta.index,
        "median_effect_loss": median_loss.values,
        "median_effect_intact": median_intact.values,
        "delta_effect": delta.values,
      })
      prefilter = prefilter.loc[
        (prefilter["target_gene"].astype(str) != lost_gene)
        & (prefilter["delta_effect"] <= -abs(args.minimum_delta))
        & (prefilter["median_effect_loss"] <= args.maximum_median_loss_effect)
      ].sort_values("delta_effect").head(args.prefilter_targets)

      for target in prefilter["target_gene"].astype(str):
        tested = test_one_target(effect, loss_ids, intact_ids, target)
        if tested is None:
          continue
        n_loss, n_intact, med_loss, med_intact, effect_delta, p_value = tested
        target = canonical_gene_name(target)
        relation = class_map.get((lost_gene, target))
        rows.append({
          "cancer": cancer,
          "lost_gene": lost_gene,
          "target_gene": target,
          "n_loss": n_loss,
          "n_intact": n_intact,
          "median_effect_loss": med_loss,
          "median_effect_intact": med_intact,
          "delta_effect": effect_delta,
          "p_value": p_value,
          "known_relation_class": relation,
          "is_preexisting_candidate": relation is not None,
          "loss_universe": args.loss_universe,
        })
      print(
        f"[{cancer} {index}/{len(lost_genes)}] {lost_gene}: "
        f"loss={len(loss_ids)}, intact={len(intact_ids)}, prefiltered={len(prefilter)}",
        flush=True,
      )

  result = pd.DataFrame(rows)
  if not result.empty:
    result["q_value_bh"] = bh_adjust(result["p_value"])
    result = result.sort_values(["q_value_bh", "delta_effect"], ascending=[True, True])
  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  result.to_csv(output, sep="\t", index=False)

  if result.empty:
    discovered = pd.DataFrame()
  else:
    supported = result.loc[
      (result["q_value_bh"] < args.fdr)
      & (result["delta_effect"] < 0)
    ].copy()
    discovered_rows: list[dict[str, object]] = []
    for record in supported.to_dict("records"):
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
        "evidence_stage": "DepMap discovery screen",
        "primary_doi": "",
        "supporting_doi": "",
        "status": "unvalidated empirical discovery",
        "discovery_n_loss": record["n_loss"],
        "discovery_n_intact": record["n_intact"],
        "discovery_delta_effect": record["delta_effect"],
        "discovery_p_value": record["p_value"],
        "discovery_q_value_bh": record["q_value_bh"],
      })
    discovered = pd.DataFrame(discovered_rows).drop_duplicates(
      ["lost_gene", "target_gene", "source_class", "colon", "stomach", "lung"]
    )

  candidate_output = resolve_path(args.candidate_output)
  candidate_output.parent.mkdir(parents=True, exist_ok=True)
  discovered.to_csv(candidate_output, sep="\t", index=False)
  print(f"Tested/prefiltered contrasts: {len(result):,}")
  print(f"FDR-supported candidate rows: {len(discovered):,}")
  print(f"Wrote {output}")
  print(f"Wrote {candidate_output}")


if __name__ == "__main__":
  main()
