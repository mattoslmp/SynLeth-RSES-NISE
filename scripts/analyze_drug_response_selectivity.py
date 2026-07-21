#!/usr/bin/env python3
"""Test whether candidate drugs preferentially affect biomarker-loss cancer models.

The input sensitivity table is produced by ``standardize_drug_sensitivity.py``.
For each candidate loss-target-drug hypothesis, response in copy-number-loss and
intact models is compared within colorectal, gastric or lung lineages. P values
are Benjamini-Hochberg adjusted within each source and cancer family.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError
from scipy.stats import mannwhitneyu

from rses_onco.depmap import cancer_model_ids, read_depmap_inputs
from rses_onco.utils import bh_adjust, canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_COLUMNS = [
  "source", "cancer", "pair_id", "lost_gene", "target_gene", "drug_name",
  "drug_id", "drug_key", "pharmacology_source", "interaction_type",
  "action_type", "n_loss", "n_intact", "median_response_loss",
  "median_response_intact", "delta_response", "supportive_delta", "p_value",
  "q_value_bh", "lower_is_more_sensitive",
]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_optional(path: Path) -> pd.DataFrame:
  if not path.exists():
    return pd.DataFrame()
  try:
    return pd.read_csv(path, sep="\t")
  except EmptyDataError:
    return pd.DataFrame()


def drug_key(value: object) -> str:
  if value is None or pd.isna(value):
    return ""
  return re.sub(r"[^A-Z0-9]+", "", str(value).upper())


def candidate_drugs(evidence: pd.DataFrame) -> pd.DataFrame:
  columns = [
    "target_gene", "drug_name", "drug_id", "drug_key",
    "pharmacology_source", "interaction_type", "action_type",
  ]
  if evidence.empty:
    return pd.DataFrame(columns=columns)
  target_series = evidence.get(
    "target_gene", pd.Series(index=evidence.index, dtype=object)
  )
  drug_name_series = evidence.get(
    "drug_name", pd.Series(index=evidence.index, dtype=object)
  )
  drug_id_series = evidence.get(
    "drug_id", pd.Series(index=evidence.index, dtype=object)
  )
  rows = evidence.loc[
    target_series.notna() & (drug_name_series.notna() | drug_id_series.notna())
  ].copy()
  if rows.empty:
    return pd.DataFrame(columns=columns)
  rows["target_gene"] = rows["target_gene"].map(canonical_gene_name)
  rows["drug_key_name"] = rows.get(
    "drug_name", pd.Series(index=rows.index, dtype=object)
  ).map(drug_key)
  rows["drug_key_id"] = rows.get(
    "drug_id", pd.Series(index=rows.index, dtype=object)
  ).map(drug_key)
  expanded = []
  for record in rows.to_dict("records"):
    for key in {record.get("drug_key_name"), record.get("drug_key_id")}:
      if not key:
        continue
      expanded.append({
        "target_gene": record["target_gene"],
        "drug_name": record.get("drug_name"),
        "drug_id": record.get("drug_id"),
        "drug_key": key,
        "pharmacology_source": record.get("source"),
        "interaction_type": record.get("interaction_type"),
        "action_type": record.get("action_type"),
      })
  return pd.DataFrame(expanded, columns=columns).drop_duplicates()


def response_test(
  values: pd.DataFrame,
  loss_ids: set[str],
  intact_ids: set[str],
  lower_is_more_sensitive: bool,
  minimum_group_size: int,
) -> dict[str, Any] | None:
  loss = pd.to_numeric(
    values.loc[values["model_id"].astype(str).isin(loss_ids), "response_value"],
    errors="coerce",
  ).dropna()
  intact = pd.to_numeric(
    values.loc[values["model_id"].astype(str).isin(intact_ids), "response_value"],
    errors="coerce",
  ).dropna()
  if len(loss) < minimum_group_size or len(intact) < minimum_group_size:
    return None
  alternative = "less" if lower_is_more_sensitive else "greater"
  p_value = float(mannwhitneyu(loss, intact, alternative=alternative).pvalue)
  median_loss = float(np.median(loss))
  median_intact = float(np.median(intact))
  delta = median_loss - median_intact
  supportive_delta = -delta if lower_is_more_sensitive else delta
  return {
    "n_loss": len(loss),
    "n_intact": len(intact),
    "median_response_loss": median_loss,
    "median_response_intact": median_intact,
    "delta_response": delta,
    "supportive_delta": supportive_delta,
    "p_value": p_value,
    "lower_is_more_sensitive": lower_is_more_sensitive,
  }


def write_empty_outputs(output: Path, evidence_output: Path, message: str) -> None:
  output.parent.mkdir(parents=True, exist_ok=True)
  evidence_output.parent.mkdir(parents=True, exist_ok=True)
  empty = pd.DataFrame(columns=OUTPUT_COLUMNS)
  empty.to_csv(output, sep="\t", index=False)
  pd.DataFrame(columns=[*OUTPUT_COLUMNS, "raw_record"]).to_csv(
    evidence_output, sep="\t", index=False
  )
  print(message)
  print(f"Wrote {output}")
  print(f"Wrote {evidence_output}")


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--gene-effect", required=True)
  parser.add_argument("--copy-number", required=True)
  parser.add_argument("--models", required=True)
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--pharmacology-evidence",
    default="data/processed/pharmacology/pharmacology_evidence_long.tsv",
  )
  parser.add_argument(
    "--sensitivity",
    default="data/processed/pharmacology/drug_sensitivity_long.tsv",
  )
  parser.add_argument("--loss-threshold", type=float, default=0.30)
  parser.add_argument("--min-group-size", type=int, default=3)
  parser.add_argument(
    "--output",
    default="data/processed/pharmacology/drug_response_selectivity.tsv",
  )
  parser.add_argument(
    "--evidence-output",
    default="data/processed/pharmacology/pharmacology_sensitivity_evidence.tsv",
  )
  args = parser.parse_args()

  output = resolve_path(args.output)
  evidence_output = resolve_path(args.evidence_output)
  ranking = pd.read_csv(resolve_path(args.ranking), sep="\t")
  evidence = read_optional(resolve_path(args.pharmacology_evidence))
  sensitivity = read_optional(resolve_path(args.sensitivity))
  if sensitivity.empty:
    write_empty_outputs(
      output,
      evidence_output,
      "No standardized PRISM/GDSC/CTRP rows were available; wrote empty outputs.",
    )
    return
  drugs = candidate_drugs(evidence)
  if drugs.empty:
    write_empty_outputs(
      output,
      evidence_output,
      "No target-linked candidate compounds were available; wrote empty outputs.",
    )
    return

  _, copy_number, models, _ = read_depmap_inputs(
    resolve_path(args.gene_effect),
    resolve_path(args.copy_number),
    resolve_path(args.models),
    None,
  )
  drug_name_series = sensitivity.get(
    "drug_name", pd.Series(index=sensitivity.index, dtype=object)
  )
  drug_id_series = sensitivity.get(
    "drug_id", pd.Series(index=sensitivity.index, dtype=object)
  )
  sensitivity["drug_key_name"] = drug_name_series.map(drug_key)
  sensitivity["drug_key_id"] = drug_id_series.map(drug_key)
  sensitivity_long = pd.concat([
    sensitivity.assign(drug_key=sensitivity["drug_key_name"]),
    sensitivity.assign(drug_key=sensitivity["drug_key_id"]),
  ], ignore_index=True)
  sensitivity_long = sensitivity_long.loc[
    sensitivity_long["drug_key"].ne("")
  ].drop_duplicates()
  response_groups = {
    (str(source), str(key)): group
    for (source, key), group in sensitivity_long.groupby(["source", "drug_key"])
  }
  drug_map: dict[str, list[dict[str, Any]]] = {}
  for record in drugs.to_dict("records"):
    drug_map.setdefault(record["target_gene"], []).append(record)

  copy_number = copy_number.set_index("ModelID")
  rows: list[dict[str, Any]] = []
  sources = sorted(set(sensitivity_long["source"].astype(str)))
  for record in ranking.to_dict("records"):
    cancer = str(record.get("cancer") or "")
    if cancer not in {"colon", "stomach", "lung"}:
      continue
    lost = canonical_gene_name(
      record.get("analysis_lost_gene") or record.get("lost_gene")
    )
    target = canonical_gene_name(
      record.get("analysis_target_gene") or record.get("target_gene")
    )
    if not lost or not target or lost not in copy_number.columns:
      continue
    cohort_ids = set(cancer_model_ids(models, cancer))
    available_ids = sorted(cohort_ids & set(copy_number.index.astype(str)))
    if not available_ids:
      continue
    cn = pd.to_numeric(copy_number.loc[available_ids, lost], errors="coerce")
    loss_ids = set(cn.index[cn < args.loss_threshold].astype(str))
    intact_ids = set(cn.index[cn >= args.loss_threshold].astype(str))
    if len(loss_ids) < args.min_group_size or len(intact_ids) < args.min_group_size:
      continue
    for drug in drug_map.get(target, []):
      for source in sources:
        values = response_groups.get((source, drug["drug_key"]))
        if values is None or values.empty:
          continue
        lower_values = values["lower_is_more_sensitive"].dropna()
        lower_is_more_sensitive = (
          bool(lower_values.iloc[0]) if not lower_values.empty else True
        )
        tested = response_test(
          values,
          loss_ids,
          intact_ids,
          lower_is_more_sensitive,
          args.min_group_size,
        )
        if tested is None:
          continue
        rows.append({
          "source": source,
          "cancer": cancer,
          "pair_id": record.get("pair_id"),
          "lost_gene": lost,
          "target_gene": target,
          "drug_name": drug.get("drug_name"),
          "drug_id": drug.get("drug_id"),
          "drug_key": drug.get("drug_key"),
          "pharmacology_source": drug.get("pharmacology_source"),
          "interaction_type": drug.get("interaction_type"),
          "action_type": drug.get("action_type"),
          **tested,
        })

  result = pd.DataFrame(rows, columns=[
    column for column in OUTPUT_COLUMNS if column != "q_value_bh"
  ])
  if not result.empty:
    result["q_value_bh"] = np.nan
    for _, indices in result.groupby(["source", "cancer"]).groups.items():
      result.loc[indices, "q_value_bh"] = bh_adjust(
        result.loc[indices, "p_value"]
      )
    result = result[OUTPUT_COLUMNS].sort_values(
      ["q_value_bh", "supportive_delta"], ascending=[True, False]
    )
  else:
    result = pd.DataFrame(columns=OUTPUT_COLUMNS)
  output.parent.mkdir(parents=True, exist_ok=True)
  evidence_output.parent.mkdir(parents=True, exist_ok=True)
  result.to_csv(output, sep="\t", index=False)
  sensitivity_evidence = result.copy()
  sensitivity_evidence["raw_record"] = sensitivity_evidence.apply(
    lambda row: row.to_json(), axis=1
  ) if not sensitivity_evidence.empty else pd.Series(dtype=object)
  sensitivity_evidence.to_csv(evidence_output, sep="\t", index=False)
  print(f"Drug-response selectivity contrasts: {len(result):,}")
  if not result.empty:
    supported = (
      (result["q_value_bh"] < 0.05)
      & (result["supportive_delta"] > 0)
    )
    print(f"FDR < 0.05 and supportive direction: {int(supported.sum()):,}")
  print(f"Wrote {output}")
  print(f"Wrote {evidence_output}")


if __name__ == "__main__":
  main()
