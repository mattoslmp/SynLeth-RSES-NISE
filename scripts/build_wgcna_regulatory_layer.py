#!/usr/bin/env python3
"""Build cancer-specific WGCNA, TF-expression and promoter-motif evidence.

WGCNA and pairwise coexpression are kept inside one expression-network domain to
avoid double counting the same expression matrix. DoRothEA TF-target associations,
TF-target expression consistency and JASPAR promoter motif predictions are combined
inside one regulatory domain. Motif predictions are never called direct binding.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Mapping

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from rses_onco.depmap import cancer_model_ids, read_depmap_inputs
from rses_onco.utils import canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
CANCERS = ("colon", "stomach", "lung")
WGCNA_SUBWEIGHTS = {
  "tom_divergence": 0.40,
  "module_divergence": 0.30,
  "kme_divergence": 0.30,
}
REGULATORY_SUBWEIGHTS = {
  "tf_association_divergence": 0.40,
  "tf_expression_profile_divergence": 0.35,
  "promoter_motif_divergence": 0.25,
}


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def clean_bool(value: object) -> bool | None:
  if value is None or (isinstance(value, float) and not math.isfinite(value)):
    return None
  text = str(value).strip().casefold()
  if text in {"1", "true", "yes", "y"}:
    return True
  if text in {"0", "false", "no", "n", "", "nan", "none"}:
    return False
  return None


def numeric(value: object) -> float | None:
  try:
    result = float(value)
  except (TypeError, ValueError):
    return None
  return result if math.isfinite(result) else None


def weighted_subscore(
  components: Mapping[str, float | None],
  weights: Mapping[str, float],
) -> dict[str, float | int]:
  total = float(sum(weights.values()))
  numerator = 0.0
  observed_weight = 0.0
  observed = 0
  for name, weight in weights.items():
    value = numeric(components.get(name))
    if value is None:
      continue
    clipped = float(np.clip(value, 0, 1))
    numerator += weight * clipped
    observed_weight += weight
    observed += 1
  raw = numerator / observed_weight if observed_weight else float("nan")
  coverage = observed_weight / total if total else float("nan")
  adjusted = raw * coverage if math.isfinite(raw) and math.isfinite(coverage) else float("nan")
  return {
    "raw": raw,
    "coverage": coverage,
    "adjusted": adjusted,
    "observed_subcomponents": observed,
  }


def jaccard_divergence(first: set[str], second: set[str]) -> float | None:
  union = first | second
  if not union:
    return None
  return 1.0 - len(first & second) / len(union)


def cosine_divergence(first: dict[str, float], second: dict[str, float]) -> float | None:
  keys = sorted(set(first) | set(second))
  if not keys:
    return None
  vector_a = np.array([first.get(key, 0.0) for key in keys], dtype=float)
  vector_b = np.array([second.get(key, 0.0) for key in keys], dtype=float)
  denominator = float(np.linalg.norm(vector_a) * np.linalg.norm(vector_b))
  if denominator == 0:
    return None
  cosine = float(np.clip(np.dot(vector_a, vector_b) / denominator, -1, 1))
  return (1.0 - cosine) / 2.0


def first_column(frame: pd.DataFrame, names: tuple[str, ...]) -> str | None:
  return next((name for name in names if name in frame), None)


def prepare_dorothea(frame: pd.DataFrame) -> pd.DataFrame:
  if frame.empty:
    return pd.DataFrame(columns=["tf", "target", "expected_sign"])
  source = first_column(frame, ("source_genesymbol", "source", "tf"))
  target = first_column(frame, ("target_genesymbol", "target", "gene"))
  if source is None or target is None:
    return pd.DataFrame(columns=["tf", "target", "expected_sign"])
  result = frame.copy()
  result["tf"] = result[source].map(canonical_gene_name)
  result["target"] = result[target].map(canonical_gene_name)
  stimulation_column = first_column(
    result,
    ("consensus_stimulation", "is_stimulation", "stimulation"),
  )
  inhibition_column = first_column(
    result,
    ("consensus_inhibition", "is_inhibition", "inhibition"),
  )
  stimulation = (
    result[stimulation_column].map(clean_bool)
    if stimulation_column
    else pd.Series(None, index=result.index)
  )
  inhibition = (
    result[inhibition_column].map(clean_bool)
    if inhibition_column
    else pd.Series(None, index=result.index)
  )
  result["expected_sign"] = np.where(
    stimulation.eq(True) & ~inhibition.eq(True),
    1,
    np.where(inhibition.eq(True) & ~stimulation.eq(True), -1, 0),
  )
  return result.loc[
    result["tf"].ne("") & result["target"].ne(""),
    ["tf", "target", "expected_sign"],
  ].drop_duplicates()


def expression_matrix_for_cancer(
  expression: pd.DataFrame,
  models: pd.DataFrame,
  cancer: str,
  candidate_genes: set[str],
  max_genes: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
  selected = set(cancer_model_ids(models, cancer).astype(str))
  matrix = expression.loc[expression["ModelID"].astype(str).isin(selected)].copy()
  gene_columns = [column for column in matrix.columns if column != "ModelID"]
  numeric_frame = matrix[gene_columns].apply(pd.to_numeric, errors="coerce")
  observed_fraction = numeric_frame.notna().mean()
  variance = numeric_frame.var(axis=0, skipna=True)
  eligible = observed_fraction.ge(0.80) & variance.gt(0)
  eligible_genes = set(eligible.index[eligible])
  candidate_present = sorted(candidate_genes & eligible_genes)
  cap = min(max_genes, max(500, len(matrix) * 50))
  variable = (
    variance.loc[list(eligible_genes)]
      .sort_values(ascending=False)
      .head(cap)
      .index
      .tolist()
  )
  selected_genes = sorted(set(variable) | set(candidate_present))
  result = pd.concat(
    [matrix[["ModelID"]].reset_index(drop=True), numeric_frame[selected_genes].reset_index(drop=True)],
    axis=1,
  )
  imputed = 0
  for column in selected_genes:
    missing = int(result[column].isna().sum())
    if missing:
      result[column] = result[column].fillna(result[column].median())
      imputed += missing
  diagnostics = {
    "cancer": cancer,
    "samples": len(result),
    "selected_genes": len(selected_genes),
    "candidate_genes_present": len(candidate_present),
    "median_imputed_cells": imputed,
    "input_expression_scale": "DepMap log2(TPM+1)",
    "gene_selection": "all candidate genes plus cancer-specific highest-variance genes",
  }
  return result, diagnostics


def run_wgcna(
  expression: pd.DataFrame,
  pairs: pd.DataFrame,
  cancer: str,
  work_dir: Path,
  rscript: str,
  script: Path,
) -> pd.DataFrame:
  cancer_dir = work_dir / cancer
  cancer_dir.mkdir(parents=True, exist_ok=True)
  expression_path = cancer_dir / "wgcna_expression_input.tsv"
  pair_path = cancer_dir / "wgcna_candidate_pairs.tsv"
  expression.to_csv(expression_path, sep="\t", index=False)
  pairs.to_csv(pair_path, sep="\t", index=False)
  subprocess.run(
    [rscript, str(script), str(expression_path), str(pair_path), str(cancer_dir), cancer],
    cwd=ROOT,
    check=True,
  )
  path = cancer_dir / "wgcna_pair_metrics.tsv"
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(path)
  return pd.read_csv(path, sep="\t", low_memory=False)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--gene-effect", required=True)
  parser.add_argument("--copy-number", required=True)
  parser.add_argument("--models", required=True)
  parser.add_argument("--expression", required=True)
  parser.add_argument("--candidates", default="data/processed/expanded_candidate_universe.tsv")
  parser.add_argument("--functional-evidence", default="data/processed/expanded_pair_functional_evidence.tsv")
  parser.add_argument("--dorothea", default="data/raw/human_functional_evidence/omnipath_dorothea.tsv")
  parser.add_argument("--promoter-motifs", default="data/processed/regulatory/jaspar_promoter_tf_summary.tsv")
  parser.add_argument("--output", default="data/processed/expanded_pair_functional_evidence.tsv")
  parser.add_argument("--work-dir", default="data/processed/regulatory/wgcna")
  parser.add_argument("--rscript", default="Rscript")
  parser.add_argument("--wgcna-script", default="scripts/run_wgcna_expression_network.R")
  parser.add_argument("--max-genes", type=int, default=3000)
  args = parser.parse_args()

  rscript = shutil.which(args.rscript)
  if rscript is None:
    raise RuntimeError("Rscript is required for the WGCNA regulatory layer")
  subprocess.run(
    [
      rscript,
      "-e",
      "stopifnot(requireNamespace('WGCNA', quietly=TRUE)); "
      "cat(as.character(packageVersion('WGCNA')))" ,
    ],
    check=True,
  )

  effect, copy_number, models, expression = read_depmap_inputs(
    resolve_path(args.gene_effect),
    resolve_path(args.copy_number),
    resolve_path(args.models),
    resolve_path(args.expression),
  )
  if expression is None:
    raise ValueError("Expression matrix is required for WGCNA")
  candidates = pd.read_csv(resolve_path(args.candidates), sep="\t", low_memory=False)
  base = pd.read_csv(resolve_path(args.functional_evidence), sep="\t", low_memory=False)
  pairs = candidates[[column for column in (
    "pair_id", "lost_gene", "target_gene", "source_class"
  ) if column in candidates]].copy()
  pairs["lost_gene"] = pairs["lost_gene"].map(canonical_gene_name)
  pairs["target_gene"] = pairs["target_gene"].map(canonical_gene_name)
  pairs = pairs.loc[pairs["lost_gene"].ne("") & pairs["target_gene"].ne("")].drop_duplicates("pair_id")
  candidate_genes = set(pairs["lost_gene"]) | set(pairs["target_gene"])

  work_dir = resolve_path(args.work_dir)
  work_dir.mkdir(parents=True, exist_ok=True)
  wgcna_frames = []
  prep_diagnostics = []
  for cancer in CANCERS:
    matrix, diagnostics = expression_matrix_for_cancer(
      expression,
      models,
      cancer,
      candidate_genes,
      args.max_genes,
    )
    if len(matrix) < 20:
      raise RuntimeError(
        f"Cancer-specific WGCNA requires at least 20 models; {cancer} has {len(matrix)}"
      )
    wgcna_frames.append(
      run_wgcna(
        matrix,
        pairs[["pair_id", "lost_gene", "target_gene"]],
        cancer,
        work_dir,
        rscript,
        resolve_path(args.wgcna_script),
      )
    )
    prep_diagnostics.append(diagnostics)
  wgcna = pd.concat(wgcna_frames, ignore_index=True)
  wgcna_components = []
  for record in wgcna.to_dict("records"):
    score = weighted_subscore(
      {
        "tom_divergence": record.get("wgcna_tom_divergence"),
        "module_divergence": record.get("wgcna_module_divergence"),
        "kme_divergence": record.get("wgcna_kME_divergence"),
      },
      WGCNA_SUBWEIGHTS,
    )
    wgcna_components.append({
      **record,
      "wgcna_expression_network_raw": score["raw"],
      "wgcna_expression_network_coverage": score["coverage"],
      "component_wgcna_expression_network": score["adjusted"],
      "wgcna_observed_subcomponents": score["observed_subcomponents"],
    })
  wgcna = pd.DataFrame(wgcna_components)
  wgcna.to_csv(work_dir / "wgcna_pair_metrics_all_cancers.tsv", sep="\t", index=False)
  pd.DataFrame(prep_diagnostics).to_csv(
    work_dir / "wgcna_input_preparation.tsv",
    sep="\t",
    index=False,
  )

  dorothea_path = resolve_path(args.dorothea)
  dorothea = prepare_dorothea(
    pd.read_csv(dorothea_path, sep="\t", low_memory=False)
    if dorothea_path.exists() and dorothea_path.stat().st_size
    else pd.DataFrame()
  )
  dorothea = dorothea.loc[dorothea["target"].isin(candidate_genes)].copy()
  regulator_sets = {
    gene: set(group["tf"].astype(str))
    for gene, group in dorothea.groupby("target")
  }
  all_regulators = set(dorothea["tf"].astype(str))

  motif_path = resolve_path(args.promoter_motifs)
  motif = (
    pd.read_csv(motif_path, sep="\t", low_memory=False)
    if motif_path.exists() and motif_path.stat().st_size
    else pd.DataFrame()
  )
  if not motif.empty:
    motif["gene"] = motif["gene"].map(canonical_gene_name)
    motif["transcription_factor"] = motif["transcription_factor"].map(canonical_gene_name)
    motif = motif.loc[motif["transcription_factor"].isin(all_regulators)].copy()
  motif_sets = {
    gene: set(group["transcription_factor"].astype(str))
    for gene, group in motif.groupby("gene")
  } if not motif.empty else {}

  regulatory_rows = []
  for cancer in CANCERS:
    selected = set(cancer_model_ids(models, cancer).astype(str))
    cancer_expression = expression.loc[
      expression["ModelID"].astype(str).isin(selected)
    ].copy()
    genes_for_correlation = sorted(
      (all_regulators | candidate_genes) & set(cancer_expression.columns)
    )
    correlation = (
      cancer_expression[genes_for_correlation]
        .apply(pd.to_numeric, errors="coerce")
        .corr(method="spearman", min_periods=3)
    )
    profiles: dict[str, dict[str, float]] = {}
    for record in dorothea.to_dict("records"):
      tf = str(record["tf"])
      target = str(record["target"])
      if tf not in correlation.index or target not in correlation.columns:
        continue
      rho = numeric(correlation.at[tf, target])
      if rho is None:
        continue
      expected_sign = int(record.get("expected_sign") or 0)
      support = rho * expected_sign if expected_sign else rho
      profiles.setdefault(target, {})[tf] = float(np.clip(support, -1, 1))

    for record in pairs.to_dict("records"):
      lost = str(record["lost_gene"])
      target = str(record["target_gene"])
      association = jaccard_divergence(
        regulator_sets.get(lost, set()),
        regulator_sets.get(target, set()),
      )
      expression_profile = cosine_divergence(
        profiles.get(lost, {}),
        profiles.get(target, {}),
      )
      promoter = jaccard_divergence(
        motif_sets.get(lost, set()),
        motif_sets.get(target, set()),
      )
      score = weighted_subscore(
        {
          "tf_association_divergence": association,
          "tf_expression_profile_divergence": expression_profile,
          "promoter_motif_divergence": promoter,
        },
        REGULATORY_SUBWEIGHTS,
      )
      regulatory_rows.append({
        "cancer": cancer,
        "pair_id": record["pair_id"],
        "regulatory_tf_association_divergence": association,
        "regulatory_tf_expression_profile_divergence": expression_profile,
        "regulatory_promoter_motif_divergence": promoter,
        "regulatory_lost_regulator_count": len(regulator_sets.get(lost, set())),
        "regulatory_target_regulator_count": len(regulator_sets.get(target, set())),
        "regulatory_lost_promoter_tf_count": len(motif_sets.get(lost, set())),
        "regulatory_target_promoter_tf_count": len(motif_sets.get(target, set())),
        "regulatory_network_raw": score["raw"],
        "regulatory_network_coverage": score["coverage"],
        "component_regulatory_network_composite": score["adjusted"],
        "regulatory_observed_subcomponents": score["observed_subcomponents"],
        "promoter_evidence_type": "JASPAR_motif_prediction_not_direct_binding",
      })
  regulatory = pd.DataFrame(regulatory_rows)
  regulatory_dir = work_dir.parent
  regulatory.to_csv(
    regulatory_dir / "promoter_tf_regulatory_pair_metrics.tsv",
    sep="\t",
    index=False,
  )

  if "component_regulatory_network" in base:
    base = base.rename(columns={
      "component_regulatory_network": "component_regulatory_network_dorothea_pair_level"
    })
  if "cancer" in base:
    expanded = base.copy()
  else:
    expanded = pd.concat(
      [base.assign(cancer=cancer) for cancer in CANCERS],
      ignore_index=True,
    )
  enriched = (
    expanded.merge(wgcna, on=["pair_id", "cancer"], how="left", suffixes=("", "_wgcna"))
      .merge(regulatory, on=["pair_id", "cancer"], how="left")
  )
  enriched["component_regulatory_network"] = enriched[
    "component_regulatory_network_composite"
  ]
  enriched["expression_network_method"] = "pairwise_Spearman_plus_signed_WGCNA"
  enriched["regulatory_network_method"] = (
    "DoRothEA_TF_target_plus_TF_expression_consistency_plus_"
    "JASPAR_promoter_motif_prediction"
  )
  enriched["regulatory_layer_version"] = "wgcna-promoter-regulatory-v1"
  enriched["regulatory_layer_generated_at_utc"] = datetime.now(timezone.utc).isoformat()

  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  temporary = output.with_suffix(output.suffix + ".tmp")
  enriched.to_csv(temporary, sep="\t", index=False)
  temporary.replace(output)
  status = {
    "version": "wgcna-promoter-regulatory-v1",
    "output_rows": len(enriched),
    "pair_count": int(enriched["pair_id"].nunique()),
    "cancers": list(CANCERS),
    "wgcna_subweights": WGCNA_SUBWEIGHTS,
    "regulatory_subweights": REGULATORY_SUBWEIGHTS,
    "direct_promoter_binding_claim": False,
    "motif_interpretation": "predicted promoter motif occurrence only",
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
  }
  (regulatory_dir / "wgcna_regulatory_layer_status.json").write_text(
    json.dumps(status, indent=2, sort_keys=True),
    encoding="utf-8",
  )
  print(f"WGCNA pair-cancer rows: {len(wgcna):,}")
  print(f"Regulatory pair-cancer rows: {len(regulatory):,}")
  print(f"Wrote enriched functional evidence: {output} ({len(enriched):,} rows)")


if __name__ == "__main__":
  main()
