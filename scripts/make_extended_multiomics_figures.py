#!/usr/bin/env python3
"""Generate Supplementary Figures S71-S78 for the extended multi-omics layer."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
  sys.path.insert(0, str(ROOT / "src"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from rses_onco.publication import (
  figure_record,
  placeholder,
  save_figure_triplet,
  set_publication_style,
  write_figure_manifest,
  write_legends_markdown,
  write_source_data,
)

SCRIPT = "scripts/make_extended_multiomics_figures.py"


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_optional(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    return pd.DataFrame()
  return pd.read_csv(path, sep="\t", low_memory=False)


def save_asset(
  fig: plt.Figure,
  item: dict[str, object],
  source: pd.DataFrame,
  inputs: list[Path],
  output_root: Path,
  strict: bool,
):
  figure_id = str(item["id"])
  file_name = str(item["file"])
  base = output_root / "figures/supplementary" / file_name
  source_path = (
    output_root
    / "source_data/figures/supplementary"
    / f"{file_name}_source_data.tsv"
  )
  write_source_data(source, source_path)
  audit = save_figure_triplet(
    fig,
    base,
    figure_id,
    strict_layout=strict,
  )
  return figure_record(
    figure_id=figure_id,
    category="supplementary",
    title=str(item["title"]),
    caption=str(item.get("caption") or item["title"]),
    base_path=base,
    source_data_path=source_path,
    input_paths=inputs,
    audit=audit,
    script=SCRIPT,
  )


def figure_source_coverage(
  inventory: pd.DataFrame,
) -> tuple[plt.Figure, pd.DataFrame]:
  if inventory.empty:
    fig, axis = plt.subplots(figsize=(9.5, 5.5), constrained_layout=True)
    placeholder(
      axis,
      "Extended multi-omics source coverage",
      "No source inventory was available.",
    )
    return fig, inventory
  data = inventory.copy()
  data["role_group"] = np.select(
    [
      data.get("direct_score_layer", False).astype(bool),
      data.get("exploratory_or_validation_layer", False).astype(bool),
    ],
    ["Scored causal/orthogonal", "Exploratory or validation"],
    default="Translation or existing layer",
  )
  summary = (
    data.groupby("role_group", dropna=False)
      .agg(
        declared_sources=("source_key", "size"),
        available_sources=("exists", "sum"),
      )
      .reset_index()
  )
  summary["availability"] = (
    summary["available_sources"]
    / summary["declared_sources"].clip(lower=1)
  )
  fig, axis = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
  x = np.arange(len(summary))
  axis.bar(x, summary["declared_sources"], label="Declared")
  axis.bar(x, summary["available_sources"], label="Available")
  axis.set_xticks(
    x,
    [str(value).replace(" ", "\n") for value in summary["role_group"]],
  )
  axis.set_ylabel("Number of sources")
  axis.set_title("Extended multi-omics source availability and declared use")
  axis.legend(frameon=False)
  axis.grid(axis="y", alpha=0.25)
  return fig, data


def figure_baseline_extended(
  ranking: pd.DataFrame,
) -> tuple[plt.Figure, pd.DataFrame]:
  columns = [
    column
    for column in (
      "pair_id",
      "cancer",
      "source_class",
      "baseline_coverage_adjusted_rses",
      "coverage_adjusted_rses",
      "extended_score_delta",
      "baseline_rank_within_cancer",
      "extended_rank_within_cancer",
    )
    if column in ranking.columns
  ]
  data = ranking[columns].copy() if columns else pd.DataFrame()
  fig, axis = plt.subplots(figsize=(7.2, 6.4), constrained_layout=True)
  required = {"baseline_coverage_adjusted_rses", "coverage_adjusted_rses"}
  if data.empty or not required.issubset(data.columns):
    placeholder(
      axis,
      "Baseline versus extended RSES-Onco",
      "No paired baseline and extended scores were available.",
    )
    return fig, data
  for cancer, subset in data.groupby("cancer"):
    axis.scatter(
      pd.to_numeric(
        subset["baseline_coverage_adjusted_rses"], errors="coerce"
      ),
      pd.to_numeric(subset["coverage_adjusted_rses"], errors="coerce"),
      s=18,
      alpha=0.55,
      label=str(cancer),
    )
  axis.plot([0, 1], [0, 1], linestyle="--", linewidth=1.0)
  axis.set_xlim(0, 1)
  axis.set_ylim(0, 1)
  axis.set_xlabel("Baseline coverage-adjusted RSES-Onco")
  axis.set_ylabel("Extended coverage-adjusted RSES-Onco")
  axis.set_title("Baseline and causal multi-omics extended scores")
  axis.legend(frameon=False)
  axis.grid(alpha=0.2)
  return fig, data


def figure_layer_heatmap(
  evidence: pd.DataFrame,
  ranking: pd.DataFrame,
) -> tuple[plt.Figure, pd.DataFrame]:
  layers = [
    "integrated_functional_loss_support",
    "dependency_probability_support",
    "protein_compensation_support",
    "rnai_orthogonal_support",
  ]
  columns = [
    "pair_id",
    "cancer",
    *[layer for layer in layers if layer in evidence.columns],
  ]
  data = evidence[columns].copy() if columns else pd.DataFrame()
  if not data.empty and "coverage_adjusted_rses" in ranking.columns:
    data = data.merge(
      ranking[["pair_id", "cancer", "coverage_adjusted_rses"]],
      on=["pair_id", "cancer"],
      how="left",
    )
    data = data.sort_values(
      "coverage_adjusted_rses", ascending=False
    ).head(25)
  fig, axis = plt.subplots(figsize=(10.8, 8.4), constrained_layout=True)
  observed_layers = [layer for layer in layers if layer in data.columns]
  if data.empty or not observed_layers:
    placeholder(
      axis,
      "Extended layer contributions",
      "No scored extended-layer evidence was available.",
    )
    return fig, data
  matrix = (
    data[observed_layers]
      .apply(pd.to_numeric, errors="coerce")
      .to_numpy()
  )
  masked = np.ma.masked_invalid(matrix)
  image = axis.imshow(masked, aspect="auto", vmin=0, vmax=1)
  labels = [f"{row.pair_id}\n{row.cancer}" for row in data.itertuples()]
  axis.set_yticks(np.arange(len(labels)), labels, fontsize=7)
  axis.set_xticks(
    np.arange(len(observed_layers)),
    [
      value.replace("_support", "").replace("_", "\n")
      for value in observed_layers
    ],
    fontsize=8,
  )
  axis.set_title("Scored causal and orthogonal multi-omics layers")
  fig.colorbar(image, ax=axis, label="Normalized support")
  return fig, data


def figure_proteomics(
  protein: pd.DataFrame,
) -> tuple[plt.Figure, pd.DataFrame]:
  status = protein.get("status", pd.Series("", index=protein.index))
  data = protein.loc[status.astype(str).eq("ok")].copy()
  if not data.empty:
    data["label"] = (
      data["pair_id"].astype(str)
      + " | "
      + data["cancer"].astype(str)
    )
    top_labels = (
      data.groupby("label")["protein_compensation_support"]
        .max()
        .sort_values(ascending=False)
        .head(20)
        .index
    )
    data = data.loc[data["label"].isin(top_labels)]
  fig, axis = plt.subplots(figsize=(11.5, 7.8), constrained_layout=True)
  if data.empty:
    placeholder(
      axis,
      "Proteomic compensation concordance",
      "No evaluable protein-compensation contrasts were available.",
    )
    return fig, data
  pivot = data.pivot_table(
    index="label",
    columns="source",
    values="protein_compensation_support",
    aggfunc="median",
  )
  pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=False).index]
  image = axis.imshow(
    np.ma.masked_invalid(pivot.to_numpy()),
    aspect="auto",
    vmin=0,
    vmax=1,
  )
  axis.set_yticks(np.arange(len(pivot.index)), pivot.index, fontsize=7)
  axis.set_xticks(
    np.arange(len(pivot.columns)),
    [
      str(value).replace("proteomics_", "").replace("_", "\n")
      for value in pivot.columns
    ],
    fontsize=8,
  )
  axis.set_title("Platform-specific target-protein compensation support")
  fig.colorbar(image, ax=axis, label="Normalized compensation support")
  return fig, data


def figure_loss_states(
  loss: pd.DataFrame,
) -> tuple[plt.Figure, pd.DataFrame]:
  if loss.empty or "functional_loss_state" not in loss.columns:
    fig, axis = plt.subplots(figsize=(8.5, 5.8), constrained_layout=True)
    placeholder(
      axis,
      "Integrated functional-loss states",
      "No functional-loss state table was available.",
    )
    return fig, loss
  summary = (
    loss.groupby("functional_loss_state")
      .agg(
        model_gene_records=("ModelID", "size"),
        genes=("gene", "nunique"),
      )
      .reset_index()
      .sort_values("model_gene_records", ascending=False)
  )
  fig, axis = plt.subplots(figsize=(9.5, 6.0), constrained_layout=True)
  x = np.arange(len(summary))
  axis.bar(x, summary["model_gene_records"])
  axis.set_xticks(
    x,
    [
      str(value).replace("_", "\n")
      for value in summary["functional_loss_state"]
    ],
  )
  axis.set_ylabel("Model–gene records")
  axis.set_title("Integrated functional-loss state composition")
  axis.grid(axis="y", alpha=0.25)
  return fig, summary


def figure_crispr_rnai(
  evidence: pd.DataFrame,
) -> tuple[plt.Figure, pd.DataFrame]:
  columns = [
    column
    for column in (
      "pair_id",
      "cancer",
      "dependency_probability_support",
      "rnai_orthogonal_support",
      "dependency_probability_n_loss",
      "rnai_n_loss",
    )
    if column in evidence.columns
  ]
  data = (
    evidence[columns].dropna(
      subset=[
        "dependency_probability_support",
        "rnai_orthogonal_support",
      ],
      how="any",
    )
    if {
      "dependency_probability_support",
      "rnai_orthogonal_support",
    }.issubset(columns)
    else pd.DataFrame()
  )
  fig, axis = plt.subplots(figsize=(7.2, 6.4), constrained_layout=True)
  if data.empty:
    placeholder(
      axis,
      "CRISPR and RNAi orthogonal support",
      "No candidate-context rows had both CRISPR dependency-probability "
      "and RNAi support.",
    )
    return fig, data
  for cancer, subset in data.groupby("cancer"):
    axis.scatter(
      subset["dependency_probability_support"],
      subset["rnai_orthogonal_support"],
      s=24,
      alpha=0.65,
      label=str(cancer),
    )
  axis.set_xlim(0, 1)
  axis.set_ylim(0, 1)
  axis.set_xlabel("CRISPR dependency-probability support")
  axis.set_ylabel("RNAi orthogonal support")
  axis.set_title("Cross-technology vulnerability support")
  axis.legend(frameon=False)
  axis.grid(alpha=0.2)
  return fig, data


def figure_context(
  covariates: pd.DataFrame,
  evidence: pd.DataFrame,
) -> tuple[plt.Figure, pd.DataFrame]:
  frames = []
  if not covariates.empty:
    subset = covariates[[
      column
      for column in (
        "pair_id",
        "cancer",
        "source",
        "global_shift_exploratory",
      )
      if column in covariates.columns
    ]].copy()
    frames.append(subset)
  if (
    not evidence.empty
    and "metabolomic_state_shift_exploratory" in evidence.columns
  ):
    subset = evidence[[
      "pair_id",
      "cancer",
      "metabolomic_state_shift_exploratory",
    ]].copy()
    subset["source"] = "metabolomics"
    subset = subset.rename(columns={
      "metabolomic_state_shift_exploratory": "global_shift_exploratory"
    })
    frames.append(subset)
  data = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
  fig, axis = plt.subplots(figsize=(11.5, 6.6), constrained_layout=True)
  if data.empty:
    placeholder(
      axis,
      "Non-causal context shifts",
      "No exploratory context-layer shifts were available.",
    )
    return fig, data
  order = (
    data.groupby("source")["global_shift_exploratory"]
      .median()
      .sort_values(ascending=False)
      .index
  )
  box_data = [
    pd.to_numeric(
      data.loc[
        data["source"].eq(source),
        "global_shift_exploratory",
      ],
      errors="coerce",
    ).dropna().to_numpy()
    for source in order
  ]
  valid = [
    (source, values)
    for source, values in zip(order, box_data)
    if len(values)
  ]
  if not valid:
    placeholder(
      axis,
      "Non-causal context shifts",
      "Context sources were present but no evaluable contrasts were available.",
    )
    return fig, data
  axis.boxplot(
    [values for _, values in valid],
    labels=[source.replace("_", "\n") for source, _ in valid],
    showfliers=False,
  )
  axis.set_ylim(0, 1)
  axis.set_ylabel("Exploratory multivariate shift")
  axis.set_title("Context layers reserved for interpretation and robustness")
  axis.tick_params(axis="x", labelsize=7)
  axis.grid(axis="y", alpha=0.2)
  return fig, data


def figure_rank_stability(
  ranking: pd.DataFrame,
) -> tuple[plt.Figure, pd.DataFrame]:
  columns = [
    column
    for column in (
      "pair_id",
      "cancer",
      "baseline_rank_within_cancer",
      "extended_rank_within_cancer",
      "extended_rank_change",
      "baseline_coverage_adjusted_rses",
      "coverage_adjusted_rses",
    )
    if column in ranking.columns
  ]
  data = ranking[columns].copy() if columns else pd.DataFrame()
  fig, axis = plt.subplots(figsize=(9.5, 6.2), constrained_layout=True)
  if data.empty or "extended_rank_change" not in data.columns:
    placeholder(
      axis,
      "Extended rank stability",
      "No baseline-to-extended rank comparison was available.",
    )
    return fig, data
  groups = []
  labels = []
  for cancer, subset in data.groupby("cancer"):
    values = pd.to_numeric(
      subset["extended_rank_change"], errors="coerce"
    ).dropna()
    if len(values):
      groups.append(values.to_numpy())
      labels.append(str(cancer))
  if not groups:
    placeholder(
      axis,
      "Extended rank stability",
      "Rank changes were unavailable after filtering.",
    )
    return fig, data
  axis.boxplot(groups, labels=labels, showfliers=False)
  axis.axhline(0, linestyle="--", linewidth=1.0)
  axis.set_ylabel("Baseline rank − extended rank\n(positive = improved rank)")
  axis.set_title("Cancer-specific ranking stability after multi-omics extension")
  axis.grid(axis="y", alpha=0.2)
  return fig, data


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--config", default="config/extended_multiomics_asset.yaml"
  )
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--processed-dir", default="data/processed/extended_multiomics"
  )
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument(
    "--strict-layout",
    action=argparse.BooleanOptionalAction,
    default=True,
  )
  args = parser.parse_args()

  config = yaml.safe_load(
    resolve(args.config).read_text(encoding="utf-8")
  ) or {}
  registry = {
    str(item["id"]): item
    for item in config.get("supplementary_figures", [])
  }
  expected = {f"Figure_S{index}" for index in range(71, 79)}
  if set(registry) != expected:
    raise RuntimeError(
      "Extended figure config must register S71-S78; "
      f"observed={sorted(registry)}"
    )

  processed = resolve(args.processed_dir)
  ranking_path = resolve(args.ranking)
  output_root = resolve(args.output_root)
  paths = {
    "inventory": processed / "extended_multiomics_source_inventory.tsv",
    "evidence": processed / "extended_pair_evidence_by_cancer.tsv",
    "protein": processed / "proteomics_pair_evidence_by_source.tsv",
    "loss": processed / "functional_loss_states.tsv",
    "covariates": processed / "extended_covariate_context.tsv",
  }
  ranking = pd.read_csv(ranking_path, sep="\t", low_memory=False)
  inventory = read_optional(paths["inventory"])
  evidence = read_optional(paths["evidence"])
  protein = read_optional(paths["protein"])
  loss = read_optional(paths["loss"])
  covariates = read_optional(paths["covariates"])

  set_publication_style()
  builders = {
    "Figure_S71": lambda: figure_source_coverage(inventory),
    "Figure_S72": lambda: figure_baseline_extended(ranking),
    "Figure_S73": lambda: figure_layer_heatmap(evidence, ranking),
    "Figure_S74": lambda: figure_proteomics(protein),
    "Figure_S75": lambda: figure_loss_states(loss),
    "Figure_S76": lambda: figure_crispr_rnai(evidence),
    "Figure_S77": lambda: figure_context(covariates, evidence),
    "Figure_S78": lambda: figure_rank_stability(ranking),
  }
  records = []
  for figure_id in sorted(
    expected,
    key=lambda value: int(value.split("S")[1]),
  ):
    fig, source = builders[figure_id]()
    records.append(save_asset(
      fig,
      registry[figure_id],
      source,
      [ranking_path, *paths.values()],
      output_root,
      args.strict_layout,
    ))
    print(f"Generated {figure_id}: {records[-1].layout_status}")
  write_figure_manifest(
    records,
    output_root / "manifests/extended_multiomics_figure_manifest.tsv",
  )
  write_legends_markdown(
    records,
    output_root
    / "manuscript_assets/extended_multiomics_figure_legends.md",
  )


if __name__ == "__main__":
  main()
