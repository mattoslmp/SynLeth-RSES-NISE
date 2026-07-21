#!/usr/bin/env python3
"""Generate every main article figure exclusively from analysis tables."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path
import string

import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError
import yaml

from rses_onco.publication import (
  FigureRecord,
  dynamic_height,
  figure_record,
  placeholder,
  save_figure_triplet,
  set_publication_style,
  wrap_label,
  write_figure_manifest,
  write_legends_markdown,
  write_source_data,
)

ROOT = Path(__file__).resolve().parents[1]
CANCER_LABELS = {"colon": "Colorectal", "stomach": "Gastric", "lung": "Lung"}


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


def registry_by_id(config: dict, section: str) -> dict[str, dict]:
  return {str(item["id"]): item for item in config.get(section, [])}


def panel_label(axis: plt.Axes, label: str) -> None:
  axis.text(
    0.01,
    0.99,
    label,
    transform=axis.transAxes,
    ha="left",
    va="top",
    fontsize=13,
    fontweight="bold",
    bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
    zorder=20,
  )


def pair_labels(frame: pd.DataFrame, width: int = 28) -> list[str]:
  lost = frame.get("analysis_lost_gene", frame.get("lost_gene", frame.get("lost_feature", "")))
  target = frame.get("analysis_target_gene", frame.get("target_gene", ""))
  return [wrap_label(f"{a} → {b}", width) for a, b in zip(lost.astype(str), target.astype(str))]


def save_record(
  *,
  fig: plt.Figure,
  item: dict,
  output_dir: Path,
  source_dir: Path,
  source_data: pd.DataFrame,
  inputs: list[Path],
  strict_layout: bool,
) -> FigureRecord:
  base = output_dir / str(item["file"])
  source_path = write_source_data(source_data, source_dir / f"{item['file']}_source_data.tsv")
  audit = save_figure_triplet(fig, base, str(item["id"]), strict_layout=strict_layout)
  return figure_record(
    figure_id=str(item["id"]),
    category="main",
    title=str(item["title"]),
    caption=str(item.get("caption") or item["title"]),
    base_path=base,
    source_data_path=source_path,
    input_paths=inputs,
    audit=audit,
    script="scripts/make_main_figures.py",
  )


def figure_1(item: dict, output_dir: Path, source_dir: Path, strict: bool) -> FigureRecord:
  set_publication_style()
  fig, axis = plt.subplots(figsize=(15.5, 8.5), constrained_layout=True)
  axis.set_axis_off()
  stages = [
    ("Candidate universe", "All NISE directions\nEnsembl paralogs\nCurated classes\nAll-target discoveries"),
    ("Human microniches", "DepMap expression\nHPA localization\nUniProt/PDB\nCRISPR phenotypes\nSTRING + DoRothEA"),
    ("Cancer evidence", "TCGA/GDC events\nDepMap dependency\nLoss selectivity\nExpression compensation"),
    ("Pharmacology", "Open Targets\nChEMBL\nDGIdb + MyChem\nPharos + CIViC\nPRISM/GDSC/CTRP"),
    ("Validation path", "Biomarker-matched models\nOrthogonal inhibition\nRescue experiments\nMechanistic assays\nIn vivo validation"),
  ]
  x_positions = np.linspace(0.03, 0.81, len(stages))
  width = 0.16
  height = 0.48
  rows = []
  for index, ((title, body), x) in enumerate(zip(stages, x_positions), start=1):
    box = patches.FancyBboxPatch(
      (x, 0.28), width, height,
      boxstyle="round,pad=0.018",
      linewidth=1.4,
      facecolor="white",
      edgecolor="0.25",
      transform=axis.transAxes,
    )
    axis.add_patch(box)
    axis.text(x + width / 2, 0.69, title, transform=axis.transAxes, ha="center", va="center", fontsize=11.5, fontweight="bold")
    axis.text(x + width / 2, 0.49, body, transform=axis.transAxes, ha="center", va="center", fontsize=9.5, linespacing=1.45)
    rows.append({"stage": index, "title": title, "contents": body.replace("\n", ";")})
    if index < len(stages):
      axis.annotate(
        "",
        xy=(x + width + 0.025, 0.52),
        xytext=(x + width + 0.005, 0.52),
        xycoords=axis.transAxes,
        arrowprops={"arrowstyle": "-|>", "linewidth": 1.5},
      )
  axis.text(
    0.5,
    0.13,
    "Coverage-aware scores preserve missing evidence as missing and prioritize experimental hypotheses, not clinical efficacy claims.",
    transform=axis.transAxes,
    ha="center",
    va="center",
    fontsize=10.5,
    fontweight="bold",
  )
  axis.set_title(str(item["title"]), pad=18)
  return save_record(
    fig=fig,
    item=item,
    output_dir=output_dir,
    source_dir=source_dir,
    source_data=pd.DataFrame(rows),
    inputs=[ROOT / "config/article_assets.yaml"],
    strict_layout=strict,
  )


def figure_2(item: dict, candidates: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path) -> FigureRecord:
  set_publication_style()
  fig = plt.figure(figsize=(15.5, 9.0), constrained_layout=True)
  grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.35])
  axis_class = fig.add_subplot(grid[0, 0])
  axis_activity = fig.add_subplot(grid[0, 1])

  class_counts = (
    candidates.assign(source_class=candidates.get("source_class", "unclassified").fillna("unclassified").astype(str))
      .groupby("source_class", as_index=False)
      .agg(directed_candidates=("pair_id", "nunique"))
      .sort_values("directed_candidates")
  )
  axis_class.barh([wrap_label(value, 25) for value in class_counts["source_class"]], class_counts["directed_candidates"])
  axis_class.set_xlabel("Directed candidate hypotheses")
  axis_class.set_title("Mechanistic classes")
  axis_class.grid(axis="x", alpha=0.25)
  panel_label(axis_class, "A")

  nise = candidates.loc[candidates.get("source_class", pd.Series(index=candidates.index, dtype=object)).astype(str).eq("NISE")].copy()
  if nise.empty or "activity" not in nise:
    placeholder(axis_activity, "Human NISE activities", "No NISE activity records were available.")
    activity_counts = pd.DataFrame(columns=["activity", "directed_candidates", "genes"])
  else:
    activity_counts = (
      nise.groupby(["group_id", "activity"], as_index=False)
        .agg(
          directed_candidates=("pair_id", "nunique"),
          lost_genes=("lost_gene", "nunique"),
          target_genes=("target_gene", "nunique"),
        )
    )
    activity_counts["genes"] = activity_counts[["lost_genes", "target_genes"]].max(axis=1)
    activity_counts = activity_counts.sort_values("directed_candidates")
    labels = [wrap_label(f"{group}: {activity}", 34) for group, activity in zip(activity_counts["group_id"], activity_counts["activity"])]
    axis_activity.barh(labels, activity_counts["directed_candidates"])
    axis_activity.set_xlabel("Directed NISE hypotheses")
    axis_activity.set_title("Complete curated NISE activity groups")
    axis_activity.grid(axis="x", alpha=0.25)
  panel_label(axis_activity, "B")
  fig.suptitle(str(item["title"]))
  source = pd.concat([
    class_counts.assign(panel="A_class_counts"),
    activity_counts.assign(panel="B_nise_activities"),
  ], ignore_index=True, sort=False)
  return save_record(
    fig=fig, item=item, output_dir=output_dir, source_dir=source_dir,
    source_data=source, inputs=[input_path], strict_layout=strict,
  )


def figure_3(item: dict, ranking: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path, top_n: int) -> FigureRecord:
  set_publication_style()
  cancers = [value for value in ("colon", "stomach", "lung") if value in set(ranking["cancer"].astype(str))]
  fig, axes = plt.subplots(1, len(cancers), figsize=(6.4 * len(cancers), 9.5), constrained_layout=True, squeeze=False)
  source_frames = []
  for index, (axis, cancer) in enumerate(zip(axes[0], cancers)):
    subset = (
      ranking.loc[ranking["cancer"].astype(str).eq(cancer)]
        .sort_values("coverage_adjusted_rses", ascending=False)
        .head(top_n)
        .sort_values("coverage_adjusted_rses")
        .copy()
    )
    axis.barh(pair_labels(subset, 24), subset["coverage_adjusted_rses"])
    axis.set_xlim(0, max(1.0, float(subset["coverage_adjusted_rses"].max()) * 1.10 if not subset.empty else 1.0))
    axis.set_xlabel("Coverage-adjusted RSES-Onco")
    axis.set_title(CANCER_LABELS[cancer])
    axis.grid(axis="x", alpha=0.25)
    panel_label(axis, string.ascii_uppercase[index])
    source_frames.append(subset.assign(panel=string.ascii_uppercase[index]))
  fig.suptitle(str(item["title"]))
  source = pd.concat(source_frames, ignore_index=True, sort=False) if source_frames else pd.DataFrame()
  return save_record(
    fig=fig, item=item, output_dir=output_dir, source_dir=source_dir,
    source_data=source, inputs=[input_path], strict_layout=strict,
  )


def figure_4(item: dict, ranking: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path) -> FigureRecord:
  set_publication_style()
  required = {"tcga_homdel_frequency", "component_selectivity", "coverage_adjusted_rses"}
  subset = ranking.dropna(subset=list(required)).copy() if required.issubset(ranking.columns) else pd.DataFrame()
  fig = plt.figure(figsize=(17.5, 9.0), constrained_layout=True)
  grid = fig.add_gridspec(2, 3, height_ratios=[3.0, 1.15])
  source_frames = []
  key_rows = []
  for index, cancer in enumerate(("colon", "stomach", "lung")):
    axis = fig.add_subplot(grid[0, index])
    group = subset.loc[subset["cancer"].astype(str).eq(cancer)].copy()
    if group.empty:
      placeholder(axis, CANCER_LABELS[cancer], "No paired TCGA event and DepMap selectivity observations were available.")
      continue
    axis.scatter(
      group["tcga_homdel_frequency"],
      group["component_selectivity"],
      s=45 + 260 * group["coverage_adjusted_rses"].fillna(0),
      alpha=0.72,
      edgecolors="0.25",
      linewidths=0.5,
    )
    axis.set_xlabel("TCGA homozygous-deletion frequency")
    axis.set_ylabel("DepMap loss-selectivity component")
    axis.set_title(CANCER_LABELS[cancer])
    axis.grid(alpha=0.25)
    panel_label(axis, string.ascii_uppercase[index])
    top = group.sort_values("coverage_adjusted_rses", ascending=False).head(5).copy()
    for number, (_, row) in enumerate(top.iterrows(), start=1):
      axis.scatter([row["tcga_homdel_frequency"]], [row["component_selectivity"]], marker=f"${number}$", s=160, color="black")
      key_rows.append({
        "panel": string.ascii_uppercase[index],
        "number": number,
        "cancer": cancer,
        "pair": f"{row.get('analysis_lost_gene', row.get('lost_gene'))} → {row.get('target_gene')}",
        "tcga_homdel_frequency": row["tcga_homdel_frequency"],
        "component_selectivity": row["component_selectivity"],
        "coverage_adjusted_rses": row["coverage_adjusted_rses"],
      })
    source_frames.append(group.assign(panel=string.ascii_uppercase[index]))

  key_axis = fig.add_subplot(grid[1, :])
  key_axis.set_axis_off()
  key_axis.set_title("Numbered high-priority points", loc="left", pad=4)
  if key_rows:
    text = "    ".join(
      f"{row['panel']}{row['number']}: {row['pair']}"
      for row in key_rows
    )
    key_axis.text(0.01, 0.55, wrap_label(text, 150), transform=key_axis.transAxes, va="center", ha="left", fontsize=9.5)
  else:
    key_axis.text(0.01, 0.55, "No eligible points.", transform=key_axis.transAxes, va="center")
  panel_label(key_axis, "D")
  fig.suptitle(str(item["title"]))
  source = pd.concat(source_frames, ignore_index=True, sort=False) if source_frames else pd.DataFrame()
  if key_rows:
    source = pd.concat([source, pd.DataFrame(key_rows).assign(record_type="key")], ignore_index=True, sort=False)
  return save_record(
    fig=fig, item=item, output_dir=output_dir, source_dir=source_dir,
    source_data=source, inputs=[input_path], strict_layout=strict,
  )


def figure_5(item: dict, ranking: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path, top_n: int) -> FigureRecord:
  set_publication_style()
  domains = [
    "microniche_expression_context",
    "microniche_localization",
    "microniche_biochemical_structural",
    "microniche_genetic_phenotype",
    "microniche_interaction_network",
    "microniche_regulatory_network",
  ]
  available = [column for column in domains if column in ranking.columns]
  pair_level = (
    ranking.sort_values("functional_microniche_adjusted", ascending=False)
      .drop_duplicates("pair_id")
      .head(top_n)
      .copy()
    if "functional_microniche_adjusted" in ranking
    else pd.DataFrame()
  )
  fig = plt.figure(figsize=(13.5, dynamic_height(len(pair_level), minimum=8.0, per_row=0.32)), constrained_layout=True)
  grid = fig.add_gridspec(1, 2, width_ratios=[24, 1])
  axis = fig.add_subplot(grid[0, 0])
  color_axis = fig.add_subplot(grid[0, 1], label="colorbar")
  if pair_level.empty or not available:
    color_axis.set_axis_off()
    placeholder(axis, str(item["title"]), "No functional-microniche component matrix was available.")
    source = pd.DataFrame()
  else:
    matrix = pair_level[available].apply(pd.to_numeric, errors="coerce")
    image = axis.imshow(matrix.to_numpy(), aspect="auto", vmin=0, vmax=1)
    axis.set_xticks(np.arange(len(available)), [wrap_label(column.removeprefix("microniche_").replace("_", " ").title(), 18) for column in available], rotation=25, ha="right")
    axis.set_yticks(np.arange(len(pair_level)), pair_labels(pair_level, 30))
    axis.set_title(str(item["title"]), pad=14)
    colorbar = fig.colorbar(image, cax=color_axis)
    colorbar.set_label("Specialization/divergence evidence")
    source = pair_level[[column for column in ["pair_id", "lost_gene", "target_gene", "source_class", *available, "functional_microniche_adjusted"] if column in pair_level]].copy()
  return save_record(
    fig=fig, item=item, output_dir=output_dir, source_dir=source_dir,
    source_data=source, inputs=[input_path], strict_layout=strict,
  )


def figure_6(item: dict, ranking: pd.DataFrame, discovery: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_paths: list[Path]) -> FigureRecord:
  set_publication_style()
  fig = plt.figure(figsize=(16.0, 8.5), constrained_layout=True)
  grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.4])
  axis_class = fig.add_subplot(grid[0, 0])
  axis_discovery = fig.add_subplot(grid[0, 1])

  class_summary = (
    ranking.assign(source_class=ranking.get("source_class", "unclassified").fillna("unclassified").astype(str))
      .groupby("source_class", as_index=False)
      .agg(
        unique_directions=("pair_id", "nunique"),
        maximum_score=("coverage_adjusted_rses", "max"),
        median_score=("coverage_adjusted_rses", "median"),
      )
      .sort_values("maximum_score")
  )
  axis_class.barh([wrap_label(value, 24) for value in class_summary["source_class"]], class_summary["maximum_score"])
  axis_class.set_xlabel("Maximum coverage-adjusted RSES-Onco")
  axis_class.set_title("Best evidence within each class")
  axis_class.grid(axis="x", alpha=0.25)
  panel_label(axis_class, "A")

  if discovery.empty or "delta_effect" not in discovery:
    placeholder(axis_discovery, "All-target discoveries", "No FDR-supported all-target discovery rows were available.")
    top_discovery = pd.DataFrame()
  else:
    q_column = "q_value_bh_within_loss_cancer" if "q_value_bh_within_loss_cancer" in discovery else "q_value_bh"
    top_discovery = discovery.sort_values([q_column, "delta_effect"], ascending=[True, True]).head(25).copy().sort_values("delta_effect", ascending=False)
    labels = [wrap_label(f"{lost} → {target} | {CANCER_LABELS.get(str(cancer), cancer)}", 34) for lost, target, cancer in zip(top_discovery["lost_gene"], top_discovery["target_gene"], top_discovery["cancer"])]
    axis_discovery.barh(labels, top_discovery["delta_effect"])
    axis_discovery.axvline(0, linewidth=1, color="black")
    axis_discovery.set_xlabel("Δ CRISPR effect (loss − intact); more negative is supportive")
    axis_discovery.set_title("All-target conditional-dependency discoveries")
    axis_discovery.grid(axis="x", alpha=0.25)
  panel_label(axis_discovery, "B")
  fig.suptitle(str(item["title"]))
  source = pd.concat([
    class_summary.assign(panel="A_class_summary"),
    top_discovery.assign(panel="B_all_target_discovery"),
  ], ignore_index=True, sort=False)
  return save_record(
    fig=fig, item=item, output_dir=output_dir, source_dir=source_dir,
    source_data=source, inputs=input_paths, strict_layout=strict,
  )


def figure_7(item: dict, pharmacology: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path, top_n: int) -> FigureRecord:
  set_publication_style()
  fig = plt.figure(figsize=(17.0, 9.0), constrained_layout=True)
  grid = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.35])
  axis_scatter = fig.add_subplot(grid[0, 0])
  axis_table = fig.add_subplot(grid[0, 1])
  axis_table.set_axis_off()

  eligible = pharmacology.dropna(subset=["coverage_adjusted_rses", "pharmacology_adjusted", "therapeutic_hypothesis_score"]).copy() if not pharmacology.empty else pd.DataFrame()
  if eligible.empty:
    placeholder(axis_scatter, "Pharmacological actionability", "No pharmacology hypotheses were available.")
    placeholder(axis_table, "Top experimental hypotheses", "Run the pharmacology acquisition and prioritization stages.")
    top = pd.DataFrame()
  else:
    axis_scatter.scatter(
      eligible["coverage_adjusted_rses"],
      eligible["pharmacology_adjusted"],
      s=45 + 280 * eligible["therapeutic_hypothesis_score"].fillna(0),
      alpha=0.70,
      edgecolors="0.25",
      linewidths=0.5,
    )
    axis_scatter.set_xlabel("Coverage-adjusted vulnerability")
    axis_scatter.set_ylabel("Coverage-adjusted pharmacology")
    axis_scatter.set_title("Concordance of vulnerability and actionability")
    axis_scatter.grid(alpha=0.25)
    top = eligible.sort_values("therapeutic_hypothesis_score", ascending=False).head(top_n).copy()
    for number, (_, row) in enumerate(top.iterrows(), start=1):
      axis_scatter.scatter(
        [row["coverage_adjusted_rses"]],
        [row["pharmacology_adjusted"]],
        marker=f"${number}$",
        s=170,
        color="black",
      )
    table_rows = []
    for number, (_, row) in enumerate(top.iterrows(), start=1):
      pair = f"{row.get('analysis_lost_gene', row.get('lost_gene'))} → {row.get('target_gene')}"
      drug = row.get("drug_name") or row.get("drug_id") or "target-only evidence"
      table_rows.append([
        number,
        wrap_label(pair, 23),
        wrap_label(drug, 20),
        f"{row['therapeutic_hypothesis_score']:.3f}",
      ])
    table = axis_table.table(
      cellText=table_rows,
      colLabels=["#", "Vulnerability", "Compound / evidence", "Score"],
      cellLoc="left",
      colLoc="left",
      loc="center",
      colWidths=[0.06, 0.34, 0.42, 0.12],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.65)
    axis_table.set_title("Top experimental target-drug hypotheses", pad=12)
  panel_label(axis_scatter, "A")
  panel_label(axis_table, "B")
  fig.suptitle(str(item["title"]))
  return save_record(
    fig=fig, item=item, output_dir=output_dir, source_dir=source_dir,
    source_data=top, inputs=[input_path], strict_layout=strict,
  )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--config", default="config/article_assets.yaml")
  parser.add_argument("--ranking", default="results/expanded_26Q1/full/expanded_rses_onco.tsv")
  parser.add_argument("--candidates", default="data/processed/expanded_candidate_universe.tsv")
  parser.add_argument("--discovery", default="results/expanded_26Q1/discovery/all_target_dependency_screen.tsv")
  parser.add_argument("--pharmacology", default="results/expanded_26Q1/pharmacology/pharmacology_ranked_hypotheses.tsv")
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument("--top-n", type=int, default=15)
  parser.add_argument("--strict-layout", action=argparse.BooleanOptionalAction, default=True)
  args = parser.parse_args()

  config_path = resolve_path(args.config)
  config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
  registry = registry_by_id(config, "main_figures")
  required_ids = {f"Figure_{index}" for index in range(1, 8)}
  missing = required_ids - set(registry)
  if missing:
    raise ValueError(f"Main figure registry is incomplete: {sorted(missing)}")

  ranking_path = resolve_path(args.ranking)
  candidates_path = resolve_path(args.candidates)
  discovery_path = resolve_path(args.discovery)
  pharmacology_path = resolve_path(args.pharmacology)
  ranking = pd.read_csv(ranking_path, sep="\t")
  candidates = pd.read_csv(candidates_path, sep="\t")
  discovery = read_optional(discovery_path)
  pharmacology = read_optional(pharmacology_path)

  output_root = resolve_path(args.output_root)
  output_dir = output_root / "figures" / "main"
  source_dir = output_root / "source_data" / "figures" / "main"
  manifest_dir = output_root / "manifests"
  legend_dir = output_root / "manuscript_assets"
  records: list[FigureRecord] = []
  records.append(figure_1(registry["Figure_1"], output_dir, source_dir, args.strict_layout))
  records.append(figure_2(registry["Figure_2"], candidates, output_dir, source_dir, args.strict_layout, candidates_path))
  records.append(figure_3(registry["Figure_3"], ranking, output_dir, source_dir, args.strict_layout, ranking_path, args.top_n))
  records.append(figure_4(registry["Figure_4"], ranking, output_dir, source_dir, args.strict_layout, ranking_path))
  records.append(figure_5(registry["Figure_5"], ranking, output_dir, source_dir, args.strict_layout, ranking_path, max(args.top_n * 2, 25)))
  records.append(figure_6(registry["Figure_6"], ranking, discovery, output_dir, source_dir, args.strict_layout, [ranking_path, discovery_path]))
  records.append(figure_7(registry["Figure_7"], pharmacology, output_dir, source_dir, args.strict_layout, pharmacology_path, args.top_n))

  write_figure_manifest(records, manifest_dir / "main_figure_manifest.tsv")
  write_legends_markdown(records, legend_dir / "main_figure_legends.md")
  print(pd.DataFrame([asdict(record) for record in records])[["figure_id", "layout_status", "base_path"]].to_string(index=False))
  print(f"Wrote all main figures to {output_dir}")


if __name__ == "__main__":
  main()
