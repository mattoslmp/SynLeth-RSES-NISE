#!/usr/bin/env python3
"""Generate every supplementary article figure exclusively from source tables."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import string

import matplotlib.pyplot as plt
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


def registry_by_id(config: dict) -> dict[str, dict]:
  return {str(item["id"]): item for item in config.get("supplementary_figures", [])}


def pair_label(frame: pd.DataFrame, width: int = 32) -> list[str]:
  lost = frame.get("analysis_lost_gene", frame.get("lost_gene", frame.get("lost_feature", "")))
  target = frame.get("analysis_target_gene", frame.get("target_gene", ""))
  return [wrap_label(f"{a} → {b}", width) for a, b in zip(lost.astype(str), target.astype(str))]


def save_record(
  fig: plt.Figure,
  item: dict,
  output_dir: Path,
  source_dir: Path,
  source: pd.DataFrame,
  inputs: list[Path],
  strict: bool,
) -> FigureRecord:
  base = output_dir / str(item["file"])
  source_path = write_source_data(source, source_dir / f"{item['file']}_source_data.tsv")
  audit = save_figure_triplet(fig, base, str(item["id"]), strict_layout=strict)
  return figure_record(
    figure_id=str(item["id"]),
    category="supplementary",
    title=str(item["title"]),
    caption=str(item.get("caption") or item["title"]),
    base_path=base,
    source_data_path=source_path,
    input_paths=inputs,
    audit=audit,
    script="scripts/make_supplementary_figures.py",
  )


def simple_bar(item: dict, frame: pd.DataFrame, label: str, value: str, output_dir: Path, source_dir: Path, inputs: list[Path], strict: bool, xlabel: str) -> FigureRecord:
  set_publication_style()
  ordered = frame.sort_values(value).copy() if not frame.empty else frame
  fig, axis = plt.subplots(figsize=(11.5, dynamic_height(len(ordered), minimum=6.0, per_row=0.30)), constrained_layout=True)
  if ordered.empty:
    placeholder(axis, str(item["title"]), "No eligible observations were available in this release.")
  else:
    axis.barh([wrap_label(item_value, 36) for item_value in ordered[label].astype(str)], ordered[value])
    axis.set_xlabel(xlabel)
    axis.set_title(str(item["title"]), pad=12)
    axis.grid(axis="x", alpha=0.25)
  return save_record(fig, item, output_dir, source_dir, ordered, inputs, strict)


def heatmap_record(
  item: dict,
  matrix: pd.DataFrame,
  row_labels: list[str],
  output_dir: Path,
  source_dir: Path,
  source: pd.DataFrame,
  inputs: list[Path],
  strict: bool,
  colorbar_label: str,
  vmin: float | None = None,
  vmax: float | None = None,
) -> FigureRecord:
  set_publication_style()
  fig = plt.figure(figsize=(12.5, dynamic_height(len(matrix), minimum=7.0, per_row=0.28)), constrained_layout=True)
  grid = fig.add_gridspec(1, 2, width_ratios=[24, 1])
  axis = fig.add_subplot(grid[0, 0])
  color_axis = fig.add_subplot(grid[0, 1], label="colorbar")
  if matrix.empty:
    color_axis.set_axis_off()
    placeholder(axis, str(item["title"]), "No eligible observations were available in this release.")
  else:
    image = axis.imshow(matrix.to_numpy(dtype=float), aspect="auto", vmin=vmin, vmax=vmax)
    axis.set_xticks(np.arange(len(matrix.columns)), [wrap_label(value, 18) for value in matrix.columns], rotation=25, ha="right")
    axis.set_yticks(np.arange(len(matrix)), row_labels)
    axis.set_title(str(item["title"]), pad=12)
    colorbar = fig.colorbar(image, cax=color_axis)
    colorbar.set_label(colorbar_label)
  return save_record(fig, item, output_dir, source_dir, source, inputs, strict)


def figure_s1(item: dict, ranking: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path) -> FigureRecord:
  domains = [
    column for column in ranking.columns
    if column.startswith("component_") or column.startswith("microniche_")
  ]
  coverage = pd.DataFrame({
    "domain": domains,
    "available_fraction": [float(ranking[column].notna().mean()) for column in domains],
    "available_rows": [int(ranking[column].notna().sum()) for column in domains],
    "total_rows": len(ranking),
  })
  return simple_bar(item, coverage, "domain", "available_fraction", output_dir, source_dir, [input_path], strict, "Fraction of scored rows with evidence")


def figure_s2(item: dict, candidates: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path) -> FigureRecord:
  nise = candidates.loc[candidates.get("source_class", pd.Series(index=candidates.index, dtype=object)).astype(str).eq("NISE")].copy()
  if nise.empty or "activity" not in nise:
    summary = pd.DataFrame(columns=["activity_label", "directed_candidates"])
  else:
    summary = (
      nise.groupby(["group_id", "activity"], as_index=False)
        .agg(
          directed_candidates=("pair_id", "nunique"),
          lost_members=("lost_gene", "nunique"),
          target_members=("target_gene", "nunique"),
        )
    )
    summary["activity_label"] = summary["group_id"].astype(str) + ": " + summary["activity"].astype(str)
  return simple_bar(item, summary, "activity_label", "directed_candidates", output_dir, source_dir, [input_path], strict, "Directed NISE hypotheses")


def figure_s3(item: dict, discovery: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path) -> FigureRecord:
  set_publication_style()
  fig, axes = plt.subplots(1, 3, figsize=(17.5, 6.5), constrained_layout=True)
  source_frames = []
  q_column = "q_value_bh_within_loss_cancer" if "q_value_bh_within_loss_cancer" in discovery else "q_value_bh"
  for index, (axis, cancer) in enumerate(zip(axes, ("colon", "stomach", "lung"))):
    group = discovery.loc[discovery.get("cancer", pd.Series(index=discovery.index, dtype=object)).astype(str).eq(cancer)].copy()
    if group.empty or q_column not in group:
      placeholder(axis, CANCER_LABELS[cancer], "No all-target discovery rows were available.")
      continue
    q = pd.to_numeric(group[q_column], errors="coerce").clip(lower=1e-300)
    group["minus_log10_q"] = -np.log10(q)
    axis.scatter(group["delta_effect"], group["minus_log10_q"], alpha=0.65, s=25)
    axis.axvline(0, color="black", linewidth=0.8)
    axis.axhline(-np.log10(0.05), color="0.4", linewidth=0.8, linestyle="--")
    axis.set_xlabel("Δ CRISPR effect")
    axis.set_ylabel("−log10 adjusted P")
    axis.set_title(CANCER_LABELS[cancer])
    axis.grid(alpha=0.2)
    axis.text(0.02, 0.98, string.ascii_uppercase[index], transform=axis.transAxes, va="top", fontweight="bold", fontsize=13)
    source_frames.append(group.assign(panel=string.ascii_uppercase[index]))
  fig.suptitle(str(item["title"]))
  source = pd.concat(source_frames, ignore_index=True, sort=False) if source_frames else pd.DataFrame()
  return save_record(fig, item, output_dir, source_dir, source, [input_path], strict)


def contrast_heatmap(item: dict, frame: pd.DataFrame, value_column: str, output_dir: Path, source_dir: Path, strict: bool, input_path: Path, colorbar_label: str, top_n: int = 60) -> FigureRecord:
  if frame.empty or value_column not in frame:
    matrix = pd.DataFrame()
    source = pd.DataFrame()
    labels = []
  else:
    source = frame.copy()
    source["pair_label"] = pair_label(source, 38)
    source["magnitude"] = pd.to_numeric(source[value_column], errors="coerce").abs()
    selected = source.sort_values("magnitude", ascending=False).drop_duplicates(["pair_label", "cancer"]).head(top_n)
    matrix = selected.pivot_table(index="pair_label", columns="cancer", values=value_column, aggfunc="first")
    matrix = matrix[[column for column in ("colon", "stomach", "lung") if column in matrix.columns]]
    matrix.columns = [CANCER_LABELS.get(column, column) for column in matrix.columns]
    labels = [wrap_label(value, 36) for value in matrix.index]
    source = selected
  limit = float(np.nanmax(np.abs(matrix.to_numpy()))) if not matrix.empty else 1.0
  return heatmap_record(item, matrix, labels, output_dir, source_dir, source, [input_path], strict, colorbar_label, vmin=-limit, vmax=limit)


def profile_heatmap(item: dict, frame: pd.DataFrame, value_column: str, output_dir: Path, source_dir: Path, strict: bool, input_path: Path, colorbar_label: str, top_n: int = 70) -> FigureRecord:
  if frame.empty or value_column not in frame:
    matrix = pd.DataFrame()
    labels = []
    source = pd.DataFrame()
  else:
    source = frame.copy()
    source["pair_label"] = [wrap_label(f"{a} → {b}", 36) for a, b in zip(source["gene_a"].astype(str), source["gene_b"].astype(str))]
    source = source.sort_values(value_column, ascending=False).head(top_n)
    matrix = source.pivot_table(index="pair_label", columns="cancer", values=value_column, aggfunc="first")
    matrix = matrix[[column for column in ("colon", "stomach", "lung") if column in matrix.columns]]
    matrix.columns = [CANCER_LABELS.get(column, column) for column in matrix.columns]
    labels = list(matrix.index)
  return heatmap_record(item, matrix, labels, output_dir, source_dir, source, [input_path], strict, colorbar_label, vmin=0, vmax=1)


def figure_s8(item: dict, evidence: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path) -> FigureRecord:
  set_publication_style()
  fig, axis = plt.subplots(figsize=(10.5, 7.5), constrained_layout=True)
  subset = evidence.dropna(subset=[column for column in ["string_neighbor_jaccard"] if column in evidence]).copy() if not evidence.empty and "string_neighbor_jaccard" in evidence else pd.DataFrame()
  if subset.empty:
    placeholder(axis, str(item["title"]), "No STRING pair-level evidence was available.")
  else:
    direct = pd.to_numeric(subset.get("string_direct_score", 0), errors="coerce").fillna(0)
    axis.scatter(subset["string_neighbor_jaccard"], direct, s=30 + 10 * pd.to_numeric(subset.get("string_shared_neighbors", 0), errors="coerce").fillna(0).clip(0, 20), alpha=0.65)
    axis.set_xlabel("STRING neighborhood Jaccard")
    axis.set_ylabel("Direct STRING score")
    axis.set_title(str(item["title"]))
    axis.grid(alpha=0.25)
  return save_record(fig, item, output_dir, source_dir, subset, [input_path], strict)


def figure_s9(item: dict, evidence: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path) -> FigureRecord:
  set_publication_style()
  subset = evidence.dropna(subset=["regulator_jaccard"]).copy() if not evidence.empty and "regulator_jaccard" in evidence else pd.DataFrame()
  fig, axis = plt.subplots(figsize=(10.5, 7.0), constrained_layout=True)
  if subset.empty:
    placeholder(axis, str(item["title"]), "No DoRothEA regulator-set evidence was available.")
  else:
    axis.hist(subset["regulator_jaccard"], bins=np.linspace(0, 1, 21), edgecolor="white")
    axis.set_xlabel("Regulator-set Jaccard")
    axis.set_ylabel("Candidate pairs")
    axis.set_title(str(item["title"]))
    axis.grid(axis="y", alpha=0.25)
  return save_record(fig, item, output_dir, source_dir, subset, [input_path], strict)


def figure_s10(item: dict, evidence: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path, top_n: int = 50) -> FigureRecord:
  subset = evidence.dropna(subset=["component_localization"]).copy() if not evidence.empty and "component_localization" in evidence else pd.DataFrame()
  if not subset.empty:
    subset["label"] = [wrap_label(f"{a} → {b}", 35) for a, b in zip(subset["lost_gene"].astype(str), subset["target_gene"].astype(str))]
    subset = subset.sort_values("component_localization", ascending=False).drop_duplicates("label").head(top_n)
  return simple_bar(item, subset, "label", "component_localization", output_dir, source_dir, [input_path], strict, "Localization divergence (1 − Jaccard)")


def figure_s11(item: dict, event_summary: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path, top_n: int = 80) -> FigureRecord:
  if event_summary.empty:
    matrix = pd.DataFrame()
    labels = []
    source = pd.DataFrame()
  else:
    frequency_column = next((column for column in ["homdel_frequency", "tcga_homdel_frequency", "frequency"] if column in event_summary), None)
    gene_column = next((column for column in ["gene", "analysis_lost_gene", "lost_gene"] if column in event_summary), None)
    cancer_column = next((column for column in ["cancer", "project_group"] if column in event_summary), None)
    if not frequency_column or not gene_column or not cancer_column:
      matrix = pd.DataFrame(); labels = []; source = pd.DataFrame()
    else:
      source = event_summary.copy()
      source[frequency_column] = pd.to_numeric(source[frequency_column], errors="coerce")
      source = source.sort_values(frequency_column, ascending=False).head(top_n)
      matrix = source.pivot_table(index=gene_column, columns=cancer_column, values=frequency_column, aggfunc="max")
      labels = [wrap_label(value, 30) for value in matrix.index]
  return heatmap_record(item, matrix, labels, output_dir, source_dir, source, [input_path], strict, "Homozygous-deletion frequency", vmin=0, vmax=float(np.nanmax(matrix.to_numpy())) if not matrix.empty else 1)


def figure_s12(item: dict, coverage: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path) -> FigureRecord:
  return simple_bar(item, coverage, "source", "evidence_rows", output_dir, source_dir, [input_path], strict, "Evidence rows")


def figure_s13(item: dict, sensitivity: pd.DataFrame, output_dir: Path, source_dir: Path, strict: bool, input_path: Path, top_n: int = 50) -> FigureRecord:
  if not sensitivity.empty and "supportive_delta" in sensitivity:
    subset = sensitivity.sort_values([column for column in ["q_value_bh", "supportive_delta"] if column in sensitivity], ascending=[True, False][:len([column for column in ["q_value_bh", "supportive_delta"] if column in sensitivity])]).head(top_n).copy()
    subset["label"] = [wrap_label(f"{lost} → {target} | {drug} | {CANCER_LABELS.get(str(cancer), cancer)}", 42) for lost, target, drug, cancer in zip(subset["lost_gene"], subset["target_gene"], subset["drug_name"], subset["cancer"])]
  else:
    subset = pd.DataFrame(columns=["label", "supportive_delta"])
  return simple_bar(item, subset, "label", "supportive_delta", output_dir, source_dir, [input_path], strict, "Biomarker-selective response effect")


def figure_s14(item: dict, output_root: Path, output_dir: Path, source_dir: Path, strict: bool) -> FigureRecord:
  audit_files = sorted((output_root / "figures").rglob("*.layout_audit.json"))
  rows = []
  for path in audit_files:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows.append({
      "figure": payload.get("figure_id"),
      "status": payload.get("status"),
      "warnings": len(payload.get("warnings") or []),
      "axes": payload.get("axes_count"),
      "text_objects": payload.get("text_objects"),
      "width_inches": payload.get("width_inches"),
      "height_inches": payload.get("height_inches"),
      "audit_file": str(path),
    })
  source = pd.DataFrame(rows)
  set_publication_style()
  fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.5), constrained_layout=True)
  if source.empty:
    placeholder(axes[0], "Layout audits", "No prior layout-audit files were available.")
    placeholder(axes[1], "Figure dimensions", "No prior layout-audit files were available.")
  else:
    status_counts = source.groupby("status", as_index=False).agg(figures=("figure", "size"))
    axes[0].bar(status_counts["status"], status_counts["figures"])
    axes[0].set_ylabel("Figures")
    axes[0].set_title("Automated layout-audit status")
    axes[0].grid(axis="y", alpha=0.25)
    axes[1].scatter(source["width_inches"], source["height_inches"], s=35)
    axes[1].set_xlabel("Width (inches)")
    axes[1].set_ylabel("Height (inches)")
    axes[1].set_title("Rendered figure dimensions")
    axes[1].grid(alpha=0.25)
  fig.suptitle(str(item["title"]))
  return save_record(fig, item, output_dir, source_dir, source, audit_files, strict)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--config", default="config/article_assets.yaml")
  parser.add_argument("--ranking", default="results/expanded_26Q1/full/expanded_rses_onco.tsv")
  parser.add_argument("--candidates", default="data/processed/expanded_candidate_universe.tsv")
  parser.add_argument("--discovery", default="results/expanded_26Q1/discovery/all_target_dependency_screen.tsv")
  parser.add_argument("--dependency", default="results/expanded_26Q1/full/expanded_dependency_contrasts.tsv")
  parser.add_argument("--expression-compensation", default="results/expanded_26Q1/full/expanded_expression_compensation.tsv")
  parser.add_argument("--phenotype-profiles", default="results/expanded_26Q1/full/expanded_crispr_phenotype_profiles.tsv")
  parser.add_argument("--expression-profiles", default="results/expanded_26Q1/full/expanded_expression_context_profiles.tsv")
  parser.add_argument("--functional-evidence", default="data/processed/expanded_pair_functional_evidence.tsv")
  parser.add_argument("--tcga-events", default="results/expanded_26Q1/full/tcga_gene_event_summary.tsv")
  parser.add_argument("--pharmacology-coverage", default="results/expanded_26Q1/pharmacology/pharmacology_source_coverage.tsv")
  parser.add_argument("--drug-sensitivity", default="data/processed/pharmacology/drug_response_selectivity.tsv")
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument("--strict-layout", action=argparse.BooleanOptionalAction, default=True)
  args = parser.parse_args()

  config_path = resolve_path(args.config)
  config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
  registry = registry_by_id(config)
  expected = {f"Figure_S{index}" for index in range(1, 15)}
  missing = expected - set(registry)
  if missing:
    raise ValueError(f"Supplementary figure registry is incomplete: {sorted(missing)}")

  paths = {
    "ranking": resolve_path(args.ranking),
    "candidates": resolve_path(args.candidates),
    "discovery": resolve_path(args.discovery),
    "dependency": resolve_path(args.dependency),
    "expression": resolve_path(args.expression_compensation),
    "phenotype": resolve_path(args.phenotype_profiles),
    "expression_profiles": resolve_path(args.expression_profiles),
    "functional": resolve_path(args.functional_evidence),
    "tcga": resolve_path(args.tcga_events),
    "pharmacology_coverage": resolve_path(args.pharmacology_coverage),
    "sensitivity": resolve_path(args.drug_sensitivity),
  }
  ranking = pd.read_csv(paths["ranking"], sep="\t")
  candidates = pd.read_csv(paths["candidates"], sep="\t")
  discovery = read_optional(paths["discovery"])
  dependency = read_optional(paths["dependency"])
  expression = read_optional(paths["expression"])
  phenotype = read_optional(paths["phenotype"])
  expression_profiles = read_optional(paths["expression_profiles"])
  functional = read_optional(paths["functional"])
  tcga = read_optional(paths["tcga"])
  pharmacology_coverage = read_optional(paths["pharmacology_coverage"])
  sensitivity = read_optional(paths["sensitivity"])

  output_root = resolve_path(args.output_root)
  output_dir = output_root / "figures" / "supplementary"
  source_dir = output_root / "source_data" / "figures" / "supplementary"
  manifest_dir = output_root / "manifests"
  legend_dir = output_root / "manuscript_assets"
  records: list[FigureRecord] = []
  records.append(figure_s1(registry["Figure_S1"], ranking, output_dir, source_dir, args.strict_layout, paths["ranking"]))
  records.append(figure_s2(registry["Figure_S2"], candidates, output_dir, source_dir, args.strict_layout, paths["candidates"]))
  records.append(figure_s3(registry["Figure_S3"], discovery, output_dir, source_dir, args.strict_layout, paths["discovery"]))
  records.append(contrast_heatmap(registry["Figure_S4"], dependency, "delta_effect", output_dir, source_dir, args.strict_layout, paths["dependency"], "Δ CRISPR effect"))
  records.append(contrast_heatmap(registry["Figure_S5"], expression, "delta_expression", output_dir, source_dir, args.strict_layout, paths["expression"], "Δ expression"))
  records.append(profile_heatmap(registry["Figure_S6"], phenotype, "divergence", output_dir, source_dir, args.strict_layout, paths["phenotype"], "CRISPR phenotype divergence"))
  records.append(profile_heatmap(registry["Figure_S7"], expression_profiles, "divergence", output_dir, source_dir, args.strict_layout, paths["expression_profiles"], "Expression-context divergence"))
  records.append(figure_s8(registry["Figure_S8"], functional, output_dir, source_dir, args.strict_layout, paths["functional"]))
  records.append(figure_s9(registry["Figure_S9"], functional, output_dir, source_dir, args.strict_layout, paths["functional"]))
  records.append(figure_s10(registry["Figure_S10"], functional, output_dir, source_dir, args.strict_layout, paths["functional"]))
  records.append(figure_s11(registry["Figure_S11"], tcga, output_dir, source_dir, args.strict_layout, paths["tcga"]))
  records.append(figure_s12(registry["Figure_S12"], pharmacology_coverage, output_dir, source_dir, args.strict_layout, paths["pharmacology_coverage"]))
  records.append(figure_s13(registry["Figure_S13"], sensitivity, output_dir, source_dir, args.strict_layout, paths["sensitivity"]))
  records.append(figure_s14(registry["Figure_S14"], output_root, output_dir, source_dir, args.strict_layout))

  write_figure_manifest(records, manifest_dir / "supplementary_figure_manifest.tsv")
  write_legends_markdown(records, legend_dir / "supplementary_figure_legends.md")
  print(pd.DataFrame([asdict(record) for record in records])[["figure_id", "layout_status", "base_path"]].to_string(index=False))
  print(f"Wrote all supplementary figures to {output_dir}")


if __name__ == "__main__":
  main()
