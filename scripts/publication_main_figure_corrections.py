#!/usr/bin/env python3
"""Scientific and layout corrections for main Figures 1, 2, 4 and 5."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
import pandas as pd

from rses_onco.publication import set_publication_style, wrap_label
from scripts.publication_scientific_semantics import add_display_pair_columns

ROOT = Path(__file__).resolve().parents[1]
CANCER_LABELS = {"colon": "Colorectal", "stomach": "Gastric", "lung": "Lung"}
DOMAIN_ORDER = [
  "expression_context", "localization", "biochemical_structural",
  "genetic_phenotype", "interaction_network", "regulatory_network",
]
DOMAIN_LABELS = {
  "expression_context": "Expression context",
  "localization": "Subcellular localization",
  "biochemical_structural": "Biochemical / structural",
  "genetic_phenotype": "CRISPR phenotype",
  "interaction_network": "Protein interaction network",
  "regulatory_network": "Regulatory network",
}


def figure_1(module: Any, item: dict, output_dir: Path, source_dir: Path, strict: bool):
  """Compact source-to-validation workflow without embedded title or footer text."""
  set_publication_style()
  fig, axis = plt.subplots(figsize=(15.5, 5.8), constrained_layout=True)
  axis.set_axis_off()
  stages = [
    ("Candidate universe", "Curated NISEs\nParalogs\nMechanistic benchmarks\nAll-target screen"),
    ("Functional microniches", "Expression\nLocalization\nBiochemistry / structure\nCRISPR phenotypes\nSTRING + DoRothEA"),
    ("Cancer evidence", "TCGA/GDC events\nDepMap dependency\nLoss selectivity\nExpression compensation"),
    ("Experimental tractability", "Open Targets\nChEMBL\nDGIdb / MyChem\nPharos / CIViC\nDrug-response datasets"),
    ("Validation", "Biomarker-matched models\nOrthogonal perturbation\nRescue\nMechanistic assays\nIn vivo follow-up"),
  ]
  x_positions = np.linspace(0.02, 0.82, len(stages))
  width = 0.155
  height = 0.67
  rows = []
  for index, ((title, body), x) in enumerate(zip(stages, x_positions), start=1):
    box = patches.FancyBboxPatch(
      (x, 0.16), width, height, boxstyle="round,pad=0.016",
      linewidth=1.4, facecolor="white", edgecolor="0.25", transform=axis.transAxes,
    )
    axis.add_patch(box)
    axis.text(x + width / 2, 0.74, title, transform=axis.transAxes, ha="center", va="center", fontsize=11.2, fontweight="bold")
    axis.text(x + width / 2, 0.47, body, transform=axis.transAxes, ha="center", va="center", fontsize=9.4, linespacing=1.35)
    rows.append({"stage": index, "stage_label": title, "essential_sources_or_steps": body.replace("\n", ";")})
    if index < len(stages):
      axis.annotate(
        "", xy=(x + width + 0.038, 0.495), xytext=(x + width + 0.008, 0.495),
        xycoords=axis.transAxes, arrowprops={"arrowstyle": "-|>", "linewidth": 1.5},
      )
  return module.save_record(
    fig=fig, item=item, output_dir=output_dir, source_dir=source_dir,
    source_data=pd.DataFrame(rows), inputs=[ROOT / "config/article_assets.yaml"],
    strict_layout=strict,
  )


def readable_class(value: object) -> str:
  text = str(value).replace("_", " ").strip()
  replacements = {
    "NISE": "Non-homologous isofunctional enzymes",
    "collateral dependency": "Collateral dependency",
    "paralog": "Homologous paralogs",
    "pathway backup": "Pathway backup",
  }
  return replacements.get(text, text[:1].upper() + text[1:] if text else "Unclassified")


def figure_2(
  module: Any,
  item: dict,
  candidates: pd.DataFrame,
  output_dir: Path,
  source_dir: Path,
  strict: bool,
  input_path: Path,
):
  """Show the full class asymmetry and a transparent low-frequency enlargement."""
  set_publication_style()
  classes = (
    candidates.assign(
      source_class=candidates.get("source_class", pd.Series("unclassified", index=candidates.index)).fillna("unclassified").astype(str)
    ).groupby("source_class", as_index=False).agg(
      directed_candidates=("pair_id", "nunique"),
      unique_targets=("target_gene", "nunique"),
      unique_origin_genes=("lost_gene", "nunique"),
    ).sort_values("directed_candidates")
  )
  classes["class_label"] = classes["source_class"].map(readable_class)
  low_threshold = max(10, int(np.ceil(classes["directed_candidates"].quantile(0.65)))) if not classes.empty else 10
  low = classes.loc[classes["directed_candidates"] <= low_threshold].copy()

  nise = candidates.loc[candidates.get("source_class", pd.Series(index=candidates.index, dtype=object)).astype(str).eq("NISE")].copy()
  if not nise.empty and "activity" in nise:
    activities = nise.groupby(["group_id", "activity"], as_index=False).agg(
      directed_candidates=("pair_id", "nunique"),
      unique_members=("lost_gene", "nunique"),
    ).sort_values("directed_candidates")
    activities["activity_label"] = activities["group_id"].astype(str) + ": " + activities["activity"].astype(str)
  else:
    activities = pd.DataFrame(columns=["group_id", "activity", "directed_candidates", "unique_members", "activity_label"])

  fig = plt.figure(figsize=(20.0, 10.5), constrained_layout=True)
  grid = fig.add_gridspec(1, 3, width_ratios=[1.1, 0.9, 1.35])
  axes = [fig.add_subplot(grid[0, index]) for index in range(3)]

  positions = np.arange(len(classes), dtype=float)
  axes[0].barh(positions, classes["directed_candidates"])
  axes[0].set_yticks(positions, [wrap_label(value, 27) for value in classes["class_label"]])
  axes[0].set_xlabel("Directed candidate hypotheses")
  axes[0].set_title("Complete mechanistic-class distribution")
  axes[0].grid(axis="x", alpha=0.25)
  for y, value in zip(positions, classes["directed_candidates"]):
    axes[0].text(float(value), y, f" {int(value)}", va="center", fontsize=8.5)
  module.panel_label(axes[0], "A")

  low_positions = np.arange(len(low), dtype=float)
  axes[1].barh(low_positions, low["directed_candidates"])
  axes[1].set_yticks(low_positions, [wrap_label(value, 25) for value in low["class_label"]])
  axes[1].set_xlabel("Directed candidate hypotheses")
  axes[1].set_title(f"Low-frequency classes (≤{low_threshold})")
  axes[1].grid(axis="x", alpha=0.25)
  for y, value in zip(low_positions, low["directed_candidates"]):
    axes[1].text(float(value), y, f" {int(value)}", va="center", fontsize=8.5)
  module.panel_label(axes[1], "B")

  activity_positions = np.arange(len(activities), dtype=float)
  axes[2].barh(activity_positions, activities["directed_candidates"])
  axes[2].set_yticks(activity_positions, [wrap_label(value, 34) for value in activities["activity_label"]])
  axes[2].set_xlabel("Directed NISE hypotheses")
  axes[2].set_title("Complete curated NISE activity groups")
  axes[2].grid(axis="x", alpha=0.25)
  for y, value in zip(activity_positions, activities["directed_candidates"]):
    axes[2].text(float(value), y, f" {int(value)}", va="center", fontsize=8.5)
  module.panel_label(axes[2], "C")

  source = pd.concat([
    classes.assign(panel="A_complete_classes", display_threshold=np.nan),
    low.assign(panel="B_low_frequency_enlargement", display_threshold=low_threshold),
    activities.assign(panel="C_nise_activities", display_threshold=np.nan),
  ], ignore_index=True, sort=False)
  return module.save_record(
    fig=fig, item=item, output_dir=output_dir, source_dir=source_dir,
    source_data=source, inputs=[input_path], strict_layout=strict,
  )


def figure_4(
  module: Any,
  item: dict,
  ranking: pd.DataFrame,
  output_dir: Path,
  source_dir: Path,
  strict: bool,
  input_path: Path,
):
  """Use biologically valid event-frequency axes and a true point-key panel."""
  set_publication_style()
  required = {"tcga_homdel_frequency", "component_selectivity", "coverage_adjusted_rses"}
  subset = ranking.dropna(subset=list(required)).copy() if required.issubset(ranking.columns) else pd.DataFrame()
  subset = add_display_pair_columns(subset)
  fig = plt.figure(figsize=(18.5, 11.0), constrained_layout=True)
  grid = fig.add_gridspec(2, 3, height_ratios=[2.35, 1.15])
  source_frames = []
  key_rows = []
  for index, cancer in enumerate(("colon", "stomach", "lung")):
    axis = fig.add_subplot(grid[0, index])
    group = subset.loc[subset["cancer"].astype(str).eq(cancer)].copy()
    if group.empty:
      module.placeholder(axis, CANCER_LABELS[cancer], "No paired TCGA event and DepMap selectivity observations were available.")
      continue
    x = pd.to_numeric(group["tcga_homdel_frequency"], errors="coerce")
    if ((x < 0) | (x > 1)).any():
      raise ValueError(f"Biologically invalid TCGA homozygous-deletion frequency in {cancer}")
    y = pd.to_numeric(group["component_selectivity"], errors="coerce")
    axis.scatter(x, y, s=38 + 220 * pd.to_numeric(group["coverage_adjusted_rses"], errors="coerce").fillna(0), alpha=0.66, edgecolors="0.25", linewidths=0.45)
    maximum = max(0.01, float(x.max()) * 1.10 if x.notna().any() else 0.01)
    axis.set_xlim(0, min(1.0, maximum))
    if np.isclose(float(x.max()), 0.0):
      axis.text(0.5, 0.04, "All eligible frequencies are 0", transform=axis.transAxes, ha="center", va="bottom", fontsize=9)
    axis.set_xlabel("TCGA homozygous-deletion frequency")
    axis.set_ylabel("DepMap loss-selectivity component")
    axis.set_title(CANCER_LABELS[cancer])
    axis.grid(alpha=0.25)
    module.panel_label(axis, chr(ord("A") + index))
    top = group.sort_values("coverage_adjusted_rses", ascending=False).head(5).copy()
    for number, (_, row) in enumerate(top.iterrows(), start=1):
      px = float(row["tcga_homdel_frequency"])
      py = float(row["component_selectivity"])
      axis.annotate(
        str(number), (px, py), xytext=(5 + 2 * (number % 2), 5 + 4 * (number % 3)),
        textcoords="offset points", fontsize=8.5, fontweight="bold",
        bbox={"boxstyle": "circle,pad=0.18", "facecolor": "white", "edgecolor": "0.25"},
      )
      key_rows.append({
        "panel": chr(ord("A") + index), "number": number, "cancer": cancer,
        "pair_id": row.get("pair_id"), "display_pair_label": row.get("display_pair_label"),
        "tcga_homdel_frequency": px, "component_dependency": row.get("component_dependency"),
        "component_selectivity": row.get("component_selectivity"),
        "n_loss": row.get("n_loss") or row.get("tcga_homdel_n"),
        "n_intact": row.get("n_intact"), "p_value": row.get("p_value"),
        "q_value_bh": row.get("q_value_bh"), "coverage": row.get("evidence_coverage"),
        "coverage_adjusted_rses": row.get("coverage_adjusted_rses"),
      })
    source_frames.append(group.assign(panel=chr(ord("A") + index), record_type="scatter_data"))

  key_axis = fig.add_subplot(grid[1, :])
  key_axis.set_axis_off()
  if key_rows:
    rows = []
    for row in key_rows:
      rows.append([
        f"{row['panel']}{row['number']}", wrap_label(row["display_pair_label"], 30),
        CANCER_LABELS.get(row["cancer"], row["cancer"]),
        f"{row['tcga_homdel_frequency']:.4f}",
        f"{float(row['component_dependency']):.3f}" if pd.notna(row.get("component_dependency")) else "NA",
        f"{float(row['component_selectivity']):.3f}" if pd.notna(row.get("component_selectivity")) else "NA",
        f"{float(row['coverage']):.3f}" if pd.notna(row.get("coverage")) else "NA",
      ])
    table = key_axis.table(
      cellText=rows, colLabels=["ID", "Pair/context", "Cancer", "Event frequency", "Dependency", "Selectivity", "Coverage"],
      cellLoc="left", colLoc="left", loc="center", colWidths=[0.06, 0.30, 0.12, 0.13, 0.12, 0.12, 0.10],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.4)
    table.scale(1.0, 1.45)
  else:
    key_axis.text(0.5, 0.5, "No eligible high-priority points.", transform=key_axis.transAxes, ha="center", va="center")
  key_axis.set_title("D  Numbered-point key and exact plotted values", loc="left", fontsize=12.5, fontweight="bold", pad=8)

  source = pd.concat([
    *source_frames,
    pd.DataFrame(key_rows).assign(panel="D_point_key", record_type="numbered_point_key") if key_rows else pd.DataFrame(),
  ], ignore_index=True, sort=False)
  return module.save_record(
    fig=fig, item=item, output_dir=output_dir, source_dir=source_dir,
    source_data=source, inputs=[input_path], strict_layout=strict,
  )


def figure_5(
  module: Any,
  item: dict,
  ranking: pd.DataFrame,
  output_dir: Path,
  source_dir: Path,
  strict: bool,
  input_path: Path,
  top_n: int,
):
  """Render observed values while explicitly marking missing and non-eligible cells."""
  audit_path = ROOT / "article_outputs" / "tables" / "qc" / "candidate_domain_evidence_audit.tsv"
  if not audit_path.exists() or audit_path.stat().st_size == 0:
    raise RuntimeError(f"Figure 5 requires the candidate-domain audit: {audit_path}")
  audit = pd.read_csv(audit_path, sep="\t", low_memory=False)
  ranked = add_display_pair_columns(ranking)
  pair_level = (
    ranked.sort_values("functional_microniche_adjusted", ascending=False)
      .drop_duplicates(["cancer", "pair_id"]).head(top_n).copy()
  )
  keys = pair_level[["cancer", "pair_id", "display_pair_label"]].copy()
  cells = audit.loc[
    audit["domain_family"].astype(str).eq("Functional microniche")
    & audit["domain"].astype(str).isin(DOMAIN_ORDER)
  ].merge(keys, left_on=["cancer", "candidate_id"], right_on=["cancer", "pair_id"], how="inner")
  cells["domain"] = pd.Categorical(cells["domain"], categories=DOMAIN_ORDER, ordered=True)
  cells = cells.sort_values(["display_pair_label", "domain"])
  pair_order = pair_level["display_pair_label"].tolist()
  value_matrix = cells.pivot_table(index="display_pair_label", columns="domain", values="component_normalized", aggfunc="first").reindex(index=pair_order, columns=DOMAIN_ORDER)
  status_matrix = cells.pivot_table(index="display_pair_label", columns="domain", values="evidence_state", aggfunc="first").reindex(index=pair_order, columns=DOMAIN_ORDER)

  set_publication_style()
  fig, axis = plt.subplots(figsize=(14.5, max(8.5, 0.46 * len(value_matrix))), constrained_layout=True)
  masked = np.ma.masked_invalid(value_matrix.to_numpy(dtype=float))
  image = axis.imshow(masked, aspect="auto", vmin=0, vmax=1)
  axis.set_xticks(np.arange(len(DOMAIN_ORDER)), [wrap_label(DOMAIN_LABELS[domain], 18) for domain in DOMAIN_ORDER], rotation=25, ha="right")
  axis.set_yticks(np.arange(len(value_matrix)), [wrap_label(value, 34) for value in value_matrix.index])
  axis.set_xlabel("Functional-microniche domain")
  axis.set_ylabel("Prioritized directed hypothesis")
  colorbar = fig.colorbar(image, ax=axis, fraction=0.025, pad=0.02)
  colorbar.set_label("Observed specialization / divergence value")

  for y in range(len(status_matrix)):
    for x in range(len(status_matrix.columns)):
      state = str(status_matrix.iloc[y, x])
      if state == "not_eligible":
        rectangle = patches.Rectangle((x - 0.5, y - 0.5), 1, 1, facecolor="0.75", edgecolor="0.35", hatch="///", linewidth=0.5)
        axis.add_patch(rectangle)
      elif state in {"missing", "technical_failure", "insufficient_sample", "nan"}:
        rectangle = patches.Rectangle((x - 0.5, y - 0.5), 1, 1, facecolor="white", edgecolor="0.55", hatch="...", linewidth=0.5)
        axis.add_patch(rectangle)
      elif pd.notna(value_matrix.iloc[y, x]):
        axis.text(x, y, f"{float(value_matrix.iloc[y, x]):.2f}", ha="center", va="center", fontsize=7.2)

  legend_handles = [
    patches.Patch(facecolor="white", edgecolor="0.55", hatch="...", label="Evidence unavailable / insufficient"),
    patches.Patch(facecolor="0.75", edgecolor="0.35", hatch="///", label="Domain not eligible"),
  ]
  axis.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=2, frameon=False)

  cells["display_value"] = cells["component_normalized"]
  cells["display_status"] = cells["evidence_state"]
  cells["calculation_rule"] = "Observed value shown; missing and non-eligible cells are not imputed."
  cells["source_and_reason"] = cells["evidence_source"].astype(str) + " | " + cells["absence_reason"].fillna("").astype(str)
  return module.save_record(
    fig=fig, item=item, output_dir=output_dir, source_dir=source_dir,
    source_data=cells, inputs=[input_path, audit_path], strict_layout=strict,
  )
