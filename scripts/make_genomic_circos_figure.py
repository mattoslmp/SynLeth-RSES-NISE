#!/usr/bin/env python3
"""Generate Supplementary Figure S70: genomic Circos for NISE and paralog hypotheses."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))
if str(ROOT / "src") not in sys.path:
  sys.path.insert(0, str(ROOT / "src"))

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd
import yaml

from rses_onco.circos import GenomeLayout
from rses_onco.publication import (
  set_publication_style,
  write_figure_manifest,
  write_legends_markdown,
)
from scripts.publication_audit_figures import panel_label, save_record

SCRIPT = "scripts/make_genomic_circos_figure.py"


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_required(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(
      f"Mandatory Circos table is absent or empty: {path}"
    )
  return pd.read_csv(path, sep="\t", low_memory=False)


def polar_point(theta: float, radius: float) -> tuple[float, float]:
  return radius * np.cos(theta), radius * np.sin(theta)


def draw_chromosomes(
  axis: plt.Axes,
  layout: GenomeLayout,
  radius: float = 1.0,
) -> None:
  colors = ("#D9DEE7", "#AEB8C6")
  for index, arc in enumerate(layout.arcs.values()):
    axis.add_patch(Wedge(
      (0, 0),
      radius,
      np.degrees(arc.theta_end),
      np.degrees(arc.theta_start),
      width=0.055,
      facecolor=colors[index % 2],
      edgecolor="#FFFFFF",
      linewidth=0.7,
    ))
    x, y = polar_point(arc.theta_mid, radius + 0.085)
    axis.text(
      x,
      y,
      arc.chromosome,
      ha="center",
      va="center",
      fontsize=7.5,
      fontweight="bold",
    )


def draw_links(
  axis: plt.Axes,
  layout: GenomeLayout,
  links: pd.DataFrame,
  radius: float,
) -> None:
  ordered = links.sort_values(
    "maximum_coverage_adjusted_rses",
    ascending=True,
  )
  for row in ordered.to_dict("records"):
    theta_a = layout.theta(
      row["lost_chromosome"],
      row["lost_position"],
    )
    theta_b = layout.theta(
      row["target_chromosome"],
      row["target_position"],
    )
    p0 = polar_point(theta_a, radius)
    p3 = polar_point(theta_b, radius)
    p1 = (p0[0] * 0.16, p0[1] * 0.16)
    p2 = (p3[0] * 0.16, p3[1] * 0.16)
    path = MplPath(
      [p0, p1, p2, p3],
      [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
      ],
    )
    axis.add_patch(PathPatch(
      path,
      facecolor="none",
      edgecolor=str(row.get("link_color") or "#777777"),
      linewidth=float(row.get("link_width") or 0.4),
      alpha=float(row.get("link_alpha") or 0.08),
      zorder=1,
    ))


def draw_gene_ticks(
  axis: plt.Axes,
  layout: GenomeLayout,
  coordinates: pd.DataFrame,
) -> None:
  color_map = {
    "NISE": "#C62828",
    "homologous_paralog": "#111111",
    "NISE_and_paralog": "#6A1B9A",
  }
  for row in coordinates.to_dict("records"):
    theta = layout.theta(
      row["chromosome"],
      row["genomic_position"],
    )
    x0, y0 = polar_point(theta, 1.005)
    x1, y1 = polar_point(theta, 1.038)
    axis.plot(
      [x0, x1],
      [y0, y1],
      color=color_map.get(
        str(row.get("gene_class")),
        "#777777",
      ),
      linewidth=0.7,
      alpha=0.85,
      zorder=8,
    )


def top_gene_labels(
  axis: plt.Axes,
  layout: GenomeLayout,
  coordinates: pd.DataFrame,
  ring_values: pd.DataFrame,
  maximum: int = 26,
) -> None:
  scores = ring_values.loc[
    ring_values["source_column"].eq("coverage_adjusted_rses")
  ].copy()
  scores["value"] = pd.to_numeric(scores["value"], errors="coerce")
  selected = (
    scores.sort_values("value", ascending=False)
    .dropna(subset=["value"])
    .head(maximum)
  )
  lookup = coordinates.set_index("gene")
  for row in selected.to_dict("records"):
    gene = str(row["gene"])
    if gene not in lookup.index:
      continue
    coordinate = lookup.loc[gene]
    theta = layout.theta(
      coordinate["chromosome"],
      coordinate["genomic_position"],
    )
    radius = 1.13
    x, y = polar_point(theta, radius)
    rotation = np.degrees(theta)
    alignment = "left"
    if np.cos(theta) < 0:
      rotation += 180
      alignment = "right"
    axis.text(
      x,
      y,
      gene,
      rotation=rotation,
      rotation_mode="anchor",
      ha=alignment,
      va="center",
      fontsize=6.2,
    )


def draw_tracks(
  axis: plt.Axes,
  layout: GenomeLayout,
  coordinates: pd.DataFrame,
  ring_values: pd.DataFrame,
  tracks: pd.DataFrame,
  *,
  outer_radius: float = 0.94,
  inner_radius: float = 0.37,
) -> None:
  ordered = tracks.sort_values("ring_order").reset_index(drop=True)
  radii = np.linspace(outer_radius, inner_radius, len(ordered))
  coordinate_lookup = coordinates.set_index("gene")
  normalize = Normalize(vmin=0, vmax=1)
  cmap = plt.get_cmap("viridis")
  for radius, track in zip(radii, ordered.to_dict("records")):
    circle = plt.Circle(
      (0, 0),
      radius,
      fill=False,
      edgecolor="#E0E0E0",
      linewidth=0.35,
      zorder=0,
    )
    axis.add_patch(circle)
    subset = ring_values.loc[
      ring_values["track_id"].eq(track["track_id"])
    ]
    for row in subset.to_dict("records"):
      gene = str(row["gene"])
      if gene not in coordinate_lookup.index:
        continue
      coordinate = coordinate_lookup.loc[gene]
      theta = layout.theta(
        coordinate["chromosome"],
        coordinate["genomic_position"],
      )
      x, y = polar_point(theta, radius)
      value = pd.to_numeric(
        pd.Series([row.get("value")]),
        errors="coerce",
      ).iloc[0]
      if pd.isna(value):
        axis.scatter(
          [x],
          [y],
          s=4.5,
          facecolors="none",
          edgecolors="#BDBDBD",
          linewidths=0.25,
          zorder=5,
        )
      else:
        axis.scatter(
          [x],
          [y],
          s=5.0,
          color=cmap(normalize(float(value))),
          linewidths=0,
          zorder=6,
        )
    lx, ly = polar_point(np.deg2rad(218), radius)
    axis.text(
      lx,
      ly,
      str(track["track_id"]),
      fontsize=5.4,
      ha="center",
      va="center",
      bbox={
        "facecolor": "white",
        "edgecolor": "none",
        "pad": 0.2,
      },
    )


def draw_panel(
  axis: plt.Axes,
  coordinates: pd.DataFrame,
  links: pd.DataFrame,
  ring_values: pd.DataFrame,
  tracks: pd.DataFrame,
  panel: str,
) -> None:
  layout = GenomeLayout()
  axis.set_aspect("equal")
  axis.set_xlim(-1.27, 1.27)
  axis.set_ylim(-1.27, 1.27)
  axis.set_xticks([])
  axis.set_yticks([])
  axis.set_axis_off()
  draw_chromosomes(axis, layout)
  draw_gene_ticks(axis, layout, coordinates)
  selected_tracks = tracks.loc[tracks["panel"].eq(panel)].copy()
  inner = 0.36 if len(selected_tracks) <= 10 else 0.30
  draw_tracks(
    axis,
    layout,
    coordinates,
    ring_values,
    selected_tracks,
    inner_radius=inner,
  )
  draw_links(
    axis,
    layout,
    links,
    radius=max(0.22, inner - 0.025),
  )
  top_gene_labels(axis, layout, coordinates, ring_values)
  axis.text(
    0,
    0.05,
    f"{len(coordinates):,} genes",
    ha="center",
    va="center",
    fontsize=9,
    fontweight="bold",
  )
  axis.text(
    0,
    -0.03,
    f"{len(links):,} NISE/paralog links",
    ha="center",
    va="center",
    fontsize=8,
  )
  axis.text(
    0,
    -0.11,
    "red = NISE; black = homologous paralog",
    ha="center",
    va="center",
    fontsize=7,
  )


def legend_panel(axis: plt.Axes, tracks: pd.DataFrame) -> None:
  axis.set_xticks([])
  axis.set_yticks([])
  axis.set_axis_off()
  y = 0.98
  axis.text(
    0.0,
    y,
    "Ring key",
    fontsize=12,
    fontweight="bold",
    va="top",
    transform=axis.transAxes,
  )
  y -= 0.045
  for panel, title in (
    ("A", "Panel A — top-level RSES-Onco"),
    ("B", "Panel B — microniche and internal layers"),
  ):
    axis.text(
      0.0,
      y,
      title,
      fontsize=9.2,
      fontweight="bold",
      va="top",
      transform=axis.transAxes,
    )
    y -= 0.032
    for row in tracks.loc[
      tracks["panel"].eq(panel)
    ].sort_values("ring_order").to_dict("records"):
      axis.text(
        0.0,
        y,
        f"{row['track_id']}  {row['track_label']}",
        fontsize=7.1,
        va="top",
        transform=axis.transAxes,
      )
      y -= 0.026
    y -= 0.015
  y -= 0.01
  axis.text(
    0.0,
    y,
    "Rendering rules",
    fontsize=9.2,
    fontweight="bold",
    va="top",
    transform=axis.transAxes,
  )
  y -= 0.033
  rules = [
    "• every eligible NISE/paralog gene is a genomic tick",
    "• every coordinate-complete pair is a chord",
    "• ring values are maxima across pair × cancer rows",
    "• hollow markers are missing/non-eligible, never zero",
    "• chromosome positions use Ensembl/GRCh38 coordinates",
    "• exact input tables are Supplementary Tables S45–S52",
  ]
  for rule in rules:
    axis.text(
      0.0,
      y,
      rule,
      fontsize=7.1,
      va="top",
      transform=axis.transAxes,
    )
    y -= 0.03


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--config", default="config/article_assets.yaml")
  parser.add_argument(
    "--coordinates",
    default=(
      "data/processed/circos/"
      "genomic_circos_gene_coordinates.tsv"
    ),
  )
  parser.add_argument(
    "--links",
    default=(
      "data/processed/circos/genomic_circos_pair_links.tsv"
    ),
  )
  parser.add_argument(
    "--ring-values",
    default=(
      "data/processed/circos/genomic_circos_ring_values.tsv"
    ),
  )
  parser.add_argument(
    "--tracks",
    default=(
      "data/processed/circos/"
      "genomic_circos_track_definitions.tsv"
    ),
  )
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument(
    "--strict-layout",
    action=argparse.BooleanOptionalAction,
    default=True,
  )
  args = parser.parse_args()

  paths = {
    name: resolve(value)
    for name, value in {
      "coordinates": args.coordinates,
      "links": args.links,
      "ring_values": args.ring_values,
      "tracks": args.tracks,
    }.items()
  }
  coordinates = read_required(paths["coordinates"])
  links = read_required(paths["links"])
  ring_values = read_required(paths["ring_values"])
  tracks = read_required(paths["tracks"])
  config = yaml.safe_load(
    resolve(args.config).read_text(encoding="utf-8")
  ) or {}
  registry = {
    str(item["id"]): item
    for item in config.get("supplementary_figures", [])
  }
  item = registry.get("Figure_S70")
  if item is None:
    raise RuntimeError(
      "Figure_S70 is not registered in config/article_assets.yaml"
    )

  set_publication_style()
  fig = plt.figure(figsize=(27.0, 13.2), constrained_layout=True)
  grid = fig.add_gridspec(
    1,
    3,
    width_ratios=[1.0, 1.0, 0.62],
  )
  axis_a = fig.add_subplot(grid[0, 0])
  axis_b = fig.add_subplot(grid[0, 1])
  axis_legend = fig.add_subplot(grid[0, 2])
  draw_panel(axis_a, coordinates, links, ring_values, tracks, "A")
  draw_panel(axis_b, coordinates, links, ring_values, tracks, "B")
  legend_panel(axis_legend, tracks)
  panel_label(axis_a, "A")
  panel_label(axis_b, "B")

  source = pd.concat([
    coordinates.assign(record_type="gene_coordinate"),
    links.assign(record_type="pair_link"),
    ring_values.assign(record_type="ring_value"),
    tracks.assign(record_type="track_definition"),
  ], ignore_index=True, sort=False)
  output_root = resolve(args.output_root)
  record = save_record(
    fig=fig,
    figure_id="Figure_S70",
    file_name=str(item["file"]),
    title=str(item["title"]),
    caption=str(item.get("caption") or item["title"]),
    output_root=output_root,
    source=source,
    inputs=list(paths.values()),
    script=SCRIPT,
    strict=args.strict_layout,
  )
  write_figure_manifest(
    [record],
    output_root / "manifests/genomic_circos_figure_manifest.tsv",
  )
  write_legends_markdown(
    [record],
    output_root
    / "manuscript_assets/genomic_circos_figure_legend.md",
  )
  summary = pd.DataFrame([asdict(record)])
  if (
    len(summary) != 1
    or summary.iloc[0]["figure_id"] != "Figure_S70"
  ):
    raise RuntimeError("Figure S70 manifest generation failed")
  print(f"Generated Figure_S70: {record.layout_status}")


if __name__ == "__main__":
  main()
