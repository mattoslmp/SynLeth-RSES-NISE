from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import textwrap
from typing import Iterable, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import pandas as pd


FORMATS = ("png", "pdf", "svg")


@dataclass(frozen=True)
class FigureAudit:
  figure_id: str
  status: str
  warnings: tuple[str, ...]
  width_inches: float
  height_inches: float
  text_objects: int
  axes_count: int


@dataclass(frozen=True)
class FigureRecord:
  figure_id: str
  category: str
  title: str
  caption: str
  base_path: str
  source_data_path: str
  input_paths: str
  formats: str
  layout_status: str
  layout_warnings: str
  script: str


def set_publication_style() -> None:
  """Apply a journal-safe vector-friendly style with legible 100% zoom text."""
  mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "axes.titlesize": 12.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 10.5,
    "xtick.labelsize": 9.0,
    "ytick.labelsize": 9.0,
    "legend.fontsize": 9.0,
    "figure.titlesize": 15.0,
    "figure.titleweight": "bold",
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
  })


def wrap_label(value: object, width: int = 28) -> str:
  text = str(value)
  return "\n".join(
    textwrap.wrap(text, width=width, break_long_words=False)
  ) or text


def dynamic_height(
  rows: int,
  minimum: float = 5.5,
  per_row: float = 0.34,
  maximum: float = 24.0,
) -> float:
  return max(minimum, min(maximum, 2.0 + per_row * max(rows, 1)))


def add_panel_label(axis: plt.Axes, label: str) -> None:
  axis.text(
    -0.12,
    1.05,
    label,
    transform=axis.transAxes,
    fontsize=13,
    fontweight="bold",
    va="top",
    ha="left",
    clip_on=False,
  )


def placeholder(axis: plt.Axes, title: str, message: str) -> None:
  axis.set_axis_off()
  axis.set_title(title, pad=12)
  axis.text(
    0.5,
    0.5,
    wrap_label(message, 52),
    transform=axis.transAxes,
    ha="center",
    va="center",
    fontsize=10.5,
    bbox={
      "boxstyle": "round,pad=0.6",
      "facecolor": "white",
      "edgecolor": "0.65",
    },
  )


def write_source_data(frame: pd.DataFrame, path: str | Path) -> Path:
  path = Path(path)
  path.parent.mkdir(parents=True, exist_ok=True)
  frame.to_csv(path, sep="\t", index=False)
  return path


def file_sha256(path: str | Path) -> str:
  digest = hashlib.sha256()
  with Path(path).open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def _bbox_intersection_area(first: Bbox, second: Bbox) -> float:
  x0 = max(first.x0, second.x0)
  y0 = max(first.y0, second.y0)
  x1 = min(first.x1, second.x1)
  y1 = min(first.y1, second.y1)
  return max(0.0, x1 - x0) * max(0.0, y1 - y0)


def _in_view_tick_bboxes(
  axis: plt.Axes,
  orientation: str,
  renderer: object,
) -> list[Bbox]:
  """Return rendered tick-label boxes only for ticks inside current view limits."""
  if orientation == "x":
    locations = axis.get_xticks()
    labels = axis.get_xticklabels()
    low, high = sorted(axis.get_xlim())
  elif orientation == "y":
    locations = axis.get_yticks()
    labels = axis.get_yticklabels()
    low, high = sorted(axis.get_ylim())
  else:
    raise ValueError(f"Unsupported tick orientation: {orientation}")

  tolerance = max(abs(high - low), 1.0) * 1e-10
  bboxes: list[Bbox] = []
  for location, label in zip(locations, labels):
    try:
      numeric_location = float(location)
    except (TypeError, ValueError):
      continue
    if numeric_location < low - tolerance or numeric_location > high + tolerance:
      continue
    if not label.get_visible() or not label.get_text().strip():
      continue
    bbox = label.get_window_extent(renderer=renderer)
    if bbox.width > 0 and bbox.height > 0:
      bboxes.append(bbox)
  return bboxes


def _outside_figure(
  bbox: Bbox,
  figure_bbox: Bbox,
  tolerance_pixels: float,
) -> bool:
  return bool(
    bbox.x0 < figure_bbox.x0 - tolerance_pixels
    or bbox.y0 < figure_bbox.y0 - tolerance_pixels
    or bbox.x1 > figure_bbox.x1 + tolerance_pixels
    or bbox.y1 > figure_bbox.y1 + tolerance_pixels
  )


def audit_figure_layout(
  fig: plt.Figure,
  figure_id: str,
  tolerance_pixels: float = 2.0,
) -> FigureAudit:
  """Audit panel overlap, clipping and visible tick-label collisions."""
  fig.canvas.draw()
  renderer = fig.canvas.get_renderer()
  figure_bbox = fig.bbox
  warnings: list[str] = []

  axes = [axis for axis in fig.axes if axis.get_visible()]
  for index, first in enumerate(axes):
    first_box = first.get_position()
    for second in axes[index + 1:]:
      second_box = second.get_position()
      area = _bbox_intersection_area(first_box, second_box)
      if area > 1e-5:
        warnings.append(
          "axes_overlap:"
          f"{first.get_label() or index}:"
          f"{second.get_label() or axes.index(second)}"
        )

  text_objects: list[mpl.text.Text] = list(fig.texts)
  for axis in axes:
    text_objects.extend([
      axis.title,
      axis.xaxis.label,
      axis.yaxis.label,
      *axis.texts,
    ])
    legend = axis.get_legend()
    if legend is not None:
      text_objects.extend(legend.get_texts())
      legend_bbox = legend.get_window_extent(renderer=renderer)
      if _outside_figure(legend_bbox, figure_bbox, tolerance_pixels):
        warnings.append(
          f"legend_outside_figure:{axis.get_label() or axes.index(axis)}"
        )

  for text in text_objects:
    if not text.get_visible() or not text.get_text().strip():
      continue
    bbox = text.get_window_extent(renderer=renderer)
    if _outside_figure(bbox, figure_bbox, tolerance_pixels):
      warnings.append(f"text_outside_figure:{text.get_text()[:60]}")

  visible_tick_count = 0
  for axis_index, axis in enumerate(axes):
    for orientation in ("x", "y"):
      bboxes = _in_view_tick_bboxes(axis, orientation, renderer)
      visible_tick_count += len(bboxes)
      for bbox in bboxes:
        if _outside_figure(bbox, figure_bbox, tolerance_pixels):
          warnings.append(f"tick_outside_figure:{axis_index}:{orientation}")
          break
      for index, first in enumerate(bboxes):
        for second in bboxes[index + 1:]:
          if _bbox_intersection_area(first, second) > 1.0:
            warnings.append(f"tick_overlap:{axis_index}:{orientation}")
            break
        if warnings and warnings[-1] == (
          f"tick_overlap:{axis_index}:{orientation}"
        ):
          break

  warnings = sorted(set(warnings))
  width, height = fig.get_size_inches()
  return FigureAudit(
    figure_id=figure_id,
    status="pass" if not warnings else "warning",
    warnings=tuple(warnings),
    width_inches=float(width),
    height_inches=float(height),
    text_objects=len(text_objects) + visible_tick_count,
    axes_count=len(axes),
  )


def save_figure_triplet(
  fig: plt.Figure,
  base_path: str | Path,
  figure_id: str,
  strict_layout: bool = True,
) -> FigureAudit:
  base_path = Path(base_path)
  base_path.parent.mkdir(parents=True, exist_ok=True)
  if getattr(fig, "_suptitle", None) is not None:
    fig._suptitle.set_text("")
  for axis in fig.axes:
    if not axis.axison:
      axis.set_xticks([])
      axis.set_yticks([])
  audit = audit_figure_layout(fig, figure_id)
  if strict_layout and audit.warnings:
    raise RuntimeError(
      f"Layout audit failed for {figure_id}: " + "; ".join(audit.warnings)
    )
  for extension in FORMATS:
    metadata = {
      "Title": figure_id,
      "Creator": "RSES-Onco scripted publication pipeline",
    }
    save_kwargs = {
      "dpi": 600 if extension == "png" else None,
      "bbox_inches": "tight",
      "pad_inches": 0.08,
    }
    if extension in {"png", "pdf", "svg"}:
      save_kwargs["metadata"] = metadata
    fig.savefig(base_path.with_suffix(f".{extension}"), **save_kwargs)
  plt.close(fig)
  audit_path = base_path.with_suffix(".layout_audit.json")
  audit_path.write_text(
    json.dumps(asdict(audit), indent=2), encoding="utf-8"
  )
  return audit


def figure_record(
  *,
  figure_id: str,
  category: str,
  title: str,
  caption: str,
  base_path: str | Path,
  source_data_path: str | Path,
  input_paths: Sequence[str | Path],
  audit: FigureAudit,
  script: str,
) -> FigureRecord:
  return FigureRecord(
    figure_id=figure_id,
    category=category,
    title=title,
    caption=caption,
    base_path=str(base_path),
    source_data_path=str(source_data_path),
    input_paths=";".join(str(path) for path in input_paths),
    formats=";".join(FORMATS),
    layout_status=audit.status,
    layout_warnings=";".join(audit.warnings),
    script=script,
  )


def write_figure_manifest(
  records: Sequence[FigureRecord],
  path: str | Path,
) -> Path:
  path = Path(path)
  path.parent.mkdir(parents=True, exist_ok=True)
  pd.DataFrame([asdict(record) for record in records]).to_csv(
    path, sep="\t", index=False
  )
  return path


def write_legends_markdown(
  records: Sequence[FigureRecord],
  path: str | Path,
) -> Path:
  path = Path(path)
  path.parent.mkdir(parents=True, exist_ok=True)
  lines = ["# Figure legends", ""]
  for record in records:
    lines.extend([
      f"## {record.figure_id}. {record.title}",
      "",
      record.caption,
      "",
    ])
  path.write_text("\n".join(lines), encoding="utf-8")
  return path


def write_sha256_manifest(
  paths: Iterable[str | Path],
  output: str | Path,
) -> Path:
  output = Path(output)
  output.parent.mkdir(parents=True, exist_ok=True)
  rows = []
  for path in sorted({Path(item) for item in paths if Path(item).exists()}):
    if path.is_file():
      rows.append(f"{file_sha256(path)}  {path}")
  output.write_text(
    "\n".join(rows) + ("\n" if rows else ""), encoding="utf-8"
  )
  return output
