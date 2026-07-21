#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)


def save_all(fig, stem: str) -> None:
  for suffix in ("png", "pdf", "svg"):
    fig.savefig(FIG / f"{stem}.{suffix}", dpi=300, bbox_inches="tight")
  plt.close(fig)


def workflow() -> None:
  fig, ax = plt.subplots(figsize=(16, 10))
  ax.set_xlim(0, 16); ax.set_ylim(0, 10); ax.axis("off")
  ax.set_title("RSES-Onco: from functional enzyme backups to experimentally testable synthetic-lethal pairs", fontsize=19, weight="bold", pad=20)
  ax.text(0.7, 8.95, "Relation validity is a separate gate", fontsize=13, weight="bold")
  ax.text(0.7, 8.55, "Strict NISE, homologous paralog, pathway backup and collateral deletion are retained as distinct biological mechanisms.", fontsize=11)

  source_boxes = [
    (0.5, 6.65, 2.7, 1.25, "Human functional-backup\natlas\nNISEs, paralogs and pathways"),
    (3.55, 6.65, 2.7, 1.25, "TCGA / GDC\nloss-of-function\nand GISTIC -2 events"),
    (6.60, 6.65, 2.7, 1.25, "DepMap\nCRISPR gene effect\ncopy number and expression"),
    (9.65, 6.65, 2.7, 1.25, "GTEx / matched normal\nexpression\nand safety context"),
    (12.70, 6.65, 2.7, 1.25, "Literature / SL priors\nmechanism, drugs\nand validation"),
  ]
  middle_boxes = [
    (1.0, 3.75, 4.0, 1.55, "Biomarker-positive versus negative\nloss frequency and target dependency shift\nselectivity and compensation"),
    (6.0, 3.75, 4.0, 1.55, "Coverage-aware RSES-Onco\nobserved-domain score\ncoverage and adjusted score"),
    (11.0, 3.75, 4.0, 1.55, "Cancer-specific prioritization\ncolorectal, gastric\nand lung cancer"),
  ]
  final_box = (4.65, 1.05, 6.7, 1.65, "Validation package\nisogenic perturbation and rescue\northogonal viability and pharmacology\nmechanism and in vivo")

  for i, (x, y, w, h, label) in enumerate(source_boxes):
    rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.08", linewidth=1.5, facecolor=plt.cm.Blues(0.16 + 0.07 * i), edgecolor="0.25")
    ax.add_patch(rect)
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=9.2, linespacing=1.25)
  for i, (x, y, w, h, label) in enumerate(middle_boxes):
    rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.08", linewidth=1.6, facecolor=plt.cm.Blues(0.24 + 0.08 * i), edgecolor="0.25")
    ax.add_patch(rect)
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=10.3, weight="bold" if i == 1 else "normal", linespacing=1.30)
  x, y, w, h, label = final_box
  ax.add_patch(patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.08", linewidth=1.8, facecolor=plt.cm.Blues(0.42), edgecolor="0.20"))
  ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=10.0, weight="bold", linespacing=1.22)

  for sx, sy, sw, sh, _ in source_boxes:
    target_x = 3.0 if sx < 6 else (8.0 if sx < 11 else 13.0)
    ax.annotate("", xy=(target_x, 5.30), xytext=(sx + sw / 2, sy), arrowprops=dict(arrowstyle="->", lw=1.5, color="0.25", shrinkA=5, shrinkB=6))
  ax.annotate("", xy=(6.0, 3.50), xytext=(5.0, 3.50), arrowprops=dict(arrowstyle="->", lw=2.0, color="0.25"))
  ax.annotate("", xy=(11.0, 3.50), xytext=(10.0, 3.50), arrowprops=dict(arrowstyle="->", lw=2.0, color="0.25"))
  ax.annotate("", xy=(7.2, 2.60), xytext=(8.0, 3.75), arrowprops=dict(arrowstyle="->", lw=1.8, color="0.25"))
  ax.annotate("", xy=(9.7, 2.60), xytext=(13.0, 3.75), arrowprops=dict(arrowstyle="->", lw=1.8, color="0.25"))
  ax.text(8.0, 0.55, "Computational prioritization is hypothesis generation, not clinical evidence.", fontsize=12, weight="bold", ha="center")
  save_all(fig, "Figure_1_RSES_Onco_workflow")


def catalog_summary() -> None:
  d = pd.read_csv(ROOT / "data/curated/human_nise_bonafide_2017.tsv", sep="\t")
  sizes = d.groupby(["group_id", "ec_number", "activity"]).size().reset_index(name="proteins").sort_values("proteins")
  fig, ax = plt.subplots(figsize=(14, 9))
  labels = [f"{r.ec_number}  |  {r.activity}" for r in sizes.itertuples()]
  labels = [label.replace(" / ", " /\n") if len(label) > 58 else label for label in labels]
  bars = ax.barh(range(len(sizes)), sizes.proteins)
  ax.set_yticks(range(len(sizes)), labels, fontsize=10)
  ax.set_xlabel("Proteins in the bona fide human analog activity group", fontsize=12)
  ax.set_title("Human intragenomic analogous enzymes: 70 proteins across 15 EC activities", fontsize=16, weight="bold", pad=15)
  ax.grid(axis="x", alpha=0.25)
  for bar, val in zip(bars, sizes.proteins):
    ax.text(val + 0.12, bar.get_y() + bar.get_height() / 2, str(int(val)), va="center", fontsize=10)
  ax.set_xlim(0, max(sizes.proteins) + 1.6)
  fig.subplots_adjust(left=0.35, bottom=0.11, right=0.96, top=0.91)
  fig.text(0.02, 0.02, "Source: Piergiorge et al. (2017), DOI 10.1093/gbe/evx119. Cluster labels describe reported structural groups, not a new phylogeny.", fontsize=9)
  save_all(fig, "Figure_2_human_NISE_catalog")


def cancer_priority() -> None:
  r = pd.read_csv(ROOT / "results/literature_anchored_candidates.tsv", sep="\t")
  top = r.head(15).copy()
  mat = top[["colon", "stomach", "lung"]].to_numpy(float) * top["coverage_adjusted_rses"].to_numpy()[:, None]
  fig, ax = plt.subplots(figsize=(10, 9))
  im = ax.imshow(mat, aspect="auto", vmin=0, vmax=1, cmap="viridis")
  ax.set_xticks([0, 1, 2], ["Colorectal", "Gastric", "Lung"], fontsize=12)
  y = [f"{a} -> {b}" for a, b in zip(top.lost_feature, top.target_gene)]
  ax.set_yticks(range(len(top)), y, fontsize=10)
  ax.set_title("Literature-anchored RSES-Onco priorities by target cancer", fontsize=16, weight="bold", pad=15)
  for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
      text = "-" if top.iloc[i][["colon", "stomach", "lung"][j]] == 0 else f"{mat[i,j]:.2f}"
      ax.text(j, i, text, ha="center", va="center", fontsize=9, color="white" if mat[i,j] > .55 else "black")
  cbar = fig.colorbar(im, ax=ax, pad=0.02)
  cbar.set_label("Coverage-adjusted literature prior", fontsize=11)
  fig.text(0.01, 0.01, "These are transparent literature priors. Cohort-specific TCGA/DepMap estimates must be generated with the empirical pipeline.", fontsize=9)
  save_all(fig, "Figure_3_cancer_specific_priorities")


def score_composition() -> None:
  r = pd.read_csv(ROOT / "results/literature_anchored_candidates.tsv", sep="\t").head(12).iloc[::-1]
  cols = ["component_tumor_event", "component_dependency", "component_selectivity", "component_functional_relation", "component_validation_tractability"]
  names = ["Tumor/lineage", "Genetic dependency", "Isogenic selectivity", "Functional relation", "Validation/tractability"]
  fig, ax = plt.subplots(figsize=(13, 8))
  left = np.zeros(len(r))
  for col, name in zip(cols, names):
    vals = r[col].fillna(0).to_numpy() / len(cols)
    ax.barh(range(len(r)), vals, left=left, label=name)
    left += vals
  labels = [f"{a} -> {b}" for a, b in zip(r.lost_feature, r.target_gene)]
  ax.set_yticks(range(len(r)), labels, fontsize=10)
  ax.set_xlabel("Equal-domain display (component / 5)", fontsize=12)
  ax.set_title("Evidence composition of leading literature-anchored candidates", fontsize=16, weight="bold", pad=15)
  ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=10, frameon=False)
  ax.grid(axis="x", alpha=.2)
  save_all(fig, "Figure_4_priority_score_composition")


def validation_roadmap() -> None:
  fig, ax = plt.subplots(figsize=(16, 10))
  ax.set_xlim(0, 16); ax.set_ylim(0, 10); ax.axis("off")
  ax.set_title("Experimental validation funnel for RSES-Onco candidates", fontsize=19, weight="bold", pad=18)
  stages = [
    (0.3, 6.2, 2.8, 2.0, "1. Computational gate", "Deletion or LOF biomarker\nTarget dependency\nSelectivity and coverage"),
    (3.45, 6.2, 2.8, 2.0, "2. Isogenic validation", "Restore or delete E1\nPerturb E2\nGenetic rescue"),
    (6.60, 6.2, 2.8, 2.0, "3. Orthogonal assays", "CRISPR and RNAi\nDegrader or inhibitor\nViability and apoptosis"),
    (9.75, 6.2, 2.8, 2.0, "4. Mechanism", "Flux or repair readout\nMetabolite or DNA damage\nTarget engagement"),
    (12.90, 6.2, 2.8, 2.0, "5. Translation", "Organoid or xenograft\nNormal-cell window\nBiomarker assay"),
  ]
  for x, y, w, h, title, body in stages:
    ax.add_patch(patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=.06", facecolor="white", edgecolor="0.2", linewidth=1.6))
    ax.text(x + w / 2, y + h * .72, title, ha="center", va="center", fontsize=11.5, weight="bold")
    ax.text(x + w / 2, y + h * .34, body, ha="center", va="center", fontsize=10, linespacing=1.35)
  for i in range(len(stages) - 1):
    x = stages[i][0] + stages[i][2]
    x2 = stages[i + 1][0]
    ax.annotate("", xy=(x2, 7.2), xytext=(x, 7.2), arrowprops=dict(arrowstyle="->", lw=2))
  panels = [
    (0.8, 1.45, 4.4, 3.1, "Colorectal cancer", ["MSI/MMR deficiency -> WRN", "MTAP loss -> PRMT5 or MAT2A", "HRD -> PARP1, POLQ or USP1", "NTHL1-centered repair hypotheses"]),
    (5.8, 1.45, 4.4, 3.1, "Gastric cancer", ["MSI/MMR deficiency -> WRN", "ARID1A loss -> ATR", "MTAP loss -> PRMT5 or MAT2A", "HRD-associated targets"]),
    (10.8, 1.45, 4.4, 3.1, "Lung cancer", ["SMARCA4 loss -> SMARCA2", "SMARCA4 loss -> CDK4/6", "MTAP loss -> PRMT5 or MAT2A", "ATM/HRD -> ATR, USP1 or POLQ"]),
  ]
  for x, y, w, h, title, lines in panels:
    ax.add_patch(patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=.06", facecolor=plt.cm.Oranges(.16), edgecolor="0.3", linewidth=1.4))
    ax.text(x + w / 2, y + h * .80, title, ha="center", va="center", fontsize=13, weight="bold")
    body = "\n".join(f"- {line}" for line in lines)
    ax.text(x + 0.25, y + h * .43, body, ha="left", va="center", fontsize=10.5, linespacing=1.55)
  ax.text(8.0, 0.55, "Advance only candidates with reproducible selectivity, mechanistic rescue and a plausible normal-tissue therapeutic window.", ha="center", fontsize=11.5, weight="bold")
  save_all(fig, "Figure_5_validation_roadmap")


def synthetic_demo() -> None:
  d = pd.read_csv(ROOT / "results/dependency_contrasts.tsv", sep="\t")
  if d.empty:
    return
  d = d.sort_values("delta_effect")
  fig, ax = plt.subplots(figsize=(11, 7))
  y = np.arange(len(d))
  ax.hlines(y, d.median_effect_intact, d.median_effect_loss, linewidth=2)
  ax.scatter(d.median_effect_intact, y, label="Intact group", s=60)
  ax.scatter(d.median_effect_loss, y, label="Loss group", s=60)
  labels = [f"{r.lineage}: {r.lost_gene} -> {r.target_gene}" for r in d.itertuples()]
  ax.set_yticks(y, labels, fontsize=9)
  ax.axvline(-0.5, linestyle="--", linewidth=1.2, label="Illustrative strong-dependency threshold")
  ax.set_xlabel("Median CRISPR gene effect (more negative = stronger dependency)", fontsize=11)
  ax.set_title("Software verification using synthetic planted dependencies", fontsize=16, weight="bold", pad=15)
  ax.legend(loc="upper right", ncol=1, frameon=False, fontsize=9)
  fig.subplots_adjust(bottom=0.15, left=0.24, right=0.97, top=0.90)
  fig.text(0.02, 0.025, "SYNTHETIC TEST FIXTURE - not TCGA, DepMap, patient or real cell-line data.", fontsize=10, weight="bold")
  save_all(fig, "Figure_S1_synthetic_software_verification")


if __name__ == "__main__":
  workflow(); catalog_summary(); cancer_priority(); score_composition(); validation_roadmap(); synthetic_demo()
  print(f"Wrote figures to {FIG}")
