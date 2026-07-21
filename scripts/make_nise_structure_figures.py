#!/usr/bin/env python3
"""Generate the complete human NISE structural atlas from PyMOL renders.

Outputs:
- Figure 8: representative high-priority NISE structural comparisons;
- Figures S15-S29: one large, legible structural atlas per NISE activity group;
- Figure S30: AlphaFold confidence and annotation coverage;
- Figure S31: catalytic, binding and ligand-residue evidence provenance;
- Figure S32: pairwise structural evidence across NISE directions.

Molecular renders contain no text labels. Gene, accession, activity, confidence and
residue descriptions are placed in dedicated figure areas to prevent overlap.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

from rses_onco.publication import (
  figure_record,
  save_figure_triplet,
  set_publication_style,
  write_figure_manifest,
  write_legends_markdown,
  write_source_data,
)

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_optional(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    return pd.DataFrame()
  return pd.read_csv(path, sep="\t")


def load_image(path: object) -> np.ndarray | None:
  if path is None or pd.isna(path):
    return None
  candidate = resolve_path(str(path))
  if not candidate.exists():
    return None
  with Image.open(candidate) as image:
    return np.asarray(image.convert("RGB"))


def image_panel(axis: plt.Axes, image: np.ndarray | None, title: str, subtitle: str = "") -> None:
  axis.set_axis_off()
  if image is None:
    axis.text(
      0.5, 0.5, "Structure render unavailable", ha="center", va="center",
      fontsize=11, transform=axis.transAxes,
      bbox={"boxstyle": "round,pad=0.5", "facecolor": "white", "edgecolor": "0.6"},
    )
  else:
    axis.imshow(image)
  axis.set_title(title, fontsize=12.5, fontweight="bold", pad=8)
  if subtitle:
    axis.text(
      0.5, -0.025, subtitle, transform=axis.transAxes,
      ha="center", va="top", fontsize=9.5, clip_on=False,
    )


def residue_summary(annotations: pd.DataFrame, accession: str, maximum: int = 12) -> str:
  if annotations.empty:
    return "No exact UniProt-numbered functional residues available"
  subset = annotations.loc[
    annotations["uniprot_accession"].astype(str).eq(str(accession))
  ].copy()
  if subset.empty:
    return "No exact UniProt-numbered functional residues available"
  subset["label"] = subset.apply(
    lambda row: f"{row.get('residue_name') or ''}{int(row['residue_number'])} "
    f"({str(row['annotation_type']).replace('_', ' ')})",
    axis=1,
  )
  labels = subset.drop_duplicates(["residue_number", "annotation_type"])["label"].tolist()
  if len(labels) > maximum:
    labels = labels[:maximum] + [f"+{len(labels) - maximum} additional residues"]
  return "; ".join(labels)


def render_lookup(render_manifest: pd.DataFrame) -> dict[tuple[str, str], str]:
  return {
    (str(record["uniprot_accession"]), str(record["view"])): str(record["render_path"])
    for record in render_manifest.to_dict("records")
    if str(record.get("status")) == "ok"
  }


def top_nise_pairs(
  ranking: pd.DataFrame,
  proteins: pd.DataFrame,
  render_manifest: pd.DataFrame,
  number: int = 3,
) -> pd.DataFrame:
  available = set(render_manifest.loc[
    render_manifest["status"].astype(str).eq("ok"), "uniprot_accession"
  ].astype(str))
  gene_to_accession = dict(zip(proteins["gene_symbol"], proteins["uniprot_accession"]))
  subset = ranking.copy()
  if "source_class" in subset:
    subset = subset.loc[subset["source_class"].astype(str).eq("NISE")]
  lost_column = "analysis_lost_gene" if "analysis_lost_gene" in subset else "lost_gene"
  target_column = "analysis_target_gene" if "analysis_target_gene" in subset else "target_gene"
  subset["lost_gene_struct"] = subset[lost_column].astype(str)
  subset["target_gene_struct"] = subset[target_column].astype(str)
  subset["lost_accession"] = subset["lost_gene_struct"].map(gene_to_accession)
  subset["target_accession"] = subset["target_gene_struct"].map(gene_to_accession)
  subset = subset.loc[
    subset["lost_accession"].astype(str).isin(available)
    & subset["target_accession"].astype(str).isin(available)
  ]
  subset = subset.sort_values(
    ["coverage_adjusted_rses", "functional_microniche_adjusted"],
    ascending=[False, False],
    na_position="last",
  ).drop_duplicates(["lost_gene_struct", "target_gene_struct"])
  return subset.head(number)


def make_main_figure(
  ranking: pd.DataFrame,
  proteins: pd.DataFrame,
  renders: pd.DataFrame,
  annotations: pd.DataFrame,
  output_root: Path,
  strict_layout: bool,
):
  pairs = top_nise_pairs(ranking, proteins, renders, 3)
  lookup = render_lookup(renders)
  rows = max(1, len(pairs))
  fig, axes = plt.subplots(
    rows, 4, figsize=(21, max(7.5, rows * 5.2)),
    constrained_layout=True, squeeze=False,
  )
  if pairs.empty:
    for axis in axes.ravel():
      image_panel(axis, None, "No ranked NISE structural pair")
  else:
    for row_index, record in enumerate(pairs.to_dict("records")):
      lost = str(record["lost_gene_struct"])
      target = str(record["target_gene_struct"])
      lost_accession = str(record["lost_accession"])
      target_accession = str(record["target_accession"])
      score = record.get("coverage_adjusted_rses")
      score_text = f"RSES-Onco adjusted = {score:.3f}" if pd.notna(score) else "RSES-Onco adjusted unavailable"
      image_panel(
        axes[row_index, 0],
        load_image(lookup.get((lost_accession, "whole"))),
        f"{lost} | AlphaFold DB",
        f"{lost_accession}; whole structure; {score_text}",
      )
      image_panel(
        axes[row_index, 1],
        load_image(lookup.get((lost_accession, "site"))),
        f"{lost} | known functional residues",
        textwrap.fill(residue_summary(annotations, lost_accession), width=55),
      )
      image_panel(
        axes[row_index, 2],
        load_image(lookup.get((target_accession, "whole"))),
        f"{target} | AlphaFold DB",
        f"{target_accession}; whole structure",
      )
      image_panel(
        axes[row_index, 3],
        load_image(lookup.get((target_accession, "site"))),
        f"{target} | known functional residues",
        textwrap.fill(residue_summary(annotations, target_accession), width=55),
      )
  fig.suptitle(
    "Structural basis of representative human NISE vulnerabilities",
    fontsize=17, fontweight="bold",
  )
  base = output_root / "figures/main/Figure_8_human_nise_structures"
  source = output_root / "source_data/figures/structures/Figure_8_human_nise_structures.tsv"
  write_source_data(pairs, source)
  audit = save_figure_triplet(fig, base, "Figure_8", strict_layout=strict_layout)
  return figure_record(
    figure_id="Figure_8",
    category="main",
    title="Structural basis of representative human NISE vulnerabilities",
    caption=(
      "AlphaFold DB models for top-ranked directional human NISE hypotheses. "
      "Whole-structure views are colored by AlphaFold confidence; enlarged site views highlight "
      "exact UniProt-numbered catalytic, binding, metal and experimentally mapped ligand-contact residues. "
      "No ligand pose is inferred from an AlphaFold model."
    ),
    base_path=base,
    source_data_path=source,
    input_paths=["AlphaFold DB", "M-CSA", "UniProtKB", "PDBe", "expanded_rses_onco.tsv"],
    audit=audit,
    script="scripts/make_nise_structure_figures.py",
  )


def make_activity_figures(
  proteins: pd.DataFrame,
  renders: pd.DataFrame,
  annotations: pd.DataFrame,
  output_root: Path,
  strict_layout: bool,
):
  lookup = render_lookup(renders)
  records = []
  groups = sorted(proteins["group_id"].astype(str).unique())
  for offset, group_id in enumerate(groups):
    figure_number = 15 + offset
    figure_id = f"Figure_S{figure_number}"
    subset = proteins.loc[proteins["group_id"].astype(str).eq(group_id)].copy()
    subset = subset.sort_values(["structural_cluster", "gene_symbol"])
    rows = max(1, len(subset))
    fig, axes = plt.subplots(
      rows, 2, figsize=(15.5, max(7, 4.65 * rows)),
      constrained_layout=True, squeeze=False,
    )
    for row_index, protein in enumerate(subset.to_dict("records")):
      gene = str(protein["gene_symbol"])
      accession = str(protein["uniprot_accession"])
      cluster = protein["structural_cluster"]
      image_panel(
        axes[row_index, 0],
        load_image(lookup.get((accession, "whole"))),
        f"{gene} | whole AlphaFold model",
        f"{accession}; structural cluster {cluster}",
      )
      image_panel(
        axes[row_index, 1],
        load_image(lookup.get((accession, "site"))),
        f"{gene} | functional-residue view",
        textwrap.fill(residue_summary(annotations, accession, maximum=16), width=75),
      )
    activity = str(subset["activity"].iloc[0]) if not subset.empty else group_id
    ec_number = str(subset["ec_number"].iloc[0]) if not subset.empty else ""
    fig.suptitle(
      f"{group_id}: {activity} | EC {ec_number}", fontsize=16, fontweight="bold"
    )
    base = output_root / f"figures/supplementary/{figure_id}_{group_id}_structure_atlas"
    source = output_root / f"source_data/figures/structures/{figure_id}_{group_id}.tsv"
    write_source_data(subset, source)
    audit = save_figure_triplet(fig, base, figure_id, strict_layout=strict_layout)
    records.append(figure_record(
      figure_id=figure_id,
      category="supplementary",
      title=f"Human NISE structural atlas for {group_id}",
      caption=(
        f"Whole-structure and functional-residue views for every curated member of {group_id} "
        f"({activity}; EC {ec_number}). Models originate from AlphaFold DB. Residue annotations "
        "are exact UniProt-numbered evidence from M-CSA, UniProtKB and, where available, PDBe."
      ),
      base_path=base,
      source_data_path=source,
      input_paths=["human_nise_bonafide_2017.tsv", "AlphaFold DB", "structural annotations"],
      audit=audit,
      script="scripts/make_nise_structure_figures.py",
    ))
  return records


def make_coverage_figure(
  structures: pd.DataFrame,
  coverage: pd.DataFrame,
  output_root: Path,
  strict_layout: bool,
):
  first = structures.sort_values(["uniprot_accession", "fragment_index"]).drop_duplicates("uniprot_accession")
  data = coverage.merge(
    first[["uniprot_accession", "mean_plddt", "confidence_class"]],
    on="uniprot_accession", how="left",
  ).sort_values(["group_id", "gene_symbol"])
  fig, axes = plt.subplots(1, 2, figsize=(17, 20), constrained_layout=True)
  y = np.arange(len(data))
  axes[0].barh(y, data["mean_plddt"].fillna(0))
  axes[0].set_yticks(y, data["gene_symbol"])
  axes[0].invert_yaxis()
  axes[0].set_xlabel("Mean AlphaFold pLDDT")
  axes[0].set_title("Model confidence")
  axes[0].axvline(70, linestyle="--", linewidth=1)
  axes[0].axvline(90, linestyle="--", linewidth=1)
  axes[1].barh(y, data["annotated_residues"].fillna(0))
  axes[1].set_yticks(y, data["gene_symbol"])
  axes[1].invert_yaxis()
  axes[1].set_xlabel("Unique exact-numbered annotated residues")
  axes[1].set_title("Functional-residue coverage")
  fig.suptitle("AlphaFold model confidence and residue-annotation coverage", fontsize=16, fontweight="bold")
  base = output_root / "figures/supplementary/Figure_S30_structure_confidence_coverage"
  source = output_root / "source_data/figures/structures/Figure_S30_structure_confidence_coverage.tsv"
  write_source_data(data, source)
  audit = save_figure_triplet(fig, base, "Figure_S30", strict_layout=strict_layout)
  return figure_record(
    figure_id="Figure_S30", category="supplementary",
    title="AlphaFold confidence and functional-residue coverage",
    caption="Mean model pLDDT and the number of exactly mapped known functional residues for all curated human NISE proteins.",
    base_path=base, source_data_path=source,
    input_paths=["alphafold_structure_manifest.tsv", "nise_structural_annotation_coverage.tsv"],
    audit=audit, script="scripts/make_nise_structure_figures.py",
  )


def make_provenance_figure(
  proteins: pd.DataFrame,
  annotations: pd.DataFrame,
  output_root: Path,
  strict_layout: bool,
):
  categories = [
    "mcsa_catalytic", "uniprot_active_site", "uniprot_binding",
    "uniprot_metal", "uniprot_site", "pdbe_ligand", "drug_binding", "curated_user",
  ]
  if annotations.empty:
    matrix = pd.DataFrame(0, index=proteins["gene_symbol"], columns=categories)
  else:
    matrix = annotations.pivot_table(
      index="gene_symbol", columns="annotation_type", values="residue_number",
      aggfunc=lambda values: pd.Series(values).nunique(), fill_value=0,
    ).reindex(columns=categories, fill_value=0)
    matrix = proteins[["group_id", "gene_symbol"]].drop_duplicates().set_index("gene_symbol").join(matrix, how="left").fillna(0)
    matrix = matrix.sort_values(["group_id", "gene_symbol"]).drop(columns="group_id")
  fig, axis = plt.subplots(figsize=(14, max(14, 0.34 * len(matrix) + 3)), constrained_layout=True)
  image = axis.imshow(matrix.to_numpy(dtype=float), aspect="auto")
  axis.set_yticks(np.arange(len(matrix)), matrix.index)
  axis.set_xticks(np.arange(len(matrix.columns)), [value.replace("_", " ").title() for value in matrix.columns], rotation=35, ha="right")
  axis.set_title("Provenance of known functional residues mapped onto AlphaFold models")
  colorbar = fig.colorbar(image, ax=axis, fraction=0.035, pad=0.02)
  colorbar.set_label("Unique residues")
  base = output_root / "figures/supplementary/Figure_S31_structural_residue_provenance"
  source = output_root / "source_data/figures/structures/Figure_S31_structural_residue_provenance.tsv"
  write_source_data(matrix.reset_index(), source)
  audit = save_figure_triplet(fig, base, "Figure_S31", strict_layout=strict_layout)
  return figure_record(
    figure_id="Figure_S31", category="supplementary",
    title="Functional-residue evidence provenance",
    caption="Counts of exact UniProt-numbered catalytic, binding, metal, site and experimentally mapped ligand-contact residues by evidence source.",
    base_path=base, source_data_path=source,
    input_paths=["nise_structural_residue_annotations.tsv"],
    audit=audit, script="scripts/make_nise_structure_figures.py",
  )


def make_pair_summary_figure(
  ranking: pd.DataFrame,
  proteins: pd.DataFrame,
  coverage: pd.DataFrame,
  output_root: Path,
  strict_layout: bool,
):
  gene_cluster = dict(zip(proteins["gene_symbol"], proteins["structural_cluster"]))
  gene_annotation = dict(zip(coverage["gene_symbol"], coverage["annotated_residues"]))
  subset = ranking.copy()
  if "source_class" in subset:
    subset = subset.loc[subset["source_class"].astype(str).eq("NISE")]
  lost_col = "analysis_lost_gene" if "analysis_lost_gene" in subset else "lost_gene"
  target_col = "analysis_target_gene" if "analysis_target_gene" in subset else "target_gene"
  subset["lost_gene_struct"] = subset[lost_col].astype(str)
  subset["target_gene_struct"] = subset[target_col].astype(str)
  subset["different_structural_cluster"] = subset.apply(
    lambda row: float(gene_cluster.get(row["lost_gene_struct"]) != gene_cluster.get(row["target_gene_struct"]))
    if row["lost_gene_struct"] in gene_cluster and row["target_gene_struct"] in gene_cluster else np.nan,
    axis=1,
  )
  subset["mapped_residue_support"] = subset.apply(
    lambda row: min(
      int(gene_annotation.get(row["lost_gene_struct"], 0)),
      int(gene_annotation.get(row["target_gene_struct"], 0)),
    ), axis=1,
  )
  subset = subset.drop_duplicates(["lost_gene_struct", "target_gene_struct", "cancer"])
  fig, axis = plt.subplots(figsize=(11, 8), constrained_layout=True)
  for cancer, group in subset.groupby("cancer"):
    axis.scatter(
      group["mapped_residue_support"], group["coverage_adjusted_rses"],
      s=45 + 80 * group["different_structural_cluster"].fillna(0),
      alpha=0.75, label=str(cancer).title(),
    )
  axis.set_xlabel("Minimum mapped functional residues across the NISE pair")
  axis.set_ylabel("Coverage-adjusted RSES-Onco")
  axis.set_title("Pairwise structural evidence across directed human NISE hypotheses")
  axis.legend(frameon=False)
  axis.grid(alpha=0.25)
  base = output_root / "figures/supplementary/Figure_S32_pairwise_structural_evidence"
  source = output_root / "source_data/figures/structures/Figure_S32_pairwise_structural_evidence.tsv"
  write_source_data(subset, source)
  audit = save_figure_triplet(fig, base, "Figure_S32", strict_layout=strict_layout)
  return figure_record(
    figure_id="Figure_S32", category="supplementary",
    title="Pairwise structural evidence across human NISE hypotheses",
    caption="Relationship between mapped functional-residue support, structural-cluster difference and the cancer-specific RSES-Onco score.",
    base_path=base, source_data_path=source,
    input_paths=["expanded_rses_onco.tsv", "human_nise_bonafide_2017.tsv", "nise_structural_annotation_coverage.tsv"],
    audit=audit, script="scripts/make_nise_structure_figures.py",
  )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--ranking", default="results/expanded_26Q1/full/expanded_rses_onco.tsv")
  parser.add_argument("--proteins", default="data/curated/human_nise_bonafide_2017.tsv")
  parser.add_argument("--structure-manifest", default="data/processed/structures/alphafold_structure_manifest.tsv")
  parser.add_argument("--render-manifest", default="data/processed/structures/nise_structure_render_manifest.tsv")
  parser.add_argument("--annotations", default="data/processed/structures/nise_structural_residue_annotations.tsv")
  parser.add_argument("--coverage", default="data/processed/structures/nise_structural_annotation_coverage.tsv")
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument("--strict-layout", action=argparse.BooleanOptionalAction, default=True)
  args = parser.parse_args()

  set_publication_style()
  ranking = pd.read_csv(resolve_path(args.ranking), sep="\t")
  proteins = pd.read_csv(resolve_path(args.proteins), sep="\t")
  structures = pd.read_csv(resolve_path(args.structure_manifest), sep="\t")
  renders = pd.read_csv(resolve_path(args.render_manifest), sep="\t")
  annotations = read_optional(resolve_path(args.annotations))
  coverage = pd.read_csv(resolve_path(args.coverage), sep="\t")
  output_root = resolve_path(args.output_root)

  main_record = make_main_figure(
    ranking, proteins, renders, annotations, output_root, args.strict_layout
  )
  supplementary_records = make_activity_figures(
    proteins, renders, annotations, output_root, args.strict_layout
  )
  supplementary_records.extend([
    make_coverage_figure(structures, coverage, output_root, args.strict_layout),
    make_provenance_figure(proteins, annotations, output_root, args.strict_layout),
    make_pair_summary_figure(ranking, proteins, coverage, output_root, args.strict_layout),
  ])

  manifest_dir = output_root / "manifests"
  manuscript_dir = output_root / "manuscript_assets"
  write_figure_manifest([main_record], manifest_dir / "structural_main_figure_manifest.tsv")
  write_figure_manifest(supplementary_records, manifest_dir / "structural_supplementary_figure_manifest.tsv")
  write_legends_markdown([main_record], manuscript_dir / "structural_main_figure_legends.md")
  write_legends_markdown(supplementary_records, manuscript_dir / "structural_supplementary_figure_legends.md")
  print(f"Generated 1 main and {len(supplementary_records)} supplementary structural figures")


if __name__ == "__main__":
  main()
