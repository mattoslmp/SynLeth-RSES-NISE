#!/usr/bin/env python3
"""Build genomic Circos inputs for every NISE and homologous-paralog hypothesis."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
  sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd

from rses_onco.circos import (
  CHROMOSOME_LENGTHS_GRCH38,
  normalize_chromosome,
)
from rses_onco.depmap import (
  cancer_model_ids,
  detect_model_id_column,
  normalize_model_id_column,
)
from rses_onco.utils import canonical_gene_name

CANCERS = ("colon", "stomach", "lung")

TRACKS = [
  ("A01", "Coverage-adjusted RSES", "coverage_adjusted_rses", "A", "RSES-Onco", "global", 1),
  ("A02", "Evidence coverage", "evidence_coverage", "A", "RSES-Onco", "coverage", 2),
  ("A03", "Tumor event", "component_tumor_event", "A", "RSES-Onco", "tumor_event", 3),
  ("A04", "Dependency", "component_dependency", "A", "RSES-Onco", "dependency", 4),
  ("A05", "Selectivity", "component_selectivity", "A", "RSES-Onco", "selectivity", 5),
  ("A06", "Expression compensation", "component_expression_compensation", "A", "RSES-Onco", "expression_compensation", 6),
  ("A07", "Functional relation", "component_functional_relation", "A", "RSES-Onco", "functional_relation", 7),
  ("A08", "Functional microniche", "component_functional_microniche", "A", "RSES-Onco", "functional_microniche", 8),
  ("A09", "Validation and tractability", "component_validation_tractability", "A", "RSES-Onco", "validation_tractability", 9),
  ("B01", "Expression context", "microniche_expression_context", "B", "Functional microniche", "expression_context", 1),
  ("B02", "Localization", "microniche_localization", "B", "Functional microniche", "localization", 2),
  ("B03", "Biochemical/structural", "microniche_biochemical_structural", "B", "Functional microniche", "biochemical_structural", 3),
  ("B04", "Genetic phenotype", "microniche_genetic_phenotype", "B", "Functional microniche", "genetic_phenotype", 4),
  ("B05", "Interaction network", "microniche_interaction_network", "B", "Functional microniche", "interaction_network", 5),
  ("B06", "Regulatory network", "microniche_regulatory_network", "B", "Functional microniche", "regulatory_network", 6),
  ("B07", "Pairwise expression context", "pairwise_expression_context", "B", "Expression context", "pairwise_expression_context", 7),
  ("B08", "WGCNA expression network", "wgcna_expression_network", "B", "Expression context", "wgcna_expression_network", 8),
  ("B09", "DoRothEA TF associations", "regulatory_tf_association_divergence", "B", "Regulatory network", "tf_association", 9),
  ("B10", "TF expression profiles", "regulatory_tf_expression_profile_divergence", "B", "Regulatory network", "tf_expression", 10),
  ("B11", "JASPAR/FIMO promoter motifs", "regulatory_promoter_motif_divergence", "B", "Regulatory network", "promoter_motif", 11),
  ("B12", "Promoter methylation", "component_promoter_methylation_context", "B", "Regulatory network", "promoter_methylation", 12),
  ("B13", "Microniche coverage", "functional_microniche_coverage", "B", "Coverage", "microniche_coverage", 13),
  ("B14", "Expression-context coverage", "expression_context_subcoverage", "B", "Coverage", "expression_context_coverage", 14),
  ("B15", "Regulatory-network coverage", "regulatory_network_subcoverage", "B", "Coverage", "regulatory_coverage", 15),
  ("B16", "Methylation coverage", "methylation_coverage", "B", "Coverage", "methylation_coverage", 16),
]


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)


def sha256(path: Path) -> str:
  digest = hashlib.sha256()
  with path.open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def relation_flags(frame: pd.DataFrame) -> pd.DataFrame:
  text = (
    frame.get("source_class", "").fillna("").astype(str)
    + " "
    + frame.get("relation_type", "").fillna("").astype(str)
    + " "
    + frame.get("ensembl_homology_type", "").fillna("").astype(str)
  ).str.casefold()
  result = frame.copy()
  result["is_nise"] = text.str.contains("nise", regex=False)
  result["is_paralog"] = text.str.contains("paralog|homolog", regex=True)
  result["pair_class"] = np.select(
    [
      result["is_nise"] & result["is_paralog"],
      result["is_nise"],
      result["is_paralog"],
    ],
    ["NISE_and_paralog", "NISE", "homologous_paralog"],
    default="other",
  )
  return result


def read_expression_subset(path: Path, genes: set[str]) -> pd.DataFrame:
  model_column, _ = detect_model_id_column(path, "expression")
  header = pd.read_csv(path, nrows=0)
  selected: list[str] = [model_column]
  rename: dict[str, str] = {model_column: "ModelID"}
  for column in header.columns:
    if str(column) == model_column:
      continue
    gene = canonical_gene_name(column)
    if gene in genes:
      selected.append(str(column))
      rename[str(column)] = gene
  frame = pd.read_csv(path, usecols=selected, low_memory=False)
  frame = frame.rename(columns=rename)
  frame = normalize_model_id_column(frame, "expression")
  gene_columns = [column for column in frame.columns if column != "ModelID"]
  if len(gene_columns) != len(set(gene_columns)):
    numeric = frame[gene_columns].apply(pd.to_numeric, errors="coerce")
    numeric.columns = gene_columns
    collapsed = numeric.T.groupby(level=0).median().T
    frame = pd.concat([frame[["ModelID"]], collapsed], axis=1)
  return frame


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--promoters",
    default="data/raw/regulatory/ensembl_promoters.tsv",
  )
  parser.add_argument(
    "--expression",
    default=(
      "data/raw/depmap/"
      "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
    ),
  )
  parser.add_argument("--models", default="data/raw/depmap/Model.csv")
  parser.add_argument("--output-dir", default="data/processed/circos")
  args = parser.parse_args()

  paths = {
    name: resolve(value)
    for name, value in {
      "ranking": args.ranking,
      "candidates": args.candidates,
      "promoters": args.promoters,
      "expression": args.expression,
      "models": args.models,
    }.items()
  }
  for name, path in paths.items():
    if not path.exists() or path.stat().st_size == 0:
      raise FileNotFoundError(
        f"Mandatory Circos input is absent or empty: {name}={path}"
      )

  ranking = pd.read_csv(paths["ranking"], sep="\t", low_memory=False)
  candidates = relation_flags(
    pd.read_csv(paths["candidates"], sep="\t", low_memory=False)
  )
  selected = candidates.loc[
    (candidates["is_nise"] | candidates["is_paralog"])
    & candidates["lost_gene"].notna()
    & candidates["target_gene"].notna()
  ].copy()
  selected["lost_gene"] = selected["lost_gene"].map(canonical_gene_name)
  selected["target_gene"] = selected["target_gene"].map(
    canonical_gene_name
  )
  selected = selected.loc[
    selected["lost_gene"].ne("") & selected["target_gene"].ne("")
  ].drop_duplicates("pair_id")
  if selected.empty:
    raise RuntimeError(
      "No simple NISE or homologous-paralog pairs were available"
    )

  genes = sorted(set(selected["lost_gene"]) | set(selected["target_gene"]))
  promoters = pd.read_csv(paths["promoters"], sep="\t", low_memory=False)
  promoters["gene"] = promoters["gene"].map(canonical_gene_name)
  promoters["chromosome"] = promoters["chromosome"].map(
    normalize_chromosome
  )
  promoters["genomic_position"] = pd.to_numeric(
    promoters.get("gene_midpoint", promoters.get("tss")),
    errors="coerce",
  )
  coordinate_columns = [
    column
    for column in (
      "gene",
      "ensembl_gene_id",
      "canonical_transcript_id",
      "assembly",
      "chromosome",
      "strand",
      "gene_start",
      "gene_end",
      "gene_midpoint",
      "tss",
      "promoter_start",
      "promoter_end",
      "status",
      "source",
      "source_url",
      "accessed_at_utc",
    )
    if column in promoters
  ]
  coordinates = promoters[
    coordinate_columns + ["genomic_position"]
  ].copy()
  coordinates = coordinates.loc[
    coordinates["gene"].isin(genes)
  ].drop_duplicates("gene")

  class_rows = []
  for gene in genes:
    membership = selected.loc[
      selected["lost_gene"].eq(gene)
      | selected["target_gene"].eq(gene)
    ]
    nise = bool(membership["is_nise"].any())
    paralog = bool(membership["is_paralog"].any())
    gene_class = (
      "NISE_and_paralog"
      if nise and paralog
      else "NISE"
      if nise
      else "homologous_paralog"
    )
    class_rows.append({
      "gene": gene,
      "gene_class": gene_class,
      "nise_pair_count": int(membership["is_nise"].sum()),
      "paralog_pair_count": int(membership["is_paralog"].sum()),
      "total_pair_count": int(len(membership)),
    })
  coordinates = pd.DataFrame(class_rows).merge(
    coordinates,
    on="gene",
    how="left",
  )
  coordinates["coordinate_status"] = np.where(
    coordinates["chromosome"].isin(CHROMOSOME_LENGTHS_GRCH38)
    & coordinates["genomic_position"].notna(),
    "available",
    "missing_or_noncanonical",
  )
  missing_coordinates = coordinates.loc[
    ~coordinates["coordinate_status"].eq("available"),
    "gene",
  ].tolist()
  if missing_coordinates:
    raise RuntimeError(
      "Every NISE/paralog gene must have a canonical genomic coordinate; "
      "missing: "
      + ", ".join(missing_coordinates[:50])
    )

  rank_selected = ranking.loc[
    ranking["pair_id"].astype(str).isin(
      set(selected["pair_id"].astype(str))
    )
  ].copy()
  selected_metadata = selected[[
    "pair_id",
    "lost_gene",
    "target_gene",
    "pair_class",
    "source_class",
    "relation_type",
  ]]
  rank_selected = rank_selected.drop(
    columns=[
      column
      for column in (
        "lost_gene",
        "target_gene",
        "source_class",
        "relation_type",
      )
      if column in rank_selected
    ],
    errors="ignore",
  ).merge(selected_metadata, on="pair_id", how="inner")

  coord_lookup = coordinates.set_index("gene")
  link_rows = []
  for (
    pair_id,
    lost,
    target,
    pair_class,
  ), group in rank_selected.groupby(
    ["pair_id", "lost_gene", "target_gene", "pair_class"],
    dropna=False,
  ):
    lost_coord = coord_lookup.loc[lost]
    target_coord = coord_lookup.loc[target]
    scores = pd.to_numeric(
      group.get("coverage_adjusted_rses"),
      errors="coerce",
    )
    coverage = pd.to_numeric(
      group.get("evidence_coverage"),
      errors="coerce",
    )
    max_score = float(scores.max()) if scores.notna().any() else np.nan
    link_rows.append({
      "pair_id": pair_id,
      "lost_gene": lost,
      "target_gene": target,
      "pair_class": pair_class,
      "lost_chromosome": lost_coord["chromosome"],
      "lost_position": lost_coord["genomic_position"],
      "target_chromosome": target_coord["chromosome"],
      "target_position": target_coord["genomic_position"],
      "cancers": (
        ";".join(sorted(set(group["cancer"].dropna().astype(str))))
        if "cancer" in group
        else ""
      ),
      "maximum_coverage_adjusted_rses": max_score,
      "median_coverage_adjusted_rses": (
        float(scores.median()) if scores.notna().any() else np.nan
      ),
      "maximum_evidence_coverage": (
        float(coverage.max()) if coverage.notna().any() else np.nan
      ),
      "link_width": 0.25
      + 2.75 * (max_score if np.isfinite(max_score) else 0.0),
      "link_alpha": 0.06
      + 0.44 * (max_score if np.isfinite(max_score) else 0.0),
      "link_color": (
        "#C62828" if "NISE" in pair_class else "#111111"
      ),
      "link_status": "available",
    })
  links = pd.DataFrame(link_rows).sort_values(
    ["maximum_coverage_adjusted_rses", "pair_id"],
    ascending=[False, True],
  )

  track_definitions = pd.DataFrame(
    TRACKS,
    columns=[
      "track_id",
      "track_label",
      "source_column",
      "panel",
      "domain_family",
      "parent_domain",
      "ring_order",
    ],
  )
  track_definitions["aggregation"] = (
    "maximum_across_all_associated_pair_cancer_rows"
  )
  track_definitions["value_range"] = "0_to_1"
  track_definitions["missing_data_rule"] = (
    "missing_remains_NA_and_is_rendered_as_hollow_marker"
  )

  gene_rows = []
  for gene in genes:
    subset = rank_selected.loc[
      rank_selected["lost_gene"].eq(gene)
      | rank_selected["target_gene"].eq(gene)
    ]
    for track in track_definitions.to_dict("records"):
      column = track["source_column"]
      values = (
        pd.to_numeric(subset[column], errors="coerce")
        if column in subset
        else pd.Series(dtype=float)
      )
      observed = values.dropna()
      gene_rows.append({
        "gene": gene,
        **{
          key: track[key]
          for key in (
            "track_id",
            "track_label",
            "source_column",
            "panel",
            "domain_family",
            "parent_domain",
            "ring_order",
          )
        },
        "value": (
          float(observed.max()) if not observed.empty else np.nan
        ),
        "median_value": (
          float(observed.median()) if not observed.empty else np.nan
        ),
        "minimum_value": (
          float(observed.min()) if not observed.empty else np.nan
        ),
        "observed_pair_cancer_rows": int(len(observed)),
        "eligible_pair_cancer_rows": int(len(subset)),
        "evidence_status": (
          "observed" if not observed.empty else "missing_or_not_eligible"
        ),
        "aggregation": track["aggregation"],
      })
  ring_values = pd.DataFrame(gene_rows).merge(
    coordinates[[
      "gene",
      "gene_class",
      "chromosome",
      "genomic_position",
    ]],
    on="gene",
    how="left",
  )

  expression = read_expression_subset(paths["expression"], set(genes))
  models = normalize_model_id_column(
    pd.read_csv(paths["models"], low_memory=False),
    "Model.csv",
  )
  expression_rows = []
  for cancer in CANCERS:
    identifiers = set(cancer_model_ids(models, cancer).astype(str))
    subset = expression.loc[
      expression["ModelID"].astype(str).isin(identifiers)
    ]
    if subset.empty:
      continue
    long = subset.melt(
      id_vars="ModelID",
      var_name="gene",
      value_name="expression_log2_tpm_plus_1",
    )
    long["cancer"] = cancer
    long["source_file"] = str(paths["expression"])
    expression_rows.append(long)
  expression_long = (
    pd.concat(expression_rows, ignore_index=True)
    if expression_rows
    else pd.DataFrame(columns=[
      "ModelID",
      "gene",
      "expression_log2_tpm_plus_1",
      "cancer",
      "source_file",
    ])
  )
  expression_long["expression_log2_tpm_plus_1"] = pd.to_numeric(
    expression_long["expression_log2_tpm_plus_1"],
    errors="coerce",
  )
  expression_summary = (
    expression_long.groupby(["cancer", "gene"], as_index=False)
    .agg(
      n_models=("ModelID", "nunique"),
      observed_values=("expression_log2_tpm_plus_1", "count"),
      median_expression=("expression_log2_tpm_plus_1", "median"),
      mean_expression=("expression_log2_tpm_plus_1", "mean"),
      q25_expression=(
        "expression_log2_tpm_plus_1",
        lambda values: values.quantile(0.25),
      ),
      q75_expression=(
        "expression_log2_tpm_plus_1",
        lambda values: values.quantile(0.75),
      ),
    )
  )
  expression_summary["unit"] = "log2(TPM+1)"
  expression_summary["source_file"] = str(paths["expression"])

  output_dir = resolve(args.output_dir)
  output_dir.mkdir(parents=True, exist_ok=True)
  outputs = {
    "gene_coordinates": (
      output_dir / "genomic_circos_gene_coordinates.tsv"
    ),
    "pair_links": output_dir / "genomic_circos_pair_links.tsv",
    "ring_values": output_dir / "genomic_circos_ring_values.tsv",
    "track_definitions": (
      output_dir / "genomic_circos_track_definitions.tsv"
    ),
    "expression_summary": (
      output_dir / "genomic_circos_expression_summary.tsv"
    ),
    "expression_model_values": (
      output_dir / "genomic_circos_expression_model_values.tsv"
    ),
    "source_provenance": (
      output_dir / "genomic_circos_source_provenance.tsv"
    ),
  }
  atomic_tsv(
    coordinates.sort_values([
      "chromosome",
      "genomic_position",
      "gene",
    ]),
    outputs["gene_coordinates"],
  )
  atomic_tsv(links, outputs["pair_links"])
  atomic_tsv(
    ring_values.sort_values([
      "panel",
      "ring_order",
      "chromosome",
      "genomic_position",
    ]),
    outputs["ring_values"],
  )
  atomic_tsv(track_definitions, outputs["track_definitions"])
  atomic_tsv(
    expression_summary.sort_values(["cancer", "gene"]),
    outputs["expression_summary"],
  )
  atomic_tsv(
    expression_long.sort_values(["cancer", "gene", "ModelID"]),
    outputs["expression_model_values"],
  )

  provenance_rows = []
  for role, path in paths.items():
    row_count = None
    try:
      row_count = (
        sum(1 for _ in path.open("rb")) - 1
        if path.suffix in {".tsv", ".csv"}
        else None
      )
    except OSError:
      row_count = None
    provenance_rows.append({
      "role": role,
      "path": str(path),
      "bytes": path.stat().st_size,
      "sha256": sha256(path),
      "approximate_data_rows": row_count,
      "used_by": "scripts/build_genomic_circos_inputs.py",
    })
  atomic_tsv(
    pd.DataFrame(provenance_rows),
    outputs["source_provenance"],
  )

  status = {
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "candidate_pairs": int(len(selected)),
    "genes": int(len(coordinates)),
    "links": int(len(links)),
    "tracks": int(len(track_definitions)),
    "ring_rows": int(len(ring_values)),
    "expression_model_rows": int(len(expression_long)),
    "coordinate_assembly": ";".join(sorted(set(
      coordinates.get(
        "assembly",
        pd.Series(dtype=str),
      ).dropna().astype(str)
    ))),
    "coordinate_missing": int((
      ~coordinates["coordinate_status"].eq("available")
    ).sum()),
    "outputs": {
      key: str(value)
      for key, value in outputs.items()
    },
  }
  (
    output_dir / "genomic_circos_status.json"
  ).write_text(
    json.dumps(status, indent=2),
    encoding="utf-8",
  )
  print(json.dumps(status, indent=2))


if __name__ == "__main__":
  main()
