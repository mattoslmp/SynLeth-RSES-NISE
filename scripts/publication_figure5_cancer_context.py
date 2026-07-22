#!/usr/bin/env python3
"""Cancer-specific Figure 5 renderer with explicit missingness semantics."""
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
  """Render one row per candidate–cancer context without zero imputation."""
  audit_path = ROOT / "article_outputs" / "tables" / "qc" / "candidate_domain_evidence_audit.tsv"
  if not audit_path.exists() or audit_path.stat().st_size == 0:
    raise RuntimeError(f"Figure 5 requires the candidate-domain audit: {audit_path}")
  audit = pd.read_csv(audit_path, sep="\t", low_memory=False)
  ranked = add_display_pair_columns(ranking)
  pair_level = (
    ranked.sort_values("functional_microniche_adjusted", ascending=False)
      .drop_duplicates(["cancer", "pair_id"])
      .head(top_n)
      .copy()
  )
  pair_level["cancer_label"] = pair_level["cancer"].astype(str).map(CANCER_LABELS).fillna(pair_level["cancer"].astype(str))
  pair_level["display_context_label"] = pair_level["display_pair_label"].astype(str) + " | " + pair_level["cancer_label"].astype(str)
  keys = pair_level[["cancer", "pair_id", "display_pair_label", "display_context_label", "cancer_label"]].copy()
  cells = audit.loc[
    audit["domain_family"].astype(str).eq("Functional microniche")
    & audit["domain"].astype(str).isin(DOMAIN_ORDER)
  ].merge(
    keys,
    left_on=["cancer", "candidate_id"],
    right_on=["cancer", "pair_id"],
    how="inner",
  )
  cells["domain"] = pd.Categorical(cells["domain"], categories=DOMAIN_ORDER, ordered=True)
  pair_order = pair_level["display_context_label"].tolist()
  value_matrix = cells.pivot_table(
    index="display_context_label", columns="domain", values="component_normalized", aggfunc="first",
  ).reindex(index=pair_order, columns=DOMAIN_ORDER)
  status_matrix = cells.pivot_table(
    index="display_context_label", columns="domain", values="evidence_state", aggfunc="first",
  ).reindex(index=pair_order, columns=DOMAIN_ORDER)

  set_publication_style()
  fig, axis = plt.subplots(
    figsize=(15.0, max(9.0, 0.52 * len(value_matrix))),
    constrained_layout=True,
  )
  masked = np.ma.masked_invalid(value_matrix.to_numpy(dtype=float))
  image = axis.imshow(masked, aspect="auto", vmin=0, vmax=1)
  axis.set_xticks(
    np.arange(len(DOMAIN_ORDER)),
    [wrap_label(DOMAIN_LABELS[domain], 18) for domain in DOMAIN_ORDER],
    rotation=25,
    ha="right",
  )
  axis.set_yticks(
    np.arange(len(value_matrix)),
    [wrap_label(value, 42) for value in value_matrix.index],
  )
  axis.set_xlabel("Functional-microniche domain")
  axis.set_ylabel("Prioritized directed hypothesis and cancer context")
  colorbar = fig.colorbar(image, ax=axis, fraction=0.025, pad=0.02)
  colorbar.set_label("Observed specialization / divergence value")

  unavailable_states = {"missing", "technical_failure", "insufficient_sample", "nan"}
  for y in range(len(status_matrix)):
    for x in range(len(status_matrix.columns)):
      state = str(status_matrix.iloc[y, x])
      if state == "not_eligible":
        axis.add_patch(patches.Rectangle(
          (x - 0.5, y - 0.5), 1, 1,
          facecolor="0.75", edgecolor="0.35", hatch="///", linewidth=0.5,
        ))
      elif state in unavailable_states:
        axis.add_patch(patches.Rectangle(
          (x - 0.5, y - 0.5), 1, 1,
          facecolor="white", edgecolor="0.55", hatch="...", linewidth=0.5,
        ))
      elif pd.notna(value_matrix.iloc[y, x]):
        axis.text(
          x, y, f"{float(value_matrix.iloc[y, x]):.2f}",
          ha="center", va="center", fontsize=7.2,
        )

  axis.legend(
    handles=[
      patches.Patch(
        facecolor="white", edgecolor="0.55", hatch="...",
        label="Evidence unavailable / insufficient",
      ),
      patches.Patch(
        facecolor="0.75", edgecolor="0.35", hatch="///",
        label="Domain not eligible",
      ),
    ],
    loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=2, frameon=False,
  )

  cells["display_value"] = cells["component_normalized"]
  cells["display_status"] = cells["evidence_state"]
  cells["calculation_rule"] = (
    "Observed value shown; missing and non-eligible cells are not imputed."
  )
  cells["source_and_reason"] = (
    cells["evidence_source"].astype(str)
    + " | "
    + cells["absence_reason"].fillna("").astype(str)
  )
  return module.save_record(
    fig=fig,
    item=item,
    output_dir=output_dir,
    source_dir=source_dir,
    source_data=cells,
    inputs=[input_path, audit_path],
    strict_layout=strict,
  )
