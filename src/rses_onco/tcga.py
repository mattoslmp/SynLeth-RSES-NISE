from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .utils import canonical_gene_name


def read_gistic_matrix(path: str | Path) -> pd.DataFrame:
  """Read a gene-by-sample discrete GISTIC matrix.

  The first column may be named Hugo_Symbol, Gene Symbol, or gene. Values are
  expected on the conventional -2, -1, 0, 1, 2 scale.
  """
  frame = pd.read_csv(path, sep=None, engine="python")
  gene_col = next((c for c in frame.columns if str(c).casefold() in {
    "hugo_symbol", "gene symbol", "gene", "gene_symbol"
  }), frame.columns[0])
  frame = frame.rename(columns={gene_col: "gene_symbol"})
  frame["gene_symbol"] = frame["gene_symbol"].map(canonical_gene_name)
  return frame.set_index("gene_symbol")


def homozygous_deletion_frequency(matrix: pd.DataFrame, gene: str) -> tuple[int, int, float] | None:
  gene = canonical_gene_name(gene)
  if gene not in matrix.index:
    return None
  values = pd.to_numeric(matrix.loc[gene], errors="coerce").dropna()
  if values.empty:
    return None
  n_deleted = int((values <= -2).sum())
  n_total = int(len(values))
  return n_deleted, n_total, n_deleted / n_total


def event_component(frequency: float, saturation_frequency: float = 0.20) -> float:
  return float(np.clip(frequency / saturation_frequency, 0.0, 1.0))
