from __future__ import annotations

import math
import re
from typing import Iterable

import numpy as np


def canonical_gene_name(value: str) -> str:
  """Return a DepMap/TCGA-compatible HGNC-like gene symbol.

  DepMap columns are commonly encoded as ``GENE (EntrezID)``. This helper
  strips the parenthesized identifier while preserving ordinary symbols.
  """
  text = str(value).strip()
  return re.sub(r"\s*\([^)]*\)\s*$", "", text).upper()


def bh_adjust(p_values: Iterable[float]) -> np.ndarray:
  """Benjamini-Hochberg adjusted P values with NaN preservation."""
  values = np.asarray(list(p_values), dtype=float)
  out = np.full(values.shape, np.nan, dtype=float)
  valid = np.isfinite(values)
  if not valid.any():
    return out
  pv = values[valid]
  order = np.argsort(pv)
  ranked = pv[order]
  n = len(ranked)
  adjusted = ranked * n / np.arange(1, n + 1)
  adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
  adjusted = np.clip(adjusted, 0.0, 1.0)
  restored = np.empty(n, dtype=float)
  restored[order] = adjusted
  out[valid] = restored
  return out


def clamp01(value: float | None) -> float | None:
  if value is None or not math.isfinite(float(value)):
    return None
  return float(max(0.0, min(1.0, float(value))))
