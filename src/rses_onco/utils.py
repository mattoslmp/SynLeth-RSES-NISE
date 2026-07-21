from __future__ import annotations

import math
import re
from typing import Iterable

import numpy as np


GENE_LIST_DELIMITER = re.compile(r"\s*(?:/|,|;|\||\+)\s*")
ATOMIC_GENE_SYMBOL = re.compile(r"^[A-Z0-9][A-Z0-9-]{0,31}$")


def canonical_gene_name(value: object) -> str:
  """Return a DepMap/TCGA-compatible HGNC-like gene symbol.

  DepMap columns are commonly encoded as ``GENE (EntrezID)``. This helper
  strips the parenthesized identifier while preserving ordinary symbols.
  Missing values remain empty and are never converted to the strings ``NAN``
  or ``NONE``.
  """
  if value is None:
    return ""
  if isinstance(value, float) and not math.isfinite(value):
    return ""
  text = str(value).strip()
  if text.casefold() in {"", "nan", "none", "<na>"}:
    return ""
  return re.sub(r"\s*\([^)]*\)\s*$", "", text).upper()


def atomic_gene_symbols(value: object) -> list[str]:
  """Parse a field explicitly declared as one gene or a delimiter-separated gene list.

  The function is intentionally strict. Every delimited component must look like an
  atomic HGNC-style symbol; prose such as ``BRCA1/BRCA2 or HRD`` therefore returns an
  empty list rather than silently inventing a partial gene set. Composite biological
  states remain represented by their original feature field and should not be passed
  here unless the field semantically denotes gene targets or members.
  """
  if value is None:
    return []
  if isinstance(value, float) and not math.isfinite(value):
    return []
  text = str(value).strip()
  if text.casefold() in {"", "nan", "none", "<na>"}:
    return []

  raw_parts = GENE_LIST_DELIMITER.split(text)
  symbols: list[str] = []
  for raw_part in raw_parts:
    symbol = canonical_gene_name(raw_part)
    if not symbol or not ATOMIC_GENE_SYMBOL.fullmatch(symbol):
      return []
    if symbol not in symbols:
      symbols.append(symbol)
  return symbols


def is_atomic_gene_symbol(value: object) -> bool:
  """Return True only when ``value`` encodes exactly one atomic gene symbol."""
  symbols = atomic_gene_symbols(value)
  return len(symbols) == 1 and canonical_gene_name(value) == symbols[0]


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