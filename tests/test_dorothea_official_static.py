from __future__ import annotations

import gzip
import io
from pathlib import Path

import pandas as pd
import pytest

from scripts.download_dorothea_official_static import (
  normalize_levels,
  parse_static_archive,
  read_valid_cache,
  validate_and_filter,
)


def example_frame() -> pd.DataFrame:
  return pd.DataFrame({
    "source": ["P04637", "P01106", "P84022"],
    "target": ["P38936", "O14746", "P05412"],
    "source_genesymbol": ["TP53", "MYC", "SMAD3"],
    "target_genesymbol": ["CDKN1A", "TERT", "JUN"],
    "dorothea_level": ["A", "D", "C"],
    "sources": ["DoRothEA", "DoRothEA", "DoRothEA"],
    "references": ["1", "2", "3"],
  })


def gzip_tsv(frame: pd.DataFrame) -> bytes:
  raw = frame.to_csv(sep="\t", index=False).encode("utf-8")
  return gzip.compress(raw)


def test_normalize_levels() -> None:
  assert normalize_levels("A,B,C") == ("A", "B", "C")
  assert normalize_levels("a, c, a") == ("A", "C")
  with pytest.raises(ValueError, match="Unsupported"):
    normalize_levels("A,Z")


def test_parse_static_archive_filters_requested_levels() -> None:
  frame = parse_static_archive(gzip_tsv(example_frame()), ("A", "B", "C"))
  assert frame["source_genesymbol"].tolist() == ["TP53", "SMAD3"]
  assert set(frame["dorothea_level"]) == {"A", "C"}


def test_validate_static_schema() -> None:
  invalid = example_frame().drop(columns=["dorothea_level"])
  with pytest.raises(ValueError, match="dorothea_level"):
    validate_and_filter(invalid, ("A", "B", "C"))


def test_invalid_gzip_is_rejected() -> None:
  with pytest.raises(ValueError, match="valid gzip"):
    parse_static_archive(b"not-gzip", ("A",))


def test_valid_cache_is_reused(tmp_path: Path) -> None:
  path = tmp_path / "dorothea.tsv"
  example_frame().to_csv(path, sep="\t", index=False)
  frame = read_valid_cache(path, ("A", "C"))
  assert len(frame) == 2


def test_html_response_is_rejected() -> None:
  with pytest.raises(ValueError, match="HTML"):
    parse_static_archive(b"<!doctype html><html></html>", ("A",))
