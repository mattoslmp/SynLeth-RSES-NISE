from __future__ import annotations

import io
import zipfile
from pathlib import Path

import pandas as pd
import pytest

from scripts.download_hpa_subcellular_current import (
  HPA_URLS,
  REQUIRED_COLUMNS,
  parse_hpa_archive,
)


def build_archive(frame: pd.DataFrame, member: str = "subcellular_location.tsv") -> bytes:
  payload = io.BytesIO()
  with zipfile.ZipFile(payload, "w", compression=zipfile.ZIP_DEFLATED) as archive:
    archive.writestr(member, frame.to_csv(sep="\t", index=False))
  return payload.getvalue()


def valid_frame() -> pd.DataFrame:
  return pd.DataFrame({
    "Gene": ["ENSG00000141510"],
    "Gene name": ["TP53"],
    "Reliability": ["Supported"],
    "Enhanced": [""],
    "Supported": ["Nucleoplasm"],
    "Approved": [""],
    "Uncertain": [""],
  })


def test_current_hpa_url_precedes_legacy_url() -> None:
  assert HPA_URLS[0].endswith("/download/tsv/subcellular_location.tsv.zip")
  assert HPA_URLS[1].endswith("/download/subcellular_location.tsv.zip")


def test_parse_hpa_archive_accepts_current_schema() -> None:
  frame, member = parse_hpa_archive(build_archive(valid_frame()))
  assert member == "subcellular_location.tsv"
  assert len(frame) == 1
  assert REQUIRED_COLUMNS.issubset(frame.columns)


def test_parse_hpa_archive_rejects_missing_required_columns() -> None:
  invalid = valid_frame().drop(columns=["Supported"])
  with pytest.raises(ValueError, match="Unexpected HPA subcellular schema"):
    parse_hpa_archive(build_archive(invalid))


def test_parse_hpa_archive_rejects_non_zip() -> None:
  with pytest.raises(ValueError, match="not a valid ZIP"):
    parse_hpa_archive(b"not a zip")
