from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str):
  path = ROOT / "scripts" / name
  specification = importlib.util.spec_from_file_location(path.stem, path)
  assert specification is not None and specification.loader is not None
  module = importlib.util.module_from_spec(specification)
  sys.modules[specification.name] = module
  specification.loader.exec_module(module)
  return module


def test_primary_tumor_is_selected_even_when_normal_is_first() -> None:
  module = load_script("aggregate_gdc_gene_cna.py")
  hit = {
    "file_id": "file-1",
    "file_name": "gene_level.tsv",
    "cases": [{
      "submitter_id": "TCGA-AA-0001",
      "samples": [
        {
          "submitter_id": "TCGA-AA-0001-11A",
          "sample_type": "Solid Tissue Normal",
        },
        {
          "submitter_id": "TCGA-AA-0001-01A",
          "sample_type": "Primary Tumor",
        },
      ],
    }],
  }
  selection = module.select_primary_tumor_sample(hit)
  assert selection.sample_id == "TCGA-AA-0001-01A"
  assert selection.sample_type == "Primary Tumor"
  assert selection.status == "primary_tumor_exact"


def test_file_name_resolves_multiple_primary_tumor_samples() -> None:
  module = load_script("aggregate_gdc_gene_cna.py")
  hit = {
    "file_id": "file-2",
    "file_name": "TCGA-AA-0001-01B.gene_level.tsv",
    "cases": [{
      "submitter_id": "TCGA-AA-0001",
      "samples": [
        {
          "submitter_id": "TCGA-AA-0001-01A",
          "sample_type": "Primary Tumor",
        },
        {
          "submitter_id": "TCGA-AA-0001-01B",
          "sample_type": "Primary Tumor",
        },
      ],
    }],
  }
  selection = module.select_primary_tumor_sample(hit)
  assert selection.sample_id == "TCGA-AA-0001-01B"
  assert selection.status == "primary_tumor_filename_match"


def test_missing_primary_tumor_record_is_rejected() -> None:
  module = load_script("aggregate_gdc_gene_cna.py")
  hit = {
    "file_id": "file-3",
    "file_name": "normal.tsv",
    "cases": [{
      "submitter_id": "TCGA-AA-0001",
      "samples": [{
        "submitter_id": "TCGA-AA-0001-11A",
        "sample_type": "Solid Tissue Normal",
      }],
    }],
  }
  with pytest.raises(ValueError, match="No Primary Tumor sample record"):
    module.select_primary_tumor_sample(hit)


def test_tcga_sample_type_validation_detects_normal_and_unparsed_labels() -> None:
  module = load_script("validate_gdc_matrices.py")
  columns = [
    "TCGA-AA-0001-01A",
    "TCGA-AA-0002-01B__12345678",
    "TCGA-AA-0003-10A",
    "TCGA-AA-0004-11A",
    "unparsed-column",
  ]
  codes, non_primary, unparsed = module.summarize_sample_types(
    module.pd.Index(columns)
  )
  assert codes["01"] == 2
  assert codes["10"] == 1
  assert codes["11"] == 1
  assert non_primary == [
    "TCGA-AA-0003-10A",
    "TCGA-AA-0004-11A",
  ]
  assert unparsed == ["unparsed-column"]
