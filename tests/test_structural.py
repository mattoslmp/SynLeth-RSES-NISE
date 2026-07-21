from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from rses_onco.structural import (
  build_pymol_script,
  parse_mcsa_residues,
  parse_pdbe_binding_sites,
  parse_uniprot_features,
  residues_to_frame,
)

ROOT = Path(__file__).resolve().parents[1]


def test_uniprot_exact_residue_features() -> None:
  payload = {
    "features": [
      {
        "type": "Active site",
        "description": "Proton donor",
        "location": {"start": {"value": 42}, "end": {"value": 42}},
      },
      {
        "type": "Binding site",
        "description": "Substrate",
        "ligand": {"name": "NAD", "id": "CHEBI:57540"},
        "location": {"start": {"value": 105}, "end": {"value": 105}},
      },
    ]
  }
  rows = parse_uniprot_features(payload, "TEST", "P00001")
  assert len(rows) == 2
  frame = residues_to_frame(rows)
  assert set(frame["residue_number"]) == {42, 105}
  assert set(frame["annotation_type"]) == {
    "uniprot_active_site", "uniprot_binding"
  }


def test_mcsa_parser_requires_residue_number_fields() -> None:
  payload = {
    "results": [
      {
        "residue_number": 17,
        "residue_name": "HIS",
        "role": "general acid/base",
      },
      {"mcsa_id": 999, "description": "not a residue record"},
    ]
  }
  rows = parse_mcsa_residues(payload, "TEST", "P00001")
  assert len(rows) == 1
  assert rows[0].residue_number == 17
  assert rows[0].annotation_type == "mcsa_catalytic"


def test_pdbe_parser_uses_only_uniprot_numbering() -> None:
  payload = {
    "sites": [
      {
        "uniprot_residue_number": 80,
        "residue_name": "ASP",
        "ligand_code": "ATP",
      },
      {
        "author_residue_number": 90,
        "residue_name": "LYS",
        "ligand_code": "DRG",
      },
    ]
  }
  rows = parse_pdbe_binding_sites(payload, "TEST", "P00001", "1abc")
  assert len(rows) == 1
  assert rows[0].residue_number == 80
  assert rows[0].pdb_id == "1ABC"


def test_pymol_script_highlights_known_residues(tmp_path: Path) -> None:
  structure = tmp_path / "model.pdb"
  structure.write_text("ATOM      1  CA  ALA A  42      0.000   0.000   0.000  1.00 95.00           C\n")
  annotations = pd.DataFrame({
    "residue_number": [42, 105],
    "annotation_type": ["mcsa_catalytic", "uniprot_binding"],
  })
  script = build_pymol_script(
    structure,
    tmp_path / "render.png",
    annotations,
    view="site",
  )
  assert "resi 42" in script
  assert "resi 105" in script
  assert "width=2800" in script
  assert "height=2200" in script
  assert "dpi=600" in script
  assert "ray=1" in script


def test_article_registry_includes_complete_structural_atlas() -> None:
  config = yaml.safe_load(
    (ROOT / "config/article_assets.yaml").read_text(encoding="utf-8")
  )
  main_ids = {record["id"] for record in config["main_figures"]}
  supplementary_ids = {record["id"] for record in config["supplementary_figures"]}
  assert main_ids == {f"Figure_{index}" for index in range(1, 9)}
  assert supplementary_ids == {f"Figure_S{index}" for index in range(1, 33)}
  assert len(config["main_tables"]) == 4
  assert len(config["supplementary_tables"]) == 18
