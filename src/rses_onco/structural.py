from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Iterable
import hashlib
import json
import math
import re

import pandas as pd

from .utils import canonical_gene_name


STRUCTURE_RESIDUE_COLORS = {
  "mcsa_catalytic": "red",
  "uniprot_active_site": "red",
  "uniprot_binding": "orange",
  "uniprot_metal": "violet",
  "uniprot_site": "yellow",
  "pdbe_ligand": "magenta",
  "drug_binding": "cyan",
  "curated_user": "marine",
}


@dataclass(frozen=True)
class StructuralResidue:
  gene_symbol: str
  uniprot_accession: str
  residue_number: int
  residue_name: str | None
  annotation_type: str
  source: str
  description: str
  ligand_code: str | None = None
  ligand_name: str | None = None
  pdb_id: str | None = None
  evidence_level: str = "exact_uniprot"
  mapping_status: str = "exact"


def file_sha256(path: str | Path) -> str:
  digest = hashlib.sha256()
  with Path(path).open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def _position(location: dict[str, Any], key: str) -> int | None:
  value = location.get(key)
  if isinstance(value, dict):
    value = value.get("value")
  try:
    return int(value)
  except (TypeError, ValueError):
    return None


def _feature_positions(location: dict[str, Any]) -> list[int]:
  start = _position(location, "start")
  end = _position(location, "end")
  if start is None:
    return []
  end = start if end is None else end
  if end < start or end - start > 500:
    return [start]
  return list(range(start, end + 1))


def parse_uniprot_features(
  payload: dict[str, Any],
  gene_symbol: str,
  accession: str,
) -> list[StructuralResidue]:
  """Parse exact UniProt residue features suitable for AlphaFold numbering."""
  mapping = {
    "Active site": "uniprot_active_site",
    "Binding site": "uniprot_binding",
    "Metal binding": "uniprot_metal",
    "Site": "uniprot_site",
    "Nucleotide binding": "uniprot_binding",
  }
  rows: list[StructuralResidue] = []
  for feature in payload.get("features", []):
    feature_type = str(feature.get("type") or "")
    annotation_type = mapping.get(feature_type)
    if annotation_type is None:
      continue
    description = str(feature.get("description") or feature_type)
    ligand = feature.get("ligand") or {}
    ligand_name = ligand.get("name") if isinstance(ligand, dict) else None
    ligand_code = ligand.get("id") if isinstance(ligand, dict) else None
    for residue_number in _feature_positions(feature.get("location") or {}):
      rows.append(StructuralResidue(
        gene_symbol=canonical_gene_name(gene_symbol),
        uniprot_accession=accession,
        residue_number=residue_number,
        residue_name=None,
        annotation_type=annotation_type,
        source="UniProtKB reviewed feature",
        description=description,
        ligand_code=str(ligand_code) if ligand_code else None,
        ligand_name=str(ligand_name) if ligand_name else None,
      ))
  return rows


def _walk_dicts(value: Any) -> Iterable[dict[str, Any]]:
  if isinstance(value, dict):
    yield value
    for child in value.values():
      yield from _walk_dicts(child)
  elif isinstance(value, list):
    for child in value:
      yield from _walk_dicts(child)


def _first_int(record: dict[str, Any], keys: Iterable[str]) -> int | None:
  for key in keys:
    value = record.get(key)
    if isinstance(value, dict):
      value = value.get("value") or value.get("residue_number")
    if isinstance(value, str):
      match = re.search(r"-?\d+", value)
      value = match.group(0) if match else None
    try:
      if value is not None:
        return int(value)
    except (TypeError, ValueError):
      continue
  return None


def parse_mcsa_residues(
  payload: Any,
  gene_symbol: str,
  accession: str,
) -> list[StructuralResidue]:
  """Parse exact-accession M-CSA residue API responses conservatively."""
  rows: list[StructuralResidue] = []
  seen: set[tuple[int, str]] = set()
  for record in _walk_dicts(payload):
    residue_number = _first_int(
      record,
      ("uniprot_residue_number", "residue_number", "resid", "residue_id"),
    )
    if residue_number is None or residue_number <= 0:
      continue
    residue_name = (
      record.get("residue_name")
      or record.get("code")
      or record.get("residue")
      or record.get("amino_acid")
    )
    text = " ".join(str(record.get(key) or "") for key in (
      "role", "roles", "function", "description", "function_location_name"
    )).strip()
    if not text and not any(key in record for key in ("residue", "resid", "residue_number")):
      continue
    key = (residue_number, text)
    if key in seen:
      continue
    seen.add(key)
    rows.append(StructuralResidue(
      gene_symbol=canonical_gene_name(gene_symbol),
      uniprot_accession=accession,
      residue_number=residue_number,
      residue_name=str(residue_name) if residue_name else None,
      annotation_type="mcsa_catalytic",
      source="M-CSA exact UniProt query",
      description=text or "Catalytic residue annotated in M-CSA",
      evidence_level="exact_or_curated_family",
    ))
  return rows


def parse_pdbe_binding_sites(
  payload: Any,
  gene_symbol: str,
  accession: str,
  pdb_id: str,
) -> list[StructuralResidue]:
  """Parse only binding residues explicitly returned in UniProt numbering.

  PDB-only numbering is intentionally not projected onto AlphaFold models without
  a residue-level SIFTS mapping.
  """
  rows: list[StructuralResidue] = []
  seen: set[tuple[int, str | None]] = set()
  for record in _walk_dicts(payload):
    residue_number = _first_int(
      record,
      ("uniprot_residue_number", "uniprot_residue", "unp_residue_number"),
    )
    if residue_number is None or residue_number <= 0:
      continue
    ligand_code = (
      record.get("ligand_code")
      or record.get("chem_comp_id")
      or record.get("het_code")
    )
    ligand_name = record.get("ligand_name") or record.get("chem_comp_name")
    key = (residue_number, str(ligand_code) if ligand_code else None)
    if key in seen:
      continue
    seen.add(key)
    rows.append(StructuralResidue(
      gene_symbol=canonical_gene_name(gene_symbol),
      uniprot_accession=accession,
      residue_number=residue_number,
      residue_name=str(record.get("residue_name") or "") or None,
      annotation_type="pdbe_ligand",
      source="PDBe Arpeggio binding-site annotation",
      description=str(record.get("description") or "Experimentally observed ligand-binding residue"),
      ligand_code=str(ligand_code) if ligand_code else None,
      ligand_name=str(ligand_name) if ligand_name else None,
      pdb_id=pdb_id.upper(),
      evidence_level="experimental_pdb_exact_uniprot_mapping",
    ))
  return rows


def residues_to_frame(rows: Iterable[StructuralResidue]) -> pd.DataFrame:
  frame = pd.DataFrame([asdict(row) for row in rows])
  if frame.empty:
    return pd.DataFrame(columns=[field.name for field in StructuralResidue.__dataclass_fields__.values()])
  return (
    frame.drop_duplicates([
      "uniprot_accession", "residue_number", "annotation_type", "source",
      "ligand_code", "pdb_id",
    ])
    .sort_values(["gene_symbol", "residue_number", "annotation_type"])
  )


def mean_plddt_from_pdb(path: str | Path) -> float | None:
  values: list[float] = []
  with Path(path).open("r", encoding="utf-8", errors="ignore") as handle:
    for line in handle:
      if not line.startswith("ATOM") or len(line) < 66:
        continue
      try:
        values.append(float(line[60:66]))
      except ValueError:
        continue
  if not values:
    return None
  return float(sum(values) / len(values))


def plddt_class(value: float | None) -> str:
  if value is None or not math.isfinite(value):
    return "unavailable"
  if value >= 90:
    return "very_high"
  if value >= 70:
    return "confident"
  if value >= 50:
    return "low"
  return "very_low"


def residue_selection(numbers: Iterable[int]) -> str:
  unique = sorted({int(number) for number in numbers if int(number) > 0})
  return "+".join(str(number) for number in unique)


def build_pymol_script(
  structure_path: str | Path,
  output_png: str | Path,
  annotations: pd.DataFrame,
  *,
  view: str = "whole",
  width: int = 2800,
  height: int = 2200,
  object_name: str = "model",
) -> str:
  """Build a deterministic publication-oriented PyMOL script."""
  structure_path = Path(structure_path).resolve()
  output_png = Path(output_png).resolve()
  lines = [
    "reinitialize",
    f"load {structure_path.as_posix()}, {object_name}",
    "hide everything, all",
    f"show cartoon, {object_name}",
    "bg_color white",
    "set ray_opaque_background, on",
    "set antialias, 2",
    "set ray_trace_mode, 1",
    "set ray_shadows, off",
    "set ambient, 0.35",
    "set direct, 0.65",
    "set specular, 0.20",
    "set cartoon_smooth_loops, on",
    "set cartoon_fancy_helices, on",
    "set stick_radius, 0.23",
    f"color gray80, {object_name}",
    f"select af_very_high, {object_name} and b >= 90",
    f"select af_confident, {object_name} and b >= 70 and b < 90",
    f"select af_low, {object_name} and b >= 50 and b < 70",
    f"select af_very_low, {object_name} and b < 50",
    "color marine, af_very_high",
    "color cyan, af_confident",
    "color yelloworange, af_low",
    "color salmon, af_very_low",
    "set cartoon_transparency, 0.05, af_low or af_very_low",
  ]
  important_numbers: list[int] = []
  if not annotations.empty:
    for annotation_type, group in annotations.groupby("annotation_type", dropna=False):
      selection = residue_selection(group["residue_number"].dropna().astype(int))
      if not selection:
        continue
      important_numbers.extend(group["residue_number"].dropna().astype(int).tolist())
      safe_type = re.sub(r"[^A-Za-z0-9_]+", "_", str(annotation_type))
      color = STRUCTURE_RESIDUE_COLORS.get(str(annotation_type), "hotpink")
      lines.extend([
        f"select site_{safe_type}, {object_name} and resi {selection}",
        f"show sticks, site_{safe_type}",
        f"show spheres, site_{safe_type} and name CA",
        f"set sphere_scale, 0.38, site_{safe_type} and name CA",
        f"color {color}, site_{safe_type}",
      ])
  lines.extend([
    f"orient {object_name}",
    f"zoom {object_name}, 6",
  ])
  if view == "site" and important_numbers:
    selection = residue_selection(important_numbers)
    lines.extend([
      f"select all_sites, {object_name} and resi {selection}",
      "orient all_sites",
      "zoom all_sites, 8",
    ])
  lines.extend([
    "set orthoscopic, on",
    f"png {output_png.as_posix()}, width={width}, height={height}, dpi=600, ray=1",
    "quit",
  ])
  return "\n".join(lines) + "\n"


def write_json(path: str | Path, payload: Any) -> Path:
  path = Path(path)
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
  return path
