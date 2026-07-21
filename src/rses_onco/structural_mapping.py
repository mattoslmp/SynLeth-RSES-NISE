from __future__ import annotations

from typing import Any, Iterable

from .structural import StructuralResidue
from .utils import canonical_gene_name


def _walk(value: Any) -> Iterable[dict[str, Any]]:
  if isinstance(value, dict):
    yield value
    for child in value.values():
      yield from _walk(child)
  elif isinstance(value, list):
    for child in value:
      yield from _walk(child)


def _explicit_uniprot_number(record: dict[str, Any]) -> int | None:
  """Return a residue number only from fields explicitly declared as UniProt.

  M-CSA's curated reference residue identifiers commonly use PDB numbering. Those
  values must not be projected directly onto an AlphaFold model. This parser accepts
  only fields whose names explicitly indicate UniProt/target-sequence numbering.
  """
  keys = (
    "uniprot_residue_number",
    "uniprot_position",
    "target_uniprot_residue_number",
    "target_sequence_residue_number",
    "aligned_uniprot_residue_number",
  )
  for key in keys:
    value = record.get(key)
    if isinstance(value, dict):
      value = value.get("value") or value.get("residue_number")
    try:
      if value is not None:
        number = int(value)
        return number if number > 0 else None
    except (TypeError, ValueError):
      continue
  return None


def parse_mcsa_uniprot_residues(
  payload: Any,
  gene_symbol: str,
  accession: str,
) -> list[StructuralResidue]:
  rows: list[StructuralResidue] = []
  seen: set[tuple[int, str]] = set()
  for record in _walk(payload):
    residue_number = _explicit_uniprot_number(record)
    if residue_number is None:
      continue
    residue_name = (
      record.get("residue_name")
      or record.get("target_residue_name")
      or record.get("amino_acid")
    )
    description = " ".join(
      str(record.get(key) or "")
      for key in ("role", "roles", "function", "description", "function_location_name")
    ).strip()
    key = (residue_number, description)
    if key in seen:
      continue
    seen.add(key)
    rows.append(StructuralResidue(
      gene_symbol=canonical_gene_name(gene_symbol),
      uniprot_accession=accession,
      residue_number=residue_number,
      residue_name=str(residue_name) if residue_name else None,
      annotation_type="mcsa_catalytic",
      source="M-CSA residue with explicit UniProt numbering",
      description=description or "Catalytic residue mapped by M-CSA",
      evidence_level="mcsa_explicit_uniprot_mapping",
      mapping_status="exact_uniprot_numbering",
    ))
  return rows
