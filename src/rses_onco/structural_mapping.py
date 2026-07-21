from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable
import re
import xml.etree.ElementTree as ET

from .structural import StructuralResidue
from .utils import canonical_gene_name


@dataclass(frozen=True)
class MCSAReferenceResidue:
  pdb_id: str
  chain_id: str
  pdb_residue_number: str
  residue_name: str | None
  description: str


def _walk(value: Any) -> Iterable[dict[str, Any]]:
  if isinstance(value, dict):
    yield value
    for child in value.values():
      yield from _walk(child)
  elif isinstance(value, list):
    for child in value:
      yield from _walk(child)


def _explicit_uniprot_number(record: dict[str, Any]) -> int | None:
  """Return a residue number only from fields explicitly declared as UniProt."""
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


def _first_text(record: dict[str, Any], keys: Iterable[str]) -> str | None:
  for key in keys:
    value = record.get(key)
    if isinstance(value, dict):
      value = value.get("id") or value.get("code") or value.get("value")
    if value is not None and str(value).strip():
      return str(value).strip()
  return None


def _pdb_id(record: dict[str, Any]) -> str | None:
  value = _first_text(record, ("pdb_id", "pdb_code", "pdb", "structure_id"))
  if value and re.fullmatch(r"[A-Za-z0-9]{4}", value):
    return value.lower()
  return None


def _description(record: dict[str, Any]) -> str:
  return " ".join(
    str(record.get(key) or "")
    for key in (
      "role", "roles", "function", "description", "function_location_name",
      "role_description", "chemical_function",
    )
  ).strip()


def parse_mcsa_uniprot_residues(
  payload: Any,
  gene_symbol: str,
  accession: str,
) -> list[StructuralResidue]:
  """Read M-CSA rows already carrying explicit UniProt-compatible numbering."""
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
    description = _description(record)
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


def extract_mcsa_reference_residues(payload: Any) -> list[MCSAReferenceResidue]:
  """Extract PDB-numbered curated M-CSA residues for subsequent SIFTS mapping.

  A record is retained only when PDB ID, chain and residue number are all explicit.
  No sequence-homology residue is accepted here.
  """
  rows: list[MCSAReferenceResidue] = []
  seen: set[tuple[str, str, str, str]] = set()
  for record in _walk(payload):
    pdb_id = _pdb_id(record)
    chain_id = _first_text(
      record,
      ("chain_id", "chain", "pdb_chain", "auth_asym_id", "struct_asym_id"),
    )
    residue_number = _first_text(
      record,
      ("pdb_residue_number", "author_residue_number", "residue_number", "resid"),
    )
    if not pdb_id or not chain_id or not residue_number:
      continue
    residue_name = _first_text(
      record,
      ("residue_name", "amino_acid", "residue", "code"),
    )
    description = _description(record) or "Catalytic residue curated by M-CSA"
    key = (pdb_id, chain_id, residue_number, description)
    if key in seen:
      continue
    seen.add(key)
    rows.append(MCSAReferenceResidue(
      pdb_id=pdb_id,
      chain_id=chain_id,
      pdb_residue_number=residue_number,
      residue_name=residue_name,
      description=description,
    ))
  return rows


def _local_name(tag: str) -> str:
  return tag.rsplit("}", 1)[-1]


def sifts_uniprot_mapping(
  xml_content: bytes,
  *,
  accession: str,
  chain_id: str,
  pdb_residue_number: str,
) -> int | None:
  """Map one exact PDB chain/residue to UniProt numbering using SIFTS XML."""
  root = ET.fromstring(xml_content)
  normalized_pdb_number = str(pdb_residue_number).strip()
  for residue in root.iter():
    if _local_name(residue.tag) != "residue":
      continue
    pdb_match = False
    uniprot_number: int | None = None
    for cross_reference in residue:
      if _local_name(cross_reference.tag) != "crossRefDb":
        continue
      attributes = cross_reference.attrib
      source = str(attributes.get("dbSource") or "").upper()
      if source == "PDB":
        mapped_chain = str(
          attributes.get("dbChainId")
          or attributes.get("dbChainID")
          or attributes.get("dbAccessionId")
          or ""
        )
        mapped_number = str(attributes.get("dbResNum") or "").strip()
        if mapped_chain == str(chain_id) and mapped_number == normalized_pdb_number:
          pdb_match = True
      elif source == "UNIPROT":
        mapped_accession = str(attributes.get("dbAccessionId") or "")
        if mapped_accession.upper() != accession.upper():
          continue
        try:
          uniprot_number = int(str(attributes.get("dbResNum") or ""))
        except ValueError:
          uniprot_number = None
    if pdb_match and uniprot_number is not None and uniprot_number > 0:
      return uniprot_number
  return None


def map_mcsa_reference_residues_with_sifts(
  reference_residues: Iterable[MCSAReferenceResidue],
  sifts_by_pdb: dict[str, bytes],
  gene_symbol: str,
  accession: str,
) -> list[StructuralResidue]:
  """Convert curated PDB-numbered M-CSA residues to exact UniProt positions."""
  rows: list[StructuralResidue] = []
  seen: set[tuple[int, str]] = set()
  for reference in reference_residues:
    xml_content = sifts_by_pdb.get(reference.pdb_id.lower())
    if not xml_content:
      continue
    mapped = sifts_uniprot_mapping(
      xml_content,
      accession=accession,
      chain_id=reference.chain_id,
      pdb_residue_number=reference.pdb_residue_number,
    )
    if mapped is None:
      continue
    key = (mapped, reference.description)
    if key in seen:
      continue
    seen.add(key)
    rows.append(StructuralResidue(
      gene_symbol=canonical_gene_name(gene_symbol),
      uniprot_accession=accession,
      residue_number=mapped,
      residue_name=reference.residue_name,
      annotation_type="mcsa_catalytic",
      source="M-CSA curated PDB residue mapped through SIFTS",
      description=reference.description,
      pdb_id=reference.pdb_id.upper(),
      evidence_level="mcsa_curated_sifts_exact",
      mapping_status="pdb_to_uniprot_sifts_exact",
    ))
  return rows
