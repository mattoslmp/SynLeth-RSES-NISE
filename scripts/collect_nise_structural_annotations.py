#!/usr/bin/env python3
"""Collect known functional residues for human NISE AlphaFold models.

Evidence hierarchy:
1. exact UniProtKB reviewed residue features;
2. M-CSA residues returned by an exact UniProt-filtered API query and carrying
   explicit UniProt/target-sequence numbering;
3. PDBe/Arpeggio ligand-binding residues only when the API explicitly returns
   UniProt residue numbering;
4. optional user-curated exact UniProt residue tables.

The script never transfers PDB residue numbers or homology-derived catalytic
residues onto AlphaFold models without an explicit mapping. AlphaFold models do
not contain ligands; known ligand/drug-contact residues are highlighted as mapped
annotations, not as predicted binding poses.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
from typing import Any

import pandas as pd
import requests

from rses_onco.structural import (
  StructuralResidue,
  parse_pdbe_binding_sites,
  parse_uniprot_features,
  residues_to_frame,
  write_json,
)
from rses_onco.structural_mapping import parse_mcsa_uniprot_residues

ROOT = Path(__file__).resolve().parents[1]
UNIPROT = "https://rest.uniprot.org/uniprotkb"
MCSA = "https://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/"
PDBE_MOLECULES = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules"
PDBE_BINDING = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites"


def resolve_path(value: str | None) -> Path | None:
  if value is None:
    return None
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def request_json(
  session: requests.Session,
  url: str,
  *,
  params: dict[str, Any] | None = None,
  retries: int = 3,
  timeout: int = 180,
) -> Any:
  error: Exception | None = None
  for attempt in range(1, retries + 1):
    try:
      response = session.get(url, params=params, timeout=timeout)
      response.raise_for_status()
      return response.json()
    except Exception as exc:
      error = exc
      if attempt == retries:
        break
      time.sleep(min(15, 2 ** attempt))
  raise RuntimeError(f"Request failed after {retries} attempts: {url}: {error}")


def uniprot_pdb_ids(payload: dict[str, Any]) -> list[str]:
  ids = []
  for reference in payload.get("uniProtKBCrossReferences", []):
    if str(reference.get("database") or "").upper() != "PDB":
      continue
    value = str(reference.get("id") or "").strip().lower()
    if value:
      ids.append(value)
  return sorted(set(ids))


def entity_ids_for_accession(payload: Any, pdb_id: str, accession: str) -> list[int]:
  entries = payload.get(pdb_id, payload) if isinstance(payload, dict) else payload
  if isinstance(entries, dict):
    entries = [entries]
  result: list[int] = []
  if not isinstance(entries, list):
    return result
  for record in entries:
    if not isinstance(record, dict):
      continue
    accessions: set[str] = set()
    for key in ("uniprot_accession", "uniprot_id", "accession"):
      value = record.get(key)
      if isinstance(value, str):
        accessions.add(value.upper())
      elif isinstance(value, list):
        accessions.update(str(item).upper() for item in value)
    for mapping in record.get("uniprot_mapping", []) or []:
      if isinstance(mapping, dict):
        value = mapping.get("accession") or mapping.get("uniprot_accession")
        if value:
          accessions.add(str(value).upper())
    if accession.upper() not in accessions:
      continue
    try:
      result.append(int(record.get("entity_id")))
    except (TypeError, ValueError):
      continue
  return sorted(set(result))


def load_curated(path: Path | None) -> list[StructuralResidue]:
  if path is None or not path.exists():
    return []
  frame = pd.read_csv(path, sep="\t")
  required = {"gene_symbol", "uniprot_accession", "residue_number", "description"}
  missing = sorted(required - set(frame.columns))
  if missing:
    raise ValueError(f"Curated residue table lacks columns: {missing}")
  rows = []
  for record in frame.to_dict("records"):
    rows.append(StructuralResidue(
      gene_symbol=str(record["gene_symbol"]),
      uniprot_accession=str(record["uniprot_accession"]),
      residue_number=int(record["residue_number"]),
      residue_name=str(record.get("residue_name") or "") or None,
      annotation_type=str(record.get("annotation_type") or "curated_user"),
      source=str(record.get("source") or "User-curated exact UniProt mapping"),
      description=str(record["description"]),
      ligand_code=str(record.get("ligand_code") or "") or None,
      ligand_name=str(record.get("ligand_name") or "") or None,
      pdb_id=str(record.get("pdb_id") or "") or None,
      evidence_level=str(record.get("evidence_level") or "user_curated_exact"),
      mapping_status=str(record.get("mapping_status") or "exact"),
    ))
  return rows


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--proteins",
    default="data/curated/human_nise_bonafide_2017.tsv",
  )
  parser.add_argument(
    "--output",
    default="data/processed/structures/nise_structural_residue_annotations.tsv",
  )
  parser.add_argument(
    "--coverage-output",
    default="data/processed/structures/nise_structural_annotation_coverage.tsv",
  )
  parser.add_argument(
    "--cache-dir",
    default="data/raw/structures/annotation_cache",
  )
  parser.add_argument("--curated-residues", default=None)
  parser.add_argument("--skip-pdbe", action="store_true")
  parser.add_argument("--max-pdb-per-protein", type=int, default=12)
  parser.add_argument("--refresh", action="store_true")
  parser.add_argument("--sleep", type=float, default=0.10)
  args = parser.parse_args()

  proteins = pd.read_csv(resolve_path(args.proteins), sep="\t")
  proteins = proteins.drop_duplicates("uniprot_accession").sort_values(
    ["group_id", "structural_cluster", "gene_symbol"]
  )
  cache_dir = resolve_path(args.cache_dir)
  assert cache_dir is not None
  cache_dir.mkdir(parents=True, exist_ok=True)
  session = requests.Session()
  session.headers.update({
    "User-Agent": "RSES-Onco/0.10 structural annotation client",
    "Accept": "application/json",
  })

  rows: list[StructuralResidue] = load_curated(resolve_path(args.curated_residues))
  status_rows: list[dict[str, Any]] = []
  for index, protein in enumerate(proteins.to_dict("records"), start=1):
    gene = str(protein["gene_symbol"])
    accession = str(protein["uniprot_accession"])
    uniprot_cache = cache_dir / "uniprot" / f"{accession}.json"
    mcsa_cache = cache_dir / "mcsa" / f"{accession}.json"
    try:
      if uniprot_cache.exists() and not args.refresh:
        uniprot_payload = json.loads(uniprot_cache.read_text(encoding="utf-8"))
      else:
        uniprot_payload = request_json(session, f"{UNIPROT}/{accession}.json")
        write_json(uniprot_cache, uniprot_payload)
      uni_rows = parse_uniprot_features(uniprot_payload, gene, accession)
      rows.extend(uni_rows)
      status_rows.append({
        "gene_symbol": gene, "uniprot_accession": accession,
        "source": "uniprot", "status": "ok", "residue_rows": len(uni_rows),
      })
    except Exception as exc:
      uniprot_payload = {}
      status_rows.append({
        "gene_symbol": gene, "uniprot_accession": accession,
        "source": "uniprot", "status": "failed", "message": str(exc),
      })

    try:
      if mcsa_cache.exists() and not args.refresh:
        mcsa_payload = json.loads(mcsa_cache.read_text(encoding="utf-8"))
      else:
        mcsa_payload = request_json(
          session,
          MCSA,
          params={
            "format": "json",
            "entries.proteins.sequences.uniprot_ids": accession,
          },
        )
        write_json(mcsa_cache, mcsa_payload)
      mcsa_rows = parse_mcsa_uniprot_residues(mcsa_payload, gene, accession)
      rows.extend(mcsa_rows)
      status_rows.append({
        "gene_symbol": gene, "uniprot_accession": accession,
        "source": "mcsa", "status": "ok", "residue_rows": len(mcsa_rows),
        "mapping_policy": "explicit_uniprot_numbering_only",
      })
    except Exception as exc:
      status_rows.append({
        "gene_symbol": gene, "uniprot_accession": accession,
        "source": "mcsa", "status": "failed", "message": str(exc),
      })

    if not args.skip_pdbe and uniprot_payload:
      pdbe_count = 0
      for pdb_id in uniprot_pdb_ids(uniprot_payload)[:args.max_pdb_per_protein]:
        try:
          molecules = request_json(session, f"{PDBE_MOLECULES}/{pdb_id}")
          entity_ids = entity_ids_for_accession(molecules, pdb_id, accession)
          for entity_id in entity_ids:
            binding = request_json(session, f"{PDBE_BINDING}/{pdb_id}/{entity_id}")
            binding_rows = parse_pdbe_binding_sites(
              binding, gene, accession, pdb_id
            )
            rows.extend(binding_rows)
            pdbe_count += len(binding_rows)
        except Exception as exc:
          status_rows.append({
            "gene_symbol": gene, "uniprot_accession": accession,
            "source": "pdbe", "pdb_id": pdb_id,
            "status": "failed", "message": str(exc),
          })
      status_rows.append({
        "gene_symbol": gene, "uniprot_accession": accession,
        "source": "pdbe", "status": "ok", "residue_rows": pdbe_count,
        "mapping_policy": "explicit_uniprot_numbering_only",
      })

    print(f"[Structural annotations {index}/{len(proteins)}] {gene} {accession}", flush=True)
    if args.sleep:
      time.sleep(args.sleep)

  frame = residues_to_frame(rows)
  output = resolve_path(args.output)
  coverage_output = resolve_path(args.coverage_output)
  assert output is not None and coverage_output is not None
  output.parent.mkdir(parents=True, exist_ok=True)
  frame.to_csv(output, sep="\t", index=False)
  status = pd.DataFrame(status_rows)
  status.to_csv(output.with_name("nise_structural_annotation_source_status.tsv"), sep="\t", index=False)

  if frame.empty:
    coverage = proteins.copy()
    coverage["annotated_residues"] = 0
  else:
    counts = (
      frame.groupby(["gene_symbol", "uniprot_accession"])
      .agg(
        annotated_residues=("residue_number", "nunique"),
        annotation_rows=("residue_number", "size"),
        annotation_sources=("source", lambda values: ";".join(sorted(set(values)))),
        catalytic_residues=("annotation_type", lambda values: int(pd.Series(values).isin(["mcsa_catalytic", "uniprot_active_site"]).sum())),
        ligand_binding_residues=("annotation_type", lambda values: int(pd.Series(values).isin(["uniprot_binding", "pdbe_ligand", "drug_binding"]).sum())),
      )
      .reset_index()
    )
    coverage = proteins.merge(counts, on=["gene_symbol", "uniprot_accession"], how="left")
    for column in ("annotated_residues", "annotation_rows", "catalytic_residues", "ligand_binding_residues"):
      coverage[column] = coverage[column].fillna(0).astype(int)
  coverage.to_csv(coverage_output, sep="\t", index=False)
  summary = {
    "proteins": int(proteins["uniprot_accession"].nunique()),
    "proteins_with_residue_annotations": int((coverage["annotated_residues"] > 0).sum()),
    "unique_annotated_residues": int(frame[["uniprot_accession", "residue_number"]].drop_duplicates().shape[0]) if not frame.empty else 0,
    "annotation_rows": int(len(frame)),
    "policy": "Only exact UniProt-numbered residues are projected onto AlphaFold models.",
  }
  write_json(output.with_suffix(".summary.json"), summary)
  print(json.dumps(summary, indent=2, sort_keys=True))
  print(f"Wrote {output}")
  print(f"Wrote {coverage_output}")


if __name__ == "__main__":
  main()
