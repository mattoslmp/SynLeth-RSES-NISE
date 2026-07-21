#!/usr/bin/env python3
"""Acquire pharmacology evidence for prioritized RSES-Onco targets.

Live, cached sources:
- Open Targets GraphQL: tractability, known drugs, indication context;
- ChEMBL REST: target mechanisms and bioactivity;
- DGIdb GraphQL: aggregated drug-gene interactions;
- MyChem.info REST: compound identifiers and annotation enrichment;
- Pharos/TCRD GraphQL: target development level and drug ligands;
- CIViC: verified gene-level precision-oncology record links.

GDSC, CTRP and PRISM are release files rather than stable per-target APIs and are
handled by ``standardize_drug_sensitivity.py`` and
``analyze_drug_response_selectivity.py``. Network/API failures are recorded in a
source-status table; they are never converted to zero evidence.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import time
from typing import Any, Iterable

import numpy as np
import pandas as pd
import requests

from rses_onco.pharmacology import (
  normalize_dgidb_interaction_score,
  normalize_open_targets_tractability,
  normalize_pharos_tdl,
)
from rses_onco.utils import canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
MYGENE = "https://mygene.info/v3/query"
OPEN_TARGETS = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
DGIDB = "https://dgidb.org/api/graphql"
MYCHEM = "https://mychem.info/v1"
PHAROS = "https://pharos-api.ncats.io/graphql"
CIVIC_LINK = "https://civicdb.org/links/entrez_name"


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def safe_name(value: str) -> str:
  return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def request_json(
  session: requests.Session,
  method: str,
  url: str,
  *,
  retries: int = 3,
  timeout: int = 180,
  **kwargs: Any,
) -> Any:
  error: Exception | None = None
  for attempt in range(1, retries + 1):
    try:
      response = session.request(method, url, timeout=timeout, **kwargs)
      response.raise_for_status()
      return response.json()
    except Exception as exc:
      error = exc
      if attempt == retries:
        break
      time.sleep(min(15.0, 2.0 ** attempt))
  raise RuntimeError(f"Request failed after {retries} attempts: {url}: {error}")


def cached_json(
  cache_path: Path,
  producer: Any,
  refresh: bool,
) -> Any:
  if cache_path.exists() and not refresh:
    return json.loads(cache_path.read_text(encoding="utf-8"))
  payload = producer()
  cache_path.parent.mkdir(parents=True, exist_ok=True)
  cache_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
  return payload


def select_targets(ranking: pd.DataFrame, minimum_score: float, maximum_targets: int) -> pd.DataFrame:
  target_column = "analysis_target_gene" if "analysis_target_gene" in ranking else "target_gene"
  selected = ranking.copy()
  selected[target_column] = selected[target_column].map(canonical_gene_name)
  selected = selected.loc[selected[target_column].ne("")]
  selected["coverage_adjusted_rses"] = pd.to_numeric(
    selected["coverage_adjusted_rses"], errors="coerce"
  )
  summary = (
    selected.groupby(target_column, as_index=False)
      .agg(
        maximum_vulnerability=("coverage_adjusted_rses", "max"),
        median_vulnerability=("coverage_adjusted_rses", "median"),
        candidate_directions=("pair_id", "nunique"),
        cancers=("cancer", lambda values: ";".join(sorted(set(map(str, values))))),
      )
      .rename(columns={target_column: "target_gene"})
      .sort_values("maximum_vulnerability", ascending=False)
  )
  summary = summary.loc[summary["maximum_vulnerability"] >= minimum_score]
  if maximum_targets > 0:
    summary = summary.head(maximum_targets)
  return summary.reset_index(drop=True)


def mygene_identifiers(session: requests.Session, gene: str) -> dict[str, Any]:
  payload = request_json(
    session,
    "GET",
    MYGENE,
    params={
      "q": f"symbol:{gene} AND species:human",
      "fields": "symbol,name,entrezgene,ensembl.gene,uniprot.Swiss-Prot",
      "size": 5,
    },
  )
  hits = payload.get("hits", []) if isinstance(payload, dict) else []
  exact = [hit for hit in hits if canonical_gene_name(hit.get("symbol")) == gene]
  hit = exact[0] if exact else (hits[0] if hits else {})
  ensembl = hit.get("ensembl")
  if isinstance(ensembl, list):
    ensembl_ids = [item.get("gene") for item in ensembl if isinstance(item, dict)]
  elif isinstance(ensembl, dict):
    ensembl_ids = [ensembl.get("gene")]
  else:
    ensembl_ids = []
  uniprot = hit.get("uniprot", {})
  swiss = uniprot.get("Swiss-Prot") if isinstance(uniprot, dict) else None
  if isinstance(swiss, str):
    swiss = [swiss]
  return {
    "symbol": canonical_gene_name(hit.get("symbol")) or gene,
    "name": hit.get("name"),
    "entrez_id": hit.get("entrezgene"),
    "ensembl_ids": [value for value in ensembl_ids if value],
    "uniprot_ids": swiss or [],
  }


def open_targets_rows(
  session: requests.Session,
  gene: str,
  ensembl_id: str | None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
  if not ensembl_id:
    return [], {"reason": "no Ensembl identifier"}
  query = """
  query TargetPharmacology($ensemblId: String!) {
    target(ensemblId: $ensemblId) {
      id
      approvedSymbol
      tractability { label modality value }
      knownDrugs(size: 200) {
        count
        rows {
          drugId
          prefName
          drugType
          mechanismOfAction
          phase
          status
          targetClass
          disease { id name }
        }
      }
    }
  }
  """
  payload = request_json(
    session,
    "POST",
    OPEN_TARGETS,
    json={"query": query, "variables": {"ensemblId": ensembl_id}},
  )
  if payload.get("errors"):
    fallback = """
    query TargetTractability($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        id approvedSymbol tractability { label modality value }
      }
    }
    """
    payload = request_json(
      session,
      "POST",
      OPEN_TARGETS,
      json={"query": fallback, "variables": {"ensemblId": ensembl_id}},
    )
  target = (payload.get("data") or {}).get("target") or {}
  tractability = target.get("tractability") or []
  tractability_score = normalize_open_targets_tractability(tractability)
  rows = []
  known_drugs = target.get("knownDrugs") or {}
  for drug in known_drugs.get("rows") or []:
    disease = drug.get("disease") or {}
    rows.append({
      "source": "open_targets",
      "target_gene": gene,
      "target_id": ensembl_id,
      "drug_id": drug.get("drugId"),
      "drug_name": drug.get("prefName"),
      "interaction_type": drug.get("mechanismOfAction"),
      "drug_type": drug.get("drugType"),
      "max_phase": drug.get("phase"),
      "status": drug.get("status"),
      "target_class": drug.get("targetClass"),
      "disease_id": disease.get("id"),
      "disease_name": disease.get("name"),
      "tractability_score": tractability_score,
      "raw_record": json.dumps(drug, sort_keys=True),
    })
  if not rows and tractability:
    rows.append({
      "source": "open_targets",
      "target_gene": gene,
      "target_id": ensembl_id,
      "drug_id": None,
      "drug_name": None,
      "tractability_score": tractability_score,
      "raw_record": json.dumps({"tractability": tractability}, sort_keys=True),
    })
  return rows, payload


def dgidb_rows(session: requests.Session, genes: list[str]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
  query = """
  query GeneInteractions($names: [String!]!) {
    genes(names: $names) {
      nodes {
        name
        longName
        conceptId
        geneCategoriesWithSources { name sourceNames }
        interactions {
          drug { name conceptId approved immunotherapy antiNeoplastic }
          interactionScore
          interactionTypes { type directionality }
          interactionAttributes { name value }
          publications { pmid }
          sources { sourceDbName }
        }
      }
    }
  }
  """
  payload = request_json(
    session,
    "POST",
    DGIDB,
    json={"query": query, "variables": {"names": genes}},
  )
  if payload.get("errors"):
    raise RuntimeError(f"DGIdb GraphQL errors: {payload['errors']}")
  rows = []
  nodes = (((payload.get("data") or {}).get("genes") or {}).get("nodes") or [])
  for node in nodes:
    gene = canonical_gene_name(node.get("name"))
    categories = node.get("geneCategoriesWithSources") or []
    for interaction in node.get("interactions") or []:
      drug = interaction.get("drug") or {}
      interaction_types = interaction.get("interactionTypes") or []
      sources = interaction.get("sources") or []
      publications = interaction.get("publications") or []
      rows.append({
        "source": "dgidb",
        "target_gene": gene,
        "target_id": node.get("conceptId"),
        "drug_id": drug.get("conceptId"),
        "drug_name": drug.get("name"),
        "interaction_score": interaction.get("interactionScore"),
        "interaction_normalized": normalize_dgidb_interaction_score(
          interaction.get("interactionScore")
        ),
        "interaction_type": ";".join(
          sorted({str(item.get("type")) for item in interaction_types if item.get("type")})
        ),
        "directionality": ";".join(
          sorted({str(item.get("directionality")) for item in interaction_types if item.get("directionality")})
        ),
        "approved": drug.get("approved"),
        "anti_neoplastic": drug.get("antiNeoplastic"),
        "immunotherapy": drug.get("immunotherapy"),
        "sources": ";".join(
          sorted({str(item.get("sourceDbName")) for item in sources if item.get("sourceDbName")})
        ),
        "pmids": ";".join(
          sorted({str(item.get("pmid")) for item in publications if item.get("pmid")})
        ),
        "gene_categories": ";".join(
          sorted({str(item.get("name")) for item in categories if item.get("name")})
        ),
        "raw_record": json.dumps(interaction, sort_keys=True),
      })
  return rows, payload


def _chembl_exact_target(targets: Iterable[dict[str, Any]], gene: str) -> dict[str, Any] | None:
  gene = canonical_gene_name(gene)
  candidates = []
  for target in targets:
    if str(target.get("organism", "")).casefold() != "homo sapiens":
      continue
    names = {canonical_gene_name(target.get("pref_name"))}
    for component in target.get("target_components") or []:
      names.add(canonical_gene_name(component.get("accession")))
      for synonym in component.get("target_component_synonyms") or []:
        names.add(canonical_gene_name(synonym.get("component_synonym")))
    priority = 0 if gene in names else 1
    candidates.append((priority, target))
  return sorted(candidates, key=lambda item: item[0])[0][1] if candidates else None


def chembl_rows(
  session: requests.Session,
  gene: str,
  activity_limit: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
  target_payload = request_json(
    session,
    "GET",
    f"{CHEMBL}/target/search.json",
    params={"q": gene, "limit": 30},
  )
  target = _chembl_exact_target(target_payload.get("targets") or [], gene)
  if not target:
    return [], target_payload
  target_id = target.get("target_chembl_id")
  mechanism_payload = request_json(
    session,
    "GET",
    f"{CHEMBL}/mechanism.json",
    params={"target_chembl_id": target_id, "limit": 1000},
  )
  activity_payload = request_json(
    session,
    "GET",
    f"{CHEMBL}/activity.json",
    params={
      "target_chembl_id": target_id,
      "pchembl_value__isnull": "false",
      "limit": activity_limit,
    },
  )
  mechanism_by_molecule: dict[str, list[dict[str, Any]]] = {}
  for mechanism in mechanism_payload.get("mechanisms") or []:
    mechanism_by_molecule.setdefault(str(mechanism.get("molecule_chembl_id")), []).append(mechanism)
  rows = []
  for activity in activity_payload.get("activities") or []:
    molecule = str(activity.get("molecule_chembl_id") or "")
    mechanisms = mechanism_by_molecule.get(molecule) or [None]
    for mechanism in mechanisms:
      mechanism = mechanism or {}
      rows.append({
        "source": "chembl",
        "target_gene": gene,
        "target_id": target_id,
        "drug_id": molecule or None,
        "drug_name": mechanism.get("molecule_name"),
        "interaction_type": mechanism.get("mechanism_of_action"),
        "action_type": mechanism.get("action_type"),
        "mechanism_of_action": mechanism.get("mechanism_of_action"),
        "max_phase": mechanism.get("max_phase"),
        "pchembl_value": activity.get("pchembl_value"),
        "standard_type": activity.get("standard_type"),
        "standard_value": activity.get("standard_value"),
        "standard_units": activity.get("standard_units"),
        "assay_chembl_id": activity.get("assay_chembl_id"),
        "document_chembl_id": activity.get("document_chembl_id"),
        "raw_record": json.dumps({"activity": activity, "mechanism": mechanism}, sort_keys=True),
      })
  for molecule, mechanisms in mechanism_by_molecule.items():
    if any(str(row.get("drug_id")) == molecule for row in rows):
      continue
    for mechanism in mechanisms:
      rows.append({
        "source": "chembl",
        "target_gene": gene,
        "target_id": target_id,
        "drug_id": molecule,
        "drug_name": mechanism.get("molecule_name"),
        "interaction_type": mechanism.get("mechanism_of_action"),
        "action_type": mechanism.get("action_type"),
        "mechanism_of_action": mechanism.get("mechanism_of_action"),
        "max_phase": mechanism.get("max_phase"),
        "raw_record": json.dumps({"mechanism": mechanism}, sort_keys=True),
      })
  return rows, {
    "target": target_payload,
    "mechanism": mechanism_payload,
    "activity": activity_payload,
  }


def pharos_rows(session: requests.Session, genes: list[str]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
  query = """
  query TargetBatch($targets: [String!]!) {
    batch(targets: $targets) {
      targetResult {
        count
        targets {
          sym
          name
          description
          tdl
          idgFamily
          ligandCounts { name value }
          ligands(skip: 0, top: 100, isdrug: true) {
            name
            isdrug
            smiles
            synonyms { name value }
          }
        }
      }
    }
  }
  """
  payload = request_json(
    session,
    "POST",
    PHAROS,
    json={"query": query, "variables": {"targets": genes}},
  )
  if payload.get("errors"):
    raise RuntimeError(f"Pharos GraphQL errors: {payload['errors']}")
  result = (((payload.get("data") or {}).get("batch") or {}).get("targetResult") or {})
  rows = []
  for target in result.get("targets") or []:
    gene = canonical_gene_name(target.get("sym"))
    tdl = target.get("tdl")
    ligands = target.get("ligands") or []
    if not ligands:
      rows.append({
        "source": "pharos",
        "target_gene": gene,
        "target_development_level": tdl,
        "target_tractability_normalized": normalize_pharos_tdl(tdl),
        "target_family": target.get("idgFamily"),
        "raw_record": json.dumps(target, sort_keys=True),
      })
    for ligand in ligands:
      rows.append({
        "source": "pharos",
        "target_gene": gene,
        "drug_name": ligand.get("name"),
        "target_development_level": tdl,
        "target_tractability_normalized": normalize_pharos_tdl(tdl),
        "target_family": target.get("idgFamily"),
        "smiles": ligand.get("smiles"),
        "raw_record": json.dumps({"target": target, "ligand": ligand}, sort_keys=True),
      })
  return rows, payload


def civic_row(session: requests.Session, gene: str) -> dict[str, Any]:
  response = session.get(
    f"{CIVIC_LINK}/{gene}",
    timeout=60,
    allow_redirects=False,
  )
  if response.status_code in {301, 302, 303, 307, 308}:
    return {
      "source": "civic",
      "target_gene": gene,
      "civic_gene_record": True,
      "civic_url": response.headers.get("Location"),
    }
  if response.status_code == 404:
    return {
      "source": "civic",
      "target_gene": gene,
      "civic_gene_record": False,
    }
  response.raise_for_status()
  return {
    "source": "civic",
    "target_gene": gene,
    "civic_gene_record": True,
    "civic_url": response.url,
  }


def mychem_annotations(
  session: requests.Session,
  compound_ids: list[str],
  batch_size: int = 500,
) -> list[dict[str, Any]]:
  rows = []
  for start in range(0, len(compound_ids), batch_size):
    batch = compound_ids[start:start + batch_size]
    if not batch:
      continue
    payload = request_json(
      session,
      "POST",
      f"{MYCHEM}/chem",
      data={
        "ids": ",".join(batch),
        "fields": "chembl,drugcentral,pubchem,chebi,ginas,pharmgkb,sider",
      },
    )
    if isinstance(payload, dict):
      payload = [payload]
    for record in payload or []:
      rows.append({
        "source": "mychem",
        "drug_id": record.get("query") or record.get("_id"),
        "drug_name": (
          (record.get("chembl") or {}).get("pref_name")
          if isinstance(record.get("chembl"), dict)
          else None
        ),
        "mychem_id": record.get("_id"),
        "found": record.get("notfound") is not True,
        "raw_record": json.dumps(record, sort_keys=True),
      })
  return rows


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--output-dir",
    default="data/processed/pharmacology",
  )
  parser.add_argument(
    "--cache-dir",
    default="data/raw/pharmacology/api_cache",
  )
  parser.add_argument("--minimum-vulnerability-score", type=float, default=0.15)
  parser.add_argument("--max-targets", type=int, default=500)
  parser.add_argument("--chembl-activity-limit", type=int, default=1000)
  parser.add_argument("--sleep", type=float, default=0.10)
  parser.add_argument("--refresh", action="store_true")
  parser.add_argument(
    "--sources",
    default="open_targets,chembl,dgidb,mychem,pharos,civic",
  )
  parser.add_argument(
    "--strict-source",
    action="append",
    default=[],
    help="Fail if a named source cannot be acquired; repeatable.",
  )
  args = parser.parse_args()

  ranking_path = resolve_path(args.ranking)
  output_dir = resolve_path(args.output_dir)
  cache_dir = resolve_path(args.cache_dir)
  output_dir.mkdir(parents=True, exist_ok=True)
  cache_dir.mkdir(parents=True, exist_ok=True)
  ranking = pd.read_csv(ranking_path, sep="\t")
  targets = select_targets(
    ranking,
    args.minimum_vulnerability_score,
    args.max_targets,
  )
  targets.to_csv(output_dir / "pharmacology_target_selection.tsv", sep="\t", index=False)
  genes = targets["target_gene"].astype(str).tolist()
  requested = {value.strip() for value in args.sources.split(",") if value.strip()}
  strict = set(args.strict_source)

  session = requests.Session()
  session.headers.update({
    "User-Agent": "RSES-Onco/0.9 pharmacology evidence client",
    "Accept": "application/json",
  })
  evidence_rows: list[dict[str, Any]] = []
  status_rows: list[dict[str, Any]] = []
  identifiers: dict[str, dict[str, Any]] = {}

  for index, gene in enumerate(genes, start=1):
    cache = cache_dir / "mygene" / f"{safe_name(gene)}.json"
    try:
      identifiers[gene] = cached_json(
        cache,
        lambda gene=gene: mygene_identifiers(session, gene),
        args.refresh,
      )
      status_rows.append({"source": "mygene", "target_gene": gene, "status": "ok"})
    except Exception as exc:
      identifiers[gene] = {"symbol": gene, "ensembl_ids": [], "uniprot_ids": []}
      status_rows.append({"source": "mygene", "target_gene": gene, "status": "failed", "message": str(exc)})
    print(f"[Identifiers {index}/{len(genes)}] {gene}", flush=True)

  if "open_targets" in requested:
    for index, gene in enumerate(genes, start=1):
      ensembl_ids = identifiers.get(gene, {}).get("ensembl_ids") or []
      ensembl_id = ensembl_ids[0] if ensembl_ids else None
      cache = cache_dir / "open_targets" / f"{safe_name(gene)}.json"
      try:
        rows, payload = open_targets_rows(session, gene, ensembl_id)
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        evidence_rows.extend(rows)
        status_rows.append({"source": "open_targets", "target_gene": gene, "status": "ok", "rows": len(rows)})
      except Exception as exc:
        status_rows.append({"source": "open_targets", "target_gene": gene, "status": "failed", "message": str(exc)})
        if "open_targets" in strict:
          raise
      print(f"[Open Targets {index}/{len(genes)}] {gene}", flush=True)
      time.sleep(args.sleep)

  if "chembl" in requested:
    for index, gene in enumerate(genes, start=1):
      cache = cache_dir / "chembl" / f"{safe_name(gene)}.json"
      try:
        rows, payload = chembl_rows(session, gene, args.chembl_activity_limit)
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        evidence_rows.extend(rows)
        status_rows.append({"source": "chembl", "target_gene": gene, "status": "ok", "rows": len(rows)})
      except Exception as exc:
        status_rows.append({"source": "chembl", "target_gene": gene, "status": "failed", "message": str(exc)})
        if "chembl" in strict:
          raise
      print(f"[ChEMBL {index}/{len(genes)}] {gene}", flush=True)
      time.sleep(args.sleep)

  if "dgidb" in requested:
    for start in range(0, len(genes), 50):
      batch = genes[start:start + 50]
      cache = cache_dir / "dgidb" / f"batch_{start:05d}.json"
      try:
        rows, payload = dgidb_rows(session, batch)
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        evidence_rows.extend(rows)
        for gene in batch:
          count = sum(1 for row in rows if row.get("target_gene") == gene)
          status_rows.append({"source": "dgidb", "target_gene": gene, "status": "ok", "rows": count})
      except Exception as exc:
        for gene in batch:
          status_rows.append({"source": "dgidb", "target_gene": gene, "status": "failed", "message": str(exc)})
        if "dgidb" in strict:
          raise
      print(f"[DGIdb {min(start + 50, len(genes))}/{len(genes)}]", flush=True)
      time.sleep(args.sleep)

  if "pharos" in requested:
    for start in range(0, len(genes), 50):
      batch = genes[start:start + 50]
      cache = cache_dir / "pharos" / f"batch_{start:05d}.json"
      try:
        rows, payload = pharos_rows(session, batch)
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        evidence_rows.extend(rows)
        for gene in batch:
          count = sum(1 for row in rows if row.get("target_gene") == gene)
          status_rows.append({"source": "pharos", "target_gene": gene, "status": "ok", "rows": count})
      except Exception as exc:
        for gene in batch:
          status_rows.append({"source": "pharos", "target_gene": gene, "status": "failed", "message": str(exc)})
        if "pharos" in strict:
          raise
      print(f"[Pharos {min(start + 50, len(genes))}/{len(genes)}]", flush=True)
      time.sleep(args.sleep)

  if "civic" in requested:
    for index, gene in enumerate(genes, start=1):
      try:
        row = civic_row(session, gene)
        evidence_rows.append(row)
        status_rows.append({"source": "civic", "target_gene": gene, "status": "ok", "rows": 1})
      except Exception as exc:
        status_rows.append({"source": "civic", "target_gene": gene, "status": "failed", "message": str(exc)})
        if "civic" in strict:
          raise
      print(f"[CIViC {index}/{len(genes)}] {gene}", flush=True)
      time.sleep(max(args.sleep, 0.35))

  evidence = pd.DataFrame(evidence_rows)
  if "mychem" in requested and not evidence.empty:
    compound_ids = sorted({
      str(value)
      for value in evidence.get("drug_id", pd.Series(dtype=object)).dropna()
      if str(value).upper().startswith("CHEMBL")
    })
    try:
      mychem_rows = mychem_annotations(session, compound_ids)
      evidence_rows.extend(mychem_rows)
      status_rows.append({"source": "mychem", "target_gene": "ALL", "status": "ok", "rows": len(mychem_rows)})
    except Exception as exc:
      status_rows.append({"source": "mychem", "target_gene": "ALL", "status": "failed", "message": str(exc)})
      if "mychem" in strict:
        raise

  evidence = pd.DataFrame(evidence_rows)
  if not evidence.empty:
    evidence["target_gene"] = evidence.get("target_gene", pd.Series(index=evidence.index, dtype=object)).map(canonical_gene_name)
    evidence = evidence.drop_duplicates()
  evidence.to_csv(output_dir / "pharmacology_evidence_long.tsv", sep="\t", index=False)
  pd.DataFrame(status_rows).to_csv(output_dir / "pharmacology_source_status.tsv", sep="\t", index=False)
  pd.DataFrame([
    {"target_gene": gene, **identifiers.get(gene, {})}
    for gene in genes
  ]).to_csv(output_dir / "pharmacology_target_identifiers.tsv", sep="\t", index=False)
  metadata = {
    "ranking": str(ranking_path),
    "target_count": len(genes),
    "minimum_vulnerability_score": args.minimum_vulnerability_score,
    "maximum_targets": args.max_targets,
    "requested_sources": sorted(requested),
    "api_endpoints": {
      "open_targets": OPEN_TARGETS,
      "chembl": CHEMBL,
      "dgidb": DGIDB,
      "mychem": MYCHEM,
      "pharos": PHAROS,
      "civic": CIVIC_LINK,
    },
    "civic_api_key_present": bool(os.environ.get("CIVIC_API_KEY")),
  }
  (output_dir / "pharmacology_acquisition_metadata.json").write_text(
    json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8"
  )
  print(f"Selected targets: {len(genes):,}")
  print(f"Evidence rows: {len(evidence):,}")
  print(f"Wrote pharmacology evidence to {output_dir}")


if __name__ == "__main__":
  main()
