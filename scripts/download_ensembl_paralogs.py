#!/usr/bin/env python3
"""Expand the analyzed human gene universe with Ensembl Compara paralogs.

This is an explicit source-bounded definition of the homologous-paralog class:
all human paralogs returned by the selected Ensembl release for the supplied seed
genes. It does not claim that all human functional backups are paralogs.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd
import requests

from rses_onco.utils import canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
ENSEMBL = "https://rest.ensembl.org"


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def get_json(session: requests.Session, url: str, retries: int = 3) -> dict:
  error: Exception | None = None
  for attempt in range(1, retries + 1):
    try:
      response = session.get(url, headers={"Content-Type": "application/json"}, timeout=120)
      response.raise_for_status()
      return response.json()
    except Exception as exc:
      error = exc
      if attempt == retries:
        break
      time.sleep(min(10, 2 ** attempt))
  raise RuntimeError(f"Ensembl request failed: {url}: {error}")


def lookup_symbol(session: requests.Session, ensembl_id: str, cache: dict[str, str]) -> str | None:
  if ensembl_id in cache:
    return cache[ensembl_id]
  payload = get_json(session, f"{ENSEMBL}/lookup/id/{ensembl_id}?expand=0")
  symbol = canonical_gene_name(payload.get("display_name"))
  if symbol:
    cache[ensembl_id] = symbol
  return symbol or None


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--output",
    default="data/raw/ensembl/human_seed_paralogs.tsv",
  )
  parser.add_argument(
    "--cache",
    default="data/raw/ensembl/ensembl_symbol_cache.json",
  )
  parser.add_argument("--sleep", type=float, default=0.15)
  args = parser.parse_args()

  candidates = pd.read_csv(resolve_path(args.candidates), sep="\t")
  genes = sorted({
    canonical_gene_name(value)
    for column in ("lost_gene", "target_gene")
    if column in candidates
    for value in candidates[column].dropna()
    if canonical_gene_name(value)
  })
  output = resolve_path(args.output)
  cache_path = resolve_path(args.cache)
  output.parent.mkdir(parents=True, exist_ok=True)
  cache_path.parent.mkdir(parents=True, exist_ok=True)
  if cache_path.exists():
    cache = json.loads(cache_path.read_text(encoding="utf-8"))
  else:
    cache = {}

  session = requests.Session()
  rows: list[dict[str, object]] = []
  for index, seed in enumerate(genes, start=1):
    payload = get_json(
      session,
      f"{ENSEMBL}/homology/symbol/human/{seed}?type=paralogues;target_species=human;sequence=none",
    )
    homologies = []
    for item in payload.get("data", []):
      homologies.extend(item.get("homologies", []))
    for homology in homologies:
      target = homology.get("target", {})
      target_id = str(target.get("id") or "")
      target_symbol = lookup_symbol(session, target_id, cache) if target_id else None
      if not target_symbol or target_symbol == seed:
        continue
      identity = target.get("perc_id")
      rows.append({
        "pair_id": f"ENSEMBL_PARALOG_{seed}_TO_{target_symbol}",
        "lost_feature": f"{seed} loss",
        "lost_gene": seed,
        "target_gene": target_symbol,
        "source_class": "homologous_paralog",
        "relation_type": "ensembl_compara_paralog",
        "mechanism": (
          f"Ensembl Compara identifies {seed} and {target_symbol} as human paralogs; "
          "the pair is evaluated as a candidate homologous backup, not assumed synthetic lethal."
        ),
        "colon": 1,
        "stomach": 1,
        "lung": 1,
        "relation_confidence": 0.85,
        "genetic_screen": 0.0,
        "isogenic_validation": 0.0,
        "in_vivo": 0.0,
        "clinical_tractability": 0.20,
        "lineage_relevance": 0.0,
        "evidence_stage": "systematic discovery",
        "primary_doi": "",
        "supporting_doi": "",
        "status": "unvalidated Ensembl paralog direction",
        "ensembl_target_id": target_id,
        "ensembl_homology_type": homology.get("type"),
        "ensembl_target_percent_identity": identity,
      })
    cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True), encoding="utf-8")
    print(f"[Ensembl {index}/{len(genes)}] {seed}: {len(homologies)} paralog records", flush=True)
    if args.sleep:
      time.sleep(args.sleep)

  result = pd.DataFrame(rows)
  if not result.empty:
    result = result.drop_duplicates(["lost_gene", "target_gene", "source_class"])
    result = result.sort_values(["lost_gene", "target_gene"])
  result.to_csv(output, sep="\t", index=False)
  print(f"Seed genes: {len(genes):,}")
  print(f"Directed Ensembl paralog candidates: {len(result):,}")
  print(f"Wrote {output}")


if __name__ == "__main__":
  main()
