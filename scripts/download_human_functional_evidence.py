#!/usr/bin/env python3
"""Acquire human functional-microniche evidence for the expanded universe.

Sources:
- STRING functional interaction partners (taxon 9606);
- OmniPath/DoRothEA transcription-factor target interactions;
- Human Protein Atlas subcellular localization;
- UniProtKB reviewed annotations for biochemical/structural traceability.

Raw responses are cached and pair-level divergence metrics are written separately.
Missing resources remain missing and never become score zero.
"""
from __future__ import annotations

import argparse
import io
import json
import time
import zipfile
from pathlib import Path
from typing import Callable

import pandas as pd
import requests

from rses_onco.networks import (
  direct_string_score,
  localization_pair_metrics,
  regulatory_pair_metrics,
  string_pair_metrics,
)
from rses_onco.utils import canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
STRING_ENDPOINT = "https://string-db.org/api/tsv/interaction_partners"
OMNIPATH_ENDPOINT = "https://omnipathdb.org/interactions"
HPA_URL = "https://www.proteinatlas.org/download/subcellular_location.tsv.zip"
UNIPROT_ENDPOINT = "https://rest.uniprot.org/uniprotkb/search"


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def request_with_retries(
  session: requests.Session,
  method: str,
  url: str,
  *,
  retries: int = 3,
  timeout: int = 180,
  **kwargs: object,
) -> requests.Response:
  error: Exception | None = None
  for attempt in range(1, retries + 1):
    try:
      response = session.request(method, url, timeout=timeout, **kwargs)
      response.raise_for_status()
      return response
    except Exception as exc:
      error = exc
      if attempt == retries:
        break
      time.sleep(min(10, 2 ** attempt))
  raise RuntimeError(f"Request failed after {retries} attempts: {url}: {error}")


def download_string(
  genes: list[str],
  output: Path,
  required_score: int,
  limit: int,
  caller_identity: str,
  sleep_seconds: float,
) -> pd.DataFrame:
  output.parent.mkdir(parents=True, exist_ok=True)
  session = requests.Session()
  frames: list[pd.DataFrame] = []
  for index, gene in enumerate(genes, start=1):
    response = request_with_retries(
      session,
      "POST",
      STRING_ENDPOINT,
      data={
        "identifiers": gene,
        "species": 9606,
        "required_score": required_score,
        "limit": limit,
        "network_type": "functional",
        "caller_identity": caller_identity,
      },
    )
    text = response.text.strip()
    if text:
      frame = pd.read_csv(io.StringIO(text), sep="\t")
      frame["query_gene"] = gene
      frames.append(frame)
    print(f"[STRING {index}/{len(genes)}] {gene}", flush=True)
    if sleep_seconds:
      time.sleep(sleep_seconds)
  result = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
  result.to_csv(output, sep="\t", index=False)
  return result


def download_dorothea(output: Path) -> pd.DataFrame:
  """Download high/medium-confidence human DoRothEA interactions.

  OmniPath deployments have accepted both ``organisms`` and ``organism`` in
  different client generations. The downloader tries the plural documented web
  service parameter first, validates the returned schema, and falls back to the
  singular spelling without silently accepting an HTML error page.
  """
  output.parent.mkdir(parents=True, exist_ok=True)
  session = requests.Session()
  common = {
    "datasets": "dorothea",
    "genesymbols": 1,
    "format": "tsv",
    "dorothea_levels": "A,B,C",
    "fields": "sources,references,dorothea_level",
  }
  errors: list[str] = []
  for organism_key in ("organisms", "organism"):
    params = {**common, organism_key: 9606}
    try:
      response = request_with_retries(
        session,
        "GET",
        OMNIPATH_ENDPOINT,
        params=params,
      )
      text = response.text.strip()
      frame = pd.read_csv(io.StringIO(text), sep="\t")
      source_present = any(
        column in frame.columns for column in ("source_genesymbol", "source", "tf")
      )
      target_present = any(
        column in frame.columns for column in ("target_genesymbol", "target", "gene")
      )
      if frame.empty or not source_present or not target_present:
        raise ValueError(
          f"unexpected OmniPath schema: rows={len(frame)}, columns={frame.columns.tolist()[:12]}"
        )
      output.write_text(text + "\n", encoding="utf-8")
      return frame
    except Exception as exc:
      errors.append(f"{organism_key}: {exc}")
  raise RuntimeError("DoRothEA acquisition failed; " + " | ".join(errors))


def download_hpa(output: Path) -> pd.DataFrame:
  output.parent.mkdir(parents=True, exist_ok=True)
  session = requests.Session()
  response = request_with_retries(session, "GET", HPA_URL)
  with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
    members = [
      name for name in archive.namelist()
      if name.lower().endswith((".tsv", ".txt"))
    ]
    if not members:
      raise RuntimeError("HPA archive does not contain a TSV file")
    content = archive.read(members[0])
  output.write_bytes(content)
  return pd.read_csv(output, sep="\t")


def download_uniprot(genes: list[str], output: Path, chunk_size: int = 40) -> pd.DataFrame:
  output.parent.mkdir(parents=True, exist_ok=True)
  session = requests.Session()
  frames: list[pd.DataFrame] = []
  fields = ",".join([
    "accession",
    "id",
    "gene_primary",
    "protein_name",
    "ec",
    "cc_function",
    "cc_catalytic_activity",
    "cc_cofactor",
    "cc_subcellular_location",
    "xref_pdb",
  ])
  for start in range(0, len(genes), chunk_size):
    chunk = genes[start:start + chunk_size]
    gene_query = " OR ".join(f"gene_exact:{gene}" for gene in chunk)
    response = request_with_retries(
      session,
      "GET",
      UNIPROT_ENDPOINT,
      params={
        "query": f"({gene_query}) AND organism_id:9606 AND reviewed:true",
        "format": "tsv",
        "fields": fields,
        "size": 500,
      },
    )
    text = response.text.strip()
    if text:
      frames.append(pd.read_csv(io.StringIO(text), sep="\t"))
    print(f"[UniProt {min(start + chunk_size, len(genes))}/{len(genes)}]", flush=True)
  result = pd.concat(frames, ignore_index=True).drop_duplicates() if frames else pd.DataFrame()
  result.to_csv(output, sep="\t", index=False)
  return result


def load_or_download(
  path: Path,
  downloader: Callable[[], pd.DataFrame],
  skip_download: bool,
) -> pd.DataFrame:
  if skip_download:
    if not path.exists():
      raise FileNotFoundError(f"--skip-download requested but file is absent: {path}")
    return pd.read_csv(path, sep="\t")
  return downloader()


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--raw-dir",
    default="data/raw/human_functional_evidence",
  )
  parser.add_argument(
    "--output",
    default="data/processed/expanded_pair_functional_evidence.tsv",
  )
  parser.add_argument("--string-required-score", type=int, default=700)
  parser.add_argument("--string-limit", type=int, default=100)
  parser.add_argument("--string-sleep", type=float, default=0.15)
  parser.add_argument("--caller-identity", default="rses-onco")
  parser.add_argument("--skip-download", action="store_true")
  args = parser.parse_args()

  candidates = pd.read_csv(resolve_path(args.candidates), sep="\t")
  required = {"pair_id", "lost_gene", "target_gene"}
  missing = sorted(required - set(candidates.columns))
  if missing:
    raise ValueError(f"Candidate universe lacks columns: {missing}")
  genes = sorted({
    canonical_gene_name(value)
    for column in ("lost_gene", "target_gene")
    for value in candidates[column].dropna()
    if canonical_gene_name(value)
  })

  raw_dir = resolve_path(args.raw_dir)
  raw_dir.mkdir(parents=True, exist_ok=True)
  string_path = raw_dir / "string_interaction_partners.tsv"
  dorothea_path = raw_dir / "omnipath_dorothea.tsv"
  hpa_path = raw_dir / "hpa_subcellular_location.tsv"
  uniprot_path = raw_dir / "uniprot_reviewed_annotations.tsv"
  metadata_path = raw_dir / "source_metadata.json"

  string_edges = load_or_download(
    string_path,
    lambda: download_string(
      genes,
      string_path,
      args.string_required_score,
      args.string_limit,
      args.caller_identity,
      args.string_sleep,
    ),
    args.skip_download,
  )
  dorothea = load_or_download(
    dorothea_path,
    lambda: download_dorothea(dorothea_path),
    args.skip_download,
  )
  hpa = load_or_download(
    hpa_path,
    lambda: download_hpa(hpa_path),
    args.skip_download,
  )
  uniprot = load_or_download(
    uniprot_path,
    lambda: download_uniprot(genes, uniprot_path),
    args.skip_download,
  )

  rows: list[dict[str, object]] = []
  for record in candidates.to_dict("records"):
    gene_a = canonical_gene_name(record["lost_gene"])
    gene_b = canonical_gene_name(record["target_gene"])
    string_metrics = string_pair_metrics(string_edges, gene_a, gene_b)
    regulatory_metrics = regulatory_pair_metrics(dorothea, gene_a, gene_b)
    localization_metrics = localization_pair_metrics(hpa, gene_a, gene_b)

    structural_divergence = None
    cluster_a = record.get("lost_structural_cluster")
    cluster_b = record.get("target_structural_cluster")
    if pd.notna(cluster_a) and pd.notna(cluster_b):
      structural_divergence = float(str(cluster_a) != str(cluster_b))

    rows.append({
      "pair_id": record["pair_id"],
      "lost_gene": gene_a,
      "target_gene": gene_b,
      "source_class": record.get("source_class"),
      "component_localization": localization_metrics.divergence,
      "localization_shared": localization_metrics.shared_count,
      "localization_exclusive_lost": localization_metrics.exclusive_a_count,
      "localization_exclusive_target": localization_metrics.exclusive_b_count,
      "component_biochemical_structural": structural_divergence,
      "component_interaction_network": string_metrics.divergence,
      "string_direct_score": direct_string_score(string_edges, gene_a, gene_b),
      "string_neighbor_jaccard": string_metrics.jaccard,
      "string_shared_neighbors": string_metrics.shared_count,
      "component_regulatory_network": regulatory_metrics.divergence,
      "regulator_jaccard": regulatory_metrics.jaccard,
      "shared_regulators": regulatory_metrics.shared_count,
    })

  result = pd.DataFrame(rows)
  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  result.to_csv(output, sep="\t", index=False)
  metadata_path.write_text(
    json.dumps(
      {
        "string_endpoint": STRING_ENDPOINT,
        "string_species": 9606,
        "string_required_score": args.string_required_score,
        "string_limit_per_query": args.string_limit,
        "omnipath_endpoint": OMNIPATH_ENDPOINT,
        "omnipath_dataset": "dorothea",
        "omnipath_confidence_levels": ["A", "B", "C"],
        "hpa_url": HPA_URL,
        "uniprot_endpoint": UNIPROT_ENDPOINT,
        "candidate_gene_count": len(genes),
        "string_edge_rows": len(string_edges),
        "dorothea_rows": len(dorothea),
        "hpa_rows": len(hpa),
        "uniprot_rows": len(uniprot),
      },
      indent=2,
      sort_keys=True,
    ),
    encoding="utf-8",
  )
  print(f"Candidate genes: {len(genes):,}")
  print(f"STRING edges: {len(string_edges):,}")
  print(f"DoRothEA interactions: {len(dorothea):,}")
  print(f"HPA localization records: {len(hpa):,}")
  print(f"UniProt reviewed records: {len(uniprot):,}")
  print(f"Wrote {len(result):,} pair evidence rows to {output}")
  print(f"Wrote {metadata_path}")


if __name__ == "__main__":
  main()
