#!/usr/bin/env python3
"""Acquire functional-microniche evidence with resilient DoRothEA fallback.

This recovery-oriented entry point reuses a complete STRING per-gene acquisition,
tries multiple official OmniPath service addresses for DoRothEA, accepts an optional
local DoRothEA TSV, and preserves a persistent OmniPath outage as explicit missing
regulatory coverage unless strict mode is requested.
"""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from rses_onco.networks import (
  direct_string_score,
  localization_pair_metrics,
  regulatory_pair_metrics,
  string_pair_metrics,
)
from rses_onco.utils import canonical_gene_name
from scripts.download_human_functional_evidence import (
  STRING_REQUIRED_COLUMNS,
  atomic_write_frame,
  download_hpa,
  download_string,
  download_uniprot,
  request_with_retries,
  resolve_path,
  write_json,
)

OMNIPATH_ENDPOINTS = (
  "https://omnipathdb.org/interactions",
  "https://omnipathdb.org/interactions/",
  "https://omnipath.org/interactions",
  "https://omnipath.org/interactions/",
)
DOROTHEA_COLUMNS = (
  "source",
  "target",
  "source_genesymbol",
  "target_genesymbol",
  "is_directed",
  "is_stimulation",
  "is_inhibition",
  "consensus_direction",
  "consensus_stimulation",
  "consensus_inhibition",
  "sources",
  "references",
  "dorothea_level",
)
DOROTHEA_SOURCE_COLUMNS = ("source_genesymbol", "source", "tf")
DOROTHEA_TARGET_COLUMNS = ("target_genesymbol", "target", "gene")


def validate_dorothea(frame: pd.DataFrame) -> None:
  """Require a non-empty TF-target table with recognized source/target columns."""
  source_present = any(column in frame.columns for column in DOROTHEA_SOURCE_COLUMNS)
  target_present = any(column in frame.columns for column in DOROTHEA_TARGET_COLUMNS)
  if frame.empty or not source_present or not target_present:
    raise ValueError(
      "unexpected DoRothEA schema: "
      f"rows={len(frame)}, columns={frame.columns.tolist()[:20]}"
    )


def read_valid_dorothea(path: Path) -> pd.DataFrame:
  frame = pd.read_csv(path, sep="\t")
  validate_dorothea(frame)
  return frame


def empty_dorothea() -> pd.DataFrame:
  return pd.DataFrame({column: pd.Series(dtype=object) for column in DOROTHEA_COLUMNS})


def acquire_dorothea(
  output: Path,
  status_output: Path,
  *,
  local_file: Path | None,
  retries: int,
  refresh: bool,
  strict: bool,
) -> tuple[pd.DataFrame, dict[str, Any]]:
  """Acquire DoRothEA from cache, local TSV or multiple official OmniPath hosts."""
  output.parent.mkdir(parents=True, exist_ok=True)
  attempts: list[dict[str, Any]] = []

  if local_file is not None:
    frame = read_valid_dorothea(local_file)
    atomic_write_frame(frame, output)
    summary = {
      "status": "local_file",
      "available": True,
      "rows": len(frame),
      "source": str(local_file),
      "attempts": attempts,
    }
    write_json(summary, status_output)
    print(f"DoRothEA: local file ({len(frame):,} rows)", flush=True)
    return frame, summary

  if output.exists() and not refresh:
    try:
      frame = read_valid_dorothea(output)
      summary = {
        "status": "cache",
        "available": True,
        "rows": len(frame),
        "source": str(output),
        "attempts": attempts,
      }
      write_json(summary, status_output)
      print(f"DoRothEA: cache ({len(frame):,} rows)", flush=True)
      return frame, summary
    except Exception as exc:
      attempts.append({
        "endpoint": str(output),
        "organism_parameter": None,
        "status": "invalid_cache",
        "detail": str(exc),
      })

  session = requests.Session()
  session.headers.update({
    "User-Agent": "RSES-Onco/0.10.3 DoRothEA acquisition",
    "Accept": "text/tab-separated-values, text/plain, application/json",
  })
  common = {
    "datasets": "dorothea",
    "genesymbols": 1,
    "format": "tsv",
    "dorothea_levels": "A,B,C",
    "fields": "sources,references,dorothea_level",
  }

  for endpoint in OMNIPATH_ENDPOINTS:
    for organism_key in ("organisms", "organism"):
      params = {**common, organism_key: 9606}
      try:
        response = request_with_retries(
          session,
          "GET",
          endpoint,
          retries=retries,
          timeout=240,
          params=params,
        )
        text = response.text.strip()
        if not text or text.lstrip().lower().startswith(("<html", "<!doctype")):
          raise ValueError("empty or HTML OmniPath response")
        frame = pd.read_csv(io.StringIO(text), sep="\t")
        validate_dorothea(frame)
        atomic_write_frame(frame, output)
        attempts.append({
          "endpoint": endpoint,
          "organism_parameter": organism_key,
          "status": "ok",
          "rows": len(frame),
        })
        summary = {
          "status": "downloaded",
          "available": True,
          "rows": len(frame),
          "source": endpoint,
          "organism_parameter": organism_key,
          "attempts": attempts,
        }
        write_json(summary, status_output)
        print(
          f"DoRothEA: downloaded from {endpoint} ({len(frame):,} rows)",
          flush=True,
        )
        return frame, summary
      except Exception as exc:
        attempts.append({
          "endpoint": endpoint,
          "organism_parameter": organism_key,
          "status": "failed",
          "detail": str(exc),
        })
        print(
          f"DoRothEA endpoint failed: {endpoint} [{organism_key}]: {exc}",
          flush=True,
        )

  frame = empty_dorothea()
  atomic_write_frame(frame, output)
  summary = {
    "status": "unavailable",
    "available": False,
    "rows": 0,
    "source": None,
    "attempts": attempts,
    "scientific_interpretation": (
      "Regulatory-network evidence is missing because all official OmniPath "
      "DoRothEA endpoints were unavailable; missingness is not scored as zero."
    ),
  }
  write_json(summary, status_output)
  if strict:
    raise RuntimeError(
      "DoRothEA remained unavailable after all official endpoint attempts; "
      f"inspect {status_output}"
    )
  print(
    "WARNING: DoRothEA unavailable; continuing with explicit missing regulatory "
    "coverage. Re-run with --refresh-dorothea when OmniPath recovers.",
    flush=True,
  )
  return frame, summary


def load_complete_string(
  genes: list[str],
  string_path: Path,
  status_path: Path,
) -> tuple[pd.DataFrame | None, dict[str, Any] | None]:
  """Reuse STRING aggregate output only when status covers all current genes."""
  if not string_path.exists() or not status_path.exists():
    return None, None
  try:
    edges = pd.read_csv(string_path, sep="\t")
    status = pd.read_csv(status_path, sep="\t")
  except Exception:
    return None, None
  if "query_gene" not in status or "status" not in status:
    return None, None
  observed = set(status["query_gene"].map(canonical_gene_name))
  expected = set(genes)
  if not expected.issubset(observed):
    return None, None
  failures = status["status"].isin(["mapping_request_failed", "request_failed"])
  if failures.any():
    return None, None
  for column in STRING_REQUIRED_COLUMNS:
    if column not in edges:
      edges[column] = pd.Series(dtype=object)
  summary = {
    "discovery_status": "complete_existing_acquisition",
    "edge_rows": len(edges),
    "status_counts": status["status"].value_counts(dropna=False).to_dict(),
    "candidate_gene_count": len(genes),
    "status_output": str(status_path),
  }
  return edges, summary


def load_or_download_table(
  path: Path,
  downloader,
  *,
  refresh: bool,
) -> tuple[pd.DataFrame, str]:
  if path.exists() and not refresh:
    try:
      return pd.read_csv(path, sep="\t"), "cache"
    except Exception:
      pass
  return downloader(), "downloaded"


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
  parser.add_argument("--string-sleep", type=float, default=1.0)
  parser.add_argument("--string-retries", type=int, default=7)
  parser.add_argument("--string-map-chunk-size", type=int, default=200)
  parser.add_argument("--caller-identity", default="rses-onco")
  parser.add_argument("--refresh-string", action="store_true")
  parser.add_argument("--dorothea-retries", type=int, default=3)
  parser.add_argument("--refresh-dorothea", action="store_true")
  parser.add_argument("--strict-dorothea", action="store_true")
  parser.add_argument("--dorothea-file", default=None)
  parser.add_argument("--refresh-hpa", action="store_true")
  parser.add_argument("--refresh-uniprot", action="store_true")
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
  string_mapping_path = raw_dir / "string_id_mapping.tsv"
  string_status_path = raw_dir / "string_acquisition_status.tsv"
  string_cache_dir = raw_dir / "string_cache"
  dorothea_path = raw_dir / "omnipath_dorothea.tsv"
  dorothea_status_path = raw_dir / "omnipath_dorothea_status.json"
  hpa_path = raw_dir / "hpa_subcellular_location.tsv"
  uniprot_path = raw_dir / "uniprot_reviewed_annotations.tsv"
  metadata_path = raw_dir / "source_metadata.json"

  string_edges: pd.DataFrame | None = None
  string_summary: dict[str, Any] | None = None
  if not args.refresh_string:
    string_edges, string_summary = load_complete_string(
      genes,
      string_path,
      string_status_path,
    )
  if string_edges is None or string_summary is None:
    string_edges, string_summary = download_string(
      genes,
      string_path,
      string_mapping_path,
      string_status_path,
      string_cache_dir,
      args.string_required_score,
      args.string_limit,
      args.caller_identity,
      args.string_sleep,
      args.string_retries,
      args.string_map_chunk_size,
      args.refresh_string,
      True,
    )
  else:
    print(
      f"STRING: reused complete acquisition ({len(string_edges):,} edge rows)",
      flush=True,
    )

  local_dorothea = resolve_path(args.dorothea_file) if args.dorothea_file else None
  dorothea, dorothea_summary = acquire_dorothea(
    dorothea_path,
    dorothea_status_path,
    local_file=local_dorothea,
    retries=args.dorothea_retries,
    refresh=args.refresh_dorothea,
    strict=args.strict_dorothea,
  )

  hpa, hpa_status = load_or_download_table(
    hpa_path,
    lambda: download_hpa(hpa_path),
    refresh=args.refresh_hpa,
  )
  uniprot, uniprot_status = load_or_download_table(
    uniprot_path,
    lambda: download_uniprot(genes, uniprot_path),
    refresh=args.refresh_uniprot,
  )

  rows: list[dict[str, object]] = []
  regulatory_available = bool(dorothea_summary.get("available"))
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
      "regulatory_source_available": regulatory_available,
      "regulatory_source_status": dorothea_summary.get("status"),
    })

  result = pd.DataFrame(rows)
  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  atomic_write_frame(result, output)

  metadata = {
    "string": string_summary,
    "string_species": 9606,
    "string_required_score": args.string_required_score,
    "string_limit_per_query": args.string_limit,
    "string_mapping_output": str(string_mapping_path),
    "string_status_output": str(string_status_path),
    "dorothea": dorothea_summary,
    "dorothea_status_output": str(dorothea_status_path),
    "dorothea_dataset": "dorothea",
    "dorothea_confidence_levels": ["A", "B", "C"],
    "hpa_status": hpa_status,
    "uniprot_status": uniprot_status,
    "candidate_gene_count": len(genes),
    "string_edge_rows": len(string_edges),
    "dorothea_rows": len(dorothea),
    "hpa_rows": len(hpa),
    "uniprot_rows": len(uniprot),
  }
  write_json(metadata, metadata_path)

  print(f"Candidate genes: {len(genes):,}")
  print(f"STRING edges: {len(string_edges):,}")
  print(
    f"DoRothEA interactions: {len(dorothea):,} "
    f"(status={dorothea_summary.get('status')})"
  )
  print(f"HPA localization records: {len(hpa):,} ({hpa_status})")
  print(f"UniProt reviewed records: {len(uniprot):,} ({uniprot_status})")
  print(f"Wrote {len(result):,} pair evidence rows to {output}")
  print(f"Wrote {metadata_path}")
  print(f"Wrote {dorothea_status_path}")


if __name__ == "__main__":
  main()
