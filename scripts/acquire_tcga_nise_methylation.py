#!/usr/bin/env python3
"""Acquire candidate-gene TCGA/GDC methylation through UCSC Xena.

The GDC is the source of record. Repbase is not used because it contains
representative repetitive-element sequences rather than sample-level methylation.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import pandas as pd

from rses_onco.methylation import (
  CANCER_PROJECTS,
  METHYLATION_SUBWEIGHTS,
  build_pair_metrics,
  choose_dataset,
  query_gene_frame,
  summarize_genes,
)
from rses_onco.utils import canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_HUB = "https://gdc.xenahubs.net"


def resolve(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)


def atomic_json(payload: dict[str, Any], path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  temporary.write_text(
    json.dumps(payload, indent=2, sort_keys=True),
    encoding="utf-8",
  )
  temporary.replace(path)


def overrides(values: list[str]) -> dict[str, str]:
  parsed = {}
  for value in values:
    project, dataset = value.split("=", 1)
    parsed[project.strip()] = dataset.strip()
  return parsed


def write_failure(
  pairs: pd.DataFrame,
  output_dir: Path,
  reason: str,
  min_samples: int,
) -> None:
  status = {
    cancer: {"status": "technical_failure", "reason": reason}
    for cancer in CANCER_PROJECTS
  }
  metrics = build_pair_metrics(
    pd.DataFrame(), pd.DataFrame(), pairs, status, min_samples
  )
  atomic_tsv(
    metrics,
    output_dir / "tcga_nise_methylation_pair_metrics.tsv",
  )
  atomic_tsv(
    pd.DataFrame(status).T.reset_index(names="cancer"),
    output_dir / "tcga_nise_methylation_source_status.tsv",
  )
  atomic_json(
    {
      "status": "technical_failure",
      "reason": reason,
      "repbase_used": False,
      "cancers": status,
    },
    output_dir / "tcga_nise_methylation_source_status.json",
  )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--output-dir",
    default="data/processed/epigenetics/methylation",
  )
  parser.add_argument("--hub", default=DEFAULT_HUB)
  parser.add_argument(
    "--dataset",
    action="append",
    default=[],
    help="Dataset override as TCGA-COAD=dataset_name",
  )
  parser.add_argument("--min-samples", type=int, default=20)
  parser.add_argument("--strict", action="store_true")
  args = parser.parse_args()

  output_dir = resolve(args.output_dir)
  output_dir.mkdir(parents=True, exist_ok=True)
  candidate_frame = pd.read_csv(
    resolve(args.candidates), sep="\t", low_memory=False
  )
  columns = [
    name
    for name in (
      "pair_id", "lost_gene", "target_gene", "source_class"
    )
    if name in candidate_frame
  ]
  pairs = candidate_frame[columns].copy()
  pairs["lost_gene"] = pairs["lost_gene"].map(canonical_gene_name)
  pairs["target_gene"] = pairs["target_gene"].map(canonical_gene_name)
  pairs = pairs.loc[
    pairs["lost_gene"].ne("") & pairs["target_gene"].ne("")
  ].drop_duplicates("pair_id")
  genes = sorted(set(pairs["lost_gene"]) | set(pairs["target_gene"]))

  try:
    import xenaPython as xena
  except ImportError:
    write_failure(
      pairs,
      output_dir,
      "xenaPython_not_installed",
      args.min_samples,
    )
    if args.strict:
      raise
    print(
      "Methylation acquisition unavailable: xenaPython not installed"
    )
    return

  dataset_overrides = overrides(args.dataset)
  accessed = datetime.now(timezone.utc).isoformat()
  long_frames = []
  status_rows = []
  cancer_status: dict[str, dict[str, Any]] = {}
  for cancer, projects in CANCER_PROJECTS.items():
    available = 0
    reasons = []
    for project in projects:
      dataset, discovery, candidate_count = choose_dataset(
        xena,
        args.hub,
        project,
        dataset_overrides.get(project),
      )
      if not dataset:
        reasons.append(f"{project}:{discovery}")
        status_rows.append({
          "cancer": cancer,
          "project": project,
          "status": "source_unavailable",
          "reason": discovery,
          "dataset": "",
          "hub": args.hub,
          "accessed_at_utc": accessed,
          "candidate_dataset_count": candidate_count,
        })
        continue
      try:
        samples = [
          str(sample)
          for sample in xena.dataset_samples(
            args.hub, dataset, None
          )
        ]
        gene_frames = []
        for gene in genes:
          position, values = xena.dataset_gene_probes_values(
            args.hub, dataset, samples, [gene]
          )
          gene_frames.append(
            query_gene_frame(position, values, samples, gene)
          )
        project_long = (
          pd.concat(gene_frames, ignore_index=True)
          if gene_frames
          else pd.DataFrame()
        )
        if not project_long.empty:
          project_long["cancer"] = cancer
          project_long["project"] = project
          project_long["dataset"] = dataset
          project_long["hub"] = args.hub
          project_long["accessed_at_utc"] = accessed
          project_long["probe_selection"] = (
            "GDC_Xena_gene_associated_CpG_probes"
          )
          long_frames.append(project_long)
        available += 1
        status_rows.append({
          "cancer": cancer,
          "project": project,
          "status": "available",
          "reason": "",
          "dataset": dataset,
          "hub": args.hub,
          "accessed_at_utc": accessed,
          "candidate_dataset_count": candidate_count,
          "dataset_discovery": discovery,
          "sample_count": len(samples),
          "candidate_gene_count": len(genes),
          "retrieved_record_count": len(project_long),
        })
      except Exception as exc:
        reason = f"{type(exc).__name__}:{exc}"
        reasons.append(f"{project}:{reason}")
        status_rows.append({
          "cancer": cancer,
          "project": project,
          "status": "technical_failure",
          "reason": reason,
          "dataset": dataset,
          "hub": args.hub,
          "accessed_at_utc": accessed,
        })
        if args.strict:
          raise
    cancer_status[cancer] = {
      "status": (
        "available" if available else "technical_failure"
      ),
      "reason": ";".join(reasons),
      "available_projects": available,
      "requested_projects": len(projects),
    }

  long = (
    pd.concat(long_frames, ignore_index=True)
    if long_frames
    else pd.DataFrame()
  )
  summaries = summarize_genes(long)
  metrics = build_pair_metrics(
    long,
    summaries,
    pairs,
    cancer_status,
    args.min_samples,
  )
  atomic_tsv(
    long,
    output_dir / "tcga_nise_gene_methylation_long.tsv",
  )
  atomic_tsv(
    summaries,
    output_dir / "tcga_nise_methylation_gene_summary.tsv",
  )
  atomic_tsv(
    metrics,
    output_dir / "tcga_nise_methylation_pair_metrics.tsv",
  )
  atomic_tsv(
    pd.DataFrame(status_rows),
    output_dir / "tcga_nise_methylation_source_status.tsv",
  )
  atomic_json(
    {
      "source_of_record": "NCI Genomic Data Commons",
      "access_layer": "UCSC Xena GDC Hub",
      "hub": args.hub,
      "gdc_beta_definition": "M/(M+U), bounded to [0,1]",
      "repbase_used": False,
      "repbase_exclusion_reason": (
        "Repbase is a repetitive DNA sequence library, not a "
        "sample-level DNA methylation resource."
      ),
      "methylation_subweights": METHYLATION_SUBWEIGHTS,
      "cancers": cancer_status,
      "accessed_at_utc": accessed,
    },
    output_dir / "tcga_nise_methylation_source_status.json",
  )
  print(f"Methylation long rows: {len(long):,}")
  print(f"Methylation gene summaries: {len(summaries):,}")
  print(f"Methylation pair rows: {len(metrics):,}")


if __name__ == "__main__":
  main()
