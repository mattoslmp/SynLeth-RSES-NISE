#!/usr/bin/env python3
"""Aggregate downloaded GDC ASCAT3 gene-level integer CN into homdel matrices.

Output is deletion-only discrete coding: 0 total copies -> -2; >0 copies -> 0;
missing remains missing. This is compatible with the RSES-Onco homozygous-
deletion frequency function but must not be described as a GISTIC matrix.

The GDC files query is restricted to ``Primary Tumor``. A case can nevertheless
contain several sample records, including matched normal samples. The aggregator
therefore selects the sample whose ``sample_type`` is exactly ``Primary Tumor``
instead of taking the first sample attached to the case. The file-to-sample
selection is exported as a provenance table.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANCER_PROJECTS = {
  "colon": ["TCGA-COAD", "TCGA-READ"],
  "stomach": ["TCGA-STAD"],
  "lung": ["TCGA-LUAD", "TCGA-LUSC"],
}


@dataclass(frozen=True)
class SampleSelection:
  sample_id: str
  sample_type: str
  status: str
  primary_tumor_candidates: str
  all_sample_records: str


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def sample_records(hit: dict) -> list[tuple[str, str]]:
  """Return unique ``(submitter_id, sample_type)`` records from a GDC hit."""
  records: list[tuple[str, str]] = []
  seen: set[tuple[str, str]] = set()
  for case in hit.get("cases") or []:
    for sample in case.get("samples") or []:
      sample_id = str(sample.get("submitter_id") or "").strip()
      sample_type = str(sample.get("sample_type") or "").strip()
      if not sample_id:
        continue
      record = (sample_id, sample_type)
      if record not in seen:
        seen.add(record)
        records.append(record)
  return records


def select_primary_tumor_sample(hit: dict) -> SampleSelection:
  """Select the Primary Tumor barcode associated with a GDC file.

  The nested ``cases.samples`` response can contain matched normals even when the
  file query itself was filtered to Primary Tumor. When more than one primary
  tumor barcode is available, prefer the barcode contained in the file name;
  otherwise choose deterministically and record the ambiguity in provenance.
  """
  records = sample_records(hit)
  primary = sorted({
    sample_id
    for sample_id, sample_type in records
    if sample_type.casefold() == "primary tumor"
  })
  file_name = str(hit.get("file_name") or "")

  if len(primary) == 1:
    selected = primary[0]
    status = "primary_tumor_exact"
  elif len(primary) > 1:
    filename_matches = [sample_id for sample_id in primary if sample_id in file_name]
    if len(filename_matches) == 1:
      selected = filename_matches[0]
      status = "primary_tumor_filename_match"
    else:
      selected = primary[0]
      status = "multiple_primary_tumors_lexicographic"
  else:
    case_ids = [
      str(case.get("submitter_id") or "").strip()
      for case in hit.get("cases") or []
      if str(case.get("submitter_id") or "").strip()
    ]
    fallback = case_ids[0] if case_ids else str(hit.get("file_id") or "").strip()
    raise ValueError(
      "No Primary Tumor sample record was available for GDC file "
      f"{file_name or hit.get('file_id')}; fallback identifier would be {fallback}"
    )

  all_records = ";".join(
    f"{sample_id}|{sample_type or 'unknown'}"
    for sample_id, sample_type in records
  )
  return SampleSelection(
    sample_id=selected,
    sample_type="Primary Tumor",
    status=status,
    primary_tumor_candidates=";".join(primary),
    all_sample_records=all_records,
  )


def sample_identifier(hit: dict) -> str:
  """Return the explicitly selected Primary Tumor TCGA sample barcode."""
  return select_primary_tumor_sample(hit).sample_id


def choose_column(columns, candidates):
  lookup = {str(c).casefold(): c for c in columns}
  for candidate in candidates:
    if candidate.casefold() in lookup:
      return lookup[candidate.casefold()]
  return None


def read_gene_copy_number(path: Path) -> pd.Series:
  try:
    frame = pd.read_csv(path, sep="\t", comment="#")
  except Exception:
    frame = pd.read_csv(path, sep=None, engine="python", comment="#")
  gene_col = choose_column(frame.columns, [
    "gene_name", "gene_symbol", "hugo_symbol", "Gene Symbol", "gene",
  ])
  cn_col = choose_column(frame.columns, ["copy_number", "Copy_Number", "copy number"])
  if gene_col is None or cn_col is None:
    raise ValueError(
      f"Could not identify gene/copy_number columns in {path}; "
      f"columns={list(frame.columns)}"
    )
  genes = (
    frame[gene_col]
      .astype(str)
      .str.replace(r"\..*$", "", regex=True)
      .str.upper()
      .str.strip()
  )
  values = pd.to_numeric(frame[cn_col], errors="coerce")
  data = pd.DataFrame({"gene": genes, "copy_number": values}).dropna()
  data = data.loc[data["gene"].ne("")]
  # If a gene is represented more than once, retain the minimum total copy number.
  return data.groupby("gene", sort=False)["copy_number"].min()


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--raw-dir", default="data/raw/gdc")
  parser.add_argument(
    "--manifest",
    default="data/raw/gdc/gdc_gene_level_copy_number_manifest.json",
  )
  parser.add_argument("--output-dir", default="data/processed")
  parser.add_argument(
    "--provenance-output",
    default="data/processed/gdc_sample_provenance.tsv",
  )
  args = parser.parse_args()
  raw_dir = resolve_path(args.raw_dir)
  manifest_path = resolve_path(args.manifest)
  manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
  outdir = resolve_path(args.output_dir)
  provenance_output = resolve_path(args.provenance_output)
  outdir.mkdir(parents=True, exist_ok=True)
  provenance_output.parent.mkdir(parents=True, exist_ok=True)

  project_matrices: dict[str, pd.DataFrame] = {}
  provenance_rows: list[dict[str, object]] = []
  for project, hits in manifest.items():
    columns: dict[str, pd.Series] = {}
    for hit in hits:
      path = raw_dir / project / hit["file_name"]
      if not path.exists():
        raise FileNotFoundError(path)
      selection = select_primary_tumor_sample(hit)
      sample = selection.sample_id
      output_sample = sample
      duplicate_disambiguated = False
      if output_sample in columns:
        output_sample = f"{sample}__{str(hit['file_id'])[:8]}"
        duplicate_disambiguated = True
      columns[output_sample] = read_gene_copy_number(path)
      provenance_rows.append({
        "project": project,
        "file_id": hit.get("file_id"),
        "file_name": hit.get("file_name"),
        "file_size": hit.get("file_size"),
        "md5sum": hit.get("md5sum"),
        **asdict(selection),
        "matrix_column": output_sample,
        "duplicate_disambiguated": duplicate_disambiguated,
        "manifest_path": str(manifest_path),
        "raw_path": str(path),
      })
    if not columns:
      raise ValueError(f"No files available for {project}")
    absolute = pd.concat(columns, axis=1).sort_index()
    discrete = pd.DataFrame(np.nan, index=absolute.index, columns=absolute.columns)
    discrete[absolute == 0] = -2
    discrete[absolute > 0] = 0
    discrete.index.name = "Hugo_Symbol"
    output = outdir / f"{project.replace('-', '_')}_homdel_discrete.tsv"
    discrete.to_csv(output, sep="\t", na_rep="NA")
    project_matrices[project] = discrete
    print(
      f"Wrote {output}: {discrete.shape[0]:,} genes x "
      f"{discrete.shape[1]:,} primary tumors",
      flush=True,
    )

  for cancer, projects in CANCER_PROJECTS.items():
    missing = [project for project in projects if project not in project_matrices]
    if missing:
      raise ValueError(f"Missing projects for {cancer}: {missing}")
    combined = pd.concat([project_matrices[p] for p in projects], axis=1)
    combined = combined.loc[:, ~combined.columns.duplicated()].sort_index()
    combined.index.name = "Hugo_Symbol"
    output = outdir / f"TCGA_{cancer.upper()}_homdel_discrete.tsv"
    combined.to_csv(output, sep="\t", na_rep="NA")
    print(
      f"Wrote {output}: {combined.shape[0]:,} genes x "
      f"{combined.shape[1]:,} primary tumors",
      flush=True,
    )

  provenance = pd.DataFrame(provenance_rows)
  provenance.to_csv(provenance_output, sep="\t", index=False)
  print(f"Wrote {provenance_output}: {len(provenance):,} file mappings", flush=True)


if __name__ == "__main__":
  main()
