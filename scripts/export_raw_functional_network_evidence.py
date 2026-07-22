#!/usr/bin/env python3
"""Export raw functional-network, localization and annotation evidence.

The exporter filters already acquired STRING, DoRothEA/OmniPath, Human Protein Atlas
and UniProtKB tables to genes present in the candidate universe. It preserves source
channels and never describes STRING combined scores as direct experiments or
DoRothEA associations as promoter binding without a separate direct source.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError

ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Record:
  evidence_family: str
  output_path: str
  source_path: str
  status: str
  rows: int
  columns: int
  sha256: str
  interpretation_boundary: str


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_optional(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    return pd.DataFrame()
  try:
    return pd.read_csv(path, sep="\t", low_memory=False)
  except EmptyDataError:
    return pd.DataFrame()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> Path:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)
  return path


def sha256(path: Path) -> str:
  digest = hashlib.sha256()
  with path.open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def first_column(frame: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
  return next((column for column in candidates if column in frame), None)


def clean_gene(series: pd.Series) -> pd.Series:
  return series.fillna("").astype(str).str.strip().str.upper()


def candidate_genes(frame: pd.DataFrame) -> set[str]:
  genes: set[str] = set()
  for column in (
    "analysis_lost_gene", "lost_gene", "analysis_target_gene", "target_gene",
  ):
    if column in frame:
      genes.update(
        clean_gene(frame[column]).loc[lambda values: values.ne("")].tolist()
      )
  return genes


def export(
  family: str,
  frame: pd.DataFrame,
  source_path: Path,
  output_path: Path,
  boundary: str,
) -> Record:
  if frame.empty:
    frame = pd.DataFrame([{
      "evidence_status": "source_absent_or_no_candidate_records",
      "source_path": str(source_path),
      "scientific_interpretation": "Unavailable evidence is not negative biological evidence.",
    }])
    status = "source_absent_or_no_candidate_records"
  else:
    status = "available"
  atomic_tsv(frame, output_path)
  return Record(
    evidence_family=family,
    output_path=str(output_path),
    source_path=str(source_path),
    status=status,
    rows=len(frame),
    columns=len(frame.columns),
    sha256=sha256(output_path),
    interpretation_boundary=boundary,
  )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--raw-dir",
    default="data/raw/human_functional_evidence",
  )
  parser.add_argument("--output-root", default="article_outputs")
  args = parser.parse_args()

  ranking_path = resolve_path(args.ranking)
  ranking = pd.read_csv(ranking_path, sep="\t", low_memory=False)
  genes = candidate_genes(ranking)
  if not genes:
    raise ValueError("No gene identifiers were available in the ranking")

  raw_dir = resolve_path(args.raw_dir)
  output_root = resolve_path(args.output_root)
  network_dir = output_root / "tables/supporting_evidence/networks/raw_sources"
  localization_dir = output_root / "tables/supporting_evidence/localization"
  structure_dir = output_root / "tables/supporting_evidence/structures"
  records: list[Record] = []
  node_frames: list[pd.DataFrame] = []

  string_path = raw_dir / "string_interaction_partners.tsv"
  string = read_optional(string_path)
  if not string.empty:
    source_column = first_column(
      string,
      ("preferredName_A", "source_genesymbol", "query_gene", "source"),
    )
    target_column = first_column(
      string,
      ("preferredName_B", "target_genesymbol", "preferredName", "target"),
    )
    if source_column and target_column:
      string["source_gene"] = clean_gene(string[source_column])
      string["target_gene"] = clean_gene(string[target_column])
      string = string.loc[
        string["source_gene"].isin(genes)
        | string["target_gene"].isin(genes)
      ].copy()
      string["source_name_column"] = source_column
      string["target_name_column"] = target_column
      string["combined_score_interpretation"] = (
        "STRING combined score integrates multiple channels and is not direct experimental evidence."
      )
      channel_columns = [
        column for column in (
          "nscore", "fscore", "pscore", "ascore", "escore", "dscore",
          "tscore", "score", "combined_score",
        ) if column in string
      ]
      string["available_evidence_channels"] = ";".join(channel_columns)
      node_frames.append(string[["source_gene"]].rename(columns={"source_gene": "gene"}))
      node_frames.append(string[["target_gene"]].rename(columns={"target_gene": "gene"}))
  records.append(export(
    "STRING_raw_candidate_edges",
    string,
    string_path,
    network_dir / "string_candidate_edges_all_channels.tsv",
    "Experimental, curated, coexpression, neighborhood, cooccurrence, fusion and text-mining channels remain distinguishable; the combined score is not called experimental evidence.",
  ))

  dorothea_path = raw_dir / "omnipath_dorothea.tsv"
  dorothea = read_optional(dorothea_path)
  if not dorothea.empty:
    source_column = first_column(
      dorothea,
      ("source_genesymbol", "source", "tf"),
    )
    target_column = first_column(
      dorothea,
      ("target_genesymbol", "target", "gene"),
    )
    if source_column and target_column:
      dorothea["transcription_factor"] = clean_gene(dorothea[source_column])
      dorothea["regulated_gene"] = clean_gene(dorothea[target_column])
      dorothea = dorothea.loc[
        dorothea["transcription_factor"].isin(genes)
        | dorothea["regulated_gene"].isin(genes)
      ].copy()
      dorothea["regulatory_evidence_type"] = "TF_target_association"
      dorothea["promoter_binding_evidence"] = "not_available_from_this_table"
      dorothea["interpretation_boundary"] = (
        "DoRothEA association is not direct promoter-binding evidence unless a separate traceable promoter source is present."
      )
      node_frames.append(dorothea[["transcription_factor"]].rename(columns={"transcription_factor": "gene"}))
      node_frames.append(dorothea[["regulated_gene"]].rename(columns={"regulated_gene": "gene"}))
  records.append(export(
    "DoRothEA_raw_candidate_edges",
    dorothea,
    dorothea_path,
    network_dir / "dorothea_candidate_regulatory_edges.tsv",
    "Regulatory direction and confidence are retained; promoter binding is not claimed from an inferred TF-target association.",
  ))

  hpa_path = raw_dir / "hpa_subcellular_location.tsv"
  hpa = read_optional(hpa_path)
  if not hpa.empty:
    gene_column = first_column(hpa, ("Gene name", "gene_name", "gene"))
    if gene_column:
      hpa["candidate_gene"] = clean_gene(hpa[gene_column])
      hpa = hpa.loc[hpa["candidate_gene"].isin(genes)].copy()
      hpa["localization_interpretation"] = (
        "Presence of an annotation does not automatically imply spatial compatibility or a maximal localization component."
      )
      node_frames.append(hpa[["candidate_gene"]].rename(columns={"candidate_gene": "gene"}))
  records.append(export(
    "Human_Protein_Atlas_candidate_localization",
    hpa,
    hpa_path,
    localization_dir / "hpa_candidate_localization.tsv",
    "Localization reliability and compartments are retained; no maximal score is assigned merely because an annotation exists.",
  ))

  uniprot_path = raw_dir / "uniprot_reviewed_annotations.tsv"
  uniprot = read_optional(uniprot_path)
  if not uniprot.empty:
    gene_column = first_column(
      uniprot,
      ("Gene Names (primary)", "Gene Names", "Gene name", "gene_name", "gene"),
    )
    if gene_column:
      uniprot["candidate_gene"] = (
        uniprot[gene_column].fillna("").astype(str).str.split().str[0].str.upper()
      )
      uniprot = uniprot.loc[uniprot["candidate_gene"].isin(genes)].copy()
      uniprot["structure_interpretation"] = (
        "Absence of an experimental structure is missing structural evidence, not demonstrated structural divergence."
      )
      node_frames.append(uniprot[["candidate_gene"]].rename(columns={"candidate_gene": "gene"}))
  records.append(export(
    "UniProtKB_candidate_annotations",
    uniprot,
    uniprot_path,
    structure_dir / "uniprot_candidate_annotations.tsv",
    "Molecular functions, reactions, domains and structure cross-references remain source-bounded.",
  ))

  nodes = (
    pd.concat(node_frames, ignore_index=True)
      .dropna()
      .drop_duplicates()
      .sort_values("gene")
    if node_frames
    else pd.DataFrame(columns=["gene"])
  )
  nodes["in_candidate_universe"] = nodes["gene"].isin(genes)
  records.append(export(
    "functional_network_nodes",
    nodes,
    ranking_path,
    network_dir / "functional_network_nodes.tsv",
    "Nodes are derived only from candidate identifiers and acquired evidence tables.",
  ))

  promoter_status = pd.DataFrame([{
    "evidence_family": "promoter_or_regulatory_element_binding",
    "status": "not_available_in_current_pipeline_sources",
    "reason": (
      "The current acquired DoRothEA table contains TF-target associations but no independently traceable promoter-region binding table."
    ),
    "scientific_rule": (
      "No promoter or direct TF-binding claim is generated without a separate direct source."
    ),
  }])
  records.append(export(
    "promoter_evidence_status",
    promoter_status,
    dorothea_path,
    network_dir / "promoter_evidence_status.tsv",
    "Absence of direct promoter evidence is explicit and is not converted into negative regulation evidence.",
  ))

  manifest = pd.DataFrame([asdict(record) for record in records])
  manifest["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
  manifest_path = network_dir / "raw_functional_evidence_manifest.tsv"
  atomic_tsv(manifest, manifest_path)
  for value in manifest["output_path"]:
    path = Path(value)
    if not path.exists() or path.stat().st_size == 0:
      raise RuntimeError(f"Missing raw functional evidence output: {path}")
  print(manifest[["evidence_family", "status", "rows", "output_path"]].to_string(index=False))
  print(f"Wrote {manifest_path}")


if __name__ == "__main__":
  main()
