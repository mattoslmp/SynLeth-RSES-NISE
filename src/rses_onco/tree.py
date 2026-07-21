from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pandas as pd
import requests


def download_uniprot_fasta(accessions: list[str], output: str | Path) -> Path:
  """Download reviewed UniProt sequences for a supplied accession list."""
  output = Path(output)
  output.parent.mkdir(parents=True, exist_ok=True)
  query = " OR ".join(f"accession:{a}" for a in accessions)
  response = requests.get(
    "https://rest.uniprot.org/uniprotkb/stream",
    params={"query": f"({query}) AND organism_id:9606", "format": "fasta"},
    timeout=120,
  )
  response.raise_for_status()
  output.write_bytes(response.content)
  return output


def build_sequence_tree(
  fasta: str | Path,
  output_prefix: str | Path,
  threads: int = 4,
) -> tuple[Path, Path]:
  """Run MAFFT and FastTree when installed, returning alignment and Newick paths."""
  if not shutil.which("mafft"):
    raise RuntimeError("MAFFT is not installed")
  if not shutil.which("FastTree") and not shutil.which("fasttree"):
    raise RuntimeError("FastTree is not installed")
  fasttree = shutil.which("FastTree") or shutil.which("fasttree")
  prefix = Path(output_prefix)
  prefix.parent.mkdir(parents=True, exist_ok=True)
  alignment = prefix.with_suffix(".aligned.fasta")
  tree = prefix.with_suffix(".nwk")
  with alignment.open("w", encoding="utf-8") as handle:
    subprocess.run(["mafft", "--thread", str(threads), "--auto", str(fasta)], check=True, stdout=handle)
  with tree.open("w", encoding="utf-8") as handle:
    subprocess.run([fasttree, "-wag", str(alignment)], check=True, stdout=handle)
  return alignment, tree


def catalog_accessions(table: str | Path) -> list[str]:
  frame = pd.read_csv(table, sep="\t")
  return sorted(frame["uniprot_accession"].dropna().astype(str).unique())
