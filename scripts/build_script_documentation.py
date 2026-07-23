#!/usr/bin/env python3
"""Generate a complete repository script/module catalogue with commands, inputs and outputs."""
from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import re

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SUPPORTED = {".py", ".sh", ".R", ".r"}


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def sha256(path: Path) -> str:
  digest = hashlib.sha256()
  with path.open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def python_docstring(text: str) -> str:
  try:
    tree = ast.parse(text)
    return (ast.get_docstring(tree) or "").strip()
  except SyntaxError:
    return ""


def shell_or_r_description(lines: list[str]) -> str:
  values = []
  for line in lines[:40]:
    stripped = line.strip()
    if stripped.startswith("#!"):
      continue
    if stripped.startswith("#"):
      value = stripped.lstrip("#").strip()
      if value:
        values.append(value)
    elif stripped:
      break
  return " ".join(values)


def purpose(path: Path, text: str) -> str:
  if path.suffix == ".py":
    value = python_docstring(text)
  else:
    value = shell_or_r_description(text.splitlines())
  return (
    re.sub(r"\s+", " ", value).strip()
    or "No module-level description was recorded."
  )


def cli_options(text: str) -> str:
  options = sorted(set(re.findall(
    r"['\"](--[A-Za-z0-9][A-Za-z0-9_-]*)['\"]",
    text,
  )))
  return ";".join(options)


def declared_paths(text: str) -> str:
  pattern = (
    r"(?:data|results|article_outputs|logs|docs|manuscript|"
    r"supplementary|config)/[A-Za-z0-9_./{}$()\-]+"
  )
  values = sorted({
    value.rstrip(".,;:'\")")
    for value in re.findall(pattern, text)
  })
  return ";".join(values[:200])


def stage(path: Path) -> str:
  name = path.name.casefold()
  mapping = [
    (("download", "acquire"), "data_acquisition"),
    (("candidate", "relation_forest"), "candidate_universe"),
    (
      ("wgcna", "regulatory", "methylation", "promoter", "jaspar"),
      "expression_and_regulation",
    ),
    (
      ("score", "rses", "dependency", "expression_compensation"),
      "scoring",
    ),
    (("structure", "alphafold", "pymol"), "structural_analysis"),
    (("pharmacology", "drug", "chembl", "civic"), "pharmacology"),
    (("circos", "genomic"), "genomic_circos"),
    (("figure", "plot"), "figures"),
    (
      ("table", "workbook", "manifest", "catalog"),
      "tables_and_reproducibility",
    ),
    (("validate", "audit", "test"), "validation"),
    (("document", "manuscript", "supplement"), "documentation"),
    (("pipeline", "resume", "run_"), "orchestration"),
  ]
  for terms, label in mapping:
    if any(term in name for term in terms):
      return label
  return "supporting_module"


def command(path: Path) -> str:
  relative = path.relative_to(ROOT).as_posix()
  if path.suffix == ".py":
    return f"python -u {relative} --help"
  if path.suffix.casefold() == ".r":
    return f"Rscript {relative} --help"
  return f"bash {relative}"


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--output-md",
    default="docs/SCRIPT_CATALOG.md",
  )
  parser.add_argument(
    "--output-tsv",
    default="docs/script_manifest.tsv",
  )
  parser.add_argument(
    "--processed-output",
    default=(
      "data/processed/documentation/"
      "pipeline_script_catalog.tsv"
    ),
  )
  args = parser.parse_args()

  paths = sorted([
    path
    for directory in (ROOT / "scripts", ROOT / "src/rses_onco")
    for path in directory.rglob("*")
    if path.is_file()
    and path.suffix in SUPPORTED
    and "__pycache__" not in path.parts
  ])
  rows = []
  generated = datetime.now(timezone.utc).isoformat()
  for path in paths:
    text = path.read_text(encoding="utf-8", errors="replace")
    relative = path.relative_to(ROOT).as_posix()
    rows.append({
      "script_path": relative,
      "language": (
        "Python"
        if path.suffix == ".py"
        else "R"
        if path.suffix.casefold() == ".r"
        else "Bash"
      ),
      "pipeline_stage": stage(path),
      "purpose": purpose(path, text),
      "entrypoint_or_reproduction_command": command(path),
      "cli_options": cli_options(text),
      "declared_input_output_paths": declared_paths(text),
      "line_count": len(text.splitlines()),
      "sha256": sha256(path),
      "documentation_status": "documented_from_source",
      "catalogued_at_utc": generated,
    })
  frame = pd.DataFrame(rows)
  if frame.empty:
    raise RuntimeError("No scripts/modules were discovered")

  outputs = [
    resolve(args.output_tsv),
    resolve(args.processed_output),
  ]
  for output in outputs:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    frame.to_csv(temporary, sep="\t", index=False)
    temporary.replace(output)

  markdown = [
    "# RSES-Onco script and module catalogue",
    "",
    (
      "This catalogue is generated directly from every Python, Bash and R "
      "source file in `scripts/` and `src/rses_onco/`. It is therefore an "
      "executable repository contract rather than a manually curated partial "
      "list."
    ),
    "",
    f"Generated entries: **{len(frame)}**",
    "",
    (
      "| Script/module | Stage | Language | Purpose | Reproduction command "
      "| Declared paths |"
    ),
    "|---|---|---|---|---|---|",
  ]
  for row in frame.to_dict("records"):
    purpose_text = str(row["purpose"]).replace("|", "\\|")
    paths_text = str(
      row["declared_input_output_paths"]
    ).replace("|", "\\|")
    markdown.append(
      f"| `{row['script_path']}` | {row['pipeline_stage']} | "
      f"{row['language']} | {purpose_text} | "
      f"`{row['entrypoint_or_reproduction_command']}` | "
      f"`{paths_text}` |"
    )
  output_md = resolve(args.output_md)
  output_md.parent.mkdir(parents=True, exist_ok=True)
  output_md.write_text(
    "\n".join(markdown) + "\n",
    encoding="utf-8",
  )
  print(
    f"Documented {len(frame)} scripts/modules in "
    f"{output_md} and {outputs[0]}"
  )


if __name__ == "__main__":
  main()
