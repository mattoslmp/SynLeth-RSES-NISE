from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]


def test_mutation_standardizer_keeps_only_lof_or_annotated_damaging(
  tmp_path: Path,
) -> None:
  models = pd.DataFrame({
    "ModelID": [
      "ACH-000001",
      "ACH-000002",
      "ACH-000003",
      "ACH-000004",
    ],
    "CCLEName": ["A", "B", "C", "D"],
  })
  models_path = tmp_path / "Model.csv"
  models.to_csv(models_path, index=False)

  mutations = pd.DataFrame({
    "ModelID": models["ModelID"],
    "HugoSymbol": ["GENEA", "GENEB", "GENEC", "GENED"],
    "Variant_Classification": [
      "Frame_Shift_Del",
      "Missense_Mutation",
      "Missense_Mutation",
      "Silent",
    ],
    "isDeleterious": [False, True, False, True],
    "ProteinChange": ["p.A10fs", "p.V20E", "p.G30A", "p.L40L"],
  })
  data_dir = tmp_path / "dmap_data"
  data_dir.mkdir()
  mutation_path = data_dir / "mutations.csv"
  mutations.to_csv(mutation_path, index=False)

  config = {
    "version": "test",
    "sources": {
      "mutation_table": {
        "filename": mutation_path.name,
        "role": "scored_variant_level_functional_loss",
        "layout": "long_mutation",
      }
    },
  }
  config_path = tmp_path / "sources.yaml"
  config_path.write_text(
    yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
  )

  matrix_path = tmp_path / "mutation_matrix.tsv"
  event_path = tmp_path / "mutation_events.tsv"
  status_path = tmp_path / "mutation_status.json"
  runtime_path = tmp_path / "runtime.yaml"
  subprocess.run([
    sys.executable,
    str(ROOT / "scripts/standardize_extended_mutation_table.py"),
    "--config",
    str(config_path),
    "--input-dir",
    str(data_dir),
    "--models",
    str(models_path),
    "--output",
    str(matrix_path),
    "--event-output",
    str(event_path),
    "--status-output",
    str(status_path),
    "--runtime-config-output",
    str(runtime_path),
  ], cwd=ROOT, check=True)

  events = pd.read_csv(event_path, sep="\t")
  assert set(events["gene"]) == {"GENEA", "GENEB"}
  assert int(events["clear_loss_of_function"].sum()) == 1
  assert int(
    events["annotated_damaging_missense_or_inframe"].sum()
  ) == 1

  matrix = pd.read_csv(matrix_path, sep="\t")
  assert {"GENEA", "GENEB"}.issubset(matrix.columns)
  assert "GENEC" not in matrix.columns
  assert "GENED" not in matrix.columns

  runtime = yaml.safe_load(runtime_path.read_text(encoding="utf-8"))
  assert runtime["sources"]["mutation_table"]["layout"] == "matrix"
  assert runtime["sources"]["mutation_table"]["gene_features"] is True
