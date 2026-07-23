from __future__ import annotations

from dataclasses import dataclass
from math import pi
from typing import Mapping

import numpy as np

CHROMOSOME_LENGTHS_GRCH38: dict[str, int] = {
  "1": 248956422,
  "2": 242193529,
  "3": 198295559,
  "4": 190214555,
  "5": 181538259,
  "6": 170805979,
  "7": 159345973,
  "8": 145138636,
  "9": 138394717,
  "10": 133797422,
  "11": 135086622,
  "12": 133275309,
  "13": 114364328,
  "14": 107043718,
  "15": 101991189,
  "16": 90338345,
  "17": 83257441,
  "18": 80373285,
  "19": 58617616,
  "20": 64444167,
  "21": 46709983,
  "22": 50818468,
  "X": 156040895,
  "Y": 57227415,
}
CHROMOSOME_ORDER = tuple(CHROMOSOME_LENGTHS_GRCH38)


def normalize_chromosome(value: object) -> str:
  text = str(value or "").strip().upper()
  if text.startswith("CHR"):
    text = text[3:]
  aliases = {"23": "X", "24": "Y"}
  return aliases.get(text, text)


@dataclass(frozen=True)
class ChromosomeArc:
  chromosome: str
  length: int
  theta_start: float
  theta_end: float

  @property
  def theta_mid(self) -> float:
    return (self.theta_start + self.theta_end) / 2.0


class GenomeLayout:
  """Map GRCh38 genomic coordinates to a clockwise circular layout."""

  def __init__(
    self,
    chromosome_lengths: Mapping[str, int] | None = None,
    gap_degrees: float = 1.35,
    start_angle: float = pi / 2.0,
  ) -> None:
    lengths = dict(
      chromosome_lengths or CHROMOSOME_LENGTHS_GRCH38
    )
    self.lengths = {
      chromosome: int(lengths[chromosome])
      for chromosome in CHROMOSOME_ORDER
      if chromosome in lengths
    }
    total = float(sum(self.lengths.values()))
    gap = np.deg2rad(float(gap_degrees))
    available = 2.0 * pi - gap * len(self.lengths)
    current = float(start_angle)
    arcs: dict[str, ChromosomeArc] = {}
    for chromosome, length in self.lengths.items():
      span = available * (float(length) / total)
      arcs[chromosome] = ChromosomeArc(
        chromosome,
        length,
        current,
        current - span,
      )
      current = current - span - gap
    self.arcs = arcs

  def theta(self, chromosome: object, position: object) -> float:
    chrom = normalize_chromosome(chromosome)
    if chrom not in self.arcs:
      raise KeyError(f"Unsupported chromosome: {chromosome}")
    arc = self.arcs[chrom]
    pos = float(position)
    fraction = float(np.clip(pos / arc.length, 0.0, 1.0))
    return arc.theta_start + fraction * (
      arc.theta_end - arc.theta_start
    )

  @staticmethod
  def point(theta: float, radius: float) -> tuple[float, float]:
    return (
      radius * float(np.cos(theta)),
      radius * float(np.sin(theta)),
    )
