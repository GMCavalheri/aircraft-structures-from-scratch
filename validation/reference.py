"""Load the cited reference values in ``validation/reference_data``."""

from __future__ import annotations

import csv
from pathlib import Path

DATA = Path(__file__).resolve().parent / "reference_data"


def load(name: str) -> list[dict[str, str]]:
    with open(DATA / name) as f:
        rows = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(rows))


def ramberg_osgood_n(material: str, direction: str, loading: str) -> tuple[float, str]:
    for row in load("ramberg_osgood.csv"):
        if (row["material"], row["direction"], row["loading"]) == (material, direction, loading):
            return float(row["n"]), row["source"]
    raise KeyError((material, direction, loading))
