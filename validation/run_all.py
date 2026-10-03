"""Regenerate every validation figure and results table.

Run with ``uv run python validation/run_all.py``.
"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = [
    "stress_strain_validation.py",
    "beam_validation.py",
    "fea_validation.py",
    "buckling_validation.py",
    "composites_validation.py",
    "fatigue_validation.py",
]

if __name__ == "__main__":
    sys.path.insert(0, str(HERE))
    for name in SCRIPTS:
        print(f"== {name}")
        runpy.run_path(str(HERE / name), run_name="__main__")
