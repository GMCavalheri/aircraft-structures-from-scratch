"""Rainflow cycle counting, ASTM E1049-85 (reapproved 2017) Sec. 5.4.4.

The history is first reduced to its reversals (turning points). Each reversal is pushed on a
stack. While the newest range X = |S[-1] - S[-2]| is at least the previous range
Y = |S[-2] - S[-3]|, Y is counted: as a half cycle if it contains the starting point (it is the
first range on the stack), otherwise as a full cycle. Ranges left on the stack at the end are
half cycles.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Cycle:
    range: float
    mean: float
    count: float  # 0.5 or 1.0
    start: int  # index of the starting reversal in the original series
    end: int

    @property
    def amplitude(self) -> float:
        return self.range / 2


def reversals(series) -> tuple[np.ndarray, np.ndarray]:
    """Turning points of ``series`` (values, indices into the original series). The end points
    are always kept; repeated values (plateaus) collapse to their first point."""
    x = np.asarray(series, float)
    keep = np.r_[True, np.diff(x) != 0] if x.size else np.zeros(0, bool)
    v, i = x[keep], np.flatnonzero(keep)
    if v.size < 3:
        return v, i
    d = np.sign(np.diff(v))
    turns = np.flatnonzero(d[1:] != d[:-1]) + 1
    sel = np.r_[0, turns, v.size - 1]
    return v[sel], i[sel]


def count_cycles(series) -> list[Cycle]:
    vals, idx = reversals(series)
    stack: list[tuple[float, int]] = []
    cycles: list[Cycle] = []

    def make(a, b, count):
        (va, ia), (vb, ib) = a, b
        return Cycle(abs(vb - va), (va + vb) / 2, count, ia, ib)

    for v, i in zip(vals, idx, strict=True):
        stack.append((float(v), int(i)))
        while len(stack) >= 3:
            X = abs(stack[-1][0] - stack[-2][0])
            Y = abs(stack[-2][0] - stack[-3][0])
            if X < Y:
                break
            if len(stack) == 3:
                cycles.append(make(stack[0], stack[1], 0.5))
                stack.pop(0)
            else:
                cycles.append(make(stack[-3], stack[-2], 1.0))
                del stack[-3:-1]
    for a, b in zip(stack[:-1], stack[1:], strict=False):
        cycles.append(make(a, b, 0.5))
    return cycles


def histogram(cycles: list[Cycle]) -> list[tuple[float, float]]:
    """(range, total count) sorted by range, the summary format of ASTM E1049 Sec. 5.4.4."""
    out: dict[float, float] = {}
    for c in cycles:
        out[c.range] = out.get(c.range, 0.0) + c.count
    return sorted(out.items())


def range_mean_matrix(cycles: list[Cycle], range_bins, mean_bins) -> np.ndarray:
    """2D histogram of counts over (range, mean) bins."""
    H, _, _ = np.histogram2d(
        [c.range for c in cycles],
        [c.mean for c in cycles],
        bins=[range_bins, mean_bins],
        weights=[c.count for c in cycles],
    )
    return H
