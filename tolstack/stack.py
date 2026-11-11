"""1D tolerance stack-up analysis: worst case, RSS and Monte Carlo.

A stack is a closed loop of dimensions that adds up to a gap (or interference).
Each dimension has a direction (+1 if it opens the gap, -1 if it closes it),
a nominal value and a possibly asymmetric tolerance (+plus / -minus).

Asymmetric tolerances are converted to an equivalent mean-shifted symmetric
tolerance: mean = nominal + (plus - minus) / 2, half = (plus + minus) / 2.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Dimension:
    name: str
    nominal: float
    plus: float
    minus: float          # magnitude, e.g. 0.05 for "-0.05"
    direction: int = 1    # +1 opens the gap, -1 closes it
    dist: str = "normal"  # "normal" (tolerance = +/-3 sigma) or "uniform"

    def __post_init__(self):
        if self.direction not in (1, -1):
            raise ValueError(f"{self.name}: direction must be +1 or -1")
        if self.plus < 0 or self.minus < 0:
            raise ValueError(f"{self.name}: give tolerances as non-negative magnitudes")
        if self.dist not in ("normal", "uniform"):
            raise ValueError(f"{self.name}: dist must be 'normal' or 'uniform'")

    @property
    def mean(self) -> float:
        return self.nominal + (self.plus - self.minus) / 2

    @property
    def half_tol(self) -> float:
        return (self.plus + self.minus) / 2

    @property
    def sigma(self) -> float:
        if self.dist == "normal":
            return self.half_tol / 3
        return self.half_tol / math.sqrt(3)   # std of a uniform distribution on +/- half_tol


@dataclass
class StackResult:
    method: str
    mean: float
    low: float
    high: float
    sigma: float | None = None
    yield_fraction: float | None = None
    cpk: float | None = None

    def __str__(self) -> str:
        s = f"{self.method:<12} mean {self.mean:9.4f}   range [{self.low:9.4f}, {self.high:9.4f}]"
        if self.yield_fraction is not None:
            s += f"   yield {self.yield_fraction:7.3%}"
        if self.cpk is not None:
            s += f"   Cpk {self.cpk:5.2f}"
        return s


class Stack:
    def __init__(self, dims: list[Dimension]):
        if not dims:
            raise ValueError("stack needs at least one dimension")
        self.dims = list(dims)

    @property
    def nominal(self) -> float:
        return sum(d.direction * d.nominal for d in self.dims)

    @property
    def mean(self) -> float:
        return sum(d.direction * d.mean for d in self.dims)

    def worst_case(self) -> StackResult:
        t = sum(d.half_tol for d in self.dims)
        return StackResult("worst_case", self.mean, self.mean - t, self.mean + t)

    def rss(self, lsl: float | None = None, usl: float | None = None) -> StackResult:
        """Root sum square, statistical +/-3 sigma result assuming independent normal parts."""
        sigma = math.sqrt(sum(d.sigma**2 for d in self.dims))
        res = StackResult("rss", self.mean, self.mean - 3 * sigma, self.mean + 3 * sigma, sigma=sigma)
        res.cpk = _cpk(self.mean, sigma, lsl, usl)
        res.yield_fraction = _normal_yield(self.mean, sigma, lsl, usl)
        return res

    def monte_carlo(self, n: int = 200_000, seed: int | None = 0,
                    lsl: float | None = None, usl: float | None = None) -> StackResult:
        rng = np.random.default_rng(seed)
        total = np.zeros(n)
        for d in self.dims:
            if d.dist == "normal":
                x = rng.normal(d.mean, d.sigma, n)
            else:
                x = rng.uniform(d.mean - d.half_tol, d.mean + d.half_tol, n)
            total += d.direction * x
        mean, sigma = float(total.mean()), float(total.std(ddof=1))
        lo, hi = np.percentile(total, [0.135, 99.865])   # +/-3 sigma equivalent coverage
        ok = np.ones(n, dtype=bool)
        if lsl is not None:
            ok &= total >= lsl
        if usl is not None:
            ok &= total <= usl
        return StackResult("monte_carlo", mean, float(lo), float(hi), sigma=sigma,
                           yield_fraction=float(ok.mean()) if (lsl is not None or usl is not None) else None,
                           cpk=_cpk(mean, sigma, lsl, usl))

    def contributions(self) -> list[tuple[str, float]]:
        """Percent contribution of each dimension to the RSS variance, largest first."""
        var = [d.sigma**2 for d in self.dims]
        total = sum(var)
        return sorted(((d.name, 100 * v / total) for d, v in zip(self.dims, var)), key=lambda t: -t[1])


def _cpk(mean: float, sigma: float, lsl: float | None, usl: float | None) -> float | None:
    if sigma == 0 or (lsl is None and usl is None):
        return None
    sides = []
    if usl is not None:
        sides.append((usl - mean) / (3 * sigma))
    if lsl is not None:
        sides.append((mean - lsl) / (3 * sigma))
    return min(sides)


def _phi(z: float) -> float:
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))


def _normal_yield(mean: float, sigma: float, lsl: float | None, usl: float | None) -> float | None:
    if lsl is None and usl is None:
        return None
    hi = _phi((usl - mean) / sigma) if usl is not None else 1.0
    lo = _phi((lsl - mean) / sigma) if lsl is not None else 0.0
    return hi - lo


def load_csv(path: str | Path) -> Stack:
    """CSV columns: name, nominal, plus, minus, direction, dist (dist optional)."""
    dims = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            dims.append(Dimension(
                name=row["name"].strip(),
                nominal=float(row["nominal"]),
                plus=float(row["plus"]),
                minus=abs(float(row["minus"])),
                direction=int(row["direction"]),
                dist=(row.get("dist") or "normal").strip(),
            ))
    return Stack(dims)
