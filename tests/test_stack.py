import math
from pathlib import Path

import pytest

from tolstack import Dimension, Stack, load_csv

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "hinge_gap.csv"


def simple_stack():
    return Stack([
        Dimension("A", 10.0, 0.1, 0.1, +1),
        Dimension("B", 4.0, 0.1, 0.1, -1),
        Dimension("C", 5.0, 0.1, 0.1, -1),
    ])


def test_nominal_and_worst_case():
    s = simple_stack()
    assert s.nominal == pytest.approx(1.0)
    wc = s.worst_case()
    assert (wc.low, wc.high) == (pytest.approx(0.7), pytest.approx(1.3))


def test_rss_three_equal_tolerances():
    rss = simple_stack().rss()
    assert rss.high - rss.mean == pytest.approx(0.1 * math.sqrt(3))


def test_rss_never_wider_than_worst_case():
    s = load_csv(EXAMPLE)
    assert s.rss().high - s.rss().low <= s.worst_case().high - s.worst_case().low


def test_asymmetric_tolerance_shifts_mean():
    d = Dimension("x", 1.5, 0.15, 0.05)
    assert d.mean == pytest.approx(1.55)
    assert d.half_tol == pytest.approx(0.10)


def test_uniform_sigma():
    d = Dimension("u", 1.0, 0.3, 0.3, dist="uniform")
    assert d.sigma == pytest.approx(0.3 / math.sqrt(3))


def test_monte_carlo_matches_rss_for_normal_parts():
    s = simple_stack()
    mc = s.monte_carlo(n=400_000, seed=1)
    assert mc.mean == pytest.approx(s.mean, abs=1e-3)
    assert mc.sigma == pytest.approx(s.rss().sigma, rel=0.01)


def test_yield_and_cpk_agree_between_methods():
    s = simple_stack()
    rss = s.rss(lsl=0.85, usl=1.15)
    mc = s.monte_carlo(n=400_000, seed=2, lsl=0.85, usl=1.15)
    assert rss.yield_fraction == pytest.approx(mc.yield_fraction, abs=2e-3)
    assert rss.cpk == pytest.approx(mc.cpk, rel=0.02)


def test_cpk_one_sided():
    rss = simple_stack().rss(lsl=0.7)
    assert rss.cpk == pytest.approx((1.0 - 0.7) / (3 * rss.sigma))


def test_contributions_sum_to_100():
    pcts = [p for _, p in load_csv(EXAMPLE).contributions()]
    assert sum(pcts) == pytest.approx(100)
    assert pcts == sorted(pcts, reverse=True)


@pytest.mark.parametrize("kwargs", [
    {"direction": 0}, {"plus": -0.1}, {"dist": "triangular"},
])
def test_invalid_dimension(kwargs):
    base = dict(name="bad", nominal=1.0, plus=0.1, minus=0.1, direction=1)
    base.update(kwargs)
    with pytest.raises(ValueError):
        Dimension(**base)


def test_example_file_loads():
    s = load_csv(EXAMPLE)
    assert len(s.dims) == 6
    assert s.nominal == pytest.approx(0.40)
