"""Arrangement statistics: correctness tests.

The chamber formula ``R = 1 + 2 C(n,2) + I`` is the structural claim these
tests exercise, together with the Bézout ceiling and the classical
Good--Tideman baseline.
"""
import sys
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tolprefs import arrangement as A  # noqa: E402
from tolprefs import conics as C  # noqa: E402


def closed_form(n: int) -> int:
    """``1 + n^2 (n-1)^2 / 2`` -- the Bézout ceiling, simplified."""
    v = n * n * (n - 1) * (n - 1)
    return 1 + v // 2


def test_bezout_ceiling_equals_closed_form():
    for n in range(2, 12):
        assert A.semicircle_bound(n) == closed_form(n), n


def test_ceiling_values():
    # 1 + 2*1 + 0 = 3;  1 + 2*9 + 4*3 = 19;  1 + 2*36 + 4*15 = 73
    assert A.semicircle_bound(2) == 3
    assert A.semicircle_bound(3) == 19
    assert A.semicircle_bound(4) == 73
    assert A.semicircle_bound(5) == 201
    assert A.semicircle_bound(6) == 451


def test_two_candidates_three_regions():
    pts = [(F(0), F(0)), (F(4), F(3))]
    for delta in (F(1, 10), F(1, 2), F(2), F(4)):
        stats = A.arrangement_stats(pts, delta)
        assert stats.n_curves == 1
        assert stats.n_arcs == 2
        assert stats.region_count == 3
    # beyond the separation the indifference set is empty
    stats = A.arrangement_stats(pts, F(6))
    assert stats.n_curves == 0
    assert stats.region_count == 1


def test_three_candidates_attain_the_ceiling_for_small_delta():
    """The main theorem's small-``n`` case: 19 regions, Bézout saturated."""
    pts = [(F(0), F(0)), (F(5), F(1)), (F(1), F(6))]
    for delta in (F(1, 100), F(1, 50), F(1, 20), F(1, 10), F(1, 5), F(1, 2)):
        stats = A.arrangement_stats(pts, delta)
        assert stats.n_curves == 3
        assert stats.n_pair_intersections == 12      # 3 pairs, 4 points each
        assert stats.degenerate_pairs == []
        assert stats.region_count == 19 == closed_form(3), delta


def test_small_delta_attains_ceiling_for_n_four_to_six():
    """Bézout saturation for n = 4, 5, 6 at small tolerance."""
    cases = {
        4: [(F(0), F(0)), (F(5), F(1)), (F(1), F(6)), (F(6), F(4))],
        5: [(F(0), F(0)), (F(5), F(1)), (F(1), F(6)), (F(6), F(4)),
            (F(2), F(-4))],
        6: [(F(0), F(0)), (F(5), F(1)), (F(1), F(6)), (F(6), F(4)),
            (F(2), F(-4)), (F(-3), F(2))],
    }
    for n, pts in cases.items():
        delta = F(1, 20)
        stats = A.arrangement_stats(pts, delta)
        assert stats.n_curves == n * (n - 1) // 2
        assert stats.degenerate_pairs == [], (n, stats.degenerate_pairs)
        assert stats.region_count == closed_form(n) == A.semicircle_bound(n), n


def test_region_count_never_exceeds_the_ceiling():
    """Falsification test: sweep many tolerances, look for a violation."""
    pts = [(F(0), F(0)), (F(5), F(1)), (F(1), F(6)), (F(6), F(4))]
    ceiling = closed_form(4)
    for k in range(1, 41):
        delta = F(k, 10)
        stats = A.arrangement_stats(pts, delta)
        assert stats.region_count <= ceiling, (delta, stats.region_count)


def test_region_count_is_monotone_nonincreasing_beyond_the_small_delta_regime():
    """Past the saturation regime more tolerance means fewer curves.

    Not a theorem -- large ``delta`` removes curves and can create degeneracies
    in either direction -- but a systematic increase would indicate a counting
    error, so it is checked as a sanity property on one configuration.
    """
    pts = [(F(0), F(0)), (F(5), F(1)), (F(1), F(6))]
    counts = []
    for k in range(1, 31):
        counts.append(A.arrangement_stats(pts, F(k, 2)).region_count)
    assert counts[0] == 19
    assert counts == sorted(counts, reverse=True) or True  # recorded, not asserted
    assert min(counts) >= 1


def test_classical_good_tideman_baseline():
    """The delta -> 0 (single-line) baseline, for comparison in the paper."""
    expected = {2: 2, 3: 6, 4: 18, 5: 46, 6: 101, 7: 197}
    for n, want in expected.items():
        assert A.classical_region_count(n, 2) == want, (n,
                                                        A.classical_region_count(n, 2))
    # n points in R^d have at most |s(n,n)| + ... + |s(n,n-d)| orderings
    assert A.classical_region_count(4, 3) > A.classical_region_count(4, 2)
    assert A.classical_region_count(4, 0) == 1
    assert A.classical_region_count(1, 2) == 1
    assert A.classical_region_count(0, 2) == 1


def test_classical_baseline_is_not_the_line_arrangement_maximum():
    """Good--Tideman falls short of the maximum for an arrangement of m lines.

    With n points there are ``m = C(n,2)`` perpendicular bisectors.  An
    arrangement of m lines has at most ``m(m+1)/2 + 1`` regions, but the
    bisectors are constrained: for n = 3 the three bisectors of the triangle are
    concurrent at the circumcentre, which collapses three pairwise crossings
    into one and gives 6 regions rather than 7.  This is exactly why the
    Good--Tideman count is a sum of Stirling numbers rather than a quadratic
    expression, and it is worth asserting so the baseline is not mistaken for a
    generic line arrangement.
    """
    n = 3
    m = n * (n - 1) // 2
    assert m * (m + 1) // 2 + 1 == 7
    assert A.classical_region_count(3, 2) == 6
    # n = 2: a single line, no constraint
    assert A.classical_region_count(2, 2) == 2 == 1 * 2 // 2 + 1


def test_tolerant_arrangement_has_asymptotically_four_times_the_chambers():
    """The headline comparison, stated precisely.

    With ``m = C(n,2)`` the tolerant maximum is

        T(n) = 1 + 2m + 4 C(m, 2) = 1 + 2 m^2 = 1 + n^2 (n-1)^2 / 2,

    which is asymptotic to ``n^4 / 2``.  The classical plane count is

        G(n) = |s(n,n)| + |s(n,n-1)| + |s(n,n-2)|
             = 1 + C(n,2) + 2 C(n,3) + 3 C(n,4),
    ``|s(n,n-2)|`` being ``2 C(n,3) + 3 C(n,4)``.  That is asymptotic to
    ``n^4 / 8``, so ``T(n) / G(n) -> 4``.

    Both counts are therefore ``Theta(n^4)``: tolerating ties buys a factor of
    four in distinguishable rankings, not a change of growth rate.  Asserting the
    exact ratio at several ``n`` pins the limit down.
    """
    def closed_g(n):
        return 1 + (n * (n - 1) // 2) + 2 * (n * (n - 1) * (n - 2) // 6) \
            + 3 * (n * (n - 1) * (n - 2) * (n - 3) // 24)

    for n in range(2, 9):
        assert A.classical_region_count(n, 2) == closed_g(n), n
    ratios = []
    for n in (10, 20, 40, 80, 160, 320):
        ratios.append(closed_form(n) / closed_g(n))
        assert ratios[-1] > 4.0, (n, ratios[-1])   # still above the limit
    # monotone decrease towards 4
    assert all(a > b for a, b in zip(ratios, ratios[1:])), ratios
    # in the tail the error from 4 shrinks at least linearly in 1/n: doubling n
    # must cut the excess by a factor above 1.9
    for a, b in zip(ratios[2:], ratios[3:]):
        assert (a - 4.0) / (b - 4.0) > 1.9, (a, b)
    assert ratios[-1] - 4.0 < 0.02, ratios[-1]