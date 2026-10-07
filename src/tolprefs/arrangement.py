"""Region counts for indifference arrangements, and the ordering interpretation.

The central object
------------------
For ``Q = {Q_0, ..., Q_{n-1}}`` and a tolerance ``delta > 0`` let

    A(Q, delta) = { H_ij : 0 <= i < j < n }

be the arrangement of the ``C(n,2)`` indifference hyperbolas, and let

    R(Q, delta) = number of chambers (connected components of the complement).

Each chamber carries a *consistent tolerant preference profile*: a sign in
``{+,-,0}`` for each pair ``{i,j}``, recording whether ``Q_i`` is preferred,
dispreferred, or neutral.  Thus ``R`` is exactly the number of distinguishable
"delta-sophisticated" rankings of the ``n`` candidates.

Why the region count is what we compute
--------------------------------------
Each ``H_ij`` is a proper embedded copy of the real line (its two branches are
two disjoint proper arcs, each running off to infinity at both ends).  If the
arrangement is in *general position* -- every pairwise intersection is
transverse, no three curves are concurrent, and no two curves share a vertical
fibre -- then the classical incremental argument for arrangements of proper
arcs gives

    R = 1 + (#arcs) + (#intersection points)                     (1)

Indeed, adding an arc meeting the previous arrangement in ``k >= 1`` transverse
points subdivides it into ``k+1`` pieces, each splitting one region; an arc
meeting nothing adds one region.  Summing over the arcs, and observing that
each pairwise intersection is counted exactly once (by the curve added later),
gives (1).  Each hyperbola contributes two arcs, so ``#arcs = 2 C(n,2)``, and
with ``I(Q, delta)`` denoting the number of intersection points,

    R(Q, delta) = 1 + 2 C(n,2) + I(Q, delta)                      (2)

in general position.  Degeneracies lower the count, and
:func:`region_count` reports ``I`` together with a degeneracy flag so that
formula (2) is never applied blindly.

Note the classical limit: at ``delta = 0`` each ``H_ij`` degenerates to a line
(the perpendicular bisector) and the arrangement is that of Good and Tideman,
whose generic region count is ``|s(n,n)| + |s(n,n-1)| + ... + |s(n,n-d)|``
(``s`` = signed Stirling number of the first kind).  Our arrangement is the
curved, "second-level" analogue that Zaslavsky identified as an open problem.
"""

from __future__ import annotations

from . import conics as C
from dataclasses import dataclass, field
from fractions import Fraction as F
import itertools



__all__ = [
    "CurvePair",
    "ArrangementStats",
    "arrangement_stats",
    "region_count",
    "classical_region_count",
    "semicircle_bound",
]


@dataclass(frozen=True)
class CurvePair:
    """A pair of indices labelling one indifference curve."""

    i: int
    j: int


@dataclass
class ArrangementStats:
    """Exact combinatorial data for one arrangement."""

    n: int
    delta: F
    n_curves: int                 # number of non-empty indifference curves
    n_arcs: int                   # 2 per curve (both branches)
    n_pair_intersections: int     # sum over curve pairs of the exact counts
    n_distinct_points: int        # after removing coincidences
    min_separation: float         # smallest gap between distinct points
    degenerate_pairs: list = field(default_factory=list)
    empty_pairs: list = field(default_factory=list)

    @property
    def region_count_general_position(self) -> int:
        """Formula (2); valid only when ``degenerate_pairs`` is empty."""
        return 1 + self.n_arcs + self.n_distinct_points

    @property
    def region_count(self) -> int:
        """Chamber count; equals formula (2) in general position."""
        return self.region_count_general_position


def arrangement_stats(pts, delta: F) -> ArrangementStats:
    """Exact region statistics for the arrangement of indifference hyperbolas.

    Parameters
    ----------
    pts:
        Sequence of distinct points with rational coordinates.
    delta:
        Positive rational tolerance.
    """
    pts = [C.as_point(p) for p in pts]
    n = len(pts)
    delta = F(delta)

    pairs: list[CurvePair] = []
    conics: list = []
    empty: list = []
    for i, j in itertools.combinations(range(n), 2):
        c = C.gamma_conic(pts[i], pts[j], delta)
        if c is None:
            empty.append(CurvePair(i, j))
        else:
            pairs.append(CurvePair(i, j))
            conics.append(c)

    total = 0
    degenerate: list = []
    coords: list[tuple[float, float]] = []
    for (e1, c1), (e2, c2) in itertools.combinations(list(zip(pairs, conics)), 2):
        k = C.conic_intersection_count(c1, c2)
        total += k
        if k != 4:  # 4 is the Bezout ceiling; fewer means a real degeneracy
            degenerate.append((e1, e2, k))
        for p in C.conic_intersections(c1, c2):
            coords.append(p.value())

    # Remove coincident points (three or more curves through one point).
    coords.sort()
    distinct: list[tuple[float, float]] = []
    min_sep = float("inf")
    for x, y in coords:
        if distinct:
            d = ((x - distinct[-1][0]) ** 2 + (y - distinct[-1][1]) ** 2) ** 0.5
            min_sep = min(min_sep, d)
            if d < 1e-7:
                continue
        distinct.append((x, y))

    return ArrangementStats(
        n=n,
        delta=delta,
        n_curves=len(pairs),
        n_arcs=2 * len(pairs),
        n_pair_intersections=total,
        n_distinct_points=len(distinct),
        min_separation=(0.0 if min_sep == float("inf") else min_sep),
        degenerate_pairs=degenerate,
        empty_pairs=empty,
    )


def region_count(pts, delta: F) -> int:
    """Number of chambers of ``A(Q, delta)`` in general position."""
    return arrangement_stats(pts, delta).region_count


def classical_region_count(n: int, d: int) -> int:
    """Good--Tideman generic region count for ``n`` points in ``R^d``.

    ``sum_{i=0}^{d} |s(n, n-i)|`` with ``s`` the signed Stirling number of the
    first kind.  Provided so that the ``delta -> 0`` limit of our arrangements
    can be compared against the classical baseline in one place.
    """
    if d < 0:
        raise ValueError("d must be non-negative")
    if n < 0:
        raise ValueError("n must be non-negative")
    total = 0
    # s(n, k) via the recurrence s(n,k) = s(n-1,k-1) - (n-1) s(n-1,k)
    row = [1]  # s(0, 0)
    for m in range(1, n + 1):
        new = [0] * (m + 1)
        for k in range(1, m + 1):
            new[k] = (row[k - 1] if k - 1 < len(row) else 0) - (m - 1) * (
                row[k] if k < len(row) else 0
            )
        new[0] = 0
        row = new
    for i in range(0, min(d, n) + 1):
        k = n - i
        total += abs(row[k] if k < len(row) else 0)
    return total


def semicircle_bound(n: int) -> int:
    """Bezout ceiling on the number of intersection points.

    Two distinct non-degenerate conics meet in at most four real points, so

        I(Q, delta) <= 4 * C(C(n,2), 2),

    and formula (2) gives ``1 + 2 C(n,2) + 4 C(C(n,2), 2)``.  Pairs of
    *coincident* curves (``delta = |Q_iQ_j|``) contribute no intersections, so
    the bound is only attained for a configuration with no such pair.
    """
    m = n * (n - 1) // 2
    return 1 + 2 * m + 4 * (m * (m - 1) // 2)