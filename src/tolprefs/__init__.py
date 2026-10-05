"""Tolerance-based ("delta-sophisticated") distance orderings of point sets.

Given ``n`` candidate points ``Q_1, ..., Q_n`` in the plane and a tolerance
``delta > 0``, a query point ``P`` prefers ``i`` over ``j`` when

    d(P, Q_i) + delta < d(P, Q_j),

is indifferent when ``|d(P,Q_i) - d(P,Q_j)| <= delta``, and otherwise prefers
``j``.  As ``P`` ranges over the plane, the induced *sign profile* varies; the
number of distinct profiles is the number of chambers of the arrangement of
indifference hyperbolas, computed exactly in :mod:`tolprefs.arrangement`.

At ``delta = 0`` this reduces to the classical Good--Tideman problem, whose
generic region count is a sum of unsigned Stirling numbers of the first kind.
"""

from .conics import (
    AlgebraicPoint,
    Conic,
    F,
    Point,
    as_point,
    circumcenter,
    conic_eval,
    conic_intersection_count,
    conic_intersections,
    dist2,
    gamma_conic,
    gamma_nonempty,
    perpendicular_bisector,
)
from .arrangement import (
    ArrangementStats,
    CurvePair,
    arrangement_stats,
    classical_region_count,
    region_count,
    semicircle_bound,
)

__all__ = [
    "F",
    "Point",
    "Conic",
    "AlgebraicPoint",
    "as_point",
    "dist2",
    "gamma_conic",
    "gamma_nonempty",
    "conic_eval",
    "conic_intersections",
    "conic_intersection_count",
    "perpendicular_bisector",
    "circumcenter",
    "CurvePair",
    "ArrangementStats",
    "arrangement_stats",
    "region_count",
    "classical_region_count",
    "semicircle_bound",
    "__version__",
]

__version__ = "0.1.0"