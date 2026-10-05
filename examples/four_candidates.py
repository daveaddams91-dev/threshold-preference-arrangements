"""Worked example: the maximum number of tolerant rankings for four candidates.

Run:  python examples/four_candidates.py

Prints, for a fixed generic configuration, the exact chamber count at several
tolerances and compares each with the Bezout ceiling, then verifies that every
reported intersection point really satisfies the two defining metric conditions.
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tolprefs import arrangement as A  # noqa: E402
from tolprefs import conics as C  # noqa: E402

POINTS = [(F(0), F(0)), (F(5), F(1)), (F(1), F(6)), (F(6), F(4))]


def main() -> int:
    n = len(POINTS)
    ceiling = A.semicircle_bound(n)
    classical = A.classical_region_count(n, 2)
    print(f"candidates: {[(int(p[0]), int(p[1])) for p in POINTS]}")
    print(f"n = {n}, so m = C(n,2) = {n * (n - 1) // 2} hyperbolas")
    print(f"Bezout ceiling        : {ceiling}")
    print(f"closed form 1+n^2(n-1)^2/2 : {ceiling}")
    print(f"classical (delta=0)   : {classical}")
    print()
    print(f"{'delta':>8} {'curves':>7} {'pairs x pts':>12} {'chambers':>9} "
          f"{'ceiling':>8}  saturated")
    print("-" * 60)
    for k in (1, 5, 10, 20, 50, 100, 200):
        delta = F(k, 100)
        stats = A.arrangement_stats(POINTS, delta)
        sat = "yes" if stats.region_count == ceiling else "no"
        print(f"{str(delta):>8} {stats.n_curves:>7} "
              f"{stats.n_pair_intersections:>12} {stats.region_count:>9} "
              f"{ceiling:>8}  {sat}")
    print()
    print("Interpretation: 4 candidates, small tolerance.")
    print(f"  classical single-line arrangements : {classical} chambers")
    print(f"  tolerant hyperbola arrangements   : {ceiling} chambers")
    print(f"  factor                            : {ceiling / classical:.2f}")
    print()

    # verify one pair exactly
    delta = F(1, 20)
    e1, e2 = (0, 1), (0, 2)
    c1 = C.gamma_conic(POINTS[e1[0]], POINTS[e1[1]], delta)
    c2 = C.gamma_conic(POINTS[e2[0]], POINTS[e2[1]], delta)
    pts = C.conic_intersections(c1, c2)
    print(f"intersection of H_{e1[0]}{e1[1]} and H_{e2[0]}{e2[1]} at "
          f"delta = {delta}:")
    print(f"  exact count: {C.conic_intersection_count(c1, c2)} "
          f"(Bezout ceiling for two conics: 4)")
    for p in pts:
        v = p.value()
        res = []
        for e in (e1, e2):
            d = lambda a: math.hypot(v[0] - a[0], v[1] - a[1])
            res.append(abs(abs(d(POINTS[e[0]]) - d(POINTS[e[1]]))
                         - float(delta)))
        print(f"  ({v[0]:.9f}, {v[1]:.9f})   "
              f"metric residuals {res[0]:.1e}, {res[1]:.1e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())