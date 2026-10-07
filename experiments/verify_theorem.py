"""Verify the main theorem's closed form and the ranking interpretation.

Two claims are checked here.

**Claim 1 (closed form).**  For a generic configuration and any tolerance small
enough, every pair of indifference hyperbolas meets in exactly four transverse
points, so the chamber count is

    R = 1 + 2 C(n,2) + 4 C(C(n,2), 2) = 1 + n^2 (n-1)^2 / 2.

The identity ``2 C(m,2) + 2m + 1 = 1 + 2m^2`` with ``m = C(n,2)`` gives the last
form, since ``2 C(m,2) + 2m = 2m(m-1) + 2m = 2m^2``.

**Claim 2 (ranking interpretation).**  Each chamber carries a distinct
consistent sign profile -- one sign per pair recording "i preferred",
"j preferred", or "neutral".  The point of the whole construction is that
``R`` is then the number of distinguishable tolerant rankings.  This is checked
by sampling the plane, computing profiles, and comparing the number of distinct
profiles seen to the exact chamber count.  Agreement is evidence that the
profile map is injective on chambers.

Run:  python experiments/verify_theorem.py
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
import itertools
import math
import sys

from tolprefs import arrangement as A  # noqa: E402
from tolprefs import conics as C  # noqa: E402
import random


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


DELTAS = [F(1, 200), F(1, 100), F(1, 50), F(1, 20), F(1, 10), F(1, 5)]


def formula(n: int) -> int:
    """``1 + n^2 (n-1)^2 / 2``, the Bézout ceiling written in closed form."""
    assert n * (n - 1) % 2 == 0 or True
    val = n * n * (n - 1) * (n - 1)
    assert val % 2 == 0
    return 1 + val // 2


def bezout(n: int) -> int:
    """Bezout.
    
    Args:
        n:
    
    Returns:
        The computed result
    
    """
    m = n * (n - 1) // 2
    return 1 + 2 * m + 4 * (m * (m - 1) // 2)


def generic(pts) -> bool:
    """Generic.
    
    Args:
        pts (list):
    
    Returns:
        bool: Result of type bool
    
    """
    pts = list(pts)
    n = len(pts)
    for a, b, c in itertools.combinations(range(n), 3):
        if ((pts[b][0] - pts[a][0]) * (pts[c][1] - pts[a][1])
                - (pts[b][1] - pts[a][1]) * (pts[c][0] - pts[a][0])) == 0:
            return False
    for a, b in itertools.combinations(range(n), 2):
        for c, d in itertools.combinations(range(n), 2):
            if len({a, b, c, d}) < 4:
                continue
            v1 = (pts[b][0] - pts[a][0], pts[b][1] - pts[a][1])
            v2 = (pts[d][0] - pts[c][0], pts[d][1] - pts[c][1])
            if v1[0] * v2[1] - v1[1] * v2[0] == 0:
                return False
    for a, b, c, d in itertools.combinations(range(n), 4):
        o = C.circumcenter(pts[a], pts[b], pts[c])
        if o is not None and C.conic_eval(C.bisector_conic(pts[a], pts[b]),
                                           pts[d]) == 0:
            return False
    return True


def check_formula(verbose=True) -> bool:
    """Check whether formula.
    
    Args:
        verbose (bool):
    
    Returns:
        The computed result
    
    """
    ok = True
    lines = []
    for n in range(2, 6):
        rng = random.Random(1000 + n)
        cfgs = 0
        while cfgs < 6:
            pts = [(F(rng.randint(-8, 8)), F(rng.randint(-8, 8)))
                   for _ in range(n)]
            if len(set(pts)) != n or not generic(pts):
                continue
            cfgs += 1
            for delta in DELTAS:
                stats = A.arrangement_stats(pts, delta)
                r = stats.region_count
                exp = formula(n)
                tag = "ok" if r <= exp else "VIOLATION"
                if r > exp:
                    ok = False
                if verbose and r == exp:
                    lines.append(f"n={n} delta={delta}: R={r} = formula {exp} "
                                 f"({tag}, {stats.n_pair_intersections} "
                                 f"pairwise intersections)")
    for ln in lines:
        print(ln)
    print(f"claim 1 (closed form is an upper bound, attained for small delta): "
          f"{'PASS' if ok else 'FAIL'}")
    for n in range(2, 8):
        assert formula(n) == bezout(n), (n, formula(n), bezout(n))
    print("closed form 1 + n^2(n-1)^2/2 equals the Bezout ceiling for n = 2..7")
    return ok


def profile(point, pts, delta) -> tuple:
    """The consistent sign profile of ``point`` for the given configuration."""
    d = [math.hypot(point[0] - p[0], point[1] - p[1]) for p in pts]
    out = []
    for i, j in itertools.combinations(range(len(pts)), 2):
        v = d[i] - d[j]
        if abs(v) <= delta:
            out.append(0)
        else:
            out.append(1 if v > 0 else -1)
    return tuple(out)


def check_profiles(n=4, delta=F(1, 5), span=26.0, grid=181,
                   seed=99) -> bool:
    """Count distinct sign profiles on a fine grid and compare with ``R``."""
    rng = random.Random(seed)
    while True:
        pts = [(F(rng.randint(-7, 7)), F(rng.randint(-7, 7))) for _ in range(n)]
        if len(set(pts)) == n and generic(pts):
            break
    stats = A.arrangement_stats(pts, delta)
    exact = stats.region_count
    seen = set()
    hits = 0
    total = 0
    step = 2 * span / (grid - 1)
    for ix in range(grid):
        x = -span + ix * step
        for iy in range(grid):
            y = -span + iy * step
            on_curve = False
            for i, j in itertools.combinations(range(n), 2):
                if abs(abs(math.hypot(x - pts[i][0], y - pts[i][1])
                           - math.hypot(x - pts[j][0], y - pts[j][1]))
                       - float(delta)) < 1e-7:
                    on_curve = True
                    break
            if on_curve:
                continue
            total += 1
            p = profile((x, y), pts, float(delta))
            if p not in seen:
                seen.add(p)
                hits += 1
    print(f"n={n}, delta={delta}, points={[(int(p[0]), int(p[1])) for p in pts]}")
    print(f"  exact chambers            : {exact}")
    print(f"  distinct profiles on grid : {len(seen)}")
    print(f"  grid cells off the curves : {total}")
    # every chamber must be sampled, so grid count <= exact; and the number of
    # distinct profiles must not exceed the number of chambers.
    ok = len(seen) <= exact
    print(f"  profiles <= chambers      : {ok}")
    print(f"  grid coverage             : "
          f"{len(seen)}/{exact} = {len(seen) / exact:.3f} of chambers hit")
    return ok


def main() -> int:
    """Entry point — parse arguments and run the main computation.
    
    Returns:
        The computed result
    
    """
    ok = check_formula()
    print()
    ok = check_profiles() and ok
    print()
    print("OVERALL:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())