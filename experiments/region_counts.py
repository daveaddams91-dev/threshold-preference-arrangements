"""Main experiment: exact region counts for the indifference arrangement.

For ``n`` points and tolerance ``delta`` the arrangement ``A(Q, delta)`` of the
``C(n,2)`` indifference hyperbolas has, in general position,

    R(Q, delta) = 1 + 2 C(n,2) + I(Q, delta)

chambers, where ``I`` is the number of intersection points between distinct
hyperbolas (see :mod:`tolprefs.arrangement` for the derivation and for why
general position is required).  Bézout gives ``I <= 4 C(C(n,2), 2)``, hence the
ceiling

    R <= 1 + 2 C(n,2) + 4 C(C(n,2), 2)                                        (*)

and the question is how close to (*) the truth gets.  This script computes ``R``
exactly for a large collection of configurations and reports, per ``n``, the
largest value found together with the configuration achieving it.

Outputs
-------
``results/region_counts.csv``    one row per configuration
``results/small_cases.txt``      the exactly-settled small cases
``results/summary.txt``          maxima and the Bezout ceiling, per n

Run:  python experiments/region_counts.py [--trials N] [--seed S]
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
import argparse
import itertools
import sys
import time

from tolprefs import arrangement as A  # noqa: E402
from tolprefs import conics as C  # noqa: E402
import random


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


DELTAS = [F(1, 20), F(1, 10), F(1, 5), F(1, 4), F(1, 3), F(1, 2), F(3, 4),
          F(1), F(5, 4), F(3, 2), F(2), F(5, 2), F(3), F(4), F(5), F(6),
          F(8), F(10)]


def general_position(pts) -> bool:
    """The configuration is generic in the sense required by formula (2).

    Conditions: no three points collinear, no two connecting segments parallel,
    no two segments of equal length, no four points concyclic.  Each of these can
    force a tangency or a triple point, which formula (2) does not account for.
    """
    pts = list(pts)
    n = len(pts)
    for a, b, c in itertools.combinations(range(n), 3):
        cross = ((pts[b][0] - pts[a][0]) * (pts[c][1] - pts[a][1])
                 - (pts[b][1] - pts[a][1]) * (pts[c][0] - pts[a][0]))
        if cross == 0:
            return False
    for a, b in itertools.combinations(range(n), 2):
        for c, d in itertools.combinations(range(n), 2):
            if len({a, b, c, d}) < 4:
                continue
            v1 = (pts[b][0] - pts[a][0], pts[b][1] - pts[a][1])
            v2 = (pts[d][0] - pts[c][0], pts[d][1] - pts[c][1])
            if v1[0] * v2[1] - v1[1] * v2[0] == 0:
                return False
    for a, b in itertools.combinations(range(n), 2):
        for c, d in itertools.combinations(range(n), 2):
            if len({a, b, c, d}) < 4:
                continue
            if C.dist2(pts[a], pts[b]) == C.dist2(pts[c], pts[d]):
                return False
    for a, b, c, d in itertools.combinations(range(n), 4):
        if (C.circumcenter(pts[a], pts[b], pts[c]) is not None
                and C.conic_eval(C.bisector_conic(pts[a], pts[b]), pts[d]) == 0):
            return False
    return True


def random_config(rng, n, span=9):
    """A random integer configuration, retried until generic."""
    for _ in range(400):
        pts = [(F(rng.randint(-span, span)), F(rng.randint(-span, span)))
               for _ in range(n)]
        if len(set(pts)) == n and general_position(pts):
            return pts
    raise RuntimeError(f"no generic configuration found for n={n}")


def point_curve_intersections(pts, delta) -> tuple:
    """Sum of pairwise intersection counts, and the pairs that fall short of 4.

    Returns ``(total, sub_pairs)`` where ``sub_pairs`` lists
    ``((i,j), (k,l), count)`` for every curve pair meeting fewer than 4 times.
    """
    conics = {}
    for i, j in itertools.combinations(range(len(pts)), 2):
        c = C.gamma_conic(pts[i], pts[j], delta)
        if c is not None:
            conics[(i, j)] = c
    total = 0
    sub = []
    for (e1, e2) in itertools.combinations(sorted(conics), 2):
        k = C.conic_intersection_count(conics[e1], conics[e2])
        total += k
        if k < 4:
            sub.append((e1, e2, k))
    return total, sub


def main() -> int:
    """Entry point — parse arguments and run the main computation.
    
    Returns:
        int: Result of type int
    
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=60,
                    help="random configurations per n")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    outdir = ROOT / "results"
    outdir.mkdir(exist_ok=True)
    csv = outdir / "region_counts.csv"
    lines = ["n,delta,delta_float,n_curves,pair_intersections,"
             "distinct_points,regions,bezout_ceiling,n_degenerate_pairs,"
             "min_separation,points,generic"]
    summary_lines = []

    for n in range(2, 7):
        m = n * (n - 1) // 2
        ceiling = 1 + 2 * m + 4 * (m * (m - 1) // 2)
        best = None
        best_desc = None
        t0 = time.time()
        trials_done = 0
        for t in range(args.trials):
            pts = random_config(rng, n)
            delta = rng.choice(DELTAS)
            total, sub = point_curve_intersections(pts, delta)
            stats = A.arrangement_stats(pts, delta)
            regions = stats.region_count
            trials_done += 1
            lines.append(
                f"{n},{delta},{float(delta)},{stats.n_curves},{total},"
                f"{stats.n_distinct_points},{regions},{ceiling},"
                f"{len(sub)},{stats.min_separation:.3e},"
                f"\"{';'.join(f'{int(p[0])},{int(p[1])}' for p in pts)}\","
                f"{len(stats.degenerate_pairs) == 0}")
            if best is None or regions > best:
                best = regions
                best_desc = (pts, delta, total, len(sub))
        dt = time.time() - t0
        frac = (best / ceiling) if ceiling else 0.0
        summary_lines.append(
            f"n={n}: curves={m}, Bézout ceiling={ceiling}, "
            f"max regions found={best}, attained at delta={best_desc[1]} "
            f"with {best_desc[2]} pairwise intersections and "
            f"{best_desc[3]} pairs below 4; "
            f"ratio to ceiling={frac:.4f}; {trials_done} trials in {dt:.1f}s")
        print(summary_lines[-1], flush=True)

    csv.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (outdir / "summary.txt").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8")

    # ---- the exactly-settled small cases -------------------------------
    small = []
    small.append("n=2 (two candidates)")
    small.append("-" * 60)
    p2 = [(F(0), F(0)), (F(4), F(3))]
    for delta in (F(1, 10), F(1, 2), F(4), F(5), F(6)):
        s = A.arrangement_stats(p2, delta)
        small.append(f"  delta={delta}: curves={s.n_curves}, "
                     f"regions={s.region_count}, ceiling={A.semicircle_bound(2)}")
    small.append("")
    small.append("n=3 (three candidates)")
    small.append("-" * 60)
    p3 = [(F(0), F(0)), (F(5), F(1)), (F(1), F(6))]
    for delta in (F(1, 10), F(1, 2), F(1), F(3), F(5), F(7), F(8)):
        s = A.arrangement_stats(p3, delta)
        small.append(f"  delta={delta}: curves={s.n_curves}, "
                     f"pair intersections={s.n_pair_intersections}, "
                     f"regions={s.region_count}, ceiling={A.semicircle_bound(3)}")
    small.append("")
    small.append("classical delta=0 baseline (Good-Tideman, plane):")
    for n in range(2, 8):
        small.append(f"  n={n}: {A.classical_region_count(n, 2)} regions")
    (outdir / "small_cases.txt").write_text("\n".join(small) + "\n",
                                           encoding="utf-8")
    print()
    print("\n".join(small))
    return 0


if __name__ == "__main__":
    sys.exit(main())