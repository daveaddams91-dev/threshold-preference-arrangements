"""Validate the exact conic-intersection counter against an independent method.

Our counter is exact (pencil elimination + Sturm isolation + exact gcd).  To
check it we use two independent references:

1. **Numerical root finding on the metric equations.**  Sample the *defining*
   conditions ``|d_i - d_j| = delta`` directly -- not the conic equations --
   and locate solutions with a fine multistart Newton iteration on
   ``f(X) = (|X-Q_i| - |X-Q_j|)^2 - delta^2 = 0`` for the two curve pairs
   simultaneously.  This shares no code with the conic algebra, so agreement is
   meaningful.  It is a *sampling* method: it can miss solutions, so it certifies
   a lower bound.

2. **Direct verification of every point we report.**  Each returned point must
   satisfy ``|d_i - d_j| = delta`` for both curve pairs.  This certifies an
   upper bound and catches spurious or mis-coordinated points outright.

Together: every reported point is real (check 2) and no point is missed relative
to the independent search (check 1).

Run:  python experiments/validate_intersections.py
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
import itertools
import math
import sys

from tolprefs import conics as C
import random


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))



def _f(p, a, b, delta):
    """Squared metric condition, smooth away from the foci."""
    da = math.hypot(p[0] - a[0], p[1] - a[1])
    db = math.hypot(p[0] - b[0], p[1] - b[1])
    return (abs(da - db) - delta) * (da + db)  # = da^2 - db^2 - delta(da+db)


def numeric_count(pts, e1, e2, delta, seed=0, starts=400, span=14.0) -> tuple:
    """Distinct real solutions of both metric conditions, by multistart Newton."""
    A, B = pts[e1[0]], pts[e1[1]]
    Cc, D = pts[e2[0]], pts[e2[1]]
    rng = random.Random(seed)
    h = 1e-7

    def g(x, y) -> tuple:
        """G.
        
        Args:
            x:
            y:
        
        Returns:
            tuple: Result of type tuple
        
        """
        return (_f((x, y), A, B, delta), _f((x, y), Cc, D, delta))

    def jac(x, y) -> tuple:
        """Jac.
        
        Args:
            x:
            y:
        
        Returns:
            tuple: Result of type tuple
        
        """
        j11 = (g(x + h, y)[0] - g(x - h, y)[0]) / (2 * h)
        j12 = (g(x, y + h)[0] - g(x, y - h)[0]) / (2 * h)
        j21 = (g(x + h, y)[1] - g(x - h, y)[1]) / (2 * h)
        j22 = (g(x, y + h)[1] - g(x, y - h)[1]) / (2 * h)
        return j11, j12, j21, j22

    sols = []
    for _ in range(starts):
        x, y = rng.uniform(-span, span), rng.uniform(-span, span)
        ok = True
        for _ in range(120):
            f1, f2 = g(x, y)
            if abs(f1) < 1e-12 and abs(f2) < 1e-12:
                break
            j11, j12, j21, j22 = jac(x, y)
            det = j11 * j22 - j12 * j21
            if abs(det) < 1e-13:
                ok = False
                break
            dx = (-f1 * j22 + f2 * j12) / det
            dy = (-j11 * f2 + j21 * f1) / det
            step = math.hypot(dx, dy)
            if step > span:
                dx, dy = dx * span / step, dy * span / step
            x += dx
            y += dy
            if step < 1e-14:
                break
        if not ok:
            continue
        if abs(g(x, y)[0]) < 1e-9 and abs(g(x, y)[1]) < 1e-9:
            if all(math.hypot(x - a, y - b) > 1e-4 for a, b in sols):
                sols.append((x, y))
    return sols


def main(trials=40, seed=20240):
    """Entry point — parse arguments and run the main computation.
    
    Args:
        trials (int):
        seed (int):
    
    Returns:
        The computed result
    
    """
    random.seed(seed)
    pairs = 0
    spurious = 0
    missed = 0
    worst = 0.0
    for _ in range(trials):
        n = random.choice([3, 3, 4])
        pts = [(F(random.randint(-6, 6)), F(random.randint(-6, 6)))
               for _ in range(n)]
        if len(set(pts)) < n:
            continue
        delta = F(random.choice([1, 2, 3, 5, 10])) / F(10)
        cs = {}
        for i, j in itertools.combinations(range(n), 2):
            c = C.gamma_conic(pts[i], pts[j], delta)
            if c is not None:
                cs[(i, j)] = c
        for e1, e2 in itertools.combinations(cs, 2):
            pairs += 1
            c1, c2 = cs[e1], cs[e2]
            ours = []
            for p in C.conic_intersections(c1, c2):
                xv, yv = p.value()
                d = lambda a: math.hypot(xv - a[0], yv - a[1])
                r1 = abs(abs(d(pts[e1[0]]) - d(pts[e1[1]])) - float(delta))
                r2 = abs(abs(d(pts[e2[0]]) - d(pts[e2[1]])) - float(delta))
                if max(r1, r2) > 1e-6:
                    spurious += 1
                    print("  SPURIOUS (%.8f,%.8f) residuals %.2e %.2e"
                          % (xv, yv, r1, r2))
                else:
                    worst = max(worst, r1, r2)
                ours.append((xv, yv))
            ref = numeric_count(pts, e1, e2, float(delta),
                                seed=pairs, starts=800, span=16.0)
            # every independently found point must be one of ours
            for rx, ry in ref:
                if all(math.hypot(rx - ox, ry - oy) > 1e-3 for ox, oy in ours):
                    missed += 1
                    print("  MISSED (%.8f,%.8f) not among our %d points"
                          % (rx, ry, len(ours)))
                    print("     pts", [(int(p[0]), int(p[1])) for p in pts],
                          "delta", delta, "pair", e1, e2)
    print(f"pairs tested            : {pairs}")
    print(f"spurious points         : {spurious}")
    print(f"points missed           : {missed}")
    print(f"worst residual (valid)  : {worst:.3e}")
    return spurious + missed


if __name__ == "__main__":
    sys.exit(1 if main() else 0)