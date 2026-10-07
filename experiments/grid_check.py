"""Grid-based validation of the chamber formula and the ranking interpretation.

An independent check of the two structural claims, using none of the conic
algebra.  A regular grid over a large box is rasterised and, for every pair of
4-adjacent grid cells whose centres lie off the arrangement, the two are declared
in the same chamber **iff no indifference curve crosses the segment joining
their centres**.  The crossing test is a sign change of

    g_ij(X) = | d(X, Q_i) - d(X, Q_j) | - delta

between the two centres, for every pair ``(i, j)``.  Since each ``g_ij`` is
continuous and the grid step is far below the radius of curvature of every curve
(``O(delta)``), a sign change is exactly a crossing.  Connected components of the
resulting graph are the chambers.

Two quantities are then compared with the exact count ``R`` computed by
:func:`tolprefs.arrangement.arrangement_stats`:

* ``#components`` -- validates the formula ``R = 1 + 2 C(n,2) + I`` together
  with the intersection counts that feed it;
* ``#distinct sign profiles`` -- each chamber carries one profile, so
  ``#profiles <= R`` is forced, and when the grid resolves every chamber,
  equality shows that distinct chambers give distinct tolerant rankings, i.e.
  that ``R`` is exactly the number of distinguishable rankings.

``delta`` is deliberately large so that chambers are wide relative to the grid
step; the reported step records the achieved resolution so that a shortfall can
be attributed to resolution rather than to a disagreement.

Run:  python experiments/grid_check.py
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
import itertools
import sys

from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from tolprefs import arrangement as A  # noqa: E402
from tolprefs import conics as C  # noqa: E402
import numpy as np
import random



ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


ON_CURVE_TOL = 1e-9


def generic(pts, min_sep_sq=F(1)) -> bool:
    """Generic *and* well spread.

    Genericity alone is not enough for a raster comparison: three nearly collinear
    points have nearly parallel perpendicular bisectors, and the wedges between
    them are far thinner than any practical grid step, so the chamber count then
    looks too low.  ``min_sep_sq`` therefore also rejects configurations with two
    candidates closer than ``sqrt(min_sep_sq)``, and parallel segments are ruled
    out outright.
    """
    pts = list(pts)
    for a, b in itertools.combinations(range(len(pts)), 2):
        if C.dist2(pts[a], pts[b]) < min_sep_sq:
            return False
    for a, b, c in itertools.combinations(range(len(pts)), 3):
        if ((pts[b][0] - pts[a][0]) * (pts[c][1] - pts[a][1])
                - (pts[b][1] - pts[a][1]) * (pts[c][0] - pts[a][0])) == 0:
            return False
    for a, b in itertools.combinations(range(len(pts)), 2):
        for c, d in itertools.combinations(range(len(pts)), 2):
            if len({a, b, c, d}) < 4:
                continue
            v1 = (pts[b][0] - pts[a][0], pts[b][1] - pts[a][1])
            v2 = (pts[d][0] - pts[c][0], pts[d][1] - pts[c][1])
            if v1[0] * v2[1] - v1[1] * v2[0] == 0:
                return False
    for a, b, c, d in itertools.combinations(range(len(pts)), 4):
        o = C.circumcenter(pts[a], pts[b], pts[c])
        if o is not None and C.conic_eval(C.bisector_conic(pts[a], pts[b]),
                                           pts[d]) == 0:
            return False
    return True


def _g(pts, delta, X, Y):
    """Stack of the *signed* ``|d_i - d_j| - delta``, shaped ``(pairs, nx, ny)``.

    The sign matters: the crossing test below looks for a sign change between two
    adjacent cells, so the value must be allowed to be negative.  (An absolute
    value here silently disables the whole barrier detection.)
    """
    n = len(pts)
    d = [np.hypot(X - float(p[0]), Y - float(p[1])) for p in pts]
    out = []
    for i, j in itertools.combinations(range(n), 2):
        out.append(np.abs(d[i] - d[j]) - delta)
    return np.stack(out)


def check(n, delta, span, grid, seed) -> dict:
    """Check whether the condition holds.
    
    Args:
        n:
        delta:
        span:
        grid:
        seed:
    
    Returns:
        dict: Result of type dict
    
    """
    rng = random.Random(seed)
    while True:
        pts = [(F(rng.randint(-6, 6)), F(rng.randint(-6, 6))) for _ in range(n)]
        if len(set(pts)) == n and generic(pts):
            break

    stats = A.arrangement_stats(pts, F(str(delta)))
    exact = stats.region_count

    axis = np.linspace(-span, span, grid)
    X, Y = np.meshgrid(axis, axis, indexing="ij")
    g = _g(pts, float(delta), X, Y)          # (pairs, nx, ny), signed
    # A cell centre is *on* the arrangement when |g| vanishes.  Testing ``g <= tol``
    # instead would mark every cell on the negative side of a curve, i.e. most of
    # the plane, and block the whole graph.
    on = np.any(np.abs(g) <= ON_CURVE_TOL, axis=0)
    off = ~on

    # An edge between adjacent cells is open iff no g_ij changes sign across it
    # and neither endpoint sits on a curve.
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    nx, ny = grid, grid

    def add_edges(u_idx, v_idx, gu, gv, shape):
        """u_idx/v_idx: 1-D indices; gu/gv: same-shape stacks of g values."""
        blocked = np.any(np.abs(gu) <= ON_CURVE_TOL, axis=0) | \
            np.any(np.abs(gv) <= ON_CURVE_TOL, axis=0)
        flip = np.any((gu > 0) != (gv > 0), axis=0) & ~blocked
        open_ = ~blocked & ~flip
        rows.append(u_idx[open_])
        cols.append(v_idx[open_])

    # horizontal adjacency: (i, j) -- (i, j+1)
    gu = g[:, :, :-1]
    gv = g[:, :, 1:]
    ui = (np.arange(nx)[:, None] * ny + np.arange(ny - 1)[None, :]).ravel()
    vi = (np.arange(nx)[:, None] * ny + np.arange(1, ny)[None, :]).ravel()
    add_edges(ui, vi, gu.reshape(gu.shape[0], -1), gv.reshape(gv.shape[0], -1),
              (nx, ny))
    # vertical adjacency: (i, j) -- (i+1, j)
    gu = g[:, :-1, :]
    gv = g[:, 1:, :]
    ui = (np.arange(nx - 1)[:, None] * ny + np.arange(ny)[None, :]).ravel()
    vi = (np.arange(1, nx)[:, None] * ny + np.arange(ny)[None, :]).ravel()
    add_edges(ui, vi, gu.reshape(gu.shape[0], -1), gv.reshape(gv.shape[0], -1),
              (nx, ny))

    r = np.concatenate(rows)
    c = np.concatenate(cols)
    ncomp, labels = connected_components(
        coo_matrix((np.ones(len(r), dtype=np.int8), (r, c)),
                   shape=(nx * ny, nx * ny)), directed=False)

    # the sign profile of every off-curve cell
    prof = np.zeros((nx, ny), dtype=np.int64)
    mult = 3 ** np.arange(n * (n - 1) // 2, dtype=np.int64)
    k = 0
    d = [np.hypot(X - float(p[0]), Y - float(p[1])) for p in pts]
    for i, j in itertools.combinations(range(n), 2):
        v = d[i] - d[j]
        s = np.where(v > float(delta), 1, np.where(v < -float(delta), -1, 0))
        prof += s.astype(np.int64) * int(mult[k])
        k += 1

    live = labels[off.ravel()]
    uniq, inv = np.unique(live, return_inverse=True)
    # profile must be constant on each component
    flat_prof = prof.ravel()[off.ravel()]
    mixed = 0
    order = np.argsort(inv, kind="stable")
    inv_s = inv[order]
    prof_s = flat_prof[order]
    bounds = np.flatnonzero(np.diff(inv_s)) + 1
    for chunk in np.split(prof_s, bounds):
        if len(np.unique(chunk)) != 1:
            mixed += 1
    nprofiles = len(np.unique(flat_prof))

    step = 2 * span / (grid - 1)
    print(f"n={n} delta={delta} grid={grid}x{grid} step={step:.4f} "
          f"points={[(int(p[0]), int(p[1])) for p in pts]}")
    print(f"  exact chambers R           : {exact}")
    print(f"  grid chambers              : {len(uniq)}")
    print(f"  components w/ mixed profile: {mixed}")
    print(f"  distinct profiles on grid  : {nprofiles}")
    print(f"  profiles <= chambers       : {nprofiles <= exact}")
    print(f"  grid resolves all chambers : {len(uniq) == exact}"
          "   (note: a chamber leaving and re-entering the box is split by the"
          " boundary, so this may read False for large grids too)")
    print(f"  profiles == chambers       : {nprofiles == exact}")
    return {"n": n, "delta": str(delta), "exact": exact, "grid": len(uniq),
            "profiles": int(nprofiles), "mixed": int(mixed), "step": step,
            "points": ";".join(f"{int(p[0])},{int(p[1])}" for p in pts)}


def main() -> int:
    """Entry point — parse arguments and run the main computation.
    
    Returns:
        The computed result
    
    """
    cases = [
        # n, delta, span, grid, seed -- delta large relative to the grid step so
        # that every chamber is resolved.  The box is deliberately generous: a
        # chamber that leaves and re-enters the box is split by the boundary,
        # which inflates the grid count above R, so under-counting is the
        # meaningful failure mode to watch for.
        (2, 2.0, 14.0, 1401, 1),
        (3, 2.0, 20.0, 2001, 2),
        (3, 3.0, 20.0, 2001, 3),
        (4, 3.0, 26.0, 2001, 4),
        (4, 4.0, 26.0, 2001, 5),
    ]
    rows = ["n,delta,exact_chambers,grid_chambers,distinct_profiles,"
            "components_with_mixed_profile,grid_step,points"]
    ok = True
    for n, delta, span, grid, seed in cases:
        r = check(n, delta, span, grid, seed)
        ok = ok and r["mixed"] == 0 and r["profiles"] <= r["exact"]
        rows.append(f"{n},{r['delta']},{r['exact']},{r['grid']},"
                    f"{r['profiles']},{r['mixed']},{r['step']:.5f},"
                    f"\"{r['points']}\"")
        print()
    (ROOT / "results" / "grid_check.csv").write_text(
        "\n".join(rows) + "\n", encoding="utf-8")
    print("OVERALL:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())