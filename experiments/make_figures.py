"""Figures for the paper and README.

All panels are drawn from the exact machinery in :mod:`tolprefs`; nothing is
plotted from hard-coded coordinates.  Each figure records the configuration and
tolerance it was produced from in its title, and the underlying numbers are
recomputed here rather than copied from the results files, so a figure cannot
silently disagree with the code.

Run:  python experiments/make_figures.py
"""
from __future__ import annotations

import sys
from fractions import Fraction as F
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tolprefs import arrangement as A  # noqa: E402
from tolprefs import conics as C  # noqa: E402

FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)


def _grid(pts, delta, span, n=900) -> tuple:
    """``Z(x, y) = min over pairs of | |d_i - d_j| - delta |`` on a grid.

    Zero level sets of ``Z`` are exactly the indifference curves.  Coordinates
    are converted to ``float`` up front: NumPy defers ``ndarray - Fraction`` to
    ``Fraction.__rsub__``, which collapses the array to a scalar, so mixing exact
    and floating arithmetic here silently destroys the grid.
    """
    axis = np.linspace(-span, span, n)
    X, Y = np.meshgrid(axis, axis, indexing="ij")
    px = [float(p[0]) for p in pts]
    py = [float(p[1]) for p in pts]
    dists = [np.hypot(X - px[i], Y - py[i]) for i in range(len(pts))]
    Z = None
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            band = np.abs(np.abs(dists[i] - dists[j]) - delta)
            Z = band if Z is None else np.minimum(Z, band)
    return axis, X, Y, Z


def fig_arrangement(pts, delta, title, fname, span=13.0, n=900, ncurves=None):
    """Fig arrangement.
    
    Args:
        pts:
        delta:
        title:
        fname:
        span (float):
        n (int):
        ncurves:
    
    """
    axis, X, Y, Z = _grid(pts, float(delta), span, n)
    stats = A.arrangement_stats(pts, F(str(delta)))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4))

    ax = axes[0]
    lev = np.array([0.0])
    cs = ax.contour(X, Y, Z, levels=lev, colors="k", linewidths=1.1)
    ax.contourf(X, Y, np.clip(Z, 0, 1), levels=[0, 1e-12, 1],
                colors=["white", "#f3f3f3"])
    ax.plot([float(p[0]) for p in pts], [float(p[1]) for p in pts],
            "o", color="#c1121f", markersize=7, zorder=5)
    for k, p in enumerate(pts):
        ax.annotate(f"$Q_{k}$", (float(p[0]), float(p[1])),
                    textcoords="offset points", xytext=(6, 6), fontsize=10)
    ax.set_title(f"arrangement $A(Q,\\delta)$, $\\delta={delta}$")
    ax.set_aspect("equal")
    ax.set_xlim(-span, span)
    ax.set_ylim(-span, span)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")

    ax = axes[1]
    for k in range(stats.n_curves):
        col = plt.cm.tab10(k % 10)
        pair = stats_pair(pts, k)
        draw_conic(ax, C.gamma_conic(pts[pair[0]], pts[pair[1]],
                                     F(str(delta))), span, col,
                   f"$H_{{{pair[0]}{pair[1]}}}$")
    ax.plot([float(p[0]) for p in pts], [float(p[1]) for p in pts],
            "k*", markersize=11, zorder=5)
    ax.set_title(f"the {stats.n_curves} indifference hyperbolas")
    ax.set_aspect("equal")
    ax.set_xlim(-span, span)
    ax.set_ylim(-span, span)

    fig.suptitle(f"{title}\nchambers $R={stats.region_count}$, "
                 f"pairwise intersections $I={stats.n_pair_intersections}$",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / fname, dpi=150)
    plt.close(fig)
    print(f"wrote {fname}  (R={stats.region_count}, "
          f"I={stats.n_pair_intersections}, curves={stats.n_curves})")


def stats_pair(pts, k):
    """Stats pair.
    
    Args:
        pts (list):
        k:
    
    Returns:
        The computed result
    
    """
    import itertools
    return list(itertools.combinations(range(len(pts)), 2))[k]


def draw_conic(ax, c, span, color, label):
    """Draw conic.
    
    Args:
        ax:
        c:
        span:
        color:
        label:
    
    """
    if c is None:
        return
    Aq, Bq, Cq, Dq, Eq, Gq = [float(v) for v in c]
    xs = np.linspace(-span, span, 1200)
    for quad, lin, const in ((Cq, lambda x:
        Bq * x + Eq,
                              lambda x: Aq * x * x + Dq * x + Gq),
                             (Aq, lambda x: Dq + Bq * x,
                              lambda x: Cq * x * x + Eq * x + Gq)):
        for other in np.linspace(-span, span, 1600):
            b = lin(other)
            q = const(other)
            if quad == 0:
                if b != 0:
                    ax.plot([other], [-q / b], ".", color=color, ms=0.6)
                continue
            disc = b * b - 4 * quad * q
            if disc < 0:
                continue
            s = np.sqrt(disc)
            for t in (1.0, -1.0):
                val = (2 * q / (-b - t * s)) if (-b - t * s) != 0 else \
                    (-b + t * s) / (2 * quad)
                if -span <= val <= span:
                    ax.plot([other], [val], ".", color=color, ms=0.4, alpha=0.5)


def fig_scaling():
    """Chamber count against n, tolerant versus the classical baseline."""
    ns = list(range(2, 8))
    tolerant = [A.semicircle_bound(n) for n in ns]
    classical = [A.classical_region_count(n, 2) for n in ns]
    wide = list(range(2, 26))
    wt = [A.semicircle_bound(n) for n in wide]
    wc = [A.classical_region_count(n, 2) for n in wide]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    ax = axes[0]
    ax.plot(wide, wt, "-", color="#c1121f", lw=2,
            label=r"tolerant, $\delta$ small: $1+\frac{n^2(n-1)^2}{2}$")
    ax.plot(wide, wc, "-", color="#1d3557", lw=2,
            label=r"classical $\delta=0$: $\sum_{i=0}^{2}|s(n,n-i)|$")
    for n, t, c in zip(ns, tolerant, classical):
        ax.annotate(f"{t}", (n, t), textcoords="offset points",
                    xytext=(0, 7), fontsize=8, ha="center", color="#c1121f")
        ax.annotate(f"{c}", (n, c), textcoords="offset points",
                    xytext=(0, -14), fontsize=8, ha="center", color="#1d3557")
    ax.set_yscale("log")
    ax.set_xlabel("$n$")
    ax.set_ylabel("chambers $R$")
    ax.set_title("exact maximum chamber counts")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(wide, [t / c for t, c in zip(wt, wc)], "-", color="#2a9d8f", lw=2)
    for n, t, c in zip(ns, tolerant, classical):
        ax.plot(n, t / c, "o", color="#2a9d8f")
        ax.annotate(f"{t / c:.2f}", (n, t / c), textcoords="offset points",
                    xytext=(4, 4), fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel("$n$")
    ax.set_ylabel("tolerant / classical")
    ax.set_title(r"ratio $\to 4$ (both are $\Theta(n^4)$)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "scaling.png", dpi=150)
    plt.close(fig)
    print("wrote scaling.png")
    for n, t, c in zip(ns, tolerant, classical):
        print(f"  n={n}: tolerant {t}, classical {c}, ratio {t / c:.2f}")
    for n in (40, 80, 160, 320):
        print(f"  n={n}: ratio {A.semicircle_bound(n) / A.classical_region_count(n, 2):.4f}")


def fig_delta_sweep():
    """Chamber count as a function of tolerance, for n = 3 and n = 4."""
    configs = {
        3: [(F(0), F(0)), (F(5), F(1)), (F(1), F(6))],
        4: [(F(0), F(0)), (F(5), F(1)), (F(1), F(6)), (F(6), F(4))],
    }
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, (n, pts) in zip(axes, configs.items()):
        deltas = [F(k, 40) for k in range(1, 61)]
        counts = [A.arrangement_stats(pts, d).region_count for d in deltas]
        ax.plot([float(d) for d in deltas], counts, "-", color="#1d3557")
        ax.axhline(A.semicircle_bound(n), ls="--", color="#e63946",
                   label="Bézout ceiling")
        ax.axhline(A.classical_region_count(n, 2), ls=":", color="#457b9d",
                   label=r"classical $\delta=0$")
        ax.set_xlabel(r"$\delta$")
        ax.set_ylabel("chambers $R$")
        ax.set_title(f"$n={n}$")
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)
    fig.suptitle("chamber count versus tolerance: saturated for small "
                 r"$\delta$, then collapsing as curves go empty", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "delta_sweep.png", dpi=150)
    plt.close(fig)
    print("wrote delta_sweep.png")


def main() -> int:
    """Entry point — parse arguments and run the main computation.
    
    Returns:
        int: Result of type int
    
    """
    pts4 = [(F(0), F(0)), (F(5), F(1)), (F(1), F(6)), (F(6), F(4))]
    fig_arrangement(pts4[:3], 0.5,
                    "$n=3$: the Bezout ceiling is saturated",
                    "arrangement_n3.png", span=12.0)
    fig_arrangement(pts4, 0.5,
                    "$n=4$: 6 hyperbolas, 15 pairs, 4 crossings each",
                    "arrangement_n4.png", span=13.0)
    fig_scaling()
    fig_delta_sweep()
    return 0


if __name__ == "__main__":
    sys.exit(main())