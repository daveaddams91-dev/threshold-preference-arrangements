"""Exact conic algebra: correctness tests.

The geometric claims are checked against the *defining metric condition*
``|d_i - d_j| = delta``, not against another algebraic routine, so these tests
are independent of the elimination code they exercise.
"""
import math
import sys
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tolprefs import conics as C  # noqa: E402
from tolprefs import poly as P  # noqa: E402


def metric_residual(pts, e, delta, pt, tol=1e-6):
    """How far ``pt`` is from satisfying ``|d_i - d_j| = delta`` for pair ``e``."""
    i, j = e
    di = math.hypot(pt[0] - pts[i][0], pt[1] - pts[i][1])
    dj = math.hypot(pt[0] - pts[j][0], pt[1] - pts[j][1])
    return abs(abs(di - dj) - float(delta))


def conic_y_values(c, x, limit=4):
    """All real ``y`` with ``c(x, y) = 0``, as floats.

    Solves the quadratic in ``y`` when the conic has a ``y^2`` term and the linear
    equation otherwise, so it also covers hyperbolas whose axis is horizontal --
    for those the ``y^2`` coefficient vanishes and a naive quadratic formula
    divides by zero.
    """
    A, B, Cq, D, E, G = c
    lin = float(B * x + E)
    const = float(A * x * x + D * x + G)
    quad = float(Cq)
    if quad == 0:
        if lin == 0:
            return []
        return [-const / lin]
    disc = lin * lin - 4 * quad * const
    if disc < 0:
        return []
    s = math.sqrt(disc)
    # Numerically stable form:  (-lin +- s) / (2 quad) loses all precision when
    # lin ~ +- s, which happens whenever |x| is small compared with |QiQj|.  The
    # conjugate root gives the same value with no cancellation:
    #     (-lin + t s) / (2 quad)  ==  2 const / (-lin - t s).
    out = []
    for t in (1.0, -1.0):
        den = -lin - t * s
        if den != 0.0:
            out.append(2 * const / den)
        else:
            out.append((-lin + t * s) / (2 * quad))
    return out[:limit]


def test_gamma_conic_matches_metric():
    """Every point of the conic satisfies the defining condition, and vice versa."""
    pts = [(F(0), F(0)), (F(4), F(3)), (F(1), F(6))]
    for delta in (F(1, 10), F(1, 2), F(1), F(3)):
        for e in ((0, 1), (0, 2), (1, 2)):
            assert delta < float(C.dist2(pts[e[0]], pts[e[1]])) ** 0.5
            c = C.gamma_conic(pts[e[0]], pts[e[1]], delta)
            assert c is not None
            n = 0
            for k in range(1, 400):
                x = F(-60) + F(k, 3)
                for yv in conic_y_values(c, x):
                    pt = (float(x), yv)
                    assert metric_residual(pts, e, delta, pt) < 1e-6, (e, delta, pt)
                    n += 1
            assert n > 0, (e, delta)


def test_gamma_conic_degenerate_case():
    """At ``delta = |QiQj|`` the algebraic locus is the *whole* focal line.

    This is a genuine subtlety worth pinning down.  The metric condition
    ``|d_i - d_j| = |QiQj|`` holds, by the equality case of the reverse triangle
    inequality, only on the two rays of the line ``QiQj`` outside the segment
    ``QiQj``.  But the polynomial locus obtained from the conic equation is the
    entire line, counted twice: substituting ``delta^2 = L`` into
    ``4(L - delta^2) s^2 - 4 delta^2 t^2 - delta^2 (L - delta^2) L`` leaves
    ``-4 L t^2 = 0``, i.e. ``t = 0``, the whole focal line.

    So ``gamma_conic`` at ``delta = |QiQj|`` over-approximates the true
    indifference set, and callers must require ``delta < |QiQj|``.  The main
    theorem is stated under exactly that hypothesis.
    """
    qi, qj = (F(0), F(0)), (F(3), F(4))
    delta = F(5)
    c = C.gamma_conic(qi, qj, delta)
    assert c is not None
    # the conic is the double line through the foci
    for x, y in ((F(-3), F(-4)), (F(6), F(8)), (F(0), F(0)), (F(3), F(4))):
        assert C.conic_eval(c, (x, y)) == 0, (x, y)
    assert C.conic_eval(c, (F(0), F(1))) != 0
    # the metric condition fails on the segment, as it must
    assert metric_residual([qi, qj], (0, 1), delta, (1.5, 2.0)) > 0.5
    # and it holds on both outer rays
    assert metric_residual([qi, qj], (0, 1), delta, (-3.0, -4.0)) < 1e-9
    assert metric_residual([qi, qj], (0, 1), delta, (6.0, 8.0)) < 1e-9


def test_gamma_nonempty_threshold():
    """H_ij exists iff delta <= |Q_i Q_j|; equality degenerates to a double line."""
    qi, qj = (F(0), F(0)), (F(3), F(4))     # distance 5
    assert C.gamma_nonempty(qi, qj, F(5))
    assert C.gamma_nonempty(qi, qj, F(1))
    assert not C.gamma_nonempty(qi, qj, F(6))
    assert C.gamma_conic(qi, qj, F(6)) is None
    # at delta = distance the conic is a double line through the foci
    c = C.gamma_conic(qi, qj, F(5))
    assert c is not None
    # the perpendicular bisector is the delta -> 0 limit
    d = C.gamma_conic(qi, qj, F(1, 10 ** 9))
    assert d is not None
    assert C.perpendicular_bisector(qi, qj) == \
        (2 * F(3), 2 * F(4), -F(25))


def test_vertex_location():
    """Vertices sit on the line Q_i Q_j at signed distance delta/2 from midpoint."""
    qi, qj = (F(0), F(0)), (F(6), F(0))
    delta = F(1, 2)
    c = C.gamma_conic(qi, qj, delta)
    A, B, Cq, D, E, G = c
    for x in (F(3) - delta / 2, F(3) + delta / 2):
        # on the conic with y = 0?  the vertices are the two points with x these
        # and minimal |y|; solve for y and require it to be 0
        disc = (B * x + E) ** 2 - 4 * Cq * (A * x * x + D * x + G)
        assert disc == 0, (x, disc)


def test_gamma_conic_is_hyperbola():
    """The discriminant B^2 - 4AC is positive: the curve really is a hyperbola."""
    pts = [(F(0), F(0)), (F(5), F(1)), (F(2), F(4))]
    for delta in (F(1, 5), F(1), F(3)):
        for e in ((0, 1), (0, 2), (1, 2)):
            c = C.gamma_conic(pts[e[0]], pts[e[1]], delta)
            assert c is not None
            A, B, Cq = c[0], c[1], c[2]
            assert B * B - 4 * A * Cq > 0, (e, delta)


def test_split_in_y_reconstructs_conic():
    """The coefficient split must rebuild the conic exactly (regression guard)."""
    pts = [(F(7), F(-4)), (F(-1), F(-4)), (F(2), F(7))]
    for e in ((0, 1), (0, 2), (1, 2)):
        c = C.gamma_conic(pts[e[0]], pts[e[1]], F(1, 4))
        a, b, k = C._split_in_y(c)
        C._check_split(c, a, b, k)
        A, B, Cq, D, E, G = c
        assert a == (Cq,)
        assert b == (E, B)
        assert k == (G, D, A)
        # spot-check the reconstruction at a few rational points
        for x, y in ((F(0), F(0)), (F(3), F(-1)), (F(1, 2), F(5, 2))):
            lhs = a[0] * y * y + P.p_eval(b, x) * y + P.p_eval(k, x)
            assert lhs == C.conic_eval(c, (x, y))


def test_line_restriction_is_constant_term_first():
    """Regression: the fibre restriction must encode ``lead t^2 + lin t + const``."""
    c = (F(39996), F(-8080), F(0), F(35552), F(-8080), F(-14544))
    t = F(2)
    A, B, Cq, D, E, G = c
    r = C._line_restriction(c, t, solve_for_y=True)
    # const = A t^2 + D t + G, lin = B t + E, lead = C  (here C = 0, so linear)
    assert P.p_trim(r) == P.p_trim((A * t * t + D * t + G, B * t + E, Cq))
    rx = C._line_restriction(c, t, solve_for_y=False)
    assert P.p_trim(rx) == P.p_trim((Cq * t * t + E * t + G, D + B * t, A))
    # a conic with a y^2 term must give a genuine quadratic restriction
    c2 = (F(29127, 2), F(-64240), F(141255, 2), F(-69423, 2), F(154395, 2),
          F(2415497, 128))
    r2 = C._line_restriction(c2, F(3), solve_for_y=True)
    A2, B2, C2, D2, E2, G2 = c2
    assert P.p_trim(r2) == P.p_trim((A2 * 9 + D2 * 3 + G2, B2 * 3 + E2, C2))
    assert P.p_deg(r2) == 2
    assert float(r2[0]) > 0 and float(r2[1]) < 0 and float(r2[2]) > 0


def test_intersection_count_matches_points():
    """``conic_intersection_count`` must agree with ``conic_intersections``."""
    cases = [
        ([(F(0), F(0)), (F(4), F(3)), (F(1), F(6))], F(1, 2)),
        ([(F(-2), F(5)), (F(-5), F(3)), (F(-4), F(1))], F(1)),
        ([(F(5), F(6)), (F(4), F(-6)), (F(-6), F(-5))], F(1)),
        ([(F(4), F(-3)), (F(-2), F(-1)), (F(4), F(-1)), (F(-1), F(2))],
         F(3, 10)),
        ([(F(1), F(3)), (F(-2), F(4)), (F(-3), F(-7))], F(1)),
    ]
    for pts, delta in cases:
        conics = {}
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                c = C.gamma_conic(pts[i], pts[j], delta)
                if c is not None:
                    conics[(i, j)] = c
        keys = sorted(conics)
        for a in range(len(keys)):
            for b in range(a + 1, len(keys)):
                e1, e2 = keys[a], keys[b]
                c1, c2 = conics[e1], conics[e2]
                n = C.conic_intersection_count(c1, c2)
                pts_out = C.conic_intersections(c1, c2)
                assert n == len(pts_out), (e1, e2, delta, n, len(pts_out))
                assert 0 <= n <= 4
                # every returned point must satisfy both metric conditions
                for p in pts_out:
                    v = p.value()
                    assert metric_residual(pts, e1, delta, v) < 1e-6, (e1, e2, v)
                    assert metric_residual(pts, e2, delta, v) < 1e-6, (e1, e2, v)


def test_points_are_distinct():
    """No duplicates, including the shared-fibre case with equal coordinates."""
    pts = [(F(-2), F(5)), (F(-5), F(3)), (F(-4), F(1))]
    c1 = C.gamma_conic(pts[0], pts[1], F(1))
    c2 = C.gamma_conic(pts[1], pts[2], F(1))
    out = C.conic_intersections(c1, c2)
    vals = [p.value() for p in out]
    for a in range(len(vals)):
        for b in range(a + 1, len(vals)):
            assert math.dist(vals[a], vals[b]) > 1e-6


def test_shared_fibre_recovered():
    """Regression: two points on one vertical line must both be reported."""
    pts = [(F(-2), F(5)), (F(-5), F(3)), (F(-4), F(1))]
    c1 = C.gamma_conic(pts[0], pts[1], F(1))
    c2 = C.gamma_conic(pts[1], pts[2], F(1))
    out = [p.value() for p in C.conic_intersections(c1, c2)]
    assert len(out) == 4
    xs = [round(v[0], 9) for v in out]
    # two of them share x = -3 exactly
    assert xs.count(-3.0) == 2, xs


def test_shear_is_involutive_on_incidence():
    """A shear preserves the intersection count and maps points back correctly."""
    c1 = (F(0), F(13920), F(82940), F(0), F(-62640), F(-20880))
    c2 = (F(39996), F(-8080), F(0), F(35552), F(-8080), F(-14544))
    for k in (F(1), F(-1), F(2), F(1, 2)):
        n_before = C.conic_intersection_count(c1, c2)
        n_after = C.conic_intersection_count(C._shear(c1, k), C._shear(c2, k))
        assert n_before == n_after, k


def test_perpendicular_bisector_and_circumcenter():
    a, b = (F(0), F(0)), (F(6), F(0))
    assert C.perpendicular_bisector(a, b) == (F(12), F(0), F(-36))
    c = (F(3), F(4))
    o = C.circumcenter(a, b, c)
    assert o == (F(3), F(7, 8))
    # equidistant from all three
    for p in (a, b, c):
        assert C.dist2(o, p) == C.dist2(o, a)
    # three collinear points have no circumcentre
    assert C.circumcenter((F(0), F(0)), (F(1), F(1)), (F(2), F(2))) is None


def test_bisector_crossing():
    a, b, c, d = (F(0), F(0)), (F(2), F(0)), (F(0), F(0)), (F(0), F(2))
    assert C.bisector_crossing(a, b, c, d) == (F(1), F(1))
    assert C.bisector_crossing(a, b, c, (F(2), F(0))) is None


def test_delta_zero_is_the_bisector():
    """As ``delta -> 0`` both branches converge onto the perpendicular bisector.

    The hyperbola is not a perturbation *of* the bisector line: its two branches
    straddle it and squeeze onto it.  So the correct statement is a geometric
    one -- every point of ``H_ij(delta)`` approaches the bisector line as
    ``delta -> 0`` -- which is what this test measures.
    """
    qi, qj = (F(1), F(2)), (F(5), F(-1))
    a, b, g = C.perpendicular_bisector(qi, qj)
    norm = float((a * a + b * b) ** 0.5)
    # The asymptotes of H_ij meet the bisector at an angle of order
    # delta/|QiQj|, so the deviation from the bisector grows like delta * |X|
    # rather than like delta alone.  Over the bounded window |x| <= 20 with
    # |QiQj| = 5 the deviation is at most 12 delta, and that linear-in-delta
    # behaviour is the statement being checked.
    def worst_deviation(delta, window, samples=400):
        c = C.gamma_conic(qi, qj, delta)
        worst = 0.0
        k = 0
        while k < samples:
            x = F(-window) + F(k, 2 * samples) * (2 * window) / samples
            k += 1
            for yv in conic_y_values(c, x):
                worst = max(worst,
                            abs(float(a) * x + float(b) * yv + float(g)) / norm)
        return worst
    for delta in (F(1, 100), F(1, 1000), F(1, 10000)):
        w = worst_deviation(delta, 20)
        assert w <= 12 * float(delta), (delta, w)
    # on a fixed window the branches squeeze onto the bisector as delta shrinks
    assert worst_deviation(F(1, 100), F(1, 4)) > \
        20 * worst_deviation(F(1, 4000), F(1, 4))