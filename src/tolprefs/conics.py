"""Exact algebra for the threshold-indifference curves.

The mathematical objects
------------------------
Let ``Q = {Q_0, ..., Q_{n-1}}`` be a finite set of points in the plane and let
``delta > 0`` be a tolerance.  For a pair ``{i, j}`` of distinct indices define
the *indifference curve*

    H_ij = { X : | d(X, Q_i) - d(X, Q_j) | = delta }.

``H_ij`` is a hyperbola with foci ``Q_i`` and ``Q_j``: semi-transverse axis
``a = delta/2`` along the line ``Q_i Q_j`` and semi-conjugate axis
``b = sqrt(|Q_i Q_j|^2 - delta^2)/2``.  Its two branches are the loci on which
one candidate is exactly ``delta`` closer than the other; together they are the
boundary of the region where the preference between ``i`` and ``j`` is
*neutral*.

Three facts drive the whole analysis.

1. **Non-emptiness.**  Since ``|d(X,Q_i) - d(X,Q_j)| <= |Q_i Q_j|`` for all ``X``,
   with equality attained on the line ``Q_i Q_j``, the curve ``H_ij`` exists iff
   ``delta <= |Q_i Q_j|``.
2. **Position.**  ``H_ij`` lies in ``{X : d(X,Q_i) != d(X,Q_j)}``, and its
   vertices lie on the line ``Q_i Q_j``, at signed distance ``delta/2`` from the
   midpoint.
3. **Degenerate limit.**  As ``delta -> 0+`` the two branches merge into the
   perpendicular bisector of ``Q_i Q_j``, which recovers the classical
   Good--Tideman arrangement of that name.

Exact representation
--------------------
Everything is computed over ``Q``.  A conic is the 6-tuple ``(A,B,C,D,E,F)``
denoting

    A x^2 + B x y + C y^2 + D x + E y + F = 0.

To obtain the equation, put ``w = Q_j - Q_i``, ``L = |w|^2``,
``m = (Q_i + Q_j)/2``, ``s = (X - m) . w``, ``t = (X - m) . w_perp`` with
``w_perp = (-w_y, w_x)``.  In the orthonormal frame aligned with ``w`` the point
has coordinates ``u = s/sqrt(L)``, ``v = t/sqrt(L)``, and ``H_ij`` is

    u^2/a^2 - v^2/b^2 = 1,    a = delta/2,   b^2 = (L - delta^2)/4.

Substituting and clearing denominators yields the rational conic

    4 (L - delta^2) s^2  -  4 delta^2 t^2  -  delta^2 (L - delta^2) L = 0.

Counting intersections
----------------------
:func:`conic_intersection_count` returns the number of *distinct real*
intersection points of two conics, exactly.  The elimination gives a quartic
``H(x)`` whose real roots are the candidate ``x``-coordinates, and one subtlety
must be handled: when the pencil relation degenerates (``Qn(x) = 0``) the two
conics share a vertical fibre and contribute *two* points with the same
``x``.  Because ``Qn`` is linear its root is rational, so this degenerate case
is always decidable in exact rational arithmetic -- no minimal polynomial of
``x`` is ever required.
"""

from __future__ import annotations

from . import poly as P
from fractions import Fraction as F



Point = tuple[F, F]
Conic = tuple[F, F, F, F, F, F]

__all__ = [
    "F",
    "Conic",
    "Point",
    "AlgebraicPoint",
    "as_point",
    "add",
    "sub",
    "dot",
    "norm2",
    "dist2",
    "gamma_nonempty",
    "gamma_conic",
    "bisector_conic",
    "perpendicular_bisector",
    "circumcenter",
    "bisector_crossing",
    "conic_eval",
    "conic_intersections",
    "conic_intersection_count",
    "triangle_areas",
]


# ------------------------------------------------------------- basics


def as_point(p) -> Point:
    """As point.
    
    Args:
        p:
    
    Returns:
        tuple: Result of type tuple
    
    """
    return (F(p[0]), F(p[1]))


def add(p: Point, q: Point) -> Point:
    """Add.
    
    Args:
        p:
        q:
    
    Returns:
        tuple: Result of type tuple
    
    """
    return (p[0] + q[0], p[1] + q[1])


def sub(p: Point, q: Point) -> Point:
    """Sub.
    
    Args:
        p:
        q:
    
    Returns:
        tuple: Result of type tuple
    
    """
    return (p[0] - q[0], p[1] - q[1])


def dot(p: Point, q: Point) -> F:
    """Dot.
    
    Args:
        p:
        q:
    
    Returns:
        The computed result
    
    """
    return p[0] * q[0] + p[1] * q[1]


def norm2(p: Point) -> F:
    """Norm2.
    
    Args:
        p:
    
    Returns:
        The computed result
    
    """
    return dot(p, p)


def dist2(p: Point, q: Point) -> F:
    """Dist2.
    
    Args:
        p:
        q:
    
    Returns:
        The computed result
    
    """
    return norm2(sub(p, q))


# ------------------------------------------------------------- conics


def _normalise(c: Conic) -> Conic:
    """Normalise.
    
    Args:
        c:
    
    Returns:
        The computed result
    
    """
    for v in c:
        if v != 0:
            return tuple(-w for w in c) if v < 0 else c  # type: ignore[return-value]
    return c


def gamma_nonempty(qi: Point, qj: Point, delta: F) -> bool:
    """``H_ij`` has a real point iff ``delta <= |Q_i Q_j|``."""
    return F(delta) * F(delta) <= dist2(qi, qj)


def gamma_conic(qi: Point, qj: Point, delta: F) -> Conic | None:
    """Implicit conic of ``H_ij``; ``None`` when the curve is empty.

    **The conic equation is faithful to the metric locus only when
    ``delta < |Q_i Q_j|``, and callers must enforce this.**  For
    ``0 < delta < |Q_i Q_j|`` the conic is a genuine hyperbola with foci
    ``Q_i, Q_j`` and is exactly ``{X : |d(X,Q_i) - d(X,Q_j)| = delta}``.

    At the endpoint ``delta = |Q_i Q_j|`` the metric locus is, by the equality
    case of the reverse triangle inequality, only the two rays of the line
    ``Q_iQ_j`` *outside* the segment; but substituting ``delta^2 = L`` into the
    equation below collapses it to ``t = 0``, the whole focal line counted
    twice.  The returned conic therefore *over*-approximates the true
    indifference set at that single value, and the discrepancy is not numerical.
    For ``delta > |Q_i Q_j|`` the locus is empty and ``None`` is returned.
    """
    delta = F(delta)
    L = dist2(qi, qj)
    if delta * delta > L:
        return None
    w = sub(qj, qi)
    m = ((qi[0] + qj[0]) / 2, (qi[1] + qj[1]) / 2)
    wp = (-w[1], w[0])

    # s = (X - m) . w  and  t = (X - m) . w_perp, as affine forms in (x, y)
    sx, sy = w[0], w[1]
    s0 = -(m[0] * w[0] + m[1] * w[1])
    tx, ty = wp[0], wp[1]
    t0 = -(m[0] * wp[0] + m[1] * wp[1])

    k1 = F(4) * (L - delta * delta)  # multiplies s^2
    k2 = -F(4) * delta * delta  # multiplies t^2
    k0 = -delta * delta * (L - delta * delta) * L

    conic = _normalise(
        (
            k1 * sx * sx + k2 * tx * tx,
            F(2) * (k1 * sx * sy + k2 * tx * ty),
            k1 * sy * sy + k2 * ty * ty,
            F(2) * (k1 * s0 * sx + k2 * t0 * tx),
            F(2) * (k1 * s0 * sy + k2 * t0 * ty),
            k1 * s0 * s0 + k2 * t0 * t0 + k0,
        )
    )
    if all(v == 0 for v in conic):
        return None
    return conic


def bisector_conic(qi: Point, qj: Point) -> Conic:
    """The ``delta = 0`` locus ``{d(X,Q_i) = d(X,Q_j)}`` as a degenerate conic."""
    w = sub(qj, qi)
    return (F(0), F(0), F(0), F(2) * w[0], F(2) * w[1], -(norm2(qj) - norm2(qi)))


def perpendicular_bisector(qi: Point, qj: Point) -> tuple[F, F, F]:
    """Perpendicular bisector.
    
    Args:
        qi:
        qj:
    
    Returns:
        tuple: Result of type tuple
    
    """
    con = bisector_conic(qi, qj)
    return (con[3], con[4], con[5])


def circumcenter(qi: Point, qj: Point, qk: Point) -> Point | None:
    """Circumcenter.
    
    Args:
        qi:
        qj:
        qk:
    
    Returns:
        tuple: Result of type tuple
    
    """
    a, b, c = perpendicular_bisector(qi, qj)
    d, e, f = perpendicular_bisector(qi, qk)
    det = a * e - b * d
    if det == 0:
        return None
    return ((b * f - c * e) / det, (c * d - a * f) / det)


def bisector_crossing(qi: Point, qj: Point, qk: Point, ql: Point) -> Point | None:
    """Bisector crossing.
    
    Args:
        qi:
        qj:
        qk:
        ql:
    
    Returns:
        tuple: Result of type tuple
    
    """
    a1, b1, c1 = perpendicular_bisector(qi, qj)
    a2, b2, c2 = perpendicular_bisector(qk, ql)
    det = a1 * b2 - b1 * a2
    if det == 0:
        return None
    return ((b1 * c2 - c1 * b2) / det, (c1 * a2 - a1 * c2) / det)


def triangle_areas(pts) -> list[F]:
    """Signed double areas of consecutive triples, cyclically.

    Used to certify general position: a configuration is in general position
    (in the sense needed here) when no three points are collinear, no two
    connecting segments are parallel, and no four points are concyclic.
    """
    n = len(pts)
    out: list[F] = []
    for i in range(n):
        a, b, c = pts[i], pts[(i + 1) % n], pts[(i + 2) % n]
        out.append(abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])))
    return out


def conic_eval(c: Conic, p: Point) -> F:
    """Conic eval.
    
    Args:
        c:
        p:
    
    Returns:
        The computed result
    
    """
    A, B, C, D, E, G = c
    x, y = p
    return A * x * x + B * x * y + C * y * y + D * x + E * y + G


# --------------------------------------------- splitting a conic in y


def _split_in_y(c: Conic) -> tuple[P.Poly, P.Poly, P.Poly]:
    """Write ``c`` as ``A(x) y^2 + B(x) y + C(x)`` with ``A, B, C`` in ``Q[x]``.

    Since ``c`` denotes ``A2 x^2 + B xy + C y^2 + D x + E y + G``, and
    polynomials are stored constant-term-first, we get

        A(x) = C,     B(x) = E + B x,     C(x) = G + D x + A2 x^2.

    The ordering of ``B`` is the subtle part: the tuple is ``(E, B)``, *not*
    ``(B, E)``.
    """
    Aq, Bq, Cq, D, E, G = c
    return ((Cq,), (E, Bq), (G, D, Aq))


def _check_split(c: Conic, a: P.Poly, b: P.Poly, k: P.Poly) -> None:
    """Assert a coefficient split really reconstructs the conic.

    Guards the constant-term-first convention.  A reversed tuple produces a
    self-consistent but wrong polynomial, so an internal identity check on the
    eliminant cannot detect it; this assertion catches it at construction.
    This corresponds to a bug found during development, where ``(B, E)`` was used
    in place of ``(E, B)`` and every eliminant was silently wrong.
    """
    A, B, C, D, E, G = c
    variants = (((C,), (E, B), (G, D, A)),)  # split in y
    for want_a, want_b, want_k in variants:
        if a == want_a and b == want_b and k == want_k:
            return
    raise AssertionError(
        f"coefficient split mis-mapped conic {c}: got a={a}, b={b}, k={k}"
    )


def _check_const_first(coeffs: P.Poly, lead: F, lin: F, const: F,
                       what: str) -> None:
    """Assert ``coeffs`` really encodes ``lead t^2 + lin t + const``."""
    want = P.p_trim((const, lin, lead))
    if P.p_trim(coeffs) != want:  # pragma: no cover - defensive
        raise AssertionError(
            f"{what}: expected constant-term-first {want}, got {coeffs}"
        )


# ------------------------------------------------------- elimination


def _eliminate(c1: Conic, c2: Conic) -> tuple[P.Poly, P.Poly, P.Poly] | None:
    """Eliminate ``y`` between two conics, or return ``None``.

    Write ``F = a1 y^2 + b1 y + k1`` and ``G = a2 y^2 + b2 y + k2`` with the
    coefficients in ``Q[x]``.  The combination ``a2 F - a1 G`` cancels the
    ``y^2`` term, leaving

        (a2 b1 - a1 b2) y + (a2 k1 - a1 k2) = 0.

    Hence ``y Qn = -const = a1 k2 - a2 k1``, i.e.

        Qn = a2 b1 - a1 b2,      Pn = a1 k2 - a2 k1,

    and substituting ``y = Pn/Qn`` back into ``F`` gives
    ``H = a1 Pn^2 + b1 Pn Qn + k1 Qn^2``, of degree at most 4 by Bezout.  The
    real roots of ``H`` are the candidate ``x``-coordinates of common points.
    """
    a1, b1, k1 = _split_in_y(c1)
    a2, b2, k2 = _split_in_y(c2)
    _check_split(c1, a1, b1, k1)
    _check_split(c2, a2, b2, k2)
    # If either conic is linear in y (no y^2 term) the pencil construction
    # degenerates: the identity a2*F - a1*G becomes trivial and the eliminant
    # collapses.  This happens whenever a hyperbola has a horizontal axis in the
    # chosen coordinates, so the caller must be able to switch variables.
    if P.p_deg(a1) < 0 or P.p_deg(a2) < 0:
        return None
    Qn = P.p_sub(P.p_mul(a2, b1), P.p_mul(a1, b2))
    if P.p_deg(Qn) < 0:
        return None
    Pn = P.p_sub(P.p_mul(a1, k2), P.p_mul(a2, k1))
    H = P.p_add(
        P.p_add(P.p_mul(a1, P.p_mul(Pn, Pn)), P.p_mul(b1, P.p_mul(Pn, Qn))),
        P.p_mul(k1, P.p_mul(Qn, Qn)),
    )
    if P.p_deg(H) < 1:
        return None
    return Pn, Qn, H


def _eliminate_x(c1: Conic, c2: Conic) -> tuple[P.Poly, P.Poly, P.Poly] | None:
    """Same elimination with ``x`` and ``y`` interchanged (fallback path)."""

    def split_in_x(c: Conic) -> tuple[P.Poly, P.Poly, P.Poly]:
        """Split in x.
        
        Args:
            c:
        
        Returns:
            tuple: Result of type tuple
        
        """
        A, B, C, D, E, G = c
        return ((A,), (D, B), (G, E, C))

    a1, b1, k1 = split_in_x(c1)
    a2, b2, k2 = split_in_x(c2)
    if P.p_deg(a1) < 0 or P.p_deg(a2) < 0:
        return None
    Qn = P.p_sub(P.p_mul(a2, b1), P.p_mul(a1, b2))
    if P.p_deg(Qn) < 0:
        return None
    Pn = P.p_sub(P.p_mul(a1, k2), P.p_mul(a2, k1))
    H = P.p_add(
        P.p_add(P.p_mul(a1, P.p_mul(Pn, Pn)), P.p_mul(b1, P.p_mul(Pn, Qn))),
        P.p_mul(k1, P.p_mul(Qn, Qn)),
    )
    if P.p_deg(H) < 1:
        return None
    return Pn, Qn, H


class AlgebraicPoint:
    """An exact point of the plane.

    ``x`` is the unique real root of ``poly`` inside ``interval``, and ``y`` is
    given either as the polynomial ``ycoeffs`` in ``x`` (an element of the field
    ``Q[x]/(poly)``) or, when ``x`` is rational, by its own minimal polynomial
    ``ypoly`` inside ``yinterval`` (or by the rational ``y0``).
    """

    __slots__ = ("poly", "ycoeffs", "interval", "ypoly", "yinterval", "x0", "y0",
                 "swapped", "solve_for_y", "shear")

    def __init__(self, poly, ycoeffs, interval, ypoly=None, yinterval=None,
                 x0=None, y0=None, swapped=False, solve_for_y=True):
        """Init.
        
        Args:
            poly:
            ycoeffs:
            interval:
            ypoly:
            yinterval:
            x0:
            y0:
            swapped (bool):
            solve_for_y (bool):
        
        """
        self.poly = poly
        self.ycoeffs = ycoeffs
        self.interval = interval
        self.ypoly = ypoly
        self.yinterval = yinterval
        self.x0 = x0
        self.y0 = y0
        self.swapped = swapped
        self.solve_for_y = solve_for_y
        self.shear = F(0)

    def key(self) -> tuple:
        """Key.
        
        Returns:
            tuple: Result of type tuple
        
        """
        return (self.poly, self.ycoeffs, self.interval, self.ypoly,
                self.yinterval, self.x0, self.y0, self.swapped,
                self.solve_for_y)

    def __eq__(self, other):
        """Eq.
        
        Args:
            other:
        
        Returns:
            The computed result
        
        """
        return isinstance(other, AlgebraicPoint) and self.key() == other.key()

    def __hash__(self):
        """Hash.
        
        Returns:
            The computed result
        
        """
        return hash(self.key())

    def value(self) -> tuple[float, float]:
        """Numeric ``(x, y)``, undoing any axis swap and variable convention."""
        # the isolating interval refers to the *rooted* coordinate
        lo, hi = self.interval
        root = float((lo + hi) / 2)
        if self.y0 is not None:
            other = float(self.y0)
        elif self.ypoly is not None and self.yinterval is not None:
            ylo, yhi = self.yinterval
            other = float((ylo + yhi) / 2)
        else:
            other = float(sum(c * root**k for k, c in enumerate(self.ycoeffs)))
        x, y = (root, other) if self.solve_for_y else (other, root)
        if self.swapped:
            x, y = y, x
        if self.shear != 0:
            # the forward map was x -> x + k y, so the inverse is x = X + k Y
            x = x + float(self.shear) * y
        return x, y

    def unshear(self, k: F) -> "AlgebraicPoint":
        """Map back from the sheared frame ``(X, Y) = (x + k y, y)``.

        The forward map sends ``(x, y)`` to ``(X, Y) = (x + k y, y)``, so the
        inverse is ``x = X + k Y``.  That mixes the two coordinates, and the
        result generally has *both* coordinates algebraic, which this
        representation cannot express exactly.  Rather than return a silently
        wrong exact object, we carry the shear parameter and apply it only when
        producing numeric values -- which is all the callers need, and the
        numeric pair then satisfies the original conics.
        """
        if k == 0:
            return self
        clone = AlgebraicPoint(self.poly, self.ycoeffs, self.interval,
                               self.ypoly, self.yinterval, self.x0, self.y0,
                               self.swapped, self.solve_for_y)
        clone.shear = k
        return clone

    def __repr__(self) -> str:  # pragma: no cover
        """Repr.
        
        Returns:
            The computed result
        
        """
        return f"AP(x~{self.value()[0]:.6g}, y~{self.value()[1]:.6g})"


def _line_restriction(c: Conic, t: F, solve_for_y: bool) -> P.Poly:
    """The conic restricted to the line ``x = t`` (or ``y = t``), in the other
    variable, stored constant-term-first.

    For ``solve_for_y`` the restriction is ``A y^2 + B y + K`` and the tuple is
    ``(K, B, A)``; for ``solve_for_x=False`` the roles of ``x`` and ``y`` are
    exchanged.  The ordering is asserted because a silent reversal here yields
    a plausible-looking but wrong gcd.
    """
    A, B, C, D, E, G = c
    if solve_for_y:
        lead = C
        lin = B * t + E
        const = A * t * t + D * t + G
    else:
        lead = A
        lin = D + B * t
        const = C * t * t + E * t + G
    out = P.p_trim((const, lin, lead))
    _check_const_first(out, lead, lin, const, "line restriction")
    return out


def _fibre_gcd(c1: Conic, c2: Conic, t: F, solve_for_y: bool) -> P.Poly:
    """Fibre gcd.
    
    Args:
        c1:
        c2:
        t:
        solve_for_y:
    
    Returns:
        The computed result
    
    """
    f1 = _line_restriction(c1, t, solve_for_y)
    f2 = _line_restriction(c2, t, solve_for_y)
    return P.p_gcd(P.p_trim(f1), P.p_trim(f2))


def _fibre(c1: Conic, c2: Conic, t: F, solve_for_y: bool = True,
           swapped: bool = False) -> list[AlgebraicPoint]:
    """Common points of the two conics on the line ``x = t`` (or ``y = t``).

    Used exactly when the pencil degenerates.  Both restrictions are univariate
    with rational coefficients, so their gcd is exact and its real roots are
    the desired remaining coordinates.
    """
    g = _fibre_gcd(c1, c2, t, solve_for_y)
    if P.p_deg(g) < 1:
        return []
    out: list[AlgebraicPoint] = []
    for lo, hi in P.real_roots(g):
        if lo == hi:
            out.append(AlgebraicPoint((F(1),), None, (t, t), x0=t, y0=lo,
                                      solve_for_y=solve_for_y, swapped=swapped))
        else:
            out.append(AlgebraicPoint((F(1),), None, (t, t), x0=t, ypoly=g,
                                      yinterval=(lo, hi),
                                      solve_for_y=solve_for_y, swapped=swapped))
    return out


def _fibre_count(c1: Conic, c2: Conic, t: F, solve_for_y: bool = True) -> int:
    """Fibre count.
    
    Args:
        c1:
        c2:
        t:
        solve_for_y (bool):
    
    Returns:
        int: Result of type int
    
    """
    g = _fibre_gcd(c1, c2, t, solve_for_y)
    if P.p_deg(g) < 1:
        return 0
    return len(P.real_roots(g))


def _lift(lo: F, hi: F, H: P.Poly, Pn: P.Poly, Qn: P.Poly,
          c1: Conic, c2: Conic, solve_for_y: bool = True,
          swapped: bool = False) -> list[AlgebraicPoint]:
    """Exact preimages of the root of ``H`` isolated by ``(lo, hi)``.

    ``H``'s root is the coordinate named ``x`` when the elimination solved for
    ``y`` (``solve_for_y=True``) and the coordinate named ``y`` otherwise.  If
    ``Qn`` does not vanish there, the other coordinate is ``Pn/Qn`` uniquely.
    When ``Qn`` vanishes, that root is rational (``Qn`` is linear) and the two
    conics share a fibre, contributing one point per real root of the gcd of the
    two restrictions to that line.
    """
    x0 = _qn_zero(Qn)
    if x0 is not None and lo <= x0 <= hi and P.p_eval(H, x0) == 0:
        return _fibre(c1, c2, x0, solve_for_y=solve_for_y, swapped=swapped)
    m = _irreducible_part_at(H, lo, hi, Qn)
    inv_q = P.invert_mod(P.p_trim(Qn), m)
    if inv_q is None:  # pragma: no cover - defensive
        return []
    other = P.p_divmod(P.p_mul(P.p_trim(Pn), inv_q), m)[1]
    if P.p_trim(other) == ():
        return []  # the recovered coordinate is 0; handled by reduction
    return [AlgebraicPoint(m, P.p_trim(other), (lo, hi),
                           solve_for_y=solve_for_y, swapped=swapped)]


def _qn_zero(Qn: P.Poly) -> F | None:
    """The rational root of the (linear or constant) ``Qn``, if it has one."""
    d = P.p_deg(Qn)
    if d < 0:
        return None
    if d == 0:
        return None
    return -Qn[0] / Qn[1]


def _irreducible_part_at(a: P.Poly, lo: F, hi: F,
                         Qn: P.Poly | None = None) -> P.Poly:
    """The irreducible factor of ``a`` carrying the root isolated by ``(lo, hi)``.

    Working modulo the *whole* square-free eliminant is wrong whenever that
    eliminant is reducible: the field ``Q[x]/(m)`` then decomposes, and an
    element that is invertible in one component can be a zero divisor overall,
    so ``invert_mod`` returns ``None`` and the point is silently lost.  Splitting
    out the single factor that has this root restores a genuine field.

    A cheap and sufficient way to obtain that factor: for a rational candidate
    ``x0`` inside the isolating interval use ``gcd(a, x - x0)``; otherwise use
    the square-free part, which is irreducible unless ``a`` factors, and in that
    case fall back to the linear rational-root split below.
    """
    sf = P.p_squarefree(a)
    rats = P.rational_roots(sf)
    for r in rats:
        if lo <= r <= hi:
            g = P.p_gcd(sf, (-r, F(1)))
            if P.p_deg(g) >= 1:
                return g
    # ``_factors_in_interval`` only splits quadratics and cubics with *cheap*
    # rational roots.  When the eliminant has a rational root whose coefficients
    # are too large for that search, the split is missed and the whole
    # eliminant is returned, reintroducing zero divisors.  Dividing by the one
    # rational linear factor we know about -- Qn's root, which is rational
    # because Qn is linear -- recovers the true factor cheaply and exactly.
    x0 = _qn_zero(Qn) if Qn is not None else None
    if x0 is not None and P.p_eval(sf, x0) == 0:
        g = P.p_gcd(sf, (-x0, F(1)))
        if P.p_deg(g) >= 1:
            rest, rem = P.p_divmod(sf, g)
            if not rem:
                if _has_root_in(rest, lo, hi):
                    return rest
                return g
    keep = [f for f in _factors_in_interval(sf, lo, hi)
            if _has_root_in(f, lo, hi)]
    if keep:
        prod = keep[0]
        for f in keep[1:]:
            prod = P.p_mul(prod, f)
        return prod
    # Last resort: split off any factor that shares a root with Qn's root, so
    # that the returned modulus is invertible there even if the split above was
    # inconclusive (e.g. a cubic that is irreducible over Q but has Qn's root
    # among its roots -- impossible, but the defensive step keeps the caller
    # from silently dropping a point).
    return sf


def _minpoly_split(sf: P.Poly, Qn: P.Poly, lo: F, hi: F) -> P.Poly:
    """A factor of ``sf`` carrying the root in ``(lo, hi)`` and coprime to ``Qn``."""
    factors = _factors_in_interval(sf, lo, hi)
    keep = [f for f in factors if _has_root_in(f, lo, hi)]
    if not keep:
        return sf
    prod = keep[0]
    for f in keep[1:]:
        prod = P.p_mul(prod, f)
    if P.p_deg(P.p_gcd(P.p_trim(prod), P.p_trim(Qn))) > 0:
        # Qn vanishes at this root; pick the factor that does not
        for f in keep:
            if P.p_deg(P.p_gcd(P.p_trim(f), P.p_trim(Qn))) == 0:
                return f
    return prod


def _has_root_in(f: P.Poly, lo: F, hi: F) -> bool:
    """Whether the square-free ``f`` has a real root in the open interval."""
    chain = P.p_sturm_chain(P.p_squarefree(f))
    if not chain:
        return False
    # Sturm counts roots in (lo, hi) as V(lo) - V(hi); guard against lo/hi
    # themselves being roots by shrinking slightly.
    eps = F(1, 10**9)
    a, b = lo + eps, hi - eps
    if not (a < b):
        a, b = lo, hi
    return P._variations(chain, a) - P._variations(chain, b) >= 1


def _factors_in_interval(a: P.Poly, lo: F, hi: F) -> list[P.Poly]:
    """Factors of ``a`` (degree <= 4, square-free) vanishing in ``(lo, hi)``.

    Uses the quadratic discriminant test; degrees 3 and 4 with no rational root
    are treated as irreducible, which is sound.
    """
    d = P.p_deg(a)
    if d <= 1:
        return [a] if d == 1 else []
    if d == 2:
        return _split_quadratic(a)
    if d == 3:
        r = P.rational_roots(a)
        if not r:
            return [a]
        out = [(-rr, F(1)) for rr in r]
        rest = a
        for rr in r:
            _, rem = P.p_divmod(rest, (-rr, F(1)))
            if not rem:
                rest = P.p_squarefree(P.p_divmod(rest, (-rr, F(1)))[0])
        if P.p_deg(rest) >= 1:
            out.append(rest)
        return out
    return [a]


def _split_quadratic(a: P.Poly) -> list[P.Poly]:
    """Split a quadratic over ``Q`` by testing the discriminant for a square."""
    ints = P._clear_denominators(a)
    c0, c1, c2 = ints
    disc = c1 * c1 - 4 * c2 * c0
    if disc < 0:
        return [a]
    s = P._isqrt(disc)
    if s * s != disc:
        return [a]
    g1 = (c1 - s, F(2) * c2)
    g2 = (c1 + s, F(2) * c2)
    q1, r1 = P.p_divmod(P.p_monic(a), P.p_monic(g1))
    if r1 == () and P.p_deg(q1) >= 1:
        q2, r2 = P.p_divmod(P.p_monic(a), P.p_monic(g2))
        if r2 == () and P.p_deg(q2) >= 1:
            return [P.p_monic(g1), P.p_monic(g2)]
    return [a]


def _rotate(c: Conic) -> Conic:
    """Apply the axis swap ``(x, y) -> (y, x)`` to a conic.

    A hyperbola whose transverse axis is horizontal has no ``y^2`` term, which
    makes the ``y``-elimination degenerate; a hyperbola whose axis is vertical
    has no ``x^2`` term, which does the same to the ``x``-elimination.  Since
    every real conic has *some* direction in which it is non-degenerate, the two
    eliminations together cover all cases, and the axis swap makes the
    fallbacks genuinely different rather than equivalent.
    """
    A, B, C, D, E, G = c
    return (C, B, A, E, D, G)


def _shear(c: Conic, k: F) -> Conic:
    """Apply the shear ``(x, y) -> (x + k y, y)`` to a conic.

    Two conics can both lack a ``y^2`` term, in which case the ``y``-elimination
    is trivial; similarly for ``x^2``.  A pair can defeat both at once, so a shear
    is used as a last resort: it adds a ``y^2`` term in general and restores a
    usable eliminant.  A shear is injective, so intersection points and their
    distinctness are preserved; only the recovered coordinates must be mapped
    back, via the inverse ``x = X + k Y``.
    """
    A, B, C, D, E, G = c
    # substituting x -> x + k y:
    #   A (x+ky)^2 + B (x+ky) y + C y^2 + D (x+ky) + E y + G
    A2 = A
    B2 = 2 * A * k + B
    C2 = A * k * k + B * k + C
    D2 = D
    E2 = k * D + E
    return _normalise((A2, B2, C2, D2, E2, G))


def _unshear(c: Conic, k: F) -> Conic:
    """Inverse of :func:`_shear`: the shear with parameter ``-k``."""
    return _shear(c, -k)


def _eliminate_any(c1: Conic, c2: Conic):
    """Try the available eliminations, returning a uniform description.

    Eliminating ``y`` needs both conics to have a ``y^2`` term; eliminating ``x``
    needs both to have an ``x^2`` term.  A pair can defeat both at once (e.g. one
    conic missing ``y^2`` while the other misses ``x^2``), so we also try the
    axis swap and a small set of shears, which preserve incidence.

    Returns ``(Pn, Qn, H, mode)`` with ``mode`` in

    * ``elim_y``, ``elim_x``, ``elim_y_swap``, ``elim_x_swap``;
    * ``elim_y_shear<k>``, ``elim_x_shear<k>`` -- as above after ``(x,y) -> (x+k y, y)``,
      so the first coordinate is recovered as ``x = X - k y``.

    The mode is kept explicit (rather than folded into a boolean) because
    mixing the two conventions up produces spurious points.
    """
    for mode in ("elim_y", "elim_x", "elim_y_swap", "elim_x_swap"):
        swap = mode.endswith("_swap")
        d1, d2 = (_rotate(c1), _rotate(c2)) if swap else (c1, c2)
        elim = _eliminate(d1, d2) if mode.startswith("elim_y") else \
            _eliminate_x(d1, d2)
        if elim is not None:
            return (*elim, mode)
    for k in (F(1), F(-1), F(2), F(-2), F(1, 2), F(-1, 2), F(1, 3), F(-1, 3)):
        d1, d2 = _shear(c1, k), _shear(c2, k)
        elim = _eliminate(d1, d2)
        if elim is not None:
            return (*elim, f"elim_y_shear{k}")
        elim = _eliminate_x(d1, d2)
        if elim is not None:
            return (*elim, f"elim_x_shear{k}")
    return None


def _parse_mode(mode: str, c1: Conic, c2: Conic):
    """Decode an elimination mode into the transformed conics and conventions.

    Returns ``(d1, d2, solve_for_y, swapped, shear)`` where ``shear`` is the
    parameter ``k`` of an applied shear (0 when none was applied), so the caller
    can map recovered coordinates back to the original frame.
    """
    swapped = mode.endswith("_swap")
    d1, d2 = (_rotate(c1), _rotate(c2)) if swapped else (c1, c2)
    solve_for_y = mode.startswith("elim_y")
    shear = F(0)
    if "_shear" in mode:
        # the mode encodes the Fraction's repr, e.g. "elim_y_shear1/2"
        token = mode.split("_shear", 1)[1]
        shear = F(token)
        d1, d2 = _shear(c1, shear), _shear(c2, shear)
        swapped = False  # already folded into the shear
    return d1, d2, solve_for_y, swapped, shear


def conic_intersections(c1: Conic, c2: Conic) -> list[AlgebraicPoint]:
    """All real intersection points of two conics, exactly and without repeats."""
    if c1 == c2:
        return []
    got = _eliminate_any(c1, c2)
    if got is None:
        return []
    Pn, Qn, H, mode = got
    d1, d2, solve_for_y, swapped, shear = _parse_mode(mode, c1, c2)
    pts: list[AlgebraicPoint] = []
    for lo, hi in P.real_roots(H):
        for p in _lift(lo, hi, H, Pn, Qn, d1, d2, solve_for_y=solve_for_y,
                       swapped=swapped):
            if shear != 0:
                p = p.unshear(shear)
            pts.append(p)
    return pts


def conic_intersection_count(c1: Conic, c2: Conic) -> int:
    """Number of *distinct* real intersection points of two conics, exactly.

    Each distinct real root of the square-free eliminant contributes one point,
    except at a root where the pencil degenerates: there the two conics share a
    fibre and contribute ``deg gcd`` points instead.  Because ``Qn`` is linear,
    such a root is rational, so this case is always resolved in exact rational
    arithmetic.
    """
    if c1 == c2:
        return 0
    got = _eliminate_any(c1, c2)
    if got is None:
        return 0
    Pn, Qn, H, mode = got
    d1, d2, solve_for_y, _swapped, _shear = _parse_mode(mode, c1, c2)
    sf = P.p_squarefree(H)
    x0 = _qn_zero(Qn)
    degenerate_at = (x0 is not None and P.p_eval(H, x0) == 0)
    total = 0
    for lo, hi in P.real_roots(sf):
        if degenerate_at and lo <= x0 <= hi:
            total += _fibre_count(d1, d2, x0, solve_for_y=solve_for_y)
        else:
            total += 1
    return total


def _vertical_fibre_count(c1: Conic, c2: Conic, x0: F) -> int:
    """Deprecated alias kept for clarity; use :func:`_fibre_count`."""
    return _fibre_count(c1, c2, x0, solve_for_y=True)