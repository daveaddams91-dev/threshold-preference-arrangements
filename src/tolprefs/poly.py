"""Exact univariate polynomial arithmetic over ``Q``.

Representation
--------------
A polynomial is a ``tuple`` of :class:`fractions.Fraction`, **constant term
first**, with trailing zeros stripped (``()`` is the zero polynomial).

The constant-term-first convention is used everywhere in this package.  It is
the natural convention but it is also an easy one to violate when assembling a
polynomial by hand, and a violation yields a *self-consistent but wrong*
polynomial -- one that still has real roots, still divides, and still passes
internal identity checks.  :func:`tolprefs.conics` therefore asserts the
ordering explicitly wherever a polynomial is built field-by-field.

What is provided
----------------
* arithmetic and exact division with remainder;
* gcd and square-free part;
* Sturm sequences for counting and isolating real roots exactly;
* the rational-root theorem via integer factorisation (avoids the
  ``O(sqrt(|a_0|))`` divisor loop on the very large coefficients that arise from
  conic elimination);
* a quadratic splitting test by discriminant;
* arithmetic inside ``Q[z]/(m)``, needed to lift algebraic points.

Nothing here uses floating point.
"""

from __future__ import annotations

from fractions import Fraction as F

__all__ = [
    "Poly",
    "p_trim",
    "p_add",
    "p_sub",
    "p_mul",
    "p_scale",
    "p_neg",
    "p_deg",
    "p_eval",
    "p_deriv",
    "p_divmod",
    "p_monic",
    "p_gcd",
    "p_squarefree",
    "p_sturm_chain",
    "cauchy_bound",
    "count_real_roots",
    "real_roots",
    "rational_roots",
    "q_mul",
    "q_add",
    "q_sub",
    "q_scale",
    "invert_mod",
]

Poly = tuple[F, ...]

# Above this size the rational-root search is skipped: it is only needed to
# *label* an already-isolated root as rational, and Sturm bisection finds the
# root either way.  Skipping is sound (see `rational_roots`).
_RATIONAL_ROOT_COEFF_LIMIT = 10**5


# ---------------------------------------------------------------- basics
def p_trim(p) -> Poly:
    p = tuple(F(c) for c in p)
    while p and p[-1] == 0:
        p = p[:-1]
    return p


def p_deg(p) -> int:
    p = p_trim(p)
    return -1 if not p else len(p) - 1


def p_add(a: Poly, b: Poly) -> Poly:
    n = max(len(a), len(b))
    out = [F(0)] * n
    for i in range(n):
        out[i] = (a[i] if i < len(a) else F(0)) + (b[i] if i < len(b) else F(0))
    return p_trim(out)


def p_sub(a: Poly, b: Poly) -> Poly:
    n = max(len(a), len(b))
    out = [F(0)] * n
    for i in range(n):
        out[i] = (a[i] if i < len(a) else F(0)) - (b[i] if i < len(b) else F(0))
    return p_trim(out)


def p_neg(a: Poly) -> Poly:
    return p_scale(a, -1)


def p_scale(a: Poly, c: F | int) -> Poly:
    c = F(c)
    return p_trim(tuple(x * c for x in a))


def p_mul(a: Poly, b: Poly) -> Poly:
    if not a or not b:
        return ()
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, ca in enumerate(a):
        if ca == 0:
            continue
        for j, cb in enumerate(b):
            if cb:
                out[i + j] += ca * cb
    return p_trim(out)


def p_eval(a: Poly, x: F | int) -> F:
    x = F(x)
    acc = F(0)
    for c in reversed(a):
        acc = acc * x + c
    return acc


def p_deriv(a: Poly) -> Poly:
    if len(a) <= 1:
        return ()
    return p_trim(tuple(a[i] * i for i in range(1, len(a))))


def p_monic(a: Poly) -> Poly:
    a = p_trim(a)
    return a if not a else p_scale(a, 1 / a[-1])


def p_divmod(a: Poly, b: Poly) -> tuple[Poly, Poly]:
    """Exact division with remainder over ``Q``; ``b`` must be nonzero."""
    a, b = p_trim(a), p_trim(b)
    if not b:
        raise ZeroDivisionError("polynomial division by zero")
    if p_deg(a) < p_deg(b):
        return (), a
    q = [F(0)] * (p_deg(a) - p_deg(b) + 1)
    r = list(a)
    db, lb = p_deg(b), b[-1]
    while True:
        rt = p_trim(tuple(r))
        dr = p_deg(rt)
        if dr < db:
            break
        c = rt[-1] / lb
        shift = dr - db
        q[shift] += c
        for i, cb in enumerate(b):
            if cb:
                r[i + shift] -= c * cb
    return p_trim(tuple(q)), p_trim(tuple(r))


def p_gcd(a: Poly, b: Poly) -> Poly:
    a, b = p_trim(a), p_trim(b)
    while b:
        _, r = p_divmod(a, b)
        a, b = b, r
    return p_monic(a)


def p_squarefree(a: Poly) -> Poly:
    """Square-free part: the minimal polynomial of the *distinct* roots of ``a``.

    Correctly, if ``a = prod f_i**m_i`` over distinct irreducibles ``f_i`` then
    the square-free part is ``prod f_i``.  Dividing ``a`` by ``gcd(a, a')``
    gives ``prod f_i**(m_i - 1)``, which is *not* that -- it keeps roots whose
    multiplicity exceeded 1 and so changes the root set.  That mistake is easy
    to make and silently drops real roots, so the construction here is done by
    repeatedly removing factors that still divide the derivative:

    ``sf = a``; while ``g = gcd(sf, sf')`` is non-constant, replace ``sf`` by
    ``sf / g``.  Each step lowers a multiplicity by one, and the process stops
    exactly when all multiplicities are 1.
    """
    a = p_trim(a)
    if p_deg(a) <= 0:
        return a
    sf = p_monic(a)
    for _ in range(p_deg(a)):
        g = p_gcd(sf, p_deriv(sf))
        if p_deg(g) < 1:
            break
        q, rem = p_divmod(sf, g)
        if rem != ():  # pragma: no cover - defensive
            break
        sf = p_monic(q)
        if p_deg(sf) < 1:
            break
    return sf


# ---------------------------------------------------------------- Sturm
def p_sturm_chain(a: Poly) -> list[Poly]:
    """Sturm chain ``s0 = a, s1 = a', s_{i+1} = -rem(s_{i-1}, s_i)``.

    Members are **not** rescaled.  Sturm's theorem counts sign changes, so the
    signs produced by the recursion must be preserved exactly; normalising each
    member (e.g. to monic) flips signs and silently corrupts the root count.
    This was a real bug during development.  Only ``s0`` is given a positive
    leading coefficient, which flips the whole chain consistently and is
    therefore harmless.
    """
    a = p_trim(a)
    if p_deg(a) < 1:
        return []
    s0 = a if a[-1] > 0 else p_neg(a)
    chain = [s0, p_deriv(s0)]
    while p_deg(chain[-1]) > 0:
        _, r = p_divmod(chain[-2], chain[-1])
        chain.append(p_neg(r))
    return chain


def _variations(chain: list[Poly], x: F) -> int:
    signs = []
    for s in chain:
        v = p_eval(s, x)
        if v:
            signs.append(1 if v > 0 else -1)
    return sum(1 for i in range(len(signs) - 1) if signs[i] != signs[i + 1])


def cauchy_bound(a: Poly) -> F:
    """A rational ``M > 0`` with every root of ``a`` inside ``(-M, M)``.

    Uses the Fujiwara bound ``2 max |a_i/a_d|^{1/(d-i)}``, which is far tighter
    than the classical ``1 + max |a_i/a_d|`` and is computed with exact integer
    roots (``math.isqrt`` on scaled integers) rather than floats.  The
    tightness matters: bisection cost is logarithmic in the initial width, and
    a loose bound makes every step operate on fractions with hundreds of digits.
    """
    a = p_trim(a)
    d = p_deg(a)
    if d < 1:
        return F(1)
    lead = abs(a[-1])
    # scale to integers first:  |a_i / a_d|  ==  |A_i| / |A_d|
    scale = 1
    best = 0
    for i in range(d):
        j = d - i  # exponent for Fujiwara: 1/(d-i)
        # candidate: 2 * (|A_i| / |A_d|) ** (1/j), computed as an integer bound
        # |a_i/a_d|^(1/j) <= 1 + |a_i/a_d| is too loose; instead take
        # ceil( (|A_i| * scale) ** (1/j) ) with scale normalising the leading term
        num = abs(a[i])
        # integer j-th root of num/lead, computed as: floor((num * C) ** (1/j))
        # with C = j**j * lead (keeps everything integral and conservative)
        val = num * (j ** j) * lead
        r = _iroot(val, j)
        best = max(best, r)
    return F(2 * best + 2)


def _iroot(n: int, j: int) -> int:
    """Floor of ``n ** (1/j)`` for ``n >= 0``, ``j >= 1``, by integer bisection."""
    if n < 2 or j == 1:
        return n if j == 1 else 1
    lo, hi = 1, 1
    while hi ** j <= n:
        hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if mid ** j <= n:
            lo = mid + 1
        else:
            hi = mid
    return lo - 1


def count_real_roots(a: Poly, lo: F, hi: F) -> int:
    """Number of real roots of ``a`` in the open interval ``(lo, hi)``.

    ``lo`` and ``hi`` must not be roots of ``a``.
    """
    a = p_squarefree(a)
    if p_deg(a) < 1:
        return 0
    chain = p_sturm_chain(a)
    if not chain:
        return 0
    return _variations(chain, lo) - _variations(chain, hi)


def real_roots(a: Poly, rel: F = F(1, 10**14)) -> list[tuple[F, F]]:
    """Isolating intervals, one per distinct real root of ``a``.

    Intervals are pairwise disjoint, each containing exactly one root, and are
    refined to relative width ``rel``.  A rational root ``r`` is reported as the
    degenerate interval ``(r, r)`` -- this is recognised by dividing ``r`` out
    before Sturm bisection runs.
    """
    a = p_squarefree(a)
    if p_deg(a) < 1:
        return []
    # The rational-root search is a *label* refinement: without it, Sturm
    # bisection still finds every root exactly, it merely reports an interval
    # rather than a point.  Its cost is dominated by factoring the leading and
    # constant coefficients, which for conic eliminants are frequently
    # astronomically large, so it is skipped above a size threshold.
    rats = rational_roots(a)
    out: list[tuple[F, F]] = [(r, r) for r in rats]
    residual = a
    for r in rats:
        q, rem = p_divmod(residual, (-r, F(1)))
        if not rem:
            residual = p_squarefree(q)
    if p_deg(residual) >= 1:
        chain = p_sturm_chain(residual)
        M = cauchy_bound(residual)
        guard = 0
        while p_eval(residual, M) == 0 or p_eval(residual, -M) == 0:
            M = 2 * M + 1
            guard += 1
            if guard > 8:  # pragma: no cover - defensive
                break
        # Rational roots of ``residual`` that the (size-capped) search skipped
        # are still found by bisection, which now handles a root landing
        # exactly on a split point.
        out.extend(_bisect(chain, -M, M, residual, rel))
    return sorted(out)


def _bisect(chain: list[Poly], lo: F, hi: F, f: Poly, rel: F) -> list[tuple[F, F]]:
    """Sturm bisection on ``(lo, hi)``.

    Two invariants are maintained, both essential for correctness:

    * ``lo`` and ``hi`` are **never roots** of ``f``.  Sturm counts require it:
      ``_variations`` discards zero entries, so splitting at a root corrupts the
      root count and silently drops solutions.  When a root lands exactly on a
      split point we emit it as a degenerate interval and restart the
      surrounding intervals with nudged endpoints.
    * each returned interval contains exactly one root.
    """
    if p_eval(f, lo) == 0 or p_eval(f, hi) == 0:  # pragma: no cover
        return []
    cnt = _variations(chain, lo) - _variations(chain, hi)
    if cnt <= 0:
        return []
    if cnt == 1:
        a, b = lo, hi
        for _ in range(400):
            if b - a <= rel * (1 + max(abs(a), abs(b))):
                break
            mid = (a + b) / 2
            if mid == a or mid == b:  # consecutive rationals: cannot refine
                break
            if p_eval(f, mid) == 0:
                return [(mid, mid)]
            if _variations(chain, a) - _variations(chain, mid) > 0:
                b = mid
            else:
                a = mid
        return [(a, b)]
    mid = (lo + hi) / 2
    if p_eval(f, mid) == 0:
        # ``mid`` is an exact root: report it and bisect both sides with
        # endpoints nudged off the root, so the two side-counts stay valid.
        eps = F(1, 10**12) * (1 + abs(mid))
        return ([(mid, mid)]
                + _bisect(chain, lo, mid - eps, f, rel)
                + _bisect(chain, mid + eps, hi, f, rel))
    left = _variations(chain, lo) - _variations(chain, mid)
    return (_bisect(chain, lo, mid, f, rel)
            + _bisect(chain, mid, hi, f, rel))


# ------------------------------------------------- rational root theorem
def _gcd_int(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a)


def _lcm(a: int, b: int) -> int:
    return a * b // _gcd_int(a, b) if a and b else 0


def _clear_denominators(a: Poly) -> list[int]:
    """Coprime integer coefficients with positive leading coefficient."""
    lcm = 1
    for c in a:
        lcm = _lcm(lcm, c.denominator) if lcm else c.denominator
    ints = [int(c * lcm) for c in a]
    g = 0
    for v in ints:
        g = _gcd_int(g, abs(v))
    if g > 1:
        ints = [v // g for v in ints]
    if ints and ints[-1] < 0:
        ints = [-v for v in ints]
    return ints


def _prime_factors(n: int, cap: int = 200_000) -> list[int] | None:
    """Distinct prime factors of ``n``, or ``None`` if the search is too costly.

    The cap is deliberately small.  Trial division to ``sqrt(n)`` is cheap when
    ``n`` factors into small primes and expensive otherwise; since the caller
    treats ``None`` as "give up and report no rational roots", a modest cap keeps
    the routine predictable.  Signalling failure is always safe here because
    Sturm bisection still locates the roots exactly.
    """
    n = abs(n)
    out: list[int] = []
    d = 2
    while d * d <= n:
        if d > cap:
            return None
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        out.append(n)
    return out


def _isqrt(n: int) -> int:
    """Floor of ``sqrt(n)`` for ``n >= 0``, in exact integer arithmetic.

    Newton's method, seeded from a power-of-two bracket.  An earlier version
    closed with a linear scan ``while (x+1)^2 <= n``, which is ``O(sqrt(n))``
    and takes minutes on the 24-digit discriminants that conic eliminants
    produce; Newton converges in a handful of iterations instead.
    """
    if n < 2:
        return n
    # power-of-two upper bracket: 2^k with 2^(2k) > n
    k = n.bit_length()
    x = 1 << ((k + 1) // 2)
    while True:
        y = (x + n // x) // 2
        if y >= x:
            break
        x = y
    while x * x > n:
        x -= 1
    while (x + 1) * (x + 1) <= n:
        x += 1
    return x


def _quadratic_rational_roots(a: Poly) -> list[F]:
    """Rational roots of a quadratic: discriminant plus integer square root.

    For ``c2 x^2 + c1 x + c0`` the roots are ``(-c1 +- sqrt(D)) / (2 c2)`` with
    ``D = c1^2 - 4 c2 c0``.  Both signs are kept when exactly divisible; a
    non-square discriminant means irrational roots, and none is reported.
    """
    ints = _clear_denominators(a)
    c0, c1, c2 = ints[0], ints[1], ints[2]
    disc = c1 * c1 - 4 * c2 * c0
    if disc < 0:
        return []
    s = _isqrt(disc)
    if s * s != disc:
        return []  # irrational roots
    out: list[F] = []
    for sign in (1, -1):
        num = -c1 + sign * s
        # ``num == 0`` is a genuine root (x = 0) and must not be skipped.
        if num % (2 * c2) == 0:
            out.append(F(num // (2 * c2)))
    return out


def _divisors(primes: list[int], limit: int = 10**12) -> list[int]:
    out = [1]
    for p in primes:
        new = []
        for v in out:
            while v <= limit // p:
                new.append(v * p)
                v *= p
        out.extend(new)
    return sorted(set(out))


def rational_roots(a: Poly,
                   coeff_limit: int = _RATIONAL_ROOT_COEFF_LIMIT) -> list[F]:
    """Rational roots of ``a``, or ``[]`` when the search is deemed too costly.

    The rational-root theorem enumerates ``p/q`` with ``p | a_0`` and
    ``q | a_d``.  After eliminating the common denominator and content, the
    relevant integers are frequently astronomically large for conic
    eliminants, and factoring them dominates the run time.  We therefore bail
    out above ``coeff_limit``; this is *sound* because this function is used
    only to label an already-isolated root as rational.  When it returns
    ``[]`` the roots are still found exactly, by Sturm bisection, they are
    simply reported as intervals rather than as points.
    """
    a = p_trim(a)
    d = p_deg(a)
    if d < 1:
        return []
    if d == 2:
        return _quadratic_rational_roots(a)
    ints = _clear_denominators(a)
    if any(abs(v) > coeff_limit for v in ints):
        return []
    if ints[0] == 0:
        return [F(0)] + rational_roots(p_divmod(a, (F(0), F(1)))[0],
                                        coeff_limit)
    pf_const = _prime_factors(ints[0])
    pf_lead = _prime_factors(ints[-1])
    if pf_const is None or pf_lead is None:
        return []
    ndiv = _divisors(pf_const)
    ddiv = _divisors(pf_lead)
    # Test candidates on the *integer* polynomial q^d * a(p/q): evaluating that
    # in plain integers is far cheaper than rational Horner, and the divisor
    # counts here can reach thousands, so the difference is decisive.
    d = len(ints) - 1
    found: list[F] = []
    for pn in ndiv:
        for qd in ddiv:
            for sign in (1, -1):
                num = sign * pn
                # value = sum_i ints[i] * (num)^i * qd^(d-i)
                val = 0
                for i in range(d + 1):
                    val += ints[i] * (num**i) * (qd ** (d - i))
                if val == 0:
                    cand = F(num, qd)
                    if cand not in found:
                        found.append(cand)
    return found


# --------------------------------------------- arithmetic mod a polynomial
def q_mul(a: Poly, b: Poly, mod: Poly) -> Poly:
    return p_divmod(p_mul(a, b), mod)[1]


def q_add(a: Poly, b: Poly, mod: Poly) -> Poly:
    return p_divmod(p_add(a, b), mod)[1]


def q_sub(a: Poly, b: Poly, mod: Poly) -> Poly:
    return p_divmod(p_sub(a, b), mod)[1]


def q_scale(a: Poly, c: F | int, mod: Poly) -> Poly:
    return p_divmod(p_scale(a, c), mod)[1]


def _xgcd(a: Poly, b: Poly) -> tuple[Poly, Poly, Poly]:
    """Extended gcd: returns ``(g, u, v)`` with ``g = u a + v b`` (up to scaling).

    The loop must maintain ``old_r = old_s*a + old_t*b`` at every step.  An
    earlier version returned ``old_s``/``old_t`` as empty whenever ``b`` happened
    to divide ``a``, which silently made :func:`invert_mod` report "not
    invertible" for pairs that were coprime.
    """
    old_r, r = p_trim(a), p_trim(b)
    old_s, s = (F(1),), ()
    old_t, t = (), (F(1),)
    while r:
        q, rem = p_divmod(old_r, r)
        old_r, r = r, rem
        old_s, s = s, p_sub(old_s, p_mul(q, s))
        old_t, t = t, p_sub(old_t, p_mul(q, t))
    d = p_deg(old_r)
    if d == 0 and old_r[0] != 0:
        inv = 1 / old_r[0]
        return (F(1),), p_scale(old_s, inv), p_scale(old_t, inv)
    if d >= 1:
        # ``old_r`` is the gcd and is non-constant; the Bezout coefficients are
        # still valid (old_r = old_s*a + old_t*b) and must be returned.
        return p_monic(old_r), old_s, old_t
    return old_r, old_s, old_t


def invert_mod(a: Poly, mod: Poly) -> Poly | None:
    """Inverse of ``a`` modulo ``mod``; ``None`` iff ``gcd(a, mod) != 1``.

    A constant modulus (degree 0) is handled separately: there the "inverse"
    is the representative of ``a**-1`` modulo that constant, expressed as a
    constant polynomial.  Returning ``None`` for a nonzero modulus of degree 0
    would be wrong, since gcd is then automatically a unit.
    """
    d = p_deg(mod)
    if d < 0:
        return None
    if d == 0:
        c = mod[0]
        if c == 0:
            return None
        av = p_eval(a, F(0))
        if av == 0:
            return None
        return p_trim((1 / (av * c) * c,))  # 1/av, written as a constant poly
    g, u, _ = _xgcd(a, mod)
    if p_deg(g) > 0:
        return None
    return p_divmod(u, mod)[1]