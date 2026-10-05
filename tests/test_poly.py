"""Exact univariate polynomial arithmetic: correctness tests.

The results checked here are all decidable by hand, so the tests are exact
rather than approximate.
"""
import sys
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tolprefs import poly as P  # noqa: E402


def test_add_sub_mul():
    a, b = (F(-6), F(5), F(-1)), (F(1), F(1))
    assert P.p_mul(a, b) == (F(-6), F(-1), F(4), F(-1))
    assert P.p_sub(P.p_add(a, b), b) == a
    assert P.p_scale(a, 3) == tuple(3 * c for c in a)


def test_divmod_exact():
    x1, x2 = (F(0), F(1)), (F(1), F(1))
    q, r = P.p_divmod(P.p_mul(x1, x2), x2)
    assert r == () and q == x1
    # (x-1)^2 (x-2)(x-3) = x^4 - 7x^3 + 17x^2 - 17x + 6
    H = (F(6), F(-17), F(17), F(-7), F(1))
    q, r = P.p_divmod(H, (F(-1), F(1)))  # divide by (x-1)
    assert r == () and P.p_deg(q) == 3


def test_gcd_and_derivative():
    assert P.p_gcd((F(1), F(2), F(1)), (F(1), F(0), F(1))) == (F(1),)
    assert P.p_gcd((F(-2), F(3)), (F(-1), F(1))) == (F(1),)
    assert P.p_deriv((F(6), F(-11), F(1))) == (F(-11), F(2))


def test_gcd_is_bezout():
    """u*a + v*b must equal the returned gcd; a regression guard for _xgcd."""
    a = (F(2), F(1))          # 2x + 1
    b = (F(1), F(1))          # x + 1
    g, u, v = P._xgcd(a, b)
    assert P.p_deg(g) == 0
    lhs = P.p_add(P.p_mul(u, a), P.p_mul(v, b))
    assert P.p_trim(lhs) == P.p_monic(g), (u, v, g)


def test_xgcd_b_divides_a():
    """Regression: when b divides a the Bezout coefficients must not be empty."""
    a = (F(4), F(2), F(1))    # x^2 + 2x + 4 = (x+1)^2 + 3
    b = (F(-1), F(1))          # x - 1
    g, u, v = P._xgcd(a, b)
    lhs = P.p_add(P.p_mul(u, a), P.p_mul(v, b))
    assert P.p_trim(lhs) == P.p_monic(g), (u, v, g)
    # a case where b divides a exactly: x^2 - 1 = (x-1)(x+1)
    a2 = (F(-1), F(0), F(1))
    b2 = (F(-1), F(1))
    g2, u2, v2 = P._xgcd(a2, b2)
    assert P.p_deg(g2) >= 1
    lhs2 = P.p_add(P.p_mul(u2, a2), P.p_mul(v2, b2))
    assert P.p_trim(lhs2) == P.p_monic(g2)


def test_invert_mod():
    assert P.invert_mod((F(2), F(1)), (F(1), F(1))) == (F(1),)     # 1/(2x+1) mod (x+1)
    assert P.invert_mod((F(1), F(1)), (F(1), F(0), F(1))) == (F(1, 2), F(-1, 2))
    assert P.invert_mod((F(2),), (F(4),)) == (F(1, 2),)           # 2^-1 mod 4
    assert P.invert_mod((F(0),), (F(5),)) is None


def test_sturm_chain_not_rescaled():
    """Members must keep their recursion signs; renormalising flips the count."""
    a = (F(2), F(0), F(1))    # x^2 + 2 has no real roots
    chain = P.p_sturm_chain(a)
    assert len(chain) == 3
    assert P._variations(chain, F(7)) - P._variations(chain, F(-7)) == 0
    b = (F(-2), F(0), F(1))    # x^2 - 2 has two real roots
    chain_b = P.p_sturm_chain(b)
    assert P._variations(chain_b, F(-7)) - P._variations(chain_b, F(7)) == 2


def test_count_real_roots():
    x2m2 = (F(-2), F(0), F(1))            # x^2 - 2 : 2 real roots
    assert P.count_real_roots(x2m2, F(-2), F(2)) == 2
    x2p2 = (F(2), F(0), F(1))             # x^2 + 2 : none
    assert P.count_real_roots(x2p2, F(-9), F(9)) == 0


def test_real_roots_positions():
    # (x-1)^2 (x-2)(x-3): distinct real roots 1, 2, 3
    H = (F(6), F(-17), F(17), F(-7), F(1))
    found = sorted(float(lo) for lo, _ in P.real_roots(H))
    assert [round(v, 9) for v in found] == [1.0, 2.0, 3.0]


def test_real_roots_handle_root_on_split_point():
    """A rational root landing exactly on a bisection point must not lose others.

    This is a regression: Sturm's variation count discards zero entries, so
    splitting at a root used to corrupt the count and drop solutions.
    """
    sf = (F(-13396319, 15480000), F(43865519, 15480000), F(-38291, 12900), F(1))
    assert P.p_deg(sf) == 3
    assert P.p_eval(sf, F(1)) == 0        # exact root sitting on a split point
    roots = sorted(round(float(lo), 6) for lo, _ in P.real_roots(sf))
    assert len(roots) == 3, roots
    for expect in (0.662977, 1.0, 1.305318):
        assert any(abs(r - expect) < 1e-4 for r in roots), (expect, roots)


def test_squarefree_keeps_distinct_roots():
    """A double root must collapse to one, not be dropped entirely."""
    # (x-1)^2 (x-2)(x-3)
    H = (F(6), F(-17), F(17), F(-7), F(1))
    sf = P.p_squarefree(H)
    assert P.p_deg(sf) == 3
    assert sorted(round(float(lo), 9) for lo, _ in P.real_roots(sf)) == [1.0, 2.0, 3.0]


def test_squarefree_triple_root():
    # (x-1)^3 (x+2)^2 = x^5 + x^4 - 5x^3 - x^2 + 8x - 4
    H = (F(-4), F(8), F(-1), F(-5), F(1), F(1))
    sf = P.p_squarefree(H)
    assert P.p_deg(sf) == 2
    assert sorted(round(float(lo), 9) for lo, _ in P.real_roots(sf)) == [-2.0, 1.0]


def test_isqrt():
    for n in (0, 1, 2, 3, 10**12, 123456789012 ** 2, 10**30):
        s = P._isqrt(n)
        assert s * s <= n < (s + 1) * (s + 1), n


def test_quadratic_rational_roots():
    # x^2 - 3x + 2 = (x-1)(x-2)
    assert sorted(P.rational_roots((F(2), F(-3), F(1)))) == [F(1), F(2)]
    # x^2 - 5x + 2 has irrational roots
    assert P.rational_roots((F(2), F(-5), F(1))) == []
    # x^2 - x = x(x-1)
    assert sorted(P.rational_roots((F(0), F(-1), F(1)))) == [F(0), F(1)]


def test_cauchy_bound_contains_roots():
    for poly in ((F(6), F(-17), F(17), F(-7), F(1)),
                 (F(-2), F(0), F(1)),
                 (F(1), F(0), F(0), F(0), F(1))):
        M = P.cauchy_bound(poly)
        assert all(abs(lo) <= M and abs(hi) <= M
                   for lo, hi in P.real_roots(poly))


def test_field_arithmetic():
    # in Q[x]/(x^2 - 2): x * x = 2
    mod = (F(-2), F(0), F(1))
    prod = P.q_mul((F(0), F(1)), (F(0), F(1)), mod)
    assert P.p_trim(prod) == (F(2),)