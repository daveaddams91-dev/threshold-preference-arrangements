# Research notes

Working notes: what was tried, what failed, and why the final claims are stated
the way they are. Kept separate from the paper, which states only what is
proved.

## How the problem was chosen

The brief was to find a genuine open problem, implement exact computations, and
prove something. The web search tool was unavailable for the whole session, so
all literature work went through `webfetch` against public bibliographic APIs:
OpenAlex, Semantic Scholar, Crossref, and the arXiv API.

Candidates considered and rejected:

| Candidate | Why rejected |
| --- | --- |
| Minimum / maximum number of orderings of `n` points in `R^d` | Already known: minimum `2n-2` (Carbonero et al. 2022, arXiv:2106.14140), maximum the Stirling sum (Good & Tideman 1977). Not open. |
| Generic special-position arrangements (Zaslavsky's Research Problems 7–9) | Requires the affine-dependence matroid of the point set; the interesting counts are known for small `n` and the general answer is close to tautological. |
| Counting regions of the bisector arrangement itself | Settled by Orlik–Solomon / Zaslavsky. |

Zaslavsky, *Abstract arrangements in the plane* (Discrete Comput. Geom. **27**
(2002) 303–351, arXiv:1001.4435), §13.3 "Hyperbolic dissections", Research
Problem 10, asks for the maximum number of preference rankings of exactly the
arrangement studied here. That was selected.

## Novelty position

What we can say: searches of OpenAlex, Semantic Scholar, Crossref and arXiv for
the phrases in `paper/main.tex` §"Novelty audit" returned the classical
bisector-arrangement literature and the recent permutation-counting work; none
of it establishes `T(n) = 1 + n^2 (n-1)^2 / 2`, the chamber formula, or the
fourfold comparison.

What we cannot say: that this is the first result of its kind. The search was
not exhaustive, we did not search non-English literature, and we did not search
every citation graph forward from the references above. The paper says exactly
this and no more.

## Why `R = 1 + 2m + I`

Each `H_ij` is a proper embedded copy of `R` with two branches, so
`m` curves give `2m` arcs. For an arrangement of proper arcs in general position,
add arcs one at a time: an arc meeting the existing arrangement in `k >= 1`
transverse points is cut into `k+1` pieces, each splitting one chamber, and an
arc meeting nothing adds one. Each crossing is counted exactly once, by whichever
of its two arcs is added later. Hence `R = 1 + 2m + I`.

Two ways this can fail, both of which *lower* the count:

- **Tangency.** A tangent branch touches the old boundary without separating, so
  it adds one chamber where two were counted.
- **Triple point.** One point in `I` arises from `C(k,2) >= k-1` pairwise
  approaches, each of which would otherwise have been a separate crossing.

Both are exactly the hypotheses of "general position", and both are automatic
for generic `Q`. `arrangement_stats` reports the number of curve pairs meeting
fewer than four times and the minimum separation between computed points, which
is enough to detect most violations but is not a certificate.

## Why Bézout is saturated for small delta

The concrete geometry, which took several attempts to see clearly:

`H_ij` has semi-transverse axis `a = delta/2` and semi-conjugate axis
`b = sqrt(L - delta^2)/2 ~ sqrt(L)/2`. Its asymptotes make an angle
`arctan(a/b) ~ delta/sqrt(L)` with the bisector `ell_ij`. So as `delta -> 0`
the hyperbola is **not** a small perturbation of the line: it is two branches
that *straddle* `ell_ij` and squeeze onto it from both sides, with a
separation of order `delta`.

Two such "doubled lines" that cross transversally give a `2 x 2` grid of
crossings near the line-crossing point, i.e. four transverse crossings. Bézout
caps two conics at four points over `C` counted with multiplicity, so those four
are all of them. Since there are finitely many pairs, one `delta_0` works
simultaneously for all of them.

The transversality argument deserves a note, because the obvious version of it
is wrong. One might try to apply the implicit function theorem at the crossing
point `p` of the two *lines*. That fails: the conic polynomial's gradient at `p`
tends to zero as `delta -> 0`, since the limiting conic is a double line, whose
gradient vanishes on it. So the IFT neighbourhood shrinks with `delta` and sees
only one branch. The working argument is instead topological: each branch runs
from `dB` to `dB`, so a Jordan-arc argument forces a crossing, and continuity of
the tangent directions forces transversality.

## The endpoint `delta = |Q_i Q_j|`

A genuine trap, found by a failing test rather than by reading. At this value
the metric locus is only the two outer rays of the line `Q_iQ_j`, by the
equality case of the reverse triangle inequality. But the conic equation
collapses to `t = 0`, the *whole* line counted twice. So `gamma_conic`
over-approximates the true indifference set there. Everything is stated under
`delta < min |Q_i Q_j|`, and there is a test pinning the behaviour down.

## The asymptotic story we got wrong first

The first guess was that tolerancing ties changes the growth rate from `Theta(n^3)`
to `Theta(n^4)`. That is wrong: the classical planar count is
`1 + C(n,2) + 2 C(n,3) + 3 C(n,4)`, which is already asymptotic to `n^4/8`.
Both are `Theta(n^4)`, and the ratio tends to exactly `4`. The paper states the
correct version, and `tests/test_arrangement.py` asserts the limit numerically
along with the rate at which the ratio approaches it.

## Validation strategy, and why there are three checks

The counting code is exact, which makes it easy to fool oneself: a wrong
elimination still returns plausible integers. Three checks that share no code
with the core:

1. **Metric verification plus independent search.** Check every returned point
   against the *defining* condition `|d_i - d_j| = delta` (not the conic
   equation), and run a multistart Newton search on the metric equations to look
   for anything missed. This is the check with the most power, because it never
   touches the algebra under test.
2. **SymPy cross-check.** Useful, but it exposed a methodological trap rather
   than a bug in our code: counting real roots of a resultant is *not* the same
   as counting distinct real intersection points. The resultant counts
   `x`-values with algebraic multiplicity, so a tangency is counted twice and a
   common vertical fibre is counted twice. SymPy's Groebner bases are not
   square-free either. Several apparent "mismatches" turned out to be the
   *reference* being wrong. Worth remembering.
3. **Grid chamber counts.** Independent of the arrangement machinery: rasterise,
   join adjacent cells unless some `|d_i - d_j| - delta` changes sign between
   their centres, take connected components. Confirms the chamber formula
   outright in the cases where the grid resolves everything.

## Bugs that mattered

Seven defects in the exact core, each found by an independent check. They are
listed in the changelog; the general lesson is worth stating separately.

**A reversed coefficient order produces a self-consistent wrong answer.**
Polynomials are stored constant-term-first. Writing a polynomial by hand
`(lead, lin, const)` instead of `(const, lin, lead)` yields something that
still evaluates consistently, still divides, still has real roots, and still
passes every internal identity check. It was caught only by comparing against
the defining metric condition. The library now asserts the ordering at every
site where a polynomial is assembled field-by-field
(`conics._check_split`, `conics._check_const_first`), and the regressions are
tests.

The other six: a square-free part computed as `a / gcd(a, a')` (which keeps
roots whose multiplicity exceeded one, changing the root set); Sturm bisection
splitting at an exact root, which corrupts the variation count and silently
drops solutions; an extended-gcd routine returning empty Bézout coefficients
when one input divides the other; a sign error in the quadratic formula plus
the root `0` being skipped; an integer square root closing with a linear scan
that took minutes on 24-digit discriminants; and mis-reading the `elim_x` result
as if it were `elim_y`, which produced 27 spurious points before the elimination
mode was made explicit.

## Performance notes

The exact path was at one point 162 seconds for a single curve pair. The causes,
in order of impact:

- `_isqrt` closing with `while (x+1)^2 <= n: x += 1` — `O(sqrt(n))`, minutes on
  the 24-digit discriminants that vertical-fibre gcds produce. Replaced with
  Newton's method seeded from a power-of-two bracket.
- The rational-root search enumerating divisors of enormous integers. Now
  skipped above a coefficient-size threshold, with a discriminant fast path for
  quadratics (the common fibre case).
- A loose root bound. `1 + max |a_i/a_d|` gave intervals hundreds of digits wide;
  the Fujiwara bound `2 max |a_i/a_d|^{1/(d-i)}`, computed with integer roots,
  is far tighter and cut the same case to under two seconds.

None of these affect correctness of the counts — Sturm theory and Bézout do —
only the time to reach them.

## What is still open

Listed in `paper/main.tex` §8. The two we would attack first:

- An **explicit** `delta_0`. Lemma 3.4 is qualitative; an explicit bound in terms
  of the minimum pairwise distance would make the edge of the plateau
  computable and turn the $\delta$-experiments into a theorem.
- **Injectivity** of the chamber-to-ranking map. Currently computational. It is
  not automatic: a path between two chambers could cross each intervening curve
  twice and restore every sign. A proof would upgrade `R` from an upper bound on
  distinguishable rankings to an exact count of them.

## Reproducibility notes

- `python experiments/run_all.py` runs every stage in a separate process and
  returns non-zero if any fails.
- Seeds are fixed in every script; the region-count sweep takes `--seed`.
- `results/` holds machine-readable output; every figure is regenerated from the
  exact machinery rather than from stored coordinates, so a figure cannot
  silently disagree with the code.
- The core needs only the standard library. `numpy`, `scipy` and `matplotlib`
  are needed only by the validation and figure scripts, and `sympy` only by the
  cross-check, which is not part of the shipped counting path.