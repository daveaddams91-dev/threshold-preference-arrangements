# v1.0.0

First public release. Exact chamber counts for tolerance-based distance orderings
of point sets — arrangements of indifference hyperbolas.

## Research question

Given `n` candidate points and a tolerance `delta > 0`, a vantage point may
declare two candidates tied whenever their distances differ by at most `delta`.
The number of distinguishable tolerant rankings is the number of chambers of the
arrangement

```
H_ij = { X : | |X - Q_i| - |X - Q_j| | = delta },    1 <= i < j <= n
```

of `C(n,2)` hyperbolas. **What is the maximum number of such chambers?**

Zaslavsky posed this as Research Problem 10 in *Abstract arrangements in the
plane* (Discrete Comput. Geom. **27** (2002) 303–351, arXiv:1001.4435). We did
not find a prior result establishing the answer; see the novelty audit below.

## Main theorem

For `n` generic candidates the maximum number of chambers is

```
T(n) = 1 + 2 C(n,2) + 4 C(C(n,2), 2) = 1 + n^2 (n-1)^2 / 2
```

attained for every sufficiently small `delta`, and an upper bound for every
`delta` and every configuration.

| n | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| tolerant `T(n)` | 3 | 19 | 73 | 201 | 451 | 883 |
| classical `G(n)` | 2 | 6 | 18 | 46 | 101 | 197 |
| ratio | 1.50 | 3.17 | 4.06 | 4.37 | 4.47 | 4.48 |

Three ingredients, each proved in the paper:

1. **Chamber formula.** For `m` proper arcs in general position the chamber count
   is `1 + #arcs + #crossings`. Here each hyperbola has two branches, so
   `R = 1 + 2m + I`.
2. **Sharp upper bound.** Two conics meet in at most four points (Bézout), so
   `R <= 1 + 2m + 4 C(m,2) = 1 + 2m^2 = 1 + n^2 (n-1)^2 / 2`.
3. **Saturation.** For small `delta` each hyperbola is a *pair* of branches
   straddling the perpendicular bisector; two such doubled lines crossing
   transversally give four transverse crossings, so every pair of curves meets
   in exactly four points and the Bézout ceiling is achieved.

A fourth result: `T(n) / G(n) -> 4` exactly, where `G(n)` is the classical
Good–Tideman count. Both are `Theta(n^4)` — tolerating ties buys a constant
factor of four in distinguishable rankings, not a change of growth rate.

## Computational contribution

`tolprefs` computes `R(Q, delta)` **exactly**, over the rationals, with no
floating-point arithmetic in the counting path. Two conics are reduced to a
quartic by pencil elimination, real roots are isolated by Sturm bisection with
rational endpoints, and coordinates are recovered inside number fields. Two
degeneracies needed care and are handled exactly: a degenerate pencil (two
conics sharing a vertical fibre, which yields two points with one `x`) and a
conic missing a quadratic term (handled by axis swap and shear).

The core needs only the Python standard library.

## Validation — three independent checks, all passing

1. **Metric verification + independent search.** For 300 random curve pairs,
   every returned point satisfies both defining conditions
   `|d(X,Q_i) - d(X,Q_j)| = delta` to within `1e-6`, and a multistart Newton
   search run on the *metric equations themselves* — never on the conic
   equations — finds nothing the exact code missed.
   **Result: 0 spurious, 0 missed, worst residual 1.2e-11.**
2. **SymPy cross-check** of exact root counts against `resultant` and Groebner
   machinery. This exposed a methodological trap worth recording: counting real
   roots of a resultant is *not* the same as counting distinct real intersection
   points, because the resultant counts `x`-values with algebraic multiplicity,
   so a tangency and a common vertical fibre are each counted twice.
3. **Grid chamber counts.** The plane is rasterised and adjacent cells joined
   unless some `|d_i - d_j| - delta` changes sign between their centres;
   connected components are compared with the exact count. The sign profile was
   constant on every chamber in every case (0 components with a mixed profile),
   profiles never exceeded chambers, and where the grid resolved all chambers
   the two agreed exactly.

**40 tests pass**, exact wherever the result is checkable by hand.

## Reproduce

```bash
git clone https://github.com/rajveersinh-is-dev/threshold-preference-arrangements
cd threshold-preference-arrangements
pip install -e ".[experiments]"
python experiments/run_all.py
```

Stages: unit tests, exact-vs-numerical validation, region counts for
`n = 2..6`, closed-form check, grid chamber counts, figures. The full run takes
about 20 minutes; `python -m pytest tests -q` takes about 2.

Requirements: Python >= 3.10. `numpy`, `scipy`, `matplotlib` only for the
validation and figure scripts.

## Known limitations

Stated in full in the README and in §8 of the paper. In brief:

- **General position is required** for `R = 1 + 2m + I`. `arrangement_stats`
  reports the number of curve pairs meeting fewer than four times and the
  minimum separation between computed points, so a caller can detect a
  violation, but it does not certify it.
- **`delta = |Q_i Q_j|` is a genuine edge case.** At that single value the conic
  equation gives the whole focal line, whereas the metric locus is only the two
  outer rays. We require `delta < min |Q_i Q_j|` throughout.
- **The chamber-to-ranking map is verified computationally, not proved.** In
  every case tested distinct chambers gave distinct rankings, but this is not
  automatic from the definition and we did not prove it.
- **The rational-root search is skipped** above coefficient size `10^5`. Sound
  (Sturm bisection still locates every root exactly) but some rational
  coordinates are held internally as isolating intervals.
- **The `10^-11` residual is a float diagnostic, not a certificate.** Exactness
  of the *counts* rests on exact rational arithmetic and Sturm theory; exactness
  of the reported *coordinates* of irrational intersection points rests on the
  number-field reconstruction, which the metric check only samples numerically.
- **`delta_0` in the saturation theorem is not explicit.** The proof is
  qualitative; an explicit bound in terms of the minimum pairwise distance would
  make the edge of the plateau computable.

## Novelty audit

Searches of OpenAlex, Semantic Scholar, Crossref and the arXiv API for
"indifference hyperbola arrangement", "tolerance distance ordering",
"second-level sophisticated voter", "hyperbolic arrangement region count", and
combinations with Zaslavsky and with Good and Tideman returned the classical
bisector-arrangement literature and recent permutation-counting work (Carbonero,
Castellano, Gordon, Kulick, Ohlinger and Schmitz, arXiv:2106.14140), none of
which establishes `T(n)`, the chamber formula, or the fourfold comparison.

**We did not find a prior result establishing the maximum number of chambers of
an arrangement of indifference hyperbolas.** We make no claim of priority beyond
that: the search was not exhaustive and we did not search non-English
mathematical literature.

**This release is not peer reviewed and is not published in a journal.** It is a
preprint-quality result with a proof and reproducible computations, and it should
be treated as such.

## Contents

```
paper/main.pdf     10-page write-up with proofs
paper/main.tex     source
src/tolprefs/      poly.py, conics.py, arrangement.py
tests/             40 tests
experiments/       validation, region counts, grid check, figures, run_all
figures/           arrangement_n3.png, arrangement_n4.png, scaling.png, delta_sweep.png
results/           machine-readable output
docs/              research notes: what was tried, what failed, why claims are as stated
examples/          worked example for four candidates
```