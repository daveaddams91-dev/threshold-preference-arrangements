# Tolerance-based distance orderings of point sets

Exact chamber counts for arrangements of indifference hyperbolas.

Given `n` candidate points and a tolerance `delta > 0`, a vantage point may
declare two candidates tied whenever their distances differ by at most `delta`.
The number of distinguishable tolerant rankings is the number of chambers of
the arrangement

```
H_ij = { X : | |X - Q_i| - |X - Q_j| | = delta },      1 <= i < j <= n
```

of `C(n,2)` hyperbolas. This package computes that number **exactly**, over the
rationals, with no floating-point arithmetic in the counting path.

## The result

For `n` generic candidates the maximum number of chambers is

```
T(n) = 1 + 2 C(n,2) + 4 C(C(n,2), 2) = 1 + n^2 (n-1)^2 / 2
```

attained for every sufficiently small `delta`, and an upper bound for every
`delta` and every configuration. Values:

| n | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| tolerant `T(n)` | 3 | 19 | 73 | 201 | 451 | 883 |
| classical `G(n)` | 2 | 6 | 18 | 46 | 101 | 197 |
| ratio | 1.50 | 3.17 | 4.06 | 4.37 | 4.47 | 4.48 |

The classical row is the Good–Tideman count for perpendicular bisectors,
`sum_{i<=2} |s(n, n-i)|`. The ratio tends to **4**: tolerating ties buys a
constant factor four in distinguishable rankings, not a change of growth rate
(both are `Theta(n^4)`).

Zaslavsky posed this as Research Problem 10 in *Abstract arrangements in the
plane* (Discrete Comput. Geom. **27** (2002) 303–351, arXiv:1001.4435). We did
not find a prior result establishing the formula. See `paper/main.pdf`.

## Install

```bash
git clone https://github.com/daveaddams91-dev/threshold-preference-arrangements
cd threshold-preference-arrangements
pip install -e .              # core needs only the standard library
pip install -e ".[experiments]"   # adds numpy, scipy, matplotlib
```

Python >= 3.10.

## Use

```python
from fractions import Fraction
from tolprefs import arrangement, semicircle_bound

pts = [(Fraction(0), Fraction(0)),
       (Fraction(5), Fraction(1)),
       (Fraction(1), Fraction(6)),
       (Fraction(6), Fraction(4))]

stats = arrangement.arrangement_stats(pts, Fraction(1, 20))
print(stats.region_count)              # 73
print(semicircle_bound(4))             # 73   <- Bezout ceiling, attained
```

Single pair, exact intersection points:

```python
from tolprefs import conics

c1 = conics.gamma_conic(pts[0], pts[1], Fraction(1, 20))
c2 = conics.gamma_conic(pts[0], pts[2], Fraction(1, 20))
print(conics.conic_intersection_count(c1, c2))   # 4, exactly
for p in conics.conic_intersections(c1, c2):
    print(p.value())                              # (x, y) as floats
```

## Reproduce

```bash
python -m pytest tests -q                          # 40 tests
python experiments/validate_intersections.py       # exact vs independent numerical
python experiments/region_counts.py --trials 12    # exact region counts, n = 2..6
python experiments/grid_check.py                   # independent grid chamber counts
python experiments/make_figures.py                 # figures/
```

Or all at once: `python experiments/run_all.py`.

## Validation

Three independent checks, all currently passing:

1. **Metric verification + independent search** — for 300 random curve pairs,
   every returned point satisfies both defining conditions
   `|d(X,Q_i) - d(X,Q_j)| = delta` to within `1e-6`, and a multistart Newton
   search run on the *metric equations themselves* finds no point the exact
   code missed. Result: `0` spurious, `0` missed, worst residual `1.2e-11`.
2. **SymPy cross-check** — counts compared against SymPy `resultant` and
   Groebner machinery. This exposed a real methodological trap: counting real
   roots of a resultant is *not* the same as counting distinct real
   intersection points, because the resultant counts `x`-values with
   algebraic multiplicity.
3. **Grid chamber counts** — the plane is rasterised, adjacent cells joined
   unless some `|d_i - d_j| - delta` changes sign between their centres, and
   connected components are compared with the exact count. The sign profile is
   constant on every chamber in every case tested.

## Repository layout

```
src/tolprefs/poly.py         exact Q[x]: Sturm chains, root isolation,
                             square-free part, rational roots, Q[x]/(m)
src/tolprefs/conics.py       the indifference curve as a conic; exact
                             intersection counts and point lifting
src/tolprefs/arrangement.py  chamber counts, Bezout ceiling, classical baseline
tests/                       40 tests, exact where checkable by hand
experiments/                 validation, region counts, grid check, figures
paper/main.tex               full write-up with proofs
results/                     machine-readable output
figures/                     PNG figures
```

## Honest limitations

- **General position is required** for `R = 1 + 2m + I`. `arrangement_stats`
  reports the number of curve pairs meeting fewer than four times and the
  minimum separation between computed points, so a caller can detect a
  violation, but it does not certify general position.
- **The rational-root search is skipped** when cleared coefficients exceed
  `10^5`. Sound — Sturm bisection still locates every root exactly, it merely
  reports an interval rather than a point — but some rational intersection
  coordinates are held internally as intervals.
- **`delta = |Q_i Q_j|` is a genuine edge case.** At that single value the
  conic equation gives the whole focal line whereas the metric locus is only
  the two outer rays. We require `delta < min |Q_i Q_j|` throughout.
- **The metric residual is a float diagnostic, not a certificate.** Exactness of
  the *counts* rests on exact rational arithmetic and Sturm theory; exactness of
  the reported *coordinates* of irrational intersection points rests on the
  number-field reconstruction, which the metric check only samples numerically.
- **Injectivity of chamber → ranking is computational, not proved.**
- **Not peer reviewed.** No claim of publication.

## Licence

MIT. See `LICENSE`.