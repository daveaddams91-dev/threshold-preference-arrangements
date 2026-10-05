# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-10-05

First public release.

### Added

- **`tolprefs.poly`** — exact univariate polynomial arithmetic over `Q`,
  constant-term-first representation. Sturm sequences (sign-preserving), tight
  real-root isolation by bisection, square-free part by repeated gcd removal,
  rational roots via integer factorisation with a discriminant fast path for
  quadratics, Newton integer square root, Fujiwara root bound, and arithmetic
  inside `Q[x]/(m)`.
- **`tolprefs.conics`** — the indifference curve
  `|d(X,Q_i) - d(X,Q_j)| = delta` as an exact rational conic; pencil
  elimination with axis-swap and shear fallbacks; exact count of *distinct* real
  intersection points, including the degenerate shared-fibre case; exact
  lifting of intersection points into number fields.
- **`tolprefs.arrangement`** — chamber counts for the arrangement of
  indifference hyperbolas, the Bézout ceiling
  `1 + 2 C(n,2) + 4 C(C(n,2), 2) = 1 + n^2 (n-1)^2 / 2`, and the classical
  Good–Tideman baseline.
- **Main theorem and proofs** (`paper/main.tex`): the chamber formula in general
  position; Bézout saturation for all sufficiently small `delta`; the exact
  maximum `T(n) = 1 + n^2 (n-1)^2 / 2`; and the fourfold comparison against the
  classical count.
- **40 tests** (`tests/`), exact wherever the result is checkable by hand.
- **Three independent validations** (`experiments/`):
  metric verification plus an independent multistart Newton search; a SymPy
  cross-check of exact root counts; and a grid flood-fill comparison of chamber
  counts.
- **Figures** (`figures/`), each regenerated from the exact machinery rather than
  from stored coordinates.
- **Reproducible pipeline** (`experiments/run_all.py`).

### Known limitations

- The chamber formula requires general position; `arrangement_stats` reports
  the data needed to detect a violation but does not certify it.
- The rational-root search is skipped above coefficient size `10^5`. Sound, but
  some rational intersection coordinates are held internally as isolating
  intervals.
- At `delta = |Q_i Q_j|` the conic equation over-approximates the metric locus
  (whole focal line versus the two outer rays); `delta < min |Q_i Q_j|` is
  required throughout.
- Injectivity of the chamber-to-ranking map is computational, not proved.
- Not peer reviewed.

### Fixed during development

Seven substantive defects in the exact core, each found by an independent check
and each now guarded by a regression test: a reversed coefficient order in the
`y`-split; hand-built polynomials assembled leading-coefficient-first; a
square-free part computed as `a/gcd(a,a')`; Sturm bisection splitting at an
exact root; an extended-gcd routine returning empty Bézout coefficients when one
input divides the other; a sign error in the quadratic formula together with the
root `0` being skipped; and an integer square root closing with a linear scan.