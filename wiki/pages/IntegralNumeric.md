# IntegralNumeric.py

**Source:** IntegralNumeric.py
**Status:** superseded by [Stable Integrand](StableIntegrand.md) for plotting; integrand definition still canonical

Direct numeric evaluation of [I(s, alpha, beta)](Integral.md) (referred to in
this script as `Z`) via `scipy.integrate.quad` on the raw integrand, gridded
over `s` and `beta` for three small values of `alpha`, surface-plotted with
`matplotlib`.

## Grid

- `s = linspace(0.1, 100, 100)` — `s_min = 0.1` is chosen deliberately: below
  it, `exp(-s p^2)` barely damps the oscillatory `J_1` term and
  `quad(limit=200)` stops converging.
- `beta = linspace(-100, 100, 100)`.
- `alpha in [1e-15, 1e-9, 1e-3]`.

## Known issue

This was the script in which the sharp peak near small `s` was first
observed, motivating [Peak Investigation](PeakInvestigation.md). The peak
turned out to be a quadrature/boundary artifact rather than a real feature of
`Z`, which is part of why this branch is named `boundary-cutoff/removal`.

## Related pages

[Integral](Integral.md), [Peak Investigation](PeakInvestigation.md), [Stable Integrand](StableIntegrand.md)
