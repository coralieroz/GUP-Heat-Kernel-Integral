# Gauss-Laguerre Integrator

**Source:** gauss_laguerre_integral.py, plot_surface_gauss_laguerre.py
**Status:** active — production integrator

Self-contained, beta=0-safe Gauss-Laguerre evaluator for
[I(s, alpha, beta)](Integral.md), built directly on the exact `p = sqrt(t/s)`
substitution. No expansions or limits are used anywhere in this evaluator —
only an exact change of variables plus the same removable-singularity rewrite
used in [Stable Integrand](StableIntegrand.md).

## Method

`p = sqrt(t/s)` folds `exp(-s p^2)` into the Laguerre weight `exp(-t)`, turning
the integral into `sum_i w_i h(t_i)` with

```
h(t) = [t / (2 s^2)] * exp(-2 alpha^2 t^2 / (3s)) * J_1(A2)/A2,
A2 = beta * sqrt(t/s) / (1 + alpha^2 t/s)
```

and `J_1(A2)/A2 -> 1/2` as `A2 -> 0`, so `h(t)` is well-defined at `beta=0` —
unlike the old `IntegralAnalytic.py` symbolic form (removed; its
`poly_piece = (t/s)*(1+a2_t_s)/beta` divided by `beta` directly and blew up
there). Node/weight generation uses `scipy.special.roots_laguerre`.

## Reliability floor

`S_MIN_RELIABLE = 0.02`. Below this, the weighted Laguerre sum suffers float64
catastrophic cancellation: each term's effective Bessel argument
`beta*sqrt(t_i/s)` grows like `1/sqrt(s)`, and the terms become large and
oscillate in sign. This was confirmed **not** to be an under-resolution
problem — `N=40` already matches `N=300` (see calibration table in
`gauss_laguerre_integral.py`'s docstring) — nor a floating-point precision
problem — `mpmath`'s arbitrary-precision `quad` fails the same way for
`s ~ 1e-3..1e-6`. The real difficulty is the highly oscillatory, weakly-damped
tail itself, which would need a dedicated oscillatory-quadrature method
(Levin-type, Filon-type, or `mpmath.quadosc` with an explicit period) — out of
scope here. `Z_laguerre` raises `UnreliableRegionError` below this floor
rather than returning a plausible-looking but meaningless number.

## Cross-validation

`plot_surface_gauss_laguerre.cross_validate()` compares `Z_laguerre` against
hardened `quad` ([Stable Integrand](StableIntegrand.md)'s `quad_safe`) and the
`alpha -> 0` closed form ([Verification](Verification.md)) across
`s in [0.02, 2]`, `beta in {-1.5, 0, 1.5}`: machine-precision agreement
(~1e-16) almost everywhere, rising to ~7e-5 only at the most extreme
exponential-decay magnitudes.

## Output

`plot_surface_gauss_laguerre.py` produces two surfaces over the same
`(s, beta)` grid as [Stable Integrand](StableIntegrand.md)
(`TestingApparentIntegralDivergences/`):
`gauss_laguerre_full.png` (reliable domain only, `s < S_MIN_RELIABLE` masked as
NaN — a documented limitation, not noise) and `gauss_laguerre_clipped.png`
(percentile-clipped, with the closed-form wireframe overlay). Both are clean,
noise-free Gaussian ridges, in contrast to the raw-`quad` plots' spikes near
`s -> 0`.

## Related pages

[Integral](Integral.md), [Stable Integrand](StableIntegrand.md), [Verification](Verification.md)
