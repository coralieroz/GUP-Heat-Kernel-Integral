# Stable Integrand

**Source:** plot_surface_stable_integrand.py, TestingApparentIntegralDivergences/
**Status:** active

Algebraic rewrite of the [I(s, alpha, beta)](Integral.md) integrand to remove
a spurious divergence as `beta -> 0`.

## The rewrite

The original integrand has `T1 = p^2(1+q)/beta` (with `q = (alpha p)^2`),
which blows up as `beta -> 0`. But using `A2 = beta p / (1+q)`:

```
T1 = p^2(1+q)/beta = p^3 / A2
```

so `T1 * J_1(A2) = p^3 * J_1(A2)/A2`, which is **entire** — `J_1(A2)/A2 -> 1/2`
as `A2 -> 0`, so the limit is `p^3/2`, no singularity. The script implements
this directly:

```python
ratio = 0.5 if abs(A2) < 1e-8 else j1(A2) / A2
return (p ** 3) * ratio * np.exp(A1)
```

This eliminates the spurious spike in the small-`(s, beta)` corner of the
surface plot that the naive `T1 * exp(A1) * j1(A2)` form produced.

## Hardened quadrature

`quad_safe()` wraps `scipy.integrate.quad` with tight `epsabs`/`epsrel`, a
raised `limit`, and `full_output=1`, reading the non-convergence flag instead
of silently trusting (or blanket-suppressing warnings on) the result. A
`_regression_check()` asserts agreement with the `alpha -> 0` closed form
(see [Verification](Verification.md)) before the script runs its main grid.

## Output — twin plots

For fixed `alpha = 1e-9` over `s in [1e-6, 2]`, `beta in [-2, 2]`, the script
produces **two** figures in `TestingApparentIntegralDivergences/`:

- `full_range_raw.png` — the raw surface, full autoscaled z-axis. This shows
  `Z`'s genuine `1/s^2`-type divergence as `s -> 0`, but also includes the
  non-converged points identified in [Verification](Verification.md) (their
  count is printed and annotated on the plot).
- `clipped_masked.png` — non-converged points (`quad`'s own flag), points
  with `Zerr` above a relative-error threshold, and percentile outliers are
  NaN-masked *before* `plot_surface` (NaN-masking avoids the boundary-wall
  artifact that `set_zlim` alone produces on out-of-range vertices), then the
  z-axis is clipped to `[p1, p99]` of what remains. This reveals the smooth
  Gaussian ridge, overlaid with the closed-form wireframe for visual
  validation.

## Resolved — the open question

The small-`s` peak from [Peak Investigation](PeakInvestigation.md) does still
appear here for some `(alpha, beta)` as the `s`-domain approaches `1e-6` —
confirmed by [Verification](Verification.md) (100% of `quad` evaluations
flagged non-converged at `s=1e-6`). The rewrite alone fixes the `beta -> 0`
divergence but does not fix `quad`'s small-`s` convergence failure; that is
addressed separately by the exact-substitution
[Gauss-Laguerre integrator](GaussLaguerre.md), down to its own calibrated
reliability floor (`S_MIN_RELIABLE = 0.02`) below which the failure mode is a
genuine open numerical-analysis problem rather than an implementation bug.

## Related pages

[Integral](Integral.md), [Peak Investigation](PeakInvestigation.md), [Verification](Verification.md), [Gauss-Laguerre](GaussLaguerre.md)
