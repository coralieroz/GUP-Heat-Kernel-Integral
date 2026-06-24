# Integral I(s, alpha, beta)

**Source:** concept, implemented in `IntegralAnalytic.py`, `IntegralNumeric.py`, `plot_surface_stable_integrand.py`
**Status:** active

The project's central object is a one-sided parametric integral over `p in [0, inf)`:

```
I(s, alpha, beta) = int_0^inf  p^2 (1 + (alpha p)^2) / beta
                     * exp( -s p^2 (1 + 2(alpha p)^2/3) )
                     * J_1( beta p / (1 + (alpha p)^2) )  dp
```

where `J_1` is the order-1 Bessel function of the first kind. `s >= 0` is a
damping/regularization parameter, `alpha` a GUP deformation parameter, and
`beta = 2 r sin(omega/2)` relates to a momentum-transfer-like variable
(`r >= 0`, `omega = tau - tau'`).

## Two implementations

- **Numeric** (`IntegralNumeric.py`, `PeakInvestigation.py`,
  `plot_surface_stable_integrand.py`): evaluate `I` directly via
  `scipy.integrate.quad` over the original `p` domain, gridded over `(s, beta)`
  for fixed `alpha`.
- **Analytic/symbolic** (`IntegralAnalytic.py`): substitute `p = sqrt(t/s)` to
  turn the integral into a Gauss-Laguerre quadrature sum, yielding a genuine
  SymPy expression in `(s, alpha, beta)` that supports `diff`/`subs`/`lambdify`
  — useful for getting closed-form-ish derivatives (e.g. `dI/ds`) for the paper.

## Known numerical hazard

The prefactor `p^2(1+q)/beta` (with `q = (alpha p)^2`) diverges as `beta -> 0`.
See [Stable Integrand](StableIntegrand.md) for the algebraic fix.

The integral is only **conditionally convergent** when `s` is very small,
since the damping `exp(-s p^2 ...)` barely suppresses the oscillatory `J_1`
tail. See [Peak Investigation](PeakInvestigation.md).

## Related pages

[IntegralAnalytic.py](IntegralAnalytic.md), [IntegralNumeric.py](IntegralNumeric.md), [Stable Integrand](StableIntegrand.md), [Peak Investigation](PeakInvestigation.md)
