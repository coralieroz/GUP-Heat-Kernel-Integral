# Integral I(s, alpha, beta)

**Source:** concept, implemented in `IntegralNumeric.py`, `plot_surface_stable_integrand.py`, `gauss_laguerre_integral.py`
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

- **Direct quadrature** (`IntegralNumeric.py`, `PeakInvestigation.py`,
  `plot_surface_stable_integrand.py`): evaluate `I` directly via
  `scipy.integrate.quad` over the original `p` domain, gridded over `(s, beta)`
  for fixed `alpha`. Correct away from `s ~ 0` (see [Verification](Verification.md))
  but genuinely fails to converge as `s -> 0` (see [Peak Investigation](PeakInvestigation.md)).
- **Gauss-Laguerre substitution** (`gauss_laguerre_integral.py`): substitute
  `p = sqrt(t/s)` to fold `exp(-s p^2)` into the Laguerre weight, turning the
  integral into a fast-converging quadrature sum that removes the small-`s`
  convergence failure (down to a calibrated floor). See
  [Gauss-Laguerre](GaussLaguerre.md).

## Validation benchmark

For `alpha -> 0`, Weber's second exponential integral gives an exact closed
form `Z(s, beta) -> (1/(4 s^2)) exp(-beta^2/4s)`, matching Expression 2.14 of
the capstone report (d=4 leading term) up to the constant `4 pi^2`. Used only
to validate the two numeric implementations above — never as a production
integrator, since it drops the alpha-dependent GUP terms. See
[Verification](Verification.md).

## Known numerical hazard

The prefactor `p^2(1+q)/beta` (with `q = (alpha p)^2`) diverges as `beta -> 0`.
See [Stable Integrand](StableIntegrand.md) for the algebraic fix.

The integral is only **conditionally convergent** when `s` is very small,
since the damping `exp(-s p^2 ...)` barely suppresses the oscillatory `J_1`
tail. See [Peak Investigation](PeakInvestigation.md).

## Related pages

[IntegralNumeric.py](IntegralNumeric.md), [Stable Integrand](StableIntegrand.md), [Gauss-Laguerre](GaussLaguerre.md), [Verification](Verification.md), [Peak Investigation](PeakInvestigation.md)
