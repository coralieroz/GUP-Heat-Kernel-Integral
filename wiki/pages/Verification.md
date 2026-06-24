# Verification

**Source:** verify_integration.py
**Status:** active — read-only diagnostics, no production code

Decisive check of whether `scipy.integrate.quad` is evaluating
[I(s, alpha, beta)](Integral.md) correctly, using the exact `alpha -> 0`
closed form as a ground truth.

## Closed form (benchmark only, never the integrator)

Weber's second exponential integral gives, for `alpha -> 0`:

```
int_0^inf p^2 J_1(beta p) exp(-s p^2) dp = beta/(4 s^2) exp(-beta^2/4s)
=>  Z(s, beta) -> (1/(4 s^2)) exp(-beta^2/4s)
```

This matches Expression 2.14 of the capstone report (d=4 leading term,
`(4 pi s)^-2 exp(-beta^2/4s)`) up to the constant `4 pi^2`. It is used **only**
to validate the numerics — it drops the alpha-dependent GUP correction, so it
is never used as a production evaluator.

## Checks performed

- **Closed-form check** — `quad` (with `full_output=1`, tight `epsabs`/`epsrel`)
  vs the closed form for `alpha=1e-12` across a `(s, beta)` grid. Worst-case
  relative error in the well-damped region: **1.24e-9**.
- **Bulk-vs-spike magnitude** — `Z(s=1, beta=0) = 0.25`, matching the closed
  form exactly. Confirms the surface is not "approximately zero" away from
  small `s`; it only looked that way once a genuine `1/s^2` divergence near
  `s=0` was on the same z-axis.
- **Convergence audit** — re-running `quad` with `full_output=1` and reading
  the non-convergence flag (see [Integral](Integral.md) for the parsing
  detail) shows **100% of evaluations flagged non-converged at s=1e-6**: the
  damping `exp(-s p^2)` is too weak there to tame the oscillatory `J_1` tail,
  so the adaptive subdivision genuinely fails, rather than returning a
  slightly-off-but-meaningful number.
- **alpha-scaling check** — the GUP correction term scales like
  `alpha^2 beta^4 / (6 s^3)` relative to the leading `1/(4s^2)`, consistent
  with blow-ups worsening for larger alpha and staying mild for alpha <~ 1e-6.

## Conclusion

`quad` is **correct** away from `s ~ 0` (sub-1e-8 relative error against an
exact benchmark) and **genuinely fails** (not just imprecise) at the smallest
`s` values used in this project's grids. Both findings carry directly into
[Stable Integrand](StableIntegrand.md) (hardened `quad`, convergence masking)
and [Gauss-Laguerre](GaussLaguerre.md) (the substitution that removes this
specific failure mode for `s` down to `S_MIN_RELIABLE`).

## Related pages

[Integral](Integral.md), [Stable Integrand](StableIntegrand.md), [Gauss-Laguerre](GaussLaguerre.md), [Peak Investigation](PeakInvestigation.md)
