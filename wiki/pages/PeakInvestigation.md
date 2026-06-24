# Peak Investigation

**Source:** PeakInvestigation.py, PeakInvestigation.png
**Status:** active — the finding that motivated the [Stable Integrand](StableIntegrand.md) and [Gauss-Laguerre](GaussLaguerre.md) work

Diagnoses a sharp positive peak seen in [IntegralNumeric.py](IntegralNumeric.md)'s
surface plot, sitting at the smallest `s` value of the grid (`s_min = 1e-6`).
At that `s`, `exp(-s p^2)` is almost switched off, so the `p`-integral is
highly oscillatory and only conditionally convergent — `quad()` with a finite
`limit` may not have actually converged there.

## Two diagnostic checks

1. **Lower-bound test** — recompute `Z` with the `s`-grid's lower bound moved
   away from 0 (`s_min = 1e-6` vs `1e-1`, same `beta` range/resolution). If the
   peak height collapses and moves off the `s`-boundary, the peak was an
   artifact of where the grid was cut off, not a property of `Z` itself.
2. **Resolution test** — recompute `Z` on the *same* domain (`s_min = 1e-6`)
   at different grid resolutions (`n = 30` vs `80`). If the peak value/location
   changes wildly (even sign) with resolution, it's unconverged numerical
   noise rather than a smooth feature.

`peak_info()` finds the largest-magnitude `(s, beta, Z)` point and reports
whether it sits on the `s`-boundary (`idx[0]==0`).

## Conclusion

The peak is a boundary-cutoff artifact of evaluating the grid down to
`s_min = 1e-6` with an under-resolved, conditionally-convergent quadrature —
not a genuine feature of `I(s, alpha, beta)`. This script also had the same
`beta -> 0` division bug as `IntegralNumeric.py`, fixed in place with the
`p^3 * J_1(A2)/A2` rewrite from [Stable Integrand](StableIntegrand.md)
(its own `meshgrid` usage was already consistent, via `Z.T`, so needed no
transpose fix).

[Verification](Verification.md) confirms this conclusion directly: re-running
`quad` with `full_output=1` shows 100% of evaluations flagged non-converged at
`s=1e-6`, and the [Gauss-Laguerre](GaussLaguerre.md) substitution — which
removes this specific failure mode — produces a clean, noise-free surface over
the same domain (down to its own calibrated reliability floor).

## Related pages

[Integral](Integral.md), [IntegralNumeric.py](IntegralNumeric.md), [Stable Integrand](StableIntegrand.md), [Verification](Verification.md), [Gauss-Laguerre](GaussLaguerre.md)
