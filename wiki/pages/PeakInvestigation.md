# Peak Investigation

**Source:** PeakInvestigation.py, PeakInvestigation.png
**Status:** active — direct motivation for the `boundary-cutoff/removal` branch

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

## Conclusion (implied by the script's framing)

The peak is a boundary-cutoff artifact of evaluating the grid down to
`s_min = 1e-6` with an under-resolved, conditionally-convergent quadrature —
not a genuine feature of `I(s, alpha, beta)`. This is why the active branch
(`boundary-cutoff/removal`) removes/avoids that boundary regime rather than
treating the peak as physical. See [Stable Integrand](StableIntegrand.md) for
the follow-up that also fixes the separate `beta -> 0` divergence in the
integrand itself.

## Related pages

[Integral](Integral.md), [IntegralNumeric.py](IntegralNumeric.md), [Stable Integrand](StableIntegrand.md)
