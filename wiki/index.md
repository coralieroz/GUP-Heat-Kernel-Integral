# Index

## Core integral

- [Integral I(s, alpha, beta)](pages/Integral.md) — the parametric Bessel integral at the center of the project; its definition, substitutions, and quadrature scheme.

## Numerics

- [IntegralAnalytic.py](pages/IntegralAnalytic.md) — builds I(s, alpha, beta) as a symbolic Gauss–Laguerre sum so it supports `diff`/`subs`; validated against `scipy.integrate.quad`.
- [IntegralNumeric.py](pages/IntegralNumeric.md) — direct numeric evaluation of Z(s, beta) over a grid via `scipy.integrate.quad`, surface-plotted for several alpha.
- [StableIntegrand.md](pages/StableIntegrand.md) — algebraic rewrite of the integrand removing the spurious beta -> 0 divergence; basis for `plot_surface_stable_integrand.py`.
- [BesselFunctionTest.py](pages/BesselFunctionTest.md) — sanity plot of J0/J1, no project-specific content.

## Investigations

- [PeakInvestigation](pages/PeakInvestigation.md) — diagnoses the sharp peak at s_min in `IntegralNumeric.py`'s grid as a boundary-cutoff/convergence artifact, not a real feature of Z. **This is the motivation for the `boundary-cutoff/removal` branch.**

## Open questions

- Is the small-s peak fully resolved by the stable-integrand rewrite, or does it persist for some (alpha, beta) once the beta -> 0 divergence is removed? See [PeakInvestigation](pages/PeakInvestigation.md) and [StableIntegrand](pages/StableIntegrand.md).
