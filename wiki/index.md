# Index

## Core integral

- [Integral I(s, alpha, beta)](pages/Integral.md) — the parametric Bessel integral at the center of the project; its definition, substitutions, and quadrature scheme.

## Numerics

- [IntegralNumeric.py](pages/IntegralNumeric.md) — direct numeric evaluation of Z(s, beta) over a grid via `scipy.integrate.quad`, surface-plotted for several alpha.
- [StableIntegrand](pages/StableIntegrand.md) — algebraic rewrite removing the spurious beta -> 0 divergence, hardened `quad`, and twin full/clipped surface plots.
- [Gauss-Laguerre](pages/GaussLaguerre.md) — exact `p=sqrt(t/s)` substitution production integrator; removes the small-s convergence failure down to a calibrated reliability floor.

## Investigations

- [PeakInvestigation](pages/PeakInvestigation.md) — diagnoses the sharp peak at s_min in `IntegralNumeric.py`'s grid as a boundary-cutoff/convergence artifact, not a real feature of Z.
- [Verification](pages/Verification.md) — validates `quad` against the exact alpha->0 closed form: correct away from s~0, genuinely non-convergent (not just imprecise) as s->0.

## Resolved questions

- Does the small-s peak persist once the beta -> 0 divergence is removed? **Yes**, for s near the grid's lower bound — confirmed by [Verification](pages/Verification.md) (100% of `quad` evaluations flagged non-converged at s=1e-6) and fixed separately by [Gauss-Laguerre](pages/GaussLaguerre.md)'s exact substitution.
