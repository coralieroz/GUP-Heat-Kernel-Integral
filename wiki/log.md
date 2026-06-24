# Log

- 2026-06-24 — Wiki initialized (per the LLM-wiki pattern from
  https://gist.github.com/karpathy/442a6bf555914893e9891c11519de94f). Ingested
  all root-level scripts on the `boundary-cutoff/removal` branch:
  `IntegralAnalytic.py`, `IntegralNumeric.py`, `PeakInvestigation.py`,
  `plot_surface_stable_integrand.py`, `BesselFunctionTest.py`. Removed scripts
  `Integral.py` / `IntegralII.py` (deleted on this branch, superseded by the
  above) were not ingested.

- 2026-06-24 — Verified the user's reported integration concerns against the
  exact `alpha -> 0` closed form (`verify_integration.py`): `quad` is correct
  away from `s ~ 0` and genuinely non-convergent (not merely imprecise) as
  `s -> 0`; the "surface looks ~ 0" report was a z-axis autoscaling artifact
  masking a real `1/s^2` divergence. Fixed the `beta -> 0` division bug in
  `IntegralNumeric.py` and `PeakInvestigation.py` (same `p^3*J_1(A2)/A2`
  rewrite already used in the stable script) and an axis-transpose meshgrid
  bug in `IntegralNumeric.py`. Hardened `plot_surface_stable_integrand.py`'s
  quadrature (explicit tolerances, `full_output=1` non-convergence flag,
  closed-form regression check) and added twin full-range/clipped surface
  plots (NaN-masked before `plot_surface` to avoid boundary-wall artifacts).
  Built an exact `p=sqrt(t/s)` Gauss-Laguerre production integrator
  (`gauss_laguerre_integral.py`) with no expansions or limits, calibrated a
  reliability floor (`S_MIN_RELIABLE = 0.02`, confirmed to be a genuine
  oscillatory-tail numerical-analysis limit rather than a fixable bug), and
  cross-validated it against both hardened `quad` and the closed form to
  machine precision. New pages: [Verification](pages/Verification.md),
  [GaussLaguerre](pages/GaussLaguerre.md); updated
  [Integral](pages/Integral.md), [IntegralNumeric](pages/IntegralNumeric.md),
  [PeakInvestigation](pages/PeakInvestigation.md),
  [StableIntegrand](pages/StableIntegrand.md).

- 2026-06-24 — Release cleanup for GitHub upload. Removed `IntegralAnalytic.py`
  (superseded by `gauss_laguerre_integral.py`, which no longer references it)
  and the scratch `BesselFunctionTest.py`, plus their wiki pages. Renamed
  `StableIntegrandPlots/` to `TestingApparentIntegralDivergences/` and removed
  the stale per-alpha plots it held from the old `IntegralNumeric.py`
  workflow, keeping only the four current result plots. Added
  `requirements.txt` and a narrative `README.md`.
