# GUP Heat-Kernel Integral

Numerical evaluation and visualisation of a parametric integral arising from
a Generalized Uncertainty Principle (GUP) heat-kernel calculation in the
accompanying capstone report.

## The integral

The central object is a one-sided integral over `p in [0, inf)`:

```
I(s, alpha, beta) = int_0^inf  p^2 (1 + (alpha p)^2) / beta
                     * exp( -s p^2 (1 + 2(alpha p)^2/3) )
                     * J_1( beta p / (1 + (alpha p)^2) )  dp
```

where `J_1` is the order-1 Bessel function of the first kind, `s >= 0` is a
damping/regularisation parameter, `alpha` is a GUP deformation parameter, and
`beta = 2 r sin(omega/2)` is a momentum-transfer-like variable. The goal is to
visualise the surface `Z(s, beta) = I(s, alpha, beta)` for fixed, small
`alpha`.

A hard constraint throughout this project: the **numerical** evaluation of
this integral must use no series expansions or limiting approximations — that
algebra is done separately, analytically, in the report. Code-level
evaluation must work with the exact integral.

## The journey

### 1. First attempt — direct quadrature

`IntegralNumeric.py` evaluates `I` directly with `scipy.integrate.quad` over
`p in [0, inf)`, gridded over `(s, beta)` for a few small `alpha`. A sharp,
resolution-sensitive peak appeared at the smallest `s` in the grid.

### 2. Diagnosing the peak

`PeakInvestigation.py` tests whether that peak is real:

- **Lower-bound test** — move the grid's `s_min` away from 0. The peak
  collapses and moves off the boundary.
- **Resolution test** — recompute on the same domain at different grid
  resolutions. The peak's value (even its sign) changes wildly.

Both point to the same conclusion: the peak is a **boundary-cutoff /
non-convergence artifact**, not a physical feature of `I`. At very small `s`,
the damping `exp(-s p^2 ...)` barely suppresses the oscillatory `J_1` tail, so
`quad`'s adaptive subdivision does not actually converge there.

### 3. A second, separate bug

The prefactor `p^2(1+q)/beta` (with `q = (alpha p)^2`) divides by `beta`
directly and diverges as `beta -> 0` — a genuine bug, independent of the
convergence issue above. Using `A2 = beta p/(1+q)`, the algebraically
identical rewrite

```
p^2(1+q)/beta * J_1(A2)  =  p^3 * J_1(A2)/A2
```

is **entire**: `J_1(A2)/A2 -> 1/2` as `A2 -> 0`, removing the singularity with
no expansion or limit involved — just an exact algebraic identity. This fix is
applied throughout (`plot_surface_stable_integrand.py`,
`IntegralNumeric.py`, `PeakInvestigation.py`, `gauss_laguerre_integral.py`).

### 4. Verification against an exact benchmark

`verify_integration.py` checks `quad` against an exact closed form available
in the `alpha -> 0` limit (Weber's second exponential integral):

```
Z(s, beta) -> (1/(4 s^2)) * exp(-beta^2 / 4s)        as alpha -> 0
```

which matches Expression 2.14 of the capstone report (d=4 leading term) up to
a constant factor. Result: `quad` agrees with this benchmark to better than
`1e-8` relative error away from `s ~ 0`, and `Z(s=1, beta=0) = 0.25` exactly
as expected — the surface is **not** "approximately zero" away from small
`s`; it only looked that way because the genuine `1/s^2` divergence near
`s = 0` dominates an autoscaled z-axis. At `s = 1e-6`, however, **100% of
`quad` evaluations are flagged non-converged** — confirming the failure mode
identified in step 2 is real, not just imprecise.

### 5. Visualisation

`plot_surface_stable_integrand.py` hardens `quad` (explicit tolerances,
`full_output=1` non-convergence flagging, a closed-form regression check) and
produces **two** views of the same surface, since a single z-axis cannot show
both the genuine divergence and the underlying smooth structure at once:

- the raw surface, full autoscaled z-axis (shows the real divergence and the
  non-converged points, separately counted);
- a convergence-masked, percentile-clipped surface, with non-converged and
  outlier points set to NaN *before* plotting (avoiding the boundary-wall
  artifact that `matplotlib`'s `set_zlim` alone produces), revealing the
  smooth Gaussian ridge underneath, overlaid against the closed-form
  benchmark.

### 6. A production integrator with no convergence failure

`gauss_laguerre_integral.py` substitutes `p = sqrt(t/s)`, which folds
`exp(-s p^2)` exactly into the Gauss-Laguerre quadrature weight `exp(-t)` —
no expansion, no limit, the same integral exactly rewritten. This removes the
specific failure mode from step 2/4 (the damping is no longer separate from
the quadrature weight) down to a calibrated reliability floor,
`S_MIN_RELIABLE = 0.02`, below which the weighted sum suffers float64
catastrophic cancellation. That floor was confirmed to be a genuine open
numerical-analysis problem — adding more nodes does not help, and even
arbitrary-precision (`mpmath`) quadrature fails the same way there, because
the real difficulty is the highly oscillatory, weakly-damped tail itself, not
floating-point precision. `plot_surface_gauss_laguerre.py` cross-validates
this integrator against both hardened `quad` and the closed form (machine
precision agreement) and produces the same twin full/clipped plots.

## Illustrative plots

Raw surface — genuine `1/s^2` divergence as `s -> 0`, with non-converged
points still included:

![Raw surface](TestingApparentIntegralDivergences/full_range_raw.png)

Convergence-masked, clipped view of the same grid — the underlying Gaussian
ridge, with the closed-form benchmark overlaid in red:

![Clipped surface](TestingApparentIntegralDivergences/clipped_masked.png)

The exact-substitution Gauss-Laguerre integrator over its reliable domain —
clean and noise-free, with no non-convergence masking required:

![Gauss-Laguerre surface](TestingApparentIntegralDivergences/gauss_laguerre_clipped.png)

## Repository layout

| Path | Role |
|---|---|
| `IntegralNumeric.py` | First direct-quadrature evaluation; surfaced the small-`s` peak. |
| `PeakInvestigation.py` | Diagnoses the peak as a boundary/convergence artifact. |
| `verify_integration.py` | Validates `quad` against the exact `alpha -> 0` closed form. |
| `plot_surface_stable_integrand.py` | Hardened `quad` + beta=0 fix + twin full/clipped plots. |
| `gauss_laguerre_integral.py` | Exact-substitution Gauss-Laguerre production integrator. |
| `plot_surface_gauss_laguerre.py` | Cross-validation and surface plots for the integrator above. |
| `TestingApparentIntegralDivergences/` | Output plots from the scripts above. |
| `wiki/` | Page-per-topic notes: definitions, findings, and how everything connects. |

## Reproduce

```
pip install -r requirements.txt
python verify_integration.py
python plot_surface_stable_integrand.py
python plot_surface_gauss_laguerre.py
```

Each plotting script prints its diagnostics to the console and saves PNGs
into `TestingApparentIntegralDivergences/`.

## Further detail

See `wiki/` for page-level notes on the integral itself, each script, and how
the findings connect — start at `wiki/index.md`.
