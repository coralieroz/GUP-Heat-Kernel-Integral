"""
Surface plot of Z(s, beta) using the exact-substitution Gauss-Laguerre
evaluator (gauss_laguerre_integral.Z_laguerre) instead of scipy.integrate.quad
applied directly to the raw oscillatory integral.

Reuses, rather than re-implementing:
  - plot_surface_stable_integrand.quad_safe           (hardened adaptive quadrature, for cross-validation)
  - plot_surface_stable_integrand.closed_form_alpha0  (exact alpha->0 reference)
  - gauss_laguerre_integral.Z_laguerre / S_MIN_RELIABLE  (the production evaluator)

For s < S_MIN_RELIABLE the Laguerre evaluator is not reliable (float64
cancellation -- see gauss_laguerre_integral.py's module docstring for the
calibration that established this); those grid points are left as NaN
(skipped, not filled with noise), and the masked region is reported. This is
documented as this project's known open numerical limitation at very small s,
not silently hidden behind a plausible-looking number.
"""

import numpy as np
import matplotlib.pyplot as plt

from gauss_laguerre_integral import Z_laguerre, S_MIN_RELIABLE
from plot_surface_stable_integrand import quad_safe, closed_form_alpha0, ALPHA


def cross_validate():
    """Compare Laguerre vs hardened quad vs the exact closed form on a spot-check grid."""
    print(f"Cross-validation (alpha={ALPHA:.1e}):")
    hdr = (f"{'s':>8} {'beta':>6}  {'Laguerre':>14} {'quad (hardened)':>16} {'closed form':>14}  "
           f"{'Lag vs CF':>10} {'quad vs CF':>10}")
    print(hdr)
    print("-" * len(hdr))
    for s in np.linspace(S_MIN_RELIABLE, 2.0, 6):
        for beta in (-1.5, 0.0, 1.5):
            lag = Z_laguerre(s, ALPHA, beta)
            q, _, flagged = quad_safe(s, ALPHA, beta)
            cf = closed_form_alpha0(s, beta)
            rel_lag = abs(lag - cf) / abs(cf) if cf else float("nan")
            rel_quad = abs(q - cf) / abs(cf) if cf else float("nan")
            flag_str = "  (quad flagged non-converged)" if flagged else ""
            print(f"{s:>8.4f} {beta:>6.2f}  {lag:>14.6e} {q:>16.6e} {cf:>14.6e}  "
                  f"{rel_lag:>10.2e} {rel_quad:>10.2e}{flag_str}")
    print()


if __name__ == "__main__":
    cross_validate()

    s_vals = np.linspace(1e-6, 2, 60)
    beta_vals = np.linspace(-2, 2, 120)

    Z = np.full((len(s_vals), len(beta_vals)), np.nan)
    n_unreliable_rows = int((s_vals < S_MIN_RELIABLE).sum())
    print(f"{n_unreliable_rows} / {len(s_vals)} s-rows below S_MIN_RELIABLE={S_MIN_RELIABLE:g} "
          f"are masked (known open numerical limitation; see gauss_laguerre_integral.py)")

    for i, s in enumerate(s_vals):
        if s < S_MIN_RELIABLE:
            continue
        for j, beta in enumerate(beta_vals):
            Z[i, j] = Z_laguerre(s, ALPHA, beta)

    S, B = np.meshgrid(s_vals, beta_vals, indexing='ij')

    # --- Figure 1: full range over the reliable domain (genuine divergence as s -> S_MIN_RELIABLE) ---
    fig1 = plt.figure(figsize=(10, 7))
    ax1 = fig1.add_subplot(111, projection='3d')
    surf1 = ax1.plot_surface(S, B, Z, cmap='viridis', linewidth=0, antialiased=True)
    fig1.colorbar(surf1, ax=ax1, shrink=0.6, pad=0.1, label='Z')
    ax1.set_xlabel('s')
    ax1.set_ylabel(r'$\beta$')
    ax1.set_zlabel('Z')
    ax1.set_title(fr'$Z(s,\beta)$ -- Gauss-Laguerre, exact substitution ($\alpha={ALPHA}$)' +
                  f'\n s < {S_MIN_RELIABLE:g} masked (open numerical limitation, not noise)')
    fig1.tight_layout()
    fig1.savefig('TestingApparentIntegralDivergences/gauss_laguerre_full.png', dpi=150)

    # --- Figure 2: clipped to [p1,p99] of the reliable values, with closed-form overlay ---
    finite = np.isfinite(Z)
    zlo, zhi = np.nanpercentile(Z[finite], [1, 99])
    Z_clip = np.clip(Z, zlo, zhi)
    Z_ref_clip = np.clip(closed_form_alpha0(S, B), zlo, zhi)

    fig2 = plt.figure(figsize=(10, 7))
    ax2 = fig2.add_subplot(111, projection='3d')
    surf2 = ax2.plot_surface(S, B, Z_clip, cmap='viridis', vmin=zlo, vmax=zhi, linewidth=0, antialiased=True)
    ax2.plot_wireframe(S, B, Z_ref_clip, color='red', linewidth=0.3, alpha=0.4, rstride=4, cstride=8)
    ax2.set_zlim(zlo, zhi)
    fig2.colorbar(surf2, ax=ax2, shrink=0.6, pad=0.1, label='Z')
    ax2.set_xlabel('s')
    ax2.set_ylabel(r'$\beta$')
    ax2.set_zlabel('Z')
    ax2.set_title(r'$Z(s,\beta)$ -- Gauss-Laguerre, clipped to [p1,p99]' +
                  '\n(red wireframe: exact $\\alpha\\to0$ closed form)')
    fig2.tight_layout()
    fig2.savefig('TestingApparentIntegralDivergences/gauss_laguerre_clipped.png', dpi=150)

    plt.show()
