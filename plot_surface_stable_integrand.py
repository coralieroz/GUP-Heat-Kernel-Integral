"""
Stable 3D surface plot of Z(s, beta) using the algebraically rewritten integrand.

The original integrand T1 * exp(A1) * j1(A2) has T1 = p^2*(1+q)/beta, which
diverges as beta -> 0.  Using A2 = beta*p/(1+q) we have p^2*(1+q)/beta = p^3/A2,
so T1 * j1(A2) = p^3 * j1(A2)/A2, which is entire (limit = p^3/2 at A2 = 0).
This eliminates the spurious spike in the small-(s, beta) corner.

Z(s, beta) GENUINELY diverges as s -> 0: at alpha -> 0 it has the closed form
(Weber's second exponential integral)

    Z(s, beta) -> (1/(4 s^2)) * exp(-beta^2 / (4s))

which matches Expression 2.14 of the capstone report up to a constant factor.
That Gaussian-in-beta has width ~ sqrt(s), so for s ~ 1e-6 the true peak is only
~0.002 wide in beta -- far narrower than this script's beta grid spacing -- and
is therefore invisible on any reasonably sized grid.  What quad actually
returns at small s, on a grid that never lands inside that sliver, is instead
dominated by non-convergence: exp(-s p^2) damps too weakly to tame the
oscillatory J_1(A2) tail, and scipy.integrate.quad's adaptive subdivision
cannot resolve it. verify_integration.py confirms this directly: at s=1e-6,
100% of quad evaluations across the beta grid are flagged non-converged, and
the values returned there are numerical noise, not the real (and much smaller
in magnitude, off the exact beta=0 line) signal.

This script therefore renders two views: the raw surface (showing the noise
dominating near s -> 0) and a convergence-masked, percentile-clipped surface
(showing the underlying smooth Gaussian ridge, validated against the closed
form). Non-convergent points are masked rather than silently trusted.
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, IntegrationWarning
from scipy.special import j1

ALPHA = 1e-9

REL_TOL = 1e-3            # relative-error threshold for flagging an unreliable quad() point
OUTLIER_PCT = 99          # percentile backstop: flag |Z| beyond this percentile of all finite |Z|
DISPLAY_PCT_LOW = 1       # lower percentile used for the masked plot's zlim/color limits
DISPLAY_PCT_HIGH = 99     # upper percentile used for the masked plot's zlim/color limits


def integrand(p, s, alpha, beta):
    q = (alpha * p) ** 2
    A1 = -s * (p ** 2) * (1.0 + 2.0 * q / 3.0)
    A2 = (beta * p) / (1.0 + q)
    ratio = 0.5 if abs(A2) < 1e-8 else j1(A2) / A2
    return (p ** 3) * ratio * np.exp(A1)


def closed_form_alpha0(s, beta):
    """
    Exact alpha -> 0 reference (Weber's second exponential integral):
    int_0^inf p^2 J_1(beta p) exp(-s p^2) dp = beta/(4 s^2) exp(-beta^2/4s).
    Used only to validate quad and as a visual overlay -- never as the
    production integrator, since it drops the alpha-dependent GUP terms.
    """
    return (1.0 / (4.0 * s ** 2)) * np.exp(-beta ** 2 / (4.0 * s))


def quad_safe(s, alpha, beta):
    """
    Hardened quad: tighter tolerances and a higher subdivision limit than the
    scipy default, plus an explicit non-convergence flag instead of silently
    trusting (or blanket-silencing warnings on) the result.
    """
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always", IntegrationWarning)
        val, err, *rest = quad(integrand, 0, np.inf, args=(s, alpha, beta),
                                limit=1000, epsabs=1e-12, epsrel=1e-10, full_output=1)
    flagged = len(rest) > 1   # quad appends a message (and sometimes explain) only when ier != 1
    return val, err, flagged


def _regression_check():
    """Quick sanity check: quad must match the alpha->0 closed form away from s~0."""
    for s, beta in [(0.5, 0.0), (1.0, 1.0), (2.0, -1.0)]:
        val, _, flagged = quad_safe(s, 1e-12, beta)
        ref = closed_form_alpha0(s, beta)
        rel = abs(val - ref) / abs(ref)
        assert not flagged and rel < 1e-6, (
            f"quad regression failed at s={s}, beta={beta}: rel={rel:.3e}, flagged={flagged}")


if __name__ == '__main__':
    _regression_check()

    s_vals = np.linspace(1e-6, 2, 60)
    beta_vals = np.linspace(-2, 2, 120)

    Z = np.empty((len(s_vals), len(beta_vals)))
    Zerr = np.empty_like(Z)
    Zflag = np.zeros_like(Z, dtype=bool)

    for i, s in enumerate(s_vals):
        for j, beta in enumerate(beta_vals):
            Z[i, j], Zerr[i, j], Zflag[i, j] = quad_safe(s, ALPHA, beta)

    n_flagged = int(Zflag.sum())
    print(f"Z min = {Z.min():.6g},  Z max = {Z.max():.6g}")
    print(f"quad flagged {n_flagged} / {Z.size} points ({n_flagged / Z.size:.2%}) as non-converged")
    bad_rows = s_vals[np.any(Zflag, axis=1)]
    if bad_rows.size:
        print(f"  non-converged rows span s in [{bad_rows.min():.3g}, {bad_rows.max():.3g}]")

    S, B = np.meshgrid(s_vals, beta_vals, indexing='ij')

    # --- Figure 1: raw surface -- shows the real near-s=0 behaviour, which is
    #     numerical noise from non-convergence, not the (far narrower, much
    #     smaller off beta=0) true divergence. ---------------------------------
    fig1 = plt.figure(figsize=(10, 7))
    ax1 = fig1.add_subplot(111, projection='3d')
    surf1 = ax1.plot_surface(S, B, Z, cmap='viridis', linewidth=0, antialiased=True)
    fig1.colorbar(surf1, ax=ax1, shrink=0.6, pad=0.1, label='Z')
    ax1.set_xlabel('s')
    ax1.set_ylabel(r'$\beta$')
    ax1.set_zlabel('Z')
    ax1.set_title(fr'$Z(s,\beta)$ -- raw, full range ($\alpha={ALPHA}$)' +
                  f'\n{n_flagged} unconverged points included (numerical noise, not signal)')
    fig1.tight_layout()
    fig1.savefig('TestingApparentIntegralDivergences/full_range_raw.png', dpi=150)

    # --- Figure 2: convergence-masked, percentile-clipped surface -- reveals
    #     the smooth bulk Gaussian ridge, validated against the closed form. ---
    abs_Z = np.abs(Z)
    finite = np.isfinite(Z)
    outlier_threshold = np.percentile(abs_Z[finite], OUTLIER_PCT)
    bad = (~finite) | Zflag | (Zerr > REL_TOL * abs_Z) | (abs_Z > outlier_threshold)

    Z_masked = Z.copy()
    Z_masked[bad] = np.nan

    n_masked = int(bad.sum())
    print(f"Masked {n_masked} / {Z.size} points ({n_masked / Z.size:.2%}) for the clipped view")

    finite_masked = Z_masked[np.isfinite(Z_masked)]
    zlo, zhi = np.percentile(finite_masked, [DISPLAY_PCT_LOW, DISPLAY_PCT_HIGH])
    print(f"Masked Z min = {finite_masked.min():.6g},  Z max = {finite_masked.max():.6g}")
    print(f"Clipped plot zlim/color limits (p{DISPLAY_PCT_LOW}/p{DISPLAY_PCT_HIGH}): [{zlo:.6g}, {zhi:.6g}]")

    fig2 = plt.figure(figsize=(10, 7))
    ax2 = fig2.add_subplot(111, projection='3d')
    surf2 = ax2.plot_surface(S, B, Z_masked, cmap='viridis', vmin=zlo, vmax=zhi,
                              linewidth=0, antialiased=True)
    ax2.set_zlim(zlo, zhi)

    Z_ref = np.clip(closed_form_alpha0(S, B), zlo, zhi)
    ax2.plot_wireframe(S, B, Z_ref, color='red', linewidth=0.3, alpha=0.4, rstride=4, cstride=8)

    fig2.colorbar(surf2, ax=ax2, shrink=0.6, pad=0.1, label='Z')
    ax2.set_xlabel('s')
    ax2.set_ylabel(r'$\beta$')
    ax2.set_zlabel('Z')
    ax2.set_title(fr'$Z(s,\beta)$ for $\alpha = ${ALPHA} -- convergence-masked, clipped to '
                  fr'[p{DISPLAY_PCT_LOW},p{DISPLAY_PCT_HIGH}]' +
                  '\n(red wireframe: exact $\\alpha\\to0$ closed form)')
    fig2.tight_layout()
    fig2.savefig('TestingApparentIntegralDivergences/clipped_masked.png', dpi=150)

    plt.show()