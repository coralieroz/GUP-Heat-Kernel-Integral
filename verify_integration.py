"""
Diagnostic script verifying claims about the numerical integration of

    f(p) = [p^2(1+alpha^2 p^2)/beta] * exp{-s p^2 (1 + (2/3) alpha^2 p^2)} * J_1(beta p/(1+alpha^2 p^2))

over p in [0, inf), as used in plot_surface_stable_integrand.py / IntegralNumeric.py.

Uses the stable rewrite f(p) = p^3 * (J_1(A2)/A2) * exp(A1), A2 = beta p/(1+q),
q = (alpha p)^2, which is algebraically identical to f(p) and well-defined at beta=0.

CLOSED FORM (alpha -> 0 reference): the Weber second exponential integral
    int_0^inf p^2 J_1(beta p) exp(-s p^2) dp = beta/(4 s^2) * exp(-beta^2/(4s))
so for alpha -> 0,
    Z(s, beta) -> (1/(4 s^2)) * exp(-beta^2/(4s))
This matches Expression 2.14 of the capstone report (d=4 leading term,
(4 pi s)^-2 * exp(-beta^2/4s)) up to the constant 4*pi^2.

This script is read-only diagnostics: no production code is changed here.
"""

import warnings
import numpy as np
from scipy.integrate import quad
from scipy.special import j1


def integrand(p, s, alpha, beta):
    q = (alpha * p) ** 2
    A1 = -s * (p ** 2) * (1.0 + 2.0 * q / 3.0)
    A2 = (beta * p) / (1.0 + q)
    ratio = 0.5 if abs(A2) < 1e-8 else j1(A2) / A2
    return (p ** 3) * ratio * np.exp(A1)


def closed_form_alpha0(s, beta):
    """Exact alpha->0 value via the Weber integral."""
    return (1.0 / (4.0 * s ** 2)) * np.exp(-beta ** 2 / (4.0 * s))


def quad_eval(s, alpha, beta, limit=400, epsabs=1.49e-8, epsrel=1.49e-8):
    """Wraps quad with full_output so we can inspect ier / convergence messages.

    scipy.integrate.quad's full_output=1 return is (y, abserr, infodict) on
    clean success (ier==1 internally) and (y, abserr, infodict, message[, explain])
    when something was flagged; len(res) > 3 is therefore the non-converged signal.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        res = quad(integrand, 0, np.inf, args=(s, alpha, beta),
                    limit=limit, epsabs=epsabs, epsrel=epsrel, full_output=1)
    val, err = res[0], res[1]
    info = res[2]
    flagged = len(res) > 3
    msg = str(res[3]) if flagged else ""
    n_warnings = len(caught)
    neval = info.get("neval", None) if isinstance(info, dict) else None
    return val, err, int(flagged), msg, neval, n_warnings


# ─────────────────────────────────────────────────────────────────────────────
# 1. CLOSED-FORM CHECK: does quad reproduce the exact alpha->0 result?
# ─────────────────────────────────────────────────────────────────────────────

def check_closed_form():
    print("=" * 78)
    print("1. CLOSED-FORM CHECK  (alpha = 1e-12, i.e. alpha -> 0)")
    print("   quad(f) vs  Z = (1/4s^2) exp(-beta^2/4s)")
    print("=" * 78)
    alpha = 1e-12
    s_vals = [0.05, 0.1, 0.5, 1.0, 2.0, 5.0]
    beta_vals = [-2.0, -0.5, 0.0, 0.5, 2.0]
    hdr = f"{'s':>6} {'beta':>6}  {'quad Z':>16} {'closed form':>16} {'rel err':>10}  {'flag':>4} {'#warn':>5}"
    print(hdr)
    print("-" * len(hdr))
    worst = 0.0
    for s in s_vals:
        for beta in beta_vals:
            val, err, ier, msg, neval, nwarn = quad_eval(s, alpha, beta)
            ref = closed_form_alpha0(s, beta)
            rel = abs(val - ref) / abs(ref) if ref != 0 else float("nan")
            worst = max(worst, rel)
            print(f"{s:>6.2f} {beta:>6.2f}  {val:>16.8e} {ref:>16.8e} {rel:>10.2e}  {ier:>4} {nwarn:>5}")
    print(f"\n  Worst-case relative error across this grid: {worst:.3e}")
    print("  -> if small (<<1), quad is CORRECTLY evaluating the integral in")
    print("     the well-damped (s not tiny) regime; the 'near-zero' appearance")
    print("     in the surface plot is a SCALING artifact, not a quadrature failure.\n")


# ─────────────────────────────────────────────────────────────────────────────
# 2. BULK-VS-SPIKE MAGNITUDE
# ─────────────────────────────────────────────────────────────────────────────

def check_bulk_vs_spike():
    print("=" * 78)
    print("2. BULK-VS-SPIKE MAGNITUDE  (alpha = 1e-9, the 'no major blow-ups' case)")
    print("=" * 78)
    alpha = 1e-9
    s_vals = np.linspace(1e-6, 2, 60)
    beta_vals = np.linspace(-2, 2, 120)
    Z = np.empty((len(s_vals), len(beta_vals)))
    for i, s in enumerate(s_vals):
        for j, beta in enumerate(beta_vals):
            Z[i, j] = quad_eval(s, alpha, beta)[0]

    bulk_val = quad_eval(1.0, alpha, 0.0)[0]
    bulk_expected = closed_form_alpha0(1.0, 0.0)
    print(f"  Representative bulk point  Z(s=1, beta=0)        = {bulk_val:.6g}")
    print(f"  Closed-form prediction at (s=1, beta=0)          = {bulk_expected:.6g}")
    print(f"  Grid Z range                                     = [{Z.min():.6g}, {Z.max():.6g}]")
    print(f"  |Z| at the s_min corner (s=1e-6, beta=0)         = {abs(Z[0, len(beta_vals)//2]):.6g}")
    print(f"  Ratio bulk / corner                              = {abs(bulk_val)/abs(Z[0, len(beta_vals)//2]):.3e}")
    print("  -> the bulk value (~0.25 near beta=0) is dwarfed by the corner spike,")
    print("     which is why an autoscaled z-axis makes the bulk LOOK like zero.\n")
    return s_vals, beta_vals, Z


# ─────────────────────────────────────────────────────────────────────────────
# 3. CONVERGENCE / IER AUDIT
# ─────────────────────────────────────────────────────────────────────────────

def check_convergence_audit():
    print("=" * 78)
    print("3. CONVERGENCE / IER AUDIT  (alpha = 1e-9, limit=400, default tol)")
    print("=" * 78)
    alpha = 1e-9
    s_vals = np.linspace(1e-6, 2, 30)
    beta_vals = np.linspace(-2, 2, 60)
    n_total = 0
    n_bad_ier = 0
    n_bad_relerr = 0
    rel_tol = 1e-3
    bad_s_hist = []
    for s in s_vals:
        bad_here = 0
        for beta in beta_vals:
            val, err, ier, msg, neval, nwarn = quad_eval(s, alpha, beta)
            n_total += 1
            if ier != 0:
                n_bad_ier += 1
            relerr = err / abs(val) if val != 0 else float("inf")
            if relerr > rel_tol:
                n_bad_relerr += 1
                bad_here += 1
        bad_s_hist.append(bad_here)
    print(f"  Total evaluations                  : {n_total}")
    print(f"  Non-zero ier (quad flags itself)   : {n_bad_ier}  ({100*n_bad_ier/n_total:.1f}%)")
    print(f"  Relative error estimate > {rel_tol:g}     : {n_bad_relerr}  ({100*n_bad_relerr/n_total:.1f}%)")
    print(f"\n  Bad-point count by s-row (first 10 rows, smallest s first):")
    for s, cnt in list(zip(s_vals, bad_s_hist))[:10]:
        print(f"    s={s:.6g}: {cnt} / {len(beta_vals)} flagged")
    print("  -> if flagged points concentrate at small s, non-convergence is")
    print("     localized there (consistent with weak damping + oscillatory tail),")
    print("     not a global quadrature failure.\n")


# ─────────────────────────────────────────────────────────────────────────────
# 4. BULLET-3 REPRODUCTION: s_min / resolution dependence
# ─────────────────────────────────────────────────────────────────────────────

def check_s_min_and_resolution():
    print("=" * 78)
    print("4. s_min AND RESOLUTION DEPENDENCE  (alpha = 1e-50, beta = 1.0 slice)")
    print("=" * 78)
    alpha = 1e-50
    beta = 1.0

    print("  (a) s_min sensitivity at fixed resolution (n=50 points from s_min to 2):")
    for s_min in [1e-6, 1e-3, 1e-1]:
        s_vals = np.linspace(s_min, 2, 50)
        vals = [quad_eval(s, alpha, beta)[0] for s in s_vals]
        refs = [closed_form_alpha0(s, beta) for s in s_vals]
        worst_rel = max(abs(v - r) / abs(r) for v, r in zip(vals, refs))
        print(f"    s_min={s_min:>8.1e}: Z(s_min)={vals[0]:>14.6g}  "
              f"closed_form(s_min)={refs[0]:>14.6g}  worst rel err over row={worst_rel:.3e}")

    print("\n  (b) resolution sensitivity at fixed s_min=1e-6 (does Z(s_min) change with n?):")
    for n in [30, 60, 120, 240]:
        s_vals = np.linspace(1e-6, 2, n)
        val = quad_eval(s_vals[0], alpha, beta)[0]
        ref = closed_form_alpha0(s_vals[0], beta)
        print(f"    n={n:>4}: Z(s_min)={val:>14.6g}  closed_form={ref:>14.6g}  "
              f"rel err={abs(val-ref)/abs(ref):.3e}")
    print("  -> Z(s_min) itself should NOT depend on the grid's resolution at all")
    print("     (it's evaluated at the same s,beta regardless of n); if it agrees")
    print("     with the closed form to high precision, quad genuinely converges")
    print("     there too, and the earlier-reported resolution-dependence of the")
    print("     PEAK LOCATION in IntegralNumeric.py/PeakInvestigation.py traces to")
    print("     their different (non-rewritten, beta-dividing) integrand instead.\n")


# ─────────────────────────────────────────────────────────────────────────────
# 5. ALPHA-SCALING OF THE GUP CORRECTION
# ─────────────────────────────────────────────────────────────────────────────

def check_alpha_scaling():
    print("=" * 78)
    print("5. ALPHA-SCALING  (does the blow-up severity track alpha as expected?)")
    print("=" * 78)
    s, beta = 0.5, 1.0
    leading = closed_form_alpha0(s, beta)
    print(f"  Leading term  1/(4s^2) exp(-beta^2/4s)  at (s={s}, beta={beta}) = {leading:.6g}")
    print(f"  {'alpha':>10}  {'quad Z':>16}  {'(Z-leading)/leading':>20}")
    for alpha in [1e-50, 1e-15, 1e-9, 1e-6, 1e-3, 1e-2, 1e-1]:
        val = quad_eval(s, alpha, beta)[0]
        dev = (val - leading) / leading
        print(f"  {alpha:>10.1e}  {val:>16.8e}  {dev:>20.3e}")
    print("  -> deviation from the alpha=0 leading term should grow with alpha")
    print("     (analytically ~ alpha^2 beta^4/(6 s^3) per Expression 2.14's GUP")
    print("     correction), confirming larger-alpha blow-ups are a real GUP effect")
    print("     amplified by the same 1/s^k divergence, not purely numerical noise.\n")


if __name__ == "__main__":
    check_closed_form()
    check_bulk_vs_spike()
    check_convergence_audit()
    check_s_min_and_resolution()
    check_alpha_scaling()
