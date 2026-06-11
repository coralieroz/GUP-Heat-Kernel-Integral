"""
Parametric integral I(s, alpha, beta) as an explicit symbolic expression via
N-point Gauss-Laguerre quadrature.

Substitution p = sqrt(t/s) maps [0, inf) -> [0, inf) and factors exp(-t)
as the Gauss-Laguerre weight, leaving a smooth, rapidly decaying integrand h(t).
The J_1 argument beta*sqrt(t/s)/(1 + alpha^2*t/s) is bounded above by
beta/(2*alpha) (maximum of sqrt(u)/(1+u) at u = 1), so the integrand is
non-oscillatory and converges fast with Laguerre nodes.  Larger beta/alpha ratios
push the effective Bessel-argument range higher and may need more nodes.

NOTE ON CONVERGENCE: the rate at which N nodes suffice depends on how concentrated
h(t) is in t-space.  After substitution, the additional damping factor
exp(-2*alpha^2*t^2/(3*s)) becomes very sharp near t = 0 when s is small (e.g.
s = 0.1), placing the integrand peak at t ~ sqrt(s/(2*alpha^2)) << 1.
Gauss-Laguerre nodes are sparse there, so small-s cases require larger N.
Similarly, large beta/alpha pushes the Bessel function into oscillatory territory,
also requiring more nodes.

Gauss-Hermite is NOT used: the integrand is one-sided (domain [0, inf)) and
symmetrising to the full line is invalid because the integrand is odd in p.
Laguerre quadrature is the correct choice.
"""

import warnings
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import j1, jv, roots_laguerre

# ── Global symbolic variables ──────────────────────────────────────────────────
s, alpha, beta = sp.symbols('s alpha beta', positive=True)


# ─────────────────────────────────────────────────────────────────────────────
# 1. PRIMARY: build a fully symbolic SymPy expression
# ─────────────────────────────────────────────────────────────────────────────

def build_I_expr(N: int = 40) -> sp.Expr:
    """
    Return a SymPy Add representing I(s, alpha, beta) via N-point
    Gauss-Laguerre quadrature.

    With p = sqrt(t/s),  dp = t^(-1/2) / (2 sqrt(s)) dt, the integrand becomes

        I = int_0^inf e^{-t} h(t) dt

    where
        h(t) =   t^(-1/2) / (2 sqrt(s))
               * (t/s)(1 + alpha^2 t/s) / beta
               * exp(-2 alpha^2 t^2 / (3s))
               * besselj(1,  beta sqrt(t/s) / (1 + alpha^2 t/s))

    The e^{-t} factor is absorbed into the Laguerre weight, giving

        I ~= sum_i  w_i * h(t_i)

    Nodes t_i and weights w_i come from numpy; each is stored as a
    high-precision sympy.Float so the returned object is a genuine symbolic Add
    that supports diff, subs, lambdify, etc.
    """
    # scipy's roots_laguerre is numerically stable up to N ~ 300 in float64;
    # numpy's laggauss overflows for N >= ~170.  The overflow warning is benign
    # because the NaN guard below discards any corrupted nodes.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        nodes, weights = roots_laguerre(N)
    terms = []
    for t_np, w_np in zip(nodes, weights):
        # skip nodes whose weight underflows to zero or NaN (can happen at very large t)
        if not (np.isfinite(t_np) and np.isfinite(w_np) and w_np > 0):
            continue
        ti = sp.Float(t_np, 50)
        wi = sp.Float(w_np, 50)

        a2_t_s     = alpha**2 * ti / s                            # alpha^2 * t_i / s
        prefactor  = ti**sp.Rational(-1, 2) / (2 * sp.sqrt(s))
        poly_piece = (ti / s) * (1 + a2_t_s) / beta
        exp_piece  = sp.exp(sp.Rational(-2, 3) * alpha**2 * ti**2 / s)
        bessel_arg = beta * sp.sqrt(ti / s) / (1 + a2_t_s)

        h_i = prefactor * poly_piece * exp_piece * sp.besselj(1, bessel_arg)
        terms.append(wi * h_i)

    return sp.Add(*terms)


# ─────────────────────────────────────────────────────────────────────────────
# 2. Build symbolic expression + fast lambdified numeric function
# ─────────────────────────────────────────────────────────────────────────────

def build_I(N: int = 40):
    """
    Return (expr, numeric_fn) where
      expr       - SymPy Add in symbols (s, alpha, beta); supports diff, subs, etc.
      numeric_fn - lambdified callable f(s_val, alpha_val, beta_val) using
                   scipy.special.jv for besselj so all Bessel evaluations are
                   numerically correct.
    """
    expr = build_I_expr(N)
    numeric_fn = sp.lambdify(
        (s, alpha, beta),
        expr,
        modules=[{'besselj': jv, 'sqrt': np.sqrt, 'exp': np.exp}, 'numpy'],
    )
    return expr, numeric_fn


# ─────────────────────────────────────────────────────────────────────────────
# Helper: reference value via scipy quad on the original integrand
# ─────────────────────────────────────────────────────────────────────────────

def _reference(s_val: float, alpha_val: float, beta_val: float):
    def integrand(p):
        q   = (alpha_val * p) ** 2
        T1  = p**2 * (1.0 + q) / beta_val
        T2  = np.exp(-s_val * p**2 * (1.0 + (2.0 / 3.0) * q))
        arg = beta_val * p / (1.0 + q)
        return T1 * T2 * j1(arg)

    result, err = quad(integrand, 0.0, np.inf, limit=500, epsabs=1e-14, epsrel=1e-14)
    return result, err


# ─────────────────────────────────────────────────────────────────────────────
# Main: validation, convergence, symbolic demonstration
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':

    # ── 1. Build N=40 expression ──────────────────────────────────────────────
    print("Building N=40 symbolic expression ...", flush=True)
    expr40, I_num40 = build_I(N=40)
    print(f"  Type            : {type(expr40).__name__}")
    print(f"  Number of terms : {len(expr40.args)}")
    print(f"  First term      : {expr40.args[0]}\n")

    # Confirm symbolic operations work
    print("  diff(expr40, beta) -- first term:")
    print(f"  {sp.diff(expr40.args[0], beta)}\n")

    # ── 3. Validation against scipy.integrate.quad ────────────────────────────
    # Each row shows the Laguerre approximation vs the quad reference.
    # The recommended N for each case appears in parentheses; use build_I(N=N_rec)
    # for production runs.  Cases with small s or large beta/alpha are harder
    # because h(t) is concentrated near t ~ 0, where Laguerre nodes are sparse.
    test_cases = [
        # (s, alpha, beta,  N_recommended)
        (1.0,  1.0,  1.0,   40),
        (0.1,  1.0,  5.0,  200),   # small s -> integrand peaked near t~0
        (10.0, 0.5,  2.0,   40),
        (2.0,  2.0,  0.3,   60),
        (0.5,  0.2, 10.0,  150),   # large beta/alpha -> wide Bessel range
    ]
    hdr = (f"{'s':>5} {'alpha':>5} {'beta':>5}  "
           f"{'N':>5}  {'Laguerre':>18}  {'quad reference':>18}  {'rel error':>11}")
    print(hdr)
    print("-" * len(hdr))
    for sv, av, bv, N_rec in test_cases:
        _, fn = build_I(N_rec)
        lag = float(fn(sv, av, bv))
        ref, _ = _reference(sv, av, bv)
        rel = abs(lag - ref) / abs(ref) if ref != 0.0 else float('nan')
        print(f"{sv:>5.1f} {av:>5.2f} {bv:>5.1f}  "
              f"{N_rec:>5}  {lag:>18.12e}  {ref:>18.12e}  {rel:>11.3e}")

    # ── 4. Convergence in N ───────────────────────────────────────────────────
    print("\nConvergence in N  (s=1, alpha=1, beta=1) -- easy case:")
    ref1, _ = _reference(1.0, 1.0, 1.0)
    print(f"  Reference (quad) : {ref1:.15e}\n")
    print(f"  {'N':>4}  {'Laguerre value':>22}  {'rel error':>12}")
    for N in [10, 20, 30, 40, 60]:
        _, fn = build_I(N)
        val = float(fn(1.0, 1.0, 1.0))
        rel = abs(val - ref1) / abs(ref1)
        print(f"  {N:>4}  {val:>22.15e}  {rel:>12.3e}")

    # Small-s hard case: after substitution, the additional factor
    # exp(-2*alpha^2*t^2/(3*s)) = exp(-6.67*t^2) concentrates the integrand
    # near t ~ 0.3, while Laguerre nodes are sparse there.  Convergence is
    # therefore algebraic in N rather than exponential; large N is needed.
    print("\nConvergence in N  (s=0.1, alpha=1, beta=5) -- small-s hard case:")
    ref2, _ = _reference(0.1, 1.0, 5.0)
    print(f"  Reference (quad) : {ref2:.15e}\n")
    print(f"  {'N':>4}  {'Laguerre value':>22}  {'rel error':>12}")
    # scipy's roots_laguerre begins to overflow for N >> 300 in float64;
    # beyond that, higher-precision arithmetic (e.g. mpmath) is needed.
    for N in [40, 80, 150, 250]:
        _, fn = build_I(N)
        val = float(fn(0.1, 1.0, 5.0))
        if not np.isfinite(val):
            print(f"  {N:>4}  {'(node overflow -- use mpmath)':>35}  {'---':>12}")
        else:
            rel = abs(val - ref2) / abs(ref2)
            print(f"  {N:>4}  {val:>22.15e}  {rel:>12.3e}")

    # ── 5. Symbolic differentiation dI/ds with N=4 ───────────────────────────
    print("\n--- Symbolic dI/ds  (N=4, readable size) ---")
    expr4 = build_I_expr(N=4)
    dI_ds = sp.diff(expr4, s)
    print(f"  dI/ds is a sum of {len(dI_ds.args)} terms.  First term:")
    sp.pprint(dI_ds.args[0], use_unicode=False)
    print("\n  Full expression as LaTeX (paste directly into paper):")
    print(sp.latex(dI_ds))

    # Numeric sanity: symbolic derivative vs. central finite difference
    _, fn4 = build_I(N=4)
    dI_ds_fn = sp.lambdify(
        (s, alpha, beta),
        dI_ds,
        modules=[{'besselj': jv, 'sqrt': np.sqrt, 'exp': np.exp}, 'numpy'],
    )
    h_fd    = 1e-5
    fd_est  = (fn4(1.0 + h_fd, 1.0, 1.0) - fn4(1.0 - h_fd, 1.0, 1.0)) / (2.0 * h_fd)
    sym_val = float(dI_ds_fn(1.0, 1.0, 1.0))
    print(f"\n  Symbolic  dI/ds at (s=1, alpha=1, beta=1) [N=4] : {sym_val:.8e}")
    print(f"  Finite-difference estimate                       : {fd_est:.8e}")
    print(f"  Relative difference                              : "
          f"{abs(sym_val - fd_est) / abs(fd_est):.3e}")
