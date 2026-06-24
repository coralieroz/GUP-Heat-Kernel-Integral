"""
Production evaluator for Z(s, alpha, beta) using the exact p = sqrt(t/s)
Gauss-Laguerre substitution, instead of scipy.integrate.quad applied directly
to the oscillatory, weakly-damped integral on [0, inf).

Why: for the alpha values used throughout this project (alpha ~ 1e-9 or
smaller), scipy.integrate.quad on the raw integrand fails to converge once s
gets small (see verify_integration.py / plot_surface_stable_integrand.py): the
damping exp(-s p^2) is too weak relative to the oscillatory J_1 tail for
quad's adaptive subdivision to resolve. Folding exp(-s p^2) into the
Gauss-Laguerre weight via p = sqrt(t/s) removes that specific failure mode and
needs no expansion or limit -- it is the exact same integral, exactly
substituted. Node/weight generation uses scipy.special.roots_laguerre.

This module uses the beta=0-safe rewrite of the integrand (the same one used
in IntegralNumeric.py, PeakInvestigation.py and
plot_surface_stable_integrand.py to remove the spurious p^2(1+q)/beta
divergence as beta -> 0):

    p^3 * J_1(A2)/A2 * exp(A1),  A2 = beta*p/(1+q), A1 = -s p^2 (1+2q/3), q = alpha^2 p^2

Under p = sqrt(t/s) (so exp(A1) = exp(-t) * exp(-2 alpha^2 t^2/(3s)) and
p^3 dp = [t/(2 s^2)] exp(-t) dt, with exp(-t) absorbed into the Laguerre
weight), this becomes:

    h(t) = [t / (2 s^2)] * exp(-2 alpha^2 t^2/(3s)) * J_1(A2)/A2,
    A2 = beta*sqrt(t/s) / (1 + alpha^2 t/s)

with J_1(A2)/A2 -> 1/2 as A2 -> 0 (beta=0 or t=0), well-defined throughout,
including at beta=0.

CALIBRATED RELIABILITY (alpha ~ 1e-9, the regime this project actually uses):
for s >= S_MIN_RELIABLE, this sum matches the alpha->0 closed form (Weber's
2nd exponential integral) to ~1e-10 or better with only N=40 nodes -- and
adding more nodes does NOT help below that threshold. That floor is therefore
NOT an under-resolution problem (more nodes would fix that); it is float64
catastrophic cancellation in the weighted sum, since each term's effective
Bessel argument beta*sqrt(t_i/s) grows like 1/sqrt(s) and the terms become
large and oscillate in sign. This is a genuine open numerical problem
(confirmed separately: even mpmath's general-purpose arbitrary-precision quad
fails the same way for s ~ 1e-3..1e-6, because the difficulty is the highly
oscillatory, weakly-damped tail itself, not floating-point precision) -- it
would need a dedicated oscillatory-quadrature method (e.g. Levin-type or
Filon-type collocation, or mpmath.quadosc with an explicit period), which is
out of scope here. Below S_MIN_RELIABLE this module raises rather than
silently returning a number that looks plausible but is numerical noise.

Calibration (alpha=1e-9, beta=1.0, against the alpha->0 closed form):
    N            s=0.005      s=0.01       s=0.02       s=0.05       s=0.1
    40           1.3e+05      2.0e-06      3.1e-11      1.4e-14      2.2e-15
    80           3.0e+05      4.7e-07      4.4e-11      9.6e-14      5.4e-15
    150          1.9e+05      1.5e-05      5.7e-11      6.9e-14      5.0e-15
    300          7.2e+05      2.0e-05      1.7e-10      8.8e-14      2.9e-14
(rel. error vs closed form.) N=40 is already as good as N=300 -- consistent
with the cancellation diagnosis above, not a sparse-node problem.
"""

from functools import lru_cache

import numpy as np
from scipy.special import j1, roots_laguerre

S_MIN_RELIABLE = 0.02   # below this, the Laguerre sum is numerical noise (see module docstring)
DEFAULT_N = 40          # sufficient for s >= S_MIN_RELIABLE at the alpha scales used in this project


class UnreliableRegionError(ValueError):
    """Raised when s falls below the calibrated reliability floor."""


@lru_cache(maxsize=None)
def _nodes_weights(N: int):
    """Laguerre nodes/weights for N, with non-finite/non-positive entries dropped
    (roots_laguerre can return spurious entries at large N)."""
    nodes, weights = roots_laguerre(N)
    keep = np.isfinite(nodes) & np.isfinite(weights) & (weights > 0)
    return nodes[keep], weights[keep]


def recommended_N(s: float, alpha: float, beta: float) -> int:
    """
    Node-count heuristic: large beta/alpha pushes the Bessel argument into
    oscillatory territory and needs more nodes. At the tiny alpha used
    throughout this project that ratio is enormous but never actually binding
    (the t-range Laguerre nodes cover never reaches alpha^2*t/s ~ O(1)), so
    DEFAULT_N=40 is calibrated to already suffice; this heuristic only
    escalates N for alpha that is *not* tiny.
    """
    if alpha < 1e-6:
        return DEFAULT_N
    ratio = abs(beta) / alpha
    if ratio > 5:
        return 200
    if ratio > 2:
        return 80
    return DEFAULT_N


def Z_laguerre(s: float, alpha: float, beta: float, N: int | None = None) -> float:
    """
    Exact-substitution, beta=0-safe Gauss-Laguerre evaluation of Z(s, alpha, beta).
    Raises UnreliableRegionError if s < S_MIN_RELIABLE (see module docstring).
    """
    if s < S_MIN_RELIABLE:
        raise UnreliableRegionError(
            f"s={s:.3g} < S_MIN_RELIABLE={S_MIN_RELIABLE:.3g}: Laguerre sum is "
            "float64 cancellation noise here, not a reliable value (see module docstring).")
    if N is None:
        N = recommended_N(s, alpha, beta)
    t, w = _nodes_weights(N)

    q = alpha ** 2 * t / s
    A2 = beta * np.sqrt(t / s) / (1.0 + q)
    nonzero = np.abs(A2) >= 1e-8
    ratio = np.full_like(A2, 0.5)
    ratio[nonzero] = j1(A2[nonzero]) / A2[nonzero]
    h = (t / (2.0 * s ** 2)) * ratio * np.exp(-2.0 * alpha ** 2 * t ** 2 / (3.0 * s))
    return float(np.sum(w * h))
