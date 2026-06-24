# IntegralAnalytic.py

**Source:** IntegralAnalytic.py
**Status:** active

Builds [I(s, alpha, beta)](Integral.md) as a fully symbolic SymPy expression
via N-point Gauss-Laguerre quadrature, instead of evaluating it numerically
point-by-point.

## Method

Substituting `p = sqrt(t/s)` maps `[0, inf) -> [0, inf)` and turns `exp(-s p^2)`
into `exp(-t)`, which is exactly the Gauss-Laguerre weight. What remains is a
smooth, non-oscillatory integrand `h(t)` (the `J_1` argument is bounded above
by `beta/(2 alpha)`), so `I ~= sum_i w_i h(t_i)` converges fast for typical
parameters. Each node/weight is stored as a 50-digit `sympy.Float`, so the
result is a true symbolic `Add` supporting `diff`, `subs`, `lambdify`.

`build_I(N)` returns `(expr, numeric_fn)`; `numeric_fn` uses `scipy.special.jv`
for the Bessel evaluation so it's numerically correct even though `expr` is
symbolic.

## Convergence caveats (documented in the module docstring)

- Small `s` concentrates `h(t)` near `t ~ sqrt(s/(2 alpha^2)) << 1`, where
  Laguerre nodes are sparse — needs larger `N` (e.g. `s=0.1, beta=5` needs
  `N~200` vs `N=40` for the easy case).
- Large `beta/alpha` pushes the Bessel argument into oscillatory territory,
  also needing more nodes.
- `scipy.special.roots_laguerre` is stable to `N ~ 300` in float64;
  `numpy.polynomial.laguerre.laggauss` overflows above `N ~ 170`.
- Gauss-Hermite is deliberately *not* used — the integrand is one-sided and
  odd in `p`, so symmetrizing to the full line would be invalid.

## What `__main__` demonstrates

1. Builds `N=40` expression, confirms it supports symbolic `diff`.
2. Validates Laguerre approximation against `scipy.integrate.quad` reference
   across 5 `(s, alpha, beta)` cases with hand-picked recommended `N`.
3. Shows convergence in `N` for an easy case (`s=1,alpha=1,beta=1`) and a hard
   small-`s` case (`s=0.1,alpha=1,beta=5`) — convergence is algebraic, not
   exponential, in the hard case.
4. Computes `dI/ds` symbolically at `N=4` and cross-checks against a central
   finite difference; prints the LaTeX form for direct use in the paper.

## Related pages

[Integral](Integral.md), [IntegralNumeric.py](IntegralNumeric.md)
