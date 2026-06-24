# Stable Integrand

**Source:** plot_surface_stable_integrand.py, StableIntegrandPlots/
**Status:** active

Algebraic rewrite of the [I(s, alpha, beta)](Integral.md) integrand to remove
a spurious divergence as `beta -> 0`.

## The rewrite

The original integrand has `T1 = p^2(1+q)/beta` (with `q = (alpha p)^2`),
which blows up as `beta -> 0`. But using `A2 = beta p / (1+q)`:

```
T1 = p^2(1+q)/beta = p^3 / A2
```

so `T1 * J_1(A2) = p^3 * J_1(A2)/A2`, which is **entire** — `J_1(A2)/A2 -> 1/2`
as `A2 -> 0`, so the limit is `p^3/2`, no singularity. The script implements
this directly:

```python
ratio = 0.5 if abs(A2) < 1e-8 else j1(A2) / A2
return (p ** 3) * ratio * np.exp(A1)
```

This eliminates the spurious spike in the small-`(s, beta)` corner of the
surface plot that the naive `T1 * exp(A1) * j1(A2)` form produced.

## Output

Surface plots `Z(s, beta)` for fixed `alpha = 1e-9` over `s in [1e-6, 2]`,
`beta in [-2, 2]`, saved per-alpha under `StableIntegrandPlots/`
(`alpha=1e-15.png` through `alpha=1e-50.png`, plus others).

A commented-out second pass masks remaining oscillatory-undersampling
artifacts using `quad`'s own relative-error estimate (`Zerr`) and a
percentile-based outlier backstop — currently disabled, not wired into the
active script.

## Open question

Whether this rewrite alone is sufficient, or whether the small-`s` peak from
[Peak Investigation](PeakInvestigation.md) can still appear here for some
`(alpha, beta)` since the `s`-domain still approaches `1e-6`. Not yet checked
explicitly against the lower-bound/resolution tests from that investigation.

## Related pages

[Integral](Integral.md), [Peak Investigation](PeakInvestigation.md)
