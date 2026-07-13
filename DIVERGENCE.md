## Investigate integral blow-up at $s\to0$ and $\beta\to0$: bounded above or divergent?
In the limit under consideration, the integral

```math
I(s,\alpha,\beta)=\int_0^\infty\frac{p^2(1+\alpha^2p^2)}{\beta}\text{exp}\bigl\lbrace -sp^2(1+\frac23\alpha^2p^2) \bigr\rbrace\mathcal{J}_1\bigl( \frac{\beta p}{1+\alpha^2p^2} \bigr)dp
```

comes down to a competition between the $`1/\beta`$ term and the Bessel function as the exponent $`\text{exp}\lbrace-s(...)\rbrace\to1`$.  For small $`x`$, $`\mathcal{J}_1(x)=\frac x2-\frac{x^3}{16}`$ and so the competing product can be expressed as

```math
\frac{p^2(1+\alpha^2p^2)}{\beta}\cdot\frac{\beta p}{2(1+\alpha^2p^2)}-\frac{p^2(1+\alpha^2p^2)}{\beta}\cdot\frac{\beta^3 p^3}{16(1+\alpha^2 p^2)^3}=\frac{p^3}2-\frac{\beta^2p^5}{16(1+\alpha^2p^2)^2}.
```

Taking the appropriate limit, the second term drops out and due to linearity of the limit and integration, $`\lim_{s,\beta\to0}I(s,\alpha,\beta)=\int_0^\infty\frac{p^3}2\ dp`$ diverges.