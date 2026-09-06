# Independent reconstruction: Gaussian deformation and confined wells

2026-09-06. Director audit of
`decisive_independent_gaussian_stability_2026_09_06.md`.
Verdict: PASS for the stated necessary condition; no convergence inference.

I reconstructed the replacement rather than using earlier audit verdicts.
For the normalized soft maximum

```math
L(W)=(\beta n)^{-1}\log\sum_{s,x}
 \exp[\beta sH_W(x)/\sqrt n],
```

the edge third derivative is a third cumulant of a variable in {-1,1}
times beta^2/n^(5/2), hence bounded by 8 beta^2/n^(5/2).
Center independent signs at a A_e, a=(1+t^2)^(-1/2); their variances are
1-a^2, matching the Gaussian replacements. Taylor's third remainder,
with absolute third moments at most 8 and 2, summed over at most n^2/2
edges, costs at most (20/3) beta^2/sqrt(n). The soft-maximum gap is at
most 2 log(2)/beta. Thus beta=n^(1/6) gives the claimed constant below
10. Every un-replaced complete sample is an actual signing, so its cap
is at least M_n; division by a gives precisely

```math
n^{-3/2}\mathbb E Q(A+tG)
\ge \sqrt{1+t^2}(m_n-10n^{-1/6}).
```

This holds for every A, and only comparison with Q(A) uses near-minimality.
The noise parameter is fixed before n; no uniform t growing with n is
asserted. The absolute objective is handled by the extra sign s throughout.

For centers z at distance r modulo reversal, exactly r(n-r) edges change
in H_G(x)-H_G(z), each with coefficient magnitude two. The increment
variance is therefore 4r(n-r), not twice that quantity. At most
2K_n binom(n,r) pairs occur. Each shell maximum has Gaussian Lipschitz
constant at most n. Gaussian concentration bounds the expected maximum
over n+1 centered shell maxima by n sqrt(2 log(n+1)); after multiplying
by fixed t and normalizing this vanishes. This argument does not assume
independence between shells, centers, or increments.

Writing eta=lim log(K_n)/n yields the source's bound

```math
c(\sqrt{1+t^2}-1)
\le 3t\sqrt\eta+e\kappa\exp[-\kappa^2/(8t^2)]
\qquad(0<t\le\kappa/4).
```

The constants follow from H(rho)<=rho log(e/rho), rho<=1/2, and
monotonicity of exp(-L)sqrt(L) for L>=2. Taking t small rules out eta=0.
The positive elementary lower bound c>=1/4 uses a bipartition, independent
random signs on one side, optimized signs on the other, Khintchine, and
global reversal of one side; it is independent of the reported .433 bound.

## What was and was not removed

The theorem excludes a particular uniformly linearly confined,
subexponential-well geometry for asymptotic minimizers. It does not assert
exponentially many exact maximizers, uniform spectral flatness, a useful
insertion row, or an all-order realization. Its current role is a falsifier
for proposed minimizer models. It is not a substitute for the limit theorem.
