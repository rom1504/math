# Appendix: sharpness of the n^(3/4) bounded-cap tail and filter scale

2026-09-16. Status: **entire appendix independently audited PASS,
including weighted variance; exact finite replay and py_compile PASS**.
This is a sharpness example
for the class of arbitrary bounded-cap full signings, not a statement
about exact or near-minimizers. It complements the bounded-cap regime in
`artifacts/bh_2026_09_16_nonlinear_filter_audit.md`.

## 1. An explicit full signing with a planted rare-energy event

Let `t>=1`, `n=2^(4t)`, and `k=2^(3t)=n^(3/4)`, so `k^2=n^(3/2)`.
Index the Sylvester matrix by binary vectors and write

```math
S_n(i,j)=(-1)^{i\cdot j},\qquad S_n^2=nI,\qquad
\operatorname{tr}S_n=0.
```

Hollow it to obtain the full signing `B`, then change every edge in
the first `k` vertices to `+1`, obtaining `A`. Since the diagonal
contributes its trace at every spin vector,

```math
H_B(x)=\tfrac12x^TS_nx,\qquad Q(B)\le\tfrac12n^{3/2}.
```

The first `k` principal block is `S_k`; its total matrix sum is `k`
and its trace is zero. Consequently the original edge-sign sum in that
clique was `k/2`. The change `Delta=H_A-H_B` has only nonnegative
coefficients, so its supremum is its coefficient sum:

```math
\|\Delta\|_\infty=\binom k2-\frac k2=\frac{k(k-2)}2,
\qquad Q(A)\le n^{3/2}-k<n^{3/2}.                         \tag{1}
```

This uses the zero trace for the exact base cap; applying the operator
norm to the hollow matrix instead would introduce an unnecessary error.

Condition on all `k` planted spins being `+1`. The remaining spins are
independent uniform signs, so every cross or exterior edge has mean zero:

```math
\mathbb E[H_A\mid x_1=\cdots=x_k=1]=\binom k2.
```

Put `Y=H_A/n^(3/2)`. Then `|Y|<=1` by (1), and its conditional mean
is `1/2-1/(2k)`. For any fixed `0<a<1/2`, writing `p` for the
conditional probability of `Y>=a` gives
`1/2-1/(2k)<=a(1-p)+p`. Thus, for `k>=1/(1/2-a)`,

```math
\Pr\{Y\ge a\}
\ge c_a2^{-k},\qquad
c_a=\frac{1/2-a}{2(1-a)}>0.                               \tag{2}
```

In particular, no uniform bounded-cap tail estimate of the form
`C exp(-gamma n^beta)` with `beta>3/4`, fixed `C,gamma>0`, can hold
for this class at such a fixed threshold. The subsequence `n=16^t`
already falsifies a purported estimate uniform in all large dimensions.
The planted clique also gives `||A||op>=k-1`, by its supported all-one
Rayleigh vector; bounded scalar cap does not enforce order-`sqrt(n)`
operator norm.

## 2. Chebyshev filters detect this family at degree Theta(n^(3/4))

Fix a proposed cap `0<c<1/2`, choose `c<a<1/2`, and put
`kappa=arccosh(a/c)>0`. The polynomial

```math
P_d(y)=T_d(y/c)
```

has supremum at most one on `[-c,c]`. On the event in (2),
`P_d(Y)>=cosh(d kappa)>=exp(d kappa)/2`, whence

```math
\|P_d(Y)\|_2^2\ge\frac{c_a}{4}
 \exp\{2\kappa d-(\log2)k\}.                             \tag{3}
```

For any fixed `delta>0`, take

```math
d=\left\lceil\left(\frac{\log2}{2\kappa}+\delta\right)k
\right\rceil.
```

Then the right side of (3) grows at least as `exp(2 kappa delta k)`
times a positive constant. This is a genuine coefficient-BH detection:
`d=Theta(n^(3/4))=o(n)`, and the actual reduced Walsh degree is exactly
`m=2d` for all sufficiently large `n`. Indeed, the leading coefficient
of `P_d` is nonzero and each top `2d`-set coefficient of `H_A^d` is
`d!` times an odd hafnian. Thus the coefficient `ell_(q_m)` norm is
at least this exponentially growing `L2` norm, overwhelming any
degree-polynomial BH constant under the hypothetical filter bound
`||P_d(Y)||infinity<=1` implied by the proposed cap `||Y||infinity<=c`.

The same is true for the nonconstant weighted certificate. Parseval,
valid for every full signing, gives

```math
\mathbb E Y^2=\frac{\binom n2}{n^3},\qquad
\Pr\{|Y|\le c\}\ge1-\frac{\binom n2}{c^2n^3}\ge\frac12
```

eventually. On this low-energy event `|P_d(Y)|<=1`. The high and low
events are disjoint, and the chosen degrees eventually have
`exp(kappa d)/2>=1`. Applying
`Var(Z)=(1/2)E(Z-Z')^2` to their two ordered cross-events yields

```math
\operatorname{Var}(P_d(Y))
\ge\frac{c_a}{2}\,2^{-k}
 \left(\tfrac12e^{\kappa d}-1\right)^2.                   \tag{4}
```

For `s>=0`, the paper-type weighted functional satisfies
`W_s(g)>=m^(-s)||g-Eg||2`. Equations (3)--(4) therefore make both
the full and weighted coefficient certificates exponentially large at
the displayed degree, despite any fixed polynomial denominator in `m`.

This proves sharpness of the **power of n** in a uniform small-degree
non-detection theorem for the bounded-cap class. It does not determine
the optimal constant multiplying `n^(3/4)`, give a coefficient lower
bound uniform over all signings, or say anything about the original
minimization over signings. In fact the planted conditional mean forces
this example's normalized scalar cap to have liminf at least `1/2`.

## 3. Reproducible finite check

Run `.venv/bin/python computations/bh_2026_09_16_planted_tail_checks.py`.
At `n=16,k=8`, exhaustive spin evaluation gives actual cap `46`,
below the proved upper bound `56`. The conditional mean is `28`, and
the exact conditional probability of `H_A>=n^(3/2)/4=16` is `115/128`.
These checks illustrate the construction; the asymptotic conclusions
are proved by (1)--(4).
