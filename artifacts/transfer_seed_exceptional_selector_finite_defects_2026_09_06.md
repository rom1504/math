# Finite necessary defects for exceptional low-cap selectors

Date: 2026-09-06. Seed-transfer track. All finite inequalities concern
actual hollow symmetric sign matrices, not an ensemble certificate.

## 1. A covariance clipping inequality without a diagonal floor

Let `A` have order `n`, put `L=A/sqrt(n)`, and let `R` be any nonzero
positive semidefinite matrix. Define

```math
v=\operatorname{tr}R/n,\quad e=\operatorname{tr}(LR)/n,\quad
s=\operatorname{tr}R^2/n,\quad
\delta^2=n^{-1}\sum_i(R_{ii}-v)^2.
```

For EVERY `b>0`,

```math
\boxed{\quad
\frac{Q(A)}{n^{3/2}}\ge\frac1\pi
\left[\frac{e-2\delta\sqrt s/b}{v+b}
       -\frac{s}{(v+b)^2\sqrt n}\right].\quad}             (1)
```

Unlike the earlier diagonal-normalization inequality, (1) has no
assumption `R_ii>=a>0`. It permits a singular covariance with zero
diagonal entries.

Proof. Keep `I={i:R_ii<=v+b}` and let `J` be its complement. Chebyshev
gives `|J|/n<=delta^2/b^2`. The matrix

```math
C=R_I/(v+b)+\operatorname{diag}_{i\in I}(1-R_{ii}/(v+b))
```

is a correlation matrix. Gaussian signs on `I`, extended by independent
unbiased signs on `J`, are actual full-order Boolean spins. Their expected
normalized energy is `sum_(i,j in I) L_ij arcsin(C_ij)/(pi n)`.
Since `|arcsin z-z|<=|z|^3<=|z|^2`, its difference from the linear term
is at most `s/[pi(v+b)^2 sqrt(n)]`; the hollow diagonal contributes
nothing. Finally, the exact row squared norms of `L` are at most one,
so Frobenius Cauchy--Schwarz gives

```math
\frac1n|\operatorname{tr}(LR)-\operatorname{tr}(L_I R_I)|
\le2\sqrt{|J|/n}\sqrt s\le2\delta\sqrt s/b.
```

This also covers an empty `I`, when the actual expected energy is zero.
These facts prove (1).

Consequently, if `Q(A)<=c n^(3/2)` with `c>0`, then

```math
e\le\pi c(v+b)+2\delta\sqrt s/b+s/[(v+b)\sqrt n]           (2)
```

for every `b>0`. Thus the covariance trace ratio cannot be high while
its diagonal is uniform and its second trace is controlled. This is a
finite necessary condition, not an assertion that such covariance
conditions hold on exceptional selectors.

## 2. An explicit two-sided cubic test using only even spectral moments

Fix the rational polynomial

```math
P(x)=x^3+\frac85x^2-\frac25x-\frac35,
\qquad R_+=P(L)^2,\quad R_-=P(-L)^2.
```

Write `mu_k=tr(L^k)/n` and put

```math
\bar v=\frac9{25}-\frac{44}{25}\mu_2
                  +\frac{44}{25}\mu_4+\mu_6,
\qquad
\bar e=\frac{12}{25}\mu_2-\frac{62}{25}\mu_4
                              +\frac{16}{5}\mu_6,         (3)
```

and

```math
S=\frac{\operatorname{tr}(R_+^2+R_-^2)}{2n},\qquad
\Delta^2=\frac12\left[\operatorname{Var}_i(R_+)_{ii}
                         +\operatorname{Var}_i(R_-)_{ii}\right]. (4)
```

Here all variances use uniform averaging over the `n` vertices.
The quantity `S` depends ONLY on even spectral moments through degree
12. The local quantity in (4) has the more concrete expression

```math
\Delta^2=
\operatorname{Var}_i\left[\frac{44}{25}(L^4)_{ii}+(L^6)_{ii}\right]
+\operatorname{Var}_i\left[-\frac{62}{25}(L^3)_{ii}
                            +\frac{16}{5}(L^5)_{ii}\right]. (5)
```

Indeed the constant and degree-two diagonal terms in `P(L)^2` are
constant across vertices, since `(L^2)_ii=(n-1)/n`, and `(L)_ii=0`.
The variance identity is `Var(X+Y)+Var(X-Y)=2 Var X+2 Var Y`.

If the actual cap is at most `c n^(3/2)`, define
`Gamma=bar e-pi c bar v`. Then the following is an EXACT finite
necessary defect:

```math
\boxed{\quad
\Delta\ge\left[
\frac{(\Gamma_+)^2}{8\pi c\sqrt S}
-\frac{\sqrt S}{2\sqrt n}\right]_+.\quad}                 (6)
```

The zero matrix case, if it occurs, is interpreted separately as a
vacuous inequality. For this particular cubic on a hollow signing,
`S>0`: simultaneous `P(L)=P(-L)=0` would require both the even and odd
parts of `P` to vanish on every eigenvalue, but they have no common
real root.

Proof. Apply (2) to `(L,R_+)` and `(-L,R_-)`. Their mean values of `v,e`
are precisely (3). Bound `(v_+ +b)^(-1),(v_-+b)^(-1)<=1/b` and use
Cauchy--Schwarz on the two diagonal-error terms. This gives

```math
\Gamma\le\pi c b+\frac{2\Delta\sqrt S+S/\sqrt n}{b}
\quad\hbox{for every }b>0.
```

The minimum of the right side is twice the geometric mean. Squaring
when `Gamma>0` proves (6); otherwise (6) is vacuous. The use of both
energy signs is essential for removing odd GLOBAL moments.

## 3. Original-problem structural consequence

Suppose a sequence of ACTUAL signings has
`Q(A_n)/n^(3/2)<=1/2+o(1)` and its even global moments through degree
12 tend to the centered variance-one semicircle moments. Then

```math
\bar v\longrightarrow\frac{178}{25},\qquad
\bar e\longrightarrow\frac{288}{25},\qquad
S\longrightarrow\frac{314747}{625}<23^2.
```

Using `pi<22/7`, the limiting gap is at least
`288/25-(11/7)(178/25)=58/175`. Equation (6) implies

```math
\liminf\Delta\ge
\frac{(58/175)^2}{8(11/7)23}
=\frac{841}{2213750}>\frac1{3000}.                         (7)
```

Thus an asymptotically minimizing family cannot have BOTH semicircle
even spectral moments through degree 12 AND vanishing variation of the
two explicit rooted-walk combinations in (5). No odd global spectral
moment hypothesis is needed. Either its even spectrum is detectably
non-semicircular or these local rooted-walk statistics stay nonuniform.

This is compatible with all currently known strict-subhalf signings:
neither parent optimality nor a cap bound forces their even spectral
moments to be semicircular. In particular (7) does NOT prove that good
exceptional selectors cannot exist.

The clipping inequality also upgrades the earlier one-sided cubic
semicircle witness, under its stated local moment conditions, from
`576/(361 pi)` to

```math
\frac{\bar e}{\pi\bar v}
=\frac{144}{89\pi}\simeq0.5150182428.
```

The added `I/10` is no longer needed: first take an order limit at
fixed `b>0` in (1), and then let `b` decrease to zero. The general
fixed-degree hierarchy already reaches `2/pi`; this cubic improvement
is useful because it needs only degrees through 12.

## 4. A completely explicit finite screening rule

Let `m_k` be the semicircle moments, and suppose `n>=10^6` and

```math
\left[\frac1n\sum_i((L^k)_{ii}-m_k)^2\right]^{1/2}
\le10^{-5}\qquad(0\le k\le12).                            (8)
```

Then

```math
\frac{Q(A)}{n^{3/2}}
\ge\frac{145071245497770289701}{288781684033417350314}
>0.50235>\frac12.                                         (9)
```

This bound uses only rational arithmetic and `pi<22/7`. To verify it,
write `R=P(L)^2`. The coefficient absolute sums of `P^2` and `P^4`
are respectively `B=276/25` and `B_2=53896/625`. With `eta=10^-5`,
(8), Jensen, and the triangle inequality imply

```math
178/25-B\eta\le v\le178/25+B\eta,\quad
e\ge288/25-B\eta,\quad\delta\le B\eta,
\quad s\le314747/625+B_2\eta<23^2.
```

Insert these values in (1) with `b=1/10` and `sqrt(n)>=1000`.
For the positive term use the upper bound on `v+b`; for the negative
term use the lower bound. The exact result is (9).

Therefore EVERY actual order-`n>=10^6` signing of normalized cap at
most one half violates at least one of the thirteen explicit tests
(8). This gives a finite check on any proposed near-optimal seed or
exceptional selector. It is not an asymptotic existence theorem.

## 5. Scope relative to selector counting and parent optimality

The bounded-cap parent theorem proves that uniformly chosen sparse
selectors satisfy the local moment conditions in probability. The
fixed-small-retention quantitative counting extension is handled in
the adversary track; it is not assumed in the finite argument above.

Combining rarity with (6) explains how an exceptional good selector
must fail the typical polynomial witness. It gives no contradiction:
a family comprising a vanishing fraction of all selectors can still
contain a good signing at every desired order. No current inequality
from parent minimality bounds the deterministic quantities in (3)--(5)
for those exceptional selectors. That missing implication is recorded
explicitly rather than replaced by a new putative transfer state.

## Preserved quantitative side calculation, not promoted to a theorem

Before this finite-defect direction was selected, a growing-operator
version of the graph proof was considered. If the normalized parent
operator bound is `C_D`, retaining its powers suggests that
`(n/D) C_D^14 ->0`, together with `C_D/sqrt(D)->0`, should suffice for
the fixed cubic moment hypotheses. For a fixed tensor seed this would
give the exponent range `alpha<1-14 log_d(||B||op/sqrt(d))`.
This is only a draft rate: denominator, odd-moment, and rooted-diagram
uniformity have NOT all been audited. No subsequent theorem here uses
it. It is preserved to avoid losing the unfinished calculation.

The exact rational constants and direct finite covariance checks are
replayed by
`computations/transfer_seed_exceptional_selector_verify_2026_09_06.py`.
