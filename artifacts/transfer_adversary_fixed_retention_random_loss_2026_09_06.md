# A quantitative macroscopic random-selector obstruction

Date: 2026-09-06. This extends the vanishing-retention random-selector
theorem to a sufficiently small FIXED positive retention, with an
elementary polynomial probability bound. No exponential rarity rate is
claimed. All cap statements refer to actual hollow signings.

## Theorem

Fix `C0<infinity` and `c<2/pi`. There exist `rho>0`, `K<infinity`,
and `n0`, depending only on `C0,c`, such that for every symmetric hollow
sign matrix `A_D` with `Q(A_D)<=C0 D^(3/2)`, and every integer
`n0<=n<=rho D`, a uniform principal `n`-subset `T` satisfies

```math
\Pr\{Q((A_D)_T)<c n^{3/2}\}\le\frac K n.                 (1)
```

Consequently at most `(K/n) binom(D,n)` selectors can have that low cap.
The constants are finite and computable from fixed walk-pattern lists,
integer polynomial coefficients, and the explicit spectral-core bound.
They are not numerically optimized here. In particular (1) applies to
any fixed retention in `(0,rho]` as the order increases.

This does NOT exclude an exceptional selector. Even a proportion `K/n`
can contain exponentially many subsets at fixed retention, and the
proof supplies no criterion showing that none of them is a good signing.

## 1. Uniform covariance control for a fixed connected graph observable

Let `B` be a full symmetric sign matrix of order `D`, with
`||B||op<=C sqrt D`, where `C>=1`. Let `G` be a fixed connected
Eulerian multigraph, allowing loops, with `v` vertices and `k` edge
occurrences. The case of one isolated vertex and `k=0` is also allowed.
For a uniform `n`-subset `T`, define

```math
X_G(T)=n^{-1-k/2}
 \sum_{\phi:V(G)\hookrightarrow T}
                 \prod_{uv\in E(G)}B_{\phi(u),\phi(v)}.   (2)
```

For fixed `r`, if both graphs have at most `r` edges and
`n>=2(2r+2)^2`, then

```math
|\operatorname{Cov}(X_G,X_{G'})|\le\frac{V_r(C)}n,         (3)
```

uniformly for `n<=D`. One admissible explicit constant is obtained by
putting

```math
B_r(C)=2(r+1)!C^r,\quad v_*=2r+2,
```

and taking

```math
V_r(C)=3v_*^2B_r(C)^2+
       v_*!(1+2^{v_*})B_{2r}(C).                         (4)
```

Only the previously audited unrestricted Eulerian graph estimate is
used: for a parity graph `F`,
`|t_F(B)|<=C^(|E(F)|) D^(-b(F))`, with
`b=|V|-components-|E|/2`. It is the precise Mingo--Speicher specialization
in `transfer_seed_sparse_bounded_cap_parent_theorem_2026_09_06.md`.

### 1.1 A normalized injective graph bound

Write `t_G^inj(B)` for the mean over injective ambient colorings; edge
multiplicities are still reduced modulo two inside the sign product.
For a connected Eulerian graph, an Euler circuit makes it a closed-walk
pattern, so `a=v-1-k/2<=b(F)`. For each partition identifying `h`
vertices, its quotient remains connected and Eulerian, and its exponent
is `a-h`. Thus

```math
|n^{a-h}t_{G/\pi}(B)|\le C^k\qquad(1<=n<=D).             (5)
```

For `a-h<=0` use the trivial sign-product bound; otherwise use the
operator graph bound and `n<=D`. Exact partition-Mobius inversion now
gives, for `D>=2v^2`,

```math
|n^a t_G^{\rm inj}(B)|\le2v!C^k\le B_k(C).               (6)
```

Indeed `D^v/(D)_v<=2`, each nontrivial quotient has the additional factor
`(n/D)^h<=1`, and the sum of absolute partition-Mobius coefficients is
`v!` (group permutations by their cycle partitions).

### 1.2 Actual sample overlaps gain one power of n

Let `G,G'` have `v_1,v_2` vertices and lengths `k_1,k_2`, and put
`v=v_1+v_2`, `a_i=v_i-1-k_i/2`. In the product of the two sums (2),
the two individually injective colorings can share `z` vertices. Such
an overlap is a partial matching between their two vertex sets; there
are at most `v!` partial matchings, since each is a permutation made
of disjoint cross-set transpositions.

If `z>=1`, the resulting graph `G_sigma` is connected and Eulerian,
with exponent `a_sigma=a_1+a_2-z+1`. Its contribution is

```math
\frac{(n)_{v-z}}{n^{2+(k_1+k_2)/2}}t_{G_\sigma}^{\rm inj}(B).
```

Equation (6) bounds its absolute value by `B_(2r)(C)/n`. The extra
`1/n` is from merging two original connected components, not from an
incorrect assumption that the two walks become independent.

### 1.3 The disjoint contribution versus the product of expectations

Let `t_disj` denote the injective average on the disjoint union. Counting
overlaps of two independent injective AMBIENT colorings gives exactly

```math
(D)_{v_1}(D)_{v_2}t_G^{\rm inj}t_{G'}^{\rm inj}
=(D)_v t_{\rm disj}
 +\sum_{\sigma:\,z\ge1}(D)_{v-z}t_{G_\sigma}^{\rm inj}.   (7)
```

The difference between the two coefficients multiplying the product of
normalized averages is at most `3v^2/n`: the sample cross-collision probability is
at most `v_1v_2/n`, and
`(D)_(v1)(D)_(v2)/(D)_v-1<=2v_1v_2/D`. Equation (6) therefore bounds
this part of the covariance by `3v^2 B_r(C)^2/n`.

Each remaining term of (7), after its sample normalization, is bounded
using

```math
\frac{(D)_{v-z}}{(D)_v}\le2^zD^{-z},\qquad
n^{a_1+a_2}D^{-z}t_{G_\sigma}^{\rm inj}
=(n/D)^z n^{a_\sigma-1}t_{G_\sigma}^{\rm inj}.
```

Its absolute value is at most `2^v B_(2r)(C)/n`. Combining these terms
with the actual sample overlaps proves (3)--(4). This proof never
approximates a whole sample by independent distinct colors.

## 2. Four polynomial statistics concentrate at fixed retention

Assume first that a hollow parent obeys `||A_D||op<=C0 sqrt D` and
complete it by `B=A_D+I`, of normalized norm at most `C=C0+1`.
Set `L=(A_D)_T/sqrt n`. Fix a finite degree `q>=1`, and use

```math
P_q=\sum_{j=0}^q U_j(x/2),\quad R=I+P_q(L)^2,\qquad
v=\frac{\operatorname{tr}R}n,\quad
w=\frac{\operatorname{tr}(LR)}n,\quad
s=\frac{\operatorname{tr}(R^2)}n,\quad
y=\frac1n\sum_iR_{ii}^2.                                (8)
```

Every trace is a finite linear combination of connected closed-walk
observables (2). A summand of `y` consists of two closed walks with a
common root, so its union is also connected Eulerian. An original hollow
loop contributes zero and is simply omitted; loops created in the
ambient injectivity calculation use the full-sign completion.

There are finitely many such patterns, with at most `4q` edge
occurrences. Formula (3) and the fixed polynomial coefficients give

```math
\max\{\operatorname{Var}v,\operatorname{Var}w,
       \operatorname{Var}s,\operatorname{Var}y\}
\le K_{q,C}/n.                                          (9)
```

The constant is explicitly bounded by the square of each coefficient
absolute sum times its number of patterns squared times `V_(4q)(C)`.
A length-`r` rooted walk or rooted two-walk pattern has at most
`(r+1)^(r+1)` canonical label patterns, so all constants here are
computable without an asymptotic oracle.

The already-proved single/rooted-double moment expansions also give
a finite constant `J_(q,C)` such that

```math
\max\{|\mathbb Ev-v_0|,|\mathbb Ew-w_0|,
       |\mathbb Es-s_0|,|\mathbb Ey-v_0^2|\}
\le J_{q,C}(\sqrt p+n^{-1/2}+D^{-1/2}),                  (10)
```

where `p=n/D`, `v_0=q+2`, `w_0=2q`, and
`s_0=E{(1+P_q(S)^2)^2}` for the variance-one semicircle law. The square-root
bound accommodates odd-length terms; even combined lengths have the
stronger earlier bound. The same finite pattern/Mobius estimates used
in Section1 provide a computable `J_(q,C)`.

## 3. Turning concentration into a fixed-p cap bound

Choose `q` so that `d_q=2q/[pi(q+2)]>c`, and set `Delta=d_q-c`.
Put `K0=(2v_0+2)(s_0+1)`, and choose

```math
0<eta<=min(1/4,pi Delta/16,(pi Delta)^2/(64K0)).
```

On the event that the four quantities in (8) are within `eta` of
their four targets in (10), one has

```math
\frac wv\ge\frac{w_0}{v_0}-2\eta,\qquad
\delta_R^2:=y-v^2\le(2v_0+2)\eta,\qquad s\le s_0+1.
```

The finite Gaussian inequality with diagonal lower bound one gives

```math
Q((A_D)_T)/n^{3/2}
\ge d_q-\frac1\pi
 \{2\eta+\sqrt{K0\eta}+(s_0+1)/\sqrt n\}>c              (11)
```

for sufficiently large `n`. Choose a fixed `p0>0` small and `n0` large
so that the right side of (10) is at most `eta/2` whenever
`n>=n0` and `n/D<=p0`. For example its two parts are controlled by
`p0<=(eta/(4J_(q,C)))^2` and `n0>=(8J_(q,C)/eta)^2`.
Increase `n0` also to satisfy the graph-size conditions and
`(s_0+1)/sqrt(n0)<=pi Delta/4`.

Chebyshev and (9), applied to all four statistics, now give

```math
\Pr\{Q((A_D)_T)<c n^{3/2}\}
\le16K_{q,C}/(eta^2 n).                                 (12)
```

This is an actual fixed-retention probability bound, uniform in the
parent within its operator-norm class. It is not an exponential estimate.

## 4. Actual bounded-cap parents and the spectral-core loss

Return to the cap hypothesis of the theorem. Choose fixed
`epsilon in (0,1/4)` and `c1<2/pi` so that
`c1(1-2epsilon)^(3/2)>c`. Principal spectral deletion supplies a
deterministic core `R_D`, of order at least `(1-epsilon)D`, satisfying

```math
\|(A_D)_{R_D}\|_{op}
\le\frac{4K_G C0}{epsilon\sqrt{1-epsilon}}\sqrt{|R_D|}.
```

This is a fixed normalized operator bound. Apply Section3 to this
class and threshold `c1`, obtaining `p0,K1,n1`. If
`rho=(1-epsilon)p0` and `n<=rho D`, then every intersection size
`m=|T intersect R_D|` satisfies `m/|R_D|<=p0`.

The hypergeometric variance bound gives

```math
\Pr\{m<(1-2epsilon)n\}\le1/(4epsilon^2 n).
```

Conditional on `m`, the intersection is uniform in the core. For
`n>=n1/(1-2epsilon)`, (12) has conditional failure probability at
most `K1/[(1-2epsilon)n]` on the likely size range. Principal cap
monotonicity and the preserved factor `(m/n)^(3/2)` now prove (1), with

```math
K=1/(4epsilon^2)+K1/(1-2epsilon).
```

## 5. Deterministic selectors are not excluded

For fixed positive `rho`, the fraction `K/n` leaves exponentially many
possible exceptional subsets. The proof gives no lower bound on the
best principal selector, no algorithm for locating an exceptional
selector, and no comparison between original minimizing sequences at
different orders. In particular it does not close the original
convergence/nonconvergence obligation.

This is stronger than a tiny-sample or vanishing-retention falsifier:
the lower failure probability is polynomially controlled on an entire
macroscopic interval `n0<=n<=rho D`, with `rho` depending only on the
actual parent cap bound and the target child coefficient.

Exact overlap and covariance decompositions are replayed by
`computations/transfer_adversary_fixed_retention_covariance_2026_09_06.py`.
