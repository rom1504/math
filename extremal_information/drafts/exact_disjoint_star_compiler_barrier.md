# Exact disjoint star compilers: an exact construction and a fan-in barrier

## Status and scope

This note studies a deliberately restricted but natural way of replacing a
same-support quadratic overlay by a disjoint future.  The future variables
are independent stars: there are no interactions among them.  This class is
wide enough to contain the obvious edge-by-edge compiler and to allow
arbitrary cancellation among all higher Fourier levels.  It is not the class
of all quadratic extensions.  In particular, a dense interacting auxiliary
block can correlate the optimizing future spins and is not covered.

The main result is stronger than the elementary statement “one edge variable
per old edge”: an exact star compiler with only `O(k)` future variables must
have fan-in `Omega(k^(2/3))`, every such compiler uses
`Omega(k^(3/2))` old--new incidences, and the all-positive cut shell needs at
least `k-1` future variables regardless of fan-in.  Thus bounded-fan-in
compilation has a vanishing `k^(3/2)` signal at the enlarged total scale.

Throughout, `x in {+-1}^k`, `E=binom(k,2)`, and

```math
H_T(x)=\sum_{i<j}T_{ij}x_ix_j,
\qquad
C_T(x)=E-H_T(x).
```

Constants in an optimized energy are immaterial, but retaining `E` makes the
exact compiler identity transparent.

## 1. The sparse edge compiler

### Proposition SC.1 (one-sided exact compilation)

For each old edge `e={i,j}` introduce an independent future spin `y_e` and
put

```math
G_T(x,y)=\sum_{i<j}y_{ij}(x_i-T_{ij}x_j).             \tag{SC.1}
```

Every used coefficient is exactly `+-1`, the future is disjoint from the old
edge support, and

```math
\max_yG_T(x,y)=C_T(x)                                \tag{SC.2}
```

for every `x`.  Consequently, for an arbitrary old energy `H_A`,

```math
\max_{x,y}\{H_A(x)+G_T(x,y)\}
=E+\max_x\{H_A(x)-H_T(x)\}.                          \tag{SC.3}
```

#### Proof

The `y_e` optimize independently and

```math
\max_{y_e=+-1}y_e(x_i-T_{ij}x_j)
=|x_i-T_{ij}x_j|=1-T_{ij}x_ix_j.
```

Summing proves (SC.2), then (SC.3). `square`

This is a one-sided response identity.  It does not by itself compile
`max_x|H_A-H_T|` into one future, because the negative channel contains
`-H_A`; two signed old channels would suffice, but switching the already
present old energy is itself part of the closure problem.

The construction has `m=E` future variables and `N=k+E=Theta(k^2)` total
variables.  A signal of order `k^(3/2)` is therefore only
`Theta(N^(3/4))`, hence is `Theta(N^(-3/4))` after normalization by
`N^(3/2)`.

## 2. Independent-star compilers

An **exact sign star compiler** for `C_T` is an identity

```math
C_T(x)=c+\max_{y in {+-1}^m}
       \sum_{a=1}^m y_a\sum_{i in S_a}\sigma_{ai}x_i
=c+\sum_{a=1}^m\left|\sum_{i in S_a}\sigma_{ai}x_i\right|,   \tag{SC.4}
```

where `sigma_ai in {+-1}` and `d_a=|S_a|`.  We allow all Fourier terms of
degree at least four from different stars to cancel; no atomwise
quadraticity is assumed.  Let `Delta=max_a d_a`.

For `d>=2` define

```math
\gamma_d=
2^{-(d-2)}{d-2\choose \lfloor(d-1)/2\rfloor}.       \tag{SC.5}
```

### Lemma SC.2 (the pair-Fourier Gram identity)

For `f_(S,sigma)(x)=|sum_(i in S)sigma_i x_i|` and distinct `i,j in S`,

```math
\widehat f_(S,sigma)({i,j})
=\gamma_|S|\,\sigma_i\sigma_j.                      \tag{SC.6}
```

It is zero if either index is outside `S`.  Moreover

```math
0<\gamma_d\le {1\over\sqrt{d-1}},
\qquad
\gamma_d{d\choose2}\le {1\over2}d^{3/2}.           \tag{SC.7}
```

#### Proof

After replacing `x_i` by `sigma_i x_i`, condition on the sum `Z` of the
other `d-2` Rademachers.  The conditional pair coefficient is

```math
{ |Z+2|-2|Z|+|Z-2| \over4}.                         \tag{SC.8}
```

If `d` is even this is one exactly when `Z=0` and zero otherwise.  If `d`
is odd it is `1/2` exactly when `|Z|=1` and zero otherwise.  Both cases give
(SC.5)--(SC.6).  The standard central-binomial estimate
`2^(-n) binom(n,floor(n/2))<=1/sqrt(n+1)` gives (SC.7). `square`

### Theorem SC.3 (fan-in, incidence, and Gram-rank barriers)

Every exact sign star compiler (SC.4) satisfies

```math
\sum_{a=1}^m\gamma_(d_a){d_a\choose2}\ge E,
\qquad
\sum_{a=1}^m d_a^{3/2}\ge k(k-1).                   \tag{SC.9}
```

Consequently

```math
m\ge {k(k-1)\over\Delta^{3/2}},
\qquad
\sum_a d_a\ge {k(k-1)\over\sqrt\Delta}.            \tag{SC.10}
```

In particular, `m=O(k)` forces `Delta=Omega(k^(2/3))`, while every compiler,
even with `Delta=k`, uses `Omega(k^(3/2))` old--new incidences.

There is also a spectral obstruction.  Let `K_T` range over positive
semidefinite matrices whose off-diagonal entries are `-T_ij`, and define

```math
g_+(T)=\min\{\operatorname{rank}K_T:K_T\succeq0,
                         (K_T)_{ij}=-T_{ij}\ (i\ne j)\}.     \tag{SC.11}
```

Then every star compiler has

```math
m\ge g_+(T).                                        \tag{SC.12}
```

For the all-positive signing `T_ij=1`, this gives the fan-in-independent bound

```math
m\ge k-1.                                           \tag{SC.13}
```

#### Proof

The target has pair Fourier coefficients `-T_ij`, all of absolute value
one.  By Lemma SC.2 and the triangle inequality,

```math
E=\sum_{i<j}|\widehat C_T({i,j})|
 \le\sum_a\gamma_(d_a){d_a\choose2}.
```

Equation (SC.7) proves (SC.9).  Since `d_a^(3/2)<=Delta^(3/2)` and
`d_a^(3/2)<=sqrt(Delta)d_a`, (SC.10) follows.

For the rank statement, let `v_a in R^k` have coordinates `sigma_ai` on
`S_a` and zero elsewhere.  Then

```math
K=\sum_a\gamma_(d_a)v_av_a^T                         \tag{SC.14}
```

is positive semidefinite, has rank at most `m`, and has off-diagonal entries
`-T_ij` by exact equality in (SC.4).  This proves (SC.12).  When `T_ij=1`,
every feasible completion is

```math
K=\operatorname{diag}(p_1,\ldots,p_k)-\mathbf1\mathbf1^T
\quad(p_i>0).
```

If `Kz=0`, then `p_i z_i=sum_j z_j` for every `i`; hence its nullspace has
dimension at most one.  Thus `rank K>=k-1`, proving (SC.13). `square`

### Corollary SC.4 (the scale-retention threshold)

If `Delta=o(k^(2/3))`, every exact sign star compiler has `N=k+m=omega(k)`.
For a `Theta(k^(3/2))` old response signal,

```math
{k^{3/2}\over N^{3/2}}=o(1).                        \tag{SC.15}
```

More quantitatively, writing `Delta=k^(2/3)/g(k)` with `g(k)->infinity`,
(SC.10) gives `N>=Omega(k g(k)^(3/2))`, so the normalized signal is
`O(g(k)^(-9/4))`.

Endpoint-local stars (`d_a<=2`) obey the stronger exact count `m>=E`:
each star has only one possible nonzero pair coefficient.  Proposition SC.1
attains equality.  Hence “one auxiliary per edge” is necessary and sufficient
inside the endpoint-local class, but not claimed outside it.

### Lemma SC.4a (the absolute Rademacher sum has every even Fourier level)

Let

```math
f_k(x)=\left|\sum_(i=1)^k x_i\right|,
\qquad alpha_(k,s)=\widehat f_k(S)\quad(|S|=s).       \tag{SC.15a}
```

Then `alpha_(k,s)!=0` for every even `s`, `2<=s<=k`.  In particular
`alpha_(k,2)=gamma_k>0`.

#### Proof

For `k=2n`, let `I_0` be the middle-slice indicator.  The cube Laplacian,
whose level-`s` eigenvalue is `2s`, obeys the pointwise identity

```math
L f_(2n)=2f_(2n)-4nI_0.                              \tag{SC.15b}
```

For `s=2r`, the middle-slice coefficient is

```math
\widehat I_0(S)=2^(-2n)(-1)^r
 { {2n\choose n}{n\choose r}\over {2n\choose2r}},  \tag{SC.15c}
```

by the middle Krawtchouk identity.  It is nonzero, and comparison of Fourier
coefficients in (SC.15b) gives

```math
alpha_(2n,2r)=-{2n\over2r-1}\widehat I_0(S).         \tag{SC.15d}
```

For `k=2n+1`, average over the last coordinate.  The resulting function on
`2n` variables is `f_(2n)+I_0`, so

```math
alpha_(2n+1,2r)
=-{2n-2r+1\over2r-1}\widehat I_0(S),                \tag{SC.15e}
```

again nonzero.  The level-two value agrees with (SC.5). `square`

### Theorem SC.4b (fully dense independent stars still need quadratic order)

Suppose every star in (SC.4) has full support `S_a=[k]`.  Exact compilation
of `C_T` then forces

```math
m\gamma_k\ge\max_x H_T(x).                          \tag{SC.15f}
```

Consequently, if both orientations `C_T` and `C_(-T)` are compiled with at
most `m` stars each, then

```math
m\ge {Q(T)\over\gamma_k}=Omega(k^2),                \tag{SC.15g}
```

where the last bound holds for every complete signing by the universal
`Q(T)=Omega(k^(3/2))` bound and `gamma_k=Theta(k^(-1/2))`.

#### Proof

Let `mu` be the empirical probability law of the star sign vectors and
symmetrize it under `sigma->-sigma`; this changes no absolute linear form.
Lemma SC.4a and the absence of target Fourier levels `4,6,...` imply

```math
E_mu\prod_(i in S)sigma_i=0\quad
(|S|>=4\hbox{ even}),                               \tag{SC.15h}
```

while level two gives

```math
E_mu sigma_i sigma_j=-{T_(ij)\over m\gamma_k}.       \tag{SC.15i}
```

All odd moments vanish by symmetrization.  Fourier inversion of the
probability mass function is therefore exact:

```math
mu(sigma)=2^(-k)\left(1-{H_T(sigma)\over m\gamma_k}\right).
                                                               \tag{SC.15j}
```

Nonnegativity for every `sigma` proves (SC.15f).  Applying it to `T` and
`-T` gives `m gamma_k>=Q(T)`.  The central-binomial estimate in (SC.5) and
Lemma SC.6 give (SC.15g). `square`

Thus neither endpoint-local nor complete-bipartite independent selectors
provide a linear-order scale-preserving compiler.  Any such escape must use
nontrivial auxiliary interaction, mixed supports in a regime not covered by
the endpoint/dense extremes, or a weaker response notion.

## 3. An interacting-selector lower bound

The independent-star hypothesis can be removed at the price of specializing
the target to the all-positive cut shell.  An arbitrary quadratic future on
`m` auxiliary spins has the form

```math
F(x)=\max_{y in {+-1}^m}\{c(y)+b(y)\mathbin\cdot x\},       \tag{SC.16}
```

where auxiliary--auxiliary interactions are absorbed into `c(y)`.  If the
old--new coefficients are unit signs or zero and old vertex `i` has future
degree `d_i`, then `|b_i(y)|<=d_i` for every selector `y`.

### Theorem SC.5 (interacting selector-covering barrier)

Let `k` be even and

```math
C_+(x)=E-\sum_{i<j}x_ix_j
      ={k^2-(\sum_i x_i)^2\over2}.                  \tag{SC.17}
```

Suppose (SC.16) obeys

```math
\max_x|F(x)-C_+(x)|\le\eta<{k^2\over4}.             \tag{SC.18}
```

Put `D^2=sum_i d_i^2` and `a=k^2/2-2eta`.  Then

```math
\sum_i d_i\ge a,                                    \tag{SC.19}
```

and, provided `D>0`,

```math
m\log2\ge {a^2\over2D^2}-\log(k+1).                 \tag{SC.20}
```

In particular, for a unit complete bipartite old--new block, `d_i=m`, so
every `o(k^2)`-accurate compiler satisfies

```math
m\ge (8\log2)^(-1/3)k-o(k)
  =0.565\ldots k-o(k).                               \tag{SC.21}
```

More quantitatively, if `m=ck+o(k)` with
`c<(8log2)^(-1/3)`, its asymptotic uniform distortion is bounded below by

```math
\liminf_{k\to\infty}{\eta\over k^2}
\ge {1-\sqrt{8c^3\log2}\over4}>0.                   \tag{SC.22}
```

#### Proof

Fix a balanced `x`, so `C_+(x)=k^2/2`, and choose a selector `y` active at
`x`.  Its affine piece lies below `F` everywhere.  Comparing it first with
the all-positive spin and then with the all-negative spin gives

```math
b(y)\mathbin\cdot(x-\mathbf1)\ge a,
\qquad
b(y)\mathbin\cdot(x+\mathbf1)\ge a.
```

Adding shows

```math
x\mathbin\cdot b(y)\ge a.                           \tag{SC.23}
```

Since `|b_i(y)|<=d_i`, (SC.19) follows.

For a fixed selector, Hoeffding's Rademacher bound gives

```math
#\{x in {+-1}^k:x\mathbin\cdot b(y)\ge a\}
\le2^k\exp\{-a^2/(2D^2)\}.                          \tag{SC.24}
```

The at most `2^m` selectors must cover every balanced `x`.  Since
`binom(k,k/2)>=2^k/(k+1)`, a union bound proves (SC.20).  With `D^2=km^2`,
(SC.20) and `eta=o(k^2)` give
`m^3 log2>=k^3/8-o(k^3)`, proving (SC.21).  Solving the same inequality for
`eta` proves (SC.22). `square`

### Lemma SC.5a (antipodal exposed-set selector bound)

The preceding covering argument has a query-relative form.  Let

```math
F(x)=\max_(q\in[K])\{c_q+b_q\mathbin\cdot x\},
\qquad ||b_q||_2<=D,                                 \tag{SC.22a}
```

and suppose `||F-f||_infinity<=eta`.  If there are `p in {+-1}^k`, a set
`X subseteq {+-1}^k`, and `a>0` such that

```math
f(x)>=max\{f(p),f(-p)\}+a+2eta\quad(x\in X),         \tag{SC.22b}
```

then

```math
\log K>=\log|X|-k\log2+{a^2\over2D^2}.              \tag{SC.22c}
```

If the slopes are `b_q=By_q`, `y_q in {+-1}^m`, one may take
`K<=2^m` and `D<=||B||_(2->2)sqrt(m)`.  Thus positive entropy of an
antipodally exposed near-top face converts a spectral bridge bound into a
linear auxiliary-state lower bound at the `k^(3/2)` scale.

Indeed, an affine piece active at `x in X` lies below the envelope at both
`p` and `-p`.  The two comparisons add to `b_q dot x>=a`.  Hoeffding says
one slope covers at most `2^k exp(-a^2/(2D^2))` cube points, and the `K`
pieces must cover `X`. `square`

Theorem SC.5 is the nearly-full-entropy instance `p=1`, `X` the balanced
slice, and `a=k^2/2-2eta`.  Lemma SC.5a also identifies the exact extra input
needed to transfer the lower bound to a flat signing: a macroscopic exposed
energy gap carried by exponentially many configurations.  A maximum value
alone does not supply that entropy.

### Corollary SC.5b (bounded-cap linear lifts cannot be universal)

For the general selector representation (SC.16), let
`b(y)=By` for a fixed old--new interaction matrix `B`.  Then

```math
osc_(x in {+-1}^k)F(x)
<=2||B||_(infinity->1).                              \tag{SC.23a}
```

If this bridge is part of a quadratic parent `P` and `Q(P)<=CN^(3/2)`, then

```math
||B||_(infinity->1)<=Q(P),
\qquad osc F<=2CN^(3/2).                             \tag{SC.23b}
```

Consequently, when `N=O(k)`, no bounded-cap exact-sign lift can approximate
the all-positive cut shell `C_+`, whose oscillation is `k^2/2`, to
`o(k^2)` uniformly.

Indeed, maxima of affine pieces are Lipschitz in their linear parts, giving
(SC.23a).  For fixed old and auxiliary spins, flipping the whole auxiliary
shore reverses the cross term and preserves both internal terms; one of the
two parent energies has absolute value at least that cross term.  This proves
(SC.23b). `square`

Unlike Theorem SC.3, this theorem permits arbitrary interaction and
correlation among the auxiliary selectors.  It is nevertheless a universal
cut-shell benchmark, not yet a lower bound tailored to the short-seed
alternating-form contexts.  It shows that any genuinely sublinear-vertex
exact-sign future fails by a quadratic, not merely `k^(3/2)`, amount; a
possible linear-overhead escape must use at least `0.565...k` auxiliary
spins and quadratic old--new incidence.  If the completed parent must also
have cap `O(N^(3/2))`, even that escape is impossible for the universal
all-positive shell.  A bounded-cap compiler must exploit a target family
whose own response oscillation is only `O(k^(3/2))`, such as the flat Gram
family, rather than simulate every signing.

## 4. Exact completion cannot be perturbative after quadratic blow-up

### Lemma SC.6 (an elementary complete-sign cap)

Every hollow complete signing `S` on `N` vertices satisfies

```math
Q(S)\ge c_0N^{3/2}                                  \tag{SC.25}
```

for an absolute `c_0>0` and all sufficiently large `N` (one may take any
fixed `c_0<2/(3sqrt(6))`).

#### Proof

Partition the vertices into `U,V` with `|U|=floor(N/3)` and randomize the
spins on `U`.  The sharp `p=1` Khintchine inequality, or its elementary
constant `1/sqrt2`, gives a choice for which

```math
\sum_{v in V}\left|\sum_{u in U}S_{uv}x_u\right|
\ge |V|\sqrt{|U|/2}.
```

Choose each `x_v` to realize this cross value.  Flipping all spins in `V`
changes the sign of the cross term and leaves both internal terms fixed, so
one of the two full quadratic energies has absolute value at least the cross
value.  Optimizing the split asymptotically at `|U|=N/3` proves (SC.25).
`square`

Thus taking the endpoint-local compiler, whose total size is `N=Theta(k^2)`,
and filling all absent pairs by arbitrary signs necessarily creates a
complete signing with cap `Omega(k^3)`.  The sparse compiled signal is only
`Theta(k^(3/2))`.  Consequently **uniform-norm perturbative completion** of
this quadratic-blow-up construction cannot preserve the signal scale.

This is not a no-go theorem for contextual cancellation atop a large common
background: response differences can sometimes survive even when every
individual cap is large.  Nor does it rule out an `O(k)`-vertex interacting
auxiliary compiler.

## 5. What remains open

For a general disjoint quadratic future

```math
G_T(x,y)=x^TB_Ty+H_(C_T)(y),                         \tag{SC.26}
```

the optimizing `y` need not factor into stars.  The effective old-spin
function is a maximum of correlated affine selectors, and the Gram identity
(SC.14) no longer applies term by term.  Proving or refuting an `O(k)`-vertex
exact-sign realization of `C_T`, with a bounded-cap parent and a single
absolute-response channel, remains the genuine closure question.

The present theorem nevertheless isolates two necessary resources for any
escape from the obvious compiler:

1. linearly many independent stars require fan-in at least `k^(2/3)` and
   superlinear old--new incidence;
2. to beat this within lower fan-in, auxiliary--auxiliary interaction must
   create cancellation before the independent absolute values are paid.

That is precisely the joint-cancellation resource absent from scalar-channel
and independently-paid decompositions.
