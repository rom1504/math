# Strict-subhalf signings can violate pointwise powered aggregation

Date: 2026-09-06. Fresh convergence track. Status: proved; Section 7 and
its imported sampling lemma independently audited PASS in
`transfer_adversary_stopped_asymmetric_power_audit_2026_09_06.md`.
This is NOT a result about
the minimum-value reverse-Fekete target, and its splits are highly
asymmetric.

Section 7 proves the stronger variant in which the small child is an
EXACT minimizing signing. It uses the director's stopped discounted
power identity, rather than prescribing a deliberately worse small cap.

Write `Q(A)=max_x |x^T A x|/2`, over Boolean signs, and
`c(A)=Q(A)/n^(3/2)` for an order-`n` hollow symmetric signing. We use the
verified all-order strict upper bound: for some fixed `U<1/2`,
`M_n<= (U+o(1))n^(3/2)`.

## 1. Exact theorem and the distinction from the convergence target

Fix `U<b<d<1/2`. For every sufficiently small fixed `s>0`, put `q=1-s`.
There are `epsilon_s,kappa_s>0`, orders `n_j->infinity`, and actual
hollow symmetric signings

```math
P_j=\begin{pmatrix}B_j&F_j\\F_j^T&D_j\end{pmatrix},
\qquad |B_j|=\lfloor qn_j\rfloor,
\quad |D_j|=n_j-\lfloor qn_j\rfloor,
```

such that ALL THREE normalized caps are uniformly below one half:

```math
\max\{c(P_j),c(B_j),c(D_j)\}\le\tfrac12-\epsilon_s,       (1)
```

but

```math
Q(B_j)^{2/3}+Q(D_j)^{2/3}-Q(P_j)^{2/3}
                         \ge\kappa_s n_j.                 (2)
```

In fact `c(D_j)->d`. The constants `s,epsilon_s,kappa_s` are fixed
BEFORE the order limit; no growing iteration depth is hidden in (2).

Thus pointwise powered aggregation with an `o(n)` defect can fail even
when parent and both children have genuine strict-subhalf Boolean caps.
The construction requires `s` sufficiently small; it does not show this
for comparable splits. Most importantly, (2) concerns the actual child
caps `Q(B_j),Q(D_j)`, not `M_|B_j|,M_|D_j|`. In particular it does NOT
refute

```math
M_{m+n}^{2/3}\ge M_m^{2/3}+M_n^{2/3}-C\sqrt{m+n}
\quad (1/2\le m/n\le2),                                  (RF)
```

the sufficient convergence criterion proved in
`transfer_fresh_reverse_fekete_2026_09_06.md`.

## 2. A primary sampling upper bound, in the exact original norm

There is an absolute constant `K_R` such that, whenever `T` is uniform
among the `h`-subsets of an order-`n` hollow symmetric signing `A`,
`h<=n/2`, and `h` is sufficiently large,

```math
\mathbb E Q(A_T)
\le K_R\left[(h/n)^2Q(A)+(h/n)^{3/2}n^{3/2}\right].        (3)
```

Here and below a principal cap is monotone: averaging independent
omitted spins proves `Q(A_S)<=Q(A_T)` for `S subset T`.

### Exact primary import

[Rudelson--Vershynin, *Sampling from large matrices: an approach through
geometric functional analysis*, arXiv:math/0503442, Theorem 1.5,
pp. 4 and 13--15](https://arxiv.org/pdf/math/0503442)
proves for Bernoulli retention `p` that

```math
\mathbb E\|A_R\|_C
\le K\{p^2\|A-\operatorname{diag}A\|_C
       +p\|\operatorname{diag}A\|_C
       +p^{3/2}(\|A\|_{\mathrm{Col}}+\|A^T\|_{\mathrm{Col}})\}.
                                                                    (4)
```

The cut norm takes the maximum absolute rectangular subset sum;
`Col` sums the Euclidean column lengths. This is a SAME-selector
principal restriction, not independent row and column sampling.

For `beta(A)=max_(x,y) |x^T A y|`, the norm comparisons are

```math
Q(A)\le\tfrac12\beta(A)\le2\|A\|_C,
\qquad \|A\|_C\le\beta(A)\le4Q(A).                       (5)
```

The last bound uses symmetric hollow polarization: for Boolean `x,y`,
put `u=(x+y)/2,v=(x-y)/2`. Then `x^T A y=u^T A u-v^T A v`, and each
quadratic expression has absolute value at most `2Q(A)` by multilinear
extension to the cube. The cut-norm comparisons also appear in the
primary proof's Section 4.1. Since `diag A=0` and
`||A||Col=||A^T||Col=n sqrt(n-1)`, (4)--(5) give

```math
\mathbb E Q(A_R)\le K'[p^2Q(A)+p^{3/2}n^{3/2}].           (6)
```

### Bernoulli to uniform exact size

Take Bernoulli retention `p=2h/n`, and on the event `E={|R|>=h}` choose
a uniform `h`-subset `T` of `R`. Conditional on `E`, this `T` is uniform
among all `h`-subsets of `[n]`, by permutation symmetry. Therefore

```math
\Pr(E)\,\mathbb E_{|T|=h}Q(A_T)
=\mathbb E[1_E Q(A_T)]\le\mathbb E Q(A_R).                (7)
```

The binomial mean is `2h` and its variance is at most `2h`, so Chebyshev
gives `Pr(E^c)<=2/h<=1/2` for `h>=4`. Inserting (6) with
`p=2h/n` proves (3). There is no assumption that conditioning a
Bernoulli selector on having exactly its mean preserves an expectation
inequality without a probability factor.

## 3. Dialing a small child to a prescribed cap

For every fixed `d>U` and all sufficiently large `h`, there exists an
order-`h` signing `D_h` with

```math
d h^{3/2}\le Q(D_h)<d h^{3/2}+2.                          (8)
```

Start from a signing with cap below `d h^(3/2)` and flip its negative
edges one at a time until reaching the all-positive signing. Its final
cap is `h(h-1)/2>d h^(3/2)`. Each edge flip changes every Boolean energy
by absolute value two, hence changes the cap by at most two. The first
crossing of the threshold proves (8). Cap monotonicity along the walk
is not needed.

This lemma deliberately chooses a child WORSE than `b`, though still
strictly below one half. It must not be read as a near-minimizer lemma.

## 4. A quantitative stopped thinning lemma

Fix `U<b<1/2`. There exist positive constants `a,rho_0,s_0`, independent
of `s`, such that, for every fixed `0<s<s_0`, there are arbitrarily
large orders `n` and actual signings `A` with `c=c(A)<b` for which a
uniform `floor((1-s)n)`-subset `S` obeys

```math
\Pr\{c(A_S)\ge c+a s\}\ge\rho_0.                         (9)
```

In addition every such child satisfies the deterministic upper bound

```math
c(A_S)\le c\,(n/|S|)^{3/2}\le b(1-s)^{-3/2}+o_n(1).      (10)
```

Proof. The uniform small-fixed-retention theorem in Section 2 of
`transfer_director_exponential_selector_cost_2026_09_06.md` supplies
`rho in (0,1)` such that, from EVERY parent of normalized cap at most
`1/2`, a uniform restriction to at most `rho` of its vertices and to
size tending to infinity has cap exceeding `1/2` with probability
tending to one. Set `L=log(1/rho)`, `K=L+2`, and

```math
r=\left\lceil\frac{L}{-\log(1-s)}\right\rceil+1
\le K/s\qquad(0<s<1).
```

Start from an all-order good signing of order `N` with normalized cap
`c_0<=U+o(1)`, and take `r` nested uniform restrictions with sizes
`n_0=N`, `n_(i+1)=floor((1-s)n_i)`. For fixed `s`, `r` is fixed and
all sizes tend to infinity. Each marginal is an exactly uniform subset
of the starting vertex set. Hence its terminal cap `C_r` is above
`1/2`, and therefore above `b`, with probability tending to one.

Let `tau` be the first index with `C_i>=b`, or `r` if none exists.
The identity

```math
C_\tau-C_0=\sum_{i=0}^{r-1}1_{\{\tau>i\}}(C_{i+1}-C_i)    (11)
```

is pointwise exact. All normalized caps are nonnegative. Since a hit
occurs with probability tending to one, writing `Delta=b-U>0` gives

```math
\mathbb E(C_\tau-C_0)\ge\Delta-o(1)\ge\Delta/2             (12)
```

for sufficiently large `N`. For some index `i`, the corresponding
summand in expectation is at least `Delta/(2r)`. Conditional on the
entire history up through `i`, the next selector is still uniform in
its parent. Consequently there is a deterministic parent realization
`A` before the stop, with `c(A)<b`, for which the mean normalized
increment `X=c(A_S)-c(A)` is at least `Delta/(2r)>=Delta s/(2K)`.

For sufficiently small fixed `s` and sufficiently large orders,
principal monotonicity gives the deterministic bound `X<=2s`; any
negative increments only help the next upper estimate. Set
`a=Delta/(4K)`. Then

```math
2as\le\mathbb EX\le as+2s\Pr\{X\ge as\},
```

so (9) holds with `rho_0=a/2`. A parent can be selected for each
sufficiently large starting order; its depth may depend on that order,
but its order is at least `floor((1-s)^r N)-O_s(1)` and tends to
infinity. This proves the lemma without a compactness-derived unknown
loss size. Equation (10) is exact principal monotonicity with the
rounding error displayed.

## 5. Cheap complement repair proves the theorem

Keep `b<d<1/2` fixed and put

```math
\gamma=d^{2/3}-b^{2/3}>0.
```

Choose any fixed positive `a_0` below the verified asymptotic universal
lower bound, for example `a_0=0.4`. All sufficiently large signings
then have normalized cap at least `a_0`. The derivative of `t^(2/3)`
on `[a_0,infinity)` is at most `L_0=(2/3)a_0^(-1/3)`. Put

```math
t=\gamma/(8L_0).
```

Take one of the arbitrarily large parents `A` from Section 4, write
`n=|A|`, and sample `S` uniformly with `|S|=floor((1-s)n)`. Its
complement `T` is an exactly uniform `h=n-|S|` subset, with
`h/n=s+O(1/n)`. By (3) and `c(A)<b`, there is a constant `K_1`,
independent of sufficiently small `s`, such that

```math
\mathbb E[Q(A_T)/n^{3/2}]\le K_1s^{3/2}+o_n(1).
```

Thus the probability that `Q(A_T)>t s n^(3/2)` is at most
`K_1 sqrt(s)/t+o_n(1)`. Choose `s` sufficiently small that this is
less than `rho_0/2` for large orders. The event in (9), of probability
at least `rho_0`, therefore overlaps the cheap-complement event.
Select one such `S,T` deterministically.

Replace ONLY the old `T`-principal block by the signing `D_h` from
(8), retaining the `S`-block and every bridge edge; call the resulting
order-`n` signing `P`. The triangle inequality for Boolean energies gives

```math
Q(P)\le Q(A)+Q(A_T)+Q(D_h).                               (13)
```

Choose `s` smaller if necessary so `d sqrt(s)<t`. Equations (8), (13)
then imply, for sufficiently large `n`,

```math
c(P)\le c+2ts+o_n(1),\qquad c=c(A)\in[a_0,b),
\quad c(A_S)\ge c+as,
\quad c(D_h)=d+o_n(1).                                   (14)
```

By concavity and the derivative bound,

```math
Q(P)^{2/3}/n\le c^{2/3}+2L_0ts+o_n(1)
                 =c^{2/3}+\gamma s/4+o_n(1).             (15)
```

The two child powers, using only `c(A_S)>=c`, satisfy

```math
\frac{Q(A_S)^{2/3}+Q(D_h)^{2/3}}n
\ge(1-s)c^{2/3}+s d^{2/3}+o_n(1)
\ge c^{2/3}+\gamma s+o_n(1).                             (16)
```

Subtracting proves (2) with `kappa_s=gamma s/2` eventually. Finally
choose `s` still smaller so both `b(1-s)^(-3/2)<1/2` and
`b+2ts<1/2`. Equations (10), (14), and `d<1/2` give (1) with a fixed
positive margin, after reducing the margin to absorb order errors.

The choices are in this order: fixed `U<b<d<1/2`; the sparse-selector
constants; then sufficiently small FIXED `s`; then the order limit.

## 6. What this has, and has not, accomplished

The theorem removes the above-half loophole in a broad pointwise
powered-aggregation proposal. Its mechanism is not a new good-signing
construction: it perturbs an existing strict-subhalf parent and uses
the established sparse-selector obstruction. The new imported upper
sampling estimate controls the perturbation cost in the original norm.

It leaves three genuinely narrower possibilities open: aggregation of
ACTUAL MINIMUM VALUES; comparable-split-only pointwise aggregation;
and an exact-minimizer structural theorem preventing the engineered
small block. No convergence or nonconvergence conclusion follows.

No finite numerical experiment is required by this existence proof.
Its immediate falsifiers are a failure of the exact-size bound (3),
the stopped expectation-to-conditional-parent extraction (11)--(12),
or the cap-repair inequality (13); all are stated explicitly above.

## 7. Stronger theorem: the small child can be an exact minimizer

Let `ell=liminf M_n/n^(3/2)`, and retain the verified all-order upper
constant `U<1/2`. For every fixed `b` with `ell<b<1/2`, and every
sufficiently small fixed RATIONAL `s>0`, there are arbitrarily large
orders `n` with `m=(1-s)n` and `h=sn` integers, and signings

```math
P=\begin{pmatrix}B&F\\F^T&D\end{pmatrix},\qquad
|B|=m,\quad |D|=h,\quad Q(D)=M_h,
```

such that all three normalized caps are at most `1/2-epsilon_s` and

```math
Q(B)^{2/3}+M_h^{2/3}-Q(P)^{2/3}\ge\kappa_s n             (17)
```

for constants `epsilon_s,kappa_s>0`. The rational restriction only
avoids floor errors in the discounted identity. Arbitrarily small
fixed splits suffice for this theorem.

### The stopped discounted identity

Put `q=1-s`, `lambda=ell^(2/3)`, and choose a fixed margin `a_1>0`
such that `b+a_1<1/2`. Let `rho` be the uniform sparse-retention
threshold for starting cap bound `1/2` and target `1/2`, as above.
Put `L=log(1/rho)` and

```math
r=\lceil L/(-\log q)\rceil+1.
```

Then `rho q^2<=q^r<rho`. If `q=u/v` in lowest terms, start with a
liminf-realizing minimizing sequence of orders `D_j`, and remove at
most `v^r-1` vertices to make its new order `N_j` divisible by `v^r`.
For fixed `s` this is a bounded deletion. Principal monotonicity and
the definition of `ell` show that these starting normalized caps
still tend to `ell`. They are at most `1/2` for large `j`.

Take `r` nested exact `q`-restrictions. Write their orders as `n_i`,
their normalized caps as `C_i`, and `u_i=C_i^(2/3)`. Let `tau` be the
first `i` with `C_i>=b`, or `r` if none exists. The terminal restriction
is marginally uniform at retention `q^r<rho`, so a hit occurs with
probability tending to one. The pointwise telescoping identity is

```math
q^\tau(u_\tau-\lambda)-(u_0-\lambda)
=\sum_{i=0}^{r-1}q^i1_{\{\tau>i\}}
                    (q u_{i+1}+s\lambda-u_i).             (18)
```

All finitely many orders tend to infinity, so `u_i>=lambda-o_j(1)`
uniformly in `i` and in the selected signing. On a hit,
`u_tau>=b^(2/3)`. Thus

```math
\liminf_j\mathbb E\left[q^\tau(u_\tau-\lambda)
                                  -(u_0-\lambda)\right]
\ge q^r(b^{2/3}-\lambda)
\ge G:=\tfrac\rho4(b^{2/3}-\lambda)>0                    (19)
```

for `q>=1/2`. No monotonicity of normalized caps is used.

### Repair each prospective small child by its actual minimum

At a pre-stop parent `A_i` of order `n_i`, the next child `B` has
order `q n_i` and complement `T` has order `h=s n_i`. Fix an arbitrary
EXACT minimizing signing `D_h` and replace only the `T` block by it.
Let `P` denote this prospective repaired parent, and define

```math
R_i=\frac{Q((A_i)_T)+M_h}{n_i^{3/2}},\qquad
\mu_h=\frac{M_h^{2/3}}h,\qquad
F_i=\frac{Q(B)^{2/3}+M_h^{2/3}-Q(P)^{2/3}}{n_i}.
```

The exact-size sampling bound (3) and the all-order upper bound give
a constant `C_R`, independent of sufficiently small fixed `s`, with

```math
\mathbb E[R_i\mid A_i]\le C_R s^{3/2}+o_j(1)             (20)
```

uniformly over all pre-stop parents. For each fixed `s`, the order
error may be absorbed by increasing `C_R` for sufficiently large `j`.
Also `mu_h>=lambda-o_j(1)`. As in (13)--(15), using the same fixed
derivative bound `L_0` on the eventual universal interval,

```math
Q(P)^{2/3}/n_i\le u_i+L_0 R_i.
```

Consequently, with `J_i=q u_(i+1)+s lambda-u_i`,

```math
J_i\le F_i+L_0R_i+s\,o_j(1).                            (21)
```

Call a prospective repair GOOD if `R_i<=a_1`. On a good repair,

```math
c(P)\le b+a_1<1/2,\qquad
c(B)\le b q^{-3/2}<1/2,\qquad
c(D_h)\le U+o_j(1)<1/2,                                 (22)
```

where the middle inequality is arranged by taking `s` small. Thus
every good repair lies in exactly the strict-subhalf scope of (17).

The useful fact on BAD repairs is the following exact inequality,
which avoids any uncontrolled high-cap tail:

```math
q u_{i+1}=Q(B)^{2/3}/n_i\le Q(A_i)^{2/3}/n_i=u_i,
\qquad J_i\le s\lambda.                                 (23)
```

Markov and (20) give conditional bad-repair probability at most
`C_R s^(3/2)/a_1`. Hence their contribution to an upper bound for the
conditional mean of `J_i` is at most
`(C_R lambda/a_1)s^(5/2)`.

### Extracting a linear defect

For a fixed starting order let `K_j` be the largest positive part
`max(F_i,0)` over all GOOD repairs at all pre-stop states of its finite
restriction process, taking `K_j=0` if there is no such repair. This
is a maximum over a finite set of actual signings. Equations
(20)--(23), conditional on each pre-stop history, give

```math
\mathbb E[J_i\mid\text{history}]
\le K_j+L_0 C_R s^{3/2}
              +(C_R\lambda/a_1)s^{5/2}+s\,o_j(1).         (24)
```

Insert (24) into (18), and use `sum_(i<r) q^i<=1/s`. By (19),

```math
G\le\liminf_j K_j/s
          +L_0 C_R\sqrt{s}+(C_R\lambda/a_1)s^{3/2}.       (25)
```

Choose fixed rational `s` small enough that the last two terms sum
to at most `G/2`. It follows that `liminf K_j>=Gs/2`. Therefore,
at every sufficiently large starting order, one good repair has
`F_i>=Gs/4`. Taking `kappa_s=Gs/4` proves (17), and (22) supplies a
uniform positive sub-half margin. Its order is at least `q^r N_j`
and tends to infinity.

The exact-minimum status applies ONLY to `D_h`. Neither `B` nor `P`
is asserted to minimize its own order. Thus the remaining obstruction
to turning this into a falsifier of the ORIGINAL minimum-value
recurrence is explicit: `Q(B)` cannot be replaced by `M_m` in (17).
The theorem is still confined to sufficiently asymmetric fixed splits.
