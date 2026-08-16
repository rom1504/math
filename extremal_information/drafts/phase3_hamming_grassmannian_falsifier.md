# Hamming Grassmannian packing: a coding obstruction and a scalable escape

**Status.** The two theorems below are proved.  The finite certificates and
small exact calculations are checked by
[`verify_phase3_hamming_grassmannian_falsifier.py`](../experiments/verify_phase3_hamming_grassmannian_falsifier.py).

This note adversarially tests the question whether packings of linear
carriers in binary Hamming space are governed by subspaces of one common
separated host.  There are two conclusions.

1. Already for one-dimensional carriers, the exact packing problem is the
   classical *nonlinear* binary coding problem, up to one codeword.  Thus a
   general closed asymptotic formula would contain the classical binary-code
   rate problem.
2. Common separated hosts are not even qualitatively complete for the **pure
   Grassmannian Hausdorff packing problem**.  An explicit seven-letter
   alphabet of two-planes, followed by an ordinary outer code, gives
   exponentially many macroscopic carrier positions in a parameter regime
   where no common separated host can contain one carrier.

The second result is only a linear-bit carrier packing, not the sought
quadratic-bit asymptotic.  Moreover, its Hausdorff gap is smaller than the
`2k` presentation toll of the mixed-profile theorem, so it does **not** give
separated multichannel response profiles.  It is nevertheless a scalable
falsifier of common-host control of the pure Grassmannian packing and
identifies a different carrier mechanism: composition of a finite alphabet
with two-sided directed separation.

## 1. The metric and the line slice

For binary linear subspaces `C,D <= F_2^N`, put

```math
\delta^\to(C,D)=\max_{c\in C}d(c,D),
\qquad
d_{\rm H}(C,D)=\max\{\delta^\to(C,D),\delta^\to(D,C)\}.       \tag{HG.1}
```

Let `A_2(N,d)` denote the largest cardinality of a not necessarily linear
binary code of length `N` and minimum Hamming distance at least `d`.  Let
`K_{N,1}(a)` be the largest number of binary lines with pairwise Hausdorff
distance strictly greater than `a`.

### Theorem HG.1 (line carriers are nonlinear codes)

For every integer `0<=a<N`,

```math
\boxed{
A_2(N,a+1)-1\le K_{N,1}(a)\le A_2(N,a+1).
}                                                               \tag{HG.2}
```

Consequently the normalized logarithmic asymptotics of line-carrier packing
and ordinary binary-code packing agree wherever either is exponentially
large.

#### Proof

Over `F_2`, a line has the unique form `L_v={0,v}` with `v!=0`.  For distinct
nonzero `v,w`, write

```math
p=\operatorname{wt}(v),\quad q=\operatorname{wt}(w),\quad
r=\operatorname{wt}(v+w).
```

Direct calculation gives

```math
d_{\rm H}(L_v,L_w)
=\max\{\min(p,r),\min(q,r)\}
=\min\{r,\max(p,q)\}.                                      \tag{HG.3}
```

Take an ordinary code of size `A_2(N,a+1)`, translate one word to zero, and
discard zero.  All remaining vectors have weight greater than `a`, and all
pairwise differences have weight greater than `a`; (HG.3) gives a line
packing of size `A_2(N,a+1)-1`.

Conversely, represent a line packing by a set `V` of nonzero vectors.  Formula
(HG.3) implies that all pairwise differences in `V` have weight greater than
`a`, and that at most one member of `V` has weight at most `a`.  If no member
is light, `{0} union V` is a distance-`a+1` code and `|V|+1<=A_2(N,a+1)`.
If `v_0` is the unique light member, translate `V` by `v_0`; the resulting
set has size `|V|`, contains zero, and all its nonzero members and pairwise
differences have weight greater than `a`.  Hence `|V|<=A_2(N,a+1)`.
This proves (HG.2). `square`

The common-host construction in this slice uses the nonzero words of one
*linear* code.  The full Grassmannian packing uses the optimal nonlinear code
rate.  Therefore the gap between common-host and unrestricted carrier
packing already contains the linear-versus-nonlinear coding gap; it cannot in
general be read from separated linear rank alone.

## 2. Directed distances compose exactly

### Lemma HG.2 (direct-sum algebra)

For subspaces `C_i,D_i` in disjoint binary coordinate blocks,

```math
\delta^\to\!\left(\bigoplus_i C_i,\bigoplus_iD_i\right)
=\sum_i\delta^\to(C_i,D_i).                                  \tag{HG.4}
```

#### Proof

For `x=(x_i)`, additivity of Hamming distance gives

```math
d\!\left(x,\bigoplus_iD_i\right)=\sum_i d(x_i,D_i).
```

Maximization over the Cartesian product `oplus_i C_i` separates into the sum
of the coordinatewise maxima. `square`

Thus a finite carrier alphabet with separation in *both* directed metrics
can be concatenated by an ordinary outer code.  This is stronger than merely
knowing a symmetric Hausdorff gap: opposing orientations cannot cancel after
composition.

## 3. A seven-letter alphabet of two-planes

In `F_2^6`, let `L_0,...,L_6` be the spans of the following pairs of rows:

```text
000111  011011
001101  100011
001111  110101
010011  100110
010110  101110
011011  101001
011100  101001
```

### Lemma HG.3 (finite carrier alphabet)

For every `i!=j`,

```math
\boxed{
\delta^\to(L_i,L_j)=\delta^\to(L_j,L_i)=3.
}                                                               \tag{HG.5}
```

#### Proof

Each span has four vectors.  Comparing its three nonzero vectors with the
four vectors of every other span gives the matrix

```math
\bigl(\delta^\to(L_i,L_j)\bigr)_{0\le i,j<7}
=3(\mathbf 1\mathbf 1^\mathsf T-I_7).                         \tag{HG.6}
```

This is a finite certificate of 42 directed comparisons; the accompanying
verifier performs all of them directly. `square`

The seven displayed planes are distinct, but no uniqueness or maximality is
claimed for this alphabet.  They span five dimensions altogether, and they
are not all subspaces of one distance-three host.  More importantly, their
two-sided metric identity persists under composition.

## 4. Scalable separation beyond every common host

For a word `u=(u_1,...,u_m) in [7]^m`, define

```math
C_u=\bigoplus_{t=1}^m L_{u_t}
\le F_2^{6m}.                                                   \tag{HG.7}
```

This is a `2m`-dimensional subspace.  Lemmas HG.2--HG.3 give the exact
isometry

```math
\delta^\to(C_u,C_v)=d_{\rm H}(C_u,C_v)
=3d_{[7]^m}(u,v).                                             \tag{HG.8}
```

Here the middle `d_H` is subspace Hausdorff distance and the last distance is
ordinary seven-letter Hamming distance.

### Theorem HG.4 (common-host packings are scalably incomplete)

There are integers `Delta_m` with

```math
{\Delta_m\over6m}\longrightarrow{3\over8}                    \tag{HG.9}
```

and families

```math
\mathcal F_m\subseteq \operatorname{Gr}_{2m}(F_2^{6m})
```

such that

```math
d_{\rm H}(C,C')>\Delta_m\quad(C\ne C'),                      \tag{HG.10}
```

and

```math
\log_2|\mathcal F_m|
\ge
\left(\log_2 7\right)
\left(1-H_7(3/4)-o(1)\right)m
= (0.0573549\ldots-o(1))m.                                  \tag{HG.11}
```

Nevertheless, no member of `F_m` lies in a common host of minimum nonzero
weight greater than `Delta_m`: every member already contains a nonzero word
of weight at most four in one block.  Independently, for all sufficiently
large `m`, the Hamming bound shows that no binary linear subspace of minimum
weight greater than `Delta_m` can have dimension `2m` at all.  Thus the
common-separated-host construction cannot produce a carrier in this regime,
while the unrestricted Grassmannian contains exponentially many mutually
separated carriers.

#### Proof

The `q`-ary Gilbert greedy bound gives seven-letter codes
`Q_m subseteq[7]^m` of minimum distance

```math
d_m=\lceil3m/4\rceil
```

and cardinality

```math
|Q_m|\ge7^{(1-H_7(3/4)-o(1))m};                              \tag{HG.12}
```

the exponent is positive because `3/4<1-1/7`.  Put
`F_m={C_u:u in Q_m}` and `Delta_m=3d_m-1`.  Equation (HG.8) proves
(HG.9)--(HG.11).

Each `C_u` contains a nonzero vector supported on one six-coordinate block,
of weight at most four.  Since `Delta_m` grows linearly, this already proves
that the displayed carriers do not lie in any host of minimum distance
greater than `Delta_m`.

For the stronger ambient statement, suppose a common host
`E_m<=F_2^{6m}` had dimension at least `2m` and minimum distance greater than
`Delta_m`.  The binary Hamming bound, applied with

```math
t_m=\left\lfloor{3d_m-1\over2}\right\rfloor,
```

would give

```math
2^{2m}\sum_{j=0}^{t_m}{6m\choose j}\le2^{6m}.                \tag{HG.13}
```

After taking logarithms and dividing by `6m`, its left exponent tends to

```math
{1\over3}+H_2(3/16)
=1.02954559\ldots>1,                                         \tag{HG.14}
```

a contradiction. `square`

## 5. When finite carrier alphabets survive presentation

The seven-plane construction has presentation radius larger than its useful
Hausdorff gap.  A different finite-alphabet theorem states exactly when this
problem disappears.

Let `(X,d_X)` be a finite metric space and give `X^m` the `ell_1` product
metric.  For `i in[q]`, let `C_i subseteq X` have a presentation
`pi_i:C_i->[0,p]`.  For a word `u in[q]^m`, use the product carrier and
additive presentation

```math
C_u=\prod_{t=1}^m C_{u_t},
\qquad
\pi_u(c_1,\ldots,c_m)=\sum_{t=1}^m\pi_{u_t}(c_t),           \tag{HG.15}
```

and let

```math
F_u(x)=\min_{c\in C_u}\{d(x,c)+\pi_u(c)\}.                  \tag{HG.16}
```

### Theorem HG.5 (joint same-sign finite-alphabet amplification)

Write `f_i` for the one-block response and define its directed response
weight

```math
\omega(i,j)=\max_{x\in X}\{f_i(x)-f_j(x)\}.                 \tag{HG.16a}
```

Then the product response metric has the exact algebra

```math
\boxed{
\|F_u-F_v\|_\infty
=\max\left\{
 \sum_t\omega(u_t,v_t),
 \sum_t\omega(v_t,u_t)
 \right\}.
}                                                            \tag{HG.16b}
```

Assume that for every `i!=j`, both directed carrier distances satisfy

```math
\delta^\to(C_i,C_j)\ge d,
\qquad
\delta^\to(C_j,C_i)\ge d.                                  \tag{HG.17}
```

If `Q subseteq[q]^m` has symbol distance at least `rho*m`, then

```math
\boxed{
\|F_u-F_v\|_\infty
\ge(d-p)\rho m
}\qquad(u\ne v\in Q).                                      \tag{HG.18}
```

In particular, a positive-rate `q`-ary Gilbert family gives a positive
linear response gap whenever

```math
\boxed{
d>p
\quad\text{and}\quad
0<\rho<1-{1\over q}.
}                                                            \tag{HG.19}
```

#### Proof

Exact factorization of (HG.16) gives

```math
F_u(x_1,\ldots,x_m)=\sum_{t=1}^m f_{u_t}(x_t).
```

Maximizing its signed difference from `F_v` separates coordinatewise, and
doing the same after reversing the sign proves (HG.16b).

If `delta^to(C_i,C_j)>=d`, choose `c in C_i` at distance at least `d` from
`C_j`.  At the query `c`,

```math
f_j(c)-f_i(c)\ge d-p.
```

The reverse directed carrier inequality gives the opposite signed response
margin.  Thus `omega(i,j),omega(j,i)>=d-p`.  In the exact algebra (HG.16b),
matching symbols contribute zero, while every differing symbol contributes
at least `d-p` to the same chosen sign.  Hence a pair at symbol distance `h`
satisfies

```math
\|F_u-F_v\|_\infty\ge(d-p)h.
```

Insert `h>=rho*m`.  The `q`-ary Gilbert exponent is positive for every
`rho<1-1/q`, proving the final assertion. `square`

This is not a reformulation of ordinary concatenated coding: the local
condition is two-sided directed Hausdorff exposure, and the theorem charges
the presentation radius on each *changed* coordinate.  A symmetric local
carrier gap without orientation can lose a factor under products, while a
local presentation toll `p>=d` can erase the certified margin.  The cruder
one-shot carrier estimate `d rho m-pm` pays even matching coordinates and is
strictly weaker by `p(1-rho)m`.

### Corollary HG.6 (binary simplex response family)

Let `S<=F_2^7` be the binary `[7,3,4]` simplex code, and use its seven nonzero
words `v_1,...,v_7`.  The line carriers

```math
L_i=\{0,v_i\}
```

obey

```math
\delta^\to(L_i,L_j)=4\qquad(i\ne j),                         \tag{HG.20}
```

because every nonzero simplex word, including `v_i+v_j`, has weight four.
Give each line its standard presentation `pi_i(0)=0`, `pi_i(v_i)=2`.

There are families of profiles on `F_2^{7m}`,

```math
F_u(x)=\min_{z\in F_2^m}
\left\{2\operatorname{wt}(z)
+\operatorname{wt}\!\left(x+
 \bigoplus_{t=1}^m z_tv_{u_t}\right)\right\},               \tag{HG.21}
```

of cardinality

```math
7^{(1-H_7(3/4)-o(1))m}                                     \tag{HG.22}
```

with pairwise uniform response distance at least

```math
{3m\over2}.                                                  \tag{HG.23}
```

#### Proof

Apply HG.5 with `q=7`, `d=4`, `p=2`, and a seven-letter outer code of relative
distance `3/4`.  Here `(d-p)rho=3/2`. `square`

This is a genuine presented-response theorem: unlike the two-plane carrier
packing, its linear separation survives the exact standard presentation
cost.

### Corollary HG.7 (rank-metric multiplication response family)

Identify `E=F_8` as a three-dimensional binary space.  For every
`a in E^*`, let `M_a` be the binary matrix of multiplication by `a` and take
the rank-metric line `R_a={0,M_a}`.  For `a!=b`, both `M_a` and
`M_a-M_b=M_(a-b)` are invertible, so

```math
\delta^\to(R_a,R_b)=3.                                      \tag{HG.24}
```

On `m` block-diagonal copies of `End_(F_2)(E)`, with additive block-rank
metric and the standard line presentation, there are

```math
7^{(1-H_7(3/4)-o(1))m}                                     \tag{HG.25}
```

profiles separated in uniform norm by at least

```math
{3m\over4}.                                                  \tag{HG.26}
```

#### Proof

Apply HG.5 with `q=7`, `d=3`, `p=2`, and `rho=3/4`.  Then
`(d-p)rho=3/4`. `square`

Thus the same response-amplification law applies in Hamming and rank
geometry.  Its actual invariant is the triple `(q,d,p)`: alphabet entropy,
two-sided directed exposure, and presentation cost.

## 6. Adversarial judgment

The conjecture that the common separated-host exponent always controls the
pure Hamming Grassmannian carrier entropy is false.  The failure is not a
small finite anomaly: a two-sided local carrier alphabet plus an outer code
creates macroscopic Hausdorff separation under arbitrarily deep direct-sum
composition, even when sphere packing rules out a host of the carrier
dimension.

There is an important non-implication.  In Theorem HG.4,

```math
k=2m,
\qquad
d_{\rm H}(C_u,C_v)\ge(9/4+o(1))m,
```

whereas the presentation radius in the mixed profile
`F_V(x)=min_z(2wt(z)+||x+Vz||)` is `2k=4m`.  The carrier-to-response lower
bound `d_H-2k` is therefore vacuous.  HG.4 cannot be cited as a multichannel
response-rate lower bound without a new presentation theorem.  What it
falsifies is the proposed common-host characterization of
`Pack(Gr_k,d_H,Delta)` itself.

This does **not** determine the quadratic-scale quantity

```math
\log_2\operatorname{Pack}
(\operatorname{Gr}_{\kappa N}(F_2^N),d_{\rm H},\delta N).
```

The construction gives `Theta(N)` carrier bits, whereas the anticode
quotient upper bound can be `Theta(N^2)`.  It therefore isolates a real middle
mechanism but does not show whether it fills a positive fraction of the
code--anticode gap on the quadratic scale, nor whether any of this carrier
entropy survives a linear presentation toll.

The clean next question is no longer “common host or quotient upper bound?”
It is:

> Can one build carrier alphabets whose logarithmic size grows quadratically
> in their block length while retaining two-sided directed Hausdorff
> separation, or prove a container theorem forcing every such alphabet into
> only exponentially many compositional types?

The line theorem warns that even the smallest-dimensional slice imports
classical coding-rate difficulty.  A plausible general invariant must retain
both nonlinear carrier packing and the orientation of directed exposure; the
separated rank alone forgets both.
