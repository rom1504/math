# Finite-alphabet amplification of extremal responses

**Status.** The theorem and both applications below are proved.  The outer
amplification step is classical coding theory.  The new response-theoretic
content is a joint-channel inequality: two signed local margins let all
channel differences be aligned before the absolute value, so matching
channels pay no presentation toll.

## 1. The amplification law

Let `(X,d)` be a finite metric space for which Cartesian products carry the
`ell_1` metric.  Let

```math
\mathcal A=\{(C_a,\pi_a):a\in[q]\}
```

be a finite alphabet of presented carriers, where `C_a subseteq X` is
nonempty and `0<=pi_a<=p`.  Assume the two directed gaps obey

```math
h^\to(C_a,C_b)\ge d_0,qquad
h^\to(C_b,C_a)\ge d_0
\quad(a\ne b).                                  \tag{FA.1}
```

For `a=(a_1,...,a_m)`, take the product carrier and additive presentation

```math
C_{\boldsymbol a}=\prod_{i=1}^m C_{a_i},qquad
\pi_{\boldsymbol a}(c)=\sum_{i=1}^m\pi_{a_i}(c_i),            \tag{FA.2}
```

with response

```math
F_{\boldsymbol a}(x)=
\min_{c\in C_{\boldsymbol a}}
\{d_1(x,c)+\pi_{\boldsymbol a}(c)\}.                          \tag{FA.3}
```

For the one-block responses, define the directed response table

```math
r(a,b)=\sup_{x\in X}\{F_a(x)-F_b(x)\}.                       \tag{FA.3a}
```

### Theorem FA.1 (directed-margin amplification)

For every two words one has the exact algebra

```math
\boxed{
\|F_{\boldsymbol a}-F_{\boldsymbol b}\|_\infty
=\max\left\{
 \sum_i r(a_i,b_i),
 \sum_i r(b_i,a_i)
 \right\}.}                                                   \tag{FA.4a}
```

If `Q subseteq[q]^m` has minimum symbol distance at least `rho m`, then the
responses indexed by `Q` obey

```math
\|F_{\boldsymbol a}-F_{\boldsymbol b}\|_\infty
\ge(d_0-p)\rho m
\qquad(\boldsymbol a\ne\boldsymbol b).                        \tag{FA.4}
```

For every `rho<1-1/q`, there are such outer codes with

```math
\log_2|Q|
\ge(1-H_q(\rho)-o(1))m\log_2q.                               \tag{FA.5}
```

Consequently a positive linear response rate is certified whenever

```math
\boxed{d_0>p.}                                               \tag{FA.6}
```

More precisely, choose any fixed `0<rho<1-1/q`.  Uniform distortion below
`(d_0-p)rho m/2` then requires at least the number of bits in (FA.5).

#### Proof

For one ordered pair, `0<=pi_a<=p` gives

```math
F_a(x)\ge d(x,C_a),\qquad F_b(x)\le d(x,C_b)+p.
```

Distance functions satisfy

```math
\sup_x\{d(x,C_a)-d(x,C_b)\}=h^\to(C_b,C_a):                 \tag{FA.7}
```

the upper bound is the triangle inequality, and equality follows by taking
`x` at a farthest point of `C_b`.  Hence (FA.1) gives both signed margins

```math
\sup_x(F_a-F_b)\ge d_0-p,qquad
\sup_x(F_b-F_a)\ge d_0-p.                                   \tag{FA.8}
```

The product response factors exactly:

```math
F_{\boldsymbol a}(x_1,\ldots,x_m)=\sum_iF_{a_i}(x_i).        \tag{FA.9}
```

Maxima over independent query coordinates therefore separate:

```math
\sup_x(F_{\boldsymbol a}-F_{\boldsymbol b})
=\sum_i\sup_{x_i}(F_{a_i}-F_{b_i}).                         \tag{FA.10}
```

Applying (FA.10) in both orientations proves (FA.4a).  Equal symbols
contribute zero and every differing symbol contributes at least `d_0-p` to
the same chosen sign.  This proves (FA.4).  The `q`-ary
Gilbert greedy bound proves (FA.5), and a separated response family requires
one decoder state per radius-half ball. `square`

The requirement of two directed gaps is essential.  A symmetric Hausdorff
gap can choose opposite witnessing directions in different blocks; then
the two directed sums may each lose a leading fraction.  Condition (FA.1)
is the finite deterministic synchronization needed for composition.

This theorem is not a new proof of the Gilbert bound.  It identifies the
signed local margin `d_0-p` that coding theory can amplify.  A coarse global
carrier comparison gives only `d_0 rho m-pm`; (FA.10) gives
`(d_0-p)rho m`, larger by `p(1-rho)m`, because matching channels pay
nothing.  This is a genuinely joint same-sign estimate rather than a sum of
independently paid scalar bounds.

For the declared query—uniform distances after arbitrary direct-product
composition—the `q x q` table `r` is an exact feature algebra.  It can be
strictly smaller than the collection of full functions `F_a:X->R`, and its
update is ordinary addition followed only at the end by the two-orientation
maximum.  It does not reconstruct the responses themselves or control an
arbitrary nonproduct coupling, so the quotient claim is deliberately
query-relative.

## 2. Binary Hamming simplex carriers

Index the seven nonzero points `x` of `F_2^3`.  For every nonzero
`a in F_2^3`, define

```math
v_a=(a\mathbin\cdot x)_{x\ne0}\in F_2^7.
```

Every `v_a` has Hamming weight four, and
`v_a+v_b=v_(a+b)` also has weight four for `a!=b`.  Thus the seven binary
lines `C_a={0,v_a}` have both directed Hausdorff distances equal to four.
Give `v_a` access cost two and zero access cost zero, so `p=2`.

Take `rho=3/4`.  Theorem FA.1 gives at least

```math
7^{(1-H_7(3/4)-o(1))m}
=2^{(0.0573549\ldots-o(1))m}                                 \tag{FA.12}
```

responses on `F_2^(7m)` separated by at least

```math
(4-2)(3/4)m={3m\over2}.                                      \tag{FA.13}
```

Concretely, for an outer word `boldsymbol a`, let `V_a:F_2^m->F_2^(7m)`
send the `i`-th basis vector to `v_(a_i)` in block `i`.  The response is the
mixed-channel profile

```math
F_{\boldsymbol a}(u)=
\min_{z\in F_2^m}
\{2\operatorname{wt}(z)+\operatorname{wt}(u+V_{\boldsymbol a}z)\}. \tag{FA.14}
```

Every carrier contains a weight-four vector supported in one block, so no
common linear host of minimum distance growing with `m` contains even one of
them.  The information comes from oriented separation between different
carriers, not internal separation inside one host.  This is a scalable
response mechanism outside the separated-host lower certificate, though its
information rate is linear rather than the unresolved quadratic rate.

## 3. Rank-metric multiplication carriers

Let `E=F_8` and represent multiplication by `a in E` as a `3 x 3` binary
matrix `M_a`.  For every nonzero `a`, `M_a` is invertible; for `a!=b`,
`M_a+M_b=M_(a+b)` is also invertible.  The seven binary lines

```math
C_a=\{0,M_a\}\subseteq M_3(F_2)
```

therefore have both directed rank-Hausdorff distances equal to three.
Again give the nonzero point access cost two.  On `m` block-diagonal matrix
blocks, rank is additive, so Theorem FA.1 with `rho=3/4` gives the same
number (FA.12) of responses separated by at least

```math
(3-2)(3/4)m={3m\over4}.                                      \tag{FA.15}
```

This second application is intrinsic to rank distance: it uses invertible
field-multiplication differences, not an embedded Hamming norm.

## 4. What the theorem does and does not establish

FA.1 is a general law for when repeated composition converts microscopic
carrier distinctions into macroscopic extremal-response information.  It
unifies three ingredients that were previously separate:

1. oriented local exposure (`d_0` in both directions);
2. access/presentation loss (`p`); and
3. outer query mass (`rho`).

It applies to Hamming and rank-metric responses and is stable under arbitrary
composition depth.  It does not compress the resulting state, determine the
quadratic Hamming Grassmannian exponent, or overcome the coding-rate problem
already present at one channel.  Its principal falsification criterion is
sharp at the local signed-margin level: if `d_0<=p`, carrier separation
alone certifies no positive local response margin, regardless of outer
distance.
