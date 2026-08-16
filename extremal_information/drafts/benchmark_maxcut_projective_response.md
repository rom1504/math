# Benchmark: pure Max-Cut through a width-`w` separator

**Status.** Working theorem draft.  Every claim below is proved, but the
Fourier calculation and the information interpretation still require an
independent audit before promotion to `theorems.md`.

## 1. Operational state, derived from future responses

Let `G` be a finite nonnegatively weighted graph with a distinguished ordered
boundary `B=[w]`.  Edges belong to one side of a gluing exactly once.  For a
boundary spin vector `sigma in {+1,-1}^w`, define

```math
h_G([\sigma])=max_{z\in\{+1,-1\}^{V(G)\setminus B}}
 \sum_{uv\in E(G)}c_{uv}{\bf1}\{s_u\ne s_v\}.       \tag{MC.1}
```

The argument is the projective class
`[sigma]={sigma,-sigma}`, because a global spin flip changes no cut.  Write

```math
X_w=\{+1,-1\}^w/\{\sigma\sim-\sigma\},
\qquad q=|X_w|=2^{w-1}.                              \tag{MC.2}
```

If a continuation `C` meets `G` only in `B`, conditional independence gives

```math
\operatorname {MaxCut}(G\cup_B C)
=\max_{s\in X_w}\{h_G(s)+h_C(s)\}.                  \tag{MC.3}
```

This derives the boundary profile as a sufficient state without assuming a
treewidth dynamic program.

### Proposition MC.1 (pure-Max-Cut response isometry)

Let the future-query family consist of all nonnegatively weighted Max-Cut
attachments through `B`.  Then

```math
\sup_C\left|
 \operatorname {MaxCut}(G\cup_B C)
-\operatorname {MaxCut}(G'\cup_B C)
\right|
=\|h_G-h_{G'}\|_\infty.                             \tag{MC.4}
```

Consequently equality of the projective boundary profiles is the coarsest
exact contextual equivalence for this fixed constraint language.  The same
identity says that its lossy response geometry is exactly the sup geometry
of the *realizable* profile class.

#### Proof

The upper bound follows from (MC.3).  For the reverse bound, fix a target
`tau in {+1,-1}^w`.  Introduce an internal anchor `a`.  If `tau_i=-1`, add a
unit edge `a i`; if `tau_i=+1`, add a two-edge unit path `a-v_i-i` through a
fresh internal vertex.  For fixed boundary `sigma` and anchor spin, an
inequality edge loses one precisely when its prescribed relation fails, and
the maximized equality path scores two when its prescribed relation holds
and one otherwise.  Maximizing the anchor therefore gives the continuation
profile

```math
p_\tau([\sigma])=c_\tau-d_w([\sigma],[\tau]),
\qquad
d_w([\sigma],[\tau])=min\{d_H(\sigma,\tau),w-d_H(\sigma,\tau)\}. \tag{MC.5}
```

Taking sufficiently many copies makes `[tau]` the unique maximizing boundary
class for both `G` and `G'`.  Their two continued optima then differ by
`h_G([tau])-h_G'([tau])`.  Choose a class attaining the sup norm and, if
necessary, exchange `G,G'`. `square`

The proof separates two issues that are easy to conflate.  Future Max-Cut
attachments can expose every profile coordinate.  That alone does not show
that the fixed Max-Cut language realizes an exponentially rich family of
profiles.  The next argument supplies that missing packing.

## 2. A general generated-profile packing lemma

Let `S` have cardinality `q`.  Suppose a compositional landscape class has
additive gadget profiles `p_t in R^S`, `t in S`: a disjoint union of gadgets
indexed by `U subset S` has profile

```math
h_U=\sum_{t\in U}p_t.                               \tag{MC.6}
```

Let `P` be the `q by q` matrix whose `t`-th column is `p_t`, and let
`s_min(P)` be its least singular value.

### Proposition MC.2 (spectral atom-packing principle)

For every fixed `0<delta<1/2`, there is a family `F` of subsets of `S`, all of
cardinality `floor(q/2)`, such that

```math
|F|\ge {2^{(1-H_2(\delta))q}\over q+1}              \tag{MC.7}
```

up to an inessential rounding change in `delta`, and every distinct `U,V` in
`F` satisfy

```math
\|h_U-h_V\|_\infty
\ge s_{\min}(P)\sqrt\delta.                         \tag{MC.8}
```

If the future-response metric is profile sup distance, any summary uniformly
answering all future queries to error less than half the right-hand side of
(MC.8) requires at least

```math
(1-H_2(\delta))q-O(\log q)                          \tag{MC.9}
```

bits on this generated subclass.

A literal finite version puts `d=ceil(delta q)` and replaces the right-hand
side of (MC.7) by

```math
{\binom q{\lfloor q/2\rfloor}\over
  \sum_{i=0}^{d-1}\binom qi};
```

the entropy form follows with `delta+O(1/q)`.

#### Proof

A greedy code in the middle layer of the Boolean `q`-cube gives (MC.7): the
middle layer has at least `2^q/(q+1)` words, while a Hamming ball of relative
radius `delta` contains at most `2^{H_2(delta)q}` words.  For two selected
incidence vectors `u,v`,

```math
\|P(u-v)\|_\infty
\ge q^{-1/2}\|P(u-v)\|_2
\ge q^{-1/2}s_{\min}(P)\|u-v\|_2
\ge s_{\min}(P)\sqrt\delta.                         \tag{MC.10}
```

The last assertion is the packing bound for deterministic lossy summaries.
`square`

This principle is not merely the metric entropy of an ambient response cube.
It turns a model-specific compositional gadget family into a lower bound via
the conditioning of its generated profile algebra.  For universal boundary
kernels `P` may be the identity.  Below, `P` is a highly structured distance
matrix generated inside pure Max-Cut.

## 3. Projective-cube Fourier calculation

Pad the gadget in (MC.5) by boundary-independent internal edges so that every
atom has the same constant `C=2w`.  Its profile is

```math
p_\tau(s)=C-d_w(s,\tau).                             \tag{MC.11}
```

Let `D_w=(d_w(s,t))_(s,t in X_w)`.  The matrix is a convolution matrix on the
projective cube.  Its characters are

```math
\chi_A([\sigma])=\prod_{i\in A}\sigma_i,
\qquad |A|\ \hbox{even}.                            \tag{MC.12}
```

### Lemma MC.3 (all projective distance channels survive)

Assume `w>=2` and put `N=2 ceil(w/2)=2m`.  For every nonempty even level
`|A|=2j`, the absolute
eigenvalue of `D_w` on `chi_A` is

```math
L_{w,j}=
\begin{cases}
\displaystyle
{(2j-2)!(2m-2j)!\over
 (m-1)!(j-1)!(m-j)!},&w=2m,\\[6pt]
\displaystyle {1\over2}
{(2j-2)!(2m-2j)!\over
 (m-1)!(j-1)!(m-j)!},&w=2m-1.
\end{cases}                                        \tag{MC.13}
```

Here `1<=j<=m` when `w=2m`, while `1<=j<=m-1` when `w=2m-1`.

Every value in (MC.13) is nonzero.  The trivial eigenvalue is the positive
row sum, so `D_w` is nonsingular.  Moreover, if

```math
\Lambda_m=
\begin{cases}
\binom{2k}{k},&m-1=2k,\\
2\binom{2k}{k},&m-1=2k+1,
\end{cases}                                        \tag{MC.14}
```

then

```math
s_{\min}(D_w)=
\begin{cases}\Lambda_m,&w=2m,\\
\Lambda_m/2,&w=2m-1,
\end{cases}
\qquad
s_{\min}(D_w)\ge {2^{m-2}\over m}.                 \tag{MC.15}
```

#### Proof

Write `S_w=sum_i sigma_i`.  Since

```math
d_w([\boldsymbol1],[\sigma])={w-|S_w|\over2},       \tag{MC.16}
```

every nontrivial even Fourier coefficient of `d_w` is minus one half the
corresponding coefficient of `|S_w|`.  For `w=2m`, the elementary
Krawtchouk transform gives

```math
c_2={\binom{2m-2}{m-1}\over2^{2m-2}},
\qquad
{c_{2j+2}\over c_{2j}}
=-{2j-1\over2m-2j-1}.                              \tag{MC.17}
```

For completeness, with
`K_r(k)=sum_l(-1)^l binom(k,l)binom(2m-k,r-l)`, the coefficient of one fixed
level-`r` character is

```math
{2^{-2m}\over\binom{2m}{r}}
\sum_{k=0}^{2m}\binom{2m}{k}|2m-2k|K_r(k).          \tag{MC.17a}
```

Pairing `k` with `2m-k` and applying Pascal's identity twice gives (MC.17).
For a nontrivial character the distance coefficient is `-c_(2j)/2`, and a
convolution eigenvalue is `|X_w|` times that coefficient.  Iterating
(MC.17) with these two normalization factors gives (MC.13).

For `w=2m-1`, condition a `2m`-spin sum on its last spin.  The first sum is
odd, and hence

```math
{1\over2}(|s+1|+|s-1|)=|s|.
```

Thus its even Fourier coefficients agree with those for `2m` spins, while
the projective cube has half as many points.  This gives the second line of
(MC.13).

Writing `a=j-1`, `b=m-j`, the first line of (MC.13) is
`(2a)!(2b)!/((a+b)!a!b!)`.  Consecutive ratios show that this is minimized
when `a,b` are as balanced as possible, giving (MC.14).  The central-binomial
bound `binom(2k,k)>=2^(2k)/(2k+1)` proves (MC.15).  The row sum dominates the
absolute value of every character sum, so it cannot be the smaller singular
value. `square`

The matrix with columns (MC.11) differs from `-D_w` only in its trivial
Fourier channel.  Choosing `C=2w` makes that channel positive and nonzero;
on every nontrivial character its singular values are exactly those of
`D_w`.  In the constant-cardinality packing below the constant channel
cancels, so only `D_w` is needed.

## 4. A restricted-language width-scale lower bound

### Theorem MC.4 (pure-Max-Cut separator response requires exponential bits)

Fix `epsilon>0`.  For all sufficiently large `w`, there is a family of
unit-weight pure Max-Cut components with width-`w` boundary such that

1. every component is the disjoint interior union of exactly
   `floor(2^(w-1)/2)` pinning gadgets;
2. distinct conditional profiles have future-response distance greater than
   `2 epsilon w`; and
3. the family has cardinality

```math
2^{(1-H_2(\delta)-o(1))2^{w-1}}                   \tag{MC.18}
```

for any fixed `0<delta<1/2`.

Therefore every deterministic summary answering **all pure-Max-Cut future
attachments** to the width-scaled raw additive error `epsilon w` requires

```math
\Omega(2^w)                                        \tag{MC.19}
```

bits in the unrestricted-size worst case.

#### Proof

For `U subset X_w`, disjoint additivity of the padded pinning gadgets gives

```math
h_U=C|U|\boldsymbol1-D_w\boldsymbol1_U.             \tag{MC.20}
```

Choose the constant-weight family in Proposition MC.2.  The first term in
(MC.20) cancels between two members, and Lemma MC.3 gives

```math
\|h_U-h_V\|_\infty
\ge s_{\min}(D_w)\sqrt\delta
\ge {\sqrt\delta\,2^{\lceil w/2\rceil-2}
       \over\lceil w/2\rceil}.                      \tag{MC.21}
```

This exceeds `2 epsilon w` for sufficiently large `w`.  Apply the isometry
(MC.4) and the packing argument. `square`

This closes a gap left by the universal-kernel benchmark: the exponential
lower bound now lives inside ordinary nonnegative Max-Cut, not in a class
allowed to install an arbitrary `q by q` boundary factor.  The price is that
the partial graphs themselves have `Theta(w2^w)` edges and total weight.  The
error `epsilon w` is therefore exponentially small relative to the total
objective scale.  The theorem does **not** imply the same lower bound under a
polynomial-size, bounded-total-weight, or constant-relative-error promise;
that is the correct next rate-distortion question.

## 5. Pure-Max-Cut lookup universality

The distance-atom construction proves a restricted-language lower bound with
unit weights.  Allowing arbitrary nonnegative weights gives a stronger and
conceptually simpler statement: pure Max-Cut realizes every projective
boundary response table, up to a controllable common offset.

### Theorem MC.5 (every projective table is a pure-Max-Cut response shape)

Let `F:X_w -> [0,W]`.  There is a nonnegatively weighted pure Max-Cut
component `G_F` of treewidth at most `w+1` whose boundary profile is

```math
h_{G_F}(s)=F(s)+(6w-2)\sum_{t\in X_w}F(t).           \tag{MC.22}
```

By adding one boundary-independent internal edge, every member of the whole
cube `[0,W]^(X_w)` can instead be realized with the same offset

```math
C_{w,W}=(6w-2)|X_w|W.                               \tag{MC.23}
```

Consequently, for `0<epsilon<=W/6`, the logarithmic covering and packing
complexities of pure-Max-Cut boundary responses, restricted to the translated
cube `C_(w,W)+[0,W]^(X_w)`, are

```math
\Theta\left(2^{w-1}\log {W\over\epsilon}\right)    \tag{MC.24}
```

bits, up to universal changes of covering/packing radius.  In particular,
with `W=Theta(w)`, additive error `epsilon_0 w` for sufficiently small fixed
`epsilon_0` requires `Theta(2^w)` bits.

#### Proof

Choose the symmetric lift `lambda_a=F([a])` for every oriented word
`a in {+1,-1}^w`.  Introduce a common internal anchor spin `z` and, for every
oriented `a`, an internal spin `t_a`.  Put

```math
y_a={1+t_az\over2},
\qquad x_i=s_i z.
```

Consider first the signed pairwise Ising energy

```math
E(s,z,t)=\sum_a\lambda_a y_a
 \left(\sum_i a_i x_i-(w-1)\right).                 \tag{MC.25}
```

For fixed `s,z`, the coefficient of `y_a` is one when `a=(s_i z)_i` and is
`1-2d_H(a,(s_i z)_i)<=-1` otherwise.  Maximizing the independent `t_a`
therefore gives exactly `lambda_(sz)=F([s])`.

Although (MC.25) looks cubic, `z^2=1` makes it pairwise:

```math
E=-{w-1\over2}\sum_a\lambda_a
 +{1\over2}\sum_{a,i}\lambda_a a_i(s_i z+t_as_i)
 -{w-1\over2}\sum_a\lambda_a t_a z.                \tag{MC.26}
```

Every signed pair term `Juv` has a pure-Max-Cut implementation up to a
constant.  If `J<=0`, a direct edge of weight `-2J` scores

```math
-J+Juv.
```

If `J>=0`, a fresh two-edge path with both weights `2J`, maximized over its
middle spin, scores

```math
3J+Juv.                                             \tag{MC.27}
```

Apply these replacements to each displayed occurrence in (MC.26), without
combining parallel signed terms.  Pairing the oriented words `a,-a` shows
that the added constant, including the removal of the first constant in
(MC.26), is `(6w-2)lambda_a` per projective class.  This proves (MC.22).

A tree decomposition has a central bag containing the boundary and `z`, one
bag containing the boundary, `z`, and each `t_a`, and three-vertex leaf bags
for the path mediators.  Its width is at most `w+1`.

The offset in (MC.22) is at most `C_(w,W)`.  An isolated internal edge of
weight `C_(w,W)-(6w-2)sum_tF(t)` is always cut at an optimum and pads every
profile to (MC.23).

Proposition MC.1 identifies contextual response distance with entrywise sup
distance.  The usual grid packing and covering of a `q`-dimensional sup cube,
where `q=2^(w-1)`, now gives (MC.24). `square`

This theorem is stronger than saying that arbitrary binary CSP tables can be
installed at a separator: it uses only positive-weight cut edges.  Global
flip symmetry is the only exact table obstruction.  Its graphs and their
common offsets are exponential in `w`; in the displayed construction each
boundary vertex has weighted load `4 sum_tF(t)`.  It therefore still does not
settle polynomial-size, unit-boundary-sensitivity, or constant-relative-error
compression.

## 6. What the benchmark predicted

The response framework predicts the projective boundary profile, not the raw
`2^w` assignment table: global flip is a contextual gauge symmetry.  It also
predicts two logically separate minimality tests:

1. query exposure, proved by the pinning continuation; and
2. metric entropy of profiles generated by the fixed language, proved by the
   projective-distance atom algebra and, more strongly, the gauge-anchored
   lookup construction.

Standard treewidth dynamic programming supplies sufficiency and gluing, but
does not by itself give Proposition MC.1 or the approximate restricted-
language lower bound MC.4.  The spectral atom-packing principle is the
general lesson: composition creates macroscopic future-response information
when a large gadget alphabet acts through a sufficiently well-conditioned
profile matrix.
