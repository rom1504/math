# Tropical defect saturation: when approximate profiles compose forever

**Status.** Convex/tropical specialist draft for Phase 3, checkpoint 2.
The min-plus Kleene closure used below is classical.  The new project-level
statement is the sharp operator defect law (TDS.8): it separates one-shot
profile approximation, bounded-size composition, and arbitrary-depth
uniform composition.  It also gives a strict, exactly composable quotient
for a nontrivial fixed-chart syndrome class.  Nothing here concerns the
original signing problem.

## 1. Kernel coarse-graining in the min-plus algebra

Let `G` be a finite abelian group and put

```math
\overline{\mathbb R}_+=[0,\infty].
```

For profiles `f,g:G -> overline R_+`, write

```math
(f\star g)(x)=\min_{y\in G}\{f(y)+g(x-y)\}.       \tag{TDS.1}
```

Let `delta_0` be zero at `0` and `+infinity` elsewhere.  Fix a finite-valued
kernel `b:G->[0,infinity)` with `b(0)=0`, and define the coarse-graining

```math
P_bf=f\star b.                                    \tag{TDS.2}
```

Write `b^(star m)` for the `m`-fold convolution power.  Since a zero-cost
step may be appended,

```math
b^{\star(m+1)}\le b^{\star m}.                    \tag{TDS.3}
```

Its pointwise limit is denoted by

```math
b_*:=\inf_{m\ge1}b^{\star m}.                    \tag{TDS.4}
```

On a finite group the infimum is attained after finitely many powers: view
`b(x-y)` as the cost of the directed edge `y -> x`.  A cheapest walk can
have every nonnegative-cost cycle deleted, and therefore uses at most
`|G|-1` nonzero steps.  Thus `b_*` is the usual shortest-path or Kleene
closure of `b`.

### Lemma TDS.1 (the exact stable kernel)

The function `b_*` is the greatest subadditive minorant of `b` satisfying
`b_*(0)=0`.  In particular,

```math
b_*\star b_*=b_*.                                 \tag{TDS.5}
```

Moreover, the following are equivalent:

1. `b\star b=b`;
2. `b` is subadditive, meaning
   `b(x)<=b(y)+b(x-y)` for every `x,y`;
3. `b=b_*`;
4. `P_bf\star P_bg=P_b(f\star g)` for all profiles `f,g`.

#### Proof

Every decomposition of `x` is a path from zero to `x`, so (TDS.4) is the
shortest-path cost.  It is subadditive by concatenating paths.  If `c<=b` is
subadditive, then `c(x)` is at most the sum of the `b`-costs along every path
to `x`; hence `c<=b_*`.  This proves the first assertion and (TDS.5).

If `b` is subadditive, every term in `(b star b)(x)` is at least `b(x)`,
while the choice `y=0` gives the opposite inequality.  Conversely,
`b star b=b` implies that `b(x)` is at most every term
`b(y)+b(x-y)`.  This proves the equivalence of 1--3.  Finally,

```math
P_bf\star P_bg=f\star g\star b\star b,
\qquad
P_b(f\star g)=f\star g\star b.                  \tag{TDS.6}
```

Thus 1 implies 4.  Taking `f=g=delta_0` proves the converse. `square`

The exact case is therefore not mysterious: stable tropical smoothing
kernels are precisely Lawvere metrics (subadditive costs).  The next theorem
is the point relevant to approximate feature algebras: exact idempotence is
not necessary for a **uniform-in-depth** error bound.

## 2. The sharp arbitrary-depth defect law

Use the sup norm on finite profiles.  Comparisons below remain meaningful
for **proper** extended profiles (finite at least at one point), because
convolution with the finite kernel then produces finite outputs.  The
identically infinite profile is excluded.  Define

```math
\Delta(b)=\|b-b_*\|_\infty.                       \tag{TDS.7}
```

### Theorem TDS.2 (tropical defect saturation)

For every integer `m>=1`,

```math
\sup_{f_1,\ldots,f_m}
\left\|
  (P_bf_1)\star\cdots\star(P_bf_m)
  -P_b(f_1\star\cdots\star f_m)
\right\|_\infty
=\|b^{\star m}-b\|_\infty.                      \tag{TDS.8}
```

Consequently,

```math
\sup_{m\ge1}\ \sup_{f_1,\ldots,f_m}(\cdots)
=\boxed{\Delta(b)}.                              \tag{TDS.9}
```

The same bound holds if `P_b` is applied again at arbitrary internal nodes
of a composition tree.  Thus:

* a one-shot use has its declared approximation error;
* after a bounded number of factors the exact additional algebra defect is
  `||b^(star m)-b||_infinity`;
* after arbitrarily many factors the defect does **not** grow linearly: it
  saturates exactly at the distance from `b` to its subadditive closure.

In particular, a non-subadditive kernel with `0<Delta(b)<=eta` gives a
genuinely nonexact approximate homomorphism whose algebra error stays at
most `eta` at every depth.

#### Proof

Associativity and commutativity give

```math
(P_bf_1)\star\cdots\star(P_bf_m)
=(f_1\star\cdots\star f_m)\star b^{\star m}.     \tag{TDS.10}
```

Min-plus convolution is sup-norm nonexpansive in either argument, so the
left side of (TDS.8) is at most the right side.  Taking every
`f_i=delta_0` makes the two compared profiles exactly `b^(star m)` and `b`,
proving equality.

The powers decrease pointwise from `b` to `b_*`, and on finite `G` they
stabilize.  Therefore

```math
\sup_m\|b^{\star m}-b\|_\infty
=\|b_*-b\|_\infty.
```

Every additional application of `P_b` merely increases the convolution
power of `b`; (TDS.9) already takes the supremum over all powers. `square`

### Corollary TDS.2a (from profile defect to extremal response)

Let `C` be a convolution-closed class of finite profiles for which

```math
\|P_bf-f\|_\infty\le\alpha\qquad(f\in C).         \tag{TDS.10a}
```

Start with the coarse profiles `P_bf_i`, combine them in any binary tree,
and optionally reapply `P_b` at any internal nodes.  If the resulting state
is decoded by its maximum, its answer differs from

```math
\max_x(f_1\star\cdots\star f_m)(x)
```

by at most `alpha+Delta(b)`, independently of `m` and of the tree depth.
The same conclusion holds for every scalar decoder which is one-Lipschitz
in profile sup norm.

#### Proof

The computed profile differs by at most `Delta(b)` from the one-blur profile
of the exact product, by Theorem TDS.2.  The latter differs from the exact
product by at most `alpha`, by convolution closure and (TDS.10a).  Both the
maximum and any declared one-Lipschitz decoder contract sup norm. `square`

For example, on `Z/3Z` put

```math
b(0)=0,\qquad b(1)=1,\qquad b(2)=2+eta
```

with `0<eta<=1`.  Then `b_*=(0,1,2)` and `Delta(b)=eta`.  This is a concrete
uniformly depth-stable approximate tropical algebra which is not an exact
homomorphic quotient.

### Why this is not merely the triangle inequality

The generic response inequality

```math
d(x_1\star\cdots\star x_m,
  y_1\star\cdots\star y_m)
\le\sum_i d(x_i,y_i)
```

would charge the approximation once per factor.  Theorem TDS.2 computes the
actual joint error before the final absolute value and proves that all those
errors collapse to one shortest-path defect.  The result is sharp for the
complete profile class.  Its classical ingredient is min-plus dynamic
programming; its new use here is a necessary-and-sufficient, arbitrary-depth
response-algebra criterion.

It is also sharply different from convex support-function uncertainty, as
the following exact calculation records.

### Proposition TDS.2b (Minkowski blur has no defect saturation)

Let `K` be a compact convex subset of a finite-dimensional normed space with
`0 in K`, and coarse-grain a compact response body `A` to `A+K`.  Put

```math
R(K)=\max_{x\in K}\|x\|.
```

For every `m>=1`, the worst-case Hausdorff discrepancy between combining
`m` coarse bodies and coarse-graining their exact Minkowski sum once is

```math
d_H(mK,K)=(m-1)R(K).                              \tag{TDS.10b}
```

Hence a compact Minkowski uncertainty body has uniformly bounded
arbitrary-depth algebra defect if and only if `K={0}`.

#### Proof

Support functions turn Minkowski addition into scalar addition, and
Hausdorff distance of compact convex bodies is the sup norm of their support
functions on the dual unit ball.  Since `0 in K`, `h_K>=0`; therefore

```math
d_H(mK,K)
=(m-1)\sup_{\|u\|_*\le1}h_K(u)
=(m-1)R(K).
```

Taking all exact bodies to be `{0}` shows this is the actual worst-case
coarse-algebra defect, not merely an upper bound. `square`

Thus bounded convex error bodies have no nontrivial analogue of shortest-
path defect saturation.  This does not rule out every convex-roof summary;
it rules out the natural scheme which pays approximation by Minkowski
uncertainty and then hopes nonexpansiveness alone prevents accumulation.

## 3. Indicator kernels: exact subgroup rigidity

We now extend the finite-kernel setup to extended-real indicator kernels.
For `K subseteq G` with `0 in K`, let `iota_K` be zero on `K` and infinity
outside it.  Then

```math
iota_K^{\star m}=iota_{mK},                       \tag{TDS.11}
```

where `mK=K+...+K`.  Its Kleene closure is the indicator of the subgroup
generated by `K`.

### Corollary TDS.3 (classification of fixed-resolution min filters)

The min filter

```math
P_Kf(x)=\min_{k\in K}f(x-k)                       \tag{TDS.12}
```

is closed under arbitrary repeated min-plus composition at the **same exact
resolution** if and only if `K` is a subgroup of `G`.

If `K` is not a subgroup, no finite uniform profile error compares the
one-blur and arbitrary-depth states on the complete extended profile class.
For Hamming balls in `F_2^w`,

```math
mB_r=B_{\min\{mr,w\}},                            \tag{TDS.13}
```

and the one-Lipschitz profile `f(x)=|x|` realizes error
`min{mr,w}` between the unblurred `f` and `P_(mB_r)f`.  (The one-blur versus
`m`-blur defect is correspondingly `min{(m-1)r,w-r}` before endpoints.)
Hence landmark or ball smoothing has a real linear-loss regime before
saturation at the full leading scale.

#### Proof

Equation (TDS.11) is immediate from convolution.  In a finite group,
`K+K=K` and `0 in K` hold exactly when `K` is a subgroup: closure under
addition makes `K` a finite submonoid, hence it contains additive inverses.
The complete-profile claim follows from Theorem TDS.2 with indicator
kernels.  Equation (TDS.13) follows by partitioning the support of a vector
of weight at most `mr` into `m` parts of size at most `r`.  Finally,

```math
P_{B_s}f(x)=\max\{|x|-s,0\},
```

so a vector of weight `w` attains error `s`. `square`

This is the promised boundary between the three notions.  A small metric
net is enough against one raw future.  A non-idempotent blur can be safe for
a bounded number of compositions.  Uniform arbitrary-depth use is governed
by the whole kernel-power orbit; for hard indicator blurs it forces an exact
subgroup quotient, while finite penalties permit a genuinely approximate
but depth-stable algebra when their triangle-inequality defect is bounded.

## 4. A strict composable syndrome quotient on a fixed chart

The subgroup case is not only a no-go.  It yields a strict quotient in a
nontrivial covering-radius model.

Let `G=F_2^w`, let `E={e_1,...,e_w}` be the fixed coordinate basis, and let
`S` range over syndrome supports satisfying `E subseteq S`.  Define the word
profile and radius

```math
\lambda_S(x)=\min\{|I|:x=\sum_{s\in I}s,\ I\subseteq S\},
\qquad
\rho(S)=\max_x\lambda_S(x).                      \tag{TDS.14}
```

Fix the coordinate subspace `H=span(e_1,...,e_r)` and store only

```math
\bar\lambda_S(x+H)=\min_{h\in H}\lambda_S(x+h).  \tag{TDS.15}
```

### Theorem TDS.4 (fixed-chart tropical quotient)

For all such supports `S,T`,

```math
\bar\lambda_{S\cup T}
=\bar\lambda_S\star_{G/H}\bar\lambda_T.         \tag{TDS.16}
```

Thus the quotient profiles form an exact associative feature algebra under
arbitrarily many support unions.  Moreover,

```math
\max_{G/H}\bar\lambda_S
\le\rho(S)
\le\max_{G/H}\bar\lambda_S+r.                   \tag{TDS.17}
```

Decoding by the midpoint gives uniform radius error at most `r/2`, before
and after every future composition.  The number of possible stored states is
at most

```math
(w+1)^{2^{w-r}},                                  \tag{TDS.18}
```

so their worst-case description length is at most

```math
2^{w-r}\log_2(w+1)                               \tag{TDS.19}
```

bits.  For `0<epsilon<1/2`, with `r=floor(2*epsilon*w)`, this is
`2^((1-2epsilon)w+o(w))` bits at additive error `epsilon*w`.

#### Proof

Because `E subseteq S`, the profile `lambda_S` is one-Lipschitz in Hamming
distance:

```math
|\lambda_S(x)-\lambda_S(y)|
\le\lambda_S(x-y)\le|x-y|.                      \tag{TDS.20}
```

The support-union identity is

```math
\lambda_{S\cup T}=\lambda_S\star_G\lambda_T:    \tag{TDS.21}
```

partition a shortest union representation according to which support
supplies each atom; the reverse inequality concatenates two
representations.  The map in (TDS.15) is convolution with the subgroup
indicator.  Lemma TDS.1, or a direct rearrangement of the two coset minima,
therefore gives (TDS.16).

Equation (TDS.20) says that the values inside one `H`-coset have oscillation
at most the Hamming diameter `r` of `H`, proving (TDS.17).  Each quotient
profile has `2^(w-r)` integer entries in `{0,...,w}`, which gives
(TDS.18)--(TDS.19). `square`

This does **not** solve the arbitrary-support syndrome question.  Its common
fixed basis is exactly what supplies a common low-diameter subgroup chart.
For an arbitrary spanning support, the known Lipschitz chart is selected
from that support; two children generally select incompatible charts, and
their subgroup minima cannot be convolved without additional information.
The theorem therefore isolates rather than hides the remaining obstruction:

> arbitrary supports need either a common low-diameter quotient chart, or a
> controlled algebra for changing charts.  A state-dependent one-shot net
> by itself supplies neither.

## 5. Cross-model reading and stopping judgment

Theorem TDS.2 applies verbatim to finite-state shortest-path kernels,
min-plus transfer matrices, exact sparse-synthesis cost profiles, and code
coset-leader profiles.  In each case, `b_*` is the cost after arbitrary path
refinement and `Delta(b)` is the exact information loss incurred by replacing
the refined transfer law with one coarse step.  The theorem therefore has a
model outside coding without changing its statement.

The result is generative but scoped:

* it disproves the claim that every nonexact composable profile incurs a
  linearly growing error;
* it proves that the relevant replacement is not generic nonexpansiveness
  but bounded distance to an algebraically closed kernel;
* it classifies indicator/landmark blurs and explains why arbitrary-depth
  stability there really does demand an exact quotient;
* it gives an actual strict composable state when a common chart exists.

The next theorem should not optimize another fixed chart.  It should decide
whether a bounded-complexity **chart-change cocycle** can transport subgroup
quotient profiles between the support-selected basis charts, with a tropical
defect bounded uniformly under repeated union.  A counterexample showing
that such transitions encode the complete basis-change orbit would close
this mechanism for unrestricted syndrome supports.
