# Recognition and sharp stability of the distance-shell algebra

**Scope.**  This note treats the finite min-plus recognition question in
isolation.  Its main input is a real square kernel, not a pre-existing
metric.  The exact classification is classical in spirit (idempotent
min-plus distance kernels are Lawvere metrics).  The useful quantitative
point is sharper: a symmetric hollow kernel with idempotence defect
`delta` is always within `delta` of a pseudometric, with the best
alphabet-dependent constant, even though its one-sided Kleene closure can
be `Theta(q delta)` away.  The note also identifies the extra invariance
test needed for a *composable twisted* shell; a bijective zero pattern by
itself is not enough.

Throughout, for kernels on a finite set `X`,

```math
(K\star L)(i,j)=\min_{k\in X}\{K(i,k)+L(k,j)\},
\qquad \|K\|_\infty=\max_{i,j}|K(i,j)|.          \tag{R.1}
```

Let `PM(X)` be the cone of finite pseudometrics on `X`: symmetric,
nonnegative, hollow kernels satisfying the triangle inequality.

## 1. Exact recognition

### Theorem R.1 (normalized idempotents are pseudometrics)

Let `K` be symmetric and hollow.  Then

```math
K\star K=K
\quad\Longleftrightarrow\quad
K\in PM(X).                                      \tag{R.2}
```

No separate nonnegativity assumption is required.

#### Proof

If `K` is a pseudometric, the triangle inequality gives
`K(i,j)<=K(i,k)+K(k,j)` for every `k`, while `k=i` gives equality in the
minimum.

Conversely, idempotence says every term in the minimum is at least
`K(i,j)`, hence gives every triangle inequality.  On the diagonal it says

```math
0=\min_k\{K(i,k)+K(k,i)\}=2\min_k K(i,k),        \tag{R.3}
```

so all entries are nonnegative. `square`

This is the exact recognition test for the *untwisted* distance kernel.  A
twist needs a further test, recorded in Section 4.

## 2. Sharp quantitative stability

For a symmetric hollow `K` on `q=|X|` points, define its triangle/idempotence
defect

```math
\tau(K)=\max_{i,j,k}\{K(i,j)-K(i,k)-K(k,j)\}.    \tag{R.4}
```

Because the choices `k=i` and `k=j` show `K star K<=K`,

```math
\tau(K)=\|K-K\star K\|_\infty.                 \tag{R.5}
```

In particular, diagonal triples include the possible negativity defect:

```math
K(i,j)\ge-\tau(K)/2.                             \tag{R.6}
```

Put, for `q>=2`,

```math
c_q=\max\left\{\frac12,\frac{q-2}{q}\right\}.
                                                               \tag{R.7}
```

### Theorem R.2 (sharp worst-case metric repair)

Every symmetric hollow `q`-point kernel satisfies

```math
\boxed{
 \inf_{d\in PM(X)}\|K-d\|_\infty
 \le c_q\,\tau(K)\le\tau(K).}                  \tag{R.8}
```

For every `q>=2`, the coefficient `c_q` is the best possible universal
coefficient.  Thus approximate min-plus idempotence has a dimension-free
two-sided repair theorem.  The optimal constant tends to one, rather than
growing with the alphabet.

#### Constructive proof of the upper bound

Write `tau=tau(K)` and `c=c_q tau`.  Make a complete undirected graph with
edge lengths

```math
W(i,j)=K(i,j)+c\quad(i\ne j),\qquad W(i,i)=0.    \tag{R.9}
```

Equation (R.6) and `c>=tau/2` make every edge nonnegative.  Let `d` be the
shortest-path pseudometric of `W`.  The direct edge gives

```math
d(i,j)\le K(i,j)+c.                              \tag{R.10}
```

For any path `v_0=i,...,v_l=j`, repeated use of the relaxed triangle
inequality gives

```math
K(i,j)\le\sum_{r=1}^l K(v_{r-1},v_r)+(l-1)\tau. \tag{R.11}
```

A shortest nonnegative path can be chosen simple, so `l<=q-1`.  Therefore
its `W`-length is at least

```math
K(i,j)+c+(l-1)(c-\tau)
\ge K(i,j)+(q-1)c-(q-2)\tau
\ge K(i,j)-c,                                    \tag{R.12}
```

where the last inequality is `q c >= (q-2)tau`.  Together with (R.10),
this proves (R.8). `square`

The construction is algorithmic: shift all off-diagonal entries by `c`,
then take all-pairs shortest paths.

#### Sharp examples

For `q=2,3`, take one entry `K(1,2)=K(2,1)=-tau/2` and every other
off-diagonal entry zero.  Its defect is exactly `tau`, while every
pseudometric is nonnegative, forcing error at least `tau/2=c_q tau`.

For `q>=4`, fix `L>=tau` and order the points `1,...,q`.  Put

```math
K(i,j)=L|i-j|-\tau\quad(i\ne j),\qquad K(i,i)=0. \tag{R.13}
```

Every triangle violation is at most `tau`, and equality holds whenever the
middle index lies between the endpoints.  If `d` is a pseudometric with
`||d-K||_infinity<=epsilon`, then

```math
L(q-1)-\tau-\epsilon
 \le d(1,q)
 \le\sum_{i=1}^{q-1}d(i,i+1)
 \le(q-1)(L-\tau+\epsilon).                     \tag{R.14}
```

Hence `epsilon >= (q-2)tau/q=c_q tau`.  This matches (R.8).

There is also a useful pointwise lower bound.  If `d` is within `epsilon`
of `K`, every violated triangle changes by at most `3epsilon`, so

```math
\inf_{d\in PM(X)}\|K-d\|_\infty
\ge\max\left\{\tau(K)/3,\ -\min_{i,j}K(i,j)\right\}.          \tag{R.15}
```

It is not always sharp, as the chain (R.13) demonstrates: mutually
incompatible local repairs force the error up to nearly `tau`.

### Robust normalization corollary

The symmetric-hollow hypothesis is numerically checkable and can be
reached stably.  For an arbitrary real square kernel `D`, let

```math
a=\max\left\{
 \max_i|D(i,i)|,\ \frac12\|D-D^T\|_\infty
\right\},
\qquad \delta=\|D\star D-D\|_\infty.            \tag{R.16}
```

Let `S` be obtained by zeroing the diagonal and averaging each off-diagonal
pair.  Then `||S-D||_infinity<=a`; min-plus multiplication is two-Lipschitz
when both factors move, so

```math
\|S\star S-S\|_\infty\le\delta+3a.             \tag{R.17}
```

Theorem R.2 therefore gives

```math
\operatorname {dist}_\infty(D,PM(X))
\le a+c_q(\delta+3a).                            \tag{R.18}
```

This constant is only a convenient robust corollary; no sharpness is
claimed outside the normalized class.

## 3. Why one-sided closure can still cost `Theta(q tau)`

The two-sided repair in Theorem R.2 must be distinguished from insisting on
the canonical minorant.  If `K>=0`, its min-plus Kleene/shortest-path
closure `K_*` obeys

```math
0\le K(i,j)-K_*(i,j)\le(q-2)\tau(K).             \tag{R.19}
```

Indeed, a simple shortest path has at most `q-1` edges, and (R.11) applies.
The chain (R.13) attains equality at `(1,q)`:

```math
K(1,q)-K_*(1,q)=(q-2)\tau.                       \tag{R.20}
```

Thus the earlier tropical defect/Kleene-closure law and Theorem R.2 answer
different questions:

* a prescribed one-sided blur pays a sharp linear-in-`q` defect;
* recognition up to changing the kernel on both sides has a sharp
  dimension-free defect.

A simple condition restores dimension-free *one-sided* repair.  If every
`K_*`-geodesic has a representative with at most `H` edges, then

```math
\|K-K_*\|_\infty\le(H-1)\tau(K).                \tag{R.21}
```

For example, if all off-diagonal entries lie in `[m,D]` with `m>0`, one
may take `H<=floor(D/m)` (or one for a direct minimizer), since a shortest
path costs at most the direct edge.  This is a genuine geometric
bounded-hop hypothesis, not a consequence of approximate idempotence.

## 4. Twisted distance shells: the missing isometry test

Suppose a kernel `K:X times X->R` has a bijective row-zero selection
`g:X->X`, meaning

```math
K(a,g(a))=0.                                     \tag{R.22}
```

Untwist the **rows**, not the columns:

```math
D(u,v)=K(g^{-1}(u),v).                            \tag{R.23}
```

Then `D` is hollow.

### Theorem R.3 (single-shell recognition)

The following are equivalent:

1. `D` is symmetric and `D star D=D`;
2. there is a pseudometric `d` on `X` such that

   ```math
   K(a,t)=d(g(a),t)\qquad(a,t\in X).             \tag{R.24}
   ```

The pseudometric is exactly `D`.  If each row of `K` has a unique zero,
then `d` is a metric.

This recognizes a single family of distance cones.  For the fixed-metric
composition law

```math
K_g\star K_h=K_{h\circ g}                        \tag{R.25}
```

one must additionally check

```math
d(g(u),g(v))=d(u,v)                              \tag{R.26}
```

for every permitted twist `g` (and similarly for `h`).  Bijective zeros
plus symmetry/idempotence after untwisting do **not** imply (R.26).

#### Minimal counterexample to omitting (R.26)

Take the three points at positions `0,1,3` on the real line and let `g`
swap the first two points.  The kernel `K(a,t)=d(g(a),t)` has a bijective
zero pattern and untwists to the exact line metric.  Nevertheless direct
calculation gives

```math
\|K_g\star K_g-K_{g^2}\|_\infty=1.              \tag{R.27}
```

The failure is precisely that the swap is not an isometry.  Thus the
candidate recognition rule is correct for one shell but incomplete for the
composable algebra.

When the twists are isometries, (R.25) follows from the triangle inequality:

```math
d(g(a),t)+d(h(t),u)\ge d(hg(a),u),               \tag{R.28}
```

and equality is attained at `t=g(a)`.  Strengths give the familiar
bottleneck law

```math
(\lambda d)_g\star(\mu d)_h
=(\min\{\lambda,\mu\}d)_{h\circ g}.             \tag{R.29}
```

## 5. Approximate recognition with exact composability

For a finite permutation group `Gamma` acting on `X`, define the full orbit
defect of a symmetric hollow `D` by

```math
\omega_\Gamma(D)=
\max_{\gamma\in\Gamma,u,v}
|D(\gamma u,\gamma v)-D(u,v)|.                  \tag{R.30}
```

### Theorem R.4 (invariant metric repair)

Let `tau=||D-D star D||_infinity`.  There is a `Gamma`-invariant
pseudometric `d` satisfying

```math
\boxed{
\|D-d\|_\infty\le\omega_\Gamma(D)+c_q\tau.}    \tag{R.31}
```

In particular, if the approximate kernel is already exactly invariant
under all permitted twists, the sharp repair construction can be made
invariant with error `c_q tau`.

#### Proof

Average over the group:

```math
\overline D(u,v)=\frac1{|\Gamma|}
 \sum_{\gamma\in\Gamma}D(\gamma u,\gamma v).    \tag{R.32}
```

This is symmetric, hollow, `Gamma`-invariant, and within `omega_Gamma(D)`
of `D`.  Averaging the relaxed triangle inequalities shows
`tau(overline D)<=tau`.  Apply the shifted-shortest-path construction of
Theorem R.2 to `overline D`, using `c_q tau`.  Both the shifted complete
graph and its shortest-path metric are `Gamma`-invariant. `square`

After row untwisting as in (R.23), Theorem R.4 is an explicit certificate
that `K` is close to an exactly composable distance-shell kernel.  Combining
it with min-plus nonexpansiveness yields cumulative error equal to the sum
of the individual recognition errors; it gives sublinear cumulative loss
exactly when that sum is sublinear at the declared response scale.

The full-orbit condition cannot in general be replaced by a dimension-free
one-generator defect.  If `g` has finite order `L` and

```math
\eta=\max_{u,v}|D(gu,gv)-D(u,v)|,                \tag{R.33}
```

then telescoping gives only

```math
\omega_{\langle g\rangle}(D)
\le\lfloor L/2\rfloor\eta.                      \tag{R.34}
```

This linear dependence is necessary.  For even `L`, take
`X=Z/LZ`, `g(i)=i+1`, and

```math
a_i=A+\frac\eta2\min\{i,L-i\},\qquad
d(i,j)=a_i+a_j\ (i\ne j),\quad d(i,i)=0,         \tag{R.35}
```

where `A>0`.  This is an exact star metric and

```math
\max_{i,j}|d(g i,g j)-d(i,j)|\le\eta.           \tag{R.36}
```

But any `g`-invariant pseudometric has constant values on the adjacent
pairs `(i,i+1)`, whereas the corresponding values of `d` have range
`eta(L-2)/2`.  Hence

```math
\inf_{\rho:\ g\text{ is a }\rho\text{-isometry}}
\|d-\rho\|_\infty\ge\frac{\eta(L-2)}4.          \tag{R.37}
```

Uniform orbit control (or bounded-order holonomy) is therefore the natural
dimension-free hypothesis for exact invariantization.

## 6. Consequences and limits for the response theory

1. **A checkable recognition certificate now exists.**  Row-zero
   untwisting, symmetry, approximate idempotence, and full-group invariance
   certify closeness to the metric-isometry/bottleneck algebra without
   first guessing the metric.

2. **The new content is not the exact Lawvere-metric equivalence.**  That is
   classical.  The useful theorem-level addition is the sharp two-sided
   stability constant, its separation from one-sided Kleene repair, and the
   invariant repair criterion.

3. **Approximate idempotence alone does not compress a generic kernel.**
   The repaired metric can still contain `q(q-1)/2` independent distances.
   Compression occurs only when the metric is fixed or belongs to a smaller
   structured family and the isometry/holonomy has a compact description.

4. **Directed rows behave correctly once recognized.**  For an isometric
   shell,

   ```math
   \sup_t\{d(g(a),t)-d(g(b),t)\}=d(a,b),          \tag{R.38}
   ```

   and the bottleneck composition preserves this directed table.  The
   holonomy remains necessary for labelled endpoint futures, so the
   directed table alone is not the complete contextual state.

5. **The sharp chain is an adversarial warning.**  Local triangle defects
   can accumulate by a factor `q` in a fixed one-sided closure even while a
   different exact metric lies within one local defect.  Any theorem using
   the Kleene representative must declare why one-sidedness is required;
   any theorem allowing model replacement should use the two-sided repair.

## 7. Independent numerical sanity check

For `q=2,...,8`, a linear program minimizing sup-norm distance to the
pseudometric cone was run on the two sharp families above.  It returned

```text
q       2     3     4     5       6       7      8
opt/tau .5    .5    .5    .6      2/3     5/7    .75
c_q     .5    .5    .5    .6      2/3     5/7    .75
```

This is only an implementation check; equations (R.9)--(R.14) are the
proof.

## Director recommendation

The exact/approximate recognition problem is solved at the natural finite
kernel level.  Promote Theorems R.2 and R.4 if the campaign wants a
recognition/stability theorem adjoining the metric-isometry bottleneck
algebra.  Do **not** promote the incomplete slogan “bijective zero pattern
plus untwisted idempotence recognizes the composable algebra”: it recognizes
distance cones, but the isometry condition is independently necessary.

The strongest next question is structural rather than another metric
repair estimate:

> Which model-generated transfer kernels have row-zero untwist defect and
> full-orbit defect `o(1)` while their repaired metric belongs to a
> subquadratic-description family?

Without the final complexity condition, recognition is exact but does not
yet explain compression.

## 8. Adversarial audit of the parallel AMR.1--AMR.2 draft

I checked `extremal_information/drafts/approximate_metric_recognition.md`
after completing the derivation above.

* **AMR.1 is correct.**  Its repair `K+delta` off the diagonal is the
  simpler, nonoptimal special case of Theorem R.2.  The claim that it is a
  metric for `delta>0` is justified by `K(i,j)>=-delta/2`, so every distinct
  pair has repaired distance at least `delta/2`.  At `delta=0` one obtains a
  pseudometric, as stated.

* **AMR.2's exact power formula is correct.**  In a `T`-factor path, let
  `p<=T` be the number of nonzero moves and `V` their total integer
  variation.  The cost is `aV-delta p`.  If `T<=|i-j|`, then `p<=T` and
  `V>=|i-j|`, so the minimum uses `p=T` monotone pieces.  If
  `T>|i-j|`, a monotone unit path with `p=|i-j|` is optimal, because any
  `p>|i-j|` has cost at least `(a-delta)p>(a-delta)|i-j|`.  Zero stays fill
  unused factors.  This proves the displayed `min{T,|i-j|}` formula.

* **The projective drift constant is correct.**  In row zero the pointwise
  difference is zero at distances zero and one and equals
  `-(T-1)delta` at every distance at least `T`.  Its oscillation is therefore
  `(T-1)delta`, and quotienting by constants divides it by two.

* **One sentence needs qualification.**  A bijective zero-centre map plus
  exact symmetry/idempotence after row untwisting recognizes a single
  twisted distance shell, but not the fixed-metric *composition algebra*.
  The zero permutation must additionally be an isometry.  Equation (R.27)
  is a three-point counterexample.

* **The two drafts complement rather than contradict one another.**  AMR.2
  shows that repeatedly composing the original approximate kernel incurs
  linear drift.  R.2 says one may instead replace it once by an exact metric
  within the sharp error `c_q delta`.  These are different operations.  If
  the application is required to use the original kernel, AMR.2 is the
  decisive no-go; if model replacement is permitted, R.2 is the stronger
  recognition theorem.
