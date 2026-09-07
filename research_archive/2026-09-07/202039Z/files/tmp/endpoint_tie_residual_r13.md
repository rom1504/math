# Endpoint ties and decrement-tolled service: exact minimax and finite wall

Status: the identities, LP duality, and finite certificates below are exact.
The exhaustive small-order table is computational but uses exact integer and
rational checks.  The result is a falsification, not a proof of the quadratic
signing limit: optimizing or averaging endpoint ties does not guarantee any
positive residual service, even for a positive centered demand.  No
large-order globally minimizing lift of the finite wall is proved.

## 1. Residual capacity has an exact matrix witness

Consider any split

```math
U=
\begin{pmatrix}
D&B\\
B^{\mathsf T}&X
\end{pmatrix}
```

and an oriented child state `(sigma,y)`, where `sigma` is `+1` or `-1` and
`y` is a sign vector on `X`.  Put

```math
h=y^{\mathsf T}Xy,
\qquad L=\lVert By\rVert_1,
\qquad \partial=Q(U)-Q(X)\ge0.
```

Its raw all-successor capacity and decrement-tolled residual are

```math
c=[2L-(Q(X)-\sigma h)]_+,
\qquad r=[c-\partial]_+.
```

The elementary identity `[[x]_+-d]_+=[x-d]_+` for nonnegative `d` gives

```math
\boxed{
r=[2\lVert By\rVert_1+\sigma y^{\mathsf T}Xy-Q(U)]_+.
}
\tag{R13.1}
```

Thus decrement tolling removes the child norm exactly.  This is not merely a
scalar comparison: every positive residual has an explicit negative-energy
sibling witness.  Choose

```math
u_i=\operatorname{sign}(\sigma(By)_i),
```

with either sign allowed when `(By)_i=0`.  Then
`sigma u^TBy=L`.  Applying the parent norm bound to the full completion
`(u,y)` gives

```math
Q(U)
\ge \sigma\left(u^{\mathsf T}Du+2u^{\mathsf T}By+y^{\mathsf T}Xy\right),
```

and hence

```math
\boxed{
\sigma u^{\mathsf T}Du
\le Q(U)-2L-\sigma h.
}
\tag{R13.2}
```

In particular, if `r>0`, the expression inside (R13.1) equals `r`, so

```math
\boxed{\sigma u^{\mathsf T}Du\le-r<0.}
\tag{R13.3}
```

This both reproves `r\le Q(D)` and sharpens it: positive residual service is
certified by a specified orientation and signing of the sibling, not only by
its norm.

For an endpoint pair at a node `U`, switch the positive endpoint to `1` and
let `U=S\sqcup T` be the positive/negative endpoint cut.  Write

```math
\mathcal R(U)=P(U)+N(U),
\qquad I(U)=|P(U)-N(U)|,
\qquad h_X=1_X^{\mathsf T}U[X]1_X.
```

The exact endpoint exposure theorem gives
`L_X=\mathcal R(U)/4` on both shores.  Since
`\mathcal R(U)/2-Q(U)=-I(U)/2`, (R13.1) specializes to

```math
\boxed{
r_X^\sigma
=\left[\sigma h_X-\frac{I(U)}2\right]_+,
\qquad
\sum_{\sigma=\pm}r_X^\sigma
=\left[|h_X|-\frac{I(U)}2\right]_+.
}
\tag{R13.4}
```

Endpoint-tie optimization is therefore exactly an optimization of the two
internal endpoint energies above the parent half-imbalance threshold.  Child
norms play no role in the residual itself.  At a balanced parent the total
root residual for a pair is simply `|h_S|+|h_T|`.

## 2. The exact joint allocation/transport LP

Fix a complete endpoint tree `tau`.  Let `b=(v,X,sigma)` run over its
positive-capacity buckets, with capacity `c_b`, and put

```math
\kappa_b=\frac{[c_b-(Q(v)-Q(X))]_+}{c_b}.
```

The allocated residual is `kappa_b a_b`.

Let `H_tau` be the full `K\le4` conserved-allocation polytope.  It has incoming
obligations `w_v`, allocations `a_b`, and the path-cover variables from
Section 10.59.  Explicitly,

```math
\begin{aligned}
w_{\rm root}&=\mathcal R(A),\\
w_v&=\sum_{b\text{ based at }v}a_b
       +\sum_{X\text{ child of }v}w_X,\\
0\le a_b&\le c_b,
\qquad 0\le w_X\le I(X),\\
\ell_{vX}&=\sum_{\sigma:c_{vX}^\sigma>0}
             \frac{a_{vX}^\sigma}{c_{vX}^\sigma},\\
z_{vX}&\ge\ell_{vX},
\qquad z_{vX}\ge\theta_X,\\
\theta_v&\ge\sum_{X\text{ child of }v}z_{vX},
\qquad \theta_{\rm leaf}=0,
\qquad \theta_{\rm root}\le4.
\end{aligned}
\tag{R13.5}
```

Existence of a point in `H_tau` for every endpoint tree is Section 10.62.
The particular half-stopping construction is one point of this larger
polytope; optimizing over `H_tau` can improve residual service while retaining
the same proved path-cover bound.

Let `I` be a prescribed finite family of centered temporal demand atoms of
nonnegative sizes `z_i`, and let `i\sim_\tau b` be any explicitly stated
matrix-derived compatibility relation.  For a prescribed tree, the exact
maximum service is the LP

```math
\begin{aligned}
S_\tau(z)=\max\quad&\sum_{i,b}\gamma_{ib}\\
\text{subject to}\quad
&\sum_{b:i\sim_\tau b}\gamma_{ib}\le z_i &&(i\in I),\\
&\sum_{i:i\sim_\tau b}\gamma_{ib}\le\kappa_ba_b &&(b\in B_\tau),\\
&(a,w,z^{\rm edge},\theta)\in H_\tau,\\
&\gamma_{ib}\ge0.
\end{aligned}
\tag{R13.6}
```

This is a genuine joint allocation/transport optimization.  It is stronger
than first fixing the Section 10.62 allocation and then applying a Hall test.

For an exact dual, define the weighted residual support function

```math
H_\tau(\omega)
=\max_{(a,\ldots)\in H_\tau}
\sum_{b\in B_\tau}\omega_b\kappa_ba_b.
\tag{R13.7}
```

For `0\le\alpha_i\le1`, put

```math
\omega_{\tau b}(\alpha)
=\max_{i:i\sim_\tau b}(1-\alpha_i),
```

with maximum zero for a bucket with no neighbor.  The ordinary bipartite
capacity dual followed by optimization over `H_tau` gives

```math
\boxed{
S_\tau(z)
=\min_{0\le\alpha\le1}
\left\{
\sum_i z_i\alpha_i+H_\tau(\omega_\tau(\alpha))
\right\}.
}
\tag{R13.8}
```

There is no hidden nonlinear optimization in `H_tau`.  If (R13.5) is written
as

```math
E_\tau x=e_\tau,
\qquad G_\tau x\le g_\tau,
\qquad x\ge0,
```

and `K_tau x` is the residual-capacity vector, standard LP duality gives

```math
\boxed{
H_\tau(\omega)
=\min_{y\ {m free},\ s\ge0}
\left{e_\tau^{\mathsf T}y+g_\tau^{\mathsf T}s:
E_\tau^{\mathsf T}y+G_\tau^{\mathsf T}s
\ge K_\tau^{\mathsf T}\omega
\right\}.
}
\tag{R13.9}
```

Equations (R13.8)--(R13.9) are the requested finite minimax/dual.  For a
fixed allocation, (R13.8) reduces to the integral cut formula

```math
S(z)=\min_{J\subseteq I}
\{z(I\setminus J)+\widehat a(N(J))\}.
```

The allocation support function in (R13.8) is the additional, genuinely
matrix/tree content.

## 3. Prescribed, optimized, and averaged endpoint claims differ

Let `\mathcal T` be the finite set of complete endpoint trees.

- A **prescribed-tree** claim concerns one value `S_tau`.

- A deterministic **optimized-tree** claim concerns

  ```math
  S_{\rm det}(z)
  =\max_{\tau\in\mathcal T}S_\tau(z)
  =\max_\tau\min_\alpha
   \{z\cdot\alpha+H_\tau(\omega_\tau(\alpha))\}.
  \tag{R13.10}
  ```

- An **averaged-tree** claim permits a random endpoint tree, with allocation
  and transport chosen conditional on the realized tree.  Its primal is

  ```math
  \begin{aligned}
  S_{\rm av}(z)=\max\quad&\sum_{i,\tau,b}\gamma_{i\tau b}\\
  \text{subject to}\quad
  &\sum_{\tau,b}\gamma_{i\tau b}\le z_i,\\
  &\sum_i\gamma_{i\tau b}\le\kappa_{\tau b}a_{\tau b},\\
  &x_\tau\in\pi_\tau H_\tau,
    \qquad \sum_\tau\pi_\tau=1,
    \qquad \pi_\tau\ge0.
  \end{aligned}
  \tag{R13.11}
  ```

  Here `x_tau in pi_tau H_tau` denotes the ordinary perspective/scaled
  version of every linear constraint in (R13.5).  Its dual is

  ```math
  \boxed{
  S_{\rm av}(z)
  =\min_{0\le\alpha\le1}
  \left\{
  z\cdot\alpha+
  \max_{\tau\in\mathcal T}H_\tau(\omega_\tau(\alpha))
  \right\}.
  }
  \tag{R13.12}
  ```

The max-min/min-max distinction between (R13.10) and (R13.12) is exact.
Thus averaging can help a multi-atom compatibility problem, but it cannot
create residual service if every pure endpoint tree has zero service.
For one demand atom adjacent to every bucket, averaging also cannot beat the
best pure tree, except that either may already saturate the demand.

If compatibility requires the endpoint child to equal the temporal retained
vertex set, only matching endpoint trees belong in the relevant set.  Merely
having the same parent order is a much weaker relaxation and must not be
confused with this set-compatible claim.

## 4. Order five: ties rescue only the relaxed compatibility

For the order-five minimizer from (10.518), exact endpoint enumeration gives
five positive and five negative projective endpoints, hence twenty-five
pairs.

- Five pairs are `1+4` splits.  For each, `h_S=h_T=0`, so (R13.4) gives zero
  residual on every root bucket.  Both child imbalances vanish, so every
  conserved allocation has zero descendant allocation as well.

- The other twenty pairs are `2+3` splits.  Up to exchanging shores their
  endpoint data are

  ```math
  (h_S,h_T)=(2,-2),
  \qquad
  (c_S^+,c_S^-)=(8,4),
  \qquad
  (c_T^+,c_T^-)=(0,4).
  ```

  The two positive residuals are two and two, for total raw residual four.
  The canonical Section 10.62 root factor is

  ```math
  t_0=\frac{16}{16+4}=\frac45,
  ```

  so its residual service is `16/5` on each such pair.  A general point of
  (R13.5) attains the full value four: allocate eight and four units in the
  two residual-positive root buckets, pass four units to the imbalanced
  order-three child, and discharge them in one capacity-four bucket there.
  The child bucket has zero residual.  The edge loads are one, one, and one,
  and the exact recursion gives `theta_root=2`.  Total raw residual four is
  the matching dual upper certificate.

Uniform averaging over all twenty-five pairs gives residual capacity

```math
\frac{20}{25}\frac{16}{5}=\frac{64}{25}
```

for the canonical allocation, and `16/5` if the allocation is also optimized.

Now prescribe the field-proportional deletion in (10.518), which deletes
vertex four deterministically.  Its exact centered demand is

```math
z_5=8-\frac{64\sqrt5}{25}>0.
\tag{R13.13}
```

Exactly one endpoint pair has the matching singleton/four-set cut, and it is
one of the five zero-service pairs.  Therefore:

- prescribed matching-tree/set-compatible service is zero;
- optimized or uniformly averaged service under the relaxed rule that every
  root bucket is compatible is positive;
- indeed `64/25>z_5`, so even the uniform canonical average serves all of the
  relaxed scalar demand.

The strict inequality is exact: it reduces to `64 sqrt(5)>136`, whose square
is `20480>18496`.  This finite rescue says nothing about an honest
vertex-set compatibility theorem.

## 5. Order eight: all endpoint optimization and averaging fail

The exact order-eight minimizer

```math
A_8=
\begin{pmatrix}
0&1&1&1&1&1&1&1\\
1&0&1&-1&1&1&-1&-1\\
1&1&0&1&-1&1&-1&-1\\
1&-1&1&0&-1&-1&-1&1\\
1&1&-1&-1&0&-1&1&-1\\
1&1&1&-1&-1&0&1&1\\
1&-1&-1&-1&1&1&0&1\\
1&-1&-1&1&-1&1&1&0
\end{pmatrix}
```

has `P=N=Q=q_8=20`, four projective positive endpoints, and four projective
negative endpoints.  Exact enumeration of all sixteen pairs gives:

```math
|S|=|T|=4,
\qquad h_S=h_T=0,
\qquad P(S)=N(S)=P(T)=N(T)=8=q_4.
\tag{R13.14}
```

Equivalently, every one of the four root buckets has capacity twelve and
decrement twelve.  Thus every raw residual is zero.  Both child imbalance
caps are also zero.  Root conservation forces all forty range units to be
allocated at the root and forces every descendant obligation to vanish.
Consequently

```math
\boxed{
S_\tau(z)=S_{\rm det}(z)=S_{\rm av}(z)=0
}
\tag{R13.15}
```

for every demand whose proposed service uses this Section 10.62 residual
allocation, even under the maximally relaxed compatibility in which the atom
is adjacent to every bucket.  Descendant tie choices cannot escape (R13.15).

Retaining either endpoint shore is a deterministic general deletion to an
exact order-four minimizer.  Its decrement and centered scale toll are

```math
d=20-8=12,
\qquad
\lambda
=\frac{20}{8^{3/2}}(8^{3/2}-4^{3/2})
=20-5\sqrt2.
```

Hence it has the strictly positive demand

```math
\boxed{\lambda-d=8-5\sqrt2>0.}
\tag{R13.16}
```

This proves that no finite constant can guarantee
`centered demand <= constant times optimized/averaged residual service`.
It defeats prescribed-tree, deterministic optimized-tree, and averaged-tree
versions simultaneously.

The scope relative to field-proportional peeling must be stated carefully.
The deterministic four-vertex deletion in (R13.16) is a valid atom for the
general centered identity (10.511), but the deletion is not deterministic
under one field-proportional Bernoulli step.  It is, however, a genuine
positive-probability **outcome** of such a step: take the deleted endpoint
shore as the heavy set.  Every one of its four probabilities `r_i/7` is
positive, so deleting exactly that shore has positive rational probability.
Consequently (R13.16) defeats an outcome-by-outcome residual charging claim.
The other partial-shore outcomes may carry compensating negative credits, so
this observation alone does not defeat a parent-averaged or grouped-credit
claim.

In every positive-ground gauge of `A_8`, the row fields are a permutation of

```math
(1,1,1,3,3,3,3,5),
```

so the singleton probabilities are `r_i/7`, never zero or one.
Nevertheless, a genuine field-proportional attempt with heavy set `{i}` has
only the empty and singleton outcomes.  Every singleton restriction has
`Q=18=q_7`, so its nonempty outcome has

```math
\lambda-d
=20\left[1-\left(\frac78\right)^{3/2}\right]-2
=18-\frac{35\sqrt{14}}8>0.
\tag{R13.17}
```

The parent-level expected atom is `(r_i/7)` times (R13.17), while the empty
outcome contributes zero.  Thus this singleton-heavy-set step has positive
parent-level field-proportional demand and no local negative outcome.  This
is still a local/random-order statement.  By itself it is not a finite
deterministic-terminal-order field-proportional process; repeated rejection
until deletion is only almost surely finite.  No stronger scope is claimed.

## 6. Exhaustive small-minimizer audit

Every switching-normalized signing was enumerated through order seven.  The
table lists the exact minimum `q_n`, the number of labelled switching classes
attaining it, the range over minimizers of the best-pair canonical
Section 10.62 residual, and the corresponding range of best raw residual.

| `n` | `q_n` | minimizing switching classes | best canonical range | best raw range |
|---:|---:|---:|---:|---:|
| 2 | 2 | 1 | 0 | 0 |
| 3 | 6 | 2 | 0 | 0 |
| 4 | 8 | 6 | `8/3` | 4 |
| 5 | 8 | 12 | `16/5` | 4 |
| 6 | 10 | 12 | `10/3` | 4 |
| 7 | 18 | 3,240 | `16/9` to `24/5` | 2 to 6 |

The uniform mean raw residual is respectively `0,0,4,16/5,10/3`, and a
range from `1` to `5/4` at order seven.  Thus endpoint optimization is
nonzero on every global minimizer at orders four through seven.  The named
order-eight minimizer is the first checked competitive example where all
endpoint pairs fail at once.  The order-eight row is a certificate for one
global minimizer, not an exhaustive classification of all order-eight
minimizers.

`tmp/endpoint_tie_residual_r13.py` verifies:

1. (R13.1)--(R13.3) on all 13,992 buckets arising from every split, child
   state, and orientation of every switching-normalized signing through
   order five;
2. (R13.4) on every endpoint bucket of every global minimizer through order
   seven and on all sixteen endpoint pairs of `A_8`;
3. every count and rational value above;
4. the `A_5` matching/optimized/averaged certificates;
5. all `A_8` child norms, imbalances, residuals, row-field profiles, and the
   exact squared inequalities proving positivity of (R13.16)--(R13.17).

## 7. Consequence for the next route

The endpoint-tie minimax is now exact, but its desired positive theorem is
false.  Equation (R13.3) identifies the structural replacement: positive
residual is equivalent to finding a cross-aligned completion with negative
oriented sibling energy.  A viable asymptotic theorem must therefore do at
least one of the following in each macroscopic order window:

1. force enough endpoint pairs with such negative sibling witnesses;
2. prove that endpoint-neutral mosaics such as (R13.14) force compensating
   negative centered temporal increments or terminal excess;
3. allow a tail-summable additive error large enough to absorb finite neutral
   blocks while still vanishing under adaptive descent.

The finite `A_8` wall rules out every error-free, positive-coefficient
endpoint-tie or averaging statement.  It does not rule out an asymptotic
grouped theorem with additive tail-summable error, because no globally
minimizing blow-up of `A_8` is known.
