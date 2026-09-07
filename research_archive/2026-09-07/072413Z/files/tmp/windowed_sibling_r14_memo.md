# Wave 14 memo: disjoint sibling energy and its causal limit

## 1. Exact completion identity

Write an endpoint edge as

```math
A[U]=\begin{pmatrix}D&B\\B^{\mathsf T}&X\end{pmatrix}
```

and fix its oriented child state `(sigma,y)`.  Put

```math
L=\lVert By\rVert_1,\qquad h=y^{\mathsf T}Xy,\qquad
\partial=Q(U)-Q(X),
```

```math
c=[2L-(Q(X)-\sigma h)]_+,
\qquad r=[c-\partial]_+.
```

Assume `r>0`, and choose `u` on the sibling so that
`u_i(\sigma By)_i=|(By)_i|`.  Define the full-state deficit

```math
e=Q(U)-\sigma (u,y)^{\mathsf T}A[U](u,y)\ge0.
```

Then the positive parts are inactive and direct expansion gives the exact
identity

```math
\boxed{-\sigma u^{\mathsf T}Du=r+e.}
```

Thus (10.533) discards a nonnegative, explicitly identified amount: the
negative sibling energy is residual plus completion deficit.  This is a
proof, not a conditional statement.

## 2. Exact disjoint-sibling aggregation

Take any subfamily of positive-residual buckets on one oriented endpoint
chain.  Their sibling blocks are disjoint.  Let `W_+` be the union of
siblings whose bucket has `sigma=+1`, and `W_-` the analogous union for
`sigma=-1`.  Independently flip the entire witness on each sibling block.
All inter-block quadratic terms have expectation zero, while internal
energies are unchanged.  Therefore some state on `W_+` has energy at most
`-sum_{sigma=+}(r_b+e_b)`, and some state on `W_-` has energy at least
`sum_{sigma=-}(r_b+e_b)`.  Consequently

```math
\boxed{
\sum_{b:\sigma_b=+}(r_b+e_b)\le N(A[W_+]),\qquad
\sum_{b:\sigma_b=-}(r_b+e_b)\le P(A[W_-]).
}
```

Here `P(M)=max_x x^T Mx` and `N(M)=-min_x x^T Mx`.  Random extension from
either shore to `W=W_+ disjoint-union W_-` gives

```math
\boxed{
\sum_b(r_b+e_b)
\le N(A[W_+])+P(A[W_-])
\le \mathcal R(A[W]).
}
```

The second inequality does not assert superadditivity of `Q`; it uses the
separate monotonicities of `P` and `N`.  This distinction is essential when
the two orientations occur on the same chain.

There is an exact hypergraph formulation.  Let buckets be the vertices and
oriented root-to-leaf chains be the hyperedges.  For loads
`ell_b=a_b/c_b` (zero when `c_b=0`), its fractional path-cover number is

```math
K_{path}(\ell)=
\min_{\nu_\pi\ge0}
\left\{\sum_\pi\nu_\pi:
\sum_{\pi\ni b}\nu_\pi\ge\ell_b\quad(\forall b)\right\}.
```

Finite LP duality gives

```math
K_{path}(\ell)=
\max_{\eta_b\ge0}
\left\{\sum_b\ell_b\eta_b:
\sum_{b\in\pi}\eta_b\le1\quad(\forall\pi)\right\}.
```

Section 10.62 proves `K_path(ell)<=4`.  Hence, for any selected bucket set
`B_0`, give a chain the exact signed-profile cost

```math
\gamma_\pi=N(A[W_{\pi,+}])+P(A[W_{\pi,-}]).
```

The strongest LP consequence of the chain theorem is the profile-weighted
cover

```math
\begin{aligned}
C_\Gamma(\ell)
&=\min_{\nu_\pi\ge0}
\left\{\sum_\pi\gamma_\pi\nu_\pi:
\sum_{\pi\ni b}\nu_\pi\ge\ell_b\quad(\forall b\in B_0)\right\}\\
&=\max_{\eta_b\ge0}
\left\{\sum_{b\in B_0}\ell_b\eta_b:
\sum_{b\in\pi\cap B_0}\eta_b\le\gamma_\pi\quad(\forall\pi)\right\}.
\end{aligned}
```

Integration against a feasible path cover and finite LP duality prove

```math
\boxed{
\sum_{b\in B_0}\ell_b(r_b+e_b)
\le C_\Gamma(\ell)
\le K_{path}(\ell)\,
\max_\pi\left[N(A[W_{\pi,+}])+P(A[W_{\pi,-}])\right]
\le4\Gamma_A(B_0).
}
```

This is the strongest unconditional window theorem obtained here.  If the
selected edges form a consecutive order window deleting in total at most
`h` vertices on every chain, then

```math
\Gamma_A(B_0)
\le\min\{\mathcal R(A),h(h-1)\}.
```

The crude `4h(h-1)` consequence already follows from (10.514) by summing
the individual sibling bounds.  What is new is the signed induced-submatrix
profile and the simultaneous control of the completion deficits `e_b`.
Without an additional structural estimate on that profile, this does not
improve the asymptotic error scale.

## 3. Why this does not yet control causal Hall deficiency

For temporal demands `z_i`, bucket capacities `ell_b r_b`, and an explicit
matrix-derived compatibility relation `i~b`, the exact service LP is

```math
\max_{f\ge0}\sum_{i,b}f_{ib}
```

subject to

```math
\sum_bf_{ib}\le z_i,\qquad
\sum_if_{ib}\le\ell_br_b,\qquad
f_{ib}=0\ \text{if }i\not\sim b.
```

Its uncovered demand is exactly

```math
\max_{J\subseteq I}
\left[z(J)-\sum_{b\in N(J)}\ell_br_b\right]_+.
```

For descendant compatibility this reduces to the laminar/antichain formula
in (10.543).  The disjoint-sibling theorem is an **upper** bound on the
second term.  Hall control needs a lower bound on compatible service, so the
new theorem alone has the wrong direction.  It becomes useful only
conditionally, after a temporal-to-endpoint matching theorem has already
been proved: it then bounds the total cost of the matched service in an
order window.

Nor can `r_b+e_b` simply be declared a localized negative temporal credit.
The temporal sign is that of `lambda-partial`, whereas the identity above
contains no `lambda` and imposes no relation between `r`, `e`, and
`lambda-partial`.

## 4. Exact finite obstruction: A5 is field-proportional, not abstract

For the order-five minimizer in (10.518), choose

```math
p=(1,-1,-1,-1,-1),\qquad n=(1,-1,-1,1,1),
```

retain `X={0,1,2}`, and discard `D={3,4}`.  For orientation `sigma=-1`,

```math
Q(U)=8,\quad Q(X)=6,\quad h=-2,\quad L=4,
\quad c=4,\quad\partial=2,\quad r=2.
```

The cross-aligned sibling is `u=(1,1)`.  It has `u^T D u=2`; its full
completion has energy `-8`, so

```math
e=0,\qquad -\sigma u^T D u=r+e=2.
```

Nevertheless this deletion has strictly positive centered temporal demand

```math
\lambda-\partial
=8\left[1-(3/5)^{3/2}\right]-2
=6-\frac{24\sqrt{15}}{25}>0.
```

This is reachable by the field-proportional rule.  In the `p` gauge the row
profile is `(0,2,0,4,2)`.  With eligible set `D`, vertices 3 and 4 are
deleted independently with probabilities `1` and `1/2`; hence deletion of
all of `D` has probability `1/2`.  The other possible outcome deletes only
vertex 3 and also has positive centered demand
`8-64 sqrt(5)/25>0`.  Thus the parent conditional average has positive
demand and no negative atom.  Even a zero-deficit/full-ground sibling
certificate is therefore not a negative temporal credit.

This does not obstruct treating the bucket as positive service in the flow
LP, and it does not defeat endpoint-tie averaging with relaxed retained-set
compatibility.

Two complementary walls remain exact:

- in A8 every endpoint pair has zero residual, while the positive-demand
  shore deletion in (10.539) is a positive-probability field-proportional
  outcome, so no positive-coefficient error-free service theorem survives;

- in the A7 deterministic suffix (10.544), the matching endpoint pair at
  the order-six child has zero residual in both orientations for retaining
  the named order-five child.  The suffix demand is positive, and its only
  compensating temporal credit occurred at the ancestor, which descendant
  causality forbids it from using.

## 5. Verdict

The completion identity and signed disjoint-sibling aggregation are
**proved**.  The path-cover/profile bound is a proved conditional costing
tool once a temporal matching exists.  Conversion into localized negative
temporal credit is **falsified**, even for a positive-probability
field-proportional A5 outcome with `e=0`.  Direct control of the maximum
antichain deficiency remains **open** and requires genuinely new lower
control of compatible residual service or a separate conditional comparison;
random relative sibling signs supply neither.

Exact checker: `tmp/check_windowed_sibling_r14.py`.
