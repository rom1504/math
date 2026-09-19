# Wave 13, route 2: exact finite causal signed flow

## Outcome

There is an exact finite theorem, but it exposes a new family of conditional
comparisons rather than deriving them from the root telescope.  Under the
primary convention below, a demand at a temporal node may use only resources
located in its descendant subtree.  Feasibility is equivalent to one cut for
every rooted subtree; on a chain these are suffix cuts.  The minimum uncovered
demand is a maximum-weight antichain of rooted-subtree deficiencies.

For centered peeling atoms, every such deficiency has an exact telescoping
formula.  If

```math
\alpha=\frac{q_n}{n^{3/2}},\qquad
V(v)=Q(A_v)-\alpha r_v^{3/2},\qquad
G=q_m-\alpha m^{3/2},
```

then the cut at a reached state `v` is exactly

```math
\pi_v\bigl(G-V(v)\bigr)\le \kappa R(T_v).
```

Here `R(T_v)` is the residual endpoint resource honestly placed in the
descendant subtree.  The root cut is the relaxed scalar comparison.  Every
proper cut is an additional conditional normalized comparison.  Telescoping
computes its two sides, but supplies no sign.  Thus the causal reformulation
is exact and tail-useful once its maximum antichain deficiency is small, but
it does not by itself close the proof.

An exact order-seven minimizer gives a realizable obstruction: its first
deletion can create an exact order-six minimizer and hence a large *early*
negative centered credit, while the later suffix has positive net demand.
The global sink uses that early credit; descendant-only service cannot.

## 1. Abstract descendant-resource theorem

Let `T` be a finite rooted tree.  At node `v` place demand `z_v\ge0` and
resource `s_v\ge0`.  Write `T_v` for the descendants of `v`, including `v`.
The causal rule is

```math
v\leadsto w\quad\Longleftrightarrow\quad w\in T_v:
```

a demand at `v` may draw resource only from `v` or a descendant.  A fractional
flow consists of `f_{vw}\ge0` on these pairs, with

```math
\sum_{w\in T_v}f_{vw}\le z_v,
\qquad
\sum_{v\preceq w}f_{vw}\le s_w.
```

Let `E_{\rm caus}` be the least total uncovered demand.  Ordinary finite
max-flow/min-cut gives

```math
E_{\rm caus}
=\max_{J\subseteq T}
\left[z(J)-s\!\left(\bigcup_{v\in J}T_v\right)\right]_+.
\tag{F1}
```

This Hall family is laminar.  Define

```math
H_v=z(T_v)-s(T_v).
```

If `A` is an antichain, its rooted subtrees are disjoint.  Then

```math
\boxed{
E_{\rm caus}
=\max_{A\ {\rm antichain}}
\sum_{v\in A}H_v,
}
\tag{F2}
```

where the empty antichain is allowed.  Nodes with zero demand may be omitted
from the maximizing antichain.

**Proof.**  Given `J`, let `A` be its ancestry-minimal nodes.  Its
neighborhood is the disjoint union of the `T_v`, `v\in A`.  Enlarging `J`
inside each of these subtrees to all demand nodes can only increase the Hall
deficiency and does not enlarge the neighborhood.  The resulting deficiency
is `\sum_{v\in A}H_v`.  Conversely, for an antichain of positive-demand
nodes, choose all demand nodes in its rooted subtrees.  This realizes the same
quantity.  If a zero-demand node has `H_v>0`, then
`H_v\le\sum_c H_c` over its children, so descending and discarding
nonpositive branches gives at least as large an antichain supported on demand
nodes.  This proves (F2).  It also gives the bottom-up recursion

```math
D_v=\max\!\left\{0,H_v,\sum_{c\text{ child of }v}D_c\right\},
```

with the `H_v` option omitted when `z_v=0`; `D_{\rm root}=E_{\rm caus}`.
This is also a constructive bottom-up greedy proof of sufficiency.

In particular, exact service is possible if and only if

```math
\boxed{z(T_v)\le s(T_v)\quad\text{for every }v.}
\tag{F3}
```

The inequalities at zero-demand nodes are redundant, so stating them for all
nodes is harmless.

The fractional dual has a useful finite ``causal envelope'' form:

```math
\boxed{
E_{\rm caus}
=\max_{0\le\theta_v\le1}
\left[
\sum_v z_v\theta_v
-\sum_w s_w\max_{v\preceq w}\theta_v
\right].
}
\tag{F4}
```

The running maximum `M_w=\max_{v\preceq w}\theta_v` is nondecreasing down
each history.  Layer-cake decomposition of `(\theta,M)` recovers (F1)--(F2).
This is the complete finite dual; no topological or continuous-time
assumptions from causal-transport theory are needed.

## 2. Exact application to centered temporal atoms

Let a finite peeling tree start at an exact order-`n` minimizer.  Every leaf
has order `m`.  A state node `v` is reached with probability `\pi_v`, has
matrix `A_v` and order `r_v`, and has children `w` with conditional
probabilities `p_{vw}`.  Put

```math
\alpha=\frac{q_n}{n^{3/2}},
\qquad
V(v)=Q(A_v)-\alpha r_v^{3/2}.
```

For an edge `v\to w`, the unweighted centered increment is

```math
\xi_{vw}=\lambda_{vw}-d_{vw}=V(w)-V(v).
\tag{F5}
```

Aggregate at the state/decision level exactly as in (10.511):

```math
x_v
=\pi_v\mathbb E_v(\lambda-d)
=\sum_w\pi_w[V(w)-V(v)].
\tag{F6}
```

Here `\pi_w=\pi_vp_{vw}`.  Set

```math
z_v=(x_v)_+,
\qquad
u_v=(-x_v)_+.
```

At a leaf `\ell`, put

```math
e_\ell=\pi_\ell\varepsilon_\ell,
\qquad
\varepsilon_\ell=Q(A_\ell)-q_m,
\qquad
G=q_m-\alpha m^{3/2}.
```

Finally, suppose decrement-tolled endpoint buckets have actually been placed
at temporal nodes.  Let `R_v` be their total `\widehat a_b` capacity at `v`
and give them dilation `\kappa`.  The honest resource is

```math
s_v=u_v+\kappa R_v
```

at internal nodes, plus `e_\ell` at leaves.

Probability flow conservation makes (F5) telescope on every descendant
subtree:

```math
\begin{aligned}
\sum_{w\in T_v\text{ internal}}x_w
&=\sum_{\ell\in L_v}\pi_\ell V(\ell)-\pi_vV(v),\\
\sum_{\ell\in L_v}\pi_\ell V(\ell)
&=\pi_vG+\sum_{\ell\in L_v}e_\ell.
\end{aligned}
\tag{F7}
```

Therefore its exact Hall deficiency is

```math
\boxed{
H_v
=\pi_v\bigl(G-V(v)\bigr)-\kappa R(T_v).
}
\tag{F8}
```

Consequently the honest causal replacement for the globally relaxed credit
in (10.516) is

```math
\boxed{
E_{\rm caus}
=\max_{A\text{ antichain}}
\sum_{v\in A}
\left[
\pi_v(G-V(v))-\kappa R(T_v)
\right],
}
\tag{F9}
```

with the empty antichain allowed and with non-demand roots omitted as in the
proof of (F2).  Equivalently, zero error requires the weakest possible cut
family

```math
\boxed{
\pi_v\bigl(G-V(v)\bigr)
\le\kappa R(T_v)
\quad\text{for every reached }v.
}
\tag{F10}
```

At the root, `V=0` and `\pi=1`; this is just
`G\le\kappa R(T)`, the full-set scalar cut.  For `R=0`, every proper cut is

```math
Q(A_v)-\alpha r_v^{3/2}\ge q_m-\alpha m^{3/2}.
\tag{F11}
```

Thus proper causal cuts are conditional normalized comparisons at the
original root slope `\alpha`.  The telescope proves equality (F8), not the
inequality (F10).  Global minimality only says `Q(A_v)\ge q_{r_v}`; using it
in (F11) asks for

```math
q_{r_v}-\alpha r_v^{3/2}\ge q_m-\alpha m^{3/2},
```

which is another normalized comparison and is circular at the present
frontier.

Equation (F9) is nevertheless tail-useful.  If a matrix theorem supplies an
honest placement of residual buckets and proves

```math
E_{\rm caus}\le E_{n,m},
```

then the root cut and (10.515) give

```math
q_m-\left(\frac mn\right)^{3/2}q_n
\le4\kappa(s-1)n+E_{n,m}.
\tag{F12}
```

For `s=o(\sqrt n)`, the first error is `o(n^{3/2})`; convergence follows if
the combined normalized errors meet the adaptive/dyadic tail condition
(10.500).  The genuinely new matrix input is therefore exactly a uniform
bound on the maximum-antichain quantity (F9), not another scalar telescope.

### Outcome-level refinement

If every deletion outcome is kept as its own atom, subdivide the temporal
tree and attach

```math
x_{vw}=\pi_w[V(w)-V(v)]
```

to the outcome edge/node.  The same theorem applies verbatim.  A subtree
whose root includes the incoming outcome atom has signed deficiency before
residual service

```math
\pi_w\bigl(G-V(v)\bigr),
```

whereas the future subtree beginning *after* that atom has
`\pi_w(G-V(w))`.  These branchwise cuts are stronger than the parent-level
conditional-expectation cuts, consistently with the Jensen warning in
(10.511).

## 3. Exact realizable suffix obstruction

Use the exact order-seven minimizer already verified in (10.522):

```math
A_7=
\begin{pmatrix}
0&1&1&1&1&1&1\\
1&0&1&1&-1&-1&1\\
1&1&0&1&-1&1&-1\\
1&1&1&0&1&-1&-1\\
1&-1&-1&1&0&-1&-1\\
1&-1&1&-1&-1&0&-1\\
1&1&-1&-1&-1&-1&0
\end{pmatrix},
\qquad Q(A_7)=q_7=18.
```

This can be made a deterministic field-proportional path, not merely an
arbitrary deletion.  The vector `x=(1,-1,-1,-1,-1,-1,-1)` has
`x^{\mathsf T}A_7x=-18`.  Orient by `-A_7` and switch `x` to
one.  The row fields are

```math
(6,0,0,0,4,4,4).
```

Choose the eligible set `H` to contain only the first vertex.  Its
field-proportional deletion probability is `6/(7-1)=1`, giving the
child `B=A_7[-0]` with

```math
Q(B)=q_6=10.
```

For `B`, the vector `y=(1,-1,-1,-1,-1,-1)` has energy
`-10`.  Orient and switch in the same way; the row fields are

```math
(1,1,1,1,1,5).
```

Choose `H` to contain only the last vertex.  It is deleted with
probability `5/(6-1)=1`.  Each of the six possible order-five
children of `B` has `Q=q_5=8`; in particular this selected
child has zero terminal excess.  With

```math
\alpha=\frac{18}{7\sqrt7},
```

the first centered increment is

```math
x_1
=10-\frac{108\sqrt{42}}{49}
=-4.284081539\ldots<0.
\tag{F13}
```

The second centered increment, and hence the complete proper suffix at the
order-six state, is

```math
\boxed{
x_2
=\frac{108\sqrt{42}-90\sqrt{35}}{49}-2
=1.417812549\ldots>0.
}
\tag{F14}
```

Indeed, `\sqrt{42}>6` and `\sqrt{35}<6`, so the fraction before
subtracting two is greater than `108/49>2`.  Yet the full two-step
comparison is

```math
x_1+x_2
=G
=8-\frac{90\sqrt{35}}{49}
=-2.866268989\ldots<0.
\tag{F15}
```

Thus the root/full-horizon cut succeeds even with no residual buckets:
`\sqrt{35}>5` makes (F15) negative.  But the proper rooted subtree at
the order-six state consists of the positive demand `x_2`, with no
negative descendant credit and zero terminal excess.  Its exact uncovered
demand is therefore `x_2>0`.  The globally relaxed sink uses the
ancestral negative credit `-x_1` to cancel it;
descendant-resource compatibility cannot.  This is an exact
field-proportional induced-submatrix counterexample, not an abstract signed
sequence.

## 4. Reverse-time convention and prefix cuts

If instead a demand may use only resources at its ancestors (equivalently,
an already available credit may pay only descendant demands), the exact Hall
conditions are

```math
z(S)\le s(S)
\quad\text{for every ancestor-closed rooted set }S.
\tag{F16}
```

Indeed, the neighborhood of `J` is its ancestral hull; enlarging `J` to that
hull only increases demand.  On a chain these are prefix inequalities.
For a stopping frontier `\tau`, telescoping gives the exact prefix equality

```math
\sum_{v\text{ before }\tau}x_v
=\mathbb E V(A_\tau)-V(A_{\rm root}).
\tag{F17}
```

Again, no sign follows.

The order-five matrix in (10.518) is an exact reverse-time obstruction.
Its first singleton deletion leaves an order-four child of norm eight.  With
`\alpha=8/(5\sqrt5)`, that first prefix has positive demand

```math
8-\frac{64\sqrt5}{25}=2.275665977\ldots>0,
```

while continuing to order two gives global net

```math
2-\frac{16\sqrt{10}}{25}=-0.023857\ldots<0.
```

The future credit makes the global sink harmless but cannot satisfy the first
prefix cut in the ancestor-resource convention.

## 5. Same-order windows

On a chain, allowing arbitrary matching within a consecutive order window is
exactly contraction of that window to one supernode.  The theorem then says
that the only cuts are suffixes of windows, and (F8) is evaluated at each
window boundary.  This can remove microscopic timing walls inside a window,
but it leaves one conditional comparison per boundary.  Taking a single
window is precisely the globally relaxed sink.

On a branching tree the same conclusion holds if each connected
history/window cell is contracted and the quotient remains a tree.  If
resources may jump between incomparable histories merely because their
orders lie in the same band, neighborhoods are no longer laminar; the exact
answer reverts to the general Hall formula (F1) for that enlarged relation.
No collection of one-node subtree cuts is then sufficient in general.

## 6. Local residual certificate and what it does not yet prove

For an endpoint split

```math
U=D\sqcup X,
\qquad
A[U]=\begin{pmatrix}D&B\\B^{\mathsf T}&X\end{pmatrix},
```

and an oriented successor state `(\sigma,y)`, its layer capacity is

```math
c_b=\left[2\lVert By\rVert_1-
\bigl(Q(X)-\sigma y^{\mathsf T}Xy\bigr)\right]_+.
```

Since `\partial_b=Q(U)-Q(X)\ge0`, the residual has the exact form

```math
\boxed{
r_b=[c_b-\partial_b]_+
=\left[2\lVert By\rVert_1+\sigma y^{\mathsf T}Xy-Q(U)\right]_+.
}
\tag{F18}
```

If `r_b>0`, choose sibling signs `u` so that
`u_i(\sigma By)_i=|(By)_i|`; at zero coordinates choose either sign.  Then

```math
Q(U)\ge
\sigma\left(y^{\mathsf T}Xy+2u^{\mathsf T}By+u^{\mathsf T}Du\right)
```

and therefore

```math
\boxed{\sigma u^{\mathsf T}Du\le-r_b.}
\tag{F19}
```

For `r_b=0`, (F19) does **not** follow; the untruncated right-hand side in
the rearrangement is merely nonnegative.  Zero-residual buckets carry no
service and should simply be discarded.

Thus every positive residual bucket has a colocated negative diagonal-energy
certificate of at least the same magnitude.  This is promising causal
structure, but it is not yet a negative centered temporal atom
`[-\pi(\lambda-d)]_+`.  Different buckets may reuse the same diagonal block
with incompatible signs, and the certificate need not enter the temporal
telescope.  The weakest missing matrix theorem is a disjoint/fractional
realization of these certificates at compatible descendant histories strong
enough to imply every cut (F10), or directly a tail-summable bound on (F9).

## Bottom line

- Telescoping proves the exact conservation identities (F7), (F8), and
  (F17).
- It gives only the root/full-horizon scalar cut after the desired normalized
  comparison is inserted.
- Every proper suffix, prefix, branch, or contracted-window inequality is new
  conditional matrix information.
- The descendant-resource error is exactly a maximum-weight antichain, so a
  tail-useful sufficient theorem has a precise target: bound (F9) uniformly
  over every order pair used in the descent.
- The order-seven chain proves that replacing the global credit by honest
  temporal atoms is genuinely stronger, even for a deterministic induced
  chain starting at an exact global minimizer.
