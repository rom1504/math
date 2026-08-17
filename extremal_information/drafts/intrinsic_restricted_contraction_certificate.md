# Intrinsic restricted contraction for max-plus quotient maps

**Status.** A sufficient certificate for the contraction hypothesis in
Theorem 17.2.  The exact width-two Ising check is in
[`../experiments/verify_width2_ising_intrinsic_contraction.py`](../experiments/verify_width2_ising_intrinsic_contraction.py).

## 1. Syndetic edge contraction

Let `G=(Q,E)` be a finite control graph for the legal switch language.  Give
each vertex `q` a metric fibre `Y_q`, and each edge `e:q->q'` a map

```math
T_e:Y_q\longrightarrow Y_(q')
```

with restricted Lipschitz coefficient at most `lambda_e<=1`.  A common
invariant set `Y_q=Y` is the literal setting of Theorem 17.2; the fibred form
is the same statement after adjoining the finite control state.

### Theorem 1 (syndetic restricted contraction)

Suppose a set `C subset E` and `theta<1` satisfy `lambda_e<=theta` for
`e in C`.  Delete the edges in `C`.  If the remaining graph is acyclic and
its longest directed path has `H` edges, then every legal composition of
length

```math
L=H+1
```

has restricted Lipschitz coefficient at most `theta`.

Indeed, every length-`H+1` path contains an edge of `C`; coefficients before
and after it are at most one, and Lipschitz coefficients multiply.  The
bounded-gap condition is exactly checkable: deleting `C` leaves an acyclic
graph if and only if legal paths avoiding `C` have bounded length.  A
topological longest-path pass costs `O(|Q|+|E|)` and does not inspect legal
`L`-products.

There is also a fractional scrambling form.  If positive vertex weights
`a_q` and `mu<1` obey

```math
lambda_e {a_(q')\over a_q}<=mu
\qquad(e:q->q'),                                             \tag{1.1}
```

then every length-`L` path has coefficient at most

```math
Cmu^L,\qquad C={max_q a_q\over min_q a_q}.                    \tag{1.2}
```

Thus any `L` with `Cmu^L<1` verifies Theorem 17.2.  In logarithmic variables,
(1.1) is a system of difference constraints.  It is feasible exactly when
every positive-coefficient directed cycle has geometric-mean coefficient at
most `mu`, so it can be checked by a maximum-cycle-mean algorithm rather
than by word enumeration.

## 2. Max-plus dominance cones give intrinsic reset edges

For the column-convention max-plus map

```math
(F_Sx)_i=max_j(S_(ij)+x_j),
```

fix an input column `k` and put

```math
D_(S,k)=\left\{x:
x_k-x_j>=max_i(S_(ij)-S_(ik))\quad\hbox{for every }j\right\}. \tag{2.1}
```

If `Y subset D_(S,k)`, column `k` wins every output maximum throughout `Y`,
possibly with ties, and

```math
F_Sx=S_(*k)+x_k\mathbf 1.                                    \tag{2.2}
```

Hence `F_S` is projectively constant on `Y` and has restricted coefficient
zero.  When `Y` is a polytope, (2.1) is checked by one linear optimization
per input column.  Combined with Theorem 1, reset-cone containment,
generator-wise invariance, and one graph acyclicity test certify the
`(L,rho)=(H+1,0)` hypothesis of Theorem 17.2.

This reset conclusion is not an artifact of using a strong sufficient
condition.  Let `Y` be convex with nonempty interior in projective space.
On a tie-free max-plus selector cell the derivative is a deterministic
row-stochastic selector matrix.  Its Dobrushin coefficient is zero when all
rows select the same input, and one otherwise.  Therefore

```math
operatorname {Lip}_H(F_S|Y)<1
\quad\Longrightarrow\quad
operatorname {Lip}_H(F_S|Y)=0.                                 \tag{2.3}
```

To prove (2.3), a nonconsensus selector cell gives a small segment inside
`Y` on which the Hilbert ratio is exactly one.  If no such cell exists, the
projective derivative vanishes on every full-dimensional cell; continuity
and convexity make the projective map constant.  Since a composition of
max-plus maps is again max-plus, the same dichotomy applies to every legal
product on a full-dimensional convex `Y`.

Consequently a common Dobrushin minorization can produce a genuinely
fractional coefficient only on a thin or otherwise restricted carrier.  On
a full-dimensional interval or polytope it either proves a reset or fails
with coefficient one.

## 3. Exact width-two Ising certificate

For two boundary spins, write

```math
(T_(J,h,V)v)(y)=h_1y_1+h_2y_2+Vy_1y_2
+max_x\{v(x)+J_1x_1y_1+J_2x_2y_2\}.                            \tag{3.1}
```

Normalize at `++` and use gaps `g_x=v(++)-v(x)`.  Let

```math
Y=\{4<=g_(-+)<=20,\ 4<=g_(+-)<=20,\ 8<=g_(--)<=40\}.           \tag{3.2}
```

This is a full three-dimensional projective box contained in the Hilbert
ball of radius `R=20`.  Consider

```math
A:T_((4,4),(6,6),0),\qquad
W:T_((1,1),(3,3),0).                                            \tag{3.3}
```

The elementary inequality

```math
min_x(a_x-b_x)<=max_xa_x-max_xb_x<=max_x(a_x-b_x)               \tag{3.4}
```

shows that `A` maps the whole projective space into `Y`: its two one-spin
gaps lie in `12+[-8,8]=[4,20]`, and its two-spin gap lies in
`24+[-16,16]=[8,40]`.

For `W`, the worst possible kernel advantages over predecessor `++` are
`2,2,4` for predecessors `-+,+-,--`.  The lower gaps in (3.2) are `4,4,8`,
so (2.1) holds strictly.  Thus `W` is a projective reset on all of `Y`; its
image has gaps `(8,8,16)` and lies back in `Y`.

The other column really is noncontractive.  For

```math
v_0=(-8,-4,-4,0),\qquad v_1=(-8,-5,-4,0),                       \tag{3.5}
```

all four selectors of `A` are distinct and

```math
\|Av_0-Av_1\|_H=\|v_0-v_1\|_H.                                \tag{3.6}
```

Let the legal control graph alternate `A` and `W`.  After deleting the reset
edge `W`, the longest path has one edge.  Theorem 1 therefore gives

```math
L=2,\qquad rho=0.                                               \tag{3.7}
```

The contraction part of Theorem 17.2 is now verified by two generator
checks and one acyclicity check.  Its simulator bound specializes to
`(1+40/h)^3` states and aggregate error `h+2(epsilon+h)`.

## 4. Rejected substitutes

1. A periodically occurring weak Ising column with small **image diameter**
   does not verify Theorem 17.2.  Every nonzero unrestricted Ising bond still
   has Hilbert coefficient one.  It supports the separate additive
   small-image reset estimate, not multiplicative contraction.  It becomes
   valid here only because the invariant gaps place the whole restricted
   set in one saturation cone.
2. A one-orbit selector/fan argument is insufficient.  Contraction compares
   two states, and different cells can follow different selectors.  Unless
   one reset cone covers the entire declared fibre, a sound verifier needs a
   paired-cell reachability graph; that is precisely the enumeration burden
   this certificate is meant to avoid.
3. Dynamic programming over all length-`L` products compresses the
   enumeration but merely restates the original hypothesis.  The accepted
   certificate instead checks generator-wise dominance/invariance and the
   absence of a reset-free control cycle.
