# Wave 33 Route 1: exact global variation against a localized block

## Status and scope

All identities and inequalities below are proved from the exact global
minimality of `Q`; the finite examples are checked by
`global_variation_r33_check.py`.  The block-deletion inequality (R33.5) is
the one-block specialization of the earlier hybrid triangle inequality
(10.565).  Its application to the Wave 32 planted block, and the
near-active minimax law (R33.9), are the useful conclusions here.

The result rules out an **essential** localized `b^2` mode at cost
`O(b^(3/2))`.  It does not by itself put a low-row point in each favorable
conditional fiber: replacement witnesses may move to a different absolute
orientation, local label, and selector.

Throughout, `A` is a complete order-`n` signing which globally minimizes

```math
Q(A)=\max_{\sigma\in\{\pm1\},\,x\in\{\pm1\}^n}
\sigma x^{\mathsf T}Ax=q_n.
```

Fix `U subset [n]`, `|U|=b`, and put `T=U^c` and
`m=binom(b,2)`.  An oriented state is `omega=(sigma,x)`, modulo the global
replacement `x -> -x`.  Define

```math
d_e(\omega)=\sigma x_ix_j,
\qquad
s_e(\omega)=a_ed_e(\omega),
\qquad
\Delta_A(\omega)=q_n-\sigma x^{\mathsf T}Ax.
```

The orientation `sigma` is retained everywhere.  In particular it is not
legal to replace an oriented local state by its unoriented cut.

## 1. One-edge and simultaneous-edge certificates

If `F` is any edge set and `A^F` is obtained by flipping exactly `F`, then
direct expansion gives

```math
\boxed{
Q(A^F)-q_n
=\max_\omega\left\{-4\sum_{e\in F}s_e(\omega)
-\Delta_A(\omega)\right\}.
}
\tag{R33.1}
```

Global minimality is therefore exactly the family of conditions

```math
\boxed{
\forall F\subseteq E(K_n)\quad\exists\omega:\quad
\Delta_A(\omega)+4\sum_{e\in F}s_e(\omega)\le0.
}
\tag{R33.2}
```

For one edge this says `s_e=-1` at some state of slack at most four.
For all edges it is the covering-radius condition (10.733).  These
certificates allow the witness to depend on `F`; summing independently
chosen one-edge witnesses is invalid.

## 2. Exact planted-block response game

Write

```math
A=\begin{pmatrix}D&C\\ C^{\mathsf T}&E\end{pmatrix}
```

in the order `U,T`.  For a local oriented state
`u=(sigma,[x_U])`, define the response of all noninternal edges by

```math
Z_U(u)=\max_{y\in\{\pm1\}^T}
\sigma\{2x_U^{\mathsf T}Cy+y^{\mathsf T}Ey\}.
\tag{R33.3}
```

The local projective class has one representative with a fixed anchor bit;
there are `2^b` choices after retaining `sigma`.  If the internal block is
replaced by a signing `B`, then exactly

```math
\boxed{
Q(A^{U\to B})
=\max_u\left\{Z_U(u)+2\sum_{e\in E(U)}b_ed_e(u)\right\}.
}
\tag{R33.4}
```

This is the absolute-orientation version of the global replacement
identity.  Since `A` is a global minimizer, the minimum of the right side
over all signings `B` is exactly `q_n`.

Let `A^{0,U}` be the weighted matrix obtained by setting the internal
`U`-edges to zero.  Then `Q(A^{0,U})=max_u Z_U(u)`.  Replace `D` by any
order-`b` exact minimizer `G`.  Global minimality and the triangle bound in
(R33.4) give

```math
\boxed{
q_n-Q(A^{0,U})\le q_b.
}
\tag{R33.5}
```

This is exactly (10.565) for the partition consisting of `U` and singleton
outside blocks.  With (10.567),

```math
q_b\le 2\sqrt{b(b-1)(b+2)\log2}=O(b^{3/2}).
\tag{R33.6}
```

Thus the internal edges of **every** `b`-vertex block can be deleted while
losing at most `O(b^(3/2))` from the norm.  At `b asymp n^(3/4)` the loss is
only `O(n^(9/8))`, far below `b^2 asymp n^(3/2)`.

### Consequence for the Wave 32 planted clique

In (10.936), `A=P_U+R`, where `P_U` is a positive `b`-clique, `R` vanishes
on `U times U`, and `||R||_op<=C sqrt(n)`.  Here

```math
Q(A^{0,U})=Q(R)\le Cn^{3/2},
```

while reverse triangle gives

```math
Q(A)\ge b(b-1)-Cn^{3/2}.
```

If this `A` were an exact minimizer, (R33.5) would force

```math
\boxed{
b(b-1)\le2Cn^{3/2}+q_b.
}
\tag{R33.7}
```

For `b=ceil(Kn^(3/4))`, this requires `K^2<=2C+o(1)`.  Therefore the
large-`K` realization used to force (10.937) is positively detected by
global block variation and cannot be an exact minimizer.  More generally,
(R33.5) forbids any localized block whose deletion creates an
`omega(b^(3/2))` norm loss.

## 3. A shared near-active law from all block replacements

The scalar deletion inequality does not say how the replacement witnesses
relate.  Convexifying the internal signs gives a stronger common-law
certificate.

For `z in [-1,1]^m`, set

```math
\Phi(z)=\max_u\left\{Z_U(u)+2\sum_ez_ed_e(u)\right\},
\qquad
V_U=\min_{z\in[-1,1]^m}\Phi(z).
```

Round a minimizer `z` independently to signs `B_e` of means `z_e`.  For a
fixed local state, the centered change is a sum of `m` variables of range
four.  Hoeffding's lemma and log-sum-exp over the `2^b` oriented local
states give

```math
\mathbb E\max_u2\sum_e(B_e-z_e)d_e(u)
\le\eta_b,
\qquad
\eta_b:=\sqrt{8mb\log2}=2b\sqrt{(b-1)\log2}.
```

Some rounding has `Phi(B)<=V_U+eta_b`.  Every rounded block is a complete
sign replacement and has norm at least `q_n`; hence

```math
\boxed{q_n-\eta_b\le V_U\le q_n.}
\tag{R33.8}
```

Finite minimax gives

```math
V_U=\max_\mu\left\{
\mathbb E_\mu Z_U(u)-2\sum_e|m_e|
\right\},
\qquad m_e=\mathbb E_\mu d_e(u).
```

Choose a maximizing law and, for every local state in its support, choose
an outside completion attaining `Z_U`.  This lifts `mu` to actual full
oriented states.  Since every original score is at most `q_n`, exact
cancellation gives

```math
\boxed{
\begin{aligned}
\mathbb E_\mu\Delta_A&\le\eta_b,\\
4\sum_{e\in E(U)}(a_em_e)_+&\le\eta_b,\\
\mathbb E_\mu\left[2\sum_{e\in E(U)}a_ed_e\right]
&\le\eta_b/2,\\
\mathbb E_\mu Z_U&\ge q_n-\eta_b.
\end{aligned}
}
\tag{R33.9}
```

Indeed the mean original score is

```math
V_U+4\sum_e(a_em_e)_+\le q_n.
```

Unlike the separate one-edge witnesses, (R33.9) gives one shared law which
is near-active on average and whose internal edge correlations are almost
anti-aligned with the actual signing.

## 4. Quantitative exclusion of a uniformly localized parent layer

Put

```math
L_U(\omega)=2\sum_{e\in E(U)}a_ed_e(\omega).
```

Always `L_U>=-2m`.  If every full state with parent slack at most `t` has
`L_U>=L_0`, Markov's inequality applied to the law in (R33.9) yields

```math
\boxed{
L_0\left(1-\frac{\eta_b}{t}\right)
-2m\frac{\eta_b}{t}\le\frac{\eta_b}{2}.
}
\tag{R33.10}
```

This is a directly falsifiable necessary condition.

For example, suppose the induced block is a switched positive clique,
`a_ij=v_iv_j`.  Then

```math
L_U(\sigma,x)=\sigma\{(v\mathbin\cdot x_U)^2-b\}.
\tag{R33.11}
```

Fix `c>0`.  Taking `t=C_c eta_b` with a sufficiently large constant
`C_c` in (R33.10) proves the asymptotic statement:

> An exact minimizer cannot have `sigma=+1` and
> `|v dot x_U|>=cb` at every `O_c(b^(3/2))`-near-ground parent state.

Thus exact global variation rules out the parent-layer version of the
localized positive mode.  If the opposite orientation has a gap much
larger than `eta_b`, it cannot carry appreciable mass in (R33.9), and the
law necessarily finds a near-ground state with small block magnetization.

## 5. Sharp scope: signed variation is not a row-square theorem

There are two remaining losses.

1. `L_U` is signed.  Opposite-orientation high magnetization makes it very
   negative and can satisfy (R33.9) while its orientation-free row-square
   contribution remains large.
2. The law in (R33.9) concerns parent near-ground states.  The center-star
   target concerns favorable completions separately for many revealed
   selectors `S`.  Nothing in the replacement game forces the chosen parent
   witnesses to be favorable for a specified child.

The orientation loss is real already at the response-game level.  Take a
switched positive clique with `m` edges and retain only its two fully
magnetized oriented profiles `u_+,u_-`.  For an arbitrary threshold `q`, set

```math
Z(u_+)=q-2m,
\qquad
Z(u_-)=q+2m.
\tag{R33.12}
```

For any internal replacement `B`, writing `S=sum_eB_e`, the two scores are

```math
q-2m+2S,
\qquad
q+2m-2S.
```

The second is always at least `q`, while for the original positive clique
both profiles score exactly `q`.  Hence **every** block replacement passes,
yet every original active profile is fully magnetized and has maximal
orientation-free local quadratic magnitude.  This is an abstract response
profile, not a claimed complete-signing realization; it shows exactly what
cannot follow from the discrete block certificate without extra structure
on `Z_U`.

Actual minimizers exhibit witness migration in milder form.  In `A_9`, for

```math
U=\{0,2,4,5,7\},
```

every exact parent ground has signed internal energy `L_U=8`, but zeroing
those ten internal edges still gives `Q(A_9^{0,U})=24=q_9`.  For
`U={0,1,2,6}`, zeroing only six edges raises the norm from `24` to `32`.
Thus even exact grounds do not identify the state which supports the
zero-block or replacement response.

## 6. Resulting frontier for this route

Global minimality does supply a useful cut-compression theorem:

- no `b`-vertex block is norm-essential beyond `q_b=O(b^(3/2))`;
- all block replacements yield one `O(b^(3/2))`-near-active parent law with
  almost anti-aligned internal correlations;
- consequently the large-`K` planted-clique construction in (10.936)--
  (10.937) is not merely "not known to minimize": it violates the exact
  block-deletion inequality.

What is still needed for the center-star/forest route is a **conditional
witness-transfer lemma**.  It must either (a) constrain the opposite-
orientation boundary shield and transfer a state from (R33.9) into many
favorable child fibers, or (b) apply a coupled family of block variations
whose common witness already records the selectors.  Without such a bridge,
global signed variation does not control the orientation-free quadratic
boundary susceptibility or the row-square cost in each fiber.

