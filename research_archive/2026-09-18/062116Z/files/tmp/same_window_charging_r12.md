# Same-window centered charging: residual harvest, exact Hall dual, and finite wall

Status: every identity and inequality below is proved exactly.  The transport
statement is **conditional**: no natural compatibility relation between a
temporal peeling tree and an endpoint tree is asserted.  The order-five
example proves that an unconditional pointwise, or even total, finite Hall
statement is false.  The surviving target is necessarily grouped and
asymptotic.

## 1. Exact temporal centering

Let a finite peeling process start from an exact order-`n` minimizer `A`, with
`Q(A)=q_n`, and stop at deterministic order `m>=rho n`.  A temporal node `v`
is a history, has reached probability `pi_v`, current order `r_v`, and a
conditional deletion outcome `H` leading to a child of order `r_{v,H}`.  Put

```math
d_{v,H}=Q(A_v)-Q(A_{v,H})\ge0,
\qquad
\lambda_{v,H}
=\frac{q_n}{n^{3/2}}
\left(r_v^{3/2}-r_{v,H}^{3/2}\right).
```

Let `epsilon_L=Q(A_L)-q_m` at a leaf.  Since both `d` and `lambda`
telescope on every realized path,

```math
\boxed{
q_m-\left(\frac mn\right)^{3/2}q_n
=\sum_v \pi_v\,\mathbb E_v(\lambda_{v,H}-d_{v,H})
-\mathbb E\epsilon_L.
}
\tag{R12.1}
```

This is the safe centered version of (10.499).  In the notation of (10.329),

```math
\boxed{
\pi_v\mathbb E_v(\lambda-d)
=\pi_v\left[\mathbb E_v\lambda+\mathbb E_v g
-(2a_v-c_v)\right].
}
\tag{R12.2}
```

The expectation is essential: generally `2a_v-c_v=E_v(d+g)`, not
`d_H+g_H` outcome by outcome.

Define the coarsest (and therefore smallest) positive demand and negative
credit at a parent by

```math
z_v=\left[\pi_v\mathbb E_v(\lambda-d)\right]_+,
\qquad
r_v=\left[-\pi_v\mathbb E_v(\lambda-d)\right]_+,
\qquad
C=\sum_v r_v+\mathbb E\epsilon_L.
```

Then (R12.1) is exactly

```math
\boxed{
q_m-\left(\frac mn\right)^{3/2}q_n
=\sum_v z_v-C.
}
\tag{R12.3}
```

If compatibility must remember individual deletion outcomes, replace the
parent atoms by

```math
z_{v,H}=\pi_v\Pr_v(H)[\lambda_{v,H}-d_{v,H}]_+,
\qquad
r_{v,H}=\pi_v\Pr_v(H)[d_{v,H}-\lambda_{v,H}]_+.
```

Their signed difference is still (R12.1), while Jensen shows that their total
positive mass is at least the parent-level mass.  Thus the parent version is
the strongest scalar reduction, and the transition version is the stronger
geometric requirement.

## 2. A new exact residualized harvesting theorem

Fix one endpoint tree and one of the allocations from Section 10.62.  For an
oriented bucket `b` on the edge

```math
e(b):\quad U_b=D_b\sqcup X_b\longrightarrow X_b,
```

write `c_b` for its capacity, `a_b` for its allocation, and

```math
\partial_b=Q(U_b)-Q(X_b).
```

Define the **decrement-tolled residual allocation**

```math
\boxed{
\widehat a_b
=\frac{a_b}{c_b}[c_b-\partial_b]_+,
}
\tag{R12.4}
```

with value zero when `c_b=0`.  Then for arbitrary nonnegative weights
`omega_b`,

```math
\boxed{
\sum_b\omega_b\widehat a_b
\le4\sup_\pi\sum_{b\in\pi}
\omega_b[c_b-\partial_b]_+
\le4\sup_\pi\sum_{b\in\pi}\omega_b Q(D_b).
}
\tag{R12.5}
```

Proof.  The path-cover measure `nu` from (10.489) has mass at most four and
`nu{pi:b in pi}>=a_b/c_b`.  Hence Tonelli gives the first inequality.  The
local bound behind (10.490) is

```math
c_b\le \mathcal L_{e(b)}
\le \partial_b+Q(D_b),
```

which gives the second.  An oriented chain contains only one oriented bucket
on each visited edge, so no factor two is lost here.

In particular, retain only edges in the order window `|U_b|>=rho n` whose
peeled siblings have `|D_b|<=s`.  Siblings on a chain are disjoint, and
`Q(D)<=|D|(|D|-1)`, so

```math
\boxed{
\sum_b\widehat a_b\le4(s-1)n,
\qquad
\sum_b\frac{\widehat a_b}{|U_b|^{3/2}}
\le\frac{4(s-1)n}{(\rho n)^{3/2}}.
}
\tag{R12.6}
```

This is the useful improvement over (10.491): after paying the full endpoint
decrement as a toll, the `Q(A)` term disappears.  Thus `s=o(sqrt(n))` makes
the residual harvest `o(n^{3/2})`.

## 3. The exact finite primal/dual

The preceding theorem supplies small sink capacities, but does not supply a
map from temporal demands to them.  Let `I` be either the parent atoms or the
finer transition atoms above.  Let `B` be the retained residual endpoint
buckets.  A proposed proof must define a genuine compatibility relation

```math
R\subseteq I\times B.
```

For example, it might require compatible vertex sets, orders in the same
window, and nested exposing states.  Nothing proved so far produces such an
`R`; order proximity alone is not compatibility.

For a fixed dilation `kappa>=0`, the strongest finite charging LP is

```math
\begin{aligned}
E_\kappa=\min\quad&\sum_{i\in I}e_i\\
\text{subject to}\quad&
\xi_i+\sum_{b:\,iRb}\gamma_{ib}+e_i=z_i &&(i\in I),\\
&\sum_i\xi_i\le C,\\
&\sum_i\gamma_{ib}\le\kappa\widehat a_b &&(b\in B),\\
&\xi_i,\gamma_{ib},e_i\ge0.
\end{aligned}
\tag{R12.7}
```

Here the one universal credit sink of capacity `C` represents all negative
centered increments and the terminal excess.  Ordinary max-flow/min-cut gives
the exact dual/cut formula

```math
\boxed{
E_\kappa
=\max_{J\subseteq I}
\left[z(J)-C-\kappa\widehat a(N(J))\right]_+,
}
\tag{R12.8}
```

where `N(J)={b: iRb for some i in J}`.  Thus all demand is chargeable if and
only if every Hall cut satisfies

```math
z(J)\le C+\kappa\widehat a(N(J)).
\tag{R12.9}
```

When `C=0` and `kappa=1`, the maximum common served fraction is

```math
\boxed{
\eta_*
=\min\left\{1,
\min_{J:\,z(J)>0}\frac{\widehat a(N(J))}{z(J)}
\right\}.
}
\tag{R12.10}
```

More generally replace the numerator by `C+widehat a(N(J))`.  The cap at one
is needed because a capacity surplus cannot serve more than the full demand.

Combining (R12.3), (R12.6), and (R12.8) proves the conditional coefficient
theorem

```math
\boxed{
q_m-\left(\frac mn\right)^{3/2}q_n
\le4\kappa(s-1)n+E_\kappa.
}
\tag{R12.11}
```

Therefore **any bounded dilation coefficient** `kappa`, not a coefficient
below `1/3`, is sufficient after decrement-tolling, provided
`s=o(sqrt(n))` and `E_kappa=o(n^{3/2})`.  To bridge an arbitrary liminf
subsequence, the normalized sum of the right-hand errors must satisfy the
tail condition (10.500).  For example, bounded `kappa` and
`s(n)<=n^{1/2-delta}` give a dyadically summable residual error.

This is a conditional reduction, not a proof of normalized descent: the
unknown content is precisely (R12.9) with a matrix-derived compatibility
relation and tail-uniform error.

## 4. Exact realizable finite obstruction

The order-five minimizer (10.350) gives a particularly clean wall.  In doubled
normalization take

```math
A=
\begin{pmatrix}
0&-1&1&-1&1\\
-1&0&-1&1&1\\
1&-1&0&1&1\\
-1&1&1&0&1\\
1&1&1&1&0
\end{pmatrix},
\qquad Q(A)=q_5=8.
```

The positive ground `p=(1,1,1,1,1)` has gauged row fields
`(0,0,2,2,4)`.  The negative ground
`n=(1,1,1,1,-1)` gives the endpoint split

```math
D=\{4\},\qquad X=\{0,1,2,3\}.
```

Here `Q(X)=q_4=8`, `Q(D)=0`, and both child imbalances vanish.  The two
oriented capacities toward `X` are zero.  The two capacities toward `D` are
both eight, but their endpoint decrement is
`Q(A)-Q(D)=8`; hence every residual (R12.4) is zero.  The Section 10.62
allocation puts eight units in each of those two saturated buckets, so its
entire decrement-tolled residual harvest is exactly zero.

Now use the field-proportional rule with heavy set `{4}`.  Since
`r_4/(5-1)=1`, vertex four is deleted deterministically.  The retained core is
the exact order-four minimizer, and

```math
d=0,\qquad g=8,\qquad a=4,\qquad c=0,\qquad\epsilon_L=0.
```

Thus the centered demand is

```math
\boxed{
z=8\left[1-\left(\frac45\right)^{3/2}\right]
=8-\frac{64\sqrt5}{25}>0,
}
\tag{R12.12}
```

while `C=0` and `widehat a(B)=0`.  This is the same positive ground, the same
`1+4` cut, a global minimizer, an optimal child, and an exact
field-proportional deletion.  Consequently `eta_*=0`, every finite exact
dilation fails, and even the total Hall cut fails.

The order-nine singleton (10.332) supplies the complementary local warning:
`d=0`, `g=16`, and the predecessor restriction has
`2||By||_1=16`, but its deficit is also sixteen, so its all-successor payoff
and residual are zero.  Thus raw or centered replenishment cannot be mapped
pointwise to the local successor bucket.

These are finite walls.  They do not rule out a large-order grouped theorem
with an additive, tail-summable error; no asymptotic blow-up preserving global
minimality is known.

## 5. Weakest missing hypothesis

For a fixed compatibility relation, (R12.8) is logically exact.  Hence the
weakest hypothesis that actually constructs a local charging map is the
Hall-deficiency estimate

```math
\max_{J\subseteq I}
\left[z(J)-C-\kappa\widehat a(N(J))\right]_+
\le E_{n,m},
\tag{R12.13}
```

for some bounded `kappa`, with the combined error
`4 kappa(s-1)n+E_{n,m}` satisfying (10.500).  For the scalar normalized
comparison alone, only the full-set instance `J=I` is needed; all subsets are
exactly what is additionally required for an honest local transport.

The order-five wall proves that (R12.13) cannot be pointwise or error-free.
The missing mathematical input must therefore be a genuinely macroscopic
exchange/thickness statement which, after grouping many temporal events in
one order window, either supplies compatible positive residual buckets or
forces compensating negative centered increments/terminal excess.  This is
strictly sharper than asking for another raw path-congestion bound.

## 6. Audit

`tmp/check_residual_small_r12.py` independently enumerates all projective
positive/negative endpoints of (10.350) and (10.445), recomputes every root
capacity, endpoint decrement, and residual, and verifies the order-five data
above.  It uses only exact integer arithmetic.
