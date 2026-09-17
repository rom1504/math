# Wave 22 Route 3: exact flip covers and the heavy-field migration wall

This note tests whether exact minimality under arbitrary simultaneous edge
flips can upgrade the pair-field bound (10.711) to the square-root cavity
reward required by (10.617).  The answer for this mechanism is negative:
all multi-edge certificates have an exact covering-code form, and their
Gibbs transforms bound coordinate affinity in the wrong direction.  This is
a structural wall, not an asymptotic counterexample to (10.617).

All algebraic statements and the finite `A_8/A_9` audits below are verified
by `heavy_field_migration_r22.py`.

## 1. Every simultaneous edge flip is one fixed-radius cover

Let `B` be an order-`r` complete signing which attains the minimum `q=q_r`,
let `E=E(K_r)`, and put `N=|E|`.  On the oriented projective state space set

```math
s_e(\omega)=\sigma b_{ij}x_ix_j,
\qquad
e_B(\omega)=2\sum_{e\in E}s_e(\omega),
\qquad
\Delta_\omega=q-e_B(\omega).
```

Write

```math
C_\omega=\{e:s_e(\omega)=-1\},
\qquad
R=\frac N2-\frac q4.
```

The complete-energy identity is exactly

```math
\boxed{
|C_\omega|=R+\frac{\Delta_\omega}{4}.
}
\tag{R22.1}
```

For an arbitrary edge set `S`, let `B^S` be obtained by flipping precisely
the signs in `S`.  The energy of the same state is

```math
e_{B^S}(\omega)
=e_B(\omega)-4\sum_{e\in S}s_e(\omega)
=2N-4|S\mathbin\triangle C_\omega|.
```

Since every `B^S` is a signing and `B` is a global minimizer,
`Q(B^S)\ge q`.  Consequently, for every `S\subseteq E` there is an
oriented state satisfying

```math
\boxed{
|S\mathbin\triangle C_\omega|\le R.
}
\tag{R22.2}
```

Conversely, (R22.2) is precisely the simultaneous-flip certificate.  Thus
the `2^r` sets `C_\omega`, which form a translate of the augmented cut code,
have radius-`R` Hamming balls covering the whole `N`-dimensional edge cube.
Equivalently, an exact minimizer is a deepest coset of that code and `R` is
its covering radius.  In particular, passing from one-edge witnesses to
*all* edge-set witnesses does not give a new kind of local certificate: it
recovers the global covering-radius definition exactly.

## 2. Restricted star and block covers

Fix any block of edges `H\subseteq E`, of size `m`.  For a state define

```math
d_H(\omega)=|C_\omega\cap H|,
\qquad
z_H(\omega)=\sum_{e\in H}s_e(\omega)=m-2d_H(\omega),
```

and

```math
\rho_H(\omega)
=R-|C_\omega\setminus H|
=d_H(\omega)-\frac{\Delta_\omega}{4}.
```

Restricting (R22.2) to targets `S\subseteq H` says exactly that the balls

```math
\boxed{
\left\{S\subseteq H:
|S\mathbin\triangle(C_\omega\cap H)|\le\rho_H(\omega)
\right\}
}
\tag{R22.3}
```

cover the whole `m`-cube.  The Gibbs exponent attached to the same restricted
ball radius is

```math
\boxed{
\Delta_\omega+2z_H(\omega)=2m-4\rho_H(\omega).
}
\tag{R22.4}
```

For a star `H=\delta(i)`, `m=r-1`, `d_H=d_i^-`, `z_H=h_i`, and the Gibbs
moment in (R22.4) is exactly the coordinate Hellinger affinity.  Hence the
full row-replacement certificate is a cover of the row cube by balls of
radius

```math
\rho_i(\omega)=d_i^-(\omega)-\Delta_\omega/4.
```

This is stronger and cleaner than the one-edge incidence description, but
it still has the unfavorable direction: states with larger restricted
radius have larger geometric mass in (R22.4).

## 3. Averaging every block competitor gives a lower moment bound

Let

```math
D_\beta(B)=\sum_\omega e^{-\beta\Delta_\omega}
```

and, retaining the original baseline `q`, define

```math
D_\beta(B^S;q)
=\sum_\omega
\exp\left[-\beta\left(
\Delta_\omega+4\sum_{e\in S}s_e(\omega)
\right)\right].
```

Because `Q(B^S)\ge q`, at least one summand is at least one, so
`D_\beta(B^S;q)\ge1`.  Averaging exactly over every `S\subseteq H` gives

```math
\boxed{
2^{-m}\sum_{S\subseteq H}D_\beta(B^S;q)
=D_\beta(B)\cosh(2\beta)^m
\mathbb E_{\nu_{\beta,B}}e^{-2\beta z_H}.
}
\tag{R22.5}
```

Indeed, state by state,

```math
2^{-m}\sum_{S\subseteq H}
e^{-4\beta\sum_{e\in S}s_e}
=\prod_{e\in H}\frac{1+e^{-4\beta s_e}}2
=e^{-2\beta z_H}\cosh(2\beta)^m.
```

Minimality therefore implies only

```math
\boxed{
\mathbb E_{\nu_{\beta,B}}e^{-2\beta z_H}
\ge
\frac{1}{D_\beta(B)\cosh(2\beta)^m}.
}
\tag{R22.6}
```

For `H=\delta(i)`, this is

```math
\boxed{
\operatorname{BC}_i
\ge
\frac{1}{D_\beta(B)\cosh(2\beta)^{r-1}},
}
\tag{R22.7}
```

which is the opposite of the exponentially small upper bound in (10.624).
For a cut `H=\delta(U)`, (R22.5) is the analogous multi-coordinate
Hellinger identity.  For an arbitrary block it is the corresponding block
field moment identity.  Thus stars, cuts, internal blocks, and arbitrary
edge blocks all retain the same wrong direction after exact Gibbs
aggregation.

When `H` is a cut, complementary perturbations are switching-equivalent and
`D_\beta(B^S;q)=D_\beta(B^{H\setminus S};q)`.  This symmetry does not reverse
(R22.6).

## 4. Why the certificates permit heavy-field migration

The covering formulation isolates both losses in the proposed argument.

First, a perturbation `S` may choose a different certifying center
`C_\omega`.  There is no congestion statement forcing one coordinate or one
Gibbs-visible state to serve many perturbations.  Second, centers far from
the original identity have large deficit by (R22.1) and are exponentially
suppressed in the original Gibbs law, but their radius-`R` balls remain fully
available to certify distant perturbations.  The existential cover does not
weight a certificate by its original Gibbs mass.

For a small perturbation there is a useful but insufficient localization:
if `|S|=k` and `\omega` certifies it, then

```math
\Delta_\omega
\le -4\sum_{e\in S}s_e(\omega)
\le4k.
\tag{R22.8}
```

So small flips cannot use arbitrarily high-deficit witnesses.  But summing
all such localized witnesses again produces a covering lower bound such as
(R22.6), not an upper bound on affinity.  A proof of (10.617) would need a
new *weighted congestion* or common-coordinate theorem: for some fixed
coordinate, most of the original Gibbs denominator would have to carry a
field of order `sqrt(r)`.  Neither the fixed-radius cover nor pairwise field
means forces that.  They allow the order-`r` field to be rare for each
coordinate and to migrate with the state, exactly the obstruction after
(10.711).

## 5. Exact `A_8/A_9` audit

The checker enumerates all oriented states and every one of the `2^{r-1}`
row replacements at every coordinate.  It verifies (R22.1)--(R22.8), the
complement symmetry for stars, and (R22.5) on stars, triangles, four-cliques,
and matchings.

- For `A_8`, `(N,q,R)=(28,20,9)`.  Every star replacement is certified with
  deficit at most four.  The eight exact grounds all have sorted field
  profile `(5,3,3,3,3,1,1,1)`, and across the grounds every coordinate
  occurs as a maximal-field coordinate.
- For `A_9`, `(N,q,R)=(36,24,12)`.  Remarkably, every star replacement is
  already certified by an exact ground.  Across its 25 oriented grounds,
  every coordinate occurs as a maximal-field coordinate; maximal fields
  range through `4,6,8` depending on the ground.  Thus actual complete
  minimizers exhibit the finite witness/field migration permitted by the
  abstract cover.
- At `beta=0.5`, the exact star affinities range from `0.323374` to
  `0.323374` on `A_8` and from `0.263762` to `0.430737` on `A_9`.  At
  `beta=1`, they range from `0.104231` to `0.104231` and from `0.172025` to
  `0.334714`, respectively.  The averaged competitor sums in (R22.5) are
  much larger than one, illustrating why the lower bound (R22.7) is slack.

The earlier exact polynomial certificate `\kappa_{\beta,i}(A_9)<3` for all
coordinates and all `\beta>0` remains compatible with every simultaneous
flip certificate above.  This is a genuine finite complete-signing wall,
but it does not supply a growing counterexample to (10.617).

## 6. Scoped conclusion

- **Verified theorem:** arbitrary multi-edge minimality is exactly the
  radius-`R` augmented-cut-code cover (R22.1)--(R22.3).
- **Verified Gibbs consequence:** averaging all block competitors gives the
  exact transform (R22.5) and only the lower moment/affinity bound
  (R22.6)--(R22.7).
- **Closed mechanism:** unweighted simultaneous-flip certificates, even for
  all stars, cuts, and edge blocks, cannot by themselves exclude rare heavy
  fields or migration; their canonical aggregation has the wrong sign.
- **Still open:** (10.617).  No asymptotic family of exact complete
  minimizers with all `\kappa_i=o(sqrt(r))` is constructed.  A successful
  continuation needs a minimizer-specific Gibbs-weighted congestion/common-
  coordinate statement not contained in covering-radius minimality.

