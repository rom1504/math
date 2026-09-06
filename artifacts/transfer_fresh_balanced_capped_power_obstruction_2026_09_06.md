# Balanced powered aggregation fails even after a sub-half truncation

Date: 2026-09-06. Fresh convergence track. Status: proved below; the director
independently reconstructed the bounded telescoping and extraction argument
at checkpoint five. This is a pointwise theorem for actual signings, not
a theorem about the original minimum values.

Let `Q(A)=max_x |x^T A x|/2`, `c(A)=Q(A)/|A|^(3/2)`, and
`ell=liminf M_n/n^(3/2)<1/2`. Fix ANY `b` with `ell<b<1/2`.
For an order-`m` signing define the truncated power

```math
T_b(A)=\min\{Q(A),b m^{3/2}\}^{2/3}
      =m\min\{c(A),b\}^{2/3}.
```

## Theorem

There are `kappa_b>0`, arbitrarily large even orders `n=2m`, and
actual hollow symmetric signings `P` with an equal-size partition
into principal children `B,D`, such that

```math
c(P)<b<1/2,
\qquad T_b(B)+T_b(D)-Q(P)^{2/3}\ge\kappa_b n.              (1)
```

Thus truncating each actual child's cap at a strict-subhalf coefficient
does not rescue universal pointwise powered aggregation even for equal
halves. It is NOT asserted that the untruncated child caps are below
one half. Nor does (1) remain proved if the child expressions are
replaced by `M_m^(2/3)`.

## Proof

The uniform small-fixed-retention theorem in Section 2 of
`transfer_director_exponential_selector_cost_2026_09_06.md` gives a
`rho>0` such that every starting parent of normalized cap at most
`1/2`, restricted uniformly to any fixed fraction at most `rho`, has
normalized cap above `1/2` with probability tending to one.

Choose a FIXED integer `r` with `2^(-r)<rho`. Take a liminf-realizing
minimizing sequence of orders `D_j`, and delete fewer than `2^r`
vertices to obtain orders `N_j` divisible by `2^r`. Principal cap
monotonicity and the definition of `ell` show the resulting actual
starting caps still converge to `ell` after normalization.

Take `r` nested uniform restrictions to exact halves, writing their
normalized caps as `C_0,...,C_r`, and put

```math
V_i=\min\{C_i,b\}^{2/3},\qquad
\Delta=b^{2/3}-\ell^{2/3}>0.
```

The terminal selector is marginally uniform of size `2^(-r)N_j`.
Thus `V_r` tends in probability, and by boundedness in mean, to
`b^(2/3)`, while `V_0->ell^(2/3)`. The exact telescoping identity gives

```math
\sum_{i=0}^{r-1}\mathbb E(V_{i+1}-V_i)=\Delta-o_j(1).
```

For every sufficiently large starting order some index therefore has
mean increment at least `Delta/(2r)`. Conditional on its deterministic
parent signing `P`, the next half is uniform. If `c(P)>=b`, its
conditional increment is nonpositive, because `V_i=b^(2/3)` is the
largest possible value. Consequently at least one deterministic parent
realization with `c(P)<b` has

```math
\mathbb E\left[\min\{c(B),b\}^{2/3}\mid P\right]
           -c(P)^{2/3}\ge\Delta/(2r).                    (2)
```

For a uniform half `B`, the complementary half `D` has the IDENTICAL
marginal law; no independence is required. Since `|B|=|D|=|P|/2`,

```math
\mathbb E\left[\frac{T_b(B)+T_b(D)}{|P|}\mid P\right]
=\mathbb E\left[\min\{c(B),b\}^{2/3}\mid P\right].
```

Choose one partition at least as large as this expectation. Equation
(2) proves (1), for example with `kappa_b=Delta/(4r)` to leave harmless
slack. The selected parent's order is at least `2^(-r)N_j` and tends
to infinity. All limits use fixed `b,r`; there is no growing-depth
argument.

## What survives this falsifier

The minimum-value reverse-Fekete target in
`transfer_fresh_reverse_fekete_2026_09_06.md` survives untouched.
The theorem instead identifies a precise limitation of trying to prove
it by substituting arbitrary principal child caps, even clipped at a
coefficient `b<1/2`. The asymmetric construction in
`transfer_fresh_subhalf_asymmetric_power_obstruction_2026_09_06.md`
additionally makes both actual children strictly sub-half, and its
Section 7 can make the small child an exact minimizer, but those
conclusions currently require highly asymmetric splits.

No new primary theorem is imported here beyond the previously audited
sparse-selector result. The proof uses only marginal uniformity,
bounded telescoping, conditional expectation, and complement symmetry.
