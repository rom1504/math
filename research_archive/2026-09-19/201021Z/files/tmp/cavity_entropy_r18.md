# Wave 18 Route 1: cavity reward is a coordinate Hellinger barrier

## Status

- **Verified:** an exact weighted cube-edge / Bhattacharyya representation of
  every soft cavity reward, including the factor of two in the ledger's
  oriented-projective convention.
- **Verified:** the extension lemma holds at every deficit layer, not only for
  child grounds.
- **Verified, scoped no-go:** the ordinary tensorized cube log-Sobolev bound
  applied to this representation yields only constant-scale reward
  information.  Antipodal Hamming balls show that cube geometry,
  antipodality, and support entropy alone cannot force an order-`sqrt(r)`
  reward.  This does **not** rule out an inequality using the quadratic Gibbs
  form or global-minimizer ancestry.
- **Numerical/exact-enumeration audit:** the identities pass on `A8`, `A9`,
  and all their principal restrictions at three temperatures.  The exact
  `A9` Bernstein certificate becomes precisely a lower bound on every
  coordinate affinity.
- **Open:** prove a near-isolation theorem for Gibbs measures of quadratic
  sign forms inherited from global minimizers.  Generic projection entropy
  does not provide it.

## 1. Exact deficit-edge identity

Fix an order-`r` matrix `B`, a vertex `i`, and `C_i=B[-i]`.  Write

```math
d_i=Q(B)-Q(C_i),
\qquad
\delta_B(\sigma,x)=Q(B)-\sigma x^{\mathsf T}Bx,
```

and use the ledger's oriented-projective deficit sum

```math
D_\beta(B)=\sum_{\omega\in\Omega_r}e^{-\beta\delta_B(\omega)}.
```

It is easiest to derive the identity on the twofold full-cube lift
`(sigma,x) in {+1,-1} x {+1,-1}^r`; both parent and child deficit sums are
then doubled, so every ratio below is unchanged.  If `y=x_{-i}` has child
deficit `delta` and

```math
h=\sigma b_i^{\mathsf T}y,
```

the two extensions have parent deficits

```math
\delta_+=d_i+\delta-2h,
\qquad
\delta_-=d_i+\delta+2h.
\tag{CE1}
```

Parent optimality implies `2|h|<=d_i+delta`.  More importantly, the geometric
mean of the two parent weights is exact:

```math
\sqrt{e^{-\beta\delta_+}e^{-\beta\delta_-}}
=e^{-\beta(d_i+\delta)}.
\tag{CE2}
```

Let `E_i(B)` denote the unordered coordinate-`i` edges in the
oriented-projective cube (equivalently, take half of the edge sum on the full
lift), and put `w_omega=e^{-beta delta_B(omega)}`.  Summing (CE2) over child
states gives

```math
\boxed{
\sum_{\{u,\tau_i u\}\in E_i(B)}\sqrt{w_u w_{\tau_i u}}
=e^{-\beta d_i}D_\beta(C_i).
}
\tag{CE3}
```

The factorization already proved in (10.614) is equivalently

```math
\kappa_{\beta,i}(B)
=d_i+\frac1\beta
\log\frac{D_\beta(B)}{2D_\beta(C_i)}.
\tag{CE4}
```

Normalize the full lifted weights to the parent Gibbs probability

```math
\nu_B(\sigma,x)
=\frac{e^{\beta\sigma x^{\mathsf T}Bx}}
{\sum_{\sigma',x'}e^{\beta\sigma' x'^{\mathsf T}Bx'}},
```

and let `tau_i` flip `x_i`.  The Bhattacharyya affinity between `nu_B` and
its coordinate flip is

```math
\operatorname{BC}_i(\nu_B)
=\sum_{\sigma,x}
\sqrt{\nu_B(\sigma,x)\nu_B(\sigma,\tau_i x)}.
```

Equations (CE3)--(CE4), with the two directed endpoints per unordered edge,
give the exact cancellation of `d_i`:

```math
\boxed{
e^{-\beta\kappa_{\beta,i}(B)}
=\operatorname{BC}_i(\nu_B)
=\frac{2e^{-\beta d_i}D_\beta(C_i)}{D_\beta(B)}.
}
\tag{CE5}
```

Thus the selected-child target at one state is exactly

```math
\kappa_{\beta,i}(B)\ge T-K
\quad\Longleftrightarrow\quad
\operatorname{BC}_i(\nu_B)\le e^{-\beta(T-K)},
\qquad
T=\alpha[r^{3/2}-(r-1)^{3/2}].
\tag{CE6}
```

In words, one coordinate of the square root of the parent Gibbs density must
be almost orthogonal to its flip, at the exponentially small scale
`exp(-Theta(sqrt(r)))`.

## 2. Every deficit layer has the two-extension injection

Let `N_B(t)` count oriented-projective states of deficit at most `t`.  From
(CE1), the better extension has deficit at most `d_i+delta`, while both
extensions have deficit at most `2(d_i+delta)`.  Projection recovers the
child state, so the maps are injective.  Hence, for every `t>=0`,

```math
\boxed{
N_B(d_i+t)\ge N_{C_i}(t),
\qquad
N_B(2(d_i+t))\ge2N_{C_i}(t).
}
\tag{CE7}
```

The second inequality at `t=0` is (10.613).  There is also an aggregate
edge-isoperimetric consequence.  On the full lift write `\widetilde N=2N`.
If `d_i<=D` for all `i`, every child state of deficit at most `t` produces an
internal edge of the parent layer of deficit at most `2(D+t)`.  The cube
edge-isoperimetric inequality therefore gives

```math
\sum_i\widetilde N_{C_i}(t)
\le e\!\left(\widetilde S_B(2(D+t))\right)
\le\frac12\widetilde N_B(2(D+t))
\log_2\widetilde N_B(2(D+t)).
\tag{CE8}
```

This is a genuine decrement-versus-degeneracy statement, but the doubled
deficit threshold and the logarithmic edge bound are too weak by themselves
to imply (CE6).

At zero temperature, (CE5) recovers the precise prefactor refinement

```math
\kappa_{\beta,i}
=d_i+\frac1\beta\log\frac{g(B)}{2g(C_i)}+o(\beta^{-1}),
\tag{CE9}
```

where the displayed form applies with the usual interpretation when the
leading ground-edge coefficient is nonzero.  In particular, a flat step has
`g(B)>=2g(C_i)` but its entropic payment vanishes as `beta` tends to infinity.

## 3. What generic cube functional inequalities actually give

Let `p_sigma` be the parent Gibbs orientation marginal and let `nu_sigma` be
the conditional law of `x`.  For the uniform cube law `u`, set

```math
\mathcal H_x(\nu_B)
=\sum_\sigma p_\sigma
D_{\mathrm{KL}}(\nu_\sigma\Vert u).
```

Apply the sharp tensorized two-point log-Sobolev inequality to
`F_sigma=sqrt(nu_sigma/u)`.  Since

```math
\frac12\mathbb E_u(F_\sigma-F_\sigma\circ\tau_i)^2
=1-\operatorname{BC}_i(\nu_\sigma),
```

averaging over `sigma` and using (CE5) gives

```math
\boxed{
\mathcal H_x(\nu_B)
\le\sum_{i=1}^r
\left(1-e^{-\beta\kappa_{\beta,i}(B)}\right).
}
\tag{CE10}
```

Consequently this scalar inequality certifies only

```math
\max_i\kappa_{\beta,i}(B)
\ge-\frac1\beta
\log\left(1-\frac{\mathcal H_x(\nu_B)}r\right).
\tag{CE11}
```

Even the formal maximal input `mathcal H_x<=r log 2` makes the right-hand
side at most the constant `-beta^{-1}log(1-log 2)`.  Thus this direct
log-Sobolev scalarization cannot reach the `Theta(sqrt(r))` target.  This is
a limitation of (CE10)--(CE11), not a proof that every entropy argument must
fail.

There is a sharper generic obstruction to using only support entropy and
cube geometry.  Let `S_{r,k}` be the union of the Hamming balls of radius
`k<r/2` around the two antipodal corners, and let `nu` be uniform on
`S_{r,k}`.  It is antipodally invariant and coordinate-transitive.  Direct
edge counting gives, for every `i`,

```math
\operatorname{BC}_i(\nu)
=\frac{2\sum_{j=0}^{k-1}\binom{r-1}{j}}
{\sum_{j=0}^{k}\binom rj}.
\tag{CE12}
```

For `k=floor(sqrt(r))`, this is asymptotic to `2/sqrt(r)`, so every
Bhattacharyya reward is only

```math
-\frac1\beta\log\operatorname{BC}_i
=\frac{\log r}{2\beta}+O_\beta(1)
=o(\sqrt r),
```

although `log|S_{r,k}|=O(sqrt(r) log r)=o(r)`.  Therefore antipodality,
subexponential support, projection/edge entropy, and coordinate symmetry do
not by themselves imply (CE6).  The Hamming-ball measure is not asserted to
be the Gibbs measure of a quadratic signing, let alone a restriction of a
global minimizer.  A successful theorem must exploit exactly that missing
algebraic/ancestral structure.

## 4. Finite audits and the `A9` wall

`tmp/check_cavity_entropy_r18.py` checks (CE1)--(CE5) by exact integer energy
enumeration and floating evaluation only after the polynomial identities are
matched.  It passes:

- 3,048 restriction-coordinate-temperature instances for all principal
  restrictions of `A8`;
- 6,885 such instances for all principal restrictions of `A9`;
- temperatures `beta=0.1,0.5,1`.

For `A8`, every decrement is `2`, and all eight rewards are respectively
`1.200847450`, `2.257889866`, and `2.261149326` at those temperatures.  For
`A9`, all decrements are zero; at `beta=0.5` its maximum reward is
`2.665414766`.

With `s=e^{-beta}`, (CE5) converts the existing exact certificate

```math
2Z_{C_i}(s)-s^3Z_{A_9}(s)>0
```

into

```math
\operatorname{BC}_i(\nu_{A_9})>s^3=e^{-3\beta},
```

which is exactly `kappa_{beta,i}(A9)<3`.  Thus the Hellinger formulation
fully preserves the all-temperature `A9` wall and merely identifies the
positive constant slack which any theorem must allow.

## 5. Quantifiers and ledger-ready conclusion

Equations (CE3)--(CE7) hold for every fixed `beta>0`, every signing `B`, and
every coordinate.  They do not select a path.  For convergence, the needed
statement remains: one universal fixed `beta`, one universal `K`, every
global minimizer root, and a deletion path whose every prefix down to every
fixed-proportion landing order satisfies (CE6) with the root-frozen
`alpha=q_n/n^(3/2)`.  A beta chosen from the root, a different beta at each
order, one favorable root/path, or one uncontrolled landing order is
insufficient.

The ledger-ready conclusion is:

> **Verified:** soft cavity reward is exactly the negative logarithm of the
> coordinate-flip Bhattacharyya affinity of the parent oriented Gibbs law.
> Every child deficit layer embeds by one good extension and by a full cube
> edge at doubled deficit.  Standard scalar cube log-Sobolev and generic
> projection/support entropy provide only constant or logarithmic reward
> information and cannot alone reach the required square-root scale;
> antipodal Hamming balls witness this scoped limitation.  This does not
> falsify the constant-shortfall theorem for quadratic Gibbs measures
> inherited from global minimizers.  The sharpened open target is a
> minimizer-specific near-isolation theorem forcing some coordinate affinity
> at most `exp(-beta[alpha Delta_r-K])` uniformly along a selected path.

