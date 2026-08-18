# A common compact-tilt certificate for the full square carrier

Status: **rigorous actual-child full-carrier response reduction**.  The
binary spike alphabet is unnecessary.  Relative entropy to the fair bridge
law replaces query-word counting and yields one common compact-tilt
edge-cavity overlap curve for every bounded-`L^2` row-product carrier.

## 1. Setup

Let `A,D` be actual contracted-temperature minimizing children of orders
`m,n`, let `N=m+n`, and put `t=beta/sqrt(N)`.  Let `L` be their exact bridge
pressure (RT.1).  Fix the linear cap from Theorem RT.2,

```math
F(B)=L(B)\wedge CN,                               \tag{CT.1}
```

where `C>3beta^2/8+2log K`.  For the full row-product carrier write

```math
\mathcal P_K=\left\{
 P=\bigotimes_{i=1}^m q_iU_n:
 E_Uq_i=1,\ \|q_i\|_2\le K
 \right\}.                                      \tag{CT.2}
```

This contains the complete exact degree-`2d` square carrier of Theorem
37.39 for its fixed `K=K_1`, not merely a finite net or a binary slice.
Every `P in P_K` obeys the dimension-free-per-row entropy bound

```math
D(P\Vert U_{mn})
=\sum_iD(q_iU_n\Vert U_n)
\le2m\log K=:D_0.                                \tag{CT.3}
```

For a bridge edge `a`, let `r_a(B_{-a})` be the exact edge-deleted Gibbs
response FI.14, and put

```math
\Gamma_F(B)=\sum_a[D_aF(B)]^2,
\qquad D_af(B)={f(B)-f(B^a)\over2}.              \tag{CT.4}
```

Clipping is one-Lipschitz, so the exact cavity identity gives

```math
\boxed{
\Gamma_F(B)\le\Gamma_L(B)
\le t^2\sum_a r_a(B_{-a})^2.}                   \tag{CT.5}
```

Define one pair of common capped tilts, independent of the queried carrier,

```math
{d\Pi_{s,F}\over dU_{mn}}
={e^{sF}\over E_Ue^{sF}},
\qquad s\in\mathbb R,                           \tag{CT.6}
```

and their normalized actual-child overlap envelope

```math
\rho_N(S)={1\over mn}
 \sup_{|s|\le S}E_{\Pi_{s,F}}\sum_a r_a^2.       \tag{CT.7}
```

This is a scalar curve on a compact interval.  It does not contain carrier
coefficients, a carrier net, or the pressure response table.

## 2. Entropy replaces carrier-word counting

The fair cube logarithmic-Sobolev inequality and the elementary exponential
gradient estimate imply, for both signs and every `s>=0`,

```math
\log E_Ue^{\pm s(F-E_UF)}
\le C_{\rm LS}s^2
 \sup_{|u|\le s}E_{\Pi_{u,F}}\Gamma_F,            \tag{CT.8}
```

where one may take an absolute constant `C_LS=1`.  Indeed, if
`psi(s)=log E exp(s(F-EF))`, then

```math
{d\over ds}{\psi(s)\over s}
\le C_{\rm LS}E_{\Pi_{s,F}}\Gamma_F,
```

and integration proves (CT.8); apply the same argument to `-F`.

For any `P in P_K`, entropy duality gives

```math
E_PF-E_UF
\le {D(P\Vert U)+\log E_Ue^{s(F-E_UF)}\over s},  \tag{CT.9}
```

and the same inequality applied to `-F` controls the opposite sign.  Combining
(CT.3), (CT.5), and (CT.8) proves the main theorem.

### Theorem CT.1 (full-carrier compact-tilt overlap theorem)

For every `0<s<=S`,

```math
\boxed{
\operatorname {range}_{P\in\mathcal P_K}E_PF
\le2\left{{D_0\over s}
 +C_{\rm LS}t^2mn\rho_N(S)s\right}.}            \tag{CT.10}
```

By the uniform response truncation RT.2,

```math
\boxed{
\operatorname {range}_{P\in\mathcal P_K}E_PL
\le O_\beta(N^{3/2}e^{-cN})
+2\left{{D_0\over s}
 +C_{\rm LS}t^2mn\rho_N(S)s\right}.}            \tag{CT.11}
```

If

```math
s_*=sqrt{{D_0\over C_{\rm LS}t^2mn\rho_N(S)}}\le S,
```

then

```math
\boxed{
\operatorname {range}_{P\in\mathcal P_K}E_PL
\le O_\beta(N^{3/2}e^{-cN})
 +4\sqrt{C_{\rm LS}D_0t^2mn\rho_N(S)}.}          \tag{CT.12}
```

On balanced splits, `D_0=O(N)` and `t^2mn=Theta(N)`.  Therefore
`rho_N(S_N)=O(N^{-alpha})` on a window containing
`s_*=O(N^(alpha/2))` gives the full-carrier power saving

```math
\boxed{
\operatorname {range}_{P\in\mathcal P_K}E_PL
=O(N^{1-\alpha/2}).}                             \tag{CT.13}
```

No dependence on the coefficient dimension
`sum_(j<=2d)binom(n,j)`, its covering number, or its `m`-fold tensor appears.

There is an exact converse at the linear scale.  If

```math
\operatorname {range}_{P\in\mathcal P_K}E_PL\ge\eta N,
```

then, for all large `N`, at least one carrier product differs from the fair
response by `eta N/3` after clipping.  Entropy duality with

```math
S_{\eta,K}={12(1+\log K)\over\eta}               \tag{CT.14}
```

forces one sign and one `|s|<=S_(eta,K)` for which

```math
\boxed{
{1\over mn}E_{\Pi_{s,F}}\sum_a r_a^2
\ge c_{\eta,K,\beta}>0}                         \tag{CT.15}
```

on balanced splits.

*Proof of the converse.*  Choose `P,Q` witnessing the range.  At least one,
say `P`, has
`|E_PF-E_UF|>=eta N/3` for large `N`, by RT.2.  For the matching sign
`sigma`, the reverse form of (CT.9) gives

```math
\log E_Ue^{\sigma S_{\eta,K}(F-E_UF)}
\ge S_{\eta,K}{\eta N\over3}-D_0
\ge2N.                                          \tag{CT.16}
```

The differential form of (CT.8) then gives some
`|s|<=S_(eta,K)` with `E_(Pi_(s,F))Gamma_F>=c_(eta,K)N`.
Equation (CT.5) and `t^2mn=Theta_beta(N)` prove (CT.15). `square`

The entropy-regularized product objective has an immediate two-scalar
version.  For

```math
\mathcal G(P)=E_PL+{1\over\lambda}D(P\Vert U_{mn}),
```

if two carrier products satisfy
`|G(P)-G(Q)|>=eta N`, then either

```math
|D(P\Vert U)-D(Q\Vert U)|\ge {\lambda\eta N\over2}, \tag{CT.17}
```

which is a directly computable sum of row entropies, or their pressure
responses differ by at least `eta N/2`, in which case (CT.15) holds with
adjusted constants.  Thus a linear full-carrier retuning is witnessed by
either one row-additive entropy scalar or the common compact-tilt overlap
curve; no carrier response table is needed for this dichotomy.

## 3. Why this escapes the finite-net obstruction

The fixed-degree square carrier contains `exp(Omega(n))` separated row
densities, so a uniform net word can have probability `exp(-Omega(mn))`.
That makes the binary-word proof of MT.3 unusable on the complete carrier.
Theorem CT.1 never chooses a prior on carrier labels.  It charges a declared
product directly by its physical relative entropy, which is only `O(m)` by
(CT.3), regardless of coefficient dimension.  This is why a bounded compact
tilt survives the full carrier packing.

The theorem is also sharp in scope.  The generic weak-bias rank-one pressure
from the finite-carrier audit has a linear response range and therefore a
positive moderate-tilt overlap, exactly as (CT.15) predicts.  CT.1 does not
claim that bounded row complexity alone makes `rho_N` vanish.

For the actual optimizing children, however, (CT.7) is one common,
low-information target for **all** recovered square factors.  The remaining
minimal lemma is now:

> prove a power-saving bound on the fair-base capped-tilt edge-cavity overlap
> `rho_N(S_N)` for actual contracted-temperature minimizers on a growing
> moderate window, or exhibit an actual minimizing sequence with a positive
> compact-tilt overlap and determine its product-phase branch.

Known zero-field child minimality verifies (RT.2), (CT.3), and the exact
cavity identity (CT.5).  It does **not** currently verify any decay of
`rho_N`; that is an external-field statement.  Thus CT.1 is a strict
response reduction and a full-carrier extension of MT.3--MT.4, but not a
Level-6 recurrence.
