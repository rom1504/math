# Wave 56 main calculation: the entropy-matched mesoscopic core

## Outcome

For `ell=o(n)` with `ell>>sqrt(n)`, the common-core transition is a moderate
deviation of the ordinary Johnson intersection.  The degree-two baseline
load, and also any fixed fraction of the full transition load, can be carried
at the same leading exponential scale by

```math
\binom nm\exp\left{-\kappa_p\frac{\ell^2}{n}
+O\left(\frac\ell{\sqrt n}+\log n\right)\right}
```

closest partners, and asymptotically fewer partners cannot carry constant
load.  Here

```math
\kappa_p=\frac{(1-p)^2}{2p^2}.
```

Thus the entropy-matched scale for the project budget
`H=n^(3/4-c)` is

```math
\boxed{\ell_*\asymp\sqrt{nH}=n^{7/8-c/2}.}
```

It satisfies `ell_*/H -> infinity`, so it is genuinely outside the direct
collision regime (10.1270), while its required partner fraction is exactly
`e^{-Theta(H)}` rather than `e^{-Theta(n)}` as at a linear core.  This does
not prove that favorable families of an exact minimizer have the required
histogram, but it identifies a quantitatively better target scale.

One-threshold compression has a second, different phase boundary.  Its
maximum possible load is `Theta(sqrt(n)/ell)`, while the required baseline is
`lambda_2=Theta(ell^2/n^2)`.  Hence it is count-feasible precisely below the
scale `ell=Theta(n^(5/6))`.  At `ell_*`, the entropy-matching and threshold-
capacity windows overlap at exponent level exactly when `c>=1/12`.  Thus for
the usable range `1/12<c<1/4`, a mesoscopic one-threshold attack is reopened;
above `n^(5/6)` the complete multilevel histogram is essential.

## 1. Exact transition parameters

Fix `m/n -> p in (1/2,1)` and a core size `ell=o(n)`.  For a fixed base
selector `S`, the ordinary uniform-slice intersection `J_0=|S intersect T|`
and the `K_ell` transition intersection `J_ell` obey

```math
J_0\sim\operatorname{Hypergeom}(n,m,m),
```

and

```math
J_\ell-\ell\sim
\operatorname{Hypergeom}(n-\ell,m-\ell,m-\ell).
```

Consequently

```math
\mu_0=\frac{m^2}{n},\qquad
\mu_\ell=\ell+\frac{(m-\ell)^2}{n-\ell},
```

with the exact shift

```math
\boxed{
\delta_\ell:=\mu_\ell-\mu_0
=\frac{\ell(n-m)^2}{n(n-\ell)}
=(1-p)^2\ell\{1+o(1)\}.}
\tag{M56.1}
```

The exact transition variance is

```math
\boxed{
\operatorname{Var}(J_\ell)
=\frac{(m-\ell)^2(n-m)^2}
{(n-\ell)^2(n-\ell-1)}
=p^2(1-p)^2n\{1+o(1)\}.}
\tag{M56.2}
```

The ordinary variance has the same leading term.  The kernel eigenvalues and
core density satisfy

```math
\rho_\ell=\frac{(m)_\ell}{(n)_\ell},\qquad
\lambda_1(\ell)=\frac{\ell(n-m)}{m(n-\ell)},
```

and

```math
\lambda_2(\ell)=
\frac{\ell(\ell-1)(n-m)(n-m-1)}
{m(m-1)(n-\ell)(n-\ell-1)}.
```

Uniformly for `ell=o(n)`,

```math
-\log\rho_\ell=\ell\log(1/p)+O(\ell^2/n),
\qquad
\lambda_1\sim\frac{1-p}{p}\frac\ell n,
\qquad
\lambda_2\sim\lambda_1^2.
\tag{M56.3}
```

The extraction denominator is therefore

```math
D_\ell=(1-\lambda_1)+n(\lambda_1-\lambda_2)
\sim\frac{1-p}{p}\ell.
\tag{M56.4}
```

Its polynomial loss is negligible on the `H` entropy scale.

## 2. Moderate-deviation partner capacity

Let `N=binom(n,m)` and

```math
A_j=\binom mj\binom{n-m}{m-j}
=N\Pr\{J_0=j\}.
```

The transition mass of shell `j` is

```math
A_jK_\ell(j)=\Pr\{J_\ell=j\}.
```

An entropy Taylor expansion at the ordinary mean, evaluated at
`mu_ell+O(sqrt(n))`, gives

```math
\log A_j
=\log N-\frac{\delta_\ell^2}
{2p^2(1-p)^2n}
+O\left(\frac\ell{\sqrt n}+\frac{\ell^3}{n^2}+\log n\right).
```

Since `ell=o(n)`, the leading moderate-deviation cost is

```math
\boxed{
\frac{\delta_\ell^2}{2p^2(1-p)^2n}
=\left\{\frac{(1-p)^2}{2p^2}+o(1)\right\}
\frac{\ell^2}{n}.}
\tag{M56.5}
```

The `K_ell` transition has width `Theta(sqrt(n))`.  Every fixed-probability
transition window therefore contains

```math
\boxed{
N\exp\left[-\kappa_p\frac{\ell^2}{n}
+O\left(\frac\ell{\sqrt n}+\frac{\ell^3}{n^2}+\log n\right)\right]
}
\tag{M56.6}
```

distinct partners.  Collecting only the polynomial degree-two mass
`lambda_2=Theta(ell^2/n^2)` changes the transition quantile by
`O(sqrt(n log(n/ell)))` and changes the logarithm above only by
`O((ell/sqrt(n))sqrt(log(n/ell))+log n)`.  Thus its leading exponential
capacity is the same.  Conversely, closest-shell rearrangement gives the
matching lower capacity boundary.  A fully uniform statement should retain
the cubic entropy error when `ell` approaches
`n^(2/3)` or larger; at `ell=sqrt(nH)` that error is
`ell^3/n^2=H^(3/2)/sqrt(n)=o(H)`.  The exact entropy rate, rather than only
the quadratic approximation, should still be used for lower-order constants.
Its ratio to the quadratic term is `ell/n=o(1)`, so (M56.5) is uniformly
correct at relative `1+o(1)` scale.

Equating the leading loss to `H` yields

```math
\ell^2/n\asymp H,
\qquad \ell\asymp\sqrt{nH}.
\tag{M56.7}
```

Because

```math
\frac{\sqrt{nH}}H=\sqrt{n/H}=n^{1/8+c/2}\to\infty,
```

this satisfies the non-direct requirement `-log rho_ell=omega(H)`.

## 3. Sharp order and phase boundary of one-threshold compression

For a threshold `s` in the transition-typical window, the ordinary Johnson
shell is in an upper moderate-deviation tail.  Consecutive shell ratios give

```math
\frac{\sum_{j\ge s}A_j}{A_s}=\Theta(n/\ell)
```

uniformly for `ell>>sqrt(n)`, while the largest `K_ell` transition atom is
`Theta(n^(-1/2))`.  Hence

```math
\boxed{
\sup_s\left\{\sum_{j=s}^{m-1}A_j\right\}K_\ell(s)
=\Theta\left(\frac{\sqrt n}{\ell}\right).}
\tag{M56.8}
```

The optimizer can shift below `mu_ell` by `Theta(n/ell)`, which is
`o(sqrt(n))`, and does not change the order.  Thresholds outside the typical
window lose exponentially in the entropy comparison or in transition mass.

The target in the spectral-excess theorem is not constant load.  From
(M56.3),

```math
\lambda_2(\ell)\sim
\left(\frac{1-p}{p}\right)^2\frac{\ell^2}{n^2}.
\tag{M56.9}
```

The exact local-limit and Mills-ratio constants sharpen (M56.8) to

```math
\sup_s\mathcal B_{n,m}(s)K_\ell(s)
\sim\frac{p}{\sqrt{2\pi}(1-p)}\frac{\sqrt n}{\ell}.
\tag{M56.10}
```

Consequently the ratio of the largest whole-slice one-threshold load to the
degree-two baseline is

```math
\boxed{
\frac{\sup_s\mathcal B_{n,m}(s)K_\ell(s)}{\lambda_2(\ell)}
\sim\frac1{\sqrt{2\pi}}
\left(\frac p{1-p}\right)^3\frac{n^{5/2}}{\ell^3}.}
\tag{M56.11}
```

Thus whole-slice threshold capacity is impossible for

```math
\ell\gg (2\pi)^{-1/6}\frac p{1-p}n^{5/6},
```

and is not ruled out by counting below that scale.  At
`ell_*=sqrt(nH)`, (M56.11) has order

```math
\frac{n}{H^{3/2}}=n^{-1/8+3c/2}.
```

It diverges exactly when `c>1/12`, is critical at `c=1/12`, and vanishes
when `c<1/12`.  Therefore a nonempty entropy-matched one-threshold window is

```math
\sqrt{nH}\lesssim\ell\lesssim n^{5/6},
\qquad 1/12<c<1/4.
\tag{M56.12}
```

This is only capacity feasibility: a signing-specific theorem must still
place the required favorable partners in the appropriate upper-overlap ball.

## Research judgment

The scale `ell~sqrt(nH)` is the first core scale whose intrinsic partner
entropy loss matches the saved selector budget.  It is more plausible than a
fixed linear core for a constructive overlap theorem.  For `c>1/12`, even a
single optimized overlap threshold is no longer defeated by whole-slice
capacity; this is a genuinely different regime from the Wave 55 linear-core
no-go.  The exact remaining lemma is signing-specific: at a controlled row
cap, prove the required `K_ell` weighted load, either through that mesoscopic
threshold ball or through the complete typical histogram.  Capacity and
kernel asymptotics alone supply neither structure.
