# Recurrence audit for the raw negative-overlap hypothesis

Date: 2026-08-18.

Status: rigorous derivation from NT.3, TR.2, and the archived
almost-subadditivity theorem.  This note records an essential scope point:
decay of the raw overlap alone does **not** imply the pressure recurrence.

## 1. Setup and exact error propagation

Fix `beta,lambda>0`, let `N=m+n`, and assume the split is comparable.  Put

```math
t={\beta\over\sqrt N}.
```

Let `A,D` be minimizing children at the contracted scaled temperatures, and
fix an orientation.  Write `L` for its bridge pressure,

```math
V_\lambda=-{1\over\lambda}\log E_Ue^{-\lambda L},
\qquad
V_\lambda^{\rm row}
=\inf_{p\text{ row product}}
 \{E_pL+\lambda^{-1}D(p\Vert U)\},
```

and take the same-temperature target

```math
T=P_m(\beta)+P_n(\beta).
```

Suppose the full soft bridge reaches this target with an upper error `E_N`:

```math
V_\lambda\le T+E_N,
\qquad E_N\ge0.                                      \tag{R.1}
```

The raw negative-overlap theorem gives

```math
0\le V_\lambda^{\rm row}-V_\lambda
={\mathcal I_\lambda^{\leftarrow}\over\lambda}
\le C_{\rm LS}\lambda t^2mn\widehat\rho_N^-(\lambda).
                                                            \tag{R.2}
```

Consequently, for `Delta_N=V_lambda^row-T`,

```math
(\Delta_N)_+
\le d_N:=E_N+C_{\rm LS}\lambda\beta^2{mn\over N}
                  \widehat\rho_N^-(\lambda).       \tag{R.3}
```

This verifies the scaling: on comparable splits, the overlap contribution
to the product target excess is `O(N rhohat_N^-)`, not a square root of that
quantity and not `O(N^2 rhohat_N^-)`.

Let `L_0=L(0)` and `h=T-L_0`.  Since `L_0<=V_lambda`, (R.1) gives
`(-h)_+<=E_N`.  In TR.2 choose

```math
a_N=2d_N+1.                                         \tag{R.4}
```

Then `a_N>Delta_N`, `h+a_N>0`, and the row-product basin satisfies

```math
U\{L\le T+a_N\}
\ge \exp\!\left(-\lambda^2\beta^2{mn\over N}\right)
     \left({a_N-\Delta_N\over h+a_N}\right)^2.
                                                            \tag{R.5}
```

Because `0<=T<=beta^2N/4` and the numerator is at least one, the last
factor is at least `c_(beta,lambda)/N^2` whenever `d_N=O(N)`.  Thus (R.5)
is a linearly rare basin.  More importantly for the scalar recurrence, its
nonemptiness gives directly

```math
\boxed{
P_N(\beta)\le P_m(\beta)+P_n(\beta)
 +2E_N+2C_{\rm LS}\lambda\beta^2{mn\over N}
       \widehat\rho_N^-(\lambda)+1.}               \tag{R.6}
```

No orientation loss occurs: if the joint two-orientation soft value obeys
(R.1), one deterministic orientation has a value no larger than the joint
soft value.  There is also no additional temperature-conversion term in
(R.6): the children are minimized at raw temperature `beta/sqrt(N)`, and
the definition of `T` already includes the complete contracted-to-same-
temperature retuning cost.

## 2. Power rates and the sharp summability requirement

If uniformly on all sufficiently large comparable splits

```math
\widehat\rho_N^-(\lambda)\le C_\rho N^{-\alpha},
\qquad
E_N\le C_E N^{1-\gamma},                           \tag{R.7}
```

with `alpha,gamma>0`, then (R.6) is

```math
P_N(\beta)\le P_m(\beta)+P_n(\beta)
 +C N^{1-\delta},
\qquad
\delta=\min\{\alpha,\gamma,1\}.                  \tag{R.8}
```

Every `alpha>0` is therefore sufficient; there is no positive threshold
larger than zero.  At merge level `j` of the balanced tree, the normalized
overlap error is `O((2^jk)^(-alpha))`, whose sum is `O(k^(-alpha))`.

For non-power rates, put

```math
q(N)={E_N\over N}+\widehat\rho_N^-(\lambda)+{1\over N}.
```

A sufficient and essentially sharp condition for this proof architecture is
the dyadic Dini condition

```math
\lim_{k\to\infty}
\sum_{j\ge0}
 \sup_{2^jk\le N\le C2^jk}q(N)=0                 \tag{R.9}
```

for a fixed shell constant `C` determined by the balanced merging scheme.
For an eventually monotone envelope this is equivalent to the usual
Hammersley condition `sum_N q(N)/N<infinity`.  Thus
`q(N)=O((log N)^(-1-epsilon))` works, while `q(N)=O(1/log N)` and an
unspecified `o(1)` do not suffice.  The archived slowly oscillating linear
sequence shows that this is a real logical issue, not just a proof artifact
of one estimate.

Under (R.8), or under (R.9), the balanced-tree theorem plus the adjacent-
order `O_beta(1)` modulus proves convergence of `P_n(beta)/n` for this fixed
`beta`.  If the hypotheses hold for every fixed `beta` (it is enough along
an unbounded sequence of fixed `beta` values), then

```math
\beta {M_n\over n^{3/2}}-\log2
\le {P_n(\beta)\over n}
\le \beta {M_n\over n^{3/2}}
```

implies

```math
\limsup_n {M_n\over n^{3/2}}
-\liminf_n {M_n\over n^{3/2}}
\le {\log2\over\beta}.
```

Sending the fixed parameter `beta` to infinity proves convergence of the
normalized ground-state sequence.  Constants and powers may deteriorate
with `beta`; no uniformity in `beta` is needed.

## 3. Logical ceiling: overlap decay alone gives no recurrence

The premise (R.1) is indispensable.  NT.3 controls the amount by which the
negative tilt or a row-product law can lower the fair bridge pressure:

```math
0\le E_UL-V_\lambda^{\rm row}
\le C_{\rm LS}\lambda t^2mn\widehat\rho_N^-(\lambda),
```

and similarly `E_UL-V_lambda` is bounded from above.  These inequalities
give a **lower** bound on `V_lambda` in terms of `E_UL`; they do not show
that `V_lambda` is below the same-temperature target `T`.  If `E_UL-T` is
linear, small overlap says precisely that the fixed negative tilt cannot
remove that linear excess.

There is an even sharper redundancy warning.  If (R.1) itself holds, then

```math
P_N(\beta)\le\min_B L(B)\le V_\lambda\le T+E_N.    \tag{R.10}
```

Thus a power-saving `E_N` already proves a sharper recurrence than (R.6),
without using overlap.  Exact target reach (`E_N=0`) already gives exact
same-temperature subadditivity.  If only `E_N=o(N)` is known without a
summable rate, adding a power-saving overlap bound cannot remove the
nonsummable `E_N` term in (R.3).

Therefore the standalone hypothesis

```math
\widehat\rho_N^-(\lambda)=O(N^{-\alpha})
```

is not by itself a Level-6 implication.  It is a sharp **no-gain/product-
approximation theorem**.  To become a nonredundant recurrence mechanism it
must be coupled to a separate statement that relates the fair-pressure
excess to the target without already asserting (R.1), or to a branch theorem
which turns failure of target reach into a constructive parent.  Neither is
contained in the current raw-overlap inequalities.
