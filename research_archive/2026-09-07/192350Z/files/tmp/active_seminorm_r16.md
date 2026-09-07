# Wave 16: active-state response versus global `Q`-distance

**Status.**  The deterministic identities and inequalities below are exact.
The `A_9` claims and the all-positive finite audit are checked exhaustively by
`tmp/check_active_seminorm_r16.py`.  The structural rectangularization lemma
is exact.  No unconditional critical-scale partition theorem is proved: the
surviving signed compatibility target is stated explicitly in Section 6.

## 1. Oriented states and three localized functionals

Use the ledger normalization

```math
Q(M)=\max_{z\in\{\pm1\}^n}|z^{\mathsf T}Mz|.
```

It is useful to retain the orientation.  Let

```math
\Omega_n=\{(\sigma,z):\sigma\in\{\pm1\},\ z\in
\{\pm1\}^n/\{\pm1\}\},
\qquad e_M(\sigma,z)=\sigma z^{\mathsf T}Mz.
```

Thus `Q(M)=max_{omega in Omega_n} e_M(omega)`.  For a fixed cross
mosaic `C`, define its oriented deficit and its `t`-near-active set by

```math
\delta_C(\omega)=Q(C)-e_C(\omega)\ge0,
\qquad
\mathcal A_C(t)=\{\omega:\delta_C(\omega)\le t\},
\quad t\ge0.
\tag{R16.1}
```

There are three closely related objects, and conflating them loses the
important sign.

First, the signed active support is

```math
h_{C,t}(K)=\max_{\omega\in\mathcal A_C(t)}e_K(\omega).
\tag{R16.2}
```

Second, its nonnegative convexification is the one-sided active sublinear
gauge

```math
\mathfrak a^+_{C,t}(K)
=\max_{\omega\in\mathcal A_C(t)}[e_K(\omega)]_+
=[h_{C,t}(K)]_+.
\tag{R16.3}
```

For fixed `(C,t)`, this functional vanishes at zero, is subadditive, and is
positively homogeneous for nonnegative scalars.  It is not symmetric, so
“one-sided sublinear gauge” or “hemi-seminorm” is more precise than norm.  Its
symmetrization is the genuine seminorm

```math
\mathfrak a_{C,t}(K)
=\max\{\mathfrak a^+_{C,t}(K),\mathfrak a^+_{C,t}(-K)\}
=\max_{\omega\in\mathcal A_C(t)}|e_K(\omega)|.
\tag{R16.4}
```

Finally, the deficit-corrected active response is

```math
\mathcal R_{C,t}(K)
=\max_{\omega\in\mathcal A_C(t)}
\{e_K(\omega)-\delta_C(\omega)\}.
\tag{R16.5}
```

This last object is the sharp one for finite perturbations.  It is convex in
`K` and satisfies `R_{C,t}(0)=0`, but the deficit destroys homogeneity, so it
is not a seminorm.  The useful lower-response defect

```math
\Gamma_{C,t}(K)=[-\mathcal R_{C,t}(K)]_+
\tag{R16.6}
```

is also not a seminorm: it is based on the **best** compatible active state,
whereas a convex seminorm majorant necessarily pays for a worst state.

The localization is genuinely weaker than global `Q`-distance:

```math
\mathfrak a^+_{C,t}(K)\le \mathfrak a_{C,t}(K)\le Q(K).
\tag{R16.7}
```

Both inequalities can be arbitrarily strict; Sections 4 and 5 give exact
examples with both one-sided active gauges zero and nonzero, or leading,
`Q(K)`.

## 2. Exact deterministic active/tail sandwich

Put

```math
\Delta_C(K)=Q(C+K)-Q(C).
```

The basic identity is

```math
\boxed{
\Delta_C(K)
=\max_{\omega\in\Omega_n}
\{e_K(\omega)-\delta_C(\omega)\}.
}
\tag{R16.8}
```

Indeed, subtracting `Q(C)` inside the finite maximum defining `Q(C+K)`
gives (R16.8) term by term.  Splitting the state space at deficit `t` gives
the exact decomposition

```math
\Delta_C(K)
=\max\{\mathcal R_{C,t}(K),\mathcal T_{C,t}(K)\},
\qquad
\mathcal T_{C,t}(K)
=\max_{\delta_C(\omega)>t}
\{e_K(\omega)-\delta_C(\omega)\},
\tag{R16.9}
```

where an empty tail has value `-infinity`.  Since `e_K(omega)<=Q(K)`, this
proves the sharp deterministic active/tail sandwich

```math
\boxed{
\mathcal R_{C,t}(K)
\le \Delta_C(K)
\le \max\{\mathcal R_{C,t}(K),\ Q(K)-t\}.
}
\tag{R16.10}
```

On the head, `0<=delta_C<=t`, so

```math
h_{C,t}(K)-t
\le\mathcal R_{C,t}(K)\le h_{C,t}(K).
```

Consequently (R16.10) has the signed-support and one-sided-gauge
corollaries

```math
\boxed{
h_{C,t}(K)-t
\le\Delta_C(K)
\le\max\{h_{C,t}(K),Q(K)-t\},
}
\tag{R16.11}
```

and

```math
\boxed{
-\mathfrak a^+_{C,t}(-K)-t
\le\Delta_C(K)
\le\max\{\mathfrak a^+_{C,t}(K),Q(K)-t\}.
}
\tag{R16.12}
```

For the lower bound in (R16.12), every active value obeys
`e_K(omega)>=-a^+_{C,t}(-K)`.  Thus no symmetry or absolute value is being
silently inserted.

The adjustable margin has a precise role.  To deduce
`Delta_C(K)<=u` from localized information, it is enough to prove

```math
\mathcal R_{C,t}(K)\le u,
\qquad t\ge Q(K)-u.
\tag{R16.13}
```

States beyond the chosen margin can otherwise overtake the old grounds.
The `A_9` capture wall shows that this is a real state switch, not a proof
artifact.

For infinitesimal perturbations no margin is needed:

```math
\lim_{\lambda\downarrow0}
\frac{Q(C+\lambda K)-Q(C)}{\lambda}
=h_{C,0}(K).
\tag{R16.14}
```

This is the standard directional derivative of a finite support function.
The finite increment can have the opposite sign because a near-active state
becomes the new ground.

## 3. What the functional says for `A=C+D`

Suppose `A=C+D` and write `q=Q(A)`.  Substitution in (R16.8) gives

```math
\Delta_C(D)=q-Q(C).
\tag{R16.15}
```

More locally, cancellation of the cross energy gives the exact compatibility
formula

```math
\boxed{
\mathcal R_{C,t}(D)
=\max_{\omega\in\mathcal A_C(t)}e_A(\omega)-Q(C).
}
\tag{R16.16}
```

Hence the cross-only overshoot satisfies

```math
\boxed{
[Q(C)-q]_+
\le\Gamma_{C,t}(D)
\le \mathfrak a^+_{C,t}(-D)+t.
}
\tag{R16.17}
```

At `t=0`, the first upper certificate is simply

```math
\Gamma_{C,0}(D)
=\left[
Q(C)-\max_{\omega:\ e_C(\omega)=Q(C)}e_A(\omega)
\right]_+.
\tag{R16.18}
```

Thus it asks whether **some** oriented `C`-ground is compatible with a high
`A`-energy.  By contrast, the one-sided sublinear majorant is

```math
\mathfrak a^+_{C,0}(-D)
=\left[
Q(C)-\min_{\omega:\ e_C(\omega)=Q(C)}e_A(\omega)
\right]_+,
\tag{R16.19}
```

which pays for the worst `C`-ground.  Equations (R16.18)--(R16.19) isolate
the convexification loss: the genuinely useful response object is signed
and best-state, not a seminorm.

For a block mosaic `D=bigoplus_i D_i` replaced by
`K=bigoplus_i G_i`, equations (R16.8)--(R16.12) apply without change.  At
block scale `s~sqrt(n)`, however,

```math
Q(K)\le\sum_iq_{m_i}=O(n^{5/4})=o(n^{3/2})
```

by (10.567).  Taking `t=Q(K)` already makes the tail nonpositive, and the
whole hybrid increment is subleading by the old triangle bound.  Active
localization therefore does not improve the leading hybrid-selection
question.  Its potentially new use is (R16.17): control the cross overshoot
while allowing the original internal `D` and its captured excess to remain
leading.

## 4. Why the range wall does not automatically transfer

For an equipartition, (10.570) uses

```math
E_{\mathcal P}=D_{\mathcal P}-pA,
\qquad
\sum_iQ(D_i)\le2[pQ(A)+Q(E_{\mathcal P})].
\tag{R16.20}
```

Its proof chooses blockwise positive or negative extrema and then controls
`E_P` at that adaptively assembled global state.  There is no reason for
that state to lie in `A_C(t)`.  Replacing `Q(E_P)` by an active seminorm in
(R16.20) is therefore not a valid inference.

The all-positive adaptive witness makes the missing implication exact.  Let
`A=J-I`, equipartition `n=ks` vertices into `k` equal blocks, assume `k` is
even and `2<=s<n/2`, and put

```math
p=\frac{s-1}{n-1},\qquad D=\bigoplus_iA[V_i],\qquad
C=A-D,\qquad E=D-pA.
```

If `r_i=sum_{v in V_i}z_v` and `T=sum_i r_i`, then

```math
z^{\mathsf T}Cz=T^2-\sum_i r_i^2,
\qquad
z^{\mathsf T}Ez=\sum_i r_i^2-pT^2-n(1-p).
\tag{R16.21}
```

It follows that

```math
Q(C)=n(n-s),
```

with a unique projective oriented ground `(sigma,z)=(+1,1)`.  At that state
`z^T E z=0`, and hence

```math
\boxed{
\mathfrak a^+_{C,0}(E)
=\mathfrak a^+_{C,0}(-E)=0.
}
\tag{R16.22}
```

In fact `C+E=(1-p)A`, so

```math
\boxed{Q(C+E)-Q(C)=0.}
\tag{R16.23}
```

Thus the zero active response is exactly right, even though the global norm
is large.  Give half the blocks each sign and make `z` constant on blocks.
Then `T=0`, every `|r_i|=s`, and

```math
z^{\mathsf T}Ez=\frac{n^2(s-1)}{n-1}.
\tag{R16.24}
```

This is in fact `Q(E)`: the upper bound follows from
`sum_i r_i^2<=ns`, while Cauchy and
`p(n/s)<1` give the matching lower absolute bound.  The state in (R16.24)
has `C`-energy `-ns`; with the orientation making `E` positive its
`C`-deficit is exactly `n^2`.  It is intentionally invisible to the active
functional.  Meanwhile

```math
\sum_iQ(D_i)=n(s-1).
\tag{R16.25}
```

At `s~sqrt(n)`, (R16.24)--(R16.25) are both of leading
`Theta(n^{3/2})` order while both one-sided active gauges are zero.  This example is
not a global minimizer and `Q(A)=Theta(n^2)`; it does **not** prove that a
useful partition exists for minimizers.  It proves the precise logical point
needed here: active response is not comparable to global `Q(E)`, so the
range proof (10.570) cannot simply be replayed with (R16.3) or (R16.6).

## 5. Exact `A_9` audits

### 5.1 The `6+1+1+1` capture wall

Use the partition from (10.579),

```math
V_1=\{0,1,2,4,5,6\},\qquad
V_2=\{3\},\quad V_3=\{7\},\quad V_4=\{8\}.
```

For the original internal matrix `D` and cross mosaic `C`, exact enumeration
gives

```math
Q(C)=22,\qquad Q(D)=22,\qquad
\Delta_C(D)=Q(A_9)-Q(C)=2.
```

The localized profiles are:

| `K` | `t` | `|A_C(t)|` | `h_{C,t}(K)` | `R_{C,t}(K)` | `a^+_{C,t}(K)` | `a^+_{C,t}(-K)` |
|---|---:|---:|---:|---:|---:|---:|
| original `D` | 0 | 3 | 2 | 2 | 2 | 6 |
| original `D` | 4 | 19 | 6 | 2 | 6 | 10 |
| original `D` | 8 | 55 | 10 | 2 | 10 | 14 |

Thus the original increment is already witnessed by an exact `C`-ground,
despite `Q(D)=22`.

There are 40 order-six minimizers `G` for which `Q(C+G)=24`.  Every one has
`Q(G)=10` and the identical localized data

| `t` | `|A_C(t)|` | `h_{C,t}(G)` | `R_{C,t}(G)` | `a^+_{C,t}(G)` | `a^+_{C,t}(-G)` | `Q(G)-t` |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 3 | -6 | -6 | 0 | 6 | 10 |
| 4 | 19 | 6 | 2 | 6 | 10 | 6 |
| 8 | 55 | 10 | 2 | 10 | 10 | 2 |

This gives three exact conclusions.

1. Exact active states alone are insufficient for finite perturbations:
   `a^+_{C,0}(G)=0` but `Delta_C(G)=2`.  This directly falsifies the tempting
   bound `Delta_C(K)<=a^+_{C,0}(K)`.

2. Deficit four is where a new state first realizes the true response
   `R=2`.  For the universal deterministic sandwich, deficit eight is where
   the tail bound also falls to two, so (R16.10) closes exactly.

3. At that closing margin the symmetric localized seminorm has become
   `a_{C,8}(G)=10=Q(G)`.  On this wall, enlarging the active set enough to
   make the generic finite-perturbation bound exact loses all seminorm gain.

Every one of the 40 completions also has 4--13 persistent signed `+24`
grounds shared with the original `A_9`.  Their `(C,D,G)` signed energy
profiles are exactly `(14,10,10)` or `(18,6,6)`.  Consequently

```math
Q(C+(1-\lambda)D+\lambda G)=24
\qquad(0\le\lambda\le1),
```

even though the common-mosaic response baseline grows by the captured
excess twelve.  This independently matches the persistent-ground audit in
`tmp/response_interpolation_r16.md`: a response-to-temporal argument based
only on variation of the global norm along replacement segments cannot see
that baseline change.

### 5.2 The `3+6` wall

For the first-three/last-six split of (10.550), the original matrices have

```math
Q(C)=24,\qquad Q(D)=20,\qquad
\Delta_C(D)=h_{C,0}(D)=0.
```

For the explicit best replacement (10.553), `Q(G)=q_3+q_6=16` and

```math
\Delta_C(G)=4=h_{C,0}(G)=\mathcal R_{C,0}(G).
```

Thus the earlier `3+6` obstruction is already visible at exact `C`-grounds,
whereas the `6+1+1+1` optimal completions require near-active states.  A
single fixed choice of margin cannot be inferred from the profile alone.

## 6. A precise structural obstruction at critical scale

The active set must remain correlated across blocks.  The following exact
lemma shows why a natural independent-block relaxation is self-defeating.

**Rectangularization lemma.**  Let `D=bigoplus_i D_i`, and write

```math
P_i=\max_x x^{\mathsf T}D_i x,
\qquad N_i=-\min_x x^{\mathsf T}D_i x.
```

For each block choose a positive extremizer `x_i^+` and a negative extremizer
`x_i^-`.  Suppose that, for one fixed orientation `sigma`, all `2^k`
concatenations obtained by independently choosing
`x_i in {x_i^+,x_i^-}` belong to `A_C(t)`.  Then

```math
\boxed{
\mathfrak a^+_{C,t}(D)+\mathfrak a^+_{C,t}(-D)
\ge\sum_i(P_i+N_i)
\ge\sum_iQ(D_i).
}
\tag{R16.26}
```

Indeed, for `sigma=+1`, the all-positive-extremizer concatenation gives
`sum_i P_i` to the first gauge and the all-negative-extremizer concatenation
gives `sum_i N_i` to the reverse gauge.  For `sigma=-1` the two roles swap.
This proves (R16.26).  In particular,

```math
\mathfrak a_{C,t}(D)\ge\frac12\sum_iQ(D_i).
\tag{R16.27}
```

At `s~sqrt(n)`, any theorem that first replaces the true active profile by
such a Cartesian product and then proves the symmetric active seminorm is
`o(n^{3/2})` also forces the entire internal norm, and hence captured excess,
to be subleading.  It recreates the range wall by a different route.  This
rules out product-cover, independent-block chaining, and rectangular hull
relaxations whenever they include the two opposite concatenated global
corners at one common orientation.  Merely containing local extrema in
different orientation-dependent pieces is insufficient.

The lemma does **not** defeat the signed best-state defect `Gamma`, because
that defect deliberately retains the correlations and orientation of the
actual active profile.  The genuinely open critical-scale partition target
is therefore:

```math
\boxed{
\begin{gathered}
s\asymp\sqrt n,\qquad
X(\mathcal P)=\sum_i[Q(D_i)-q_{m_i}]
=\Omega(n^{3/2}),\\
t=o(n^{3/2}),\qquad
\Gamma_{C_{\mathcal P},t}(D_{\mathcal P})
=o(n^{3/2}).
\end{gathered}
}
\tag{R16.28}
```

By (R16.17), (R16.28) would imply
`[Q(C_P)-q_n]_+=o(n^{3/2})` without invoking the global centered norm and
without forcing `X(P)` small.  Neither the exact algebra nor the finite
audits prove existence of such a partition for global minimizers.  The
remaining difficulty is precisely to control a **best compatible actual
near-ground**, not a uniform norm and not a Cartesianized family.

## 7. Checker

Run

```bash
.venv/bin/python tmp/check_active_seminorm_r16.py
```

The script:

- verifies (R16.8)--(R16.12) at every deficit breakpoint and neighboring
  integer margin used in the finite examples;
- enumerates all 384 labelled order-six minimizers and all 40 optimal
  `6+1+1+1` completions;
- checks the active profiles, the persistent-ground interpolation data, and
  the explicit `3+6` replacement; and
- verifies the all-positive active/global separation and (10.572) exactly
  after multiplying the centered mask by `n-1`.
