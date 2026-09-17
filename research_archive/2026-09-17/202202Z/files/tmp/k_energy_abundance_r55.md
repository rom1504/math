# Wave 55 main audit: far K-profile abundance is universally too small

## Outcome

The diverging-depth abundance lemmas (10.1277), (10.1285), and the K-branch
of (10.1297) are impossible.  The obstruction is not coefficient
localization: before the exact K-profile is considered, its necessary cross
energy already has uniform local-state probability
`exp{-Omega(omega_n^2 H)}` at depth `r>=omega_n T_n`.

This closes the K-functional tail as an independent annealed implementation.
At bounded depth `r=O(T_n)`, saved local-state mass is already
recurrence-strong by inverse Hanson--Wright.  At diverging depth the required
local states are superexponentially rare on the `H` scale.  The conditional
Montgomery--Smith completion theorem remains correct; what fails is the
proposed saved population of states on which to apply it.

## 1. Uniform cross-energy upper tail

Fix a selector `S`, put `T=S^c`, `m=|S|`, `k=|T|`, and let

```math
B=A[T,S],\qquad R=B^{\mathsf T}B,\qquad
G=\lVert By\rVert_2^2=y^{\mathsf T}Ry.
```

Under a uniform projective local spin `y`, `G` has the same distribution as
under a fully independent Rademacher vector, because it is invariant under
the global sign flip.  Exactly,

```math
\mathbb E_yG=\operatorname{tr}R=mk.
```

For an exact order-`n` minimizer, the verified operator estimate gives

```math
\lVert R\rVert_{\rm op}=\lVert B\rVert_{\rm op}^2
\le\lVert A\rVert_{\rm op}^2\le2q_n,
```

and positivity of `R` gives

```math
\lVert R\rVert_F^2=\operatorname{tr}(R^2)
\le\lVert R\rVert_{\rm op}\operatorname{tr}R
\le2q_nmk.
```

Hanson--Wright therefore proves, uniformly in `S`,

```math
\boxed{
\Pr_y\{G\ge L\}
\le2\exp\left[-c\min\left\{
\frac{(L-mk)^2}{2q_nmk},
\frac{L-mk}{2q_n}
\right\}\right]
}
\tag{R55.1}
```

for `L>=mk`.  No local-energy conditioning or independence between `e` and
`G` is used.

## 2. The far profile has `exp{-Omega(omega^2 H)}` mass

Fix `0<c<1/4` and

```math
H=\lceil n^{3/4-c}\rceil,\qquad T_n=n^{3/2-c}.
```

The exact dual bound `K_H(u)<=sqrt(H)||u||_2` shows

```math
K_H(u)\ge r+b_H
\quad\Longrightarrow\quad
G\ge\frac{(r+b_H)^2}{H}\ge\frac{r^2}{H}.
\tag{R55.2}
```

Let `a_n>=1` and suppose `r>=a_nT_n`.  With

```math
L=\frac{a_n^2T_n^2}{H}
=\Theta(a_n^2n^{9/4-c}),
```

one has `L/(mk)->infinity`, uniformly on a compact fixed-density window.
Using `q_n=O(n^(3/2))`, the two exponents in (R55.1) satisfy

```math
\frac{L}{q_n}=\Omega(a_n^2n^{3/4-c})
=\Omega(a_n^2H),
```

and

```math
\frac{L^2}{q_nmk}
=\Omega(a_n^4n^{1-2c})
=\Omega(a_n^2H)\,a_n^2n^{1/4-c}.
```

Consequently the **Verified uniform mass theorem** is

```math
\boxed{
\nu_m\{g<0,\ r\ge a_nT_n,\ K_H(u)\ge r+b_H\}
\le2e^{-c_1a_n^2H}.
}
\tag{R55.3}
```

The selector average and orientation bit do not change the bound.

For every `omega_n->infinity`, (R55.3) is `e^{-omega(H)}`.  Hence it
contradicts the proposed lower bound `e^{-C_0H}` in (10.1277).  The scalar
class (10.1285) is a subset of the same exact K-success event by (10.1284),
so it is also impossible.

## 3. Consequence for the orientation-relaxed lemma

The conditional implication (10.1297) is algebraically correct, but its K
branch has mass at most `2e^{-c omega_n^2H}` by (R55.3).  Therefore, if the
union in (10.1297) had mass `e^{-C_0H}`, then for all sufficiently large
`n` essentially all of that saved mass would lie in the opposite-orientation
local-band branch.  That branch already proves the recurrence directly.

There is no intermediate asymptotic scale left for this implementation:

- if `r/T_n=O(1)`, saved local-state mass with `g=-p_2r>=-O(T_n)` proves
  the restriction recurrence by the Wave 51 inverse Hanson--Wright argument;
- if `r/T_n->infinity`, (R55.3) makes the necessary energy-qualified state
  population `e^{-omega(H)}`.

Thus the exact completion lower tail (10.1276) remains a valid conditional
probability statement, but cannot provide the missing annealed saved mass in
the proposed far-state framework.

## 4. Fixed-depth truncation makes every entropy-H K population local

The obstruction does not depend on having chosen a diverging threshold in
advance.  For every fixed `D>0`, the same proof, now with

```math
L_D=\frac{D^2T_n^2}{H},
```

gives, uniformly for all sufficiently large `n`,

```math
\boxed{
\nu_m\{g<0,\ r\ge DT_n,\ K_H(u)\ge r+b_H\}
\le 2e^{-c_2D^2H}.
}
\tag{R55.4}
```

Indeed the two Hanson--Wright exponents are respectively at least constant
multiples of `D^4 n^(1-2c)` and `D^2 H`; their ratio is
`D^2 n^(1/4-c)`, so the latter controls the minimum for fixed `D` and large
`n`.

Now suppose, for some fixed `C_0`, that the entire K-success population has
mass at least `e^{-C_0H}`.  Choose a fixed `D` so large that
`c_2D^2>C_0+1`.  Equation (R55.4) leaves mass at least
`(1/2)e^{-C_0H}` on states with

```math
-p_2DT_n\le g<0.
```

That is already a saved local-band population.  Conditioning on a selector
and applying the established inverse Hanson--Wright argument proves
`q_m<=p^(3/2)q_n+O(T_n)` directly, without using the completion theorem.
Consequently:

```math
\boxed{
\text{Every }e^{-O(H)}\text{-abundant K-success population is
recurrence-local after a fixed-depth truncation.}
}
\tag{R55.5}
```

This retires the K-functional condition not merely at the particular
diverging threshold in (10.1277), but as any independent entropy-`H`
annealed implementation based on uniform local states.  Its conditional
pointwise completion theorem can still be used inside a future argument whose
averaging law or entropy bookkeeping is genuinely different.

## Scope and successor

These theorems retire the orientation/K abundance implementation, not the
bare arbitrary-cut incidence `Z_t`, the box/spectral route, or a completion
mechanism that does not first select uniform local states through the event
`G>=r^2/H`.

A successor would have to use a different, selected local-state law whose
entropy cost is explicitly repaid, or return to the exact-minimizer box and
weighted spectral alternatives.  Merely weakening the clipped profile cannot
evade (R55.3), because its Euclidean necessary condition is already the wall.
Section 5 below shows that averaging the completion first under the same
uniform local-state law does not evade localization either.

## 5. Direct joint Hanson--Wright closes the full uniform-completion detour

There is a stronger argument which does average over the completion before
imposing a cross-energy event.  For fixed `S`, let `P_S` denote coordinate
projection, put

```math
C_S=A-P_SAP_S=
\begin{pmatrix}0&A[S,T]\\A[T,S]&A[T]\end{pmatrix},
```

and write `x=(y,w)`.  The centered completion polynomial is exactly

```math
Z(w)=\sigma x^{\mathsf T}C_Sx.
```

For either orientation, a uniform projective `y` followed by uniform `w`
has the same law for this globally sign-invariant quadratic form as a full
independent Rademacher vector.  Moreover

```math
\operatorname{tr}C_S=0,qquad
\lVert C_S\rVert_F^2\le n^2,qquad
\lVert C_S\rVert_{\rm op}
\le\lVert A\rVert_{\rm op}+\lVert P_SAP_S\rVert_{\rm op}
\le2\sqrt{2q_n}.
```

Hanson--Wright therefore gives the **Verified joint completion upper tail**

```math
\boxed{
\Pr_{a,w}\{Z(w)\le-u\}
\le2\exp\left[-c_3\min\left\{
\frac{u^2}{n^2},\frac{u}{\sqrt{q_n}}
\right\}\right].
}
\tag{R55.6}
```

At `u=a_nT_n`, `a_n>=1`, the two exponents are respectively
`Omega(a_n^2n^(1-2c))` and `Omega(a_nH)`.  Their ratio tends to infinity,
so

```math
\Pr_{a,w}\{Z(w)\le-a_nT_n\}
\le2e^{-c_4a_nH}.
\tag{R55.7}
```

Now use the exact normal form (10.1244), and denote its selector-averaged
incidence by `Z_t` as in the ledger:

```math
Z_t=\mathbb E_{S,a}\Pr_w\{Z(w)\le g/p_2\}.
```

For every fixed `D>0`, split according to the local margin.  On
`g<-p_2DT_n`, a favorable completion requires `Z(w)<-DT_n`, so (R55.7),
uniformly in `S`, proves

```math
\boxed{
Z_t\le\nu_m\{g\ge-p_2DT_n\}+2e^{-c_4DH}.
}
\tag{R55.8}
```

If `Z_t>=e^{-C_0H}`, choose fixed `D` with `c_4D>C_0+1`.  Then the first
term in (R55.8) is at least `(1/2)e^{-C_0H}` for all large `n`.  This is
again saved local-band mass and proves the restriction recurrence directly
by inverse Hanson--Wright.  Hence the full uniform-local-state annealed
completion incidence is not an independent far-tail mechanism either:

```math
\boxed{
Z_t\ge e^{-O(H)}\quad\Longrightarrow\quad
\text{saved }g\ge-O(T_n)\text{ mass}\quad\Longrightarrow\quad
q_m\le p^{3/2}q_n+O(T_n).
}
\tag{R55.9}
```

This does **not** retire the bare fixed-cut selector tail (10.795): there the
global cut is selected rather than averaged uniformly, and its restrictions
do not have the product local-state law used in (R55.6).  Nor does it retire
box/coarea extraction of such a selected low-row cut.  It retires the direct
uniform-cut completion detour as a genuinely nonlocal implementation.
