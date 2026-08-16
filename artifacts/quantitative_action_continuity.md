# Quantitative action continuity under a spectral bound

## Status

**Verified supporting theorem.**  This sharpens the qualitative continuity
statement in `concentration_compactness_boolean_profiles.md`.  It does not
provide the missing all-order realization of signed action limits.

Let `S,T` be `P`-operators on probability spaces, and write

```math
\Phi(T)=\sup_{|f|\le1}|\langle f,Tf\rangle|.
```

For a hollow finite matrix `A` of order `n`, viewed as
`T_A=A/sqrt(n)` on the uniform probability space, coordinatewise affinity
gives

```math
\Phi(T_A)={1\over n^{3/2}}\max_{x\in\{\pm1\}^n}|x^TAx|
={2Q(A)\over n^{3/2}}.
```

Let `S_1(T)` be the one-profile: the set of laws of `(f,Tf)` over all
measurable `|f|<=1`.  Closures are taken in the weak topology and distances
between probability laws use the Levy--Prokhorov metric.
The profile and action metric are those of
[Backhausz--Szegedy, *Action convergence of operators and graphs*](https://arxiv.org/abs/1811.00626).

## Theorem

Suppose

```math
\|S\|_{2\mathbin\to2},\|T\|_{2\mathbin\to2}\le C
```

and put

```math
\delta=d_H^{\rm LP}\bigl(\overline{S_1(S)},
                          \overline{S_1(T)}\bigr).
```

For `0<delta<=1`,

```math
\boxed{|\Phi(S)-\Phi(T)|\le5C\sqrt\delta+\delta.}       \tag{1}
```

For the standard action metric `d_M`, whose `k=1` profile term has weight
`1/2`, this implies

```math
\boxed{
|\Phi(S)-\Phi(T)|
\le5C\sqrt{2d_M(S,T)}+2d_M(S,T).}                    \tag{2}
```

## Proof

Take laws `mu,nu` in the two closed one-profiles with
`d_LP(mu,nu)<=delta`.  Write `(X,Y)~mu` and `(X',Y')~nu`.  Every law in a
one-profile has `|X|<=1` and

```math
\mathbb E Y^2\le C^2;
```

the latter inequality passes to weak limits by lower semicontinuity.
Strassen's theorem supplies a coupling for which the Euclidean distance
between `(X,Y)` and `(X',Y')` exceeds `delta` with probability at most
`delta` (first use `delta+o(1)` if the infimum is not attained).

For `R>0`, truncate by

```math
\theta_R(y)=\max(-R,\min(y,R)),
\qquad g_R(x,y)=x\theta_R(y).
```

On the good coupling event, `g_R` changes by at most `(R+1)delta`; on
the exceptional event it changes by at most `2R`.  Also

```math
\mathbb E|XY-g_R(X,Y)|
\le\mathbb E[|Y|1_{\{|Y|>R\}}]
\le{C^2\over R},
```

and the same holds for `nu`.  Consequently

```math
\left|\int xy\,d\mu-\int xy\,d\nu\right|
\le {2C^2\over R}+(3R+1)\delta.                    \tag{3}
```

Taking `R=C/sqrt(delta)` yields (1); `C=0` is immediate.  Hausdorff
matching works in both directions, so taking the suprema of the absolute
energy integral proves the claim.  Finally
`d_H^{LP}(S_1(S),S_1(T))<=2d_M(S,T)`, which gives (2).

## Boundary

A common `2 -> 2` bound is essential.  Under a pure `infinity -> 1` bound,
an `n^(3/4)`-vertex Boolean spike can vanish in every fixed action profile
while carrying a constant normalized quadratic energy.  The theorem therefore
closes the continuity step only after spectral regularization.  It supplies
neither a fixed-bound `o(n^(3/2))` regularization theorem nor exact sign
realizers at every sufficiently large order.
