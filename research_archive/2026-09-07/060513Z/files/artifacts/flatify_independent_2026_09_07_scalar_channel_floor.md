# Exact parameter mechanism and a scalar-criterion floor

Date: 2026-09-07. Status: analytic derivation; numerical values are labeled.
This is NOT a new lower bound for the original M_n.

Write

```math
C(p,t)={t+p\log2+E_t(\nu_p)\over2t\sqrt p},\qquad
0<p\le1,\quad t>0.
```

## 1. An explicit informative channel

For q in [0,1], let a nonzero input X reveal its sign as a label +/-.
When X=0, assign label 0 with probability 1-q, and each sign label with
probability q/2. Put r=p+(1-p)q, and let h denote binary entropy in nats.
Then

```math
I(X;L)=h(r)-(1-p)h(q)+p\log2,
\qquad E\operatorname{Var}(X\mid L)={r-p\over r}.
```

Indeed the sign labels each have mass r/2 and posterior mean +/-sqrt(p)/r;
the zero label has posterior mean zero. The mean posterior-mean square is
p/r. The conditional label entropy is (1-p)[h(q)+q log2], giving the
information formula. Thus the definition of E gives the exact bound

```math
E_t(\nu_p)\ge
g_t((r-p)/r)-h(r)+(1-p)h(q)-p\log2.                       (1)
```

At q=1 this is

```math
E_t(\nu_p)\ge g_t(1-p)-p\log2.                          (2)
```

At q=0 it is the full-revelation bound E_t>=-h(p)-p log2. Neither endpoint
requires solving the rate-distortion problem. Intermediate q expresses
partial revelation of the zero input and explains the second informative
branch in the numerical pilot. This interpretation does not assert exact
optimality of these channels: an optimal finite-temperature channel can
also make rare sign errors.

At p=.96,t=4.85, the q=1 lower witness gives E>=-0.8241705351707043
(numerical), within about 2.8e-6 of the first optimized informative branch.
The q-parametrized family is therefore a mechanism for the observed
transition, not a numerical optimizer-completeness argument.

## 2. Two exact scalar lower bounds

The uninformative channel gives E_t(nu_p)>=g_t(1). Minimizing its consequent
bound in t exactly yields

```math
C(p,t)\ge f(p):={\sqrt{1-2^{-4p}}\over2\sqrt p}.          (3)
```

For completeness, put rho=(sqrt(1+16t^2)-1)/(4t). Then
g_t(1)=-t(1-rho)+log(1-rho^2)/4, and differentiation of the Gaussian
criterion shows its unique minimum occurs at 1-rho^2=2^(-4p). The minimum
value is rho/(2sqrt(p)), which proves (3). This is the archived Gaussian
certificate floor, retained here with its p-dependence.

Independence is an admissible Gaussian self-coupling, so g_t(v)>=-tv.
Equation (2) consequently implies

```math
C(p,t)\ge {\sqrt p\over2}.                               (4)
```

There is separately an ACTUAL rank-one-weave floor sqrt(p)/2, proved by
matching-fibre spin states in `transfer_adversary_seed_loss_2026_09_06.md`.
Equation (4) is a direct scalar-channel proof. Equation (3) is only a
certificate floor; combining them must NOT promote (3) into an actual
weave lower bound.

## 3. A uniform exact combined barrier

The function f in (3) strictly decreases, since

```math
{d\over dp}f(p)^2
 ={(1+4p\log2)2^{-4p}-1\over4p^2}<0.
```

The right side of (4) strictly increases. They cross at the unique
positive solution p_* in (0,1) of

```math
p_*^2+2^{-4p_*}=1.
```

Uniqueness follows because the left side minus one is strictly convex,
vanishes at zero, has negative derivative there, and is positive at one.
Therefore, for every p,t,

```math
\boxed{C(p,t)\ge {\sqrt{p_*}\over2}.}                     (5)
```

Numerically p_*=0.9649459054049875... and the barrier is
0.4911583007048205..., improving the old uniform scalar floor
sqrt(15)/8=0.4841229182759271.... These decimals are illustrative; the
implicit exact value in (5) is the theorem.

The criterion has stronger constraints than the coarse barrier (5): (2)
retains g_t(1-p), and (1) retains all partial-zero channels. A numerical
minimization using only (2) and the Gaussian branch gave about
0.493595582557 at p about .959339,t about 4.90715. This last number is
not a certified global lower bound and is not used in any original result.

## 4. Original-limit scope

The newly certified original upper 0.493608093589... lies close to this
specific criterion's numerical floor. Further parameter precision would
not provide an all-seed closure. The archived column-sign gauge makes the
independently signed rank-one-weave ensemble seed blind; its scalar E
recursion is not a realization theorem for arbitrary original minimizers.
Thus neither this barrier nor the certified upper settles convergence.
