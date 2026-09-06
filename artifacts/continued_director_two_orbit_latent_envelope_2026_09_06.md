# The ternary latent lower envelope needs only two posterior orbits

Date: 2026-09-06. Exact finite-dimensional reduction of a LOWER bound
on the deep Bellman value, not an upper certificate or a new cap.

Use the latent theorem from
`continued_convergence_deep_latent_envelope_2026_09_06.md`. Write

```math
 \nu_p=(1-p)\delta_0+\frac p2\delta_{-1/\sqrt p}
                       +\frac p2\delta_{1/\sqrt p},
 \phi_G(u)=\sup_{0\le r<1}
        \{-u(1-r)+\tfrac14\log(1-r^2)\}.                (1)
```

The last expression is the exact Gaussian terminal potential, with
r/(1-r^2)=2u at u>0. Its convexity in u is important below.
The channel lower envelope relevant to the weave exponent is

```math
 E(p,t)=p\log2+t(1-\sqrt p)
       +\sup_L\{\phi_G(t\,E\operatorname{Var}(X|L))-I(X;L)\},
 \qquad X\sim\nu_p.                                    (2)
```

It lower-bounds the limiting Bellman criterion. A positive value would
exclude every finite depth for this p,t. A negative value does NOT
prove that the Bellman criterion is negative.

## 1. Posterior coordinates and exact channel realization

A posterior distribution on the three input values is described by

```math
 z=\Pr(X\ne0|L),\qquad b=\Pr(X=1/\sqrt p|L)
                         -\Pr(X=-1/\sqrt p|L),
 \quad 0\le z\le1,\quad |b|\le z.
```

Its entropy is

```math
 h(z)+z h\!\left(\frac{1+b/z}{2}\right),                 (3)
```

with the continuous convention at z=0. Pair every posterior with its
reflection b->-b. This preserves entropy and b² and enforces the zero
signed marginal. The only remaining input-marginal condition is Ez=p.
Every such posterior mixture is an actual channel: define the joint
law as mixture weight times posterior probability and then divide by
the positive input probabilities to obtain Pr(L|X).

Writing b=zs, 0<=s<=1 after reflection,

```math
 V=E\operatorname{Var}(X|L)=1-E(z^2s^2)/p,
 \quad I(X;L)=h(p)+p\log2-E[h(z)+z h((1+s)/2)].           (4)
```

These formulas include zero-mass posteriors and fully revealing labels.

## 2. Two posterior reflection-orbits suffice

The supremum in (2) is attained by a mixture of at most TWO posterior
reflection-orbits. Thus at most four latent labels suffice, with opposite
labels combined when their posterior is already symmetric.

To prove this, maximize

```math
 E[h(z)+z h((1+s)/2)]+\phi_G(t[1-E z^2s^2/p])             (5)
```

over probability measures on the compact square (z,s), subject to Ez=p.
The integrands extend continuously to its boundary. The objective is
continuous and CONVEX in the mixing measure, because phi_G is convex.
Any measure's three moments (Ez,E z²s²,E entropy) are reproduced by a
measure on at most four points, by finite-dimensional convex hull
reduction. In the polytope of weights on those points with Ez=p, each
extreme point has at most two nonzero weights. Decompose the original
weight vector into such extreme points; convexity implies at least
one of them has objective no smaller. Compactness then gives attainment.

Explicitly, choose z1<=p<=z2 and s1,s2 in [0,1], with

```math
 a=\frac{z_2-p}{z_2-z_1},\qquad
 V=1-\frac{a z_1^2s_1^2+(1-a)z_2^2s_2^2}{p}.            (6)
```

Substitute (6) and the weighted posterior entropies in (4). The case
z1=z2=p is the one-orbit limit. This is an exact four-real-parameter
description of (2), not an assumption about a binary sign channel.
Fixing V separately would generally require three posterior orbits;
the two-orbit conclusion uses maximization of the FULL convex Gaussian
potential, not a fixed-distortion rate function.

## 3. Equivalent one-dimensional entropy-envelope formula

There is a second useful exact reduction. For kappa>=0 define

```math
 g_\kappa(z)=h(z)+z\max_{0\le s\le1}
       \{h((1+s)/2)+\kappa z s^2\}.
```

Let cav(g_kappa) be its least concave majorant on [0,1]. Then

```math
 E(p,t)=-h(p)+t(1-\sqrt p)
 +\sup_{0\le r<1}\left\{-t(1-r)+\tfrac14\log(1-r^2)
             +\operatorname{cav}(g_{t(1-r)/p})(p)\right\}.       (7)
```

Indeed expose (1) inside (5), interchange two suprema, and fix r.
The remaining objective is linear in the posterior mixing measure;
maximize s pointwise and mix z values with mean p. This is exactly
the concave-envelope value. No minimax interchange of opposing extrema
is being used.

The inner maximization is elementary. Put a=kappa z. Its derivative is
`-atanh(s)+2as`. If 2a<=1, the unique maximizer is s=0. If 2a>1, the
unique positive maximizer solves

```math
 \operatorname{atanh}s=2as.
```

The ratio atanh(s)/s is strictly increasing from1 to infinity, as seen
from its positive power series. The derivative changes from positive
to negative at that root. Thus (7) is a scalar maximization and a
one-dimensional concave-envelope operation, with no hidden growing
alphabet. It is structurally the familiar binary mean-field entropy
maximization inside a ternary density constraint; no external phase-
diagram theorem or analogy is needed for the identity.

## 4. Current numerical role and limits

`computations/continued_director_two_orbit_latent_envelope_2026_09_06.py`
tests actual two-orbit channels by numerical optimization. Every
feasible channel gives only a LOWER witness for (2). Small floating
negative mutual informations at an uninformative channel are roundoff,
not negative information. Neither agreement of random starts nor a
negative selected-channel value is an upper certificate for (2), much
less for the Bellman maximum.

The point of the reduction is to test an infinite-depth obstruction
without assuming that its only possibilities are a Gaussian bulk and
complete entropy-paying condensation. Those two cases do not exhaust
general finite input laws, as the companion four-atom counterexample
proves. For the actual ternary root, its exact envelope is now a small
and independently falsifiable variational calculation.
