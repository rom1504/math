# Partial random flips do not close the actual-child carrier

Status: **rigorous equivalence and ceiling audit**.  Randomly flipping only
part of an optimizing child does give exact temperature-contraction and
one-vertex deletion identities.  The unrestricted family is equivalent to
the original sign-flip minimization; its exchangeable versions are radial;
and cut-correlated flips are universal gauge identities.  None determines
the norm-bounded low-degree product carrier, the sector--Gram direction, or
the target orientation without one additional nonradial closure theorem.

## 1. Arbitrary random flips are convex consequences

Fix an order-`m` child `A` which minimizes pressure at temperature `t`, and
write

```math
R_A(S)={Z_{A^S}(t)\over Z_A(t)}\ge1
\qquad(S\subseteq E(K_m)).                           \tag{PFR.1}
```

For every probability law `Pi` on edge subsets, independent of the child
spins,

```math
\boxed{\mathbb E_{S\sim\Pi}R_A(S)\ge1.}             \tag{PFR.2}
```

This contains no stronger inequality than (PFR.1): it is a convex
combination of the vertex inequalities, and allowing point masses `Pi`
recovers every one of them.  Thus an unrestricted random-flip theorem is
exactly full sign minimization in averaged notation.

For independent edge flips of probabilities `p_e`, FC.1 gives the more
useful temperature form

```math
\mathbb E_SR_A(S)
={P_A(((1-2p_e)\tanh t)_e)\over P_A(\tanh t\,\mathbf1)}\ge1.       \tag{PFR.3}
```

Because `P_A` is multiaffine, knowing (PFR.3) on the full probability cube
is again equivalent to knowing all vertex inequalities (PFR.1).  Any strict
reduction must restrict the probability family and prove that the omitted
directions are irrelevant; randomization alone supplies no such theorem.

## 2. Star contraction is exact one-vertex minimality

There is a useful nonradial specialization, but it exposes precisely the
already missing extension landscape.  Fix a vertex `v`, put

```math
L_v(x)=\sum_{j\ne v}a_{vj}x_j,
```

and let `nu_(A-v,t)` be the augmented Gibbs law of the deleted child.  Exact
summation over `x_v` gives

```math
{Z_A(t)\over Z_{A-v}(t)}
=\mathbb E_{\nu_{A-v,t}}\cosh(tL_v).                \tag{PFR.4}
```

Taking `p_e=1/2` on the star of `v` and zero elsewhere in (PFR.3) yields

```math
\boxed{
\mathbb E_{\nu_{A-v,t}}\cosh(tL_v)
\le(\cosh t)^{m-1}.}                                \tag{PFR.5}
```

More generally, replacing the incident row signs by any
`b in {+-1}^{m-1}` is a sign flip of `A`; hence exact child minimality says

```math
\boxed{
\mathbb E_{\nu_{A-v,t}}\cosh(t\langle a_v,X\rangle)
\le
\mathbb E_{\nu_{A-v,t}}\cosh(t\langle b,X\rangle)
\quad\hbox{for every }b.}                            \tag{PFR.6}
```

Thus the row of an actual minimizer solves the complete one-vertex
extension problem against its deleted-child Gibbs law.  Independent
partial flips within this star produce the inhomogeneous contraction of
the same extension response.  Since its boundary is the entire Boolean
table in (PFR.6), the full star box is not a lower-information statistic.

At physical scale `t=beta/sqrt(N)`, (PFR.5) gives only

```math
\mathbb E L_v^2
\le {2\{(\cosh t)^{m-1}-1\}\over t^2}=O_\beta(N).    \tag{PFR.7}
```

This controls one signed local-field quadratic form, not the Frobenius
mass of the child correlation matrix.  Expanding (PFR.6) to second order
would identify the quadratic extension form, but at physical scale the
higher even cumulants have no uniform negligible remainder.  Keeping all
of them reconstructs (PFR.6).  This is the same information-versus-tail
barrier as the bridge problem, now in a one-vertex slice.

## 3. Correlated cut flips are universal Ward identities

One might hope that correlated partial flips avoid the full edge cube.
The most natural family does not use minimality at all.  Given vertex signs
`sigma`, flip exactly the cut

```math
S_\sigma=\{ij:\sigma_i\sigma_j=-1\}.
```

Then

```math
A^{S_\sigma}=\operatorname {diag}(\sigma)
 A\operatorname {diag}(\sigma),
```

so the change of variables `x->diag(sigma)x` gives

```math
\boxed{R_A(S_\sigma)=1}                              \tag{PFR.8}
```

for every signing `A`, minimizing or not.  Averaging biased or correlated
vertex signs and differentiating their bias parameters therefore produces
only gauge/Ward identities.  Such identities cannot distinguish optimizing
children or constrain their retuning phase.

## 4. Orientation blindness and the radial ceiling

All child flip ratios are blind to global signing reversal:

```math
R_{-A}(S)=R_A(S),                                    \tag{PFR.9}
```

because `Z_{-C}(t)=Z_C(t)` after `tau->-tau`.  The two bridge orientations
replace one child by its global negative.  Therefore no statistic made only
from the random-flip ratios can select the bias-canceling orientation or
prove that it is target-relevant.

If the flip law is exchangeable over edges, (PFR.2) retains only the
fixed-size averages.  FC.18--FC.21 show that the complete collection of
these averages is exactly the radial absolute-energy histogram.  The two
certified order-eight actual minimizers in FC.22--FC.26 have the same such
histogram and every exchangeable partial-flip value, but different
zero-temperature two-replica correlation Frobenius masses (`14` versus
`10`, hence also different values at all sufficiently large finite
temperature) and different one-vertex response tables.  Thus exchangeable
partial flips determine
neither sector--Gram data nor even the smallest nontrivial extension
response.

## 5. Consequence for the low-degree carrier

The carrier in LDC.1 is determined by bounded-degree moments of the
*two-child optimal product factors*.  Equations (PFR.3) and (PFR.6) provide
no identity for those moments:

- the full inhomogeneous data are equivalent to all sign-flip or extension
  values;
- exchangeable data are falsified as sufficient by the actual order-eight
  pair;
- cut-correlated data are universal;
- all data are orientation-blind.

The only potentially useful intermediate statement left by partial flips
is quantitative deletion/extension rigidity, for example a theorem that
the gap in (PFR.5) is small enough to suppress the higher extension
cumulants uniformly at `t=beta/sqrt(N)`.  No such estimate follows from
minimality itself.  Without it, randomized partial-flip minimality neither
evaluates the carrier value `V^(d,K)` nor bounds the physical sector--Gram
direction, and it does not advance the target orientation question.
