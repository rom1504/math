# A direct-E all-order upper bound below 0.493609

Date: 2026-09-07. Status: proved by analytic reduction and directed interval
coverage; pending independent audit of the new mean-value rectangle bound.

## Statement and scope

For the original hollow symmetric full-sign problem,

```math
\limsup_{n\to\infty}{M_n\over n^{3/2}}
\le {\frac{97}{20}+\frac{24}{25}\log2-\frac{5151}{6250}
          \over 2\frac{97}{20}\sqrt{\frac{24}{25}}}
<0.493609.
```

The middle expression is about 0.49360809358874863. This improves only
the upper endpoint. It neither establishes a limit nor identifies the
original optimum with the ternary construction's variational value.

## 1. Exact scalar reduction

Let p=24/25, t=97/20 and

```math
\nu_p=(1-p)\delta_0+{p\over2}
 (\delta_{-1/\sqrt p}+\delta_{1/\sqrt p}),\qquad
c_t(\lambda)={1\over4}\log{\lambda(2t-\lambda)\over t^2}.
```

The previously proved direct-E realization gives, for every fixed p,t,

```math
\limsup_n {M_n\over n^{3/2}}
\le {t+p\log2+E_t(\nu_p)\over2t\sqrt p}.                 (1)
```

Its sufficient dependencies are reconstructed in
`decisive_audit_standalone_direct_E_upper_2026_09_07.md`, Sections 1--3.
The closure and orbital estimates are for arbitrary fixed t>0; only the
old final numerical substitution specialized them to t=4. Fixed recursion
depth precedes all-order limits. The same terminal H2/H12 orders and
principal restriction fill all orders for the present rational p. No
new all-order sampling or depth-uniform theorem is needed.

The exact ternary source reduction, proved in Sections 1--2 of
`decisive_bridge_exact_ternary_gaussian_phase_2026_09_07.md`, is valid for
every 0<p<1 and t>0, not just its title's old parameter pair. It gives

```math
E_t(\nu_p)=\sup_{0<\lambda\le t,\,0\le z\le1}
\left[c_t(\lambda)-\lambda+
 \max_{0\le v\le1}\{p\log(1+vB)-\log(1+vA)\}\right],    (2)
```

where A=exp(lambda z^2/p)-1, B=cosh(2lambda z/p)-1. The optimal reproduction
measure has support {0,+a,-a}; the source KKT witness's positive local
maximum is unique because the coefficient ratios in its critical-point
quotient are increasing. The source endpoints are not missing reproduction
maxima. The change z=a sqrt(p) gives the displayed compact range.

For A,B>0 the maximizing weight is exactly

```math
v_* = \operatorname{clip}_{[0,1]}
       {pB-A\over(1-p)AB}.                                (3)
```

At a zero denominator the original continuous maximization is used.
The three gain branches are 0, p log(1+B)-log(1+A), and, in the interior,

```math
p\log p+(1-p)\log(1-p)-p\log A+\log B-(1-p)\log(B-A).     (4)
```

This reconstruction does not assume that there is just one local maximum
of (2) in the joint precision/reproduction plane. In fact the numerical
pilot found two distinct informative local maxima at the present parameters.

## 2. Low precision is controlled analytically

Since p>=1/3, the normalized ternary source is unit subgaussian: its even
moment p^(1-j) is at most (2j-1)!!, by induction. Entropy duality implies
I(X;L)>=E[E(X|L)^2]/2. Thus J_lambda(nu_p)=lambda for lambda<=1/2.

The maximum of c_t(lambda)-lambda over this interval is

```math
g_t(1)=-t(1-\rho)+{1\over4}\log(1-\rho^2),\qquad
\rho={\sqrt{1+16t^2}-1\over4t}.                           (5)
```

Its maximizer lambda=t(1-rho) lies in (0,1/4), so no endpoint is lost.
Directed interval evaluation in the verifier proves g_t(1)<-5151/6250.
Numerically g_t(1) is about -0.8244733448121. Unlike the old p=31/32,t=4
point, the present point is NOT asserted to be in Gaussian phase.

## 3. High-precision coverage and the new gradient bound

Put F(lambda,z)=c_t(lambda)-lambda+G(A,B), with G the optimized gain.
The closed rectangle [1/2,t] times [0,1] is covered by rational dyadic
subdivisions. Every accepted rectangle receives one of two rigorous bounds.

First, G increases with B and decreases with A. Moreover

```math
{d\over d\lambda}(c_t(\lambda)-\lambda)
={1\over4}\left({1\over\lambda}-{1\over2t-\lambda}\right)-1
\le -{1\over2}\quad(1/2\le\lambda\le t).
```

Consequently a rectangle [l,u] times [a,b] has the monotone upper bound

```math
c_t(l)-l+G(exp(l a^2/p)-1,\ cosh(2u b/p)-1).              (6)
```

Second, if a>0, A and B are positive throughout the rectangle. The maximizer
(3) is unique, including at its clipping transitions. Envelope
differentiation is valid, and the first derivatives agree at those
transitions because the v-derivative is zero there. Write
q=2lambda z/p and R=exp(lambda z^2/p). Then

```math
F_\lambda={1\over4}(1/\lambda-1/(2t-\lambda))-1
 +{2v_* z\sinh q\over1+v_*B}
 -{v_* z^2 R\over p(1+v_*A)},                             (7)
F_z={2v_*\lambda\sinh q\over1+v_*B}
 -{2v_*\lambda zR\over p(1+v_*A)}.                       (8)
```

The interval evaluation of (3), clipped monotonically to [0,1], encloses
every optimizing weight on the rectangle. Substitution into (7)--(8)
therefore gives rigorous derivative intervals. If their maximum absolute
endpoints are L_lambda,L_z, and (lambda_0,z_0) is the rectangle center,
the fundamental theorem of calculus on two axis-parallel segments gives

```math
\sup F\le F(\lambda_0,z_0)
       +L_\lambda(u-l)/2+L_z(b-a)/2.                     (9)
```

There is no convexity, Hessian-sign, or numerical stationary-point premise
in (9). Rectangles meeting z=0 use only (6). The point-center evaluation
is itself an upper bound obtained by downward rounding A, upward rounding
B, and selecting the gain branch with exact rational comparisons.

## 4. Complete certificate

Run from the repository root:

```bash
.venv/bin/python -B computations/flatify_independent_2026_09_07_ternary_interval.py \
  --output computations/results/flatify_independent_2026_09_07_ternary_interval.json
```

The default p,t,target are exactly 24/25,97/20,-5151/6250. The completed
50-decimal-digit directed interval run checked 86,041 tree nodes and
accepted 43,021 leaves, 42,895 using the new mean-value bound. Every leaf's
upper bound is strictly below the target. The split tree exhausts the
initial closed rectangle; the program never prunes using floating point.
The reported run took about 170 seconds.

A separate `--monotone-only` run stopped at its declared 2,000,000-node
budget after about 1188 seconds, with 19 rectangles still pending. Its
record is `computations/results/flatify_independent_2026_09_07_ternary_monotone_interval.json`.
It is explicitly NOT a second certificate. The completed certificate above
retains the independently audited mean-value bound; the unfinished run is
preserved as a failed bounded diagnostic, not silently discarded.

All endpoint arithmetic and branch comparisons use Python Fractions
extracted exactly from mpmath interval endpoints. All exp/log/sqrt calls
are directed interval calls. Rectangular coordinates are exact rational
numbers throughout. The floating cap in the JSON is explicitly display-only.
Its exact rational upper endpoint is

```math
{738722853867255520612371821466032639475585245768797
 \over1496577676626844588240573268701473812127674924007424}
< {493609\over10^6}.                                     (10)
```

Equations (2), (5), and the high-precision cover prove E_t(nu_p)<=-5151/6250.
Substitution in (1) and (10) proves the stated original upper bound.

## 5. Limits of this improvement

The scalar recursion is closed under its information-theoretic Bellman
operation, but this does not make it closed under arbitrary original
minimizing seeds. The old seed-loss artifact proves that arbitrary
independent Hadamard-column signs erase the off-diagonal weave seed:
each edge sign is absorbed by a separate pair of column signs.
Consequently this particular scalar/orbit construction cannot establish
liminf-to-all-order value transfer merely by inserting a better seed.
A genuinely correlated seed-dependent source class would be needed, and
the current E proof does not supply its closure. This is a scoped limitation
of the ensemble, not a nonconvergence theorem.
