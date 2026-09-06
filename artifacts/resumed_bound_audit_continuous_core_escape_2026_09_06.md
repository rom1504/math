# Independent reconstruction: continuous cores and purification escape

Date: 2026-09-06. This note independently checks the analytic-core,
continuous-core and fractional-tie claims of the response agent. It
also records the compactness strengthening developed jointly after
those checks. All frames below are finite and actual in the fixed
canonical creation space. No uniformity over expanding frames follows.

## 1. Analytic cyclic core plus measurable causal appendages

Let G be a finite orthonormal Gaussian creation-closed frame. Its first
m inverse features are globally analytic functions of G_1,...,G_m;
each later inverse is measurable in the preceding coordinates. For
odd f, even q and |f|<=q<=1 use

    p=Eq, a=EfG, K=a dot h, t=p-|a|^2,
    J=E(1-q)Psi(K,t),
    s(K)=2Phi(K/sqrt(t))-1,
    B=E[(1-q)phi(K/sqrt(t))]/sqrt(t),
    w=(1-q)s(K)-2BK, V=Uw.

At positive value t,B>0. A zero unrestricted first-order gap forces
V=A dot G within the old frame and the exact threshold rule
f=sign(V)q, with q=1 on |V|+B>Psi and q=0 on the opposite strict set.

The largest-index argument in the triangular audit removes every
appended coefficient of both a and A when V is nonzero. It uses no
regularity of the core. Thus K,w are analytic core functions. On each
V-sign halfspace, the analytic score cannot vanish identically: analytic
continuation would equate the positive function Psi to an affine
function taking negative values. Ties therefore have zero measure.
Positive p supplies a nonempty open strict-support set, where
(1-q)s(K)=w+2BK=0. Analytic continuation makes this zero globally,
contradicting EfV=-2B|a|^2<0 and EfV=Eq|V|>=0.

The only remaining branch is V=0. The distributional calculation in
the earlier triangular audit, which uses no causality, gives
J<=phi(0)=1/sqrt(2pi). Thus all higher-value pairs in this analytic-core
class have a positive unrestricted first-order gap. The full proof in
`resumed_response_analytic_core_causal_extension_stationarity_2026_09_06.md`
passes independently, including its fixed-frame weak-* compactness
upgrade to a uniform first-order gap.

## 2. Continuity suffices when the response is ternary

Now assume only that the core inverse functions are globally continuous,
and q=|f| is binary. After the same elimination of appended coefficients,
K,w and R=w+2BK are continuous core functions, while V=A dot G is a
nondegenerate Gaussian linear form. Set b_0=Psi(0,t)-B.

On K!=0, R/s(K) is continuous and equals the binary function H=1-q
almost surely. Full Gaussian support forces it to be binary everywhere
on that open set, hence locally constant. At K=0 away from |V|=b_0,
the strict score is nonzero and continuity makes q locally constant
on a neighborhood. These representatives agree on their overlaps.
Consequently q is constant on the connected components outside the
at most two hyperplanes |V|=b_0.

For b_0<=0, parity forces q constant almost surely, incompatible with
positive J and p. For b_0>0, the central slab has negative score since
Psi(K,t)>=Psi(0,t); hence q=0 there. Parity and p>0 give q=1 on both
outer components. Therefore a=C A with C>0 by Gaussian regression,
and K=Cw. The inverse identity reduces to

    (1/C+2B)K=H s(K).

Its possible K values are the finite set {0,k_*,-k_*}, since s(k)/k
is strictly decreasing for positive k. Connectedness and continuity
make K constant; the outer support has H=0 and hence K=0. This
contradicts nonzero V. The V=0 cap remains phi(0).

This reconstruction verifies all details of Section 1 of
`resumed_response_continuous_core_and_fractional_tie_escape_2026_09_06.md`.

## 3. Fractional stationary ties: flat motion then improvement

For an arbitrary actual finite closed frame, suppose a stationary pair
has V!=0 and positive mass of 0<q<1. Fractionality is confined to exact
ties. Its intersection with K=0 lies in the null fixed Gaussian level
|V|=Psi(0,t)-B. Hence

    beta^2=E[q(1-q)s(K)^2]>0.

Append an independent Gaussian gate W=U h_W, with h_W old-measurable
and orthogonal to the old inverse span. This is possible because the
even L2 space of a nontrivial finite Gaussian frame is infinite
dimensional. An even transform of W supplies independent uniform U_0.
Put q'=1{U_0<q} and f'=sign(V)q'. Their conditional old expectations
are q,f and their gate first-chaos coefficient is zero. Thus all of
p,a,K,t,B,J remain unchanged along the segment to this pair.

The new inverse gradient is w_theta=w-theta(q'-q)s(K). Its added term
is orthogonal to every old inverse and to h_W, because its conditional
old mean is zero. It therefore creates an independent Gaussian
innovation of norm theta beta. On the fractional tie set the increase
of the convex positive-part term is at least theta beta phi(0).
Thus the new full gap is at least

    theta beta phi(0) P(0<q<1).

A sufficiently small ordinary line-search step gives a positive
quadratic-in-theta objective increase at arbitrarily small distance
from the original pair. This verifies the response artifact's explicit
quadratic escape, but does not turn an originally zero first-order
gap into a positive first-order gap at that original point.

## 4. Uniform two-stage escape on a fixed continuous-core frame

There is nonetheless a fixed-frame compactness theorem. Fix c>phi(0)
and the old continuous-core/causal frame. On its feasible superlevel
set J>=c, use the product weak-* topology for (f,q). The set is compact,
and p,a,K,t,J are continuous; the Gaussian moment bound gives a uniform
t>0 floor. The quantities B and A_j=<w,h_j> are continuous, while
sigma^2=||w||^2-|A|^2 is lower semicontinuous.

For a GENERAL pair, choose h_W additionally orthogonal to w. Use the
lossless feasible purification

    q'=1{U_0<q},
    f'=(f/q)q', with ratio zero where q=0.

This preserves all moments and J and makes q' binary; it need not make
f' ternary. It is enough for the argument: if its full gap were zero,
stationarity would force |f'|=q' when V' is nonzero, and the continuous
core theorem would contradict J>c. The V'=0 branch is already capped.

More explicitly, put v=E[q(1-q)s(K)^2]. The post-purification full
gradient has old coefficients A, zero coefficient on the gate, and
independent residual variance

    rho^2=sigma^2+v.

This variance is CONTINUOUS under weak-* convergence, not merely lower
semicontinuous. Indeed the q^2 terms cancel exactly:

    ||w||^2+v
      =E[(1-q)s(K)^2]
         -4B E[(1-q)s(K)K]+4B^2||K||^2.

Each term on the right is continuous, as is |A|^2. The post-purification
full gap is therefore a continuous, strictly positive function on the
compact old-frame superlevel set. It has a positive minimum depending
on that frame and c. A uniform bounded-Hessian line search after the
lossless gate yields a uniform two-stage objective improvement there.

For any fixed partial purification fraction 0<theta<1, the residual
variance is sigma^2+theta^2 v, which is lower semicontinuous since it
equals theta^2 rho^2+(1-theta^2)sigma^2. The resulting gap remains
strictly positive: a previously positive gap cannot decrease under
independent Gaussian variance increase, while a zero old gap must be
fractional and has v>0. This yields the analogous uniform improvement
after any prescribed positive partial purification fraction.

These are FIXED-frame constants. The construction does not bound them
away from zero after successively adjoining new gradient coordinates,
and it supplies no convergence theorem for the original signing minima.
