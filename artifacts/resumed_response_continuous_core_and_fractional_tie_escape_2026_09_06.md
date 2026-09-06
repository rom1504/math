# Continuous cyclic cores and an exact quadratic escape from fractional ties

Date: 2026-09-06. Status: Sections 1--3 and the quantitative tie-escape
formulas independently reconstructed by the bound-audit agent. Section 4
was proposed by that agent and independently reconstructed here.
This note strengthens the analytic-core theorem in two different ways.
Neither assertion assumes that a large derivative, by itself, prevents
same-space Gaussian creation.

## 1. A continuous core suffices for genuinely ternary stationary states

Use the finite-frame setup of
`resumed_response_analytic_core_causal_extension_stationarity_2026_09_06.md`,
but replace global real analyticity of the core inverse features by
global CONTINUITY. Appended inverse features remain arbitrary measurable
strictly causal functions of earlier Gaussian coordinates. Suppose now
that the actual response is ternary:

    q=|f| in {0,1} almost surely.

Then every unrestricted stationary point still has

    J<=phi(0)=1/sqrt(2 pi).

The positive-value nonzero-gradient branch is impossible. The proof
eliminating all active appended a and A coefficients uses no regularity
of the core, so it applies verbatim. Thus, if V!=0, the functions

    K=sum_(j<=m) a_j h_j,
    w=sum_(j<=m) A_j h_j,
    R=w+2BK=H s(K)

are continuous functions of the core coordinates x in R^m, while
V=A dot x is a nonzero linear function. Set

    b_0=Psi(0,t)-B,   D(x)=|V(x)|+B-Psi(K(x),t).

On the open set K!=0, the quotient H=R/s(K) has a continuous
representative taking values in {0,1}: the almost-sure binary identity
forces this everywhere there, by continuity and the full support of
Gaussian measure. Hence H, and q, are locally constant on this set.

At a point with K=0 and |V|!=b_0, the score D is nonzero. Continuity
makes its sign fixed on a neighborhood, and the strict stationary
threshold rule makes q constant on that neighborhood as well. These
local representatives agree on overlaps because they agree almost
surely. Therefore q has a locally constant representative on

    R^m minus {|V|=b_0}.                           (1)

If b_0<0, the removed set is empty, so q is constant almost surely;
p>0 would imply q=1 and J=0. If b_0=0, the two halfspaces have the
same constant value by global evenness of q, yielding the same
contradiction. Thus b_0>0.

The three connected regions in (1) are the central slab |V|<b_0 and
the two tails. The central slab has D<0, since Psi(K,t)>=Psi(0,t),
so q=0 there. The two tail values agree by parity, and p>0 forces
both to equal one. Consequently

    q=1{|V|>b_0},   f=sign(V)q.                   (2)

Gaussian regression gives

    a=C A,   C=2phi(b_0/|A|)/|A|>0,
    K=Cw.

Hence the stationary inverse identity becomes

    H s(K)=(1/C+2B)K.                             (3)

On the tails H=0, so K=0. On the center H=1, the strict decrease of
s(k)/k for k>0 shows that the solutions of (3) are either just zero
or the finite set {0,k_*,-k_*}. Thus the continuous function K maps
the connected space R^m into a finite set and is constant. Its zero
tail value makes it identically zero, contradicting V!=0. The only
remaining stationary branch has V=0, whose phi(0) bound needs no
regularity or causality assumptions.

The binary hypothesis matters. A continuous fractional q can vary
within a positive-measure exact tie region. Analytic continuation
excluded such tie regions in the preceding analytic-core theorem;
mere continuity does not.

## 2. Fractional stationary ties have a nearby improving perturbation

The following result applies to ANY actual finite creation-closed frame,
without a continuous or causal core assumption. It does not say that
fractional stationary points fail the first-order test; they can pass
it. It says that a nonzero-gradient fractional stationary point cannot
be a local maximum, and gives an explicit quadratic escape.

Suppose (f,q) is stationary, V!=0, and

    kappa=P(0<q<1)>0.

Then stationarity gives f=sign(V)q almost surely. Fractional q can
occur only on the tie set D=|V|+B-Psi(K,t)=0. On K=0 that would be
the fixed level event |V|=Psi(0,t)-B, which has Gaussian probability
zero. Thus

    beta^2=E[q(1-q)s(K)^2]>0.                     (4)

Adjoin one independent standard Gaussian gate coordinate whose inverse
feature is old-measurable. It exists by choosing an even old L2 function
outside the finite inverse span, projecting off that span, normalizing,
and applying U. Let U_gate=2Phi(|G_gate|)-1, independent uniform on (0,1),
and set

    q'=1{U_gate<q},   f'=sign(V)q'.

These are even/odd, and their conditional old expectations are q and f.
Their first-chaos coefficients on the old frame are unchanged. Their
coefficient on the gate coordinate is zero by evenness in that coordinate.
No other coefficient is created because the functions are measurable in
the enlarged finite Gaussian frame. Hence p,a,K,t and J are unchanged.

For 0<theta<1 form the feasible flat segment

    q_theta=(1-theta)q+theta q',
    f_theta=(1-theta)f+theta f'.                    (5)

All its values have EXACTLY the same objective J, not merely zero
initial derivative. B is also unchanged. Its inverse full gradient is

    w_theta=w-theta r,
    r=(q'-q)s(K).

The function r has conditional old mean zero and squared norm beta^2.
It is orthogonal to every old inverse feature and to the new gate
inverse, since those are all old-measurable. Therefore U r/beta is
a standard Gaussian independent of the ENTIRE enlarged old frame.
The gradient at (5), conditional on that frame, is distributed as

    V_theta=V+theta beta Z,

with Z independent standard Gaussian. This is an exact creation-space
innovation identity, not a claimed conditional law of a matrix response.

All subtraction terms in the full gap remain unchanged along (5):
a,A,B,p and E q_theta Psi equal their original values. Since the original
gap is zero, the new gap is the increase of the convex Gaussian-averaged
positive-part term. Jensen makes the contribution off the old tie set
nonnegative. On the fractional tie set, D=0 and V!=0 almost surely, and

    (|V+sigma Z|-|V|)_+
      >=(sigma sign(V)Z)_+.

Its Gaussian expectation is at least sigma phi(0). Taking
sigma=theta beta gives the explicit bound

    gap(f_theta,q_theta)>=theta beta phi(0) kappa. (6)

There is a uniform Hessian bound M=M(t)<infinity for the line toward
the new full linearized maximizer, restricted to a step at most 1/2.
For example use the M_0 displayed in Section 7 of the analytic-core
note with t_0=t: all moments along the flat segment are unchanged, and
the subsequent half-line has residual variance at least t/2. Put

    g_theta=theta beta phi(0) kappa.

For sufficiently small theta, choose the second line-search step
lambda=g_theta/M<=1/2. Taylor's inequality and (6) give a final feasible
pair with objective at least

    J+g_theta^2/(2M)
      =J+theta^2 beta^2 phi(0)^2 kappa^2/(2M).     (7)

Its L-infinity distance from the original pair is at most
theta+2lambda in the f coordinate and theta+lambda in q. Thus strictly
larger values occur arbitrarily close to the original stationary pair.
The perturbation uses at most two new independent Gaussian coordinates:
the gate, then the gradient innovation. Their inverse features are
successively old-measurable, so both additions are causal.

## 3. Consequence and precise limits

Combining Sections 1 and 2: in a finite continuous creation-closed cyclic
core followed by finitely many arbitrary measurable causal appendages,
NO feasible pair with J>phi(0) is a local maximum against arbitrary
finite Gaussian extensions. If its first-order gap is positive, ordinary
line search improves it. If the gap is zero, V=0 is ruled out at that
value; a binary q is ruled out by Section 1; and a fractional q has the
explicit nearby improvement (7).

For an actual purified ternary trajectory in this class, the stronger
conclusion holds that the full first-order gap itself is always positive
above phi(0). For arbitrary fractional pairs, a zero first-order gap is
not excluded by continuity alone, so the uniform first-order gap theorem
proved for analytic cores must NOT be carried over without a new argument.

This does not yet classify stationary ternary states in genuinely rough
cyclic cores. The exact matrix-height kernel criterion remains available
to test their realizability. No uniform gain over expanding frames, no
sharp variational value, and no convergence statement for the original
matrix problem follow from the local-escape results alone.

## 4. Uniform two-stage escape on every fixed continuous-core frame

The bound-audit agent identified a cancellation which strengthens the
continuous-core result to a uniform TWO-STAGE statement. Fix one frame
of the continuous-core/causal-extension class, and c>phi(0). On its
nonempty compact superlevel set J>=c, there is a constant

    delta_2(frame,c)>0                             (8)

such that an objective-preserving one-gate support purification gives a
full gradient gap at least delta_2. A subsequent uniformly controlled
line-search step therefore increases the original objective by a
positive amount depending only on this fixed frame and c.

This is not a uniform bound on the ORIGINAL first-order gap; fractional
ties can make that gap zero. It is also not uniform over changing frames.

### Exact one-gate transformation for an arbitrary feasible pair

The old pair need not be stationary. Write

    w=(1-q)s(K)-2BK,
    A_j=<w,h_j>,
    sigma^2=||w||_2^2-|A|^2,
    v=E[q(1-q)s(K)^2].

Choose a new normalized even old-measurable inverse h_gate orthogonal
to the finite old inverse span AND to w. There is always such a function:
the even L2 space of a nondegenerate finite Gaussian frame is infinite
dimensional. Then G_gate=U h_gate is independent of the old frame and
of the old full gradient. Its even uniform U_gate gives

    q'=1{U_gate<q},
    f'=(f/q)q', with f/q=0 when q=0.               (9)

This is feasible in the conditional domain |f'|<=q'<=1. It preserves
parity, conditional means f,q, all first-chaos moments, p,K,t,B, and
the objective J. The amplitude f' need not itself be ternary; the
conditional-ternary realization theorem can be applied afterward if
an actual ternary response is desired. One gate is sufficient for (9).

Put r=(q'-q)s(K). Then r has zero conditional old mean, is orthogonal
to all old inverse features and h_gate, has norm squared v, and is
orthogonal to w. The new inverse gradient is w'=w-r. Its old projection
is still A; its coefficient on G_gate is zero by the choice of h_gate.
Therefore its conditional independent Gaussian residual has variance

    rho^2=sigma^2+v.                               (10)

In particular the new full gap is EXACTLY

    gap_post=E_G,Z[(|A dot G+rho Z|+B-Psi(K,t))_+]
               -a dot A-Bp+E[q Psi(K,t)].          (11)

Although the chosen gate inverse can depend on the old pair, the scalar
formula (11) does not depend on that choice. The variance in (10) and
the conditional Gaussian representation are exact isometry statements.

### Strict positivity at every high-value old pair

If the old gap is positive, rho>=sigma and convex Gaussian averaging
give gap_post>=gap_old>0. If the old gap is zero, the old pair is fully
stationary. Since J>phi(0), its gradient is nonzero. By Section 1 it
cannot have binary q; by Section 2 its fractional mass gives v>0 and
the gate transformation has strictly positive gap. Thus gap_post>0
at EVERY pair in the fixed-frame superlevel set.

Equivalently, one can apply the continuous-core theorem directly after
(9): a zero post-gap would imply stationarity, whose nonzero gradient
forces f'=sign(V')q'. Since q' is binary, this would be a genuinely
ternary stationary pair in a causal extension of the continuous core,
which Section 1 forbids at this value.

### Weak-* continuity and compactness

Use the same compact weak-* feasible set as in Section 7 of the analytic
core note. The moments p,a,K,t and J are continuous, and J>=c supplies
a uniform t>=t_0>0. B and A are continuous as before. The separate
sigma^2 is only lower semicontinuous, and v need not be continuous.
Their SUM, however, is continuous because the q^2 terms cancel:

    rho^2
      =E[(1-q)s(K)^2]
         -4B E[(1-q)K s(K)]
         +4B^2 |a|^2-|A|^2.                     (12)

Each term on the right is weak-* continuous: K varies strongly in L2,
s(K) strongly in L2 and boundedly, and q is tested against fixed L1
functions after this strong convergence. The expression is nonnegative
by its exact variance interpretation (10), so rho is continuous too.

Formula (11) is consequently continuous on the compact superlevel set.
Its strict positivity just proved implies that it attains a positive
minimum, establishing (8). This avoids treating a weak-* limit of binary
support functions as still binary; the larger conditional domain is
retained throughout the compactness argument.

Finally use the uniform Hessian bound M_0 from the analytic-core note,
with the same t_0. The gate step leaves the value and moments unchanged.
The line-search step

    lambda=min(1/2,delta_2/M_0)

then improves J by at least lambda delta_2/2. All added inverse features
are successively old-measurable; each fixed finite stage remains a valid
finite response construction. The uniform gain belongs to the specified
old frame. Repeating the process changes that frame, so (8) does not
by itself establish a global iteration rate or identify its limit.
