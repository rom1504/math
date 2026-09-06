# No high-value stationary point after finitely many births from an analytic cyclic core

Date: 2026-09-06. Status: independently reconstructed by the bound-audit
agent; submitted to the director for a further reconstruction.
This extends the strictly causal theorem in
`resumed_bound_audit_triangular_stationarity_tie_obstruction_2026_09_06.md`.
The initial core may be genuinely cyclic and countably supported in the
canonical Gaussian coordinates. All appended inverse features may be
arbitrary measurable functions, including aligned discontinuous threshold
and softsign-threshold functions. They are NOT assumed polynomial,
analytic, Sobolev, or to have uncancellable jumps.

## 1. The class and the theorem

Let G=(G_1,...,G_d) be an actual finite standard creation-closed Gaussian
frame with orthonormal jointly even inverse features h_j=U*G_j. For some
1<=m<=d suppose:

* For j<=m, h_j is a globally real-analytic function of G_1,...,G_m.
  These core features may depend on their own coordinates and on each
  other; strict causality of the core is not required.
* For j>m, h_j is measurable in G_1,...,G_(j-1).

Let odd f and even q be measurable in this frame with |f|<=q<=1. Set

    p=E q, a_j=E f G_j, K=sum_j a_j h_j,
    t=p-|a|^2, H=1-q,
    J=E H Psi(K,t),   Psi(k,t)=E|k+sqrt(t)N|.

If this pair is stationary against ALL feasible first variations of
the purified full-response functional, including Gaussian coordinate
extensions, then

    J <= phi(0) = 1/sqrt(2 pi) < 0.4.             (1)

In fact, any positive-value stationary point in this class has identically
zero full Gaussian gradient V. The nonzero-gradient branch is impossible.
The zero-gradient branch really can occur below the bound, as the explicit
two-coordinate causal example in the cited note demonstrates.

Consequently every pair in this class above phi(0), in particular above
the present 0.43146 lower certificate, has a strictly improving feasible
first-order direction. The direction is the full-gradient threshold
followed by a sufficiently small convex line-search step. It adds at
most one independent Gaussian coordinate, with an old-measurable inverse,
and hence stays inside this same class after each finite number of steps.

This statement supplies no uniform step size or increment, and does not
identify the limiting value of the variational iteration. It does not
apply to arbitrary finite cyclic frames with nonanalytic core inverses,
nor to an infinite-depth limiting frame without further work.

## 2. Full stationarity and its pointwise rule

At positive J the Gaussian moment-body inequality gives p>0 and t>0.
Write

    s(k)=2Phi(k/sqrt(t))-1,
    B=E[H phi(K/sqrt(t))]/sqrt(t),
    w=H s(K)-2B K,   V=U w.

The unrestricted first derivative is

    DJ[delta f,delta q]=E[V delta f+(B-Psi(K,t))delta q].

Thus stationarity requires the pointwise maximizer

    f=sign(V)q,
    q=1 if |V|+B>Psi(K,t),
    q=0 if |V|+B<Psi(K,t),                         (2)

with ties allowed. If V has a nonzero independent Gaussian component
outside the old frame, conditioning on the old frame gives a strictly
positive improvement: the two-sided Gaussian tails make Jensen strict
for (|v|+B-Psi)_+. Therefore stationarity forces

    V=sum_j A_j G_j,   w=sum_j A_j h_j.             (3)

We first suppose V is not identically zero and derive a contradiction.
Since J>0, one has EH>0 and consequently B>0.

## 3. All active appended coefficients can be eliminated

If a=0, then K=0, w=0, and V=0, contrary to the current case. Let k be
the largest index with a_k!=0. We claim k<=m.

Suppose k>m. By the stated core/causal structure, K is measurable in
G_<k. Let ell be the largest index with A_ell!=0.

### Case ell>=k

Conditional on G_<k, V has nondegenerate Gaussian tail variance
sigma^2=sum_(j>=k) A_j^2, while the threshold Psi(K,t)-B is fixed.
Ties have zero conditional probability. Direct one-variable Gaussian
integration in (2) gives, for all j>=k,

    a_j=A_j C,                                    (4)

where C is the same strictly positive finite scalar. Explicitly, with
conditional mean v_0 and theta=max(Psi(K,t)-B,0), its integrand is

    [phi((theta-v_0)/sigma)+phi((theta+v_0)/sigma)]/sigma.

Since a_j=0 for j>k, (4) forces A_j=0 for j>k and A_k!=0. Thus both
w=sum_(j<=k) A_j h_j and K are G_<k-measurable. So Hs(K)=w+2BK is also
G_<k-measurable.

On an old event where K!=0 and Psi(K,t)>B, however, (2) makes H the
indicator of a finite nonempty interval in the nondegenerate Gaussian
G_k. Its conditional probability is strictly between zero and one.
It cannot have the required earlier measurability. Such old events
therefore have probability zero. On the complementary events either
s(K)=0 or H=0 almost surely, and hence Hs(K)=0 everywhere. It follows
that w=-2BK and

    E fV=<K,w>=-2B|a|^2<0,

contradicting E fV=E q|V|>=0.

### Case ell<k

Here w,K,V are all G_<k-measurable. On K!=0 the relation
Hs(K)=w+2BK determines H, hence q, from these earlier coordinates.
On K=0 the threshold rule determines q except on the event

    |V|=Psi(0,t)-B.

This is a fixed level set of the nondegenerate Gaussian V, and has
probability zero. Also V=0 has probability zero. Thus f=sign(V)q is
G_<k-measurable, which forces a_k=0, another contradiction. This proves
that a has no appended-coordinate component.

Now K is measurable in the core G_<=m. If A had any nonzero appended
component, condition on the core and repeat (4) for all j>m. It would
force a_j=A_j C, with C>0, contradicting that all these a_j are zero.
Therefore A also has no appended-coordinate component.

This elimination is the reason aligned jumps in different appended
inverse features do not defeat the theorem. At stationarity their
coefficients in BOTH K and w vanish. No smooth-complement assumption
or jump-rank calculation is used.

## 4. Global analyticity of the core excludes the remaining branch

We now have

    K=sum_(j<=m) a_j h_j,
    w=sum_(j<=m) A_j h_j,
    V=sum_(j<=m) A_j G_j != 0.

Thus K,w, and R=w+2BK are globally real-analytic functions on R^m.
The Gaussian-smoothed absolute value Psi(k,t) and s(k) are globally
real analytic in k because t>0. On each of the two halfspaces V>0
and V<0, the threshold score

    D=|V|+B-Psi(K,t)

is therefore real analytic. Neither sign-specific analytic expression
can vanish identically on a nonempty open set. For example, if
V+B-Psi(K,t) vanished there, the real-analytic identity theorem would
force Psi(K,t)=V+B everywhere on R^m. The right side is negative on
an open halfspace, whereas the left side is strictly positive. The
same argument applies to -V+B-Psi(K,t).

The zero set of a nonzero real-analytic function on a connected open
set has Lebesgue measure zero, hence Gaussian measure zero. One can
see this without an imported probabilistic theorem by partitioning
zeros according to the first nonvanishing partial derivative and
applying the implicit-function theorem to that derivative; analyticity
excludes points where all derivatives vanish unless the function is
identically zero. Consequently threshold ties in (2) have probability
zero.

Since p=Eq>0, the support set D>0 has positive measure and is a nonempty
open set in the core variables. On this set H=0 almost surely, so

    R=w+2BK=H s(K)=0

almost everywhere there. Continuity first makes R identically zero
on that open set, and global real analyticity then gives R=0 everywhere.
Thus w=-2BK. As in Section 3, this contradicts

    E fV=E q|V|>=0,   E fV=-2B|a|^2<0.

We have proved that a positive-value stationary point in the stated
core/extension class must have V=0.

## 5. The zero-gradient bound uses no causality hypothesis

For completeness, the phi(0) bound for V=0 is the distributional
calculation in Section 4 of the cited triangular-stationarity note;
none of its steps uses strict causality or analyticity. Here is the
reduction, making that scope explicit.

If a=0 then J=(1-p)sqrt(p)sqrt(2/pi)<=sqrt(8/(27pi))<phi(0).
Otherwise w=0 and the pointwise rule force B=Psi(0,t): if B is smaller,
q=0 everywhere; if B is larger, substituting Hs(K)=2BK into the defining
equation for B gives the impossible identity

    1=E_(K!=0)[2z phi(z)/(2Phi(z)-1)]<1,
    z=|K|/sqrt(t).

With B=Psi(0,t), H=1 on K!=0 and |K| there equals a constant
k=z sqrt(t)>0, because (2Phi(z)-1)/z is strictly decreasing. Put
r^2=|a|^2, p=t+r^2, e=exp(-z^2/2). The exact identities are

    t=(2Phi(z)-1)/(4phi(0)z),
    1=3t+r^2[1+(1-e)/(t z^2)],
    J=4phi(0)p sqrt(t).

In particular t<=1/3. The elementary inequality
(integral_0^z exp(-u^2/2)du)^2<=2(1-e) gives

    (1-e)/(t z^2)>=2t,
    p<=(1-2t+2t^2)/(1+2t)<=1/(1+4t).

Finally 1+4t>=4sqrt(t), yielding J<=phi(0). This proves (1).

## 6. Direct relevance to the current constructive response class

The canonical finite-resolvent anchor-plus-innovation frame has polynomial
inverse features in its finite Gaussian frame, despite its countably
supported canonical coefficients. It is therefore an analytic cyclic
core in the sense above. Appending W=UH(V) is strictly causal relative
to that core, even though its inverse is discontinuous. Each full-gradient
birth subsequently appends exactly one inverse feature measurable in
the previous frame. Randomized purification gates can likewise be
appended causally: choose an old even L2 function outside the finite
inverse-feature span, project off that span, normalize, and apply U.
This produces a new independent standard Gaussian with an old-measurable
inverse. An even function of it gives the independent uniform gate.
The same construction supplies any additional independent noise
coordinates needed for the conditional-ternary realization. Hence every
fixed finite stage of these constructive updates is covered by the
theorem.

The monotone feasible variational mechanism therefore cannot stop at a
finite-stage stationary point above phi(0). This is stronger than a
piecewise-polynomial nonstationarity observation: after the first update
the appended inverse generally contains softsigns and is not polynomial,
but it remains covered without any regularity assumption.

Important limits remain. The theorem supplies no lower bound on the
unrestricted gradient gap uniform over these finite stages, so it does
not force a particular limiting value. It does not justify running an
unbounded number of response stages at the same matrix dimension. It
does not rule out a high-value stationary frame with genuinely rough
cyclic core inverses, where the exact PSD-kernel recovery criterion is
still the applicable feasibility test.

## 7. Uniform escape from each fixed finite frame

There is a useful compactness strengthening, proposed by the director
and independently checked here. Fix ONE frame satisfying Section 1,
and fix c>phi(0). If the set of its feasible pairs with J>=c is nonempty,
then the full first-order gap has a strictly positive uniform lower bound

    gap(f,q)>=delta(frame,c)>0.                    (5)

The constant may depend on the actual inverse features, not merely their
number. No uniform bound over an expanding sequence of frames is claimed.

To prove this, put the bounded feasible set |f|<=q<=1, with the parity
constraints, in the product weak-* topology of L-infinity over the fixed
finite Gaussian space. It is weak-* closed and compact. Gaussian L1 is
separable, so this bounded compact topology is metrizable; sequential
arguments suffice.

The moments p=Eq and a=EfG are continuous, since 1 and each G_j belong
to L1. Thus K=a dot h varies strongly in L2 and t=p-|a|^2 is continuous.
The inequality

    |Psi(k,t)-Psi(k',t')|
      <=|k-k'|+E|N| |sqrt(t)-sqrt(t')|

shows that J is continuous, including at t=0. Hence the superlevel
set J>=c is compact. Moreover J<=sqrt(p), so p>=c^2 on this set. The
exact Gaussian moment body gives |a|<=m(p), where

    m(p)=2phi(Phi^(-1)(1-p/2)),
    t>=t_0:=min_(p in [c^2,1]) [p-m(p)^2]>0.       (6)

Strict positivity follows from m(p)^2<p for p>0 and compactness. This
uniform residual-variance margin makes all derivative expressions below
well behaved.

For a pair in this compact set let

    w=(1-q)s(K)-2BK,
    A_j=<w,h_j>,
    sigma^2=||w||_2^2-|A|^2>=0.

The scalar B is continuous: its bounded integrand varies strongly in
L1 as K,t vary, while q varies weak-*. The vector w varies weakly in
L2. Indeed for a fixed L2 test z, the factor z s(K) belongs to L1;
the changes in s(K) converge strongly in L2 by (6), and q remains
bounded. Thus A is continuous, while sigma^2 is lower semicontinuous
by weak lower semicontinuity of the L2 norm. All these quantities are
uniformly bounded on the compact superlevel set.

The full gradient has the conditional representation

    V=A dot G+sigma Z,

where Z is standard Gaussian independent of the fixed old frame. Its
pointwise full gap is therefore

    gap = E_G,Z[(|A dot G+sigma Z|+B-Psi(K,t))_+]
          -a dot A-Bp+E[q Psi(K,t)].              (7)

The second line is continuous in the weak-* topology. The first term
is continuous in A,B,sigma,K,t in their stated finite/strong-L1
topologies, by the Lipschitz positive-part and absolute-value bounds.
It is nondecreasing in sigma: adding independent centered Gaussian
noise increases the expectation of the convex function
(|v|+B-Psi)_+. Consequently its composition with the lower-semicontinuous
sigma is lower semicontinuous. For example, along a sequence realizing
the liminf, take a subsequence on which the bounded sigma_n converges;
its limit is at least sigma, and monotonicity gives the required bound.

Thus gap is lower semicontinuous and attains a minimum on the compact
superlevel set. A zero gap is exactly full first-order stationarity:
the pointwise optimum already ranges over every feasible value of
(f,q), and depends only on the old variables and this one independent
Gaussian. Sections 1--5 rule out zero at J>=c>phi(0). The attained
minimum is therefore positive, proving (5).

This also gives a uniform positive ONE-STEP energy gain from the fixed
frame into a possibly enlarged frame. Let the full linearized maximizer
be the other endpoint of a convex line search. Along theta<=1/2, its
residual variance is at least t_0/2. Using |a_new-a|<=2 and
|p_new-p|+2max(|a|,|a_new|)|a_new-a|<=5, the explicit Hessian bound in
the variational-stationarity note is uniformly at most

    M_0=4+26phi(0)/sqrt(t_0/2)
           +20phi(1)/(t_0/2)
           +25phi(0)/[2(t_0/2)^(3/2)].

The step theta=min(1/2,delta/M_0) improves the full response value by
at least theta delta/2. Once a new Gaussian inverse is appended, the
frame has changed; neither this delta nor a resulting iteration rate
is automatically uniform for later stages.
