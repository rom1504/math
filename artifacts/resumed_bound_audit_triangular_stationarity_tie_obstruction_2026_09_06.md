# Finite causal frames have no high-value stationary point; the zero-gradient tie exception is real

Date: 2026-09-06. Joint reconstruction with the director. This note
corrects the initially proposed pure-noise stationary ceiling: ties with
an identically zero gradient allow larger stationary values. The correct
rigorous obstruction proved here is 1/sqrt(2 pi), which is already below
the banked lower bound.

## 1. Statement and exact scope

Let G_1,...,G_d be independent standard Gaussian first-chaos coordinates
with orthonormal even inverse features h_j=U^{-1}G_j. Assume the strictly
causal property

    h_j is measurable in G_1,...,G_(j-1).             (1)

The features may be bounded measurable, discontinuous, or nonpolynomial;
no canonical-tree assumption is needed. Let odd f and even q satisfy
|f|<=q<=1 and be measurable in this finite frame. Use the conditional
ternary functional, notation, and realizability of
`resumed_bound_audit_variational_stationarity_dual_2026_09_06.md`:

    p=Eq, a_j=EfG_j, K=sum a_j h_j, t=p-|a|^2,
    J=E(1-q) Psi(K,t), Psi(k,t)=E|k+sqrt(t)N|.

Suppose the pair is stationary against ALL feasible first variations,
including the single-Gaussian coordinate extension of that note. Then

    J<=1/sqrt(2 pi)<0.4.                            (2)

Consequently every finite causal-frame response above this value has a
strictly improving feasible first-order direction, implemented by the
full-gradient threshold followed by a sufficiently small line-search
step. This is not an assertion of a uniform improvement size. It does
not apply to finite frames whose inverse features depend on their own
coordinates through an infinite-tree fixed point.

## 2. First variation and an innovation outside the frame

Write H=1-q, s(K)=2Phi(K/sqrt(t))-1, and

    B=E[H phi(K/sqrt(t))]/sqrt(t),
    w=H s(K)-2B K,   V=U w.

At positive value, p>0 and t>0. The latter follows from the exact
Gaussian moment-body bound |a|<=2phi(Phi^{-1}(1-p/2))<sqrt(p).
The first variation is E[V delta f+(B-Psi)delta q]. Its pointwise
maximum obeys

    f=sign(V)q,
    q=1 where |V|+B>Psi,
    q=0 where |V|+B<Psi,                            (3)

with arbitrary tie mixtures. If V has nonzero independent Gaussian
innovation outside the old frame, conditioning on the old frame makes
the optimal linearized value strictly greater than its value at the old
conditional mean. The reason is strict Jensen for
(|v|+B-Psi)_+ under a nondegenerate Gaussian: its two-sided unbounded
tails cross an interval on which its slopes differ. The current old-
measurable pair cannot attain that larger value. Thus stationarity
forces V=sum_j A_j G_j, where A_j=E[w h_j].

## 3. A nonzero old-frame gradient is impossible at positive value

If a=0, then K=0, s(K)=0, and V=0; moreover

    J=(1-p)sqrt(p)sqrt(2/pi)<=sqrt(8/(27 pi)).

Otherwise let k be the largest index with a_k!=0, and let ell be the
largest index with A_ell!=0. First assume V is not identically zero.
By causality K is measurable in G_<k.

If ell>=k, conditional on G_<k the tail part of V is a nondegenerate
Gaussian. The boundary ties in (3) have zero conditional probability.
Gaussian integration for the signed threshold gives, for every j>=k,

    a_j=A_j C,  C>0,                               (4)

where C is the same scalar for all these indices. Explicitly C is the
expectation of the sum of the two normal boundary densities divided by
the tail standard deviation, with threshold max(Psi-B,0). It is finite
and strictly positive. Since a_j=0 for j>k, (4) forces A_j=0 for j>k
and A_k!=0. Then w=sum_(j<=k) A_j h_j and K are both G_<k-measurable.
So H s(K) is G_<k-measurable as well.

But on any old-frame event where K!=0 and Psi(K)>B, the threshold H in
(3) is a genuinely nonconstant function of the nondegenerate G_k: its
conditional probability lies strictly between zero and one. This is
incompatible with measurability of H s(K). On Psi<=B, H=0 except a
zero-probability Gaussian boundary; on K=0, s(K)=0. Hence H s(K)=0
almost surely. Therefore w=-2B K, and

    E fV=-2B |a|^2<0,

contradicting E fV=E q|V|>=0. Here B>0 because J>0 implies EH>0 and
the normal density is everywhere positive.

It remains to justify the often-missed case ell<k. Now w,K,V are all
G_<k-measurable. On K!=0, the identity Hs(K)=w+2BK makes q earlier-
measurable. On K=0, (3) makes q earlier-measurable except on

    {|V|=Psi(0,t)-B}.

This is a fixed level set of a nondegenerate Gaussian and has probability
zero. Also V=0 has probability zero, so f=sign(V)q is earlier-measurable
everywhere relevant. That forces a_k=E fG_k=0, a contradiction.

Thus positive-value stationarity in a finite causal frame requires V=0.

## 4. Zero-gradient ties: exact reduction and the rigorous 0.4 bound

Assume a!=0 and w=0. Let Psi0=Psi(0,t)=2sqrt(t)phi(0).
If B<Psi0, (3) forces q=0 everywhere, impossible. If B>Psi0, then
q=1 and H=0 wherever K=0. On K!=0, w=0 implies

    H=2B sqrt(t) z/s(z), z=|K|/sqrt(t), s(z)=2Phi(z)-1.

Insert this into the defining equation for B and divide by B:

    1=E_(K!=0)[2z phi(z)/s(z)].                    (5)

Every integrand is strictly less than one because
s(z)=2 integral_0^z phi(u)du>2zphi(z), and the right side is strictly
less than one. This contradiction proves B=Psi0.

Consequently q=0 and H=1 wherever K!=0, because Psi(K)>Psi0 there.
The equation w=0 then forces |K| to have one constant nonzero magnitude
k on that event: s(z)/z is strictly decreasing for z>0. Write

    z=k/sqrt(t), e=exp(-z^2/2), r^2=|a|^2,
    alpha=P(K!=0)=r^2/(t z^2).

The exact equations are

    t=s(z)/(4 phi(0) z),
    1=3t+r^2[1+(1-e)/(t z^2)],
    J=2Bp=4phi(0)p sqrt(t), p=t+r^2.               (6)

In particular t<=1/3. Here is a simple exact bound, avoiding a numerical
maximization of the scalar family. Put I(z)=integral_0^z exp(-u^2/2)du.
The derivative of 2(1-e)-I(z)^2 equals

    2 exp(-z^2/2)(z-I(z))>=0.

Thus s(z)^2<=(4/pi)(1-e), and (6) gives

    (1-e)/(t z^2)>=2t,
    r^2<=(1-3t)/(1+2t),
    p<=(1-2t+2t^2)/(1+2t)<=1/(1+4t).

The last inequality follows by expanding the difference, whose numerator
is -2t^2(3-4t)<=0 on [0,1/3]. Finally 1+4t>=4sqrt(t), so

    J<=4phi(0) sqrt(t)/(1+4t)<=phi(0).

This proves (2), including the degenerate-gradient tie branch.

## 5. Why the stronger pure-noise stationary claim was false

The zero-gradient branch is not a merely formal concern. The following
finite causal two-coordinate construction realizes it above the pure-
noise value. Set z=2 and define EXACTLY

    t=(2Phi(2)-1)/(8phi(0)),
    r^2=(1-3t)/(1+(1-exp(-2))/(4t)),
    p=t+r^2, alpha=r^2/(4t).

Choose two disjoint even events A_plus,A_minus in G_1, each of probability
alpha/2. Let

    h_1=1,
    h_2=alpha^(-1/2)(1_A_plus-1_A_minus),
    G_2=U h_2.

Then h_2 is mean zero, even, norm one; G_2 is standard Gaussian independent
of G_1, and h_2 is G_1-measurable. Choose 0<=l<u<=infinity so that

    P(l<|G_2|<u)=p/(1-alpha),
    E[|G_2|1{l<|G_2|<u}]=r/(1-alpha).

Such an annulus exists by continuity: at fixed probability its first
absolute moment ranges from the centered interval's moment to the tail
interval's moment. The exact rational replay
`computations/resumed_bound_audit_tie_parameter_box.py` proves that the
requested moment lies between 0.256 and 0.259 and the requested
probability lies between 0.377 and 0.379. It encloses t by the alternating
integrated exponential sums through orders 16 and 17, and exp(-2) by the
alternating exponential sums through orders 18 and 19; all subsequent
interval arithmetic is rational, including squared bounds for r.

For clarity the annulus-existence comparison itself needs no numerical
normal quantiles. The centered interval of the required probability has
endpoint below 1/2, because

    P(|N|<1/2) >=(23/24)phi(0)>(23/24)*.398>.379.

Its absolute first moment is therefore at most .5*.379<.256. The tail
interval of that probability has endpoint greater than 3/4, because

    P(|N|>3/4)>=1-(3/2)phi(0)>1-(3/2)*.4=.4>.379.

Its absolute first moment is at least .75*.377>.259. The classical
rational bounds on pi imply .398<phi(0)<.4. Continuous movement of a
fixed-probability annulus between these two endpoints supplies exactly
the desired first moment.

Set F=sign(G_2) on this annulus intersected with the complement of
A_plus union A_minus, and set F=0 elsewhere. Then H=1-|F|, a=(0,r),
K=r h_2, and all equations (6) hold. On K!=0, F=0 and Psi>B; on K=0,
Psi=B and the entire support of F is a zero-gradient tie. Thus every
feasible first derivative is nonpositive even though

    J=0.3130615374... >sqrt(8/(27pi))=0.3071059106... .

The same rational replay proves J>.31, while the pure-noise cap is below
.308, so the strict separation does not rely on the displayed decimal
evaluations. They are not needed for the rigorous high-value obstruction.
This stationary construction is
not asserted to be a local maximum. The full theorem says that no such
stationarity, degenerate or otherwise, can persist above phi(0) in a
finite causal frame.
