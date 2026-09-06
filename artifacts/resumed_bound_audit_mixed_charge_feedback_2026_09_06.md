# Independent mixed-charge feedback audit and quantitative scope

Date: 2026-09-06. Checkpoint: 06:32 UTC.

Outcome: PASS for the finite-depth Haar paired-query recursion,
nondegeneracy, physical hard-sign passage, and retained-energy theorem
in `resumed_response_mixed_charge_involution_feedback_2026_09_06.md`.
The two-step distance certificate was also independently replayed.
Sections 4--5 below add an explicit bounded-depth lower bound and a
constructive obstruction to a depth-uniform gain.

## 1. Conditional Haar law and the population recursion

After exposing queries v_j and replies Bv_j, the span K of all these
vectors is B-invariant and B restricted to K is known. Inductively,
conditional on the history and independent seed data, the remaining
restriction is Haar on the involution orbit with its remaining positive
and negative multiplicities. An adaptive next query is already a fixed
vector under this conditioning; adaptivity imposes no further matrix
constraint.

Remove the known projection of the next query, and write its residual
as L e, where ||e||=1 and e lies in K-perp. If r and s eigendirections
remain and d=r+s, the exact longitudinal coefficient is

    z=<e,Be> ~ 2 Beta(r/2,s/2)-1,
    E z=(r-s)/d,
    Var z=2(1-((r-s)/d)^2)/(d+2).

Conditional on z, Be=z e+sqrt(1-z^2) xi, where xi is uniform on the
sphere in K-perp intersect e-perp. On the newly exposed plane, the exact
matrix in basis (e,xi) is

    [[z,sqrt(1-z^2)],[sqrt(1-z^2),-z]].

This removes one positive and one negative eigendirection. Hence an
exactly balanced starting complement stays exactly balanced; a fixed
number of removals also preserves asymptotic balance for the broader
asymptotically balanced version. The posterior on the new complement
is again its Haar orbit by the stabilizer fixing this plane. It is not
necessary, and would be incorrect, to assert e and Be are exactly
orthogonal at finite n.

Couple xi by projecting a fresh standard Gaussian vector onto this
fixed-codimension complement and normalizing. The removed Gaussian
squared norm has conditional expectation O(depth). At fixed depth,
the normalized-L2 projection and normalization errors vanish. Also
z->0. Thus the plane block becomes a swap in the limit. Induction gives
the displayed orthonormal pairs (q_j,p_j), new independent Gaussian
response, deterministic empirical moments, and empirical W2 law.
Finite-rank projection costs do not require coordinate incoherence of
the exposed vectors. Positive limiting residual variance is needed when
normalizing each new residual.

For a Boolean query u, projection coefficients a_j=<u,q_j>,
b_j=<u,p_j> give

    y=k+sigma N,  k=sum_j(a_j p_j+b_j q_j),
    E k^2=1-sigma^2,  E uy/2=sum_j a_j b_j.

The factor in the half-energy is correct. These are analysis
coefficients, not additional population-oracle operations in the actual
coordinate algorithm. The floating Sobol code follows this recursion;
its finite-sample orthogonalization does not make it an exact integral
certificate.

## 2. Hard signs and retained energy

For u_next=sign(k+sigma N+alpha u), condition on all row history before
the newest N. Every vector in the newly exposed finite span is affine
in N. The next spin is a nonconstant threshold of N whenever sigma>0
and alpha is finite. Such a threshold cannot equal an affine function
almost surely. Thus its residual variance is positive, starting from
sigma_0=1. This is the needed induction, not just a generic assertion
that nonlinear updates have positive innovation.

The score has conditional density bounded by 1/(sigma sqrt(2pi)).
Consequently its probability in a width-2h threshold window is at most
2h/(sigma sqrt(2pi)). This controls normalized-L2 hard/soft-sign error
and makes finite-n exact-zero conventions irrelevant. At each fixed
depth, choose approximation errors before taking n to infinity, and
remove them afterward. Bounded operator norm propagates field errors;
the no-atom bound controls the following discontinuous update. The
audited WZF module transfers these fixed equivariant approximations to
exact incoherent involutions and the explicitly operator-close hollow
class. No arbitrary-signing inference follows.

For v=sign(Bu), g_n=n^-1 sum(|(Bu)_i|-u_i(Bu)_i), and
m_eta=(1-eta)u+eta v, the exact quadratic expansion and ||B||op<=1 give

    [Q_B(m_eta)-Q_B(u)]/n >= eta g_n-2 eta^2.

The limiting g is strictly positive when sigma>0. Writing
Gamma_sigma(z)=E|z+sigma N|-z for z>=0, the function
x -> Gamma_sigma(sqrt(x)) is convex, since its second derivative has
numerator z Gamma_sigma''(z)-Gamma_sigma'(z)>0. Therefore

    g >= Gamma_sigma(sqrt(1-sigma^2)) > 0.

Also g<=2. The fixed choice eta=g/4 gives gain at least g^2/8.
Independent coordinate rounding preserves normalized energy for hollow
B, or for this incoherent involution class where the maximum diagonal
vanishes. For a general nonhollow unit-op matrix there would be a
diagonal mean correction, so that broader rounding assertion is not
used. Conditional off-diagonal rounding variance is O(n).

## 3. Exact two-step departure from the whole old frame

The source's decomposition at alpha_0=alpha_1=1/2 was reconstructed
including the coupling at B and -B. Let H be the central mask,
F=sign(Z)(1-H), mu=EH, b=E ZF, and tau^2=EF^2-b^2. Then

    Bu_1=W+bS+tau N,
    W=U H,  Var W=mu,  Cov(W,Z)=mu,

where W is in the old Gaussian frame and N is independent of that whole
frame and S. Independence from the countable old sigma-field follows
from independence with every finite subfamily; W itself is its L2
first-chaos limit. Matrix-sign charge keeps W odd and BF=bS+tau N even
in the actual common-seed coupling.

With A=W+alpha_1 F and c=b+alpha_1 H, the odd half of the second update
is exactly

    D_2=sign(A) 1{|A|>|c+tau R|}, R=SN.

R remains independent of the old frame and S. Thus conditional variance
is p(1-p), where p=P(|c+tau R|<|A| | old frame). It is positive almost
surely because Var(W|Z)=mu-mu^2>0 and tau>0.

The exact rectangle calculation is correct: on |Z|<=1/4 and
1/2<=W<=3/5, its mass is >3/640; the gate probability is >1/240 and
<1/2. Hence the squared distance from the complete old sigma-field is
strictly greater than 1/102400. Independent replay of
`computations/resumed_response_mixed_charge_distance_certificate_2026_09_06.py`
is byte-identical to the canonical result, saved at
`tmp/resumed_bound_audit_mixed_charge_distance_replay_2026_09_06.json`.
This is a structural distance certificate, not an energy bound.

## 4. Explicit uniform innovation for bounded depth and inertia

There is a quantitative strengthening of the nondegeneracy argument.
Put s(z)=2 Phi(z)-1 and

    d(z)=1-s(z)^2-4 phi(z)^2.

This is the squared L2 error of the best affine approximation to
sign(N+z) using 1 and N. It is positive for every finite z. It is even
and strictly decreasing for z>0, since

    d'(z)=4 phi(z) [2z phi(z)-s(z)] < 0.

At the next update, enlarge the already exposed linear space to all
functions a(history)+b(history)N of the newest normal. Enlargement can
only reduce projection error, so

    sigma_(t+1)^2 >= E d((k_t+alpha_t u_t)/sigma_t).            (10)

Assume |alpha_t|<=A. Since ||k_t||_2<=1 and |u_t|=1,
E(k_t+alpha_t u_t)^2<=(1+A)^2. Chebyshev gives probability at least
one half to |k_t+alpha_tu_t|<=sqrt(2)(1+A). Therefore the explicit
deterministic recursion

    s_0=1,
    s_(t+1)=sqrt(d(sqrt(2)(1+A)/s_t)/2)                       (11)

satisfies 0<s_t<=sigma_t for every schedule in [-A,A]^t. In particular
the retained gain after step t is uniformly at least

    Gamma_(s_t)(sqrt(1-s_t^2))^2/8.                          (12)

Here Gamma increases with its variance parameter and decreases with its
nonnegative argument. The bound is extremely small after even modest
depth, but is explicit and positive. This is uniform for the stated
bounded-depth, bounded-inertia class; it is not uniform over all finite
depths or all arbitrary coordinate algorithms.

## 5. Why no positive depth-uniform stability gap is possible

The general finite paired-query class has a constructive obstruction,
not merely an appeal to a performance ceiling. Fix 0<epsilon<1 and
eta=epsilon/4. Starting from the initial Boolean seed, repeatedly query
y=Bu, put v=sign(y), and independently at each coordinate choose the
next spin to be v with probability eta and u otherwise. The conditional
mean is m_eta above. Implement the choices by thresholding an unused
even Gaussian gate, so the rule remains signed-permutation equivariant.

Each fixed step is legitimate. If the preceding query has positive
innovation, then g>0 and disagreement has positive probability. The
fresh independent selector gives conditional variance
4 eta(1-eta) on that disagreement event. Since the old exposed span
does not use this selector, the next query again has positive residual
variance. All threshold atoms vanish; finite Lipschitz approximation and
the same deterministic-matrix transfer apply.

Let e_t be its deterministic limiting half-energy. Conditional rounding
and the retained inequality imply

    e_(t+1)-e_t >= eta g_t-2 eta^2.

If g_t>=epsilon this increment is at least epsilon^2/8. Since e_0=0
and e_t<=1/2 from the elementary unit-operator-norm bound, among the
first floor(4/epsilon^2)+1 queried states at least one has g_t<epsilon.
The index is chosen from the deterministic population sequence and is
independent of n. Thus there are genuine fixed finite-depth Boolean
feedback states with arbitrarily small positive gap. Their innovations
must tend to zero by the positive bound in Section 2.

Combining this construction afterward with the separate fixed-GFOM
Haar ceiling e_t<=sqrt(15)/8 shows that arbitrarily small stability gaps
can coexist with a fixed positive deficit from 1/2. Hence no positive
uniform gain depending only on that deficit can hold for this general
finite feedback class. This does not prove that the specific synchronous
inertial schedule has gaps tending to zero, and it says nothing about
the optimum over all signings or all Boolean vectors.

## 6. Audit of strict extension and finite-algorithm nonattainment

The subsequent theorem in
`resumed_response_paired_history_actual_energy_extension_2026_09_06.md`
also passes independent reconstruction, including its Section 5.
For a fixed paired history H, the moment body
{E[Hf]: |f|<=1, f odd} is compact by weak-* compactness and H in L1.
The actual energy is c^T Jc/2. At a maximizer, the directional derivative
forces f=sign(H^T Jc) away from ties.

If ||c||=1, Bessel equality gives f=H^T c Boolean almost surely.
Descending through the newest independent Gaussian p_j=N_j, boundedness
first forces its coefficient b_j=0. Then f is measurable before N_j,
whereas a nonzero coefficient a_j in H^T Jc would make its conditional
sign take both values. Thus a_j=0 also. Repetition eliminates all
coefficients, a contradiction. Ties do not invalidate this argument:
when a_j!=0 the conditional Gaussian has no zero atom. Compactness of
the argmax set supplies the claimed positive residual floor, and also
the uniform floor for sufficiently near-optimal moment vectors.

The fresh even-gate/odd-seed Boolean rounding preserves c exactly.
One genuine additional query has variance 1-||c||^2 bounded below, so
the retained-energy inequality gives a fixed positive improvement.
Arbitrary bounded measurable moment-body maximizers should be understood
in the ordered-L2 variational closure. Choose fixed Lipschitz and gate
approximants with error below a fraction of the strict gain to obtain a
genuine finite algorithm exceeding the old optimum. This distinction
does not create an attainment or supremum gap.

For a fixed globally Lipschitz GFOM, a zero-residual query is an L2
linear combination of previously exposed query/response vectors.
Bounded operator norm identifies its response with the swapped
combination. Replacing that query and propagating through the remaining
finite Lipschitz computation preserves normalized-L2 limits. Repeating
leaves one fixed nondegenerate limiting frame, not a sequence of frames
with potentially vanishing strict-extension constants. Therefore a
finite GFOM attaining the global fixed-rule supremum would yield one
finite frame attaining it, and its strict extension is a contradiction.

This verifies nonattainment of the fixed-rule supremum, not its value,
not a uniform gain over histories, and not convergence of the original
signing optimum.
