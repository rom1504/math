# An explicitly evaluated mixed-charge feedback mechanism

Date: 2026-09-06. Status: exact algebra and finite integer replay; the
adaptive Haar-conditioning proof, nondegeneracy, and retained-energy
increment were independently reconstructed by the convergence agent,
the bound-audit agent, and the director. See
`resumed_convergence_haar_adaptive_query_kernel_2026_09_06.md`.
Floating state-evolution integrals are explicitly
not numerical certificates. Matrix scope: balanced Haar involutions,
then the fixed-rule universality class in
`resumed_convergence_involution_local_ceiling_2026_09_06.md`.

## 1. The rule, and why this is outside the terminal architecture

Start with iid Rademacher spins u_0=S and iterate a FIXED number of times

    y_t=B u_t,
    u_(t+1)=sign(y_t+alpha_t u_t).                 (1)

At an exact zero score keep u_t. The inertia schedule used below is

    alpha_t=1/2 for 0<=t<4,
            3/10 for 4<=t<8,
            1/5 for 8<=t<16,
            1/10 for 16<=t<40.                  (2)

This is a finite signed-permutation-equivariant feedback rule. Under
B->-B, u_0 is unchanged whereas y_0 changes sign; already u_1 mixes both
matrix-sign charges. Later queries use this whole mixed-charge output.
Its charge-odd part need not be an unmarked function of the old Gaussian
creation frame; the next subsection proves this explicitly at two
updates. Thus the terminal functional's <.45 ceiling does not apply.
Section 6 of
`resumed_response_flat_involution_terminal_endpoint_2026_09_06.md`
gives the precise charge normal form and its limit of applicability.

The rule (1) need NOT increase the energy at every step for every
matrix. Neither synchronous sign updates nor small inertia have such a
general guarantee. Section 4 gives a separate rigorous damped update.

### 1.1 An explicit failure of old-frame measurability at two updates

Take alpha_0>0 and write Z=BS. Set

    H=1{|Z|<=alpha_0}, F=sign(Z)(1-H),
    u_1=F+S H,
    mu=E H in (0,1), b=E[ZF]=2 phi(alpha_0),
    tau^2=E F^2-b^2>0.

For the original normalized Hadamard coupling, the proved full local
channel and the actual creation relation give

    B u_1=W+b S+tau N, W=U H,                   (2a)

in the joint limiting local law. Here W belongs to the complete old
canonical Gaussian frame, Var W=mu, Cov(W,Z)=mu, and N is an independent
standard Gaussian, also independent of S and every old canonical
coordinate. Independence from the whole countable sigma-field follows
from independence with each finite subfamily. The two inputs on the
right of (2a) must not be merged into a single Gaussian if matrix-sign
charges are to be tracked: W is odd-charge and bS+tau N is even-charge.

The analogous Haar statement additionally uses the elementary unmarked
Haar response module in
`resumed_convergence_haar_response_traffic_2026_09_06.md`, Section 4;
it does not follow from only the old-tree recursion or the finite paired
query theorem. The original Hadamard architecture-separation claim does
not require that additional dependency.

Couple the computations at B and -B using the same seeds. The first
update is respectively F+S H and -F+S H, so their next threshold scores
are

    A+T and -A+T,
    A=W+alpha_1 F,
    T=S(b+alpha_1 H)+tau N.

Thus the charge-odd half D_2 of u_2 satisfies

    D_2=sign(A) 1{|A|>|b+alpha_1 H+tau R|},
    R=S N.                                    (2b)

The Gaussian R is independent of S and of the complete old frame.
Conditional on Z, W has variance mu-mu^2>0; therefore A is nonzero
almost surely. Given the old frame, the gate in (2b) has probability
strictly between zero and one. Consequently

    Var(D_2 | complete old canonical frame)=p_*(1-p_*)>0 almost surely,
    p_*=Phi((|A|-c)/tau)-Phi((-|A|-c)/tau),
    c=b+alpha_1 H.

In particular the exact squared L2 distance from that whole old sigma-
field is E[p_*(1-p_*)], not merely a qualitative nonmeasurability claim.

At alpha_0=alpha_1=1/2 there is the elementary exact bound

    E[p_*(1-p_*)] > 1/102400.                   (2c)

Here are all ingredients. Exact outward Gaussian intervals give
mu in (.38,.39), b in (.70,.71), tau in (.32,.37), and
sqrt(mu-mu^2) in (.48,.49). On the rectangle
|Z|<=1/4, 1/2<=W<=3/5, one has H=1, A=W and c=b+1/2 in (1.20,1.21).
The conditional W standardized coordinate has absolute value below
1.46, so its density is >2 phi(1.46)>1/4. Also phi(1/4)>3/8.
The rectangle therefore has probability >3/640.
The gate interval contains [-5/2,-9/4], giving
p_*>phi(5/2)/4>1/240. Its upper endpoint is negative, so p_*<1/2
and p_*(1-p_*)>1/480. Their product is 1/102400.
Every displayed transcendental enclosure and rectangle comparison is
replayed with exact Fractions in
`computations/resumed_response_mixed_charge_distance_certificate_2026_09_06.py`;
the canonical output has the same basename under `computations/results/`.

This proves that D_2 is not an unmarked old-frame function in the
limiting actual coupling. One could manufacture an unrelated even gate
with the same conditional marginal distribution from unused old
coordinates, but that would not preserve its correlation with the
already exposed field BF. It is precisely that correlation which the
feedback uses. Equality of one-root marginals does not put the actual
coupled response back into the terminal architecture.

## 2. Paired-query state evolution on the Haar comparator

Expose adaptive matrix queries one at a time. The span of every exposed
query and its response is B-invariant. Conditional on all exposed data,
the restriction of B to the orthogonal complement is still a Haar
involution, with one positive and one negative eigendimension removed
by each nondegenerate query/response pair.

For clarity the limiting recursion is expressed in orthonormal pairs
(q_j,p_j) on a scalar probability space, with B swapping q_j and p_j.
The initial query is q_0=S and its response is p_0=N_0, an independent
standard Gaussian. Given a subsequent scalar query u, set

    a_j=E[u q_j], b_j=E[u p_j],
    r=u-sum_j(a_j q_j+b_j p_j), sigma^2=E r^2,
    k=sum_j(a_j p_j+b_j q_j).                    (3)

If sigma>0, append

    q_new=r/sigma, p_new=N_new,
    y=k+sigma N_new,                            (4)

where N_new is an independent standard Gaussian. Then

    E u y/2=sum_j a_j b_j.                      (5)

These formulas are population identities, not empirical fitted
coefficients in the actual algorithm (1).

Here is the finite-dimensional conditioning justification. If a
residual query v lies in the unexplored complement, revealing Bv fixes
the invariant plane span(v,Bv). With probability one v is not an
eigenvector, so this plane contains one positive and one negative
eigenvector. The remaining conditional distribution is invariant under
all orthogonal transformations of the new complement; hence it is the
Haar orbit of the remaining involution spectrum. The longitudinal
overlap <v,Bv>/||v||^2 has the balanced beta law up to a fixed-rank
imbalance and tends to zero. Conditional on that overlap, the transverse
direction of Bv is uniform on the appropriate sphere.

For any fixed number of queries, replace that uniform sphere by a fresh
Gaussian vector projected off the finite exposed span and normalized.
The finite-rank projection and normalization have vanishing normalized
L2 effect. This yields (3)--(4), by induction, for fixed Lipschitz
updates and their empirical second moments. It also explains why one
must not assert exact orthogonality of v and Bv at finite n.

The hard rule (1) has positive innovation variance at every FIXED step.
Indeed, condition on the history before the immediately preceding
fresh Gaussian N. The next spin is a nonconstant sign of an affine
function with nonzero coefficient of N. Every vector in the previously
exposed linear span is, conditionally, an affine function of N. A
nonconstant sign cannot equal an affine Gaussian function almost surely.
Thus its residual variance in (3) is strictly positive. Induction starts
at sigma_0=1. The same positive variance gives a conditional continuous
density for each threshold score, so hard signs are normalized-L2 limits
of globally Lipschitz equivariant soft rules. The existing universality
module then transfers each fixed finite rule to exact flat involutions
and their stated operator-close hollow perturbations. No depth-growing
statement or quantitative finite-n rate is asserted.

## 3. Explicit evaluations and their different statuses

The floating scalar recursion is implemented in

    computations/resumed_response_mixed_charge_involution_feedback_2026_09_06.py

using scrambled Sobol Gaussian integration. Finite-sample Gram correction
is used to reduce spurious accumulated overlaps; it disappears in the
population recursion. Consequently these computations are diagnostics,
not exact integration bounds.

For schedule (2), 2^20 points and scramble seed 831 give energy

    t=20: .4462669993,
    t=30: .4506474710,
    t=40: .4516227082.

An independent 2^18-point run with seed 717 gives .45171239 at t=40.
The canonical 2^20-point output is
`computations/results/resumed_response_mixed_charge_involution_feedback_schedule_2026_09_06.json`.
These support, but do NOT prove, an asymptotic crossing of .45.

There is also an EXACT INTEGER finite-matrix replay:

    computations/resumed_response_mixed_charge_hadamard_replay_2026_09_06.py

For the symmetric Sylvester matrix H at n=2^18, B=H/512, seed 617,
and schedule (2), the final Boolean vector has

    Q_B(u_40)/n=30324543/67108864
                =.45187090337276458740234375.    (6)

All Walsh-Hadamard products, inertia comparisons, and the final energy
are exact integers; there is no floating comparison in this replay.
The exact-zero convention preserves the current Boolean spin. Sylvester
H has trace zero, so deleting its diagonal leaves every Boolean energy
unchanged. Thus (6) also describes an actual hollow sign matrix in the
original class. It is an example at a fixed size, not a universal lower
bound and not an asymptotic certificate.

Independent finite examples from the same source are

    n=2^20, seed103: 121260049/268435456=.4517288841307163...,
    n=2^20, seed717: 242447589/536870912=.4515938255935907....

The main replay JSON is
`computations/results/resumed_response_mixed_charge_hadamard_replay_2026_09_06.json`.

## 4. A genuine energy-increment invariant after any nondegenerate query

There is an exact dynamic inequality which does not rely on treating a
new mixed-charge output as an old terminal rule. For symmetric ||B||op<=1,
a Boolean u, y=Bu, v=sign(y) (an arbitrary Boolean choice at zero), put

    g_n=(1/n) sum_i (|y_i|-u_i y_i),
    m_eta=(1-eta)u+eta v, 0<=eta<=1.

Then m_eta is feasible and the quadratic identity gives

    [Q_B(m_eta)-Q_B(u)]/n
       >=eta g_n-(eta^2/2n)||v-u||^2
       >=eta g_n-2eta^2.                        (7)

For the nondegenerate population query (3)--(4),

    g=E[Psi(k,sigma^2)-u k]>0.                   (8)

The strict positivity follows pointwise from
Psi(k,sigma^2)>|k|>=u k. Moreover E k^2=1-sigma^2, and the function

    x -> Gamma_sigma(sqrt x),
    Gamma_sigma(z)=Psi(z,sigma^2)-z for z>=0,

is convex: its second derivative has numerator
z Gamma_sigma''(z)-Gamma_sigma'(z)>0. Jensen therefore gives

    g>=Gamma_sigma(sqrt(1-sigma^2))>0.           (9)

Because g<=2, choosing the FIXED eta=g/4 is feasible and (7) gives
an asymptotic energy gain at least g^2/8. Independent Boolean rounding
of m_eta preserves its normalized limiting energy (conditional error
variance O(n)), so this is an actual further signing improvement on
the flat-involution class. If desired, softsign first gives an arbitrarily
close gain without relying on a convention at zero.

This is a strict improvement after EVERY fixed finite iterate (1), since
Section 2 proved sigma_t>0. It does not produce a uniform increment over
depth: sigma_t, g_t, and hence the bound may tend to zero. Nor does it
prove that indefinite feedback reaches 1/2, or even crosses .45; the
fixed-GFOM Haar upper ceiling .4841229... remains in force. The concrete
distinction from terminal mask optimization is that (7) improves the
actual retained endpoint energy and queries its entire mixed-charge
state. It does not reset the proof to a single final J(F,H).

## 5. A quantitative, but depth-dependent, innovation floor

The director supplied the following useful strengthening of Section 2.
For N standard Gaussian, the squared distance of sign(N+h) from
span{1,N} is exactly

    v(h)=4 Phi(h) Phi(-h)-4 phi(h)^2>0.

It is even and nonincreasing in |h|, because for h>=0

    v'(h)=-4 phi(h)[2 Phi(h)-1-2h phi(h)]<=0.

At the next update, every old exposed vector is conditionally affine
in the immediately preceding fresh N. Enlarging that span to ALL
conditionally affine functions can only decrease the residual distance.
Therefore

    sigma_next^2 >= E v((k+alpha u)/sigma)
       >= (1/2) v(sqrt(2 E(k+alpha u)^2)/sigma)>0. (10)

The last inequality is Markov plus monotonicity; when the second moment
is zero it remains valid directly. The moment in (10) is

    E(k+alpha u)^2=1-sigma^2+alpha^2+4 alpha e,
    e=E[u Bu]/2,

and can alternatively be bounded by
(sqrt(1-sigma^2)+|alpha|)^2. This makes nondegeneracy reproducible at any
fixed depth, but repeated use can be extremely weak. It gives no
depth-uniform positive innovation and no asymptotic energy guarantee
beyond those explicitly stated above.
