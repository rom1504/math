# Restricted nonlinear-channel theorem and the center update

Date: 2026-09-06. Status: proved, with independent director and response-track
reconstructions in `resumed_nonlinear_channel_director_audit_2026_09_06.md`
and `resumed_response_center_theorem_independent_audit_2026_09_06.md`.
This is NEW follow-up mathematics, not a dependency of the separately
reconstructed 0.4306581794055286 bound. No extra decimal is asserted.

The detailed deterministic forest lemma, first-chaos contraction split,
and input-transfer bounds are in
`resumed_bound_audit_nonlinear_forest_flattening_2026_09_06.md`, Sections
1--8. This note packages their precise consequence and the remaining
bounded-response/sign limits, including zero local variances.

## 1. Precise fixed-polynomial statement

Let B=A/sqrt(n-1), where A is an actual symmetric hollow signing, and
assume ||B||op<=L with FIXED L. Write Q=B^2. Let X_i be a fixed finite
ancestor-closed family of normalized injective marked odd-tree fields,
including G_i=(BS)_i. Its one-root Gaussian limit is denoted by X, with
independent coordinates G_T and edge coordinate G_0.

Let F(X) be a fixed odd polynomial. Define

    a_T=E[F(X)G_T],
    K_F(X)=sum_T a_T h_T(X).

This is an even polynomial of finitely many old coordinates. Let R be
a finite set of odd positive integers at least three, and put

    f_r=E[F(X)h_r(G_0)],
    s^2=sum_{r in R} f_r^2 > 0,
    h=sum_{r in R} (f_r/s) h_r,
    Z=B h(BS),
    v_i=(B [sum_{r in R}(f_r^2/s^2) Q^(circ r)] B)_ii.

Here h_r=He_r/sqrt(r!). The signs of f_r are retained in h, and only
their squares enter v_i. One has 0<=v_i<=L^2. Let d_i be any bounded
deterministic weights, uniformly in n. For every fixed even polynomial
M and fixed one-variable polynomial psi of the specified parity,

    (U)  n^-1 sum_i d_i E[M(X_i) psi_odd(Z_i) (BF(X))_i]
         = s E M(X) n^-1 sum_i d_i
                 E[(sqrt(v_i)N) psi_odd(sqrt(v_i)N)] + o(1),

    (S)  n^-1 sum_i d_i E[S_i M(X_i) psi_even(Z_i) (BF(X))_i]
         = E[M(X)K_F(X)] n^-1 sum_i d_i
                 E[psi_even(sqrt(v_i)N)] + o(1).

The errors are uniform over bounded deterministic weights and the
bounded-operator signing class. No convergence of the v_i is required.
N is standard Gaussian independent of X. These are averaged test
identities, NOT an unqualified statement of per-root conditional means.

The version with an arbitrary finite odd h (not aligned to F) has
covariance coefficient c_i=sum_r h_r f_r (B Q^(circ r) B)_ii in place
of s v_i in the Gaussian integration-by-parts form of (U). This version
is useful when approximating a bounded F while keeping h fixed.

## 2. Decomposition: separate the nongeneric first-chaos part

Expand the odd polynomial F in the independent local Gaussian Hermite
basis, and write

    F=sum_T a_T G_T+R_F.

Every monomial in R_F has odd LOCAL Hermite degree k>=3. This local
degree is not its original input-spin chaos degree. Confusing these two
degrees would invalidate the following argument.

Local Wick/forest replacement represents each such monomial by the
symmetrized tensor product of k old branch kernels, up to rootwise
o(1) L2 error. Whole equal-branch contractions are the Hermite
subtractions, distinct whole-branch contractions vanish, and proper
contractions vanish by old fixed-root flattening bounds. The Rademacher
version follows from the same restricted leading-pairing argument.
At fixed L these errors transport through B in averaged L2.

The forest flattening lemma proves that the transported main kernel

    sum_j B_ij tensor_{a=1}^k K_Ta,j

has bounded Hilbert norm and all proper flattenings O(n^-1/2),
uniformly in i. The proof classifies a slot cut by the number of
branches split by it; it does not assume the branches are independent
after a common root has been summed. Therefore each transported
monomial has vanishing fourth cumulant. Keep monomials of different
input-chaos degrees as separate vector components, use the multivariate
Gaussian-chaos theorem, and then add them.

Consequently B R_F is jointly Gaussian in the averaged asymptotic
sense with the fixed local old family and the finitely many unmarked
star channels. Its covariance with old variables can be retained as a
possible odd linear drift. Nothing below requires it to be zero.
The star channels are independent of the local old Gaussian family.
Output-own-input removal costs o(1), so the root spin is also
independent of this nonlinear Gaussian vector.

## 3. Identifying the nonlinear component's tested projection

The existing complete deterministic-weighted one-Z theorem gives, for
every r in R,

    n^-1 sum_i d_i E[(B R_F)_i (B h_r(G))_i]
      = f_r n^-1 sum_i d_i (B Q^(circ r) B)_ii + o(1).

The first-chaos local part has no nonlinear edge Hermite coefficient,
so removing it does not alter f_r. Uniformity over all bounded d_i
implies an averaged absolute covariance error: apply the statement to
the deterministic signs of those errors. Therefore the covariance of
B R_F with the chosen Z equals s v_i up to an averaged o(1) error.

For a joint Gaussian vector (X,Z,Y) with Cov(X,Z)=0, Var(Z)=v and
Cov(Y,Z)=c, Gaussian integration by parts gives

    E[Y M(X) psi(Z)]
      = sum_T Cov(Y,G_T) E[G_T M(X)] E psi(Z)
          + c E M(X) E psi'(Z).

The first term is zero when psi is odd (and also when M is even).
This identity remains valid at v=0 and does not divide by v. Gaussian
Stein gives v E psi'(sqrt(v)N)=E[sqrt(v)N psi(sqrt(v)N)]. Thus the
nonlinear part gives (U). Its independent root spin makes its
contribution to (S) vanish.

Compactness of bounded covariance matrices and subsequence contradiction
justify this Gaussian calculation without assuming the variances or
old drift covariances converge. All polynomial moments needed for the
passage are uniformly bounded for the forest main kernels.

## 4. First-chaos part: the two different mechanisms

For one old tree T, the one-level identity gives

    (B X_T)_i=sum_j Q_ij S_j h_T(X_j)+o(1)

in averaged L2. In general this is NOT Gaussian. Coherent entries of
Q can retain finite Rademacher-spin atoms, and no invariance-in-law for
this vector is used.

For (S), remove the own S_i from Z_i and use the exact Boolean identity

    E[S_i S_j f_i h_T(X_j)]
       = E[(D_j f_i)(D_i h_T(X_j))],  j!=i,
    f_i=M(X_i)psi(Z_i).

Both factors are own-free in the needed coordinates. Their derivatives
have Lp size O(n^-1/2), so each expectation is O(1/n). Since
sum_j |Q_ij|=O(sqrt(n)), the off-diagonal contribution vanishes.
The j=i term gives E[M h_T] E psi(sqrt(v_i)N) by the local old/Z CLT.
Sum with coefficients a_T to obtain (S).

For (U), work first with Gaussian input, retaining B X_T as a
potentially non-Gaussian chaos. Every mixed contraction of its kernel
with an unmarked-star kernel vanishes in averaged squared Hilbert norm.
Proper-star contractions follow from the star's small self-contractions.
For a full r-slot star contraction into B X_T, write the latter as
Q[N h_T(X)] and split derivatives according to whether they hit N_a.

The root-not-hit part is sum_a Q_ia N_a U_a,i, where
U_a,i=sum_l B_il D_bl^r h_T(X_a). Uniform joint moment independence of
old fields at ANY root a and the star channel at root i, followed by
the positive covariance-of-squares contraction formula, gives
max_{a,i}||U_a,i||_2=o(1). The U_a,i exclude N_a. The degree-bounded
Gaussian creation inequality

    ||sum_a q_a N_a P_a||_2^2
       <= (d+1) sum_a q_a^2 ||P_a||_2^2

then proves this class vanishes. To see that inequality, expand P_a
in the normalized multivariate Hermite basis. Multiplying by N_a is
pure creation because P_a excludes N_a. Each output multiindex can
receive at most its total degree <=d+1 contributions. Weighted
Cauchy--Schwarz on these contributions, followed by summation of
Hermite coefficients, gives the bound. For a polynomial with several
degrees, apply the same bound in each orthogonal chaos and add.

The root-hit part is, up to a fixed factor,

    diag[Q (B circ K_(r-1)) B],
    K_(r-1)(a,l)=D_bl^(r-1)h_T(X_a).

The inherited directional-derivative matrix theorem gives fixed Lp
operator norms bounded by polylog(n). Flatness gives
||B circ K_(r-1)||_F<=polylog(n), so its displayed diagonal has averaged
squared L2 norm O(polylog(n)^2/n). This is where an averaged conclusion,
rather than an unproved uniform-root estimate, is deliberately used.

Vanishing mixed contractions make the star vector asymptotically
independent in polynomial moments of the joint collection (old X_i,
B X_T). Its odd psi mean is zero. Hence the first-chaos part contributes
zero to (U), proving both polynomial identities.

## 5. Input transfer and own-root handling

For nonlinear R_F, collision deletion and Lindeberg apply to its
low-influence forest kernels. The per-local-root repeated-marked-slot
Hilbert error is O(n^-1/2), because coefficients are O(n^-D/2) and
only O(n^(D-1)) tuples are forbidden. Transport by B preserves this
error in averaged Hilbert norm. The remaining squarefree kernel has
vanishing averaged maximum influence, by the proper-flattening lemma.
Joint Gaussian replacement and removal of the own output input follow.

The first-chaos unmarked test uses a different, stronger whole-functional
transfer. Put Y_i=(B X_T)_i and f_i=M(X_i)psi(Z_i). Y_i is exactly
squarefree and enters only linearly. After the usual squarefree star
replacement, f_i is a fixed polynomial of low-influence squarefree
fields. Its k-th repeated coordinate derivative has fixed Lp norm
O(n^-k/2), under all hybrid product laws. Fourth-order Lindeberg gives
only Y_i D_j^4 f_i and 4(D_jY_i)D_j^3 f_i because D_j^2Y_i=0.
The total influence bound and Cauchy--Schwarz over j give total error
O(n^-1)||Y_i||_2, and therefore o(1) after averaging roots. Fixed-degree
hybrid hypercontractivity controls the Taylor remainders. This does not
require low influences of Y_i, nor a Gaussian law for its coherent atoms.

For the marked identity the direct Boolean two-spin proof already
retains S_i exactly; there is no root-spin Gaussian replacement.

## 6. Bounded F and H; keep the channel finite

Let now F be a bounded odd Gaussian-a.e.-continuous response of finitely
many old coordinates, and H an even response in [0,1], with
|F(x)|+H(x)<=1 pointwise. Keep a FINITE set R of odd orders and define
f_r,s,h,Z,v as in Section 1 from this actual F. No infinite nonlinear
Hermite series is needed. Define K_F from its finitely many local
first-chaos coefficients as before.

Approximate F in Gaussian L2 by odd polynomials P. Under fixed L,

    n^-1 E||B(F(X)-P(X))||_2^2
       <= L^2 n^-1 sum_i E|F(X_i)-P(X_i)|^2.

Uniform root CLT and higher polynomial moments make the right side tend
to L^2||F-P||_2^2. The first-chaos isometry gives
||K_F-K_P||_2<=||F-P||_2, and the finitely many edge Hermite
coefficients converge too. Apply the arbitrary-channel polynomial
version, keeping h fixed while P improves. Thus bounded test versions
of (U),(S) pass to the actual F. Polynomial approximants are analysis
devices and are never used as infeasible Boolean means.

## 7. A careful bounded-test extension at zero channel variance

Approximating sign(z) uniformly in L2 over N(0,v), 0<=v<=L^2, is NOT
valid: as v decreases to zero any odd polynomial goes to zero while
sign(sqrt(v)N) still has unit norm. The proof must not use this shortcut.

Instead add an independent standard Gaussian auxiliary variable xi_i
and fix delta>0. Let Z_i^delta=Z_i+delta xi_i. Its limiting variance
u_i=v_i+delta^2 lies in the compact interval [delta^2,L^2+delta^2].
The density N(0,u_i) is bounded by
sqrt(L^2+delta^2)/delta times the density N(0,L^2+delta^2). Therefore
ordinary multivariate Gaussian L2 approximation controls bounded
joint tests of (old X,Z^delta) UNIFORMLY in i. Symmetrization preserves
old parity and channel parity; polynomials are finite sums of the
separable tests already proved. The matrix-side error is controlled
by averaged ||BF||_2 and the root CLT, exactly as in Section 6.

For these smoothed-channel tests the identities give the tested
regression

    S K_F(X) + s [v_i/(v_i+delta^2)] Z^delta.

No assertion of the remaining odd-old drift is needed: the permitted
test functions are even in old X. Equivalently this is the expectation
of S K_F(X)+s Z conditional on (X,S,Z+delta xi), with independent
Z=sqrt(v_i)N. The denominator is harmless because delta is fixed
throughout this approximation and its matrix-size limit.

## 8. Exact center update and removal of the auxiliary smoothing

For any actual bounded feasible F,H, choose the conditional spin means

    mu_plus  =  F(X) + H(X) sign(BF(X)),
    mu_minus = -F(X) + H(X) sign(BF(X)).

They lie in the cube because |F|+H<=1. Independent coordinate rounding
preserves both expected quadratic energies. Symmetry and hollowness give
the EXACT inequality

    Q(A)/(n sqrt(n-1))
      >= n^-1 sum_i E[H(X_i)|(BF(X))_i|].             (3)

There is no clipping loss proportional to ||B||op in this exact update.

For a lower bound on the right side, use the bounded odd soft sign
tau_eta(t)=max(-1,min(1,t/eta)), with eta>0 fixed. Test against

    H(X_i) tau_eta(S_i K_F(X_i)+s Z_i^delta).

Its part even in S_i is even in old X and odd in Z^delta; its part odd
in S_i is S_i times a function even in old X and in Z^delta. Thus it
is exactly within the bounded joint-test extension of Section 7.
After the matrix limit, the result is bounded below by

    n^-1 sum_i E H(X) |S K_F(X)+s sqrt(v_i)N|
       - eta - 2s delta E|xi| + o(1).                (4)

Indeed, for real y,e,

    y tau_eta(y+e) >= |y|-eta-2|e|.

Apply this with y=S K_F+s Z and e=s delta xi, and use 0<=H<=1.
This uniform regret inequality is why zero v_i cause no sign-limit
problem. First take n->infinity for fixed L,delta,eta and all finite
polynomial approximants. Then improve those approximants, and finally
let eta,delta decrease to zero. No discontinuous sign is justified by
an invalid uniform small-variance polynomial approximation.

## 9. Removing the matrix-dependent variances and the operator cap

For the chosen normalized odd Hermite mixture, the exact Schur
mean-standard-deviation inequality gives

    n^-1 sum_i sqrt(v_i) >= 1.

Its hypotheses are Q>=0, Q_ii=1, B symmetric with B^2=Q, and nonnegative
mixture weights f_r^2/s^2 summing to one. Those all hold here. For every
fixed real k, t -> E|k+s t N| is convex for t>=0, being an expectation
of absolute values of affine functions of t; it is also nondecreasing,
since its derivative is nonnegative for this symmetric Gaussian law.
Jensen and the mean-standard-deviation inequality therefore turn (4)
into

    liminf Q(A_n)/n^(3/2)
       >= E[H(X) E_N |K_F(X)+s N|]                  (5)

for EVERY fixed bounded-operator signing family, with the SAME right
side independent of L. The root spin was removed using the symmetry
of N and evenness of absolute value, not by replacing it by a Gaussian.

The diagonal-Grothendieck principal deletion and exact principal
monotonicity from the frozen minimal proof now extend (5) to arbitrary
signings. Fix the deletion fraction, apply the fixed-L theorem to the
retained growing orders, incur (1-epsilon)^(3/2), then let epsilon
decrease to zero. No parameter in the finite response or the sign
smoothing grows during a matrix-size limit.

Thus the finite-response center-update theorem is

    liminf_n M_n/n^(3/2)
      >= E[H(X) E_N |K_F(X)+s N|],

with bounded pointwise-feasible odd F, even 0<=H<=1 on a finite
ancestor-closed old family, and any fixed finite nonzero selection of
their nonlinear odd edge Hermite coefficients defining s. The independent
reconstructions checked the forest Wick replacement, uniform any-root
old/star independence, and full-contraction product formula. In particular,
the root-not-hit covariance-square argument first replaces the child Wick
polynomial by its pure top input chaos; it must not be applied directly
to an arbitrary sum of chaoses. Numerical evaluation is a separate task.
