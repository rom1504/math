# Arbitrary frequency phases do not rescue a small family of independent-cloud matchings

Date: 2026-09-06. Status: director derivation; the bound-audit researcher
independently reconstructed the concentration argument and Gaussian gap.
This is a scalable obstruction to independent-cloud realization, not a
bound on jointly designed profiles, R, or the original signing minimum.

## 1. Exact theorem

Fix a centrally symmetric law X on {+1,-1}^k with covariance C, and a
fixed real matrix A satisfying A^T C A=C. For N=2^m take independent
samples X_1,...,X_N and Y_1,...,Y_N from this law. Put their rows in
matrices F,G and let H be normalized Walsh. For a permutation sigma of
the N Fourier rows define

    Delta_sigma=min_(s_i in {+1,-1})
        N^(-1/2)||P_sigma H F-diag(s) H G A||_F.       (1)

The signs are chosen after seeing BOTH complete samples. Any complex
unit phases give the same minimum here, since the row inner products
are real. Let Sigma_N be any nonempty deterministic permutation family
with log|Sigma_N|=o(N). Then

    min_(sigma in Sigma_N) Delta_sigma^2
        -> d(C):=2tr(C)-2E|Z dot Y| >0                (2)

in probability, where Z,Y are independent centered Gaussians of
covariance C. In particular

    d(C)>=2(1-2/pi)tr(C)=2(1-2/pi)k.                 (3)

Thus even exponentially many freely chosen common frequency signs do
not remove the positive defect of a subexponential permutation family.
For every fixed epsilon<d(C), the probability that the left side of (2)
is at most epsilon is bounded by exp[-cN+o(N)], for some c>0 depending
only on the fixed law, A, and epsilon.

## 2. The expected squared defect is permutation-independent

Rowwise elimination of the signs gives

    Delta_sigma^2=(1/N)sum_i
       (||Z_sigma(i)||^2+||W_i||^2-2|Z_sigma(i) dot W_i|),
    Z=HF, W=HGA.                                   (4)

Central symmetry makes the marginal law of EVERY Walsh row the same
normalized iid sum. The F and G samples are independent. Hence the
expectation of (4) is exactly independent of sigma. The ordinary
bounded-vector CLT sends the pair of rows to independent N(0,C)
vectors, because A^T C A=C.

The cross inner products are uniformly integrable: their second moment
is exactly tr(C^2). Also E||Z_i||^2=E||W_i||^2=tr(C)=k. Consequently

    E Delta_sigma^2 -> d(C), uniformly in sigma.      (5)

To prove (3), diagonalize C with nonnegative eigenvalues lambda_j.
Then Z dot Y has the law sum_j lambda_j g_j h_j, with independent
standard Gaussian pairs. The triangle inequality and
E|g_j h_j|=2/pi give E|Z dot Y|<=(2/pi)sum_j lambda_j.
This bound is sharp when C has rank one. It is not an independence
claim about ALL transformed rows.

## 3. Product concentration after optimizing all phases

Set a=||A||op. All Delta_sigma are bounded by L=sqrt(k)(1+a). If two
input sample arrays differ in h of their 2N independent rows, then
Walsh orthogonality and the Lipschitz property of distance to a set give

    |Delta_sigma(z)-Delta_sigma(z')|
       <=C_0 sqrt(h/N),  C_0=2sqrt(k(1+a^2)).        (6)

This bound already minimizes over every sign choice. It is uniform
in sigma and uses only bounded input rows, not transformed independence.

Here is an elementary concentration proof adequate for the exponential
union bound. For any event E with probability at least 1/2 in a product
of R=2N independent row spaces, its Hamming distance D_E has bounded
differences one. McDiarmid's lower-tail bound at zero yields

    1/2<=P(D_E=0)<=exp[-2(E D_E)^2/R],
    E D_E<=sqrt(R log(2)/2)=sqrt(N log(2)).

The upper-tail bound then gives

    P(D_E>=r)<=exp[-(r-sqrt(N log(2)))_+^2/N].

Apply this to each of the two median half-events of Delta_sigma. If
its value differs from its median by at least t, (6) forces Hamming
distance at least Nt^2/C_0^2 from the opposite half-event. Thus

    P(|Delta_sigma-med Delta_sigma|>=t)
      <=2 exp[-(sqrt(N)t^2/C_0^2-sqrt(log(2)))_+^2]. (7)

In particular this is exp(-c_t N) for every fixed t>0. Together with
the uniform bound L, (7) implies

    E Delta_sigma^2-(med Delta_sigma)^2 ->0

uniformly over all permutations. Equation (5) therefore places every
median at sqrt(d(C))+o(1), uniformly. The exponential lower tail and
log|Sigma_N|=o(N) prove the lower side of (2); one fixed member and
(5)--(7) prove its upper side. This proves the theorem.

## 4. Application and precise scope

The quadratic dual-index counting theorem supplies exp(O(m^4))
permutations in its injective-polar version and exp(O(m^6)) in its
broader affine/function-injective version. Both are exp(o(2^m)).
Theorem (2) therefore remains applicable even if EVERY index map is
given an arbitrary data-dependent common Fourier phase pattern. No
separate count of Arf signs or their representations is required.

For exactly balanced profiles on a doubled Boolean cube, whose Fourier
support is one affine hyperplane, permutations preserving that active
hyperplane reduce to the theorem with N equal to its size. The inactive
extension is irrelevant. Permutations which move active rows to inactive
rows can only reduce the expected cross term in (4); the positive lower
bound and the same concentration scale still apply, though the exact
limit d(C) need not be the same for all such permutations.

The arbitrary-matching norm theorem is not contradicted: the number of
all permutations is exp(Theta(N log N)), outside the hypothesis. The
bilinear-compatible cloud-matching theorem is not contradicted either:
its large compatible family is not the much smaller realizable quadratic
index family.

The result is specifically about INDEPENDENT bounded profile samples.
It leaves deliberately correlated or jointly constructed profiles,
nonquadratic index families of sufficient richness, and genuinely growing
block rotations open. It does not permit independent phases for each
feature column; the common row phase is the one supplied by a signed
profile module. No assertion R<T or c_liminf<c_limsup follows.
