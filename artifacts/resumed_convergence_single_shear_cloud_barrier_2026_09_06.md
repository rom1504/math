# One arbitrary controlled translation still cannot match independent clouds

Date: 2026-09-06. Status: a scoped discriminator for the new positive
controlled-index family. It does not bound its norm on correlated
Boolean profiles, rule out networks, or imply R<T.

## 1. Setup and the exact block decomposition

Let N=2^m and split the active Fourier indices as (a,c), where
c ranges over an L=2^k point binary group and a over N/L fibers.
Let Z=H_N X and W=H_N Y A be the independent Walsh clouds from
the majorant matching construction. The Boolean rows of X,Y
are independent samples of a symmetric law with covariance C,
and A^T C A=C. Thus every fixed collection of distinct cloud
rows converges jointly to independent N(0,C) vectors.

For a controlled translation sigma_f(a,c)=(a,c+f(a)), allow
also completely arbitrary signs at the matched frequencies.
Its best normalized squared matching defect is

    D_(N,L)^2=(1/N) sum_a min_(t in F_2^k)
               sum_c min_(epsilon_c=+/-1)
                   ||Z_(a,c+t)-epsilon_c W_(a,c)||^2.    (1)

The control truth table can be optimized independently in
each a-fiber; this is why (1) remains exact even though there
are L^(N/L) possible controlled translations.

## 2. Fixed fiber size has a positive Gaussian matching limit

For fixed L, let Z_1,...,Z_L,Y_1,...,Y_L be independent
N(0,C) vectors and define

    d_L(C)=E min_t (1/L) sum_c min_(epsilon_c=+/-1)
                        ||Z_(c+t)-epsilon_c Y_c||^2.     (2)

Then, as N tends to infinity through powers of two with L fixed,

    D_(N,L)^2 -> d_L(C) in probability.                   (3)

For the mean, every one-fiber joint Walsh CLT uses only 2L
rows. Its error is uniform in the fiber: orthogonality gives
block-diagonal covariance, and the bounded-summand
characteristic-function remainder is O_L(N^(-1/2)).
Uniform fourth moments make the quadratic minimum in (1)
uniformly integrable. For the variance, apply the same argument
to two distinct fibers, involving 4L rows. The limiting fibers
are independent; uniform eighth moments justify convergence
of their product moments. Averaging the fiber costs proves (3).

If C is nonzero, then

    d_L(C)>0.                                             (4)

Indeed, for each of the finitely many shifts and sign choices,
exact equality of every paired Gaussian vector has probability
zero on the nonzero range of C. The nonnegative finite minimum
is therefore strictly positive almost surely, and so is its
expectation. This includes singular nonzero covariance matrices.

For fixed feature dimension, d_L is continuous on the compact
set {C>=0: tr C=1}: couple all Gaussians using C^(1/2), and
use dominated convergence. Hence its minimum on that set is
positive for every fixed L. This is a positive normalized
defect, not merely failure of exact finite matching.

## 3. A simple explicit bound when no extra frequency signs are used

Remove the inner minimum over epsilon in (2), and call the
result d_L^unsigned(C). Conditional on the Y vectors, the
L cross inner products

    S_t=(1/L) sum_c Z_(c+t) dot Y_c

are centered Gaussians with the same variance

    V=(1/L^2) sum_c Y_c^T C Y_c.

The Gaussian exponential bound and Jensen give, for L>=2,

    E max_t S_t <=sqrt(2 log L) E sqrt(V)
                <=sqrt(2 log L/L) ||C||_F.

Since ||C||_F<=tr C and L is a power of two,

    d_L^unsigned(C)
       >=2[1-sqrt(log 2)] tr C >0.                       (5)

For L=1 the defect is exactly 2 tr C. Bound (5) is not
claimed with arbitrary signs; Section 2 supplies the signed
positive-gap statement at every fixed L.

## 4. Growing fibers return to the subexponential-family barrier

If L=L_N tends to infinity, the number of possible shears is

    #Sigma_(N,L)=L^(N/L),
    log #Sigma_(N,L)=N log(L)/L=o(N).                     (6)

The independently audited phase-robust cloud obstruction in
`resumed_director_phase_robust_cloud_obstruction_2026_09_06.md`
therefore applies: arbitrary Fourier signs have already been
optimized, and a subexponential number of index choices cannot
make the independent-cloud defect vanish.

Together, Sections 2 and 4 exclude vanishing independent-cloud
matching defect for any sequence of single controlled
translations with a fixed standard coordinate splitting:
every bounded-fiber subsequence has a fixed-L subsubsequence,
and every unbounded-fiber sequence has a subsubsequence with
L tending to infinity. No classification of the exact control
truth tables is needed.

## 5. Precise consequence for the positive realization theorem

The arbitrary-control promotion genuinely escapes the earlier
quadratic inverse-pencil counting bound: scalar controls already
have exponential cardinality in N. But those choices only
swap pairs independently; they do not provide arbitrary
Gaussian point-cloud transport.

This result does NOT exclude a deliberately correlated Boolean
profile construction for the same shears, a compatible network
of many shears, other high-capacity index families, or a
different proof of R=T. It is a test of the independent-cloud
landing strategy, not an original-problem upper or lower bound.
