# Boolean-eigenbasis Hamming rigidity for actual signings

Date: 2026-09-06. Status: exact independently reconstructed theorem.
This extracts and sharpens the scope of Section 5 of
`resumed_bound_audit_original_lower_sampling_and_optimizer_laws_2026_09_06.md`.
It is a structural necessary condition for a sub-1/2 minimizing sequence,
not a convergence theorem or a critical-window ground-state law.

## 1. Exact theorem

Let H be a symmetric n-by-n matrix with entries in {+1,-1} and H^2=nI.
Assume there are n pairwise orthogonal Boolean eigenvectors x^(j), with

    H x^(j)=sigma_j sqrt(n) x^(j),  sigma_j in {+1,-1}.

Let A be any hollow symmetric signing, and let d count the UNORDERED
off-diagonal edges where A and H differ. With

    H_A(x)=sum_(i<j) a_ij x_i x_j,
    Q(A)=max_(x Boolean) |H_A(x)|,

one has

    Q(A) >= [n(n-1)/2-2d]/sqrt(n),                    (1)
    Q(A)/n^(3/2) >= 1/2-1/(2n)-2d/n^2.              (2)

Proof: the matrix with columns x^(j)/sqrt(n) is orthogonal. Therefore

    (1/n) sum_j sigma_j x^(j)(x^(j))^T=H/sqrt(n).

Choose the oriented Boolean spin (sigma_j,x^(j)) uniformly. Its expected
oriented A-energy is

    (1/sqrt(n)) sum_(i<j) a_ij H_ij
      =[n(n-1)/2-2d]/sqrt(n).

Every oriented energy is at most Q(A), proving (1). The argument does
not assert that this law is supported on A-ground states.

In particular, if Q(A)/n^(3/2)<=1/2-delta, then

    d >= (delta/2)n^2-n/4.                          (3)

Thus a sequence of actual minimizers staying delta below 1/2 cannot be
o(n^2) edge edits from this Hadamard class. The same conclusion holds
for any signing sequence with that cap, independently of exact optimality.

## 2. A nonempty infinite class and its exact cap window

For every a>=1, let n=4^a and

    H=(J_4-2I_4)^(tensor a).

The four Walsh vectors diagonalize J_4-2I_4: the constant vector has
eigenvalue +2, and the three balanced vectors have eigenvalue -2. Their
tensor products give the required n orthogonal Boolean eigenvectors.
Simultaneously permuting rows and columns, or switching H to DHD with
D a diagonal sign matrix, preserves all hypotheses.

For the displayed tensor H, its diagonal is the constant delta=(-1)^a.
The hollow signing A=H-delta I has the exact cap

    Q(A)=(n sqrt(n)+n)/2.

Indeed its spectral norm is sqrt(n)+1, and a Boolean eigenvector with
sigma=-delta attains the resulting quadratic upper bound. For any
Boolean eigenvector, its naturally oriented energy is

    sigma H_A(x)=(n sqrt(n)-sigma delta n)/2.

Hence the two orientations of this eigenbasis law have gaps 0 and n
from Q(A). A positive fraction have the latter gap. The law proves the
1/2 limiting scale but is not concentrated in an O(sqrt(n)) cap window.
This is why (1) does not close the fixed-block insertion obligation.

## 3. Scope checks

The hypothesis is a COMPLETE ORTHOGONAL Boolean eigenbasis, not merely
one Boolean eigenvector, a symmetric Hadamard matrix, or a bounded
normalized operator norm. No claim is made that every symmetric
Hadamard or conference matrix possesses such a basis.

The edge count d in (1) excludes the diagonal and counts each symmetric
edge once. Counting both matrix entries would replace the term 2d by
the corresponding ordered-entry mismatch count. The lower bound is
allowed to become negative when d is large; its content is near the
displayed structured class.

The result blocks a particular proposed o(n^2)-edit recovery of a
genuinely sub-1/2 minimizing subsequence. It neither establishes such
a subsequence nor rules out a recovery by larger edits, a different
structured family, or an amplification that changes the order.
