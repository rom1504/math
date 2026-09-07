# Rank-one variance completion at small Boolean-cap cost

Date: 2026-09-07. Status: **proved; independently audited PASS by the
director**. The audit checked the aggregate mass, product-deficit row
error, simultaneous concentration, and separate small-cap estimate.
This strengthens the director's variance completion lemma by removing
the original bounded-entry hypothesis from the COMPLETION operation.
The output is weighted, not a sign matrix.

## 1. Statement

Let W_N be real, symmetric, and hollow, with uniform row square sums
N-1+o(N), and Q(W_N)<=C N^(3/2), where

    H_W(x)=x^T W x/2,
    Q(W)=max_(x Boolean) |H_W(x)|,
    width(W)=(max H_W-min H_W)/2.

For every sufficiently small fixed epsilon>0, there is a principal set
I of size k>=(1-epsilon)N and a symmetric hollow real matrix V_k with

    ||V_k||op <= O_C,epsilon(sqrt(N)),
    max_i |sum_j V_ij^2-(k-1)| = O_epsilon(sqrt(N log N)),
    width(V_k) <= (1+O(epsilon)+o(1)) width(W_N)
                         +O(sqrt(epsilon)N^(3/2)).       (1)

The same estimate holds for Q. The completion entries are bounded by
O(epsilon^(-1/2)), regardless of the amplitudes of W. In particular,
if max|W_ij|/sqrt(N)->0, then max|V_ij|/sqrt(k)->0 as well.

No claim is made here that the marked lower theorem is already proved
for this larger delocalized class. The regularization implication is
separated from the unresolved fixed-operator full-overlap diagrams.

## 2. Retained mass is controlled in aggregate

Use the same simultaneous Grothendieck diagonal-majorant deletion as in
`principle_director_variance_completion_2026_09_07.md`. It gives I of
size k=N-r, r<=epsilon N, whose principal matrix W_0 satisfies

    ||W_0||op <= O_C(epsilon^(-1) sqrt(N)).

Let q_i=sum_(j in I) W_ij^2, and q_+=max_(i in I) q_i. Certainly
q_+<=N+o(N). More importantly, deleting r rows and the corresponding
columns removes at most the sum of their row and column square masses:

    sum_(i in I) q_i >= N(N-1)-2r(N-1)+o(N^2).

Thus q_+>=N-O(epsilon N)+o(N). These aggregate estimates require no
individual entry bound, and remain valid for arbitrary variance profiles.

Set

    gamma^2=(k-1)/(q_+ + epsilon N),
    d_i=(k-1)-gamma^2 q_i,
    S=sum_(i in I) d_i.

Then gamma=1+O(epsilon)+o(1), and, for all sufficiently large N at fixed
epsilon,

    c epsilon N <= d_i <= k-1,
    c epsilon N^2 <= S <= C_0 epsilon N^2,             (2)

where c,C_0 are absolute for small epsilon. The lower bounds come from
the common epsilon N added in gamma's denominator. For the upper bound
on S, insert the retained-mass estimate and gamma^2=1+O(epsilon)+o(1).
The original uniform o(N) row error is absorbed at fixed epsilon.

## 3. Product deficits give an approximately exact hollow completion

For distinct i,j put

    u_ij=d_i d_j/S,        u_ii=0.

Then u is symmetric and nonnegative, and

    sum_j u_ij=d_i-d_i^2/S,
    max_ij u_ij=O(1/epsilon),
    max_i d_i^2/S=O(1/epsilon),
    sum_(i<j) u_ij <= S/2=O(epsilon N^2).             (3)

Exact row sums are unnecessary: the missing hollow diagonal contribution
is only O(1/epsilon), uniformly, which is negligible compared with N at
fixed epsilon. This avoids an exact fractional-degree realization problem.

Choose independent signs on unordered pairs and put

    R_ij=sign_ij sqrt(u_ij),        V=gamma W_0+R.

The diagonal remains zero. The row square sum is exactly

    sum_j V_ij^2=(k-1)-d_i^2/S
                   +2 gamma sum_j W_ij R_ij.         (4)

The cross term has subgaussian variance proxy at most
(max u_ij)q_i=O(N/epsilon), independent of the largest individual W entry.
A union bound over rows gives the row-error statement of (1) with
probability tending to one.

## 4. Separate operator control from objective control

The operator norm estimate uses U=max u_ij=O(1/epsilon). For any real
unit vector x, x^T R x has subgaussian variance proxy at most 2U. A
1/4-net of the sphere of cardinality at most 9^k gives

    ||R||op <= O(sqrt(k/epsilon))

with probability at least 1-exp(-c_1 k), for a suitable fixed constant.
Together with the principal-core bound this proves the required fixed
normalized operator bound for V at fixed epsilon.

The operator bound alone is NOT small enough for the objective. For a
fixed Boolean x, however, H_R(x) has variance proxy

    sum_(i<j) u_ij <= S/2.

The elementary sign exponential-moment bound gives

    Pr(|H_R(x)|>=t) <= 2 exp(-t^2/S).

Union bound over 2^k Boolean vectors, with t=A sqrt(kS), gives

    Q(R)<=A sqrt(kS)=O(sqrt(epsilon)N^(3/2))          (5)

with probability at least 1-2 exp(-(A^2-log 2)k). Choose A fixed and
large. The row, operator, and objective events have nonempty simultaneous
intersection. Fix any realization in that intersection.

Half-width and cap obey the triangle inequality, and principal averaging
does not increase either objective. Therefore

    width(V)<=gamma width(W_0)+width(R)
             <=gamma width(W_N)+Q(R),

which proves (1). The corresponding Q inequality is identical.

## 5. Consequences and exact boundary

This supplies a fixed-operator core with approximately restored row
variance and small Boolean objective cost using only aggregate lost
variance. Large deficits at a few rows need not be individually small.
The noise can have normalized operator norm O(epsilon^(-1/2)) while its
normalized Boolean cap is O(sqrt(epsilon)); these are different estimates
and are deliberately not conflated.

For the already established bounded-amplitude weighted lower theorem,
the lemma offers an alternative regularization proof, since its added
entries are bounded at every fixed epsilon. For the prospective merely
delocalized class, it removes spectral regularization as an additional
obligation: only the fixed-operator old-frame/full-overlap extension
would remain.

It does not include localized fixed-size blocks scaled by sqrt(N) in
that prospective class, because their normalized maximum entry does not
tend to zero. It does not flatten weights to signs or identify optimized
cap values between different orders.
