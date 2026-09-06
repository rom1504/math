# A sparse-flip lower bound on the number of near-active states

Date: 2026-09-06. Status: exact finite theorem derived from global edge
optimality. It is a genuine critical-window constraint on minimizing
signings, but it does not give the joint geometry needed for insertion
or imply convergence. No exhaustive novelty claim is made.

Let A be an exact order-n minimizing hollow signing, with cap M=M_n.
Represent the absolute objective by oriented projective spins v, where

    v_(ij)=sigma x_i x_j, sigma in {+1,-1}, x modulo global negation.

The oriented energy is a.v<=M and its gap is g(v)=M-a.v. Let N_A(W)
count these distinct states with gap at most W.

## 1. Exact theorem

For an integer k>=1 put

    D_k=ceil((k+1)^((k+1)/2)).

If

    M>(2^k+2)D_k,

then

    N_A(2D_k)>k.                                      (1)

The proof only uses the finite signed-feature minimax structure, not
a Gaussian limit, bounded operator norm, or a hypothesized common
ground-state law. The quadratic problem supplies its growing cap.

## 2. An elementary bounded integer positive-combination lemma

Let S be a subset of {+1,-1}^l such that conv(S) contains a vector whose
coordinates are all strictly positive. There are nonnegative integers
t_z, z in S, such that

    1<=sum_z t_z<=D_l,
    sum_z t_z z_j>=1 for every j=1,...,l.              (2)

To prove this, maximize delta subject to

    p_z>=0, sum_z p_z=1, sum_z p_z z_j>=delta.

The optimum is positive. Choose an optimal basic solution. If delta=1,
the all-one vector occurs in S and proves the claim immediately. In
the other case, if r probabilities are positive, r<=l, and r independent
tight coordinate constraints together with sum p=1 determine those
r probabilities and delta. Their coefficient matrix has integer entries;
the r coordinate rows have Euclidean norm sqrt(r+1), and the final row
has norm sqrt(r). Its nonzero determinant q has absolute value at most

    sqrt(r)(r+1)^(r/2)<=D_l

by Hadamard's determinant inequality. Cramer's rule makes q p_z and
q delta integers. Take t_z=q p_z, using the absolute denominator and
the actual rational solution. Their sum is q, and q delta is a positive
integer. This proves (2). Reducing denominators could only improve the
bound; it is not needed.

## 3. Proof by a bounded edge flip

Suppose instead that at most k states have gap <=2D_k, and enumerate
them as v^(1),...,v^(l), l<=k. There is at least one exact ground state,
so l>=1. To each edge e assign the alignment pattern

    z_e=(a_e v_e^(1),...,a_e v_e^(l)) in {+1,-1}^l.

Discard every pattern occurring fewer than D_k times. Fewer than
2^l D_k<=2^k D_k edges are discarded. For every coordinate j, the sum
of all alignment signs is M-g(v^(j))>=M-2D_k. After discarding the rare
patterns, its sum is still strictly positive by the assumed bound on M.
Thus the convex hull of the remaining, abundant patterns meets the
strictly positive orthant.

Apply (2). It selects a multiset of at most D_l<=D_k abundant patterns
whose coordinate sums are at least one. Each pattern is available on
at least D_k distinct edges, so this multiset can be realized by an
ACTUAL subset F of distinct edges, |F|<=D_k. Flip exactly those signs.

For every listed near-active state its oriented energy changes by

    -2 sum_(e in F) a_e v_e^(j)<=-2.

Every other oriented state had gap >2D_k, and any |F|-edge flip can
increase its energy by at most 2|F|<=2D_k. It therefore also remains
strictly below M. The new signing has cap strictly less than M, a
contradiction to exact global minimality. This proves (1).

Notice why independent edge witnesses alone were insufficient: the
bounded integer combination simultaneously defeats the ENTIRE assumed
small near-active family. Abundant pattern multiplicities are what
make its fractional positive combination a genuine edge subset.

## 4. Consequence at the insertion-relevant scale

For every fixed epsilon>0, choose

    k=floor((1-epsilon)log(n)/log(log(n))).

Then

    log D_k=(1/2+o(1))k log k,
    2D_k=n^((1-epsilon)/2+o(1))=o(sqrt(n)),
    (2^k+2)D_k=n^((1-epsilon)/2+o(1)).

Even the elementary variance lower bound
M_n>=sqrt(binomial(n,2)) suffices to make the hypothesis of (1) hold
for all large n. Hence, for every fixed C>0,

    N_A(C sqrt(n)) >= (1-o(1)) log(n)/log(log(n))        (3)

uniformly over exact order-n minimizers. The stronger current universal
n^(3/2) lower bound is not needed for this critical-window consequence.

The count includes oriented states and identifies global spin negation.
Without orientation, the count is the same in any window W<M, because
only one of the two orientations can then be near-active.

## 5. Scope and remaining limitation

The theorem requires EXACT minimization. An asymptotic near-minimizer
can tolerate the bounded decrease produced here, so the statement does
not automatically survive an o(n^(3/2)) cap error or regularization.

The bound is logarithmic, not exponential, and it gives no control of
overlap, Gaussian width, or discrepancy of the near-active spins. In
particular it supplies neither a good bridge signing nor the desired
fixed-block cap increment. It removes only the possibility that an
exact minimizer's critical window is eventually a bounded list of
exceptional witnesses.

A bounded archive search found the earlier randomized repair entropy
and separate-edge stationarity arguments, but no determinant/abundant-
pattern version of this finite-window lemma. This is a provenance note,
not a claim of newness in the external discrepancy literature.
