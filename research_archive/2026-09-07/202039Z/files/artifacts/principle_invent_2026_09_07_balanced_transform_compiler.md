# Exact balanced-column repair and a seed/bulk compiler checkpoint

2026-09-07. **Balanced-transform theorem proved; mixed-fibre cap transfer
remains open.** This is a positive actual-sign operation investigated
after the hierarchy-complexity theorem. It is not a claim that a
factorization witness alone solves flatification or seed recovery.

## 1. Exact repair with a quantitative operator error

Let H be a q by m sign matrix, q even, with ||H||op<=sqrt(M), where
q,m,M are comparable. Put P=I-J/q. Suppose every column sum s_j obeys
|s_j|<=S. Independently for each column, flip exactly |s_j|/2 entries
chosen uniformly from its majority sign. The resulting matrix T is
again FULL signs and every column sums EXACTLY to zero.

If S=o(M/log^2 M), then there are realizations with

    ||T-PH||op
       <=sqrt(M)S/q+O(sqrt(S)log M)=o(sqrt(M)).      (1)

For the customary S=O(sqrt(M log M)), the right side is
O(M^(1/4)log^(5/4)M). In particular the operation can remove the physical
constant mode without replacing sign entries by fractional projections.

Here is the full calculation. For a column h with sum s, direct sampling
gives

    E T_j=d_j Ph,  d_j=q/(q+|s|).                  (2)

For example if s>0, a positive entry is flipped with probability
s/(q+s), so E T_j=q h/(q+s)-s1/(q+s), which is (2). The negative case
is identical after sign reversal. Thus

    ||E T-PH||op<=||PH||op max_j(1-d_j)
                 <=sqrt(M)S/q.                    (3)

Let xi_j=T_j-E T_j. Flipping |s|/2 entries changes h by a vector of
squared norm2|s|, so ||xi_j||<=2sqrt(2S). The covariance of a uniform
fixed-size subset of its majority positions is a scalar multiple of
the projection perpendicular to the all-ones vector on those positions.
Explicitly its nonzero eigenvalue before the factor4 from flipping is
p(1-p)N/(N-1), where N=(q+|s|)/2 and p=|s|/(q+|s|). Consequently

    E xi_j xi_j^T <= (8S/q)I,
    E||xi_j||^2<=2S.                               (4)

For independent rectangular summands Z_j=xi_j e_j^T, their operator
bound is2sqrt(2S), and their rectangular matrix-variance parameter is
at most max(8mS/q,2S)=O(S). The rectangular matrix Bernstein inequality
therefore gives

    ||sum_j Z_j||op=O(sqrt(S log M)+sqrt(S)log M)

with probability tending to one. The stated hypotheses are exactly
those of [Tropp, Theorem1.6](https://arxiv.org/pdf/1004.4389), which I
read directly: independent mean-zero rectangular summands, a uniform
operator bound, and the maximum of the two covariance sums. Combining
this with (3) proves (1).

For a partial Hadamard, first dephase a selected column to all ones by
a physical row gauge and omit that column. The remaining full-order
columns are balanced before row selection. If their selected column
excesses satisfy the displayed S bound, (1) applies; the omitted DC
column projects to zero. In this case HH^T=MI-J and

    TT^T=MP+o_op(M).                               (5)

The column-excess hypothesis is a deterministic part of the theorem,
not silently asserted for arbitrary physical selectors.

## 2. Exact transfer through a signed weave

For each fibre i let H_i be a q by m row transform. A masked signed
weave can be written W(H)=F_H^T J F_H, where F_H is the block-diagonal
map x_i ->H_i^T x_i, and J is the signed swap of the directed edge
coordinates, set to zero on any deleted edges. Thus ||J||op<=1.

Apply the repair independently in every fibre, with a uniform error
bound (1), and write P_all for the projection off every physical fibre
constant. Then

    ||W(T)-P_all W(H)P_all||op
       <=(max_i||T_i||+max_i||PH_i||)
                          max_i||T_i-PH_i||
        =o(M).                                    (6)

This estimate is uniform over the edge mask and signed-swap phases.
It gives an o(n^(3/2)) energy error when n=mq and m,q,M are comparable.
The repaired weave has EXACTLY zero row and column sums in every
remaining cross-fibre block. Its constant-fibre subspace is therefore
invariant and receives zero energy from those blocks.

An edge mask need not be claimed negligible in cap. Instead, if the
deleted edge phases were independent fair signs in the original
construction, convexity gives the exact conditional inequality

    Q(W_remaining)<=E_deleted Q(W_remaining+W_deleted). (7)

The same applies to independent diagonal-fibre signs after physical
hollowing. Averaging (7) over the remaining randomness preserves any
proved expected-cap upper bound for the original signing. It does not
assert that every realization of an arbitrary deletion has small cap.

For physical words BALANCED WITHIN EVERY FIBRE, P_all x=x. Thus (6)--(7)
show that any expected-cap certificate for the original weave transfers,
up to o(n^(3/2)), to the repaired masked weave on this balanced Boolean
face. In particular the actual strict-subhalf certificate is retained
on that face whenever the selector/excess event and repair are included
with their vanishing conditioning cost.

Deleted off-diagonal blocks may now be filled by signed J_q blocks,
which act only on the constant-fibre subspace. This is an actual way to
insert a chosen sparse seed code while keeping the known subhalf bulk
certificate on its centered Boolean face. Diagonal fibres may be filled
by signed J_q-I_q blocks, with their constant and bounded diagonal
corrections explicitly retained. Every physical off-diagonal entry of
the resulting matrix is a sign.

## 3. The precise unresolved mixed response

This operation does NOT yet transfer an arbitrary seed cap. For a general
Boolean fibre word of mean a, its centered component has coordinates
1-a and -1-a. Those are not in the balanced Boolean cube, nor even in
[-1,1] when a is nonzero. Therefore projecting the old weave and invoking
its ordinary Boolean cap would be invalid.

The missing quantity is a concrete constrained response of the repaired
bulk at specified fibre means, not a renamed unconstrained control law:

    max {|H_bulk(x)| : x_i Boolean, mean(x_i)=a_i}.

Its values at a_i=0 are controlled by the transferred strict certificate;
its values when every |a_i|=1 are zero. A norm bound gives a linear
penalty in sum(1-a_i^2), but that relaxation has the already diagnosed
spectral half floor and cannot recover a strict-subhalf seed. A successful
continuation must prove a genuinely nonlinear Boolean bound between
these faces, with the seed energy included, or find an actual mixed-word
counterexample for the proposed completion.

There is no entropy-free convex decomposition of a nonbalanced Boolean
word into balanced Boolean words and constant words: it is an extreme
point of the full cube. This elementary obstruction is why knowing the
two endpoint faces alone does not finish the compiler.

## 4. Relation to the hierarchy witness

The hierarchy theorem forces a sufficiently rich near-level phase code
in liminf minimizers. Averaging gauges from such a code produces
A Hadamard-product E(xx^T), but independently rounding that weighted
matrix pays a square-root variance cost. The factorization lower bound
alone has not been shown to remove that cost. The balanced-column
operation above instead supplies an exact physical invariant mode and
a controlled actual-sign error on its complement. Whether a near-level
phase witness can choose its mixed response successfully is still open.

This checkpoint is preserved as a proved operation and an explicit
remaining test, not counted as another cap theorem or convergence proof.
