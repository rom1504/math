# Independent audit: sparse restrictions of actual bounded-cap parents

Date: 2026-09-06. Status: the theorem below is verified by reconstruction;
the alternate calculation at the end is an unfinished research note.

Throughout this note Q(A)=max_x |x^T A x|/2. Some older spectral-deletion
files call the quantity without the factor 1/2 Q; their constants must not
be copied without conversion.

## Verified statement and dependency map

For every sequence of hollow symmetric sign matrices with
Q(A_D)=O(D^(3/2)), and every simultaneous n->infinity, n/D->0, a uniformly
random exact-size principal restriction satisfies

    Q((A_D)_T)/n^(3/2) >= 2/pi - o_P(1).

This concerns the actual Boolean cap, not just a failed certificate. It
applies to actual minimizing parents. It says nothing about exceptional
deterministic selectors, and is not a lower bound on M_n.

The reconstruction uses the following mathematical steps, not earlier audit
verdicts.

1. **Spectral core from the actual cap.** Reconstruct Section 1 of
   `resumed_bound_audit_minimal_proof_2026_09_06.md`. Cube polarization gives
   beta(A)<=4Q(A). Odd-tensor Gaussian rounding proves the real
   Grothendieck bound with K=pi/(2 asinh(1)). The SDP simultaneous majorant
   diag(d)>=+/-A has dual X,Y>=0, diag(X+Y)=1. Gram vectors
   (a_i,b_i),(a_i,-b_i) are unit, so its value is <=K beta(A), with no
   extra factor two. Strict finite SDP feasibility gives a nonnegative
   majorant of this trace. Removing d_i>K beta(A)/(epsilon D) leaves a
   (1-epsilon)D core with operator norm <=4K Q(A)/(epsilon D).

2. **Uniform local moments on that core.** The primary imported theorem is
   [Mingo--Speicher, Theorem 6, pp. 6--7](https://arxiv.org/pdf/0909.4277).
   Its graph sum is over unrestricted colorings; loops and multiple edges
   are permitted. An Eulerian component has no cut edge, so its exponent
   in that theorem is one. Thus an Eulerian graph with v vertices, e
   edges and c components has normalized sum at most
   D^(c-v)||B||op^e. Complete the hollow parent as B=A+I before parity
   cancellation. For b=v-c-e/2, this is O_C(D^(-b)).

   In a closed-walk diagram let a=v-1-k/2. Edges of positive even
   multiplicity connect the parity components, proving a<=b. Positive
   excess is O_C((n/D)^a), negative excess is O(n^a), and nonempty
   zero-excess parity graphs have a nonloop edge. Conditioning on all
   other colors leaves B_ij f(i)g(j), bounded by beta(B)/D^2=O_C(D^-1/2).
   The parity graph has no multiple edges after cancellation; hence this
   conditioning does leave exactly one edge between those two variables.

   A partition identifying h ambient colors contributes exactly
   (n/D)^h times the bounded quotient-walk expression, including loops.
   This is why no n^2/D collision hypothesis is needed. Doubled trees are
   the only leading diagrams. For two walks with a shared root, their
   doubled-tree edge sets cannot overlap and meet only at the root.
   Consequently the leading rooted second moment factors. This gives
   E avg_i[((A_T/sqrt n)^k)_ii-m_k]^2
   =O_{k,C}(n/D+1/n+D^-1/2).

3. **Actual Boolean witness.** With semicircle orthogonal polynomials
   V_0=1,V_1=x,V_(j+1)=xV_j-V_(j-1), take
   R=I+(sum_(j=0)^q V_j(L))^2. At FIXED q, local moments give
   tr(R)/n->q+2, tr(LR)/n->2q, diagonal variance->0, and bounded
   tr(R^2)/n. The finite covariance inequality in
   `transfer_seed_finite_spectral_gaussian_cap_2026_09_06.md` follows
   from Gaussian signs, |arcsin z-z|<=|z|^3, and the flat row norms of
   L=A_T/sqrt n. It yields 2q/[pi(q+2)]-o_P(1). No operator bound for
   the sampled matrix or growing-degree limit is assumed.

4. **Remove the core deletion.** Conditional on m=|T intersect core|,
   the intersection is uniform in the core. Hypergeometric variance
   <=n/4 gives m/n>=1-epsilon-o_P(1). The fixed-degree estimates are
   uniform over the likely m range. Principal monotonicity follows by
   averaging omitted unbiased spins, without any bridge penalty.
   The resulting bound is (1-epsilon)^(3/2)2q/[pi(q+2)]-o_P(1).
   For a desired fixed accuracy choose epsilon and q first, then take
   the order limit. This proves the assertion for all bounded-cap parents.

The full proof is
`transfer_seed_sparse_bounded_cap_parent_theorem_2026_09_06.md`; its
elementary Hadamard predecessor is
`transfer_seed_sparse_hadamard_local_moments_2026_09_06.md`.

## Why the earlier free-projection theorem was not enough

I separately checked [Anderson--Farrell, *Asymptotically liberating
sequences of random unitary matrices*, Corollaries 3.5 and 3.7](https://arxiv.org/pdf/1302.5688v4).
The common-permutation formulation handles exact-size principal sampling
of a Hadamard conjugation. The diagonal formulation permits the two
diagonal random variables to be dependent, including identical variables.
It gives fixed-retention free-projection limits. It does not by itself
justify simultaneous retention p_D->0 or the needed diagonal concentration.

For symmetric U=H/sqrt D and P a projection, Q=UPU, even compressed
moments can be expressed using Tr(PQ)^k. At fixed retention p, their
normalized limiting first even moments are

    1, 2-p, 5-6p+2p^2, 14-28p+20p^2-5p^3.

Sending p to zero after this fixed-p theorem suggests the semicircle,
but exchanging those limits would be a gap. The graph proof above removes
that gap rather than assuming uniformity of the imported result.

## Unfinished alternate calculation, preserved rather than discarded

An independent attempted route centered the sampling projection:
P0=P-pI. In an expansion of Tr(P0 U P0 U ... P0 U) with 2k factors,
singleton diagonal colors have zero expectation under independent
Bernoulli sampling, and surviving patterns use at most k colors.
Flat Hadamard entries suggest an unnormalized O_k(p^k) bound when pD
is large. Expanding back through the two-projection algebra might then
recover a simultaneous error O_k(1/(pD)). I did not complete the
coefficient bookkeeping, exact-size correction, odd moments, or local
diagonal estimate in this alternate argument. It is NOT an additional
proof and is superseded operationally by the verified graph argument.

## Original question still open

The theorem eliminates a random-selection transfer, not all-order
recovery. Any useful principal-selector route must use exceptional,
globally correlated selections. There is presently no verified estimate
showing such selections exist at enough target orders with normalized cap
approaching the parent's cap. Nor has their impossibility been proved.
