# Same-order repair: linearly balanced near-minimizers

2026-09-06. Director theorem; independently verified in
`transfer_adversary_same_order_orientation_repair_audit_2026_09_06.md`.
This strengthens the earlier insertion repair. It changes actual edge signs
at the SAME order, not a weighted relaxation, and preserves the original
normalized optimum up to a power-saving error. It does not prove convergence.

Write `P(A)=max H_A`, `R(A)=max(-H_A)`, `Q(A)=max(P,R)` and
`K=pi/(2 asinh(1))`. The norm convention is `H_A(x)=x^T A x/2`.

## 1. A finite repair lemma

Reverse every edge sign if necessary so `P=Q>=R`, and set `Delta=P-R`.
If `Delta=0`, do nothing. Otherwise let `k` be the smallest integer with
`binom(k,2)>=Delta`. Assume `k<=floor(n/4)`.
There exists a hollow signing `A'` at the same order satisfying

```math
Q(A')\le Q+k(4KQ/n+1),
\qquad |P(A')-R(A')|\le k((8+8K)Q/n+1).             (1)
```

Only edges of one k-vertex principal block are changed. The construction
uses actual positive and negative maximizing spins; no algorithmic speed
claim is made.

### Choosing the block

The diagonal-majorant lemma gives a principal set `U` of at least `n/2`
vertices with `||A_U||op<=8KQ/n`. Its self-contained proof is Section 1 of
`resumed_bound_audit_minimal_proof_2026_09_06.md`: Boolean polarization,
Gaussian-hyperplane Grothendieck rounding, SDP duality and diagonal deletion.

Let x be a positive maximizing spin. Exact single-vertex optimality gives
`l_i=x_i(Ax)_i>=0`, and `sum_i l_i=2P`. Fewer than `n/4` indices have
`l_i>8P/n`. Thus at least `floor(n/4)` indices of U obey this bound.
Choose any k of them, called S. Put `u=Q(A_S)` and `L=sum_(i in S) l_i`.
Then

```math
u\le4K kQ/n,\qquad L\le8kP/n.                       (2)
```

### The actual signing modification

Let y maximize `-H_A`. Replace the S-block by the switched negative clique
whose off-diagonal entries are `-y_i y_j`; leave all other edges unchanged.
Put `c=binom(k,2)` and `b=floor(k/2)`. Its energy on S lies in `[-c,b]`.
Pointwise comparison and evaluation at the old negative maximizer yield

```math
P'\le P+u+b,\qquad R+c-u\le R'\le R+c+u.            (3)
```

There is also a positive lower bound without keeping the old S-spins.
The old positive spin restricted outside S has energy
`P-L+H_(A_S)(x)>=P-L-u`. Keep these outside spins. On S choose uniformly
among spins with `sum_i y_i z_i=0` if k is even; if k is odd mix the two
uniform slices with weighted sums +1 and -1 equally. Every S-spin has
mean zero; the bridge expectation vanishes and the new clique energy
is always b. Consequently

```math
P'\ge P-L-u+b.                                      (4)
```

Since `0<=c-Delta<k`, equations (3)--(4) give
`Q(A')<=P+u+k` and `abs(P'-R')<=L+2u+k`, proving (1).
Neither an involution nor spectral flatness is used. The only spectral
bound is the selected principal core, not a bound on the full modified matrix.

## 2. Iteration improves orientation error to O(n)

**Theorem.** Fix `C0>0`. For every sufficiently large n and every hollow
signing A with `Q(A)<=C0 n^(3/2)`, there is a hollow signing B of the SAME
order with

```math
Q(B)\le Q(A)+O_{C0}(n^{5/4}),
\qquad |P(B)-R(B)|=O_{C0}(n).                        (5)
```

For integer `Delta>=1`, `k<=sqrt(2Delta)+2<=4sqrt(Delta)`. Under the
bootstrap `Q<=2C0 n^(3/2)`, (1) therefore implies

```math
Delta'\le a\sqrt{n Delta},\qquad
Q'\le Q+b_0\sqrt{n Delta},
a=4[2C0(8+8K)+1],\quad b_0=4[8KC0+1].               (6)
```

Apply (1), reversing overall orientation as needed, while
`Delta>4a^2 n`. Each step then has `Delta'<=Delta/2`. Starting from
`Delta_0<=Q(A)<=C0 n^(3/2)`, the total cap increase before stopping is at most

```math
\frac{b_0}{1-2^{-1/2}}\sqrt{n Delta_0}
\le\frac{b_0\sqrt{C0}}{1-2^{-1/2}} n^{5/4}.          (7)
```

For all sufficiently large n this is at most `C0 n^(3/2)`, closing the
bootstrap inductively. At every step `k=O_(C0)(n^(3/4))<=floor(n/4)`.
The iteration stops in `O_(C0)(log n)` steps, with `Delta<=4a^2 n`.
All comparisons are finite; n then tends to infinity with C0 fixed.

## 3. Consequence for actual optimal values

Take an EXACT minimizing signing at each n, using only the established
uniform bound `M_n<=C0 n^(3/2)`. Equation (5) yields

```math
M_n\le Q(B_n)\le M_n+O(n^{5/4}),
\qquad |P(B_n)-R(B_n)|=O(n).                         (8)
```

Thus the liminf and limsup of the original normalized minimum are unchanged
if at each order one selects these asymptotically balanced near-minimizers.
This is a genuine selectable structural property, not a claim that ALL
near-minimizers, or even one exact minimizer, have this orientation gap.

The previous clique counterexamples show why this quantifier distinction
matters. Balanced orientations still do not imply an isotropic near-ground
law, common optimizer response, low response complexity or an all-order
extension from a fixed seed. The full signing norm may grow faster than
sqrt(n). No original convergence recurrence follows from (8) alone.

## 4. Status and reproducibility

The director read the complete independent audit and reconstructed its
inequalities. The exact checker
`computations/transfer_adversary_same_order_orientation_repair_2026_09_06.py`
has 57 raw actual-signing comparisons, 33 satisfying its finite core/size
conditions, including 12 nonzero-gap repairs. Its results JSON separates
these counts. The analytical proof above is the all-order argument; finite
tests cannot establish (8).

The related theorem
`transfer_adversary_balanced_nearmin_isotropic_defect_2026_09_06.md` proves
that some all-order balanced near-minimizers still have isotropic mean
slack of order n^(4/3). It uses the exact additional estimate
`|I(B)-I(A)|<=k sqrt(k-1)` for a k-block edit. This demonstrates the
remaining distinction between orientation balance and isotropic response.
