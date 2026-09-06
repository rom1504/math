# Proper cuts of transported tensor powers

Date: 2026-09-06. Independently derived structural lemma for the next
feedback stage. This is not a complete covariance or energy theorem for
canonical old-tree feedback; its collision and coherent-root-map obligations
remain separate.

## 1. Statement

Let K_1,...,K_n be order-q tensors on a fixed finite tensor-slot pattern,
with q>=2 fixed. The dimension of each slot may depend on n. Assume:

1. The global root map with rows K_i has bounded operator norm.
2. Every proper fixed-root flattening of K_i has operator norm at most
   epsilon_n, uniformly in i and the finitely many cuts.

In particular each row Hilbert norm is uniformly bounded. Let B have
unit row Euclidean norms, entries O(n^-1/2), and bounded operator norm.
For a fixed integer p>=2 define the ordered tensor

```math
L_j=\sum_{i=1}^n B_{ji}K_i^{\otimes p}.
```

Then the global root map of L is bounded and every proper flattening is

```math
O(n^{-1/2}+\epsilon_n+\sqrt n\,\epsilon_n^2). \tag{1}
```

Consequently epsilon_n=O(n^-1/2) implies the same small proper-cut scale
for L. Fixed finite symmetrizations and finite tensor families satisfy
the corresponding statement, with constants depending on that finite data.

## 2. Global root map

Let Gamma_ik=<K_i,K_k>. The row Gram of K_i^{tensor p} is Gamma^{circ p}.
Its operator norm is bounded by the Schur multiplication bound
`||Gamma^{circ p}||op<=max_i Gamma_ii^{p-1} ||Gamma||op`.
Left multiplication by B preserves a bounded root map. This estimate
alone does not prove proper-cut smallness; the latter uses the next
classification.

## 3. Classify the straddling original branches

Fix a left/right cut of the pq formal slots. A branch K_i straddles the
cut if it has slots on both sides. All other branches lie wholly on one
side. Grouping complete branches produces row tensors with bounded root
maps by the same Schur Gram estimate.

If there are no straddling branches, properness puts at least one whole
branch on each side. The flattened L_j is the product of the two bounded
group root maps with middle diagonal diag(B_j). Its norm is O(n^-1/2).

If there are at least two straddling branches, each summand is a tensor
product of at least two proper K_i flattenings and bounded remaining row
tensors. The triangle inequality therefore gives

```math
\|L_j^{\rm cut}\|_{\rm op}
\le C\sum_i|B_{ji}|\epsilon_n^2
\le C\sqrt n\,\epsilon_n^2.
```

Suppose exactly one branch straddles. Its flattening at row i is A_i,
with ||A_i||op<=epsilon_n. Let U_i and V_i be the whole-branch tensors on
the left and right. If both groups are present, the flattening factors as

```math
(I\otimes U^T D_{B_j})
\,\operatorname{blockdiag}(A_i)
\,(I\otimes V),
```

up to permutation of tensor slots. The outer maps have bounded norm and
the middle map has norm epsilon_n; in fact the displayed coefficient
placement gives an extra n^-1/2 when both groups are present.

If only the right group is present, use V to map the right whole tensor
to the direct sum indexed by i, apply blockdiag(A_i), then contract that
root index with row B_j. The last contraction has norm ||B_j||_2=1.
Thus the total norm is O(epsilon_n). The case of only the left group is
the transpose. Since p>=2, at least one whole group is present in this
exactly-one-straddler case. This proves (1).

## 4. Why the degree-one transport is excluded

The condition p>=2 is essential. Let B be an actual normalized symmetric
conference signing, so B²=I. Take three independent Boolean seed colors
and put D_i=S_i^(1)S_i^(2)S_i^(3). Its degree-three tensor is the product
of the three coordinate vectors at vertex i. The degree-three row tensor
K_i of Y=BD has bounded global root map and all proper cuts of norm
max_j|B_ij|=1/sqrt(n-1): the relevant matrices are diagonal after the
color-slot identifications.

But BY=B²D=D exactly. The returned tensor at root i is now a single
coordinate product, with proper-cut norm one and high own-coordinate
influence. Its law is a Boolean sign, not Gaussian. This is an actual
signing example with genuinely independent Boolean colors. The small-cut
class is not closed under an unaccompanied linear transport.

## 5. Deliberate limitation

The lemma concerns the displayed ordered tensors. In a Boolean feedback
application, internal source Walsh projections must remain exact or their
errors must be charged separately. Small proper cuts do not justify
restoring an unrestricted pre-Walsh source tensor; the conference cube
falsifier in the coherent-Walsh audit still applies.

Products of p copies may require deleting marked-label coincidences
between copies. Their local Hilbert error can be small from influences,
and a bounded B transports an averaged L² error safely, but that is a
separate argument from the fixed-root cut estimate proved here. Nor does
(1) establish bounded positive-Walsh root maps for arbitrary polynomials
of the newly coherent return BY. Those are necessary checks before using
this module in a full next-stage energy theorem.

## 6. Exact distinctness preserves proper cuts up to a fixed factor

The feedback agent supplied a useful sharpening, independently checked
here. In the flattening matrix of a fixed tensor, imposing equality
between one left slot and one right slot is rectangular block pinching:
group rows and columns by that slot label and retain only matching blocks.
Its operator norm as a map on matrices is at most one, since each block
is a row/column compression of the original matrix and the resulting
blocks act between orthogonal subspaces. Excluding that equality is the
identity minus pinching, whose norm is at most two.

An equality or inequality involving only slots on one side is simply a
row or column coordinate projection, with norm at most one. Thus imposing
all distinctness constraints among a fixed number of marked slots enlarges
any given proper-cut norm by at most 2 to the number of cross-cut pairs.
This constant is independent of n. The global root-map norm cannot
increase, because the same restriction is an orthogonal projection on
the common marked tensor coordinate space.

The distinctness projection commutes with the outer root transport.
Therefore (1) also holds, up to fixed degree-dependent constants, for
the EXACT diagonal-free projected tensors L_j. This avoids any need to
claim that an unrestricted source tensor has small Hilbert distance from
its Walsh projection. The conference cube example still disproves that
different claim: its full row Hilbert mass is on the deleted diagonal,
even though its proper-cut norm is small.

## 7. Stronger director lemma from all global cuts

The director supplied a stronger complementary argument, independently
verified here. For each of p>=2 rooted tensor families K_j^b, assume every
GLOBAL cut with the root on one side and a NONEMPTY marked subset on the
other has operator norm at most C_b. No small local-cut assumption is
needed. Consider

```math
L_i=\sum_j B_{ij}\bigotimes_{b=1}^p K_j^b.
```

Every proper marked-slot flattening satisfies

```math
\boxed{\|L_i^{\rm cut}\|_{op}
\le\max_j|B_{ij}|\prod_b C_b.} \tag{2}
```

Here is the complete orientation check. Assign each wholly-right branch
to group A and each wholly-left branch to group C. Assign straddling
branches so both groups are nonempty. This is always possible: if there
is only one straddler, at least one other whole branch exists; if there
are none, properness supplies whole branches on both sides; if there
are at least two, assign one to each group.

For every A branch use the global cut mapping its right marked slots
to (root,left slots). Its right marked set is nonempty. Tensor these
maps and compress their output root indices to the common value j.
This is a coisometry, hence gives a bounded map
`A: R_A -> (j,L_A)` of norm at most the corresponding product of C_b.

For every C branch use the transpose global cut mapping (root,right
slots) to its left marked slots, which are nonempty. Embed one shared
input root j into the repeated root indices of their tensor product.
This is an isometry, giving `C:(j,R_C)->L_C` with the same product bound.
The desired flattened tensor is exactly

```math
(I_{L_A}\otimes C)
\,[D_{B_i}\otimes I]\,
(A\otimes I_{R_C}),
```

up to slot permutations. The middle diagonal has norm max_j|B_ij|,
proving (2). Root summation is the displayed composition, not an
uncontrolled partial trace. Empty left or right pieces of individual
branches are harmless because the chosen opposite marked part is
always nonempty.

The distinctness projection in Section 6 preserves (2) up to a fixed
degree-dependent factor. Thus all-global-cut bounds that are fixed
powers of log(n) give proper local cuts n^-1/2 times a fixed logarithmic
power after any transport of at least two factors. A single factor is
still excluded, consistently with the exact conference-return example.

## 8. Consequence for the centered local noise coefficient

In the first marked-history setting, suppose A(W) is a fixed EVEN
coherent polynomial and subtract its row mean. Its positive Boolean
Walsh degrees are a>=2. Let Z_q be one exact homogeneous residual-noise
channel, q>=3. The established local noise/coherent contraction argument
charges the difference between (A-EA)Z_q and its squarefree tensor-product
main in averaged L²: expand A into the primitive old coherent fields
of degree at most three, merge a noise-touching coincidence first, and
use its small proper cut. The equal-degree-three complete pair is the
already small h3(G)--marked-input covariance. All internal coherent
collisions remain exact.

The all-global-cut bounds from higher Boolean derivative matrices apply
to both exact factors. Applying (2) to B[(A_a)Z_q] gives a homogeneous
output of degree a+q>=5 with proper cuts O(n^-1/2 polylog(n)). At a
fixed root, pair this with any fixed polynomial c0(W), expanded into
the degree-at-most-three coherent primitives. The first-merge bound
gains the small cut; all remaining merges are Hilbert-contractive.
Consequently the normalized local pairing tends to zero. Symmetry of
B rewrites it as

```math
\frac1n\mathbb E\sum_{ij}B_{ij}
  (A_i-\mathbb EA_i)Z_i c_j^0=o(1). \tag{3}
```

The source-error transfer uses only the fixed operator cap and averaged
L², followed by bounded fixed-degree moments of c0. This is precisely
the centered-coefficient obstruction isolated in the marked-star working
note. It does not, by itself, certify every other mixed-noise term in
the complete marked feedback energy theorem.
