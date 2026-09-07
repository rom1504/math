# Primary discrepancy literature and the exact polarity interface

Date: 2026-09-07. Outcome: an exact row-free graph-intersection mapping,
but no imported sharp Pareto theorem, polarity balance of exact minimizers,
or reverse-value inequality. The scope distinctions below are essential.

Write `H_A(x)=sum_(i<j) A_ij x_i x_j`, `P=max H_A`, `R=max(-H_A)`.
Let G contain the positive edges of A, `m=binom(n,2)`, `S0=sum_(i<j)A_ij`,
and `p=e(G)/m=1/2+S0/(2m)`.

## 1. Induced-subgraph discrepancy is not the same functional

For `x=2 1_S-1`, direct expansion gives

```math
e_G(S)-p\binom{|S|}{2}
=\frac{H_A(x)+(A\mathbf1)\cdot x+S_0}{8}
 -\frac{S_0}{2m}\binom{|S|}{2}.                       (1)
```

The row-linear term need not be negligible at scale n^(3/2). Switching
A preserves P and R but changes the row term and the graph density.
In particular, switching a positive ground state to the all-one vector
does not eliminate the term; its total local-field sum is `2P`.

[Bollobas--Scott, Discrepancy in graphs and hypergraphs](https://people.maths.ox.ac.uk/~scott/Papers/disc.pdf),
Theorem 1, assumes the ACTUAL density p and `p(1-p)>=1/n`, and proves
`disc_p^+ disc_p^- >= p(1-p)n^3/6400`. The quantities optimize (1) over
induced subsets. The weighted extension, Theorem 8, still uses induced
subset sums, with zero total weight and fixed total absolute weight.
It does not remove the row-linear term. The paper leaves leading
constants unresolved, so neither statement supplies an exact Pareto
frontier near balanced discrepancy.

## 2. A row-free exact mapping via graph intersections

Assume n is even, and let `H=K_(n/2) union K_(n/2)`, with
`h=e(H)=n(n-2)/4` and `q=h/m`. Relabeling H is equivalent to choosing
a balanced spin x. Then

```math
e(G\cap H_x)-\frac{e(G)e(H)}m
=\frac14 H_A(x)+\frac{S_0}{4(n-1)}.                 (2)
```

Indeed the sum of A over edges inside the two shores is
`(S0+H_A(x))/2`, so the overlap count is `h/2+(S0+H_A(x))/4`.
Subtracting `ph` gives (2), using `1-2h/m=1/(n-1)`.
Consequently, with P_bal and R_bal restricted to balanced spins,

```math
\operatorname{disc}^+(G,H)
 =\frac{P_{bal}}4+\frac{S_0}{4(n-1)},\qquad
\operatorname{disc}^-(G,H)
 =\frac{R_{bal}}4-\frac{S_0}{4(n-1)}.                (3)
```

No degree regularity, switching convention, or approximate row cancellation
is assumed. For `Q(A)<=C n^(3/2)`, `p=1/2+O_C(n^(-1/2))`, and the
centering term in (3) is only `O_C(sqrt(n))`.

[Bollobas--Scott, Intersections of graphs](https://people.maths.ox.ac.uk/scott/Papers/graphint.pdf),
Theorem 1, printed p.5, proves for `16/n<=p,q<=1-16/n` that

```math
\operatorname{disc}^+(G,H)\operatorname{disc}^-(G,H)
\ge \frac{p^4(1-p)^4q^4(1-q)^4}{10^{20}}n^3.         (4)
```

Thus this primary theorem genuinely applies to our quadratic problem,
through (3). It only bounds the product at a very small constant; it
does not force the two factors to become equal. The proof uses disjoint
transpositions and controls their quadratic interaction; Theorem 4
bounds the product by a squared average absolute transposition response.
No global minimization over signings occurs. Section 6.2 explicitly
discusses unresolved sharp constants and density dependence.

For example, the direct asymptotic consequence for a bounded-cap sequence is
`liminf P_bal R_bal/n^3 >= 1/(4096*10^20)`. Since P>=P_bal and R>=R_bal
and both restricted extrema eventually have positive order n^(3/2), the
same weak product lower bound follows for P R. This is much weaker than
the original-problem lower chain and is NOT proposed as progress on its
constant. Its value here is the correct, degree-free normalization.

## 3. What the older halves theorem actually says

[Erdos--Goldberg--Pach--Spencer, Cutting a graph into two dissimilar halves](https://www.cs.rpi.edu/~goldberg/publications/discrep.pdf),
Theorem 1, treats `n<e(G)<n(n-1)/4` and finds two n/2-subsets whose
induced edge counts differ by order `sqrt(n e(G))`. These subsets need
not be disjoint: the paper expressly notes that regular graphs have
equal induced counts on every pair of complementary equal shores.
Section 3 defines bipartite discrepancy and conjectures its n^(3/2)
lower order; the graph-intersection result above supplies the later
theorem. Neither the dissimilar-halves title nor its induced-count
conclusion licenses a parent-child cap comparison.

## 4. Exact judicious partitions have the wrong universal quantifier

[Bollobas--Scott, Exact bounds for judicious partitions of graphs](https://people.maths.ox.ac.uk/scott/Papers/bipartitions.pdf),
Theorem 1, guarantees a partition of a graph with e edges for which

```math
e(V_i)\le e/4+\sqrt{e/32+1/256}-1/16,
\qquad
e(V_1,V_2)\ge e/2+\sqrt{e/8+1/64}-1/8.
```

These bounds are exact on suitable odd complete graphs. They control
raw counts, not maximized signed quadratic energies, and do not prescribe
shore sizes. Applied separately to switches of G, the partition can
depend on the switch. A child cap would require a SINGLE partition
for all exponentially many internal switches. Interchanging these
quantifiers is precisely the absent lemma, not a corollary of the
judicious theorem. Its deterministic vertex-move proof also uses
nonnegative edge counts and does not survive replacing each count by
a separately optimized signed cap.

## 5. Recent positive-discrepancy bounds do not cover the transition

[Raty--Sudakov--Tomon, Positive discrepancy, MaxCut, and eigenvalues of graphs](https://people.math.ethz.ch/~sudakovb/positive-discrepancy-and-eigenvalues.pdf)
gives refined lower orders in several average-degree regimes. Its dense
positive-discrepancy regime stops at `d<=(1/2-epsilon)n` for fixed epsilon;
its complementary MaxCut consequence assumes regularity and a density
bounded away from 1/2. Actual bounded-cap signings correspond to
`d=(n-1)/2+O(sqrt(n))`, exactly the excluded transition. The results
therefore cannot silently be used with epsilon tending to zero, or
with a switching-dependent nonregular graph.

## 6. The remaining precise obligation

To force polarity balance from a Pareto method one needs a statement
specific to the actual minimizing signing, for example: every exact
Q-minimizer with P>=R has `P-R=O(n)`, or a cap-preserving change to an
equally minimizing signing with such balance. None of the cited
theorems supplies either. An order-sharp but constant-nonsharp bound
`PR>=c n^3` is not enough: with `max(P,R)<=C n^(3/2)` it only forces
`min(P,R)>= (c/C)n^(3/2)`, not an o(n^(3/2)) gap.

Likewise, an existential child-cap partition inequality would need to
hold for a selected EXACT minimizing parent and prescribed shore sizes;
raw judicious edge counts do not establish it. The separately audited
balanced near-minimizer planting theorem reinforces why exactness
cannot be dropped from such an assertion.

The accompanying exact rational checker verifies (1)--(3) on all
order-4 signings and bounded random samples at orders 6 and 8.
