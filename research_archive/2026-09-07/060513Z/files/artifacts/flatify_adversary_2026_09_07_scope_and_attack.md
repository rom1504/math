# Selected-optimal-child flatification: adversarial scope and first-principles attack

Status: elementary identities proved below; proposed mechanisms remain open.
Independent work started 2026-09-07 during the six-hour campaign.

## 1. The existential target is exactly a scalar recurrence

Write `P(A)=max H_A`, `R(A)=max(-H_A)`, `W(A)=(P+R)/2`, and
`I(A)=(P-R)/2`. For positive weights a,b and disjoint blocks,

```math
Q(\operatorname{diag}(aA,bC))
=aW(A)+bW(C)+|aI(A)+bI(C)|.
```

This follows by choosing each block spin independently: the positive endpoint
is aP(A)+bP(C), the negative endpoint is aR(A)+bR(C). No intervening energy
values or polarity-balance assumption is needed.

For actual cap minimizers of orders m,l, choose their global edge polarities
so that P(A)=M_m and P(C)=M_l. Global edge negation preserves optimality.
Then, with a=sqrt((N-1)/(m-1)), b=sqrt((N-1)/(l-1)),

```math
Q(B)=aM_m+bM_l.
```

Consequently the existential selected-child replacement with any prescribed
nonnegative error E(N) is **equivalent**, not merely sufficient, to

```math
M_N\le aM_m+bM_l+E(N).
```

Forward: Q(B)<=aM_m+bM_l for every choice. Reverse: orient arbitrary exact
minimizers as above and choose an exact optimal parent. The formulation does
not demand that the parent be algorithmically obtained from the children.
In particular, classifying optimal-child isomorphism types or proving their
rigidity has no effect on the bare existential target. Such information may
still be essential to a concrete construction proving the scalar inequality.

For identical optimal children with opposite polarities, instead,

```math
Q(\operatorname{diag}(aA,-aA))=2aW(A).
```

Demanding the theorem for every optimal-child pair, or specifically for these
opposite polarities, is a genuinely stronger width-sensitive conjecture. A
counterexample to it does not disprove selectable-child flatification.

## 2. What convergence alone would give

If M_n/n^(3/2) converges to c, then uniformly for comparable m,l,

```math
M_N-aM_m-bM_l=o(N^{3/2}).
```

Indeed divide by N^(3/2); the weighted sum is
`sqrt(1-1/N)[(m/N)sqrt(m/(m-1)) m_m+(l/N)sqrt(l/(l-1)) m_l]`,
which tends uniformly to c. This supplies no power saving. Conversely the
recorded summable-error balanced-tree proof is substantive: a mere vanishing
one-step defect is not automatically enough. Thus the selected-child
formulation with nonquantitative error should not be mistaken for an easier
structural rigidity statement.

## 3. Actual full-sign finite census

`tmp/flatify_adversary_2026_09_07_exact_minimizers.py` independently enumerates
all full signings modulo first-row switching and all spins modulo reversal.
All energies are integers. Its JSON contains full endpoint histograms and
representative signings, not only claimed optima.

| n | M_n | minimum W | endpoint pairs among cap minimizers |
|---|---:|---:|---|
| 3 | 3 | 2 | (3,1), (1,3) |
| 4 | 4 | 4 | (4,4) |
| 5 | 4 | 4 | (4,4) |
| 6 | 5 | 5 | (5,5) |
| 7 | 9 | 8 | (9,7), (7,9) |

These are exact exhaustive finite checks. They show that midpoint balance is
not automatic even for *every* exact minimizer at an order. They neither
establish an asymptotic midpoint obstruction nor refute the selected-child
target. The order-four width minimizers additionally include endpoint pairs
(6,2) and (2,6), illustrating why width minimization cannot silently replace
cap minimization in a structural premise.

## 4. Independent first-principles attack

The original optimization is over the vertices of an edge cube, with objective
the support function of the signed spin-correlation polytope. An exact
minimizer A satisfies the following exact global exchange covering property:
for every edge subset F there are a spin x and polarity s such that

```math
sH_A(x)-2s\sum_{e\in F}a_e x_i x_j\ge M_n.
```

This is just optimality of the sign matrix obtained by flipping F, and it is
valid for arbitrarily large simultaneous exchanges. It is stronger than
single-edge stationarity, which can be vacuous on an integer cap plateau.

A concrete prospective mechanism is to apply this covering property to a
distribution of **global** parent edge exchanges whose covariance preserves
the row-square budget. One would need to show that failure of favorable
rounding produces a subset F improving one actual child; a generic bounded-cap
counterexample does not satisfy the required global-optimality premise.
The missing inequality is a transfer from the parent's failed rounded-spin
witnesses to a common child exchange F. Witnesses depend on the entire random
rounding, so exchanging max and expectation is not valid.

This proposal is not yet a reduction: without that transfer it simply
reexpresses optimality. A bounded diagnostic separately optimizes all bridge
edges for several small exact optimal children while leaving old edges fixed;
its purpose is to measure the additional price of freezing old edges, not to
test the fully global target. Solver lower conclusions are numerical unless
replaced by exhaustive or exact certificates.

### Archive collision and diagnostic outcome

The simultaneous exchange identity was subsequently found already recorded in
`decisive_independent_reverse_value_restart_2026_09_07.md`. It is not a new
structural result. Independent edge-noise entropy arguments also collide with
`decisive_independent_gaussian_stability_2026_09_06.md`, which gives the stronger
correlated-Gaussian-width necessary condition. The small endpoint census
independently replays a subset of the earlier minimizing-orbit screens; it is
verification, not new primary progress.

The frozen-child bridge diagnostic is saved as
`tmp/flatify_adversary_2026_09_07_frozen_children.{py,json}`. For both relative
polarities, it attains caps 10,13,18,25 at parent orders 8,10,12,14. HiGHS
reports optimality at orders 8 and 10; the runs at 12 and 14 stop after 45
seconds each, with lower bounds far below the returned caps. Every returned
bridge has its cap exhaustively rechecked with integer spin energies. These
runs do not establish an obstruction. In particular, the archive's balanced
partition screen already supplies order-14 optimal parents with two exact
order-seven children, though not necessarily the particular switching class
and relative orientation frozen by this diagnostic.
