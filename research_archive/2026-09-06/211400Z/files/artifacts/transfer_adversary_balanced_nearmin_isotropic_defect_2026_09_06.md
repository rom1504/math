# Balanced near-minimizers can still have superlinear isotropic response defect

Date: 2026-09-06. Actual same-order signing theorem, deduced from the
audited switched-clique repair and planted-clique construction. This is
a counterexample for ALL balanced near-minimizers, not for exact
minimizers or the existence of a selected good sequence.

Write `I(A)=max_(E xx^T=I) E|H_A(x)|` and
`D(A)=Q(A)-I(A)`, the minimum isotropic mean ground slack.

## 1. Exact block-edit continuity of I

Suppose A and B are hollow signings of the same order, differing only
inside a k-vertex principal block. For `E=B-A` and ANY isotropic law,

```math
E_\mu|H_E(x)|\le\tfrac12\operatorname{tr}|E|
\le\tfrac12\sqrt{k}\|E\|_F
\le k\sqrt{k-1}.
```

The first inequality uses spectral absolute value and the raw identity
second moment. The last uses zero diagonal and entries bounded by two
in that block. Triangle inequalities in both directions give

```math
|I(B)-I(A)|\le k\sqrt{k-1}.                            (1)
```

This is valid for every actual signing edit and requires neither
independence nor a bounded operator norm for the whole parent.

## 2. The same-order orientation repair approximately preserves I

Use the algorithm in
`transfer_director_same_order_orientation_repair_2026_09_06.md`, whose
independent audit is
`transfer_adversary_same_order_orientation_repair_audit_2026_09_06.md`.
For input `Q(A)<=C0 n^(3/2)`, let `Delta_0=|P(A)-R(A)|`.
Its performed steps have gaps `Delta_j<=2^(-j)Delta_0` and block sizes
`k_j<=4sqrt(Delta_j)`. Global orientation reversals preserve I exactly.
Consequently (1) yields

```math
|I(B)-I(A)|
\le\frac{8}{1-2^{-3/4}}\Delta_0^{3/4}
=O_{C0}(n^{9/8}).                                     (2)
```

The finite repair lemma also gives `Q_next>=Q_old-u`, where
`u<=4K_G k Q_old/n`. With the proved bootstrap
`Q_old<=2C0 n^(3/2)`, summing the same geometric sequence gives

```math
Q(B)\ge Q(A)
-\frac{32K_G C0}{1-2^{-1/2}}\sqrt{n\Delta_0}.
```

Together with the root's upper bound, the selected repaired B obeys

```math
|Q(B)-Q(A)|=O_{C0}(n^{5/4}),\qquad
|P(B)-R(B)|=O_{C0}(n).                                 (3)
```

Thus this specific orientation repair preserves the leading isotropic
response as well as the leading cap. It cannot automatically repair a
larger response defect. Equations (2)--(3) are properties of the
constructed output, not every arbitrary orientation-balancing edit.

## 3. An actual all-order near-minimizer family

The theorem in
`transfer_adversary_nearmin_clique_orientation_gap_2026_09_06.md`
starts at an exact minimizer of each order and plants a positive clique
of size `r=floor(n^(2/3))` inside a bounded-operator principal core.
It supplies actual signings A_n with

```math
Q(A_n)=M_n+O(n^{4/3}),\qquad
D(A_n)\ge(1/4-o(1))n^{4/3}.                            (4)
```

For clarity, the second claim follows from that theorem's orientation
gap at least `(1/2-o(1))r^2` and its exact chord estimate
`D(A)>=P(P-R)/(P+R)>=(P-R)/2` in the dominant positive orientation.
No computation of a large-order minimizer is required for this
existence construction.

Apply the same-order repair to A_n. Its normalized cap is uniformly
bounded, so (2)--(3) apply with a fixed C0. We only use
`Delta_0<=Q(A_n)=O(n^(3/2))`; we do NOT assume that the original exact
minimizer had balanced extrema. Combining (2)--(4) proves

```math
\boxed{\begin{aligned}
M_n\le Q(B_n)&\le M_n+O(n^{4/3}),\\
|P(B_n)-R(B_n)|&=O(n),\\
Q(B_n)-I(B_n)&\ge(1/4-o(1))n^{4/3}.
\end{aligned}}                                       (5)
```

The potential cap decrease is O(n^(5/4)) and the potential increase
in I is O(n^(9/8)); both are smaller than the defect in (4).
This is why a lower bound on the repaired cap is needed in addition
to its upper near-minimizer estimate.

## 4. Scope and the exact-minimizer question

Equation (5) disproves a claim that EVERY asymptotically minimizing
sequence with O(n)-balanced oriented extrema has
`I(A)>=Q(A)-O(n)`. Such a sequence can instead require mean isotropic
slack of order n^(4/3), despite having both oriented extrema within
O(n) of its cap. In particular, no isotropic law is supported on an
absolute near-ground shell of slack o(n^(4/3)) for these examples.

The normalized defect in (5) still tends to zero; no macroscopic
order-n^(3/2) defect for actual minimizers is asserted. The matrices
B_n are not claimed to be EXACT minimizers. Therefore (5) leaves open
both root's precise targets: whether EVERY exact minimizer, or at least
SOME exact minimizer at each order, has `I>=Q-O(n)`.
No original convergence or nonconvergence conclusion is drawn.

The finite block identities underlying this corollary were replayed by
the same-order repair checker. Equations (1)--(5) above are analytical
all-order estimates and do not depend on an extrapolation of that
finite experiment.
