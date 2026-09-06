# Independent audit of same-order orientation repair

Date: 2026-09-06. Audited theorem:
`transfer_director_same_order_orientation_repair_2026_09_06.md`.
Verdict: PASS for both the finite lemma and its iterated corollary.
It gives selectable same-order near-minimizers, not a property of every
exact minimizer and not an original convergence recurrence.

## Finite signs, local fields, and the two lower bounds

Orient the parent so `P=Q>=R`, and write `Delta=P-R`.
If Delta is zero no operation is needed. Otherwise let k be minimal
with `c=binom(k,2)>=Delta`, and require the explicit finite condition
`k<=floor(n/4)`.

At an actual positive maximizing spin x, every
`l_i=x_i(Ax)_i` is nonnegative: changing that spin alone cannot increase
the positive maximum. Also `sum l_i=2P`. Therefore fewer than n/4
indices exceed `8P/n`. The already audited simultaneous diagonal-majorant
core has at least n/2 vertices and operator norm at most `8 K_G Q/n`.
Its intersection with the low-field set has at least floor(n/4)
vertices, so the required S exists. For `u=Q(A_S)` and
`L=sum_(i in S)l_i`,

```math
u<=4K_G kQ/n,\qquad L<=8kP/n.
```

Let y be an actual negative maximizing spin. Replacing exactly the
S-principal block by `-y_S y_S^T+I` gives a hollow sign block with
energy interval `[-c,b]`, `b=floor(k/2)`. All other entries, including
the bridge, stay unchanged. Directly,

```math
P-L-u+b<=P'<=P+u+b,
\qquad R+c-u<=R'<=R+c+u.                              (1)
```

The negative lower bound evaluates y itself; its old block energy is
at least -u. For the positive lower bound the outside portion of x has
energy `P-L+H_(A_S)(x)>=P-L-u`. Choose the weighted S-spins uniformly
on the zero-sum slice when k is even, and on the equal mixture of the
sum-plus-one and sum-minus-one slices when k is odd. Both choices give
coordinate means zero and constant new-block energy b. The bridge
expectation vanishes with the outside spins fixed. Some actual spin
therefore attains at least the claimed positive lower bound. Choosing
only ONE of the odd slices would not give this mean-zero argument.

Minimality of k gives `0<=c-Delta<k-1` when Delta>0; the theorem's
weaker `<k` bound is safe. Consequently

```math
Q-u<=Q'<=Q+u+k,
\qquad |P'-R'|<=L+2u+k.                               (2)
```

This also records a useful cap LOWER bound, not needed for the root's
upper theorem but relevant when tracking another response functional.

## Iteration and bootstrap

For integer Delta>=1, `k<=4sqrt(Delta)`. Under
`Q<=2C0 n^(3/2)`, (2) gives exactly the root's constants

```math
Delta'<=a sqrt(n Delta),\qquad Q'<=Q+b_0 sqrt(n Delta),
a=4[2C0(8+8K_G)+1],\quad b_0=4[8K_G C0+1].
```

Above the stopping threshold `Delta>4a^2 n`, the next gap is at most
Delta/2. Summing sqrt(Delta) over the resulting geometric sequence
bounds the total cap increase by
`b_0 sqrt(n Delta_0)/(1-2^(-1/2))=O_(C0)(n^(5/4))`.
The bootstrap is noncircular: at a hypothetical first failed step,
all earlier bounds and the candidate step obey the same total estimate;
for sufficiently large n it is at most `C0 n^(3/2)`, a contradiction.

The block-size hypothesis holds uniformly at every performed step,
since `Delta<=Q<=2C0 n^(3/2)` implies
`k<=4sqrt(2C0)n^(3/4)<=floor(n/4)` for all sufficiently large n.
The number of steps is O(log n), but no growing-order limit is taken
inside an unproved finite lemma. Constants depend only on fixed C0.
Global orientation reversals between steps preserve Q and do not
change any cap estimates.

Thus every bounded-cap input admits a same-order signing with cap
increase O(n^(5/4)) and orientation gap O(n). Starting from an exact
minimizer proves the stated selectable near-minimizer consequence.
It does not preserve exact optimality and does not prove that any exact
minimizer was originally balanced.

## Finite integer replay

Source: `computations/transfer_adversary_same_order_orientation_repair_2026_09_06.py`.
Data: `computations/results/transfer_adversary_same_order_orientation_repair_2026_09_06.json`.
The final replay passed 57 actual-signing tests of every inequality in
(1)--(2). Of these, 33 meet the explicit root-theorem size/core conditions,
including 12 nonzero-gap repairs. The other records test only the raw
block inequalities and are explicitly marked as such.

For these small matrices the full vertex set is a certified spectral
core whenever `n^3(n-1)<=144Q^2`: Frobenius gives its operator bound,
and `K_G>3/2` makes this sufficient. The checker verifies the selected
row mass and block cap using integer inequalities. It exhausts the full
projective spin cube of each parent and modified matrix. The probabilistic
slice argument and asymptotic iteration remain analytical proofs, not
claims certified by a finite experiment. No external download is needed
for replay; the two external-origin upper witnesses are read from the
project's preserved mathematical JSON.
