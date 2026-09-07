# Complete finite Hadamard-bridge obstruction for optimal order-eight children

2026-09-07. Exact finite theorem, checked by two different algorithms.

Let A,D range over ALL hollow symmetric full signings of order eight with
Q(A)=Q(D)=M_8=10, and let B range over ALL order-eight real Hadamard matrices.
Then

    min Q([[A,B],[B^T,D]]) = 30.

The cap-30 witness is the actual matrix saved and independently checked in
`flatify_adversary_2026_09_07_actual_children_hadamard_bridge.md` and its
linked computation files. The new content is the complete lower bound,
including a FREE choice of the second optimal child, not merely sampled
bridge phases, sampled permutations, or frozen-child local trades.

The normalized parent cap is 30/64=.46875, so this remains a positive
sub-half construction. But it misses the favorable equal-child weighted
target 20 sqrt(15/7)=29.277... . Since every order-sixteen full-sign energy
is an even integer, attaining that target would require cap at most 28.
That is impossible in this complete finite Hadamard-bridge class.

## 1. Exact shell intervals

Use projective spins x,y with their first coordinate fixed to +1; independent
global reversals leave child energies unchanged and reverse bridge energy.
Consequently

    Q(parent)=max_{x,y} (|H_A(x)+H_D(y)|+|x^T B y|).

For fixed A,B, a cap at most28 is equivalent to the 128 inequalities

    L_y <= H_D(y) <= U_y,
    L_y=max_x(-28+|x^T B y|-H_A(x)),
    U_y=min_x( 28-|x^T B y|-H_A(x)).

Because D must itself be optimal, intersect these intervals with [-10,10].
This last requirement is an explicit hypothesis of the theorem; it does not
follow just from parent cap28. The interval reduction is exact, not an
entropy or ground-shell approximation.

## 2. Complete child enumeration, independent of catalogue completeness

A vertex switching uniquely makes all seven edges in the first row positive.
Enumerate all 2^21 choices for the remaining edges. Encode a signing by the
28-bit mask a of its negative edges, and encode a spin by its 28-bit cut mask
c_x. Then, exactly,

    H_A(x)=28-2 popcount(a XOR c_x).

Checking every projective spin leaves exactly4200 normalized signings with
cap at most10, and none with cap at most8. Thus M_8=10 is independently
reproved, not simply imported from the given small-order table.

The verifier next forms every permutation fixing vertex0 of the four source
representatives (atlas indices580,722,870,1005). Their resulting set of
normalized signings equals EXACTLY the directly enumerated4200 matrices.
Thus these representatives are proved complete for this computation.
Expanding the128 distinct switching gauges gives exactly537600 distinct
full optimal second children, with no use of external completeness claims.

Switching and permutation of the first child preserve the Hadamard property
of its bridge. It is therefore sufficient to test the four representatives
as the first child.

## 3. Complete Hadamard-frame enumeration

Column signs and permutations of B are absorbed by the FREE second child D
and preserve its cap. Thus each column of B may be projectively normalized
to have first coordinate+1, and columns may be sorted.

Every such B is exactly an eight-clique in the orthogonality graph on the128
projective Boolean vectors. The solver-free checker enumerates these cliques
directly using integer bitsets and finds exactly480 distinct frames.

The first, independent algorithm used a different complete generation: after
switching by one chosen column, a frame contains the all-ones vector; the
remaining seven columns form an orthogonal clique among35 balanced projective
vectors. There are30 normalized frames. Switching them by all128 projective
vectors and deduplicating gives480 frames, consistently with30*128=480*8.
The adversarial researcher separately checked this quotient argument.

## 4. Two complete checks

First algorithm: for each first-child representative and frame, solve the
28 Boolean edge variables of D subject to all shell intervals. There are
1728 distinct interval systems after exact deduplication. CP-SAT returned
INFEASIBLE for all of them, with no UNKNOWN or timeout cases; runtime228s.

Second algorithm: NO SAT, MILP, LP, floating-point, or numerical tolerance.
For every one of the4*480 representative/frame pairs, begin with the directly
enumerated537600 optimal second children. Intersect the shell intervals by
exact cut-mask popcounts, stopping when no candidates remain. Every pair
has zero surviving children. Runtime about25s. This independently verifies
the solver's global infeasibility conclusions rather than merely replaying
its statuses.

The matching cap30 witness plus parity proves the stated minimum exactly.

## 5. An algebraic shortcut also fails in this finite class

The proposed correlated completion D=-B^T A B/8 would make the off-diagonal
blocks of the parent's square vanish. It cannot produce a hollow full
signing D for any of the four optimal children. Each child has exactly12
projective zero-energy spins, but none of the480 Hadamard frames consists
entirely of them. The exact checker records the strongest zero-column count
and the least squared energy sum over frames.

Even asking only that every OFF-diagonal entry of B^T A B/8 be a full sign
does not avoid this obstruction. Orthogonal conjugation preserves squared
Frobenius norm56; its56 off-diagonal entries already consume this entire
norm if they have magnitude one, forcing every diagonal entry to vanish.

This finite algebraic obstruction is not a statement about approximate
conjugation, other bridge matrices, or asymptotic minimizing children.

## Reproduction and scope

- `computations/flatify_independent_2026_09_07_hadamard_free_child.py` and
  its same-stem result JSON preserve the first CP-SAT check.
- `computations/flatify_independent_2026_09_07_hadamard_free_child_exact.py`
  and its same-stem result JSON give the solver-free exhaustive replay.

No claim is made about arbitrary non-Hadamard bridges at order16, the actual
value of M_16, or a limiting obstruction for larger children. In particular,
this does not disprove an asymptotic recurrence with an o(n^(3/2)) error.
