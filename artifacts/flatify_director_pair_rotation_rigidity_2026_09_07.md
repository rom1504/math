# Pair-local orthogonal flatification: exact rigidity

Date: 2026-09-07. Status: proved elementary classification; not a
flatification theorem for unrestricted operations or a convergence result.

## Statement

Let A,D be hollow symmetric sign matrices of common order k>=3. Pair
vertex i of A with vertex i of D. In this ordering the off-pair tiles of
B=sqrt(2) diag(A,D) are sqrt(2) diag(a_ij,d_ij). Allow an independent real
rotation R(theta_i) in each pair and conjugate B by their direct sum.

Every entry of every off-pair 2 by 2 tile has absolute value one if and
only if there are signs t_i such that

    d_ij = -t_i t_j a_ij   (i != j).

The within-pair tiles remain zero; filling those k missing matching edges
costs at most k in Q. The scale sqrt(2), rather than the exact row-normalizing
sqrt((2k-1)/(k-1)), changes Q by O(sqrt(k)) for bounded-normalized-cap
children. Thus that minor normalization does not alter the leading issue.

This statement covers SO(2) pair rotations only as stated; no assertion
about unrestricted orthogonal changes of basis is needed.

## Proof

Put sigma_ij=a_ij d_ij. For sigma=+1, the tile is, up to its sign,
sqrt(2) times a rotation of angle theta_i-theta_j. It is flat precisely
when theta_j-theta_i=pi/4 modulo pi/2. For sigma=-1, conjugation of a
reflection gives the condition theta_j+theta_i=pi/4 modulo pi/2.
With phi_i=4 theta_i/pi, both conditions are

    phi_j = sigma_ij phi_i + 1   modulo 2.

Compose this identity along a triangle. If its product of three sigma's
is +1, the result is phi_i=phi_i+1 modulo 2, impossible: the sum of three
signed ones is odd. Hence every triangle has product -1.

Conversely, this triangle condition implies sigma_ij=-t_i t_j: fix t_1=1
and t_i=-sigma_1i, and use each triangle containing vertex 1. Choose
phi_i=t_i/2 modulo 2. If t_i=t_j then sigma=-1 and the equation holds;
if t_i=-t_j then sigma=+1 and the equation also holds. This supplies the
rotations and proves the assertion.

## Consequence and precise limitation

Selectable equal children can always meet this condition by taking D=-A.
The resulting operation is the H2 lift, up to switches, vertex order,
and its missing matching. Its Boolean cap is not controlled by its
orthogonal spectrum: it sees all clique-flipped child objectives.
That exact formula and its scoped failures are already in
`decisive_independent_h2_exact_profiles_and_algebra_2026_09_07.md` and
the subsequent H2 counterexample files. This note explains why varying
the pair rotation angles does not produce an independent escape.

In particular, it neither disproves H2 stability for suitably selected
exact minimizers nor excludes nonlocal sign replacement. It removes only
the proposed additional freedom of exact pair-local rotation angles.

## A different, unbounded-amplitude caution

For N divisible by 4, take N/4 disjoint hollow bipartite H2 blocks, each
scaled by sqrt((N-1)/2). Every row has squared norm N-1 and the total cap
is (N/4)sqrt(2(N-1)), asymptotic to N^(3/2)/(2 sqrt(2)). Thus row-square
normalization alone cannot give lossless full-sign realization: its cap
is below the established full-sign lower bound. Here amplitudes grow as
sqrt(N); this example says nothing against the bounded-amplitude,
two-optimal-child target. The block cap follows from
max_{x,y in {+1,-1}^2}|x^T H2 y|=2, checked directly from H2 y.
