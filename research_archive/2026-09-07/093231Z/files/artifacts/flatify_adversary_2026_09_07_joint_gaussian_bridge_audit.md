# Independent reconstruction of the joint-Gaussian bridge sector

2026-09-07. Verdict: PASS, with the sector and compatible-order restrictions
in the source retained. Source:
`flatify_independent_2026_09_07_joint_gaussian_bridge_sector.md`.
The source's exact rational checker was read and independently replayed.

## Construction-law coupling

For one physical two-replica word v and a conditioned Gaussian word g,
first leave g unchanged and choose a signed permutation h(g) minimizing
the quotient matching cost between h(g)v and g. Apply a new independent
uniform element q of the signed-permutation group to both arrays. Conditional
on g, q h(g) is uniform, so the physical marginal is exactly the prescribed
common signed-permutation orbit. The Gaussian marginal remains its original
conditioned law: averaging a group-invariant law over independent q leaves
it invariant. This argument does not require h(g) to be independent of g.
Both replica coordinates receive the SAME sign. Repeated physical columns
and stabilizers cause no problem because uniform group action gives the
uniform orbit measure with its correct multiplicities.

The Gaussian event is invariant under this group because it uses the
empirical antipodal-quotient law. Equal-size empirical optimal transport
can be attained by a permutation, with the cheaper common sign for each
matched pair. Triangle inequality and rowwise Minkowski give the stated
global Frobenius error. In particular this is not an illicit coupling of
two independently signed replicas.

## Exponent and determinant

The identity C/sqrt(n)=UKV^T has K orthogonal: each reciprocal two-slot
pair carries H2/sqrt(2), and every slot is used exactly once. Thus the
bridge bilinear energies divided by n^(3/2) are the corresponding diagonal
entries of A^T K B divided by n. Physical Grams equal nG by orthogonality
and the global overlap constraint; no fibrewise exact overlap is needed.

Splitting A^T K B-a^T K b into (A-a)^T K B+a^T K(B-b) yields the stated
error e n, e=eta(2 sqrt(2)+eta). The same bound applies to each Gram.
The penalty-error coefficient is consequently exactly
||T||F+(||L||F+||R||F)/2. Dropping all row conditioning costs at most
2^(2m), before rotating b by K; the rotation need not preserve the
conditioning events. After they have been dropped, orthogonal invariance
does apply, and the resulting n Gaussian four-vectors are independent.

Their original covariance is diag(G,G), whose square-root determinant
is det(G). Integration therefore contributes -log det(G), not half that
quantity, per coordinate. This reconstructs the source's F exactly.
At the rational point the inverse precision is the displayed block
matrix Sigma. Its positive leading minors certify positive definiteness;
its determinant ratio is 25/108. The trace term satisfies
Tr(LG)=t1*(37/50), because the center cross covariance is zero. Hence
t1*(37/50)-F=(1/2)log(108/25), with no missing factor of two.
Changing T to -T is block-sign conjugation of the precision and handles
the other candidate-energy polarity while retaining the absolute center
constraint.

## Uniformity and limitations

There are at most exp(2n h(1/10)) physical cloud pairs. The quotient
profile sector is determined before randomization and is invariant under
the random column operations. Separate couplings for separate candidates
are legitimate for their marginal probability bounds, followed by this
ordinary union bound. No common coupling of all candidates is asserted.
The rate margin is greater than .08, so all o(n) conditioning and W2
errors are absorbed uniformly. Each center feature vector has covariance
I: signs cancel distinct-coordinate correlations and permutation spreads
the fixed row norm. Independence of the two sides gives E h0^2=1/n.

Thus the actual sign bridge exists with center o(n^(3/2)) and the claimed
.74 threshold uniformly on this sector. It does not cover the rest of
either Hamming sphere, nor other energy shells, and does not establish a
child-preserving recurrence. The comparison with the active upper
constant is a genuine restricted-sector saving only above the explicitly
stated seed threshold. No inference for smaller seed constants follows.

## Arbitrary-center extension

The subsequent `flatify_independent_2026_09_07_joint_noisy_center_bridge_sector.md`
also passes independent reconstruction. Here the reference first replica
is an exact signed permutation of the actual center word, and the second
is alpha times it plus independent Gaussian noise. This row law remains
invariant under the common signed-permutation group, so the same matching
proof applies. Matching may change the physical center relative to the
reference center; this is allowed because the Frobenius comparison already
pays the resulting Gram and cross-energy errors.

Direct inversion gives D^-1=[[1/5,1/10],[1/10,1/5]]. Consequently
A-16D^-1=I exactly. The unperturbed conditional replica precision is
(25/9)I, giving determinant ratio (625/81)/(100/3)=25/108. Integration
therefore leaves a pointwise constant in both center arrays, including
after rotating the right array by K. Gaussianity or independence of
center coordinates is genuinely unnecessary.

The uniform transport proof retains each large center coordinate exactly
in the tail coupling, paying only its Gaussian-noise cost. Thus a tail
mass bound R^-2 suffices; an unjustified uniform center second-moment
tail bound is NOT being used. In bounded-center bins, the weighted sum
of empirical Gaussian errors is bounded by log(k+1) sqrt(J/k). The stated
choice R=k^(1/6), h=k^(-1/6) yields O(k^-1/3 log k).

For actual independent bit noise, Hadamard orthogonality gives zero
distinct-coordinate noise covariance. Fixed-dimensional Lindeberg
replacement is uniform because every coefficient is O(k^-1/2), even for
an arbitrarily spiky deterministic center word. Fixed center bins, bounded
noise tests, and uniform fourth moments then justify the claimed uniform
expected W2 convergence. The Hamming-window Gram penalty is exactly
sqrt(2)||L||F omega per n, since both marginal discrepancies contribute
half that trace norm bound each. Typicality is only a 1-o(1) statement;
the source correctly does not discard the exceptional exponential family.
