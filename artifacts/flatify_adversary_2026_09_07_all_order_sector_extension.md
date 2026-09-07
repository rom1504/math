# All-order extension of the actual noisy-center bridge SECTOR theorem

2026-09-07. This removes an order restriction from a sector operation,
not from an original-value recurrence. It uses the strict c=47/100
certificate in `flatify_adversary_2026_09_07_conditional_shell_rate.md`.

## Dense available orders

Hadamard orders 2 and 12, with Kronecker products, supply orders
k=2^a 12^b. Since log_2(12) is irrational, finitely many of its nonnegative
integer multiples form an arbitrarily fine net modulo one. For any fixed
eta>0 and all sufficiently large targets t, choose b from such a finite
set and then a>=1 so that t<=k<=(1+eta)t. Consequently there is a choice
of even Hadamard order k with k>=sqrt(2n) and k/sqrt(2n)->1 for EVERY
integer n tending to infinity. Put n_+=k^2/2 and r=n_+-n=o(n).
No quantitative density rate is needed. Order 12 can, for example, be
constructed from the quadratic-character matrix Q on F_11 as
[[1,1^T],[1,Q-I]]; Q^T=-Q and QQ^T=11I-J verify the Hadamard identity.

## Padding and restriction are uniform operations

Pad arbitrary physical centers x0,y0 of length n to length n_+ by fixed
signs. Pad every candidate x,y by those SAME signs. Construct the full
rank-two sign bridge C_+ of order n_+ using the sector theorem, then let
C be its n-by-n principal rectangular restriction (the two sides need
not use identical coordinate labels). Its entries remain full signs.

The exact full construction has ||C_+||op=sqrt(n_+). Splitting a bilinear
form into original and padded coordinate parts yields, uniformly in all
physical spins,

    |x_+^T C_+ y_+-x^T C y| <= 2 n_+ sqrt(r)
                               =o(n^(3/2)).

The same estimate applies to the centers. Thus a full center value
o(n_+^(3/2)) remains o(n^(3/2)) after restriction, with no probabilistic
assumption on the padded spins. The explicit n^(5/4) center rate need
not survive; the vanishing normalized rate does.

For an exact common Hamming distance j on the original two sides, the
padded overlap is alpha_+=1-2j/n_+, whereas alpha=1-2j/n. Their difference
is at most 2r/n_+=o(1), uniformly j. When alpha>=0 the padded overlap is
also nonnegative. Negative overlaps are first reduced by global spin
flips on the original side and then padded. This preserves child
quadratic energies and absolute bridge energies. The original child
energy window contributes 2c alpha^2 n^(3/2), which differs from its
padded counterpart 2c alpha_+^2 n_+^(3/2) by o(n^(3/2)).

## Precisely defined inherited profile sector

For each original pair use its padded two-state Hadamard feature words
and the reference law (Z,alpha_+ Z+sqrt(1-alpha_+^2)G), with Z the EXACT
padded center coefficient law in that fibre. This defines a genuine
sector of original physical pairs. It is not silently identified with
an unpadded square-Hadamard profile. Apply the all-noise theorem on the
full order n_+ to this sector. It already controls every exact common
noise level, so no new entropy union is required by padding.

For all c>=47/100, the strictly positive uniform certified cap gap absorbs
the above restriction, overlap, and normalization errors. Thus there
exists, for every sufficiently large integer n, a full n-by-n sign bridge
whose center is o(n^(3/2)) and whose parent objective in the inherited
typical-energy/profile sectors is strictly below 2sqrt(2)c n^(3/2) by a
fixed leading gap. Endpoint overlap intervals continue to use the full
construction's deterministic operator estimate before restriction.

## Typicality survives, uniformly in arbitrary centers

Consider independent physical bit noise of a common parameter delta in
the middle compact overlap interval. Generate full-order iid noise with
the deterministic parameter delta_+=delta*n/n_+. Couple the original
coordinates by their common uniforms. Changing
the noise parameter costs O(r) expected mismatches; fixing all padded
coordinates costs at most r further mismatches. The squared Euclidean
distance between the two candidate spin words is therefore O(r) in
expectation, uniformly in the centers. Block-normalized Hadamard feature
maps are orthogonal, so the averaged squared quotient-W2 distance between
their feature profiles is O(r/n_+)=o(1).

The full iid-noise profiles converge uniformly to the retained-center
reference law by the previously proved bivariate replacement argument.
Combining that result with the padding comparison gives typicality of
the inherited sector for original iid noise, uniformly in all centers.
The actual overlap concentrates about the deterministic reference overlap.
For Hamming windows use their actual overlap parameter; its o(1) changes
in the middle compact interval alter the reference law by o(1) in W2.
This statement does not claim iid typicality transfers to an exact
Hamming sphere by dividing by its probability; no such unjustified
conditioning step is used.

## Scope

Both children may be actual minimizers at the original arbitrary order n;
none of their edges is changed by this bridge restriction. Both energy
polarities are covered by the absolute center constraint and the two
source signs in the original theorem. Nevertheless only the stated
same-noise-level profile and energy sectors are controlled. Exceptional
pairs, unequal noise levels, and a paid covering by centers remain open.
This all-order result therefore supplies no full-parent bound and no
convergence conclusion.
