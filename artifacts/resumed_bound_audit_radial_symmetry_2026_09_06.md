# Independent audit of signed edge-transitive radiality

Date: 2026-09-06. Outcome: PASS for
`resumed_director_radial_symmetry_principle_2026_09_06.md`.
Assume n>=2, so the hollow signing is nonzero.

For a group of signed permutations satisfying
P_g^T A P_g=epsilon_g A, the sign epsilon is uniquely determined and
multiplicative. The reversed conjugation P_g A P_g^T=epsilon_g A also
holds. Starting from any oriented ground (x,sigma), the samples
(P_g x,epsilon_g sigma) retain oriented energy Q(A).

For M=E[epsilon_g sigma P_g xx^T P_g^T], a change of group variable
g -> hg gives P_h M P_h^T=epsilon_h M. Comparing with the identical
transformation law of A shows that M_ij/A_ij is constant on underlying
unordered-pair orbits: the two switching signs cancel. Pair transitivity
therefore gives M_ij=rho A_ij for i!=j. Summing the expected oriented
energy fixes rho=2Q(A)/(n(n-1)), with no additional factor of two.

The diagonal is M_ii=sigma E epsilon, not necessarily zero. It vanishes
for a nontrivial character. Hollowness of the tested B makes it irrelevant
without that extra symmetry; a nonhollow B would require a separate
diagonal term.

The averaged oriented B energy is rho sum_(i<j) A_ij B_ij. Taking its
absolute value proves the stated lower bound. For d signing edits,
sum A_ij B_ij=n(n-1)/2-2d, giving exactly

    Q(B)>=Q(A) |1-4d/(n(n-1))|.

All order-n signing caps have the parity of n(n-1)/2. Thus a strict cap
decrease is at least two, whereas d<n(n-1)/(2Q(A)) permits a decrease
strictly smaller than two. The local stability radius follows.

If Q(A_n)/n^(3/2)->c>0 and Q(B_n)/n^(3/2)<=c-delta+o(1), solving the
absolute-value inequality bounds d/n^2 below by delta/(4c)-o(1) and
above by 1/2-delta/(4c)+o(1). The upper bound is precisely the same
lower bound for distance from -A_n. Thus the coefficient delta/(4c)
and its two-orientation scope are correct.

For symmetric Paley cores, the affine action over the field is pair
transitive and has epsilon equal to the quadratic character of the
scale. Arbitrary minimizers need not have this symmetry. Averaging over
all switchings or relabelings that move A does not create a radial
ground-state law for the fixed A.
