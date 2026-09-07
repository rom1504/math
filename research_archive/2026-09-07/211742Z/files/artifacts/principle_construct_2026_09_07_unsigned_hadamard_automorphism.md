# Exact unsigned transitivity for the normalized Paley12 and tensor24

2026-09-07. **Exact integer certificate PASS.** This checks the finite
symmetry property proposed for the stratified-node balanced-face argument;
it does not replace the entropy/Schur proof using that property. In fact
the output ROW PERMUTATIONS mix distinct child laws under averaging, so
unsigned transitivity alone is insufficient to complete that argument.
The finite automorphism certificate remains valid, but L24 balanced-face
inheritance is not established here.

Use the normalized Paley12 matrix in the director's selector checker.
Form its incidence graph with12 column points and22 blocks, the positive
and negative supports of each nonconstant row. A graph isomorphism
individualizing point0 to point1 supplies the unsigned column permutation

    [1,0,11,2,7,3,8,6,4,10,9,5].

Direct integer multiplication verifies that H[:,p] H^T is12 times a
signed permutation matrix and that its (0,0) entry is+12. Thus the DC
row is fixed positively, and no input sign switches are used.
The finite-field translation [0,2,3,4,5,6,7,8,9,10,11,1] has the same
property. The generated orbit of point0 is all12 points.

For H24=H2 tensor H12, tensor these two generators with identity and add
the H2 column swap. Every resulting H24[:,p] H24^T is24 times a signed
permutation fixing positive DC; their point orbit has size24. These are
the required transitive unsigned input symmetries.

Reproducible code:
`../computations/principle_construct_2026_09_07_unsigned_hadamard_automorphism.py`.
The corresponding results JSON contains the full matrix, all generators,
the signed output permutation, and both exact point orbits. Graph search
discovers the witness; the finite integer matrix equalities independently
certify it and can be replayed without relying on automorphism software.

## The stronger diagonal-output condition forces a Sylvester order

Suppose a transitive unsigned input permutation group satisfies H P=D_P H
with every D_P diagonal signs. Conjugation by H identifies this group
faithfully with a subgroup of the diagonal sign group. It is elementary
abelian2, so each point orbit has power-of-two size. In particular L must
be a power of2. More precisely a faithful transitive abelian action is
regular: an element fixing one point fixes every point by commutation.
Thus the group has order L, and each output row is a character on this
group. Orthogonality supplies all distinct characters, giving a Sylvester
character table up to permutations and row/column signs.

Consequently no alternative search at order24 can furnish the stronger
diagonal-output transitivity premise. Using Paley24 requires an actual
additional envelope/entropy argument; the checked unsigned automorphisms
with output row permutations cannot be silently upgraded.
