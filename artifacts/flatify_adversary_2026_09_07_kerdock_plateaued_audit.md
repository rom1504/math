# Independent quadratic-phase MUB pentagon obstruction audit

Status: **PASS**. Independently reconstructed the growing-dimension argument
in `flatify_construct_2026_09_07_kerdock_plateaued_witness.md`, before relying
on its finite verifier. The unused-reference scope is valid as stated.

The result is exact: with the specified actual optimal pentagon seed, every
five-quadratic-phase MUB cross construction, and the stated subsets of the
explicit Kerdock family including Id, has cap at least `5t^(3/2)` against
every arbitrary hollow diagonal-fibre completion. It is not an obstruction
for arbitrary selected large minimizing children.

## 1. The explicit field family and its nonsingular pair differences

For odd m and K=F_(2^m), write W=K direct-sum F2. The displayed B_a is
alternating because `Tr((ax)²)=Tr(ax)` and binary squaring fixes the trace.
For a!=b, substitute c=a+b, alpha=a/c, u=cx, v=cy. Expanding the two products
of traces gives exactly equation (2) of the source; no cross term is omitted.

If (u,s) is in its radical, the coefficient of r gives Tr(u)=0. The remaining
coefficient of v is the field trace pairing with
`u+[Tr(alpha u)+s]1`. Nondegeneracy of the trace pairing sets this to zero.
Taking its trace, using Tr(1)=1 for odd m, gives Tr(alpha u)+s=0, hence u=0
and then s=0. Thus every pair difference is nonsingular. The ambient binary
dimension m+1 is even, as required for nonsingular alternating forms.

For a quadratic with nonsingular polar, squaring its Walsh transform and
changing variables to the input difference makes every nonzero difference
vanish by character orthogonality. All transform magnitudes are sqrt(t).
Consequently the specified H D_Q/sqrt(t) bases, plus Id, really are mutually
unbiased. This reconstructs their existence for all the claimed dimensions.

## 2. The actual Boolean witness, not a spherical relaxation

The source's displayed five unit vectors have Gram energy exactly 5 for its
displayed switching of A5. Their first two flat sign patterns on coordinates
0,u,v,u+v have binary polar bits 0 and 1 respectively; the second flat pattern
is repeated for fibres 1 and 2. The final two vectors are the same coordinate
vector and require no flatness condition.

Fix any nonzero u. The functionals `(B0+B1)(u,.)` and `(B0+B2)(u,.)` are
nonzero, and their difference `(B1+B2)(u,.)` is nonzero. Over F2 they are
therefore linearly independent. Setting both equal to one has exactly t/4
solutions v. Neither 0 nor u solves it, since all polars are alternating.
Thus the two-dimensional support is genuine. With c=B1(u,v), the three
restricted polar bits are c+1,c,c.

Multiplying all five vectors by the common phase c alpha beta preserves every
Gram entry. After multiplication by their own source quadratic phases, the
first three flat patterns all have polar bit one. Their unnormalized Walsh
sums on four points are therefore +-2, and the vector magnitudes are 1/2;
the complete inverse-Walsh outputs have entries +-1. For a coordinate vector,
the same conclusion follows immediately from a Walsh column. This proves
the literal Boolean equation x_i=sqrt(t)U_i w_i.

It follows that U_i^T x_i=sqrt(t)w_i, and the cross energy is
`sqrt(t) sum_{i<j}a_ij (sqrt(t)w_i)·(sqrt(t)w_j)=5t^(3/2)`.
The powers of t and the hollow-energy factor are correct.

## 3. Why every Heisenberg state stays Boolean and why D_i cancels

The common map `T_(p,l)w(z)=(-1)^(l·z)w(z+p)` is an orthogonal signed
permutation, so all pair Gram products are unchanged. Translation turns the
support into an affine plane; restricting a quadratic to its translate changes
only linear and constant coefficients, not its polar bit. Modulation also
adds only a linear term. Thus every orbit vector still meets the exact
two-dimensional bent condition, or remains a coordinate vector.

For each unit w, averaging the modulation eliminates every off-diagonal entry
of `(Tw)(Tw)^T`, by binary character orthogonality. Averaging translations
then replaces the diagonal by the constant 1/t. Hence its covariance is I/t,
and the corresponding Boolean x_i covariance is I. Every fixed hollow D_i
has zero mean quadratic energy. All orbit states have the same cross energy,
so at least one state has total full energy at least 5t^(3/2).

This uses an actual finite Boolean distribution with an exact covariance;
there is no unsupported cancellation assumption about child or bridge maxima.

## 4. Unused-reference factorization, including the identity basis

If all five selected bases already have form H D_Q/sqrt(t), no reference
change is needed. If Id is selected, the argument assumes an unused quadratic
basis mutually unbiased to every selected basis. Such a reference exists in
the displayed complete family at t>=16, since it has t/2+1>=9 bases. The proof
does not claim it for an arbitrary unextendible five-basis family containing Id.

Right multiplication by U_b^T leaves every cross product U_i U_j^T unchanged.
For a selected quadratic basis, the resulting matrix is `H D_Q H/t`, with
Q=Q_a+Q_b having nonsingular polar B. Completing the finite quadratic square
shows that its entries are

```math
{\gamma\over\sqrt t}(-1)^{Q^*(r+s)},\qquad
Q^*(z)=Q(B^{-1}z),\quad\gamma\in\{+1,-1\},
```

up to harmless linear/constant conventions. Its polar is B^(-1), which is
invertible. Expanding Q^*(r+s) separates a row phase, a column phase, and the
invertible pairing `r^T B^(-1)s`. The latter is absorbed by an invertible row
permutation. Therefore this is a signed-row-permutation copy of H D_R/sqrt(t).
For selected Id the relative matrix is D_b H/sqrt(t), immediately a row-switched
Walsh matrix. These physical cube automorphisms preserve the full cap and
merely conjugate arbitrary D_i to other hollow full signings.

The resulting five quadratic column phases still have nonsingular pair
differences: row signed permutations preserve mutual-unbiasedness, and the
quadratic Walsh criterion above then applies. Thus Section 2 is valid for
all these transformed bases. This checks the all-selected-subsets scope.

## 5. Fresh exact replay

`tmp/flatify_adversary_2026_09_07_kerdock_plateaued_replay.py` imports only the
candidate basis generator, not the proposing witness program. At t=16 it:

- checks all 630 relative Hadamards have quadratic column phases and that
  their row quotients are precisely a permutation of all sixteen signed
  characters, explicitly checking the factorization used above;
- constructs witnesses for all 126 subsets using u=3 and the largest valid
  v, different from the proposing program's choices;
- verifies divisibility before every integer normalization, all Boolean
  coordinates, and every exact energy 320;
- checks all 256 Heisenberg states for three subsets, including the earlier
  heuristic miss, for Booleanity, Gram preservation, energy, and covariance;
- fills arbitrary random hollow full D_i in thirty cases and verifies the
  exact orbit average is still 320 and an attaining state exists.

All checks passed. The complete parameters, witnesses, and completion results
are saved in the adjacent JSON. The result is independent of the separate
dimension-16 affine-equivalence screen and genuinely applies at growing t.

At N=5t the forced normalized cap is 1/sqrt(5), above the pentagon's
row-normalized coefficient 2/5. This invalidates the specified finite-seed
quadratic-phase transfer, but does not settle original convergence or the
selected-large-child flatification conjecture.
