# Adversarial audit: multichannel response entropy

**Verdict.**  The binary `Omega(Dk)` response-packing theorem in
[`phase3_multichannel_response_entropy.md`](phase3_multichannel_response_entropy.md)
is valid, including its claim that no channel-labelled queries are used.
The parallel finite-field theorem in
[`phase3_qary_multichannel_holonomy.md`](phase3_qary_multichannel_holonomy.md)
is also valid with its stated constants.  This audit distinguishes the
actual invariances from a tempting but false `GL(k)` quotient.

## 1. Binary theorem: independent proof audit

For

```math
F_V(u)=\min_z(2|z|+|u+Vz|)
```

and `C=im V`, dropping the first term gives `F_V>=d_C`.  Choosing an
independent subset of the input columns shows that a nearest codeword has a
representation with at most `rank(V)` letters, so

```math
d_C\le F_V\le d_C+2\operatorname{rank}(V).    \tag{A.1}
```

No factor four is needed in the packing argument.  If `C,C'` are distinct
`k`-subspaces of a linear host `C_0`, choose `c in C\C'`.  Host minimum
distance gives `d(c,C')>=d(C_0)`.  At this *one-sided* witness,

```math
F_C(c)\le2k,
\qquad
F_{C'}(c)\ge d(C_0),                           \tag{A.2}
```

and hence separation `d(C_0)-2k`.  This is robust under every choice of
basis for the two subspaces.  The product formula indeed gives

```math
{r\brack k}_2\ge2^{k(r-k)}:                   \tag{A.3}
```

after extracting `2^{r-k}` from each factor, the remaining numerator
factor is at least the denominator factor because `r>=k`.  A random
`r=floor(rho D)` dimensional binary subspace with
`rho<1-H_2(delta)` supplies the required host by the usual union bound.
Thus `2k+2epsilon D<delta D-o(D)` forces distinct summary states and gives
`k(r-k)=Omega(Dk)` bits.

The exact collision classification is also correct.  A column of weight at
most two is replaced by coordinate letters at cost at most two, and
duplicates are useless.  Conversely, for `|u|>=3`, a value `F_V(u)=2` can
only use one cost-two channel with zero residual.  Therefore the exact
profile recovers the set of distinct columns of weight at least three.

This shows why dividing matrices by right `GL(k,2)` is invalid.  Coefficient
Hamming weight is not invariant under a general basis change.  Only channel
permutations preserve the paired atom set.  At scale `o(D)`, the span is an
approximate quotient when `k=o(D)` by (A.1), but it is not an exact quotient.
The lower bound deliberately chooses only one basis per subspace, so it
survives even this coarser `GL(k)` forgetting.

Finally, the two-fragment interpretation has no hidden labelled query.  If
the quotient basis contains a zero lift and a `v_j` lift in channel `j`, a
quotient-zero selection uses both or neither, giving exactly the displayed
profile.  Both individual fragments are shear-trivial.  The query is only
the kernel endpoint `(u,0)`; the quotient basis supplies structural pairing
of atoms, not a queried channel label.

## 2. Finite-field theorem

The scalar-fibre reduction in
`phase3_qary_multichannel_holonomy.md` is sound.  After summing all letters
in channel `j`, a nonzero total coefficient on the lifted side forces its
opposite nonzero coefficient on the zero-lift side, costing at least two;
a zero total has zero kernel holonomy and all such letters can be deleted.
This proves the exact normal form over every finite field, including
characteristic two.  The local shear fixes `W` and removes all offsets.

The random-host estimate is conservative but correct: the low-weight ball
has size at most

```math
q^{(H_2(1/8)+1/8)D},
```

while membership of a fixed nonzero vector in a random `floor(D/4)`-space
is at most `q^{-3D/4+O(1)}`; since
`H_2(1/8)+1/8<3/4`, the expected number tends to zero.  For
`k<=D/32`, `r-k>=3D/16` for all sufficiently large `D`, so the Gaussian
count is at least `q^{3Dk/16}`.  Minimum distance `>D/8` minus coefficient
cost `2k<=D/16` gives the strict `>D/16` profile gap.  Error below `D/32`
therefore decodes the packing, and the deterministic and Fano bounds have
the stated constants.

The scope statement is accurate: this is a response-entropy lower bound for
scalar-closed Cayley/coset-leader profiles.  It does not assert that every
mixed-holonomy family is incompressible, nor that the weighted generator
metric is `GL(k,q)` invariant.
