# Archive collision: FR.5 versus the universal affine shell

**Verdict.**  The packing conclusion of FR.5 is strictly subsumed by
Theorem 36.7 / Corollary LA.4 already in the repository.  The fractional
reservoir mechanism may remain independently useful as a finite-anchor
extension principle, but it is not the first proof of a growing
energy-scale shell packing.

Let `A` be any hollow signing with `Q(A)<=C_+n^(3/2)`, and choose an absolute
ground state with the independent global energy orientation required in
Theorem 36.7.  In the notation of LA.1, take the largest even

```math
k\le c\sqrt n
```

for any fixed `c>0`.  The `2^k` affine masks all have absolute deficit at
most

```math
\Delta_k={4kQ(A)\over n}+2k(k-1)=O_(c,C_+)(n)=o(n^{3/2}).
                                                               \tag{FC.1}
```

Color every mask by a sign orienting its energy.  One color class has at
least `2^(k-1)` elements.  Greedily deleting mask-Hamming balls of radius
`floor(k/4)-1` leaves

```math
\exp(\Omega_c(k))=\exp(\Omega_c(\sqrt n))          \tag{FC.2}
```

masks in one common energy orientation, at mutual mask distance
`d in [k/4,k]`.  Since `k<n/2`, the corresponding oriented augmented cuts
have exact pairwise Hamming, and projective Hamming, distance

```math
d(n-d)=\Theta_c(n^{3/2}).                           \tag{FC.3}
```

Using the established constants

```math
0<C_-n^{3/2}\le M_n\le C_+n^{3/2},
```

(FC.3) is `Theta(M_n)` with fixed constants.  This applies to every
bounded-cap signing, hence in particular every exact minimizer.  The
orientation, projective cardinality, and coding argument were independently
audited in
`drafts/local_field_affine_shell_algebra_audit.md`, Sections 1 and 4.

By comparison, FR.5 gives only

```math
L_n=\left\lfloor {\alpha\log n\over\log\log n}\right\rfloor
```

positive shell words at deficit `o(M_n)` and pairwise distance
`(1/4-o(1))M_n`.  Its fixed constant `1/4` is aesthetically cleaner, but
LA.4 already has a fixed positive constant after choosing `c`; the latter
also has exponentially larger cardinality and the sharper deficit `O(n)`.

Therefore the correct classification is:

```text
FR.1--FR.3: distinct fractional finite-anchor reservoir mechanism.
FR.5 as a growing energy-scale packing: ARCHIVE COLLISION / no frontier gain.
Any checkpoint calling FR.5 the first such packing: must be repaired.
```

This collision does not affect the fixed-edge-scale question.  Both
constructions live at projective edge distance `Theta(M_n)=Theta(n^(3/2))`,
which is a vanishing fraction of `binom(n,2)`, and neither by itself gives a
scalar bounded-cap child-response packing.
