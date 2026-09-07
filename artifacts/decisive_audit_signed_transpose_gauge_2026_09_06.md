# Signed-transpose gauge: exact ambient invariance, not Boolean-source invariance

Date: 2026-09-06. A scoped algebraic observation for the actual grouped
weave of `decisive_bridge_vector_seed_weave_unitary_obstruction_2026_09_06.md`.

For a full symmetric sign seed `A` of order `k`, define
`(R_A z)_ij=A_ij z_ji` on real `k`-by-`k` arrays. For `i<j` choose
`d_ij=1,d_ji=A_ij`, and on the diagonal choose `d_ii=1`. Let
`D_A z=(d_ij z_ij)`. Then

```math
 (D_A R_A D_A z)_{ij}=z_{ji}\quad(i\ne j),\qquad
 (D_A R_A D_A z)_{ii}=A_{ii}z_{ii}.                        (1)
```

Thus all such reflections with the SAME sign diagonal are conjugate
by coordinate sign gauges. Consequently a physical kernel problem whose
admissible source class is invariant under every such gauge has the same
optimized value for every offdiagonal seed. This is an exact statement
about that source class; it is not a claim about the literal weave source.

The invariance hypothesis is substantial. Let

```math
 H=\begin{pmatrix}
 1&1&1&1\\1&-1&1&-1\\1&1&-1&-1\\1&-1&-1&1
 \end{pmatrix},\qquad x=(1,1,1,-1)^T.
```

Then `h=H^T x=(2,2,2,-2)^T`. Flipping only its first spectral coordinate
gives `h'=(-2,2,2,-2)^T`, and

```math
                     Hh'/4=(0,0,0,-2)^T,                 (2)
```

which is not Boolean. This single-coordinate row gauge occurs in the
triangular construction (1) for a seed with just one negative offdiagonal
edge. Hence even a tiny gauge needed to remove the signed transpose may
take a literal Hadamard spectrum out of its allowed source set.

The local Hilbert-norm contraction really is seed blind, as proved in
the bridge note, because unitary edge absorption preserves norms without
requiring source invariance. The positive transport expression keeps the
source domain and therefore is NOT shown seed blind by (1). Any argument
using (1) there must first prove the missing domain invariance; (2) refutes
that claim for fixed actual Hadamard Boolean spectra in general.
