# A hidden third return despite small off-diagonal square entries

Date: 2026-09-06. Status: probabilistic construction submitted for audit;
the finite fixed-path identity is exactly replayed. This is an actual
signing counterexample, not a minimizing construction.

There exist symmetric hollow signings `B_n=A_n/sqrt(n-1)` such that

```math
\limsup\|B_n\|_{\rm op}\le2,\qquad
\max_{i\ne j}|(B_n^2)_{ij}|\to0,\qquad
\frac1n\operatorname{Tr}(B_n^3)\to0,
```

but, for a positive linear number of disjoint off-diagonal pairs,

```math
(B_n^3)_{ij}\to\frac18.
```

Thus small off-diagonal entries of `Q=B²` do not by themselves justify
higher-power or general open-cactus delocalization.

## Construction

Let `m=4^r`, let `W=H_m/sqrt(m)` be normalized Sylvester Hadamard, and
let `D` be the diagonal sign matrix of the Boolean quadratic form
`q(x)=sum_l x_(2l)x_(2l+1)` on `F_2^{2r}`. Its elementary product
Fourier transform is bent, so `W D W` is again orthogonal and has all
entries of absolute value `1/sqrt(m)`.

Make a symmetric four-by-four block matrix `B0` of order `4m`, with
each block scaled by `1/2`. Fix its three upper blocks

```math
A_{12}=W,\qquad A_{23}=DW,\qquad A_{34}=WDW.
```

Then `A12 A23=WDW`, `A23 A34=W`, but

```math
A_{12}A_{23}A_{34}=I.
```

For every other off-diagonal block independently choose `P W R` with
independent random diagonal sign matrices `P,R`, and use its transpose
below the diagonal. For each diagonal block choose an independent
`P W P`, which is symmetric. Every block is orthogonal, and every full
entry of `B0` has magnitude `1/(2sqrt(m))`. The block norm bound gives
`||B0||op<=2` deterministically.

Multiply by `2sqrt(m)`, hollow the resulting full sign matrix, and
normalize by `sqrt(4m-1)` to obtain `B`. Its operator distance from `B0`
is `O(m^{-1/2})`. Matrix powers through three have the same vanishing
operator error, hence the same uniform entrywise limits.

## Delocalization of the other short paths

The diagonal blocks of `B0²` equal `I` exactly: each is one quarter
of four orthogonal block products. Every off-diagonal block is a sum
of four two-step products. The two products consisting only of the
fixed consecutive edges are flat as displayed above. Every other
two-step product has at least one independent sign gauge in its
interior. Each fixed entry is a Rademacher linear or quadratic chaos
with second moment `O(1/m)`, apart from flat deterministic terms.
Fixed-degree hypercontractivity and a union bound therefore make all
of these entries `o(1)` simultaneously with probability tending to one.

For the `(1,4)` block of `B0³`, the path `1→2→3→4` contributes `I/8`.
Every other three-step path is either a backtrack reducing to a flat
block, or contains an independent gauge. For a path containing a
random off-diagonal block `P W R`, condition on the other blocks and
use one of its two gauges: one adjacent factor is a flat `W` row or
column, and the opposite factor is orthogonal, giving entry variance
`O(1/m)`. If a diagonal signed Hadamard is interior, the same entry is
a degree-two sign chaos; its coefficient-square sum is `O(1/m)`.
Paths with two repeated adjacent diagonal blocks reduce using
`(PWP)²=I`. The remaining finite cases obey the same bounded-degree
moment estimate. Uniformly, an entry's non-flat random part has
`L^p` norm at most `C p^2/sqrt(m)` for `p>=2`. Choosing `p` a suitable
multiple of `log m` proves simultaneous `o(1)` bounds over all entries
and finitely many path types.

The diagonal blocks of `B0³` have no fixed all-distinct triangle:
every triangle among three different block vertices contains one of
the independent random blocks. Repeated-vertex terms reduce to flat
blocks or the same sign-chaos bound. Their diagonal entries therefore
vanish uniformly, proving the cubic-trace assertion.

With positive probability all these events hold, so choosing one good
gauge realization for each `m=4^r` gives the claimed deterministic
sequence. The code
`computations/continued_feedback_hidden_third_return_2026_09_06.py`
checks the fixed three-step identity by integer arithmetic and supplies
finite stress tests. Its floating outputs do not replace the moment
argument above.

The construction rules out this particular covariance-to-higher-return
shortcut. It does not refute a theorem imposing delocalization of every
required higher return separately, or a minimizer-specific argument
that excludes this signing family.
