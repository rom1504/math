# Deterministic good children can be planted in near-minimizing parents

Date: 2026-09-06. This is an actual finite modification theorem and a
counterpoint to the random-selector obstructions. It does not show that
an EXACT global minimizer contains a prescribed child, and does not
prove convergence or supply a lossless upward construction.

## 1. Finite planting lemma

Let `A` be any symmetric hollow sign matrix of order `D`, and let `B`
be any prescribed symmetric hollow sign matrix of order `n<D`. Write
`p=n/D` and `K=pi/(2 asinh(1))`. There is a principal set `S` of size
`n` and a hollow signing `A'`, identical to `A` outside `S`-by-`S`,
such that `(A')_S=B` exactly and

```math
Q(A')\le\left(1+\frac{2Kp}{1-p}\right)Q(A)+Q(B).          (1)
```

For `n<=D/2`, the simpler estimate is

```math
Q(A')\le(1+4Kp)Q(A)+Q(B).                               (2)
```

Proof. Use the independently proved finite simultaneous diagonal
majorant, whose trace is at most `K beta(A)<=4K Q(A)`.
Spectral deletion with `epsilon=1-p` leaves a core of size at least
`pD=n`, with operator norm at most

```math
4K Q(A)/[(1-p)D].
```

Choose any `n` coordinates `S` inside it. Then

```math
Q(A_S)\le\frac n2\|A_S\|_{op}
       \le\frac{2Kp}{1-p}Q(A).                          (3)
```

Replace that principal block by the prescribed `B`. For every Boolean
spin on the full order, the energy change is exactly the new block
energy minus the old block energy. The triangle inequality gives
`Q(A')<=Q(A)+Q(A_S)+Q(B)`, proving (1). Alternatively use fixed
`epsilon=1/2` to obtain (2).

The vertices outside the spectral core are NOT deleted from the final
parent. The core is used only to locate the block to overwrite. This
is why the retained-core fraction need not approach one, and the loss
is `O(p)`, not `O(sqrt p)`.

The majorant and its constant are proved in Section 1 of
`resumed_bound_audit_minimal_proof_2026_09_06.md`. Its proof uses
polarization, an explicit Gaussian/tensor Grothendieck rounding, finite
SDP duality, and deletion of the largest diagonal weights. It has no
minimizer or convergence hypothesis.

## 2. Every sublinear good child can occur in near-minimizing parents

For every `D`, choose an actual minimizer `A_D`, so `Q(A_D)=M_D`.
Let `n=n(D)=o(D)` and prescribe any children `B_n` satisfying
`Q(B_n)<=C1 n^(3/2)` with fixed `C1`. The elementary universal upper
bound gives `M_D=O(D^(3/2))`. Applying (2) yields

```math
0\le Q(A'_D)-M_D
\le4K(n/D)M_D+Q(B_n)=o(D^{3/2}).                         (4)
```

Thus `A'_D` is asymptotically minimizing in the precise additive sense
`Q(A'_D)=M_D+o(D^(3/2))`, while containing the prescribed child EXACTLY.
In particular one may prescribe an actual minimizer at every child
order. The construction works for every sublinear sequence, with no
relation between its growth rate and a power or logarithm of `D`.

This conclusion is compatible with the proven typical-selector floor
`2/pi`: the specially planted selector is exceptional. It rules out
inferring a universal deterministic-selector obstruction from the random
one, even if parent quality is required to be asymptotically optimal.

## 3. A fixed positive retention can preserve strict-subhalf examples

Let `c_*<1/2` be the proved all-order upper coefficient, and choose
actual parent and child minimizers at orders `D` and `floor(pD)`.
For fixed `p in (0,1)`, (1) gives planted parents satisfying

```math
\limsup_{D\to\infty}\frac{Q(A'_D)}{D^{3/2}}
\le c_*\left(1+\frac{2Kp}{1-p}+p^{3/2}\right).           (5)
```

The right side remains strictly below one half for all sufficiently
small fixed positive `p`. Hence both the parent family and its prescribed
macroscopic child family can have strict-subhalf caps. No numerical
choice of `p` is needed for this continuity consequence, and none is
asserted optimal.

At fixed positive `p`, (5) is NOT an asymptotic-minimality assertion:
the allowed additive normalized parent loss is a positive constant.
Equation (4), where the loss tends to zero, requires `p->0`.

## 4. The original gap remains

The construction starts with an already-existing large minimizer.
It does not transport a given small minimizer upward at its original
coefficient, and supplies no comparison forcing the original liminf
and limsup to agree. It identifies a genuine deterministic escape
from the new selector rarity theorems, with the finite cost (1).
