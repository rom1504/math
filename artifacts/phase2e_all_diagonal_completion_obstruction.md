# No diagonal completion repairs the common Sylvester lift

Date: 2026-07-31. This is an agent-authored research report. It extends the
balanced-diagonal audit to every diagonal sign choice.

## Exact finite-family theorem

Let `A` be the saved exact order-14 conference minimizer. For an arbitrary
diagonal sign matrix

```math
D=\operatorname{diag}(d_1,\ldots,d_{14}),\qquad d_i\in\{-1,1\},             \tag{AC1}
```

define the common-Sylvester lift

```math
S_D(k)=A\otimes H_k+D\otimes
       \bigl(H_k-\operatorname{diag}(H_k)\bigr).                            \tag{AC2}
```

There are `2^14=16384` labeled choices of `D`. Exact exhaustive finite-family
verification establishes:

```math
\boxed{
 \text{For every }D\in\{\pm1\}^{14},\text{ there is }z_D\in\{\pm1\}^{56}
 \text{ with }|H_{S_D(4)}(z_D)|\ge210.}                                    \tag{AC3}
```

No choice was left as a heuristic holdout. The discovery procedure was
heuristic single-spin ascent, but (AC3) follows solely from the 16,384 saved
witnesses and direct exact integer evaluation. Every diagonal mask occurs
exactly once. The minimum saved absolute energy is 210.

The exact comparison

```math
4\cdot210^2=176400>175616=56^3                                      \tag{AC4}
```

shows that every member is strictly above normalized cap `1/2`.

## Why no balance assumption is needed

Put `Delta_k=diag(H_k)` and `Q=A+D`, so

```math
S_D(k)=Q\otimes H_k-D\otimes\Delta_k.                                 \tag{AC5}
```

For the order-four Sylvester matrix,

```math
\operatorname{tr}(\Delta_4)=\operatorname{tr}(H_4)=0.                  \tag{AC6}
```

Because a Boolean vector has coordinate squares equal to one, (AC6) gives,
for every `D` and every `z`,

```math
z^{\mathsf T}(D\otimes\Delta_4)z
=\operatorname{tr}(D)\operatorname{tr}(\Delta_4)=0.                   \tag{AC7}
```

Thus the earlier use of `tr(D)=0` was unnecessary; the balanced Hadamard
diagonal alone kills the correction.

Let `L=4^t`, `v=(-1,-1,-1,1)`, and `y=v^(tensor t)`. Since `H_4v=2v`,
`H_Ly=sqrt(L)y`. Using
`H_(4L)=H_4 tensor H_L` and
`Delta_(4L)=Delta_4 tensor Delta_L`, equations (AC5)--(AC7) give

```math
|H_{S_D(4L)}(z_D\otimes y)|
=|H_{S_D(4)}(z_D)|L^{3/2}.                                             \tag{AC8}
```

Consequently, uniformly over all diagonal choices and all `t>=0`, at
`N=56L`,

```math
{\operatorname{cap}(S_D(4^{t+1}))\over N^{3/2}}
\ge {210\over56^{3/2}}
=0.501114828585795\ldots>\frac12.                                     \tag{AC9}
```

Even choosing a different one of the 16,384 macro diagonals at each scale
cannot evade (AC9). Relative to the project's all-order upper construction
`M_N<=(1/2+o(1))N^(3/2)`, this is the same certified linear landing gap

```math
u_N-b_N\ge(0.000936049705\ldots-o(1))N.                               \tag{AC10}
```

## Certificate and verification

The generalized search program is
`computations/phase2d_audit_balanced_diagonals.py`, invoked with
`--diagonal-family all`. The 2.3 MB certificate is
`computations/results/phase2e_all_diagonal_family_audit.json`.

The independent exact verifier is
`computations/phase2d_verify_balanced_diagonal_family_audit.py`. It checks
the source hash and conference identity, all 16,384 masks, every lifted
signing and witness energy, the strict threshold comparison, the zero-trace
and Boolean-eigenvector tensor identities, and the canonical sorted-record
hash

```text
760daa2b18a0e355f857f9e65b355f3e9bd0f7b33cb8b3c731357377120ad9e0
```

The exact saved-energy distribution is

```text
210: 15, 212: 46, 214: 252, 216: 1322, 218: 2983, 220: 1982,
222: 1871, 224: 3336, 226: 2813, 228: 305, 230: 105,
232: 955, 234: 297, 240: 102.
```

## Research judgment

This decides the entire diagonal-completion freedom of the common
Sylvester microblock: **none of its diagonal completions is a near-optimal
landing family**. The obstruction is scalable and uniform, not a collection
of unrelated finite caps.

The remaining escape routes must change something structural: the micro
algebra, the common-fiber Kronecker rule, or the fixed order-14 macro seed.
Searching further over diagonal completions of this same lift has no
mathematical leverage.
