# All balanced macro diagonals fail the fixed Sylvester lift

Date: 2026-07-31. This is an agent-authored research report. The finite
search that found the witnesses was heuristic; the theorem below uses only
the saved witnesses and exact integer verification.

## Result

Let `A` be the saved order-14 symmetric conference signing and let `D` range
over **every** diagonal sign matrix with `tr(D)=0`. For the symmetric
Sylvester matrix `H_k`, define

```math
S_D(k)=A\otimes H_k+D\otimes
       \bigl(H_k-\operatorname{diag}(H_k)\bigr).       \tag{AD1}
```

There are exactly

```math
\binom{14}{7}=3432                                      \tag{AD2}
```

such labeled diagonals. Exhaustive enumeration produced, for every one of
them, an explicit Boolean vector `z_D` of length 56 satisfying

```math
\left|H_{S_D(4)}(z_D)\right|\ge 210.                   \tag{AD3}
```

All 3,432 inequalities in (AD3) were independently recomputed from compact
56-bit witnesses by exact integer matrix multiplication. The smallest saved
absolute energy is exactly 210. Since

```math
4\cdot210^2=176400>175616=56^3,                        \tag{AD4}
```

(AD3) is strictly above normalized cap `1/2` without relying on a
floating-point comparison.

## Exact tensor amplification

Write `Delta_k=diag(H_k)` and `Q=A+D`. Then

```math
S_D(k)=Q\otimes H_k-D\otimes\Delta_k.                 \tag{AD5}
```

Let `L=4^t` and let

```math
v=(-1,-1,-1,1),\qquad H_4v=2v,\qquad
y=v^{\otimes t}.                                      \tag{AD6}
```

Thus `y` is Boolean and `H_L y=sqrt(L)y`. Because both `D` and the diagonal
of `H_4` are balanced,

```math
z_D^{\mathsf T}(D\otimes\Delta_4)z_D
=\operatorname{tr}(D)\operatorname{tr}(\Delta_4)=0.  \tag{AD7}
```

Apply (AD5) to `z_D tensor y`, using
`H_(4L)=H_4 tensor H_L` and
`Delta_(4L)=Delta_4 tensor Delta_L`. Equation (AD7) kills the diagonal
correction and gives the exact identity

```math
\left|H_{S_D(4L)}(z_D\otimes y)\right|
=\left|H_{S_D(4)}(z_D)\right|L^{3/2}.                 \tag{AD8}
```

Consequently, for **every balanced `D`** and every `t>=0`, at order
`N=56L`,

```math
\boxed{
 {\operatorname{cap}(S_D(4^{t+1}))\over N^{3/2}}
 \ge {210\over56^{3/2}}
 =0.501114828585795\ldots >\frac12.}                 \tag{AD9}
```

The current all-order construction bound
`M_N <= (1/2+o(1))N^(3/2)` therefore implies a linear landing gap on these
native orders. If `u_N=cap(S_D)^(2/3)` and `b_N=M_N^(2/3)`, then

```math
u_N-b_N\ge
\left[
 \left({210\over56^{3/2}}\right)^{2/3}
 -2^{-2/3}-o(1)
\right]N
=(0.000936049705\ldots-o(1))N.                       \tag{AD10}
```

## Reproducible certificate

The search program is
`computations/phase2d_audit_balanced_diagonals.py`. It enumerates every
seven-element positive support, runs batched single-spin ascent on both
energy signs, and saves the first threshold witness. Coordinate ascent is
only a discovery method and is not used as a certificate.

The certificate
`computations/results/phase2d_balanced_diagonal_family_audit.json` contains
all 3,432 diagonal masks and Boolean witnesses. The canonical sorted record
hash is

```text
d86922dc01741ecf21729a4c924b2dddc09c449e779e62a0b2a5626eec4a0475
```

The independent verifier
`computations/phase2d_verify_balanced_diagonal_family_audit.py` checks:

1. the source conference identity and source hash;
2. that every balanced diagonal occurs exactly once;
3. that every lifted matrix is a zero-diagonal signing;
4. every saved energy by direct integer multiplication;
5. the strict exact comparison (AD4), record hash, and tensor channel (AD6).

The exact saved-energy distribution is

```text
210: 8, 212: 1, 214: 31, 216: 248, 218: 337, 220: 551,
222: 234, 224: 1024, 226: 324, 230: 96, 232: 313, 234: 265.
```

## Research judgment and scope

This upgrades the earlier one-diagonal example to a complete finite-family
obstruction. Choosing the balanced macro diagonal cannot repair the common
Sylvester fiber construction: every possible choice admits an entangled
order-56 witness that tensorizes to a positive leading-constant excess.

The result does not reject scale-dependent micro algebras, non-Kronecker
fusion, or a construction whose state changes with the lifted order. It does
establish a sharp screening rule for fixed-fiber proposals: enumerate all
bounded state choices at the first nontrivial fiber level and test for a
Boolean witness above the target leading constant before attempting a
composition theorem.
