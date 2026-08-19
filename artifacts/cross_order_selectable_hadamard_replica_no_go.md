# The bounded selectable-Hadamard core test already fails at order five

Status: **certified finite falsifier, with an exact amplification theorem for
one explicit outer matrix at every dyadic replica order**.  This note concerns
only the unfilled Kronecker core in
[`cross_order_growing_hadamard_replica.md`](cross_order_growing_hadamard_replica.md).
It is not a lower bound on the best completed parent.

Let `H_2` and `H_4=H_2 tensor H_2` be the Sylvester Hadamard matrices.  For an
exact order-`r` cap minimizer `A`, put

```math
R_{k,r}(A)=M(H_k\otimes A)-k^{3/2}M_r.              \tag{1}
```

The proposed bounded-selection test was

```math
\min_{A\in\operatorname{Argmin}M_r/\!\sim}
\min_{k\in\{2,4\}}R_{k,r}(A)\le0.                 \tag{2}
```

Here the first minimum may select the most favorable exact-minimizer class;
`sim` is signed permutation and global sign.

## 1. Exact small-order test

Both signs of every displayed Boolean quadratic form were optimized exactly.
The minimizer-class lists are themselves exhaustive.

| `r` | exact classes | `M_r` | `M(H_2 tensor A)` | `R_(2,r)` | `M(H_4 tensor A)` | `R_(4,r)` |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 1 | 3 | 6 | `6-6 sqrt(2)` | 24 | 0 |
| 4 | 1 | 4 | 12 | `12-8 sqrt(2)` | 32 | 0 |
| 5 | 1 | 4 | 12 | `12-8 sqrt(2)` | 40 | 8 |

Consequently the unique order-five class gives

```math
\boxed{
\min_{k\in\{2,4\}}R_{k,5}(A)
=12-8\sqrt2
=0.686291501015\ldots>0.}                          \tag{3}
```

Thus (2) is false, and order five is the smallest falsifier: it holds at
orders three and four.  Selection among exact-minimizer classes cannot repair
the order-five failure because there is only one class.

The decisive two caps also have a short exact proof independent of the
solver.  A switching-equivalent representative of the unique order-five
class is

```math
A=\begin{pmatrix}
0&-1&1&-1&1\\
-1&0&-1&1&1\\
1&-1&0&1&-1\\
-1&1&1&0&-1\\
1&1&-1&-1&0
\end{pmatrix},
\qquad A^2=5I-J.                                   \tag{3a}
```

For `K=H_k tensor A`, write a Boolean vector `z` in `k` blocks of length
five and let `s_a` be the sum of block `a`.  Every `s_a` is odd.  Since
`H_k^T H_k=kI`, (3a) gives

```math
\|Kz\|_2^2
=k\sum_{a=1}^k(25-s_a^2)
\le24k^2.                                          \tag{3b}
```

Consequently

```math
|H_K(z)|={1\over2}|z^T Kz|
\le {1\over2}\sqrt{5k}\sqrt{24k^2}
={1\over2}\sqrt{120k^3}.                          \tag{3c}
```

The nonzero-support graph of `K` has degree `4k`, which is even.  A single
spin flip therefore changes its energy by a multiple of four.  Moreover the
all-one energy vanishes (`A1=0`), so every energy is in `4 Z`.  At `k=2`,
(3c) and this congruence give `M(K)<=12`; at `k=4` they
give `M(K)<=40`.  The saved Boolean witnesses attain `12` and `40`, proving
both values exactly by elementary arithmetic.

The exact core-only defects in the convergence normalization are

```math
12^{2/3}-2\,4^{2/3}=0.201798588838\ldots,           \tag{4}
```

and

```math
40^{2/3}-4\,4^{2/3}=1.616702553692\ldots.           \tag{5}
```

These positive numbers are not lower bounds on `b_10-2b_5` or
`b_20-4b_5`: filling the zero matching-coordinate edges can cancel part of
the core.

## 2. What the proved completion inequality actually certifies

Equation (18b) of the growing-replica note uses

```math
X=k^{3/2}M_r,
\qquad
C_{k,r}=\sqrt{k(k-1)r(kr+2)\log2}
```

and proves

```math
b_{kr}-kb_r
\le {2\over3}X^{-1/3}(R_{k,r}(A)+C_{k,r})_+.       \tag{6}
```

For the decisive order-five child this gives, without hiding the completion
payment,

```math
b_{10}-2b_5
\le {2\over3}(8\sqrt2)^{-1/3}
 \left(12-8\sqrt2+\sqrt{120\log2}\right)
=2.912190612148\ldots,                              \tag{7}
```

while `k=4` gives

```math
b_{20}-4b_5
\le {2\over3}32^{-1/3}
 \left(8+\sqrt{1320\log2}\right)
=8.031621060076\ldots.                              \tag{8}
```

So the favorable choice in this two-element family is still `k=2`, but (7)
has a fixed positive finite defect.  More importantly, (3) shows that the
core term itself cannot be discarded or signed favorably by selecting only
between these two bounded replica orders.

## 3. Dyadic amplification of the order-five obstruction

The obstruction is not confined to the two seed sizes.  Let

```math
R_q=(J_4-2I_4)^{\otimes d},\qquad q=4^d.
```

This is a symmetric regular Hadamard matrix and
`R_q 1=sqrt(q) 1`.  If `S` is either seed `H_2` or `H_4`, then
`T=R_q tensor S` is a symmetric Hadamard matrix.  Testing `T tensor A` on
the Boolean vector `1 tensor z`, where `z` maximizes `S tensor A`, proves

```math
M(T\otimes A)\ge q^{3/2}M(S\otimes A).              \tag{9}
```

For `S=H_2`, total replica order `K=2q`, (9) yields

```math
M(T\otimes A)-K^{3/2}M(A)
\ge q^{3/2}(12-8\sqrt2),                            \tag{10}
```

with relative excess

```math
\eta_2={3\over2\sqrt2}-1=0.060660171780\ldots.
```

For `S=H_4`, total replica order `K=4q`, it yields

```math
M(T\otimes A)-K^{3/2}M(A)\ge8q^{3/2},              \tag{11}
```

with relative excess `eta_4=1/4`.  Equations (10)--(11) supply one explicit
bad symmetric-Hadamard outer matrix at every dyadic replica order `K>=2`.
They do **not** say that every symmetric Hadamard of that order is bad, so
they do not falsify selection over all outer matrices.

The pressure amplification (21) in the growing-replica note is also
quantitative.  Along the `K=2q` family its linear lower bound becomes positive
when

```math
\beta>{5^{3/2}\log2\over3\sqrt2-4}
=31.938670973644\ldots,                             \tag{12}
```

and along `K=4q` when

```math
\beta>5^{3/2}\log2
=7.749621070722\ldots.                              \tag{13}
```

This is a scalable obstruction to those particular core certificates at
fixed child order five.  It is not an obstruction along selected actual
minimizers with `r` tending to infinity.

## 4. Reproduction and scope

The exact minimizer representatives come from
`computations/results/m{3,4,5}_minimizer_orbits.json`.  Recompute every core
entry and both objective signs with

```text
.venv/bin/python computations/audit_hadamard_replica_cap_selection.py \
  --orders 3 4 5 --ks 2 4 --time-limit 120 --workers 8 \
  --output computations/results/hadamard_replica_cap_selection.json
```

All twelve CP-SAT optimizations return `OPTIMAL`; the saved spin witnesses
independently attain the recorded objectives.  OR-Tools supplies a solver
certificate, not a standalone proof object.

The exact conclusion is narrow but decisive: **bounded selection between the
two canonical replica orders does not make the actual-minimizer core excess
nonpositive**.  The remaining viable forms of the replica route must either
select a genuinely different outer Hadamard using child structure or prove
cancellation in the sign completion rather than require a favorable unfilled
core.
