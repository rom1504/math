# Two-fiber quotient escape and exponential response obstruction

Date: 2026-07-31. This is an agent-authored audit of the smallest genuinely
nonhomogeneous coherent-configuration candidate.

## 1. Exact two-fiber signing and quotient

Partition the vertices into fibers of sizes `m,n` and write

```math
S=\begin{pmatrix}A&C\\C^{\mathsf T}&B\end{pmatrix}.    \tag{TF1}
```

Assume the blocks are relation-regular:

```math
\begin{aligned}
A\mathbf1_m&=a\mathbf1_m,&
B\mathbf1_n&=b\mathbf1_n,\\
C\mathbf1_n&=c\mathbf1_m,&
C^{\mathsf T}\mathbf1_m&=d\mathbf1_n.
\end{aligned}                                         \tag{TF2}
```

Here `A,B` are zero-diagonal signings and `C` is a rectangular sign matrix.
Counting the entries of `C` gives

```math
mc=nd.                                                 \tag{TF3}
```

On the orthonormal fiber-constant basis, the quotient of `S` is

```math
T=\begin{pmatrix}
a&c\sqrt{m/n}\\
d\sqrt{n/m}&b
\end{pmatrix}
=\begin{pmatrix}a&\gamma\\\gamma&b\end{pmatrix},
\qquad \gamma=c\sqrt{m/n}.                            \tag{TF4}
```

For a block-constant Boolean spin `(sigma 1_m,tau 1_n)`, the energy is

```math
H={ma+nb\over2}+\sigma\tau\,mc.                     \tag{TF5}
```

Maximizing the absolute value over the relative sign proves the exact
quotient lower bound

```math
\boxed{
\operatorname{cap}(S)
\ge {1\over2}|ma+nb|+|mc|.}                           \tag{TF6}
```

This is the complete trivial-channel obligation for two fibers.

## 2. Two fibers can genuinely cancel the `1/2` channel

If `S` is conference,

```math
S^2=(m+n-1)I,                                         \tag{TF7}
```

then the invariant quotient satisfies

```math
T^2=(m+n-1)I_2.                                       \tag{TF8}
```

Unlike the homogeneous case, (TF8) does not force either fiber-constant
Boolean vector to be an eigenvector.

This is not merely a real-valued possibility. Let

```math
m=n=s=2k^2+1,qquad
a=2k,quad b=-2k,quad c=d=1.                         \tag{TF9}
```

All row-sum parities are compatible with sign blocks: `a,b` have the parity
of `s-1`, while `c,d` have the parity of `s`. The quotient is

```math
T_k=\begin{pmatrix}2k&1\\1&-2k\end{pmatrix},qquad
T_k^2=(4k^2+1)I=(2s-1)I.                              \tag{TF10}
```

But (TF6) is only

```math
\operatorname{cap}(S)\ge s=O(m+n),                   \tag{TF11}
```

because the two internal row sums cancel and the cross row sum is one.
Thus a balanced two-fiber coherent configuration can make the entire trivial
Boolean quotient negligible at the project scale. This is the smallest
genuine escape from the homogeneous `1/2` obstruction.

An exact scalable construction target is now explicit. Find sign blocks
satisfying (TF2), (TF9), and

```math
\begin{aligned}
A^2+CC^{\mathsf T}&=(2s-1)I,\\
B^2+C^{\mathsf T}C&=(2s-1)I,\\
AC+CB&=0.
\end{aligned}                                         \tag{TF12}
```

Then (TF1) is a conference signing of order `2s`. The identities are
polynomial-time checkable and use a two-dimensional quotient state. They
supply the ordinary cap upper bound `s sqrt(2s-1)`, but by themselves do not
prove Boolean slack below the `1/2` constant.

The parameter system is nonempty at `k=1`, not merely arithmetically
consistent. With `s=3`, take

```math
A=J_3-I_3,qquad B=-A,qquad C=J_3-2I_3.               \tag{TF12a}
```

These blocks have row sums `2,-2,1`, respectively, and direct multiplication
verifies (TF12). The resulting order-six matrix is a switching/permutation of
the Paley conference signing and has exact cap five. This finite base does
not establish scalable existence, but it verifies the proposed block algebra
and its normalization exactly.

## 3. Fixed row and column sums leave an exponential bridge response

The quotient escape does not make the bridge response low-dimensional. Take
the balanced case `m=n=s`, with `s` odd and every row and column sum of `C`
equal to one. A column is a sign vector `t in {+1,-1}^s` with

```math
\mathbf1^{\mathsf T}t=1.                              \tag{TF13}
```

There are

```math
h_s={s\choose(s+1)/2}
=2^{s-O(\log s)}                                      \tag{TF14}
```

such oriented types. Let `q_t` count the columns of each type. After the
right-hand Boolean signs are optimized, the complete left response is

```math
F_q(x)=\sum_tq_t|x^{\mathsf T}t|,qquad
x\in\{+1,-1\}^s.                                     \tag{TF15}
```

Bowlin's exact bipartite response theorem says that the matrix taking all
column-type counts modulo negation to the optimized responses, after the
total-count row is adjoined, is invertible and has rank `2^(s-1)`. Restrict
that matrix to the types (TF13). Its columns remain linearly independent.
Consequently

```math
q\longmapsto\bigl(\sum_tq_t,(F_q(x))_x\bigr)           \tag{TF16}
```

has rank `h_s`.

Fixing the `s` row sums and the total number of columns imposes at most
`s+1` independent linear conditions on `q`. Therefore the affine family of
fractional type-count vectors with the **same two-fiber quotient** still has
response dimension at least

```math
\boxed{h_s-s-1=2^{s-O(\log s)}.}                      \tag{TF17}
```

The dimension assertion is genuine for the convex relaxation: the uniform
vector `q_t=s/h_s` is strictly positive, has total `s`, and by symmetry has
every row sum equal to one. Thus a neighborhood in all kernel directions is
feasible.

There is also an exact integer obstruction at square size. For every type
`t` in (TF13), take the circulant matrix whose first column is `t`. Every row
and column sum is one. Two types in different cyclic-shift orbits give
different type-count vectors, and (TF16) gives different exact response
profiles. Hence bridges with the same quotient realize at least

```math
{1\over s}{s\choose(s+1)/2}=2^{s-O(\log s)}           \tag{TF18}
```

distinct optimized responses. This count does not rely on fractional or
large-multiplicity bridges.

Equations (TF17)--(TF18) are the precise obstruction: balancing all four quotient row
sums removes the trivial channel but does not compress the optimized bridge
response by more than polynomially many dimensions.

## 4. Research judgment

The two-fiber quotient (TF9)--(TF12) is a genuine algebraic candidate and is
strictly different from fixed tensor continuation. It shows that a
bounded-fiber no-go theorem based only on block-constant spins is false.

However, a proposed composition theorem retaining only

```text
(m,n,a,b,c,d)
```

cannot control the parent Boolean cap: bridges with exactly those data carry
the exponential response freedom (TF17). Supplying the full response (TF15)
is full bridge optimization in compressed notation. The conference
identities (TF12) avoid that optimization only for the spectral upper bound,
which returns the existing `1/2` constant.

Therefore a useful two-fiber mechanism needs one additional theorem-level
ingredient:

> Prove, for an explicitly generated solution of (TF12), a uniform Boolean
> inequality on the nontrivial modules giving `o(N^(3/2))` distance from the
> desired landing constant, without retaining (TF15).

The falsifier is equally concrete: an explicit solution family of (TF12)
with a finite nontrivial-module witness that amplifies to a constant above
the desired one, or two solutions with the same bounded algebraic parameters
and a leading cap separation.

At present, two fibers successfully evade the regular Boolean channel but do
not yet give a bounded-complexity cap state. The exact surviving candidate is
(TF12); the exact obstruction to quotient-only composition is (TF17).
