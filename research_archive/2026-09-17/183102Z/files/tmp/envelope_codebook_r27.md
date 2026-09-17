# Wave 27 Route 2: a row-regular block-coset codebook

## Status

The random-diagonal row-regularization theorem, the exact compressed-envelope
identity, the greedy coset-cover lemma, and the flip formulas below are
**proved**.  The `A_6/A_8/A_9` block-cover census is an **exact finite
enumeration**; displayed decimals only evaluate explicit radicals.  The
compressed Boolean-norm lower bound aligned with `Y_A` is **open**.  Thus this
memo does not prove (10.808), (10.795), or convergence.

The useful advance is that the row-square and entropy parts of (10.808) can be
met simultaneously by a whole structured codebook, for every exact minimizer.
The sole remaining issue for this codebook is a clean, selector-aligned
compressed norm inequality.

## 1. Exact block-coset notation

Partition `[n]` into `k` nonempty blocks `B_1,...,B_k`, of maximum size `b`,
and let `P` be their `n x k` incidence matrix.  For a vertex signing
`epsilon in {+-1}^n`, put `D=diag(epsilon)`.  The associated oriented
block-coset codebook is

```math
\mathcal C(D,P)
=\{d_{\sigma,z}:x=DPz,\ \sigma\in\{\pm1\},\ z\in\{\pm1\}^k\}/(z\sim-z).
\tag{R27.1}
```

It has at most `2^k` oriented cuts.  Every row square in this codebook is

```math
R_2(d_{\sigma,z})=\lVert ADPz\rVert_2^2.
\tag{R27.2}
```

For a selector `S`, with indicator `xi_S`, define the centered selector
matrix

```math
H_S=A\circ(\xi_S\xi_S^{\mathsf T}-p_2\mathbf1\mathbf1^{\mathsf T}).
\tag{R27.3}
```

The zero diagonal of `A` makes the convention on the diagonal of the second
factor irrelevant.  Direct substitution gives the exact aligned-envelope
identity

```math
\boxed{
\max_{d\in\mathcal C(D,P)}X_d(S)
=Q\!\left(P^{\mathsf T}DH_SDP\right),
}
\tag{R27.4}
```

where for the possibly nonzero-diagonal compressed matrix
`Q(B)=max_z |z^T Bz|`.  This is genuinely coupled to `Y_A`; it is not a
scalar sign-tail statement.

## 2. A complete row-regular codebook theorem

**Theorem (proved).**  For every symmetric zero-diagonal sign matrix `A` and
every fixed block partition as above, some diagonal sign matrix `D` satisfies

```math
\boxed{
\max_z\lVert ADPz\rVert_2^2
\le 2k\log(4(n+k))
\max\{\lVert A\rVert_{\rm op}^2,b(n-1)\}.
}
\tag{R27.5}
```

Consequently, let `A` be an exact order-`n` minimizer and fix
`0<c<1/4`.  Choose a balanced partition with

```math
k=\left\lfloor\frac{n^{3/4-c}}{C_0\log n}\right\rfloor,
\qquad
b=O(n^{1/4+c}\log n).
\tag{R27.6}
```

For all sufficiently large `n`, `b=o(sqrt(n))`.  Since exact minimality gives
`||A||_op^2<=2q_n=O(n^(3/2))`, (R27.5) yields one codebook for which

```math
|\mathcal C|\le2^k=\exp\{O(n^{3/4-c})\},
\qquad
\max_{d\in\mathcal C}R_2(d)=O(n^{9/4-c}).
\tag{R27.7}
```

Thus no extra row pruning or costly-cut extraction is needed.

**Proof.**  Write `g(i)` for the block containing `i`, and let `e_a` denote a
coordinate vector in `R^k`.  With independent Rademacher `epsilon_i`,

```math
ADP=\sum_{i=1}^n\epsilon_i Z_i,
\qquad
Z_i=(Ae_i)e_{g(i)}^{\mathsf T}.
```

The two rectangular matrix-series variances are exactly

```math
\sum_iZ_iZ_i^{\mathsf T}=A^2,
\qquad
\sum_iZ_i^{\mathsf T}Z_i
=\operatorname{diag}_{a}\{|B_a|(n-1)\}.
\tag{R27.8}
```

Hence the self-adjoint dilation has variance

```math
v=\max\{\lVert A\rVert_{\rm op}^2,b(n-1)\}.
```

The matrix Rademacher-series bound gives

```math
\Pr\{\lVert ADP\rVert_{\rm op}\ge t\}
\le2(n+k)e^{-t^2/(2v)}.
```

Taking `t^2=2v log(4(n+k))` leaves probability at most `1/2`, so some `D`
obeys that bound.  Finally `||z||_2^2=k` for every block word, proving
(R27.5).  Substitution of (R27.6) proves (R27.7).  `square`

## 3. The exact remaining alignment lemma

For the `D,P` furnished above, (10.808) reduces to precisely

```math
\boxed{
Q\!\left(P^{\mathsf T}DH_SDP\right)
\ge Y_A(S)-t
\quad\text{for every }S\in\binom{[n]}m,
\qquad t=O(n^{3/2-c}).
}
\tag{R27.9}
```

The row theorem proves that the same `D` may be chosen row-good, but it does
not prove (R27.9).  Establishing (R27.9) for an `A`-adapted partition and one
row-good diagonal is now an exact sufficient lemma with the correct entropy
and row scales.  It is substantially more specific than asking for an
unspecified small codebook.

There is also a fractional/greedy version which may be easier.

**Greedy coset-cover lemma (proved).**  Let `mu` be any law supported on
diagonals `D` for which (R27.5) holds.  If, uniformly for every selector,

```math
\mu\left\{D:
Q(P^{\mathsf T}DH_SDP)\ge Y_A(S)-t
\right\}\ge\delta,
\tag{R27.10}
```

then at most

```math
L=\left\lceil\frac{\log\binom nm+1}{\delta}\right\rceil
\tag{R27.11}
```

such cosets cover the whole slice.  Their union has at most `L 2^k` cuts,
all satisfying (R27.5).  In particular, if
`-log(delta)=O(n^(3/4-c))`, then the union has the size required in (10.808).

Indeed, `L` independent cosets miss a fixed selector with probability at most
`exp(-delta L)`; the union bound is below one.  Equivalently, the usual greedy
set-cover averaging proves the deterministic statement.  The missing input
is (R27.10), not the cover argument.

## 4. Why arbitrary block rounding is insufficient

Let `x'=x odot h`, where `H={i:h_i=-1}`, and use unordered shores

```math
W_H(T)=\sum_{i\in H\cap T,\ j\in T\setminus H}A_{ij}x_ix_j.
```

For either orientation `sigma`, exact expansion gives

```math
\boxed{
X_{d'}(S)-X_d(S)
=-4\sigma[W_H(S)-p_2W_H([n])].
}
\tag{R27.12}
```

Also

```math
\boxed{
R_2(x')-R_2(x)
=-4\sum_{i\in H}x_i(A^2x)_i
+4\lVert A[:,H]x_H\rVert_2^2.
}
\tag{R27.13}
```

In particular, the worst-case estimate from (R27.12) is
`|Delta X|<=8|H|n`.  A Hamming projection argument alone therefore needs
every relevant child word within `O(n^(1/2-c))` flips of the block coset to
stay inside the target tolerance.  No such simultaneous approximation of all
child grounds is known.  Formula (R27.13) also shows that row square is not
monotone under nearest-block projection.  The random-diagonal theorem avoids
this row problem, but not the aligned norm problem (R27.9).

## 5. Exact finite census

The checker `tmp/envelope_block_cover_r27.py` enumerates every oriented
projective cut, every selector, and every partition of the indicated type.
It compares partitions by the exact integer numerator of

```math
Q(A[S])-\max_{d\in\mathcal C}X_d(S),
```

and tests the remaining radical inequality by exact squaring.  It also checks
(R27.12)--(R27.13) on random instances.

- `A_6`: among all 15 pair partitions, the partition
  `(01)(23)(45)` gives an eight-word coset.  Every word has `R_2=30`, and the
  coset covers with `t=0` for `m=3,4,5`.

- `A_8`: among all 105 pair partitions, 16-word cosets at cap `R_2<=64`
  cover with `t=0` for `m=4,6,7` (the optimizing partition may depend on
  `m`).

- `A_9`: among all 945 `2+2+2+2+1` partitions, subsets of a 32-word coset at
  cap `R_2<=80` cover with `t=0` for `m=5,7,8`; optimizing codebooks retain
  respectively 26, 30, and 20 eligible words.

- A finite compression wall occurs for the more compressed `A_9` triple
  class.  Exhausting all 280 `3+3+3` partitions at `m=5`, even with cap
  `R_2<=112`, the smallest possible worst exact base gap is
  `736/72=92/9`.  After subtracting
  `p^(3/2)q_9=40sqrt(5)/9`, the necessary tolerance is at least

  ```math
  \boxed{\frac{92-40\sqrt5}{9}=0.284142322223\ldots>0.}
  \tag{R27.14}
  ```

  Thus no eligible eight-word triple-block coset (or its low-row subset)
  proves the zero-tolerance envelope at this finite point.  This is scoped to
  `A_9,m=5,t=0` and this fixed block-coset architecture; it is not an
  asymptotic falsifier and is harmless for a tolerance growing like
  `n^(3/2-c)`.

The output is saved in `tmp/envelope_block_cover_r27.out` and ends with
`PASS`.

## 6. Frontier

What is proved is a nontrivial metric-entropy construction: an entire
`exp(O(n^(3/4-c)))` block-coset family is simultaneously low-row at the
required scale, and a genuine greedy theorem turns a per-selector coset hit
probability into a full envelope codebook.  Finite exact minimizers support
pair-block alignment but already reject overcompression to three blocks at
zero tolerance.

The open step is exactly (R27.9), or its fractional form (R27.10), for an
`A`-adapted balanced partition.  Neither scalar tails of a preselected
`X_d`, random signs without the quotient maximum, nor the tautological full
cut set establishes that inequality.
