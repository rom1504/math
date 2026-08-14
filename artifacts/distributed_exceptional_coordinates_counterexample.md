# Distributed exceptional coordinates under conference Frobenius control

Status: **proved scalable counterexample to `o(n)`-vertex puncturing**.  A
bounded-operator-norm signing can have the full averaged
diagonal-monomial/conference structure and conference limiting pressure,
while one fixed power has large off-diagonal entries on a linear matching.
Consequently those violations cannot, in general, be covered by `o(n)`
vertices.

## 1. Statement

Fix

```math
0<\beta<{1\over2},
\qquad 0<\eta<{1\over2}.
\tag{1.1}
```

Along any infinite sequence of symmetric conference signings `C_n`, there
are symmetric zero-diagonal signings `A_n` such that, with

```math
X_n={\beta A_n\over\sqrt n},
\qquad c_n^2=\beta^2(1-1/n),
\tag{1.2}
```

the following all hold:

```math
\|X_n\|_{\rm op}\le\beta+o(1)<{1\over2},
\tag{1.3}
```

```math
{1\over\sqrt n}\|X_n^2-c_n^2I\|_{\rm F}=o(1),
\tag{1.4}
```

and the averaged diagonal-monomial distribution of `X_n` is the same as
the conference/Bernoulli distribution.  Nevertheless, there are
`floor(n/10)` vertex-disjoint pairs `(i,j_i)` for which

```math
\boxed{
 |(X_n^2)_{i,j_i}|
 =2\beta^2 n^{-1/2+\eta}(1+o(1)).}
\tag{1.5}
```

Hence, for every fixed `0<epsilon<eta`, the graph of entries satisfying

```math
|(X_n^2)_{ij}|\ge n^{-1/2+\varepsilon}
\tag{1.6}
```

has vertex-cover number at least `floor(n/10)` for all sufficiently large
`n`.  No `o(n)` set of exceptional coordinates covers the violations.

The example is thermodynamically invisible at leading order:

```math
\log\overline Z(X_n)
=\log\overline Z\left({\beta C_n\over\sqrt n}\right)+o(n)
=n\psi(\beta)+o(n).
\tag{1.7}
```

Thus even conference-level pressure does not force an exceptional set of
sublinear cardinality.

## 2. Construction of a linear exceptional matching

Let `m=floor(n/10)`.  Choose pairwise disjoint sets `I,J` of size `m`, a
bijection `i -> j_i` from `I` to `J`, and put

```math
L=[n]\setminus(I\cup J).
\tag{2.1}
```

For each `i in I`, define the eligible leaf set

```math
E_i=\{k\in L:(C_n)_{ik}(C_n)_{kj_i}=-1\}.
\tag{2.2}
```

Because `C_n^2=(n-1)I`, for distinct `i,j_i` exactly `(n-2)/2` of the
indices other than the endpoints have product `-1`.  Removing `I\cup J`
therefore leaves

```math
|E_i|\ge(3/10)n-O(1).
\tag{2.3}
```

Put `r=floor(n^(1/2+eta))`.  Independently for every `i in I`, choose a
uniform `r`-subset `S_i` of `E_i`, and flip the conference edges `{i,k}`
for `k in S_i`.  Denote the perturbation by

```math
D_n=A_n-C_n.
\tag{2.4}
```

It is supported only between `I` and `L`, with

```math
(D_n)_{ik}=-2(C_n)_{ik}
\quad(k\in S_i),
\qquad
\|D_n\|_{\rm F}^2=8mr.
\tag{2.5}
```

For a matched pair `(i,j_i)`, columns indexed by `J` vanish identically in
`D_n`.  Consequently

```math
(C_nD_n)_{i,j_i}=(D_n^2)_{i,j_i}=0,
\tag{2.6}
```

while every selected leaf contributes exactly two to the other cross term:

```math
(D_nC_n)_{i,j_i}
=\sum_{k\in S_i}-2(C_n)_{ik}(C_n)_{kj_i}=2r.
\tag{2.7}
```

Since `(C_n^2)_{i,j_i}=0`, equations (1.2) and (2.6)--(2.7) prove (1.5).
The pairs form a matching of size `m`, proving the vertex-cover assertion.

## 3. The perturbation has `o(sqrt(n))` operator norm

Write `H` for the `I` by `L` block of `D_n`; then
`||D_n||_op=||H||_op`.  Let `N_i=|E_i|`.  On eligible entries,

```math
-2(C_n)_{ik}{\bf1}_{\{(C_n)_{ik}(C_n)_{kj_i}=-1\}}
=(C_n)_{kj_i}-(C_n)_{ik}.
\tag{3.1}
```

It follows that

```math
\mathbb EH
=R\{(C_n)_{J,L}-(C_n)_{I,L}\},
\qquad
R=\operatorname{diag}(r/N_i).
\tag{3.2}
```

Every submatrix of `C_n` has operator norm at most `sqrt(n-1)`, so (2.3)
gives

```math
\|\mathbb EH\|_{\rm op}=O(r/\sqrt n).
\tag{3.3}
```

The centered rows of `H` are independent.  A centered row has Euclidean
norm at most `4 sqrt(r)`.  Sampling without replacement from `N_i`
eligible coordinates gives row-covariance norm at most

```math
8r/N_i=O(r/n).
\tag{3.4}
```

Indeed, before multiplication by the fixed coordinate signs, its covariance
on the eligible set is

```math
4p_i(1-p_i)\left{
 {N_i\over N_i-1}I-{1\over N_i-1}ss^{\mathsf T}
 \right},
\qquad p_i=r/N_i,
\tag{3.5}
```

where `s` is the vector of coordinate signs.  Thus both rectangular matrix-
Bernstein variance parameters are `O(r)`: the row-side one is bounded by
the maximum centered row second moment, and the column-side one by summing
(3.4) over `m=O(n)` rows.  Rectangular matrix Bernstein now gives, with
probability tending to one,

```math
\|H-\mathbb EH\|_{\rm op}=O(\sqrt r\log n).
\tag{3.6}
```

In particular, the probabilistic method supplies deterministic choices of
the sets `S_i` for which

```math
\boxed{
\|D_n\|_{\rm op}
=O(r/\sqrt n+\sqrt r\log n)=o(\sqrt n).}
\tag{3.7}
```

This proves (1.3).

## 4. Frobenius and averaged diagonal-monomial structure

Using `C_n^2=(n-1)I`,

```math
X_n^2-c_n^2I
={\beta^2\over n}(C_nD_n+D_nC_n+D_n^2).
\tag{4.1}
```

Moreover,

```math
\|C_nD_n\|_{\rm F}=\sqrt{n-1}\,\|D_n\|_{\rm F},
\qquad
\|D_n^2\|_{\rm F}
\le\|D_n\|_{\rm op}\|D_n\|_{\rm F}.
\tag{4.2}
```

Combining (2.5), (3.7), and (4.1)--(4.2) yields

```math
{1\over\sqrt n}\|X_n^2-c_n^2I\|_{\rm F}
=O\left(\sqrt{r/n}+{r\log n\over n}+(r/n)^{3/2}\right)=o(1),
\tag{4.3}
```

because `r=n^(1/2+eta)(1+o(1))` and `eta<1/2`.

For completeness, (4.3), `Delta(X_n)=0`, and the uniform operator-norm
bound imply convergence of every averaged diagonal monomial.  In the
normalized Frobenius seminorm `||M||_(2,n)=n^(-1/2)||M||_F`, multiplication
by a uniformly bounded matrix is continuous and the diagonal projection
`Delta` is a contraction.  Induction over products and nested applications
of `Delta` therefore permits replacing `X_n^2` by `c_n^2I` and
`Delta(X_n)` by zero in every fixed diagonal monomial, with `o(1)` error.
This is precisely the conference/Bernoulli averaged diagonal distribution.

## 5. Pressure and the consequence for regularization

For symmetric matrices `X,Y`, the function `log cosh` is one-Lipschitz and
`||x||_2^2=n` on the Boolean cube.  Termwise comparison of the partition
functions gives the exact bound

```math
|\log\overline Z(X)-\log\overline Z(Y)|
\le {n\over2}\|X-Y\|_{\rm op}.
\tag{5.1}
```

Taking `X=beta A_n/sqrt(n)`, `Y=beta C_n/sqrt(n)`, and using (3.7) proves
the first equality in (1.7).  The second is the conference
high-temperature free-energy limit.

This counterexample separates two notions that an exceptional-coordinate
argument would otherwise conflate:

1. **Sublinear endpoint cover is false.**  Deleting endpoints so that every
   entry in (1.6) disappears requires at least `n/10+o(n)` vertices.
2. **Pressure-cheap global editing remains possible.**  Reversing all of the
   distributed flips touches a linear set of center vertices and
   `Theta(nr)` edges, yet costs only `o(n)` in log pressure by (5.1).

Thus bounded norm, Frobenius-near conference structure, averaged diagonal-
monomial convergence, and even conference-level asymptotic pressure do
not imply an `o(n)` exceptional-coordinate cover.  Any viable regularization
theorem must allow a distributed operator-small edge edit or use a
universality theorem based on an averaged norm; puncturing `o(n)` bad
vertices cannot follow from the stated hypotheses.
