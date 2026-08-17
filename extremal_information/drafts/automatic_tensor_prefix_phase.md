# Automatic phase laws for regular-Hadamard prefix hierarchies

Status: rigorous theorem draft.  This abstracts the order-four Walsh prefix
example.  It concerns one explicit coherent family of dense signings, not the
minimizing values in the motivating problem.

## Theorem ATP.1 (continuous tensor-prefix phase law)

Let `H` be a symmetric Hadamard matrix of order `h>1`, normalized so that its
top-left entry is `1`.  Assume that there is a Boolean vector `u` with

```math
Hu=\sqrt h\,u.                                             \tag{ATP.1}
```

Thus `H` is regular in the only sense needed below.  Put `H_r=H^(tensor r)`
in the lexicographic ordering `H_(r+1)=H tensor H_r`.
Because the top-left entry of `H` is one, the matrices `H_r` are compatible
leading principal blocks of one infinite symmetric sign matrix `S`.  Let
`A_n` be the hollow leading `n by n` block of `S`, and, for `R_r=h^r`, set

```math
F_r(t)={Q(A_(floor(tR_r)))\over R_r^(3/2)},
\qquad 1\le t\le h.                                      \tag{ATP.2}
```

Then `F_r` converges uniformly on `[1,h]` to a continuous nondecreasing
function `F`.  If

```math
G_r(t)={Q(A_(floor(th^r)))\over floor(th^r)^(3/2)},
```

then `G_r` converges uniformly on `[1,h]` to
`L(t)=F(t)/t^(3/2)`.  In particular,

```math
{Q(A_(floor(th^r)))\over floor(th^r)^(3/2)}
\longrightarrow L(t):={F(t)\over t^(3/2)}                 \tag{ATP.3}
```

Equivalently, if

```math
r(n)=floor(log_h n),\qquad t_n={n\over h^(r(n))}\in[1,h),
```

then

```math
{Q(A_n)\over n^(3/2)}-L(t_n)\longrightarrow0.             \tag{ATP.4}
```

In particular this all-order normalized sequence converges if and only if
the continuous phase profile `L` is constant.

### Proof

Fix a base-`h` rational `t=p/h^k` in `[1,h]`, where
`h^k<=p<=h^(k+1)` is an integer.  For `r>=k`, the leading
`p h^(r-k)` block of

```math
H_(r+1)=H_(k+1) tensor H_(r-k)
```

is

```math
B_(p,k) tensor H_(r-k),                               \tag{ATP.5}
```

where `B_(p,k)` is the leading `p by p` block of `H_(k+1)`.  The outer
template is fixed.  The Boolean lift `x mapsto x tensor u` preserves the
normalized full quadratic response, while the diagonal deletion has size
`O((p h^(r-k))^(-1/2))`.  The regular-Hadamard amplification theorem
therefore makes `F_r(t)` converge at every base-`h` rational `t`.

It remains to promote convergence on this dense set.  For `1<=t<=s<=h`,
write

```math
n=floor(tR_r),\qquad m=floor(sR_r),\qquad d=m-n
```

and decompose the hollow prefix at order `m` as

```math
A_m=\begin{pmatrix}A_n&C\\C^T&D\end{pmatrix}.          \tag{ATP.6}
```

Both `C` and the unhollowed matrix underlying `D` are coordinate
compressions of `H_(r+1)`, so

```math
||C||_(2->2)\le\sqrt{hR_r},\qquad
||D||_(2->2)\le\sqrt{hR_r}+1.                         \tag{ATP.7}
```

For Boolean block vectors `x,y`, the newly exposed half-quadratic energy is
bounded by

```math
|x^TCy+{1\over2}y^TDy|
\le \sqrt{hR_rnd}+{d\over2}(\sqrt{hR_r}+1).            \tag{ATP.8}
```

Principal deletion gives `Q(A_n)<=Q(A_m)`: extend a maximizing spin on the
first block by independent unbiased missing spins, average the signed
quadratic, and select an extension with the required sign.  Dividing
(ATP.8) by `R_r^(3/2)` yields

```math
0\le F_r(s)-F_r(t)
\le \sqrt h\sqrt{(n/R_r)(d/R_r)}
   +{\sqrt h\over2}{d\over R_r}
   +{d\over2R_r^(3/2)}.                               \tag{ATP.9}
```

Since `n/R_r<=h` and `d/R_r<=s-t+R_r^(-1)`, this is a common modulus
`O_h(sqrt(s-t)+s-t)+o_r(1)`.  The functions are uniformly bounded and
asymptotically equicontinuous.  Convergence on a finite base-`h` rational
net therefore makes them uniformly Cauchy.  Their uniform limit `F` is
continuous and nondecreasing.  Equations (ATP.3)--(ATP.4) and the uniform
convergence of `G_r` now follow by division by the uniformly positive factor
`(floor(th^r)/h^r)^(3/2)`, which converges uniformly to `t^(3/2)`.  `square`

## Corollary ATP.2 (endpoint normalization)

At the geometric orders `h^r`, the spectral upper bound and the Boolean
eigenvector in (ATP.1) give

```math
{Q(A_(h^r))\over(h^r)^(3/2)}={1\over2}+O(h^(-r/2)).   \tag{ATP.10}
```

The same holds at the right endpoint because `h h^r=h^(r+1)`.  Hence
`L(1)=L(h)=1/2`.  A single fixed outer prefix template with limiting ratio
strictly above `1/2` proves that `L` is nonconstant and gives genuine
all-order nonconvergence for the coherent prefix sequence.

## Scope and lesson

The theorem identifies the exact residual all-order state of a broad tensor
prefix hierarchy: a continuous mantissa/scale phase.  Exact recovery at
every geometric refinement and convergence at every fixed outer template do
not synchronize those phases.  The compact phase is much smaller than the
Boolean landscape, but it is not automatically a singleton.

This is a strict cross-scale mechanism, not a theorem about `M_n`: nothing
here says that the explicit prefixes are near-minimizers, or that minimization
over signings preserves one coherent tensor hierarchy.
