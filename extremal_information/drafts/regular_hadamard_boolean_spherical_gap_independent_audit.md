# Independent audit: regular-Hadamard Boolean--spherical gap

**Audited files.** `regular_hadamard_boolean_spherical_gap.md` and
`verify_regular_hadamard_boolean_spherical_gap.py`.

## Verdict

**PASS WITH SCOPE REPAIRS.**  BG.6--BG.20 are exact, with the stated
normalizations.  The tensor family is valid, and exact-sign completion
changes each optimum by at most `Q(C)`, so it cannot erase the leading gap.
The result decisively falsifies uniform rounding of the SA.3 spherical trust
value to the Boolean old-spin value at bounded total port mass.

Three claims should be worded more narrowly:

1. In the graph model `m` is an integer number of repeated spins.  “For
   every `m>=0`” is exact for nonnegative integer `m`, or for real `m` only
   after reinterpreting it as a field weight.
2. BG.1's displayed algebraic hypotheses do not imply that `A` is a signing.
   The “exact-sign completion” language additionally needs `H` to have
   off-diagonal entries in `{+-1}` (and the Hadamard normalization
   `r=sqrt(n)`).  The explicit tensor family has these properties.
3. A nonzero spherical integrality gap does not prove that `(G,R)` lacks
   enough information to determine or approximate the Boolean response by
   some different functional.  It proves that the **spherical trust
   functional itself**, and any rounding theorem seeking to attain it, fail
   uniformly.  A true insufficiency theorem for the Gram state would require
   two states with the same `(G,R)` and separated Boolean responses.

## 1. Exact Boolean formula

For a Boolean old spin `u`,

```math
{1\over2}u^TAu
={1\over2}u^THu-{1\over2}\sum_iH_{ii}u_i^2
={1\over2}u^THu                                   \tag{ABG.1}
```

because `u_i^2=1` and `tr H=0`.  This confirms the normalization in BG.4.
For fixed `u`, the two uncompleted auxiliary shores can independently align
their spins with both the outer sign and their incident field.  Optimizing
them and the outer absolute channel therefore gives BG.3 exactly.

If `a,b` are orthogonal Boolean words, their agreement and disagreement
sets each have size `n/2`.  Hence for all endpoint signs

```math
\|\epsilon_1a+\epsilon_2b\|_1=n.                    \tag{ABG.2}
```

It follows that

```math
|a^Tu|+|b^Tu|\le n                                  \tag{ABG.3}
```

on the Boolean cube.  The spectral bound gives
`|u^THu|/2<=rn/2`, so BG.6 is an upper bound.  At `u=a`, both upper bounds
are attained simultaneously: the quadratic is `rn/2`, while the two fields
are `n` and zero.  Thus

```math
\mathcal B_m={rn\over2}+mn                           \tag{ABG.4}
```

with no hidden incompatibility between the maximizers.

The statement is algebraically valid for every real field coefficient
`m>=0`; the literal repeated-vertex construction requires integer `m`.

## 2. Exact spherical formula and SA.3/GE normalization

On the sphere `||u||^2=n`, duality and orthogonality give

```math
|a^Tu|+|b^Tu|
\le\max_{\epsilon_1,\epsilon_2}
 \|\epsilon_1a+\epsilon_2b\|_2\sqrt n
=\sqrt2,n.                                          \tag{ABG.5}
```

Together with the spectral quadratic bound this proves the upper half of
BG.7.  The vector `(a+b)/sqrt(2)` has norm `sqrt(n)`, lies entirely in the
`+r` eigenspace, and pairs with each port by `n/sqrt(2)`.  It attains both
bounds, proving

```math
\mathcal S_m={rn\over2}+\sqrt2,mn.                  \tag{ABG.6}
```

The gap `(sqrt(2)-1)mn` is therefore exact.

The same constants follow directly from SA.3 and GE.  Here

```math
G=R=I_2,
\qquad g=h=2                                         \tag{ABG.7}
```

in the positive outer channel.  With `mu=m/r`, GE.13 has `a=4,b=0`, so
the trust formula gives

```math
{\mathcal S_m\over rn}={1\over2}+\sqrt2\,{m\over r}. \tag{ABG.8}
```

The Boolean value is `1/2+m/r`.  Since the total mass is
`c=2m/r`, the normalized gap is

```math
{c(\sqrt2-1)\over2}.                                 \tag{ABG.9}
```

At `m=r/2`, this is BG.15.  The corresponding trust minimizer is
`t=1/sqrt(2)`, bounded away from the hard boundary.  Thus the statement that
a trust margin does not cure this particular integrality gap is justified.

The spherical relaxation here is representation-dependent: it relaxes
`u^THu/2`, whose diagonal is harmless on the Boolean cube because of
ABG.1.  It is not the literal spherical relaxation of the hollow matrix
`A`, for which the deleted diagonal need not remain constant off the cube.
The draft explicitly records this important scope.

## 3. Tensor family and total-order scale

For the order-16 base matrix, `H_16^2=16I`, `tr H_16=0`, and both
`1` and `v_0` are orthogonal Boolean `+4` eigenvectors.  Tensor powers give

```math
H_j^2=16^jI=r_j^2I,
\qquad \operatorname{tr}H_j=(\operatorname{tr}H_{16})^j=0. \tag{ABG.10}
```

The words in BG.18 are Boolean and satisfy

```math
H_ja_j=r_ja_j,
\qquad H_jb_j=r_jb_j,
\qquad a_j^Tb_j=0.                                   \tag{ABG.11}
```

Here `r_j=4^j=sqrt(n_j)` and `m_j=r_j/2` is integral.  Consequently

```math
{(\sqrt2-1)m_jn_j\over n_j^{3/2}}
={\sqrt2-1\over2}.                                   \tag{ABG.12}
```

The completed parent order before auxiliary interactions is
`N_j=n_j+r_j`, so

```math
{n_j^{3/2}\over N_j^{3/2}}
=\left(1+{1\over\sqrt{n_j}}\right)^{-3/2}\longrightarrow1. \tag{ABG.13}
```

BG.19--BG.20 are therefore correctly normalized.  The canonical verifier
checks the tensor identities through order 4096; the independent verifier
checks an additional order-256 construction without relying on the closed
form alone.

For a general abstract matrix satisfying only BG.1--BG.2, `r/2` need not be
an integer and `A` need not have sign entries.  Exact graph-family language
should be attached to the regular-Hadamard tensor specialization, or the
general setup should explicitly assume a symmetric sign Hadamard matrix.

## 4. Exact-sign completion Lipschitzness

Let `E_0(u,y)` be the uncompleted energy, where the auxiliary vector `y`
remains Boolean, and let `E_C(y)` be the hollow auxiliary completion.  Its
uniform norm is

```math
\sup_y|E_C(y)|=Q(C).                                 \tag{ABG.14}
```

For any common feasible domain, including Boolean old spins or spherical
old spins with Boolean auxiliaries,

```math
\left|\sup|E_0+E_C|-\sup|E_0|\right|\le Q(C).        \tag{ABG.15}
```

This proves both inequalities in BG.21.  Combining them gives the more
explicit completed-gap bound

```math
S_C-B_C
\ge(\mathcal S_m-\mathcal B_m)-2Q(C).               \tag{ABG.16}
```

Since a hollow signing on `2m` vertices has
`Q(C)<=binom(2m,2)`, at `m=r/2` the loss is `O(r^2)=O(n)`, whereas the gap
is `Theta(rn)=Theta(n^(3/2))`.  BG.23's `-O(n)` notation correctly absorbs
the factor two.

The independent verifier exhausts the order-16 old cube and all 16
auxiliary words for two different complete signings on four auxiliary
vertices.  It computes both the completed Boolean optimum and the exact
old-spin spherical optimum and confirms ABG.15--ABG.16.  Thus no assumption
about a shared optimizer is hidden in the Lipschitz argument.

This exact-sign claim requires the old hollow matrix `A` itself to be a sign
matrix.  BG.1 alone does not ensure that; the explicit regular-Walsh family
does.

## 5. What is and is not falsified

The example proves all of the following:

- bounded `mp/r` does not imply an `o(rn)` Boolean--spherical gap;
- low spectral rank of the exposed field does not make the sphere optimizer
  Boolean-close;
- an interior trust minimizer does not eliminate this `l_1/l_2` gap;
- arbitrary exact-sign auxiliary completion cannot remove it at leading
  scale.

It therefore falsifies a uniform theorem asserting that the SA.3 spherical
trust response approximates, or can always be rounded to, the Boolean
response under the sole assumption `mp/r=O(1)`.

It does **not** by itself prove that the Gram pair `(G,R)` is an insufficient
statistic for the Boolean response.  In this family the common state
`G=R=I_2` could in principle be assigned the corrected Boolean value by a
different response functional.  To prove information insufficiency one
would need a collision: equal Gram data with macroscopically different
Boolean responses.  Accordingly the claims that a Boolean-net property is
“strictly extra information beyond `(G,R)`” and that the PSD/Gram carrier
cannot be promoted by *any* functional are stronger than BG.1 proves.  They
should be narrowed to rounding/recovery of the spherical trust value unless
such a collision is supplied from another theorem.

## 6. Verifier assessment

The canonical verifier passes 59 checks.  It exhausts the base Boolean old
spins for `m=r/2` and `m=r`, and verifies the tensor eigenvectors,
orthogonality, support norms, and normalized formulas through order 4096.
It evaluates the spherical side from the proved closed formula rather than
an independent optimizer, which is adequate for the algebra but leaves the
completion claim untested.

The independent verifier adds:

- five direct SA.3/GE normalization checks;
- two exact complete-signing Lipschitz tests;
- an independent order-256 tensor check.

Run:

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_regular_hadamard_boolean_spherical_gap.py

./.venv/bin/python \
  extremal_information/experiments/verify_regular_hadamard_boolean_spherical_gap_independent_audit.py
```

Both pass.
