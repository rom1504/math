# The prescribed dense Hadamard stabilization has an exact full-sign half floor

Status: complete elementary proof; independently audited PASS by
`transfer_seeds` (including H144, same-spin conversion, stabilization,
rectangular witness, dimension factors, and Walsh-only scope).
This closes the weaker full-sign landing statement left open in Section4
of `fresh_fixed_walsh_factorization_final_audit_2026_09_05.md`, WITHOUT
proving the stronger identity R=T. It does not prove original convergence.

## 1. Precisely the prescribed outer family

Let H4=J4-2I4. Fix the symmetric regular Hadamard H144 constructed in
Section1.1 of `fresh_limit_algebra_2026_09_05.md`. For clarity, it uses an
order12 Hadamard F, the symmetric matrix
`K_((i,j),(p,q))=F_iq F_pj`, and diagonal `D_((i,j),(i,j))=F_ij`;
then H144=DKD. Directly `K^2=144I`, and K vec(F)=12 vec(F), so H144
is symmetric, regular with row sum12, and has sign entries.

For `s=4^a144^b`, a,b>=0, use only

```math
H_s=H_4^{\otimes a}\otimes H_{144}^{\otimes b},\qquad
H_s^2=sI,\qquad H_s\mathbf1=\sqrt s\mathbf1.
```

For a real symmetric matrix B of order n, INCLUDING its diagonal, define

```math
q(B)=\tfrac12\max_{z\in\{\pm1\}^n}|z^TBz|,
\quad\beta(B)=\max_{x,y\in\{\pm1\}^n}|x^TBy|,
\quad R(B)=\sup_{s=4^a144^b}{q(H_s\otimes B)\over s^{3/2}}.
```

No arbitrary orthogonal basis or additional Hadamard matrix is allowed.
The tensor factors commute up to a simultaneous coordinate permutation,
which does not alter q or beta.

## 2. Exact same-spin conversion, not ordinary polarization

For any symmetric C and Boolean x,y, take the Boolean vector
`z=(x,x,y,y)` in four blocks. Because H4=J4-2I4,

```math
\tfrac12 z^T(H_4\otimes C)z=4x^TCy.
```

Thus

```math
q(H_4\otimes C)\ge4\beta(C).                         (1)
```

This remains true when C has nonzero diagonal. Conversely
`q(C)<=beta(C)/2`. Applying these two inequalities at each allowed outer
gives the EXACT identity

```math
R(B)=\tfrac12\sup_s {\beta(H_s\otimes B)\over s^{3/2}}. (2)
```

For the lower direction use outer4s in (1), whose normalization is
`(4s)^(3/2)=8s^(3/2)`. Hence there is no hidden factor-two loss in converting
the bilinear witness below to the actual same-spin objective.

Regularity gives monotonicity in each outer exponent: tensor any witness
with an all-one vector in a further factor. Therefore the supremum is
also the cofinal limit, and for every allowed t,

```math
R(H_t\otimes B)=t^{3/2}R(B).                           (3)
```

Indeed the left side divided by t^(3/2) is the supremum over a shifted
cofinal quadrant of the two nondecreasing outer exponents. Finiteness follows
from `R(B)<=n ||B||op/2`.

## 3. Rectangular Hadamard-column witness

Let s be an allowed outer order with s>=n. Let E be the s-by-n matrix
consisting of the first n columns of I_s, and let

```math
M=H_s^T E.
```

EVERY entry of M is a sign. Under the standard column-stacking convention,

```math
(B\otimes H_s)\operatorname{vec}(M)
=\operatorname{vec}(H_s M B^T)
=s\operatorname{vec}(E B^T).                          (4)
```

Taking the coordinatewise sign of the right-hand side as the other
Boolean vector (arbitrary signs at zeros), one obtains

```math
\beta(H_s\otimes B)\ge s\sum_{i,j}|B_{ij}|,
\qquad R(B)\ge{\sum_{i,j}|B_{ij}|\over2\sqrt s}.         (5)
```

Equation (4) uses only orthogonality of the ACTUAL prescribed H_s and its
flat entries. It does not replace a polar matrix by an arbitrary orthogonal
catalyst. For s=n it is simply
`(B tensor H_n) vec(H_n^T)=n vec(B^T)`.

## 4. Dense prescribed orders remove the dimension mismatch

The allowed s=4^a144^b are multiplicatively relatively dense: successive
ratios tend to1. This follows by the elementary irrational-rotation argument
from `log144/log4` irrational, as already proved in the outer construction.
Let s_+(m) be the smallest allowed order at least m. Then

```math
s_+(m)/m\longrightarrow1.                              (6)
```

Apply (5) not only to B but to `C_t=H_t tensor B` for any allowed t tending
to infinity. Its order is nt, and its entrywise absolute sum is
`t^2 sum_ij |B_ij|`. Choose the allowed outer s=s_+(nt), and use (3):

```math
t^{3/2}R(B)=R(C_t)
\ge {t^2\sum_{i,j}|B_{ij}|\over2\sqrt{s_+(nt)}}.
```

After division and passage to the limit using (6), this proves the exact
finite-seed inequality

```math
\boxed{R(B)\ge{1\over2\sqrt n}\sum_{i,j}|B_{ij}|.}       (7)
```

For every symmetric FULL sign seed, (7) is the desired exact half floor:

```math
\boxed{R(B)\ge\tfrac12 n^{3/2}.}                       (8)
```

For every hollow sign seed,

```math
R(A)\ge\tfrac12(n-1)\sqrt n.                          (9)
```

There is no claim R=T for general real seeds. Equation (7) is weaker than
that equality and is obtained by an explicit Boolean witness and exact
same-spin conversion.

## 5. Walsh-only specialization and its exact scope

If the outer family is restricted to H4 powers ONLY, (1)--(5) and (3)
remain valid. At seed orders n=4^r, take s=n directly in (5). Hence the
same full-sign floor (8), and hollow floor (9), hold at those seed orders
without using H144 or any density argument.

For general fixed n, the two-generator density argument in Section4 is
not available for Walsh-only outers; (7) is NOT asserted for that smaller
family at arbitrary n on the basis of this proof. The power-four seed
subsequence already suffices to obstruct universal lossless Walsh transfer
on actual asymptotically minimizing seeds.

## 6. Original-signing consequence and hollowing

For each order n choose an actual hollow minimizer A_n and any sign diagonal
D_n. Then B_n=A_n+D_n is a full symmetric sign seed and

```math
q(B_n)\le M_n+n/2,\qquad R(B_n)\ge n^{3/2}/2.
```

The banked original all-order upper bound `limsup m_n<=c_U<1/2` therefore
gives an actual-minimizer stabilization loss at least
`(1/2-c_U-o(1)) n^(3/2)`. This is not just a PSD-majorant certificate gap.
By monotonicity, for every fixed seed B_n, all sufficiently far cofinal
prescribed outer exponents have normalized SAME-SPIN cap arbitrarily close
to R(B_n), hence at least the half floor minus any chosen tolerance.

The matrix `H_s tensor B_n` is a full actual sign matrix of order N=sn.
Deleting its diagonal changes the quadratic energy by at most N/2. Its
normalized original cap therefore has liminf at least1/2 along such
cofinal outer orders, for every fixed seed. There are no missing off-diagonal
fibres when the seed is full.

For growing seeds, choose the outer exponents after the seed and the desired
tolerance; no uniform finite-depth convergence is needed. On n=4^r the
Walsh-only witness already occurs at outer4n, and all further regular
Walsh amplifications retain its normalized value; hollowing again costs
only O(N).

Consequently the minimum asymptotic stabilized full-sign coefficient equals
1/2: the lower bound is (8), and regular symmetric Hadamard seeds attain it
by the spectral upper bound. No full sign seed can satisfy the previously
proposed strict certificate `R(B)<n^(3/2)/2` in this prescribed two-generator
family. This rules out that seed-transfer architecture; it neither proves
nor disproves convergence of the original M_n/n^(3/2).

### Exact finite Walsh witness with no diagonal error

There is an especially clean choice when n=4^r, r>=1. Complete ANY hollow
sign seed A by a sign diagonal D having trace zero. Then q(A+D)=Q(A)
exactly. Set B=A+D, use outer H_n in the column witness, and the additional
H4 in the same-spin conversion. The resulting full sign matrix
K=H_(4n) tensor B has order N=4n^2 and admits a Boolean vector z with

```math
\tfrac12 z^TKz=4n^3=\tfrac12N^{3/2}.
```

Its trace is tr(H_(4n)) tr(B)=0. Therefore deleting its diagonal changes
NO same-spin quadratic energy, and the resulting actual hollow signing
satisfies Q>=N^(3/2)/2 EXACTLY. All further regular Walsh amplifications
retain this exact normalized witness, and retain zero trace. Thus the
Walsh-only obstruction on power-four actual minimizing seeds has neither
an asymptotic diagonal error nor an unspecified stabilization depth.
