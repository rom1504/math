# Wave 14 memo: balanced endpoint residual is a cross-Gram entry

Status: all identities, linear-algebra statements, and finite computations in
this memo are **Verified** by direct derivation and by
`tmp/check_endpoint_crossgram_r14.py`.  The final temporal interpretation is
only a possible successor target.  No asymptotic comparison is proved.

## 1. Exact endpoint identity

Let `A` be a symmetric zero-diagonal sign matrix.  Write

```math
P=P(A),\qquad N=N(A),\qquad Q=\max\{P,N\},\qquad I=|P-N|.
```

Choose a positive ground `p` and a negative ground `n`, so

```math
p^{\mathsf T}Ap=P,\qquad n^{\mathsf T}An=-N.
```

Switch by `p`.  Thus `p=\mathbf 1`, and write the negative endpoint as
`n=\mathbf 1_S\oplus(-\mathbf 1_T)`.  Put

```math
h_S=\mathbf1_S^{\mathsf T}A[S]\mathbf1_S,
\qquad
h_T=\mathbf1_T^{\mathsf T}A[T]\mathbf1_T,
\qquad
b=p^{\mathsf T}An.
```

If `c=\mathbf1_S^{\mathsf T}A[S,T]\mathbf1_T`, direct expansion gives

```math
P=h_S+h_T+2c,\qquad -N=h_S+h_T-2c.
```

Consequently, with `H=(P-N)/2`,

```math
\boxed{h_S+h_T=H,\qquad h_S-h_T=b.}
\tag{R14.1}
```

Equation (10.534) says that the decrement-tolled residual in shore `X`
and orientation `\sigma` is

```math
r_X^\sigma=[\sigma h_X-I/2]_+.
```

The total of the four root buckets is therefore

```math
\mathfrak r(p,n)
=[|h_S|-|H|]_+ + [|h_T|-|H|]_+.
\tag{R14.2}
```

Let `a=|H|` and `t=|b|`.  From (R14.1), the unordered pair
`\{|H+b|,|H-b|\}` is `\{a+t,|a-t|\}`.  Hence

```math
\begin{aligned}
\mathfrak r(p,n)
&=\frac12\bigl([|H+b|-2a]_+ + [|H-b|-2a]_+\bigr)\\
&=\frac12\bigl([t-a]_+ + [t-3a]_+\bigr).
\end{aligned}
```

Thus the general exact threshold formula is

```math
\boxed{
\mathfrak r(p,n)
=\frac12\left(
[|p^{\mathsf T}An|-I/2]_+
+[|p^{\mathsf T}An|-3I/2]_+
\right).
}
\tag{R14.3}
```

This formula also audits all threshold constants in (10.534).  In
particular, the conclusion requested here is specifically a **balanced-parent
statement**:

```math
\boxed{
P(A)=N(A)=Q(A)
\quad\Longrightarrow\quad
\mathfrak r(p,n)=|p^{\mathsf T}An|.
}
\tag{R14.4}
```

Because `h_S` is twice an integer edge sum, balance gives
`b=2h_S\in4\mathbb Z`.  Every nonzero balanced root residual is therefore at
least four.  Away from balance, zero residual means only
`|p^{\mathsf T}An|\le I/2`; it does **not** imply cross annihilation.

## 2. Cross-Gram neutrality and the sharp dimension theorem

Let

```math
L_+=\operatorname{span}\mathcal G_+(A),\qquad
L_-=\operatorname{span}\mathcal G_-(A),
```

and let `d_\pm=\dim L_\pm`.  For any spanning column matrices `G_+` and
`G_-`, define the endpoint cross-Gram matrix

```math
C_{+-}=G_+^{\mathsf T}AG_-.
```

At balance, (R14.4) and bilinearity give the exact equivalence

```math
\boxed{
\text{every endpoint pair is residual-neutral}
\iff C_{+-}=0
\iff L_+^{\mathsf T}AL_-=0.
}
\tag{R14.5}
```

There is a sharp rank version.  Put

```math
K=\ker A,\quad k_\pm=\dim(L_\pm\cap K),\quad r=\operatorname{rank}A,
```

and let `\rho` be the rank of the restricted pairing
`(u,v)\mapsto u^{\mathsf T}Av` on `L_+\times L_-`; equivalently,
`\rho=\operatorname{rank}C_{+-}` for any spanning ground matrices.  Then

```math
\boxed{
\rho\ge
(d_+-k_+)+(d_--k_-)-r.
}
\tag{R14.6}
```

The right side may of course be replaced by its positive part.  To prove
(R14.6), quotient `\mathbb R^n` by `K`.  The form induced by `A` is
nondegenerate on the `r`-dimensional quotient, while the images of `L_+`
and `L_-` have dimensions `d_+-k_+` and `d_--k_-`.  The kernel of the
map from the second image to the dual of the first has dimension
`d_--k_--\rho` and lies in an orthogonal complement of dimension
`r-(d_+-k_+)`.  Rearrangement proves the claim.

For a neutral balanced matrix this becomes

```math
\boxed{
(d_+-k_+)+(d_--k_-)\le\operatorname{rank}A.
}
\tag{R14.7}
```

In particular,

```math
\boxed{A\text{ invertible and neutral}\quad\Longrightarrow\quad
d_++d_-\le n.}
\tag{R14.8}
```

This is the correct invertible statement.  With a kernel, the exact quotient
form (R14.7) is preferable to the weaker `d_++d_-\le n+\nullity A`.
The small minimizers below show equality in (R14.6) at ranks zero, four, and
six, so the abstract rank inequality cannot be improved.

There is also a small arithmetic payment.  At balance every entry of
`C_{+-}` lies in `4\mathbb Z`.  A rank-`\rho` matrix has a nonzero
`\rho\times\rho` minor, whose support contains a matching.  Hence one can
choose `\rho` distinct positive grounds and `\rho` distinct negative grounds
and pair them so that

```math
\sum_{j=1}^{\rho}\mathfrak r(p_j,n_j)\ge4\rho.
\tag{R14.9}
```

This is only a linear-size certificate when `\rho=O(n)`, and one endpoint
tree chooses only one root pair.  It is therefore structural information,
not yet a valid temporal harvest.

## 3. Stable frame form, and why it repeats the old scale wall

Let independent laws on the two exact ground clouds have covariance matrices

```math
R_+=\mathbb E(pp^{\mathsf T})\succeq\kappa_+P_+,
\qquad
R_-=\mathbb E(nn^{\mathsf T})\succeq\kappa_-P_-,
```

where `P_\pm` are the Euclidean projections onto `L_\pm`.  At balance,
(R14.4) gives

```math
\begin{aligned}
\mathbb E\,\mathfrak r(p,n)^2
&=\operatorname{tr}(R_+AR_-A)\\
&\ge\kappa_+\kappa_-\|P_+AP_-\|_F^2.
\end{aligned}
\tag{R14.10}
```

If `s_1(A)\ge\cdots\ge s_n(A)` and
`a=(n-d_+)+(n-d_-)=2n-d_+-d_-`, then

```math
\boxed{
\mathbb E\,\mathfrak r(p,n)^2
\ge\kappa_+\kappa_-
\sum_{j>a}s_j(A)^2.
}
\tag{R14.11}
```

Indeed,
`\operatorname{rank}(A-P_+AP_-)\le a`, so Eckart--Young bounds the
Frobenius error of this particular rank-`a` approximation.  Formula
(R14.11) is a stable version of the rank theorem: for invertible `A`, it sees
the smallest `d_++d_--n` singular directions whenever the two ground spans
have dimension surplus.

This does **not** bypass the stopped exact-ground frame route in Section
10.53.2.  With thick linear-dimensional frames and `\Theta(n^2)` captured
spectral mass, it forces at best one endpoint residual of order `n`, not the
`n^{3/2}` leading scale.  More decisively, when `d_++d_-\le n` its spectral
tail is empty.  The `A_8` wall below has `\kappa_+=\kappa_-=1`, is uniformly
invertible, and nevertheless has `d_++d_-=n` and zero residual for every
pair.  Thus “well-conditioned exact-ground frames” alone is not a viable
successor to Section 10.53.

There is a possible, but incomplete, temporal bridge.  For a positive ground
put

```math
f_i^p=p_i(Ap)_i\ge0,\qquad \sum_i f_i^p=Q,
```

and for a negative ground put `f_i^n=-n_i(An)_i\ge0`.  If
`z=p\circ n`, then

```math
p^{\mathsf T}An=\langle f^p,z\rangle=-\langle f^n,z\rangle.
\tag{R14.12}
```

Writing `\Delta_p=\|f^p-(Q/n)\mathbf1\|_2`, balance implies

```math
\boxed{
\mathfrak r(p,n)+\sqrt n\,\Delta_p
\ge\frac Qn|p^{\mathsf T}n|,
}
\tag{R14.13}
```

and likewise with `\Delta_n`.  Thus a Hamming-unbalanced endpoint pair
either has a large residual or a large row-field variance, which is visible
in the row-square term of (10.330).  This does not close the argument:
replenishment remains uncontrolled, and `A_8` has
`p^{\mathsf T}n=0` for every endpoint pair, so it evades (R14.13) completely.

## 4. Exact small-minimizer audit

The checker exhausts every first-row-positive switching representative
through order seven.  The global-minimizer counts and balanced counts are

| `n` | `q_n` | normalized minimizers | balanced minimizers |
|---:|---:|---:|---:|
| 2 | 2 | 1 | 1 |
| 3 | 6 | 2 | 0 |
| 4 | 8 | 6 | 6 |
| 5 | 8 | 12 | 12 |
| 6 | 10 | 12 | 12 |
| 7 | 18 | 3240 | 0 |

All balanced minimizers of each listed order have the same data below.  The
order-eight row refers only to the named matrix (10.445), not to an exhaustive
order-eight classification.

| `n` | projective grounds `g_+/g_-` | `d_+/d_-` | nullity `A` | singular values of `A` | cross-Gram data |
|---:|---:|---:|---:|---|---|
| 2 | `1/1` | `1/1` | 0 | `1^2` | `C=(0)` |
| 4 | `1/1` | `1/1` | 0 | `(\sqrt5)^2,1^2` | `C=(\pm4)`, rank 1 |
| 5 | `5/5` | `5/5` | 1 | `(\sqrt5)^4,0` | rank 4; singular values `(4\sqrt5)^4,0`; 5 zero and 20 magnitude-4 entries |
| 6 | `6/6` | `6/6` | 0 | `(\sqrt5)^6` | rank 6; singular values `(4\sqrt5)^6`; 6 zero and 30 magnitude-4 entries |
| 8 | `4/4` | `4/4` | 0 | `3^6,1^2` | `C=0_{4\times4}` |

For reproducibility, using the ground order emitted by the checker, the
nontrivial order-five and order-six matrices are

```math
C_5=4
\begin{pmatrix}
-1&-1&0&1&-1\\
-1&-1&1&0&1\\
-1&0&-1&-1&1\\
-1&1&-1&1&0\\
0&1&1&1&1
\end{pmatrix},
```

```math
C_6=4
\begin{pmatrix}
-1&-1&0&1&-1&1\\
-1&-1&1&0&1&-1\\
-1&0&-1&-1&1&1\\
-1&1&-1&1&0&-1\\
-1&1&1&-1&-1&0\\
0&1&1&1&1&1
\end{pmatrix}.
```

For the named `A_8`, exact enumeration gives four projective grounds of each
sign and

```math
G_+^{\mathsf T}AG_-=0,
\qquad
G_+^{\mathsf T}G_-=0.
\tag{R14.14}
```

The nonzero eigenvalues of each uniform endpoint covariance are
`3,3,1,1`, so both clouds are quantitatively good frames with lower frame
constant one.  Also

```math
\det A_8=729,
\qquad
\operatorname{spec}(A_8)=(-3)^3,-1,1,3^3.
```

The two four-dimensional endpoint spans are Euclidean complements.  By
(R14.14) they are `A_8`-invariant, and the restricted spectra are

```math
\operatorname{spec}(A_8|_{L_+})=(-1,3,3,3),
\qquad
\operatorname{spec}(A_8|_{L_-})=(-3,-3,-3,1).
\tag{R14.15}
```

Thus `A_8` realizes the sharp neutral equality `d_++d_-=n` with no kernel,
no small singular value, and no thin-frame defect.  It is a precise finite
wall to every argument claiming that balanced endpoint neutrality by itself
must produce degeneracy or poor exact-ground frames.

## 5. Disposition

The useful new exact statement is the identification of the root residual
matrix with the absolute endpoint cross-Gram matrix at balance, together
with the quotient-rank theorem (R14.6).  It creates a clean dichotomy:

1. dimension surplus beyond the nondegenerate ambient rank forces genuinely
   nonneutral endpoint channels (and a rank-sized arithmetic matching); or
2. all-pair neutrality forces the two ground spans to fit as orthogonal
   sectors for the bilinear form of `A`.

The second outcome is real and robust at order eight.  To turn the first
outcome into a temporal comparison one still needs a compatibility theorem
that can use several endpoint pairs, or a mechanism accumulating the
order-`n` frame payment across compatible windows.  Reapplying a bare
second-moment frame inequality would repeat the stopped scale wall.
