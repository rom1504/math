# A regular-Hadamard obstruction to terminal-stability bounds

## Status and scope

This note gives an explicit infinite family showing that the following data do
**not** imply a small row-sign recoupling defect:

1. a full symmetric signing with `Q(A)=O(n^(3/2))`;
2. the row-sign shore identities at a particular input `X`;
3. both anchored shores jointly; and
4. strict one-coordinate stability of the two augmented-shore outputs.

The resulting possible local-search defect is

```math
\left({1\over4\sqrt2}+o(1)\right)n^{3/2},
```

whose leading coefficient `0.176776...` exceeds the available doubled-scale
budget `sqrt(2/pi)-c_*=0.124897...`.  This is a rigorous obstruction to an
argument based only on project-scale cap and the terminal local-field
identity.  It is **not** a counterexample to the expectation-over-`X`
row-sign local-repair lemma: the construction below singles out one `X`, and
the prescribed deterministic ascent is not proved to enter the bad local
optimum at every order.

All energies use doubled normalization

```math
Q(A)=\max_{z\in\{\pm1\}^n}|z^{\mathsf T}Az|.
```

## 1. A regular symmetric Hadamard matrix

Fix `k>=2`, put `s=2^k` and `m=s^2`, and identify the index set with

```math
V=\mathbb F_2^k\times\mathbb F_2^k.
```

For `u=(a,b)` define `q(u)=a\mathbin\cdot b` in `F_2`, and set

```math
K_{u,v}=(-1)^{q(u+v)}.                                  \tag{1}
```

Then `K` is symmetric, has diagonal `+1`, and

```math
K^2=mI,
\qquad K\mathbf1=s\mathbf1.                            \tag{2}
```

Indeed, polarization gives

```math
q(u+v)=q(u)+q(v)+\beta(u,v),
\qquad
\beta((a,b),(c,d))=a\mathbin\cdot d+c\mathbin\cdot b.
```

Thus `K=DHD`, where `D_uu=(-1)^{q(u)}` and
`H_uv=(-1)^{beta(u,v)}` is the character table of the nondegenerate
symplectic pairing `beta`.  Character orthogonality gives `H^2=mI` and
hence `K^2=mI`.  Translation in (1) gives a constant row sum, equal to the
quadratic Gauss sum

```math
\sum_{a,b}(-1)^{a\cdot b}
=\sum_a\sum_b(-1)^{a\cdot b}=2^k=s.
```

Consequently

```math
C=K-I                                                     \tag{3}
```

is a symmetric zero-diagonal `+-1` signing satisfying

```math
C\mathbf1=(s-1)\mathbf1,
\qquad \|C\|_{\rm op}=s+1.                              \tag{4}
```

## 2. An explicit strict local optimum

Write coordinates of `a,b` as `a_0,a_1,...` and `b_0,b_1,...`.  Let

```math
S=\{(a,b): b_1=a_1,\ b_0=a_0+a_1\}.                    \tag{5}
```

This is a codimension-two subspace, so `|S|=m/4`.  Define

```math
r=\mathbf1-2\mathbf1_S.                                 \tag{6}
```

We claim the exact Walsh identity

```math
K\mathbf1_S={s\over2}r.                                 \tag{7}
```

To prove it, fix `u=(a,b)`.  Summing the coordinates numbered at least two
in

```math
(K\mathbf1_S)_u=\sum_{v\in S}(-1)^{q(u+v)}
```

contributes `2^(k-2)=s/4`.  After translating the two remaining free
coordinates to `x=a_0+c_0` and `y=a_1+c_1`, their contribution is

```math
\sum_{x,y\in\mathbb F_2}
(-1)^{xy+Ax+By}=2(-1)^{AB},                              \tag{8}
```

where

```math
A=b_0+a_0+a_1+1,
\qquad B=b_1+a_1+1.
```

Now `AB=1` exactly when both defining equations in (5) hold.  Hence (8) is
`2r_u`, proving (7).

Since `K1=s1`, equations (6)--(7) imply

```math
Kr=s\mathbf1-2K\mathbf1_S=2s\mathbf1_S.                \tag{9}
```

Using `C=K-I`, we therefore obtain the coordinatewise identity

```math
r_u(Cr)_u=
\begin{cases}
-1,&u\notin S,\\
-(2s+1),&u\in S.
\end{cases}                                             \tag{10}
```

In particular, `r` is a strict one-flip local maximum of
`F(r)=-r^T C r`: flipping coordinate `u` changes `F` by
`4r_u(Cr)_u<0`.  Summing (10) gives its exact value

```math
F(r)=-r^{\mathsf T}Cr=m\left({s\over2}+1\right).        \tag{11}
```

Thus even the full terminal identity
`F(r)=sum_u |(-Cr)_u|`, together with strict stability, permits a value only
about half of the natural `m^(3/2)` shore scale.

## 3. Embedding both bad shores in a full project-scale signing

Let `H_0` be a Sylvester Hadamard matrix of order `m/2`, let

```math
L=\begin{pmatrix}1&-1\\-1&1\end{pmatrix},
\qquad B=H_0\otimes L,                                  \tag{12}
```

and define the order-`N=2m` full signing

```math
A=\begin{pmatrix}C&B\\B^{\mathsf T}&-C\end{pmatrix}.    \tag{13}
```

Every off-diagonal entry of `A` is `+-1`, its diagonal is zero, and it is
symmetric.  Since `L1=0`,

```math
B\mathbf1=B^{\mathsf T}\mathbf1=0.                     \tag{14}
```

Moreover `||H_0||_op=sqrt(m/2)` and `||L||_op=2`, so

```math
\|B\|_{\rm op}=\sqrt{2m}=\sqrt2s.
```

Splitting (13) into its block-diagonal and off-diagonal parts gives

```math
\|A\|_{\rm op}\le (s+1)+\sqrt2s.
```

Consequently

```math
Q(A)\le N\|A\|_{\rm op}
\le2m\bigl((1+\sqrt2)s+1\bigr)=O(N^{3/2}).              \tag{15}
```

Thus this is genuinely a project-scale family, not a hidden quadratic-cap
example.

Take `X=1`.  By (4) and (14),

```math
AX=((s-1)\mathbf1,-(s-1)\mathbf1).
```

The row-sign law therefore chooses the first block as the agreement shore
`I` and the second as the disagreement shore `J`.  With `p=q=1`,

```math
P=p^{\mathsf T}A[I]p=m(s-1),
\qquad
R=q^{\mathsf T}A[J]q=-m(s-1),                           \tag{16}
```

and both anchored cross fields vanish by (14).  The two sign-specific
augmented-shore objectives are consequently identical: each asks for a
one-sided lower witness for `-r^T C r` (plus an isolated collapsed spin).
The vector (6) is a strict stable point for both objectives.

If the two local searches terminate at these stable points, their joint
certificate and clipped defect are

```math
K=P+m\left({s\over2}+1\right),
\qquad
\Delta=(P-R)-K
=m\left({s\over2}-2\right).                             \tag{17}
```

Since `N=2m` and `s=sqrt(m)`,

```math
{\Delta\over N^{3/2}}
={1\over4\sqrt2}-{1\over\sqrt2s}
\longrightarrow {1\over4\sqrt2}=0.176776695\ldots .   \tag{18}
```

This exceeds `0.124897...`, so no universal terminal-stability estimate at
the required constant follows from (15)--(16), even after taking the better
of the two anchors.

## 4. What is and is not falsified

The algebra above proves existence of a bad strict local optimum at every
order `m=4^k`, and proves that both anchored problems possess the same bad
point.  It falsifies any proposed proof using only:

- `Q(A)=O(N^(3/2))` or an operator-norm substitute;
- row-sign shore signs and energies at a fixed input;
- the one-flip terminal inequalities; and
- selection of the better anchored shore.

It does not prove that the deterministic best-improvement algorithm with its
specified tie rule reaches this point for every `k`.  Direct exact replay
from the all-ones initialization does reach (6), with terminal values
`48,320,2304` and flip counts `4,16,64`, for `m=16,64,256`, respectively.
That finite observation is reproducible with
`computations/check_regular_hadamard_local_obstruction.py`; it is deliberately
not promoted to an asymptotic theorem.

Nor does one exceptional input `X=1` control the uniform average over `X`.
A proof of the expected row-sign local-repair lemma could still exploit the
distribution of switched shores or a basin-of-attraction invariant of the
particular greedy dynamics.  The present theorem says that such additional
information is necessary.
