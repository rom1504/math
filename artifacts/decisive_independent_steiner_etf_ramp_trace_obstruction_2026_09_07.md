# Actual full-sign Steiner matrices obstruct every sublinear-trace PSD ramp

Date: 2026-09-07. Status: proved from the primary Steiner-frame existence
theorem and the direct calculations below. This is a scoped falsifier of
ramp availability for arbitrary bounded-cap/full-sign parents. It is not
a statement about actual width minimizers or impossibility of other sign
surgeries, and no literature novelty claim is made.

## 1. The infinite real family and the exact imported hypotheses

For an integer j>=1 put

```math
r=16^j-1,\qquad v=4r+1,\qquad d={vr\over5},
\qquad n=v(r+1)=4r^2+5r+1.
```

The primary source is Fickus--Mixon--Tremain, *Steiner equiangular tight
frames*, [arXiv:1009.5730](https://arxiv.org/pdf/1009.5730), Theorems 1--2.
Theorem 1 is on PDF page 4. Theorem 2 provides a Steiner system S(2,k,v)
for fixed k and all sufficiently large admissible v; Theorem 1 gives a
real frame when a real Hadamard of order 1+(v-1)/(k-1) exists. Here k=5,
(v-1)/4=r is integral, and vr/5 is integral because 16^j-1 is divisible
by 5. Indeed v is 1 modulo 20. The required Hadamard order r+1=16^j is
an unconditional Sylvester order. Thus the source supplies real frames
at every sufficiently large member of this particular infinite sequence.
No Hadamard conjecture, variable-k design theorem, or density assertion is
being imported.

In fact Section 3.1.4 (PDF page 9) explicitly supplies S(2,5,v) for every
admissible v congruent to 1 or 5 modulo 20. This stronger fixed-block-size
statement removes the sufficiently-large qualifier for our sequence; the
weaker Theorem 2 already suffices for the asymptotic obstruction.

For clarity, the frame algebra can be reconstructed directly. A Steiner
system has d blocks, each of size 5, and each point belongs to r blocks.
At each point, insert the r nonconstant rows of a normalized-sign Sylvester
matrix of order r+1 into its incident block positions. Concatenate these
v blocks of columns and divide by sqrt(r), obtaining a real d by n matrix F.
Distinct frame rows are orthogonal: two design blocks share at most one
point, and at that point distinct Hadamard rows were used. Each frame row
has squared norm 5(r+1)/r. Each column has norm one. Two columns at the same
point have inner product -1/r, and columns at distinct points have inner
product +1/r or -1/r, since exactly one design block contains that pair.
Hence, with G=F^T F,

```math
G_{ii}=1,\quad G_{ij}\in\{-1/r,1/r\}\ (i\ne j),\qquad
FF^T={5(r+1)\over r}I_d.
```

## 2. Actual signs and endpoint bounds

Define A=r(G-I_n). This is an actual hollow symmetric full-sign matrix.
Its exact eigenvalues are

```math
4r+5\quad\text{(multiplicity d)},\qquad
-r\quad\text{(multiplicity n-d)}.
```

Write H_A(x)=x^T A x/2, P=max H_A, R=max(-H_A), w=(P+R)/2,
and I=(P-R)/2. The spectral bounds give

```math
P\le {n(4r+5)\over2},\quad R\le {nr\over2},\quad
w\le {n(5r+5)\over4},\quad I\le {n(4r+5)\over4}.
```

Since sqrt(n)/r tends to 2, this is a bounded-operator and bounded-cap
family: ||A||op/sqrt(n)->2 and Q(A)/n^(3/2)<=1+o(1). Its width upper bound
is (5/8+o(1))n^(3/2).

## 3. A bounded-covariance Boolean witness above that width bound

Take a centered Gaussian Z with covariance G, and X_i=sign(Z_i).
Every coordinate has variance one, so no ambiguity at zero occurs.
The two-dimensional Gaussian angle formula gives
E X_i X_j=(2/pi)arcsin(G_ij). Put

```math
a_r={2\over\pi}\arcsin(1/r).
```

Because all nonzero off-diagonal magnitudes are exactly equal, the entire
Boolean covariance is the exact affine matrix

```math
K:=\mathbb E XX^T=I_n+a_r A,\qquad
\|K\|_{\rm op}=1+a_r(4r+5)\longrightarrow1+{8\over\pi}.
```

No approximation by entrywise powers or independence assertion is needed.
Using Tr A=0 and Tr A^2=n(n-1),

```math
\mathbb E H_A(X)={a_r n(n-1)\over2}
               =\left({2\over\pi}+o(1)\right)n^{3/2}.
```

This is strictly larger than the width upper bound because
2/pi>5/8. In particular P>=E H_A(X) and R<=nr/2 ensure I>0 eventually,
so the positive orientation is the one required by the ramp formulation.

## 4. Every PSD ramp has linear trace

Let T be ANY positive-semidefinite matrix, not necessarily a projector or
a contraction. Suppose that for every Boolean x,

```math
H_A(x)\le w+{I\over n}x^T T x+e_n,
\qquad e_n=o(n^{3/2}).
```

Averaging the actual Boolean law in Section 3 and using positivity gives

```math
\mathbb E H_A(X)-w-e_n
\le {I\over n}\operatorname{Tr}(TK)
\le {I\over n}\|K\|_{\rm op}\operatorname{Tr}T.
```

The numerator on the left is positive eventually. Insert the endpoint
and covariance bounds from Sections 2--3. It follows that

```math
\liminf_j{\operatorname{Tr}T\over n}
\ge {2/\pi-5/8\over(1/2)(1+8/\pi)}
= {16-5\pi\over4(\pi+8)}>{1\over156}.
```

The last inequality follows already from pi<22/7. Thus this family has
no o(n^(3/2))-error ramp with trace o(n), even before imposing compatibility,
entry bounds, masking, idempotence, or an operator bound on T. In particular
it escapes all three proposed low-complexity ramp surgery regimes for
arbitrary bounded-cap parents, including the trace-only sign-compatible
variant.

This does NOT refute the conditional surgery theorems: their certificate
premise is absent here. It also does not establish that these matrices
minimize width, so it does not refute a width-minimizer-specific certificate
or the conjectured comparison M_n-W_n=o(n^(3/2)). Nor does it rule out an
actual sign repair not represented by such a ramp.

## 5. An explicit finite instance for reproducibility

At r=15, v=61, a concrete cyclic S(2,5,61) is obtained by translating
each of the following three blocks through Z/61Z:

```text
{1,9,20,34,58}, {4,14,19,36,49}, {13,15,16,22,56}.
```

Their ordered differences partition all 60 nonzero residues, so every
unordered point pair lies in exactly one of the 183 translated blocks.
Using H_16 produces a real 183 by 976 frame and an actual order-976 sign
matrix with eigenvalues 65 and -15. These claims can all be checked by
integer arithmetic before any Gaussian calculation.

For this instance, elementary bounds arcsin(1/r)>=1/r, 3<pi<22/7, and
arcsin(1/r)<=1/(r-1) give E H_A(X)>=222040/11, w<=19520, and
||K||op<=86/21. Every zero-error PSD ramp therefore has
Tr T>=61488/6149. The finite check is supplementary; the linear-trace
obstruction is the infinite-sequence proof above.

The program `computations/decisive_independent_steiner_etf_ramp_check_2026_09_07.py`
passed: all 1830 point pairs covered exactly once, exact 183 by 976 row
Gram, unit-sign off-diagonal Gram entries, both eigenvalue multiplicity
identities, and the displayed rational finite trace bound.
