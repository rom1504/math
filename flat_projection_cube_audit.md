# Flat-diagonal cube projections: exact audit

## 1. Frame identity and a universal lower bound

Let \(P\) be an orthogonal projection of order \(n=2k\), rank \(k\),
with

\[
P_{ii}=\frac12.
\]

Write \(P=WW^\top\), where \(W\) has orthonormal columns.  If \(w_i\)
are the rows of \(W\), put \(u_i=\sqrt2\,w_i\).  Then

\[
\|u_i\|_2=1,\qquad \sum_{i=1}^{2k}u_i u_i^\top=2I_k.
\]

Duality gives the exact identity

\[
\boxed{
\max_{x\in\{\pm1\}^{n}}x^\top Px
=
\left(
\max_{\|z\|_2=1}\|Wz\|_1
\right)^2
=
\frac12
\left(
\max_{\|z\|_2=1}\sum_i|u_i\cdot z|
\right)^2.
}                                                            \tag{1.1}
\]

For \(z\) uniform on \(S^{k-1}\),

\[
\mathbb E|u_i\cdot z|
=
a_k
:=
\frac{\Gamma(k/2)}
{\sqrt\pi\,\Gamma((k+1)/2)}.
\]

Therefore

\[
\boxed{
\frac1n\max_xx^\top Px
\ge
k a_k^2
=
\frac2\pi+O(k^{-1}).
}                                                            \tag{1.2}
\]

This is the strongest bound obtained from the first spherical moment.
It is not known here to be the sharp universal constant.

---

## 2. A rank-two counterexample to the Haar threshold

Let

\[
u_j=
\left(
\cos\frac{j\pi}{4},
\sin\frac{j\pi}{4}
\right),
\qquad j=0,1,2,3,
\]

and let \(W\) have rows \(u_j/\sqrt2\).  Then \(P=WW^\top\) is a
rank-\(2\) projection of order \(4\) with diagonal \(1/2\).  Directly,

\[
P=
\begin{pmatrix}
\frac12&\frac1{2\sqrt2}&0&-\frac1{2\sqrt2}\\
\frac1{2\sqrt2}&\frac12&\frac1{2\sqrt2}&0\\
0&\frac1{2\sqrt2}&\frac12&\frac1{2\sqrt2}\\
-\frac1{2\sqrt2}&0&\frac1{2\sqrt2}&\frac12
\end{pmatrix}.
\]

Exhausting the eight sign rays, or applying (1.1), gives

\[
\boxed{
\max_xx^\top Px=2+\sqrt2,
\qquad
\frac14\max_xx^\top Px
=
\frac{2+\sqrt2}{4}
=0.8535533905\ldots.
}                                                            \tag{2.1}
\]

This is strictly below the Haar threshold

\[
\beta_*=
\frac12\left(1+\frac{\sqrt{15}}4\right)
=0.9841229182\ldots.
\]

Direct sums of this projection give an infinite balanced-diagonal
counterfamily with the same ratio.  Thus \(\beta_*\) is not universal
for flat-diagonal half-dimensional subspaces.

---

## 3. A coherent order-\(16\) counterexample

Put

\[
U_0=2P-I_4
=
\frac1{\sqrt2}
\begin{pmatrix}
0&1&0&-1\\
1&0&1&0\\
0&1&0&1\\
-1&0&1&0
\end{pmatrix},
\]

so \(U_0^\top=U_0\), \(U_0^2=I\), and \(\operatorname{diag}U_0=0\).
Let \(H_4\) be a symmetric Hadamard matrix and \(V=H_4/2\).  Define

\[
U=U_0\otimes V,\qquad P_{16}=\frac{I_{16}+U}{2}.
\]

Then \(P_{16}\) is a rank-\(8\) projection with diagonal \(1/2\), and

\[
\max_{i\ne j}|U_{ij}|=\frac1{\sqrt8}
=\frac{\sqrt2}{\sqrt{16}}.
\]

Thus it has the natural \(O(n^{-1/2})\) coherence scale, although it
has structural zeros and is not an exact conference projection.

To compute its cube norm, write a sign vector as four sign vectors
\(x_1,\ldots,x_4\in\{\pm1\}^4\).  Optimizing \(x_2,x_4\) first gives

\[
\max x^\top Ux
=
2\sqrt2
\max_{a,c\in\{\pm1\}^4}
\left[
\|V p\|_1+\|V r\|_1
\right],
\]

where

\[
p=\frac{a+c}{2},\qquad r=\frac{c-a}{2}
\]

have complementary supports.  If \(|\operatorname{supp}p|=s\), the
four possibilities are:

\[
\begin{array}{c|c}
s&\max(\|Vp\|_1+\|Vr\|_1)\\ \hline
0,4&4\\
1,3&5\\
2&4.
\end{array}
\]

For \(s=1,3\), a three-supported vector has one Walsh coefficient of
magnitude \(3/2\) and three of magnitude \(1/2\), while the singleton
has \(\ell_1\)-transform norm \(2\).  Hence the maximum is \(5\), and

\[
\boxed{
\max_xx^\top Ux=10\sqrt2,
\qquad
\frac1{16}\max_xx^\top P_{16}x
=
\frac12+\frac{5\sqrt2}{16}
=0.9419417382\ldots<\beta_*.
}                                                            \tag{3.1}
\]

So even \(O(n^{-1/2})\) coherence at a finite order does not force the
Haar threshold.  This construction does not give an asymptotic
counterexample with exact entry-flatness.

---

## 4. Exact conference specialization

Let \(C\) be a symmetric conference matrix:

\[
C^2=(n-1)I,\qquad U=\frac C{\sqrt{n-1}},
\qquad P=\frac{I+U}{2}.
\]

Then \(P\) is a rank-\(n/2\) projection with

\[
P_{ii}=\frac12,\qquad
|2P_{ij}|=\frac1{\sqrt{n-1}}\quad(i\ne j).
\]

The frame in Section 1 is now an equiangular tight frame with

\[
|u_i\cdot u_j|=\rho:=\frac1{\sqrt{n-1}}.
\]

For \(z\) uniform on \(S^{k-1}\),

\[
\mathbb E|u_i\cdot z||u_j\cdot z|
=
\frac2{\pi k}
\left(
\sqrt{1-\rho^2}+\rho\arcsin\rho
\right).
\]

Consequently, the exact second spherical moment yields only

\[
\boxed{
\frac1n\max_xx^\top Px
\ge
\frac1n+
\frac{2(n-1)}{\pi n}
\left(
\sqrt{1-\rho^2}+\rho\arcsin\rho
\right)
=
\frac2\pi+O(n^{-1}).
}                                                            \tag{4.1}
\]

This is much weaker than the already proved conference rounding bound

\[
\frac1n\max_xx^\top Px
\ge
0.8364933644\ldots.
\]

Fixed Boolean moments are also too small.  For

\[
q(x)=x^\top Ux,
\]

conference orthogonality gives

\[
\boxed{
\mathbb E_xq(x)^2=2n,
}
                                                               \tag{4.2}
\]

and

\[
\boxed{
\mathbb E_xq(x)^4
=
\frac{4n(3n^2-15n+20)}{n-1}
=12n^2+O(n).
}                                                            \tag{4.3}
\]

For (4.3), if

\[
\mathcal Q_4(C)
=
\sum_{\text{undirected four-cycles}}
\prod_{e\in C_4}c_e,
\]

then \(\operatorname{tr}C^4=n(n-1)^2\) gives

\[
\mathcal Q_4(C)=-\frac{n(n-1)(n-2)}8.
\]

Substitution into

\[
\mathbb E
\left(\sum_{i<j}c_{ij}x_ix_j\right)^4
=
3E_n^2-2E_n+24\mathcal Q_4(C)
\]

gives (4.3).  Thus every fixed second/fourth-moment argument sees
\(\max q\) only at order \(\sqrt n\), whereas the desired projection
threshold requires \(q=\Theta(n)\).  Exact entry-flatness enters only
perturbatively at fixed moment order; the \(\beta_*\) question is a
linear-order, nonperturbative tail problem.

No proof of

\[
\liminf
\frac1n\max_xx^\top Px\ge\beta_*
\]

for all exact conference projections was found, nor was an asymptotic
exact-conference counterfamily found.

---

## 5. Connection to the original signing problem and stability needed

For a symmetric conference matrix,

\[
\frac1n\max_xx^\top Px
=
\frac12+
\frac1{2n\sqrt{n-1}}\max_xx^\top Cx.
                                                               \tag{5.1}
\]

Applying the same statement to \(-C\),

\[
\frac{\max_x|x^\top Cx|}
{n\sqrt{n-1}}
=
2\max\{\beta(P_C),\beta(P_{-C})\}-1.
                                                               \tag{5.2}
\]

Therefore a universal exact-conference theorem with threshold
\(\beta_*\) would imply the one-copy conference lower constant

\[
\frac12(2\beta_*-1)=\frac{\sqrt{15}}8.
\]

It would not by itself prove the same bound for \(F(n)\), because a
minimizing signing need not be conference.

A direct robust extension could start with

\[
P_A=
\frac12
\left(
I+\frac A{\sqrt{n-1}}
\right).
\]

This matrix has the exact desired diagonal and off-diagonal
entry-flatness, but it is a projection only when
\(A^2=(n-1)I\).  A stability theorem would need to cover approximate
projections under a condition such as

\[
\boxed{
\|A^2-(n-1)I\|_{\rm op}=o(n).
}                                                            \tag{5.3}
\]

Then

\[
\|P_A^2-P_A\|_{\rm op}=o(1),
\]

and a robust projection inequality could transfer directly to the
quadratic form of \(A\).

Current near-minimizer information is far weaker.  The bound

\[
\max_x|x^\top Ax|
\ge\frac{\|A\|_{\rm op}^3}{n}
\]

only gives \(\|A\|_{\rm op}=O(n^{5/6})\) for competitors of order
\(n^{3/2}\).  It implies mean-square row pseudorthogonality, but not
(5.3).  Wigner sign matrices and conference matrices both pass that
mean-square test while having very different Boolean maxima.

Hence the exact-conference projection problem is only one branch of a
full lower-bound proof.  Away from (5.3), a separate entropic or
spin-glass lower bound is still required.
