# Independent audit: relative Lindeberg defect transfer

## Verdict

The relative-invariance theorem in
`orientation_even_stability_audit.md`, Sections 7--9, is correct for
every fixed \(t\in\mathbb R\) and \(\tau>0\), with constants allowed to
depend on \((t,\tau)\).  I reconstructed the argument without using the
original derivation and found no missing factor.

There is one notation point that should be made explicit.  Before
gauging a pair \(i\ne j\), the scalar entering the one-pair law is

\[
d_{ij}=a_{ij}(A^2)_{ij},
\]

not the ungauged \((A^2)_{ij}\).  After switching the pair so that
\(a_{ij}=1\), this becomes the new \((A^2)_{ij}\).  Since
\(d_{ij}^2=(A^2)_{ij}^2\), the final Frobenius-defect theorem is
unchanged.

Precisely, if \(m=n-1\ge2\),

\[
X_i^\sigma=\operatorname{sign}\left(
\sigma\frac{(A\xi)_i}{\sqrt m}+t\xi_i+\tau Z_i
\right),
\qquad \sigma\in\{\pm1\},
\]

and

\[
f(u)=\frac{\psi_\tau(u+t)+\psi_\tau(u-t)}2,\qquad
K(q)=\mathbb E f(G)f(H),
\]

where \(\psi_\tau(u)=2\Phi(u/\tau)-1\) and
\((G,H)\) are standard Gaussian with correlation \(q\), then there is
\(C_{t,\tau}<\infty\) such that every admissible pair parameter \(d\)
satisfies

\[
\boxed{
\left|C_m(d)-K(d/m)\right|
\le C_{t,\tau}\frac{|d|}{m^{3/2}}.
}
\tag{A}
\]

Consequently, writing
\(\phi_v(t)=(2\pi v)^{-1/2}e^{-t^2/(2v)}\),

\[
\boxed{
\frac12\sum_{\sigma=\pm1}
\mathbb E (X^\sigma)^TA^2X^\sigma
\ge nm+
\left(
4\phi_{1+\tau^2}(t)^2-\frac{C_{t,\tau}}{\sqrt m}
\right)
\frac{\|A^2-mI\|_F^2}{m}.
}
\tag{B}
\]

The coefficient in parentheses can of course be negative at small
\(m\); the assertion of interest is its uniform asymptotic form.

## 1. Pair gauge, including the missing ungauged sign

Fix \(i\ne j\).  Choose a diagonal sign matrix \(D\) with

\[
D_{ii}=1,\qquad D_{jj}=a_{ij},\qquad D_{kk}=a_{ik}\quad(k\ne i,j),
\]

and put \(B=DAD\).  Then

\[
b_{ij}=1,\qquad b_{ik}=1,\qquad
b_{jk}=a_{ij}a_{ik}a_{jk}=:c_k.
\]

Thus

\[
(B^2)_{ij}=\sum_{k\ne i,j}c_k
=a_{ij}(A^2)_{ij}=:d.
\tag{1}
\]

The witness law is switching covariant: under
\((\eta,W)=(D\xi,DZ)\), the witness for \(B\) is \(DX^\sigma\).
Therefore, if

\[
\overline C_{ij}(A)
=\frac12\sum_\sigma\mathbb E X_i^\sigma X_j^\sigma,
\]

then

\[
\overline C_{ij}(B)=a_{ij}\overline C_{ij}(A).
\tag{2}
\]

There are \(m-1=n-2\) common coordinates \(k\ne i,j\).  In this gauge
the two pre-threshold fields are

\[
\sigma m^{-1/2}
\left(\xi_j+\sum_k\eta_k\right)+t\xi_i+\tau Z_i,
\]

\[
\sigma m^{-1/2}
\left(\xi_i+\sum_kc_k\eta_k\right)+t\xi_j+\tau Z_j.
\tag{3}
\]

Permutation invariance of the common Rademachers shows that the
orientation-averaged correlation depends on the coefficient pattern
only through

\[
d=\sum_kc_k\in
\{-(m-1),-(m-3),\ldots,m-1\}.
\]

Call this correlation \(C_m(d)\).

## 2. Exact oddness and parity anchors

In the law with coefficients \(-c_k\), apply

\[
(\eta_k,\sigma,\xi_j,Z_j)
\longmapsto
(-\eta_k,-\sigma,-\xi_j,-Z_j),
\]

leaving \((\xi_i,Z_i)\) fixed.  The first field in (3) is unchanged
and the second is negated.  The transformation preserves the joint
measure, so

\[
\boxed{C_m(-d)=-C_m(d).}
\tag{4}
\]

This proves \(C_m(0)=0\) when \(m-1\) is even.  When \(m-1\) is odd,
the admissible lattice misses zero and oddness instead gives

\[
2C_m(1)=C_m(1)-C_m(-1).
\tag{5}
\]

These are exactly the two anchors needed for telescoping.

## 3. Exact one-coordinate increment and its factor

Flip one common coefficient from \(-1\) to \(+1\).  If the sum of the
other \(m-2\) coefficients is \(r\), this changes the total parameter
from \(r-1\) to \(r+1\).  Put \(s=m^{-1/2}\), integrate out the
independent \(Z_i,Z_j\), and define

\[
\Delta_s(w)=\psi_\tau(w+s)-\psi_\tau(w-s).
\]

After conditioning on all variables except the flipped Rademacher,
write

\[
W_i^\sigma=
\sigma s\left(\xi_j+\sum_{k\ne\ell}\eta_k\right)+t\xi_i,
\]

\[
W_j^\sigma=
\sigma s\left(\xi_i+\sum_{k\ne\ell}c_k\eta_k\right)+t\xi_j.
\]

For a fixed orientation, the difference between the \(+1\) and
\(-1\) coefficient laws is

\[
\frac12
\mathbb E\left[
\Delta_s(W_i^\sigma)\Delta_s(W_j^\sigma)
\right].
\]

The extra \(1/2\) from averaging the two orientations gives the exact
identity

\[
\boxed{
C_m(r+1)-C_m(r-1)
=\frac14\sum_{\sigma=\pm1}
\mathbb E\left[
\Delta_s(W_i^\sigma)\Delta_s(W_j^\sigma)
\right].
}
\tag{6}
\]

Thus the \(1/4\) in the claimed formula is correct.

Because \(\psi_\tau'>0\), (6) also proves monotonicity.  The finite
lower bound in Section 7 of the source note follows as written:
Chebyshev at radius \(2\sqrt{1+t^2}\) puts both residual fields in the
good interval with probability at least \(1/2\), and then each
increment is at least \(d_{t,\tau}^2/m\).  Telescoping from (4)--(5)
gives

\[
d\,C_m(d)\ge \frac{d_{t,\tau}^2}{2m}d^2.
\]

## 4. Uniform two-dimensional Lindeberg estimate

For fixed \((t,\tau)\), all derivatives of \(\psi_\tau\) are bounded.
For \(0<s\le1\) and \(0\le a\le3\),

\[
\|\Delta_s^{(a)}\|_\infty
\le 2s\|\psi_\tau^{(a+1)}\|_\infty.
\tag{7}
\]

Hence every partial derivative of total order at most three of

\[
(u,v)\longmapsto
\Delta_s(u+\alpha)\Delta_s(v+\beta)
\]

is \(O_{t,\tau}(s^2)\), uniformly in the shifts
\(\alpha,\beta\).

The \(m-2\) remaining common Rademachers enter as independent
two-dimensional vectors

\[
s\eta_k(1,c_k).
\]

Replace them one at a time by
\(sg_k(1,c_k)\), \(g_k\sim N(0,1)\).  The first two moments match.
Taylor's theorem through order two, (7), and

\[
\mathbb E\|s\eta_k(1,c_k)\|_2^3+
\mathbb E\|sg_k(1,c_k)\|_2^3=O(s^3)
\]

show that one replacement costs \(O_{t,\tau}(s^5)\).  Therefore all
\(m-2\) replacements cost

\[
O_{t,\tau}(m s^5)
=O_{t,\tau}(m^{-3/2}).
\tag{8}
\]

This proof uses only bounded derivatives; it never divides by a
covariance determinant.  It is therefore uniform when the resulting
two-dimensional Gaussian covariance is singular.

## 5. The Gaussianized increment, including endpoint edges

After replacement, let \((U,V)\) be the remaining common Gaussian
fields.  They have

\[
\operatorname{Var}U=\operatorname{Var}V=1-\frac2m,
\qquad
\operatorname{Cov}(U,V)=\frac rm.
\tag{9}
\]

The direct \(ij\)-edge shifts the two smooth arguments by
\(\sigma s\xi_j\) and \(\sigma s\xi_i\).  The gradient of the product
test in (7) is \(O_{t,\tau}(s^2)\), so deleting these two shifts costs
only

\[
O_{t,\tau}(s^3)=O_{t,\tau}(m^{-3/2}).
\tag{10}
\]

Average the remaining endpoint spins and define

\[
h_s(u)
=\frac{\Delta_s(u+t)+\Delta_s(u-t)}2.
\]

Directly from the definition of \(f\),

\[
h_s(u)=f(u+s)-f(u-s).
\tag{11}
\]

The orientation average in (6) then becomes

\[
\frac12\mathbb E h_s(U)h_s(V)
=2s^2\mathbb E f'(U)f'(V)+O_{t,\tau}(s^4).
\tag{12}
\]

Add independent \(N(0,2/m)\) noise to \(U\) and \(V\).  The resulting
\((G,H)\) have unit variances and covariance \(r/m\).  Bounded
derivatives give

\[
\mathbb E f'(U)f'(V)
=\mathbb E f'(G)f'(H)+O_{t,\tau}(m^{-1/2}).
\tag{13}
\]

This coupling remains valid when (9) is singular.  Gaussian
integration by parts, or Price's identity, gives

\[
K'(q)=\mathbb E f'(G)f'(H).
\tag{14}
\]

Combining (8)--(14) proves the uniform step formula

\[
\boxed{
C_m(r+1)-C_m(r-1)
=\frac2mK'(r/m)+O_{t,\tau}(m^{-3/2}).
}
\tag{15}
\]

The endpoint/direct-edge terms and the variance normalization are
therefore genuinely lower-order at the relative scale claimed.

## 6. Telescoping and the relative error

The function \(f\) is odd, hence \(K\) is odd and \(K(0)=0\).
Moreover \(K''\) is uniformly bounded on \([-1,1]\), since Price's
identity gives it in terms of bounded derivatives of \(f\).  Therefore

\[
K((r+1)/m)-K((r-1)/m)
=\frac2mK'(r/m)+O_{t,\tau}(m^{-2}).
\tag{16}
\]

If \(d\) is even, telescope (15)--(16) from zero.  If \(d\) is odd,
first use (5) and the step \(-1\to1\), then telescope from \(1\).
Oddness handles negative \(d\).  There are \(O(|d|)\) steps, so

\[
|C_m(d)-K(d/m)|
\le C_{t,\tau}|d|m^{-3/2}
+O_{t,\tau}(|d|m^{-2}),
\]

which is (A).  In particular, the estimate remains relative when
\(|d|=1\), and remains uniform when \(|d|=\Theta(m)\).

## 7. From the pair law to \(A^2\)-energy

Expand the orientation-averaged squared-field energy:

\[
\frac12\sum_\sigma
\mathbb E(X^\sigma)^TA^2X^\sigma
=nm+\sum_{i\ne j}(A^2)_{ij}\overline C_{ij}(A).
\tag{17}
\]

By (1)--(2), each off-diagonal summand is

\[
(A^2)_{ij}\overline C_{ij}(A)
=d_{ij}C_m(d_{ij}).
\tag{18}
\]

The Gaussian noise stability of \(f\) has Hermite expansion

\[
K(q)=\sum_{\substack{k\ge1\\k\ {\rm odd}}}w_kq^k,
\qquad w_k\ge0.
\]

Its first coefficient is

\[
w_1=(\mathbb Ef'(G))^2
=4\phi_{1+\tau^2}(t)^2.
\tag{19}
\]

For every \(q\in[-1,1]\),

\[
qK(q)=\sum_{k\ {\rm odd}}w_kq^{k+1}
\ge w_1q^2.
\tag{20}
\]

Apply (A), multiply by \(d\), and use (20):

\[
dC_m(d)
\ge
\left(w_1-\frac{C_{t,\tau}}{\sqrt m}\right)\frac{d^2}{m}.
\tag{21}
\]

Finally,

\[
\sum_{i\ne j}d_{ij}^2
=\sum_{i\ne j}(A^2)_{ij}^2
=\|A^2-mI\|_F^2.
\tag{22}
\]

Substitution of (21)--(22) into (17) proves (B), with every ordered-pair
factor accounted for.

## 8. Exact numerical checks

All checks below used \(t=0.7,\tau=0.8\).  Gaussian expectations were
evaluated by 70-node Gauss--Hermite quadrature; discrete expectations
were exhaustively enumerated.

First, for every admissible \(d\), define

\[
E_m(d)=
\frac{m^{3/2}}{|d|}
\left|C_m(d)-K(d/m)\right|
\quad(d\ne0).
\]

The maximum over admissible \(d\) was:

| \(m\) | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| \(\max_d E_m(d)\) | .18823 | .08412 | .07959 | .06636 | .05999 | .05473 | .05075 | .04753 | .04485 |

This is consistent with the uniform relative estimate (A).

Second, random sign matrices generated with NumPy seed `20260725`
were checked in two independent ways:

1. direct enumeration of every \(\xi\) and both orientations, with the
   Gaussian dithers integrated analytically through \(\psi_\tau\);
2. the gauged pair sum
   \(nm+\sum_{i\ne j}d_{ij}C_m(d_{ij})\).

| \(n\) | defect \(\|A^2-mI\|_F^2\) | direct energy | pair-sum energy | difference |
|---:|---:|---:|---:|---:|
| 3 | 6 | 7.269760267934 | 7.269760267934 | \(8.9\cdot10^{-16}\) |
| 4 | 16 | 13.816009738042 | 13.816009738042 | \(3.6\cdot10^{-15}\) |
| 5 | 52 | 24.190438157496 | 24.190438157496 | \(3.6\cdot10^{-15}\) |
| 6 | 64 | 33.927554564423 | 33.927554564423 | \(7.1\cdot10^{-15}\) |
| 7 | 218 | 53.217573348779 | 53.217573348779 | \(2.1\cdot10^{-14}\) |
| 8 | 256 | 67.104252792437 | 67.104252792437 | \(9.9\cdot10^{-14}\) |

Finally, the exact increment identity (6) was checked for
\(m=2,\ldots,7\) and agreed to at worst \(1.7\cdot10^{-16}\).

## 9. Scope of the certification

The audit certifies the relative pair theorem and its conversion to
the leading-scale \(A^2\)-energy defect coefficient.  It does not
supply the missing positive-heavy-field tail estimate needed to turn
that defect energy into a strictly larger universal bound for
\(Q(A)\).  The downstream bottleneck identified in Section 9 of the
source note is therefore real, but it is not caused by an error in the
relative Lindeberg transfer.
