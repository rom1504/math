# Thick near-cap rank and local-field thresholding

## 1. Purpose and outcome

Continue with the top/bottom decomposition

\[
A=c+h,\qquad \widetilde A=c-h
\]

from `dual_cut_cone_diffuse_midpoint.md`.  The cross-supported part
\(c\) has rectangular matrix \(C\), the internal part is
\(h=B\oplus D\), and

\[
\|C\|_{\infty\to1}=W,\qquad h\cdot\mathbf1=d.
\]

The exact capped-cone inequality is

\[
\boxed{
|d-H_B(p)-H_D(q)|
\le W-|p^\mathsf TCq|.
}
\tag{1.1}
\]

The exact cap layer can have almost zero feature rank even for global
width minimizers.  This note therefore studies Hamming balls around a
cap state.

The positive result is sharp:

* singleton and double flips of low-field vertices enter an explicit
  thick cap layer;
* their internal cut features span the entire edge space of the
  low-field core;
* the corresponding feature operator has smallest singular value
  exactly \(2\), independent of the core size.

The obstruction is also sharp: the allowed slacks on these spanning
features are of order of the local fields, typically \(\Theta(\sqrt
n)\).  Full rank therefore does not make the midpoint tangent small.
Thresholding high-field vertices does not repair the scale because
their number/mass tradeoff is exactly critical.

## 2. Local fields at a cross cap

Switch within \(U,V\) so that

\[
\mathbf1^\mathsf TC\mathbf1=W.
\tag{2.1}
\]

Define the cross row and column fields

\[
r_i=\sum_{j\in V}C_{ij}\quad(i\in U),
\qquad
s_j=\sum_{i\in U}C_{ij}\quad(j\in V).
\tag{2.2}
\]

One-spin optimality for the rectangular norm gives

\[
r_i,s_j\ge0,
\qquad
\sum_{i\in U}r_i=\sum_{j\in V}s_j=W.
\tag{2.3}
\]

For \(X\subseteq U,Y\subseteq V\), flip precisely those cross-spin
coordinates.  Put

\[
\Gamma(X,Y)
=r(X)+s(Y)-2C(X,Y),
\tag{2.4}
\]

where \(C(X,Y)=\sum_{i\in X,j\in Y}C_{ij}\).  Then

\[
p^\mathsf TCq=W-2\Gamma(X,Y).
\tag{2.5}
\]

Since \(|p^\mathsf TCq|\le W\),

\[
0\le\Gamma(X,Y)\le W.
\tag{2.6}
\]

Writing

\[
\partial_B(X)
=\sum_{i\in X,\ i'\in U\setminus X}B_{ii'},
\qquad
\partial_D(Y)
=\sum_{j\in Y,\ j'\in V\setminus Y}D_{jj'},
\]

equation (1.1) becomes

\[
\boxed{
|\partial_B(X)+\partial_D(Y)|
\le
\min\{\Gamma(X,Y),W-\Gamma(X,Y)\}.
}
\tag{2.7}
\]

For \(Y=\varnothing\), this specializes to

\[
\boxed{
|\partial_B(X)|
\le
\min\{r(X),W-r(X)\}.
}
\tag{2.8}
\]

There is an identical column-side inequality for \(D\).

## 3. Singleton and double-flip layers

Let

\[
R_i=\sum_{\substack{i'\in U\\i'\ne i}}B_{ii'}
\]

be the internal row sum at \(i\).  If \(r_i\le W/2\), (2.8) gives

\[
\boxed{|R_i|\le r_i.}
\tag{3.1}
\]

For distinct \(i,j\in U\), put

\[
R_{ij}=\partial_B(\{i,j\})=R_i+R_j-2B_{ij}.
\tag{3.2}
\]

If \(r_i+r_j\le W/2\), then

\[
\boxed{|R_{ij}|\le r_i+r_j.}
\tag{3.3}
\]

Fix a threshold \(T\le W/4\) and define the low-field core

\[
L_T=\{i\in U:r_i\le T\}.
\tag{3.4}
\]

Every singleton cut in \(L_T\) has absolute cross slack at most
\(2T\), and every double cut has slack at most \(4T\).  Therefore the
\(4T\)-near-cap layer contains all their internal feature vectors.

## 4. Exact feature rank and singular values

Let \(k=|L_T|\), let \(E_k=E(K_k)\), and work in
\(\mathbb R^{E_k}\).  For \(i\in[k]\), let \(s_i\) be the incidence
vector of the full star at \(i\).  For \(ij\in E_k\), let

\[
t_{ij}=\delta(\{i,j\})=s_i+s_j-2e_{ij}.
\tag{4.1}
\]

Thus

\[
\boxed{
e_{ij}=\frac{s_i+s_j-t_{ij}}2.
}
\tag{4.2}
\]

The singleton and double-cut features span all of
\(\mathbb R^{E_k}\).

There is a dimension-free quantitative statement.  Let \(S\) be the
matrix with rows \(s_i\), let \(T_2\) be the matrix with rows
\(t_{ij}\), and put

\[
\mathcal F=
\begin{pmatrix}
S\\T_2
\end{pmatrix}.
\tag{4.3}
\]

If \(L(K_k)\) denotes the line graph of \(K_k\), with adjacency matrix
\(J_k^{\rm line}\), then

\[
S^\mathsf TS=2I+J_k^{\rm line},
\qquad
T_2=J_k^{\rm line}.
\]

The three eigenvalues of \(J_k^{\rm line}\) are

\[
2(k-2),\qquad k-4,\qquad -2
\]

on the constant, vertex-standard, and cycle subspaces.  Consequently

\[
\mathcal F^\mathsf T\mathcal F
=2I+J_k^{\rm line}+(J_k^{\rm line})^2
\tag{4.4}
\]

has eigenvalues

\[
\begin{aligned}
&2+2(k-2)+4(k-2)^2,\\
&2+(k-4)+(k-4)^2,\\
&4.
\end{aligned}
\tag{4.5}
\]

In particular,

\[
\boxed{
\sigma_{\min}(\mathcal F)\ge\sqrt2,
\qquad
\sigma_{\min}(\mathcal F)=2\ \text{for }k\ge5.
}
\tag{4.6}
\]

Thus failure of the thick-cap route is not a rank loss or a poor
restricted-invertibility constant.

## 5. Why the full-rank constraints remain too weak

Equations (3.1)--(3.3), together with (4.2), give

\[
|B_{ij}|
\le\frac{|R_i|+|R_j|+|R_{ij}|}{2}
\le r_i+r_j.
\tag{5.1}
\]

But \(B_{ij}=\pm1\), while a competitive cross cap has typical field

\[
r_i\asymp\frac W{|U|}=\Theta(\sqrt n).
\tag{5.2}
\]

Hence even the exact coordinate inverse is weaker than coefficient
flatness by a factor \(\sqrt n\).

The global \(\ell_2\) version has the same loss.  If all \(r_i\le T\)
on a \(k\)-vertex core, the vector of singleton and double-cut
evaluations has squared norm at most

\[
kT^2+4\binom k2T^2=O(k^2T^2).
\tag{5.3}
\]

Even when \(L_T\) has a boundary in \(U\), the exact inverse (4.2)
cancels all boundary-edge contributions:

\[
B_{ij}=\frac{R_i+R_j-R_{ij}}2.
\]

Summing the square of this identity over \(ij\in E(L_T)\), together
with (5.3), yields only

\[
\|B[L_T]\|_F=O(kT),
\tag{5.4}
\]

where the true flat norm is already
\(\Theta(k)\).  It becomes informative only at \(T=O(1)\), while the
mean field is \(\Theta(\sqrt n)\).

## 6. Threshold decomposition

From (2.3),

\[
\boxed{
|\{i:r_i>T\}|\le\frac WT.
}
\tag{6.1}
\]

If \(W=O(n^{3/2})\) and

\[
T=K\sqrt n,
\]

then the exceptional high-field set has size \(O(n/K)\).  Taking
\(K\to\infty\) makes it \(o(n)\), so the global block-replacement
theorem can remove a truly localized spike there.

However, the low-field margin in (5.3) becomes

\[
T=K\sqrt n,
\]

which is even larger than the natural \(\sqrt n\) scale.  In the other
direction, taking \(T=o(\sqrt n)\) would make the thick-cap constraints
quantitatively sharp, but (6.1) no longer makes the exceptional set
small: it may contain a linear fraction of all vertices.

This is the exact critical tradeoff:

\[
\boxed{
\text{small exceptional set}\ \Longleftrightarrow\
\text{non-informative thick-layer margins}.
}
\tag{6.2}
\]

The induced-block exchange inequality

\[
W(A[S])\le W_{|S|}
+2\|A_{S,S^c}\|_{\infty\to1}
\tag{6.3}
\]

rules out a high-field set only when its induced width dominates its
boundary norm.  The scalar mass bound (6.1) supplies no such boundary
control.

## 7. Exact small-order face audit

The program `cut_cap_face_rank_audit.cpp` exhaustively enumerates every
switching class through \(n=8\), selects all global half-range
minimizers, runs over every top/bottom endpoint pair, and computes the
real rank of the internal cut features at the exact absolute cross cap
\(|p^\mathsf TCq|=W\).

| \(n\) | \(W_n\) | endpoint profile | internal dimension | exact-cap rank |
|---:|---:|---|---:|---:|
| 3 | 2 | \(1+2,\ |d|=1\) | 1 | 0 |
| 4 | 4 | \(2+2,\ |d|=0\) | 2 | 0 |
| 4 | 4 | \(2+2,\ |d|=2\) | 2 | 0 |
| 5 | 4 | \(1+4,\ d=0\) | 6 | 0 |
| 5 | 4 | \(2+3,\ d=0\) | 4 | 1 |
| 6 | 5 | \(1+5,\ d=0\) | 10 | 0 |
| 6 | 5 | \(3+3,\ d=0\) | 6 | 2 |
| 7 | 8 | \(2+5,\ |d|=1\) | 11 | 1 |
| 7 | 8 | \(3+4,\ |d|=1\) | 9 | 0 or 1 |
| 8 | 10 | \(4+4,\ d=0\) | 12 | 0 |

At \(n=8\), all \(67{,}200\) endpoint pairs have only the two trivial
absolute-cap cut representatives and rank zero.

Therefore exact extremality of the capped cut polytope is not merely
unproved; it is decisively absent on the known optimizers, including
centered ones.  A deficient exact face cannot by itself imply a
width-lowering replacement.

## 8. Revised frontier

The near-cap program now has a precise viable form:

1. use a *weighted* thick layer, not exact cap rank;
2. exploit the actual joint distribution of the margins
   \(W-|p^\mathsf TCq|\), not only a threshold \(T\);
3. construct a correlated internal signing \(g=B'\oplus D'\) obeying

   \[
   |g\cdot\delta(S)|
   \le
   \min\{c\cdot\delta(S),W-c\cdot\delta(S)\}
   \tag{8.1}
   \]

   for all cuts, while \(g\cdot\mathbf1=o(n^{3/2})\).

Such a \(g\) would make \(c+g\) another width-\(W_n\) signing with a
centered midpoint.  Random internal signs turn (8.1) into the exact
traffic sum

\[
\sum_{p,q}
\Pr\left\{
|S_{k(p,q)}|
>
\frac{W-|p^\mathsf TCq|}{2}
\right\},
\tag{8.2}
\]

where

\[
k(p,q)
=|\delta_U(p)|+|\delta_V(q)|
\]

is the number of active internal edges and \(S_k\) is a \(k\)-term
Rademacher sum.  If (8.2) is less than one (with room for the
high-probability event \(g\cdot\mathbf1=O(n)\)), the desired centered
minimizer exists.

The singleton/double rank theorem proves that no hidden linear
dimension is missing.  What remains is exactly the weighted
margin-versus-traffic inequality (8.2).
