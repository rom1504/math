# Induced hafnian anti-cancellation audit

## Setup

For a symmetric zero-diagonal signing \(A=(a_{ij})\), define

\[
h_k(A)=\sum_{\substack{S\subset[n]\\ |S|=2k}}\operatorname{haf}(A_S)^2,
\qquad
\operatorname{haf}(A_S)=
\sum_{M\in\mathsf{PM}(S)}\prod_{e\in M}a_e .
\]

This is the squared norm of the top Fourier layer of \(H_A^k\):

\[
\|\operatorname{Proj}_{2k}H_A^k\|_2^2=(k!)^2h_k(A).
\]

## 1. Exact pair-of-matchings expansion

An ordered pair \((M,N)\) of perfect matchings on the same vertex set
has a union whose connected components are:

* a doubled edge, when an edge belongs to both matchings; or
* a simple alternating cycle of even length at least four.

A doubled edge has weight \(a_e^2=1\).  A simple even cycle \(C\)
has weight \(\chi_A(C)=\prod_{e\in C}a_e\), and its two alternating
two-colourings correspond to the two possible ordered pairs of matchings.
Consequently

\[
h_k(A)=
\sum_{\substack{\mathcal F:\text{ vertex-disjoint}\\
                 \text{even cycles covering }2k\text{ vertices}}}
\prod_{\substack{C\in\mathcal F\\|C|=2}}1
\prod_{\substack{C\in\mathcal F\\|C|\ge4}}2\chi_A(C).
\tag{1}
\]

This formula makes the cancellation mechanism exact: cancellation takes
place between even-cycle packings.

Let \(z_i\) be commuting square-zero variables.  In the algebra
\(\mathbb R[z_1,\ldots,z_n]/(z_1^2,\ldots,z_n^2)\),

\[
\boxed{
\sum_{S\subset[n]}\operatorname{haf}(A_S)^2z_S
=
\exp\left[
\sum_{\substack{m\ge2\\m\ {\rm even}}}
\frac1m
\sum_{\substack{i_1,\ldots,i_m\\\text{all distinct}}}
\left(\prod_{r=1}^m a_{i_ri_{r+1}}\right)
\prod_{r=1}^m z_{i_r}
\right],
}
\tag{2}
\]

where \(i_{m+1}=i_1\).  Equivalently this is the square-free part of

\[
\boxed{
\det\!\left(I-(ZA)^2\right)^{-1/2},\qquad
Z=\operatorname{diag}(z_1,\ldots,z_n).
}
\tag{3}
\]

Indeed, the logarithm in (3) is
\(\sum_{m\ {\rm even}}\operatorname{tr}(ZA)^m/m\); a simple undirected
cycle of length \(m\ge4\) occurs \(2m\) times in the trace and hence gets
weight \(2\), while the length-two term gives one copy of each doubled
edge.  Square-free projection deletes every repeated-vertex walk.

## 2. Exact \(k=2\) formula and spectral minimum

For a four-set \(S\),

\[
\operatorname{haf}(A_S)^2
=3+2\sum_{\text{three 4-cycles on }S}\chi_A(C).
\]

Separating simple cycles from the repeated-vertex walks in
\(\operatorname{tr}A^4\) gives

\[
\boxed{
h_2(A)
=3\binom n4+
\frac14\left[
\operatorname{tr}A^4-n(n-1)(2n-3)
\right].
}
\tag{4}
\]

Since \(\operatorname{tr}A^2=n(n-1)\),
Cauchy--Schwarz yields
\(\operatorname{tr}A^4\ge n(n-1)^2\).  Therefore

\[
h_2(A)\ge
\frac{n(n-1)(n-2)(n-5)}8.
\tag{5}
\]

Equality in the spectral step occurs for a conference matrix
\(A^2=(n-1)I\).  At order six this gives \(h_2=15=\binom64\), the
pointwise parity floor: every four-vertex principal hafnian is \(\pm1\).
For large \(n\), however, (5) is only a \(1-O(1/n)\) loss from the random
sign benchmark \(3\binom n4\).

## 3. Exact full-order counterfamily

Let \(B_k\) have vertices
\[
L_1,R_1,\ldots,L_k,R_k.
\]
Every off-diagonal entry is \(+1\), except
\[
B_k(R_r,L_s)=B_k(L_s,R_r)=-1\qquad(r<s).
\tag{6}
\]

Then

\[
\boxed{\operatorname{haf}(B_k)=1\quad\text{for every }k\ge1.}
\tag{7}
\]

Proof.  Expand the hafnian at the last two vertices.  Their mutual match
contributes \(\operatorname{haf}(B_{k-1})\).  On the old vertices,
\(R_k\) has all \(+1\) incidences, while \(L_k\) is \(+1\) on old
\(L\)-vertices and \(-1\) on old \(R\)-vertices.  Hence all cross matches
contribute

\[
2\sum_{i<j}
\operatorname{haf}\!\left(B_{k-1}\setminus\{L_i,L_j\}\right)
-
2\sum_{i<j}
\operatorname{haf}\!\left(B_{k-1}\setminus\{R_i,R_j\}\right).
\tag{8}
\]

The map \(L_i\leftrightarrow R_{k-i}\) (indices \(1,\ldots,k-1\))
is an automorphism of \(B_{k-1}\), so the two sums agree.  Thus
\(\operatorname{haf}(B_k)=\operatorname{haf}(B_{k-1})\), and the base
case is \(1\).

Exact dynamic programming independently verified (7) for \(k\le10\).

## 4. Consequence

At the endpoint \(n=2k\),

\[
h_k(B_k)=1,
\]

whereas the random-sign mean is

\[
\mathbb E_A h_k(A)=(2k-1)!!.
\]

Therefore no universal estimate of the form

\[
h_k(A)\ge (2k-1)!!\,e^{-Ck}
\]

can hold with an absolute \(C\) through \(\alpha=k/n=1/2\):
the ratio \(1/(2k-1)!!\) is \(e^{-k\log k+O(k)}\).
The proposed universal RMS perfect-matching scale is false.

For fixed \(\alpha<1/2\), the construction does not by itself settle the
question.  Exact values up to \(n=20\) show that its induced-minor sum is
often at or above the random benchmark in the bulk and only collapses
near full order.  Thus a theorem restricted to
\(\alpha\le1/2-\varepsilon\) remains logically possible, but it would
still need many Fourier boundary layers to compete with the current
\(0.33649\) lower bound.
