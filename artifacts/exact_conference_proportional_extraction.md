# Exact-conference proportional extraction

Let \(C\) be a symmetric conference matrix of order \(N\):
\[
C^2=(N-1)I,\qquad c_{ii}=0,\quad c_{ij}\in\{\pm1\}.
\]
For \(S\subset[N]\), \(|S|=m\), write
\[
A=C[S,S],\qquad B=C[S,S^c].
\]
The desired universal extraction statement is
\[
\max_{x\in\{\pm1\}^m}|x^\top A x|
\le (1+o(1))m^{3/2}.
\tag{1}
\]

No proof of (1) for all conference matrices was obtained.  This note records
the exact reductions and the barriers found.

## 1. Sparse-cube formulation

Put \(U=C/\sqrt{N-1}\), a symmetric orthogonal involution.  A pair
\((S,x)\) is the same as a ternary vector
\[
z\in\{0,\pm1\}^N,\qquad |\operatorname{supp}z|=m.
\]
Then
\[
x^\top A x=z^\top Cz
=\sqrt{N-1}\,\langle z,Uz\rangle.
\]
Writing \(m/N\to\alpha\), (1) asks for one coordinate \(m\)-plane whose
entire Boolean cube avoids the two cones
\[
\boxed{
\frac{|\langle z,Uz\rangle|}{\|z\|_2^2}
>\sqrt\alpha+o(1).}
\tag{2}
\]
This is a coordinate-subspace avoidance problem, not an ordinary spectral
compression problem.

## 2. Exact cross-Gram identity and joint bad-witness condition

The \(S\times S\) block of \(C^2=(N-1)I\) gives
\[
\boxed{A^2+BB^\top=(N-1)I_m.}
\tag{3}
\]
For every Boolean \(x\),
\[
\|B^\top x\|_2^2=(N-1)m-\|Ax\|_2^2
\tag{4}
\]
and Cauchy--Schwarz gives
\[
(x^\top Ax)^2\le m\|Ax\|_2^2.
\tag{5}
\]
Consequently, if
\[
|x^\top Ax|>(1+\eta)m^{3/2},
\]
then
\[
\boxed{
\|B^\top x\|_2^2
<
(N-1)m-(1+\eta)^2m^2.}
\tag{6}
\]
Since the sign-average of the left side is
\[
\mathbb E_x\|B^\top x\|_2^2=\|B\|_F^2=m(N-m),
\]
a bad internal Rayleigh quotient forces an anomalously small cross norm.
Thus the exact missing object is the **joint** profile
\[
\left(x^\top Ax,\ \|B^\top x\|_2^2\right).
\]

Equivalently, define the even defect
\[
D_S=A^2-(m-1)I.
\]
Then
\[
\boxed{
(x^\top Ax)^2
\le m^2(m-1)+m\,x^\top D_Sx.}
\tag{7}
\]
Conference orthogonality identifies, for \(i\ne j\in S\),
\[
(D_S)_{ij}
=(A^2)_{ij}
=-\sum_{k\notin S}c_{ik}c_{kj}.
\tag{8}
\]
The strong condition
\[
\max_xx^\top D_Sx=o(N^2)
\tag{9}
\]
would prove (1).  But (9) is too strong: for a random proportional support,
the entries in (8) naturally have size \(\sqrt N\), and their Boolean
quadratic norm is naturally of order \(N^2\).  A successful proof must use
the fact that the vector maximizing the defect also has to have a large
Rayleigh quotient for \(A\).

## 3. Exact field-square simplification

For a full Boolean vector \(x\), define the switched fields
\[
r_i=x_i(Cx)_i.
\]
Conference orthogonality gives the exact identity
\[
\boxed{
\sum_i r_i^2=\|Cx\|_2^2=N(N-1),}
\tag{10}
\]
independently of \(x\).  Therefore the random-principal-restriction
second-moment formula has no uncontrolled field-profile term for a
conference matrix.  This removes the principal variance obstruction, but a
second moment cannot control the maximum over \(2^m\) signs.

## 4. Why direct bad-pair counting is too weak

There are \(\binom Nm2^m\) ternary vectors of support \(m\).  To prove that
some support is good by counting bad pairs \((S,x)\), the bad probability
for a uniformly random ternary vector must be smaller than \(2^{-m}\).

Even in the ideal Gaussian/Haar flat-involution model, the correlation tail
at level \(t\) has rate
\[
I_{\rm flat}(t)=\frac14\log\frac1{1-t^2}.
\]
At the target \(t=\sqrt\alpha\), the exponent is
\[
NI_{\rm flat}(\sqrt\alpha)
=\frac N4\log\frac1{1-\alpha}.
\tag{11}
\]
Paying for all \(2^m\) sign vectors requires
\[
\frac14\log\frac1{1-\alpha}>\alpha\log2.
\tag{12}
\]
The nonzero root is
\[
\alpha_0=0.9225232669\ldots.
\]
Thus even the ideal flat rate is insufficient for a first-moment proof at
\(\alpha=1/2\) and throughout most of \((0,1)\).  Moreover, arithmetic
conference resonances prevent using the flat tail uniformly without an
additional argument.

Support-choice entropy must therefore be used at the level of the
**lower tail of the maximum within a support**, not by union-bounding its
\(2^m\) corners.

## 5. Spectral/interlacing wall

Principal compression of a balanced orthogonal involution has the
free-compression edge
\[
2\sqrt{\alpha(1-\alpha)}
\]
for \(\alpha\le1/2\), and an edge at \(1\) once the coordinate subspace has
forced intersection with an eigenspace.  Characteristic-polynomial
interlacing sees this spectral edge, not Boolean corners.  After restoring
the \(\sqrt N\) scale it is generally larger than \(\sqrt m\), so
\[
\frac m2\|C[S]\|_{\rm op}
\]
cannot yield the target \(\frac12m^{3/2}\).

## 6. Verdict

The conference task is reduced sharply to:

> Find a support \(S\) for which no Boolean vector has both a large internal
> Rayleigh quotient \(x^\top A x\) and the cross-norm deficit forced by
> (6).

Pure spectral paving, Grothendieck--Pietsch weights, a defect norm bound, and
first-moment counting of bad ternary vectors all discard essential joint
information and do not prove proportional extraction.

