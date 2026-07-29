# Orientation-even cancellation via Grothendieck factorization

Let \(A\) be a real symmetric zero-diagonal \(n\times n\) matrix with
\(|a_{ij}|\le 1\), and put
\[
Q(A)=\max_{x\in\{\pm1\}^n}|x^\top A x|.
\]

## 1. Mixed-sign spectral modes: exact vector witness

Let
\[
Au_r=\lambda_r u_r
\]
be any orthonormal collection of eigenpairs, of either sign, and let
\(w_r\ge 0\). Define
\[
D=\max_i\sum_r w_r u_r(i)^2.
\]
If \(D>0\), set
\[
p_i=D^{-1/2}(\sqrt{w_r}u_r(i))_r,\qquad
q_i=D^{-1/2}(\operatorname{sgn}\lambda_r\sqrt{w_r}u_r(i))_r.
\]
Then \(\|p_i\|_2,\|q_i\|_2\le1\), and
\[
\sum_{i,j}a_{ij}\langle p_i,q_j\rangle
=\frac1D\sum_r w_r|\lambda_r|.
\]
The real Grothendieck inequality therefore gives
\[
\frac1D\sum_r w_r|\lambda_r|
\le K_G\|A\|_{\infty\to1}.
\]
The cube-polarization identity gives
\[
\|A\|_{\infty\to1}\le2Q(A),
\]
so
\[
\boxed{
Q(A)\ge
\frac{\sum_r w_r|\lambda_r|}
{2K_G\max_i\sum_r w_r u_r(i)^2}.}
\tag{1}
\]

There is no same-sign hypothesis.  The spectral signs are placed in the
second vector family \(q_i\), so positive/negative cancellation disappears
exactly.  Equation (1) is an orientation-even, multi-mode version of the
single-eigenvector localization argument.

## 2. Grothendieck--Pietsch common-support deletion

Put \(B=\|A\|_{\infty\to1}\).  The Grothendieck factorization theorem
provides probability vectors \(\mu,\nu\) such that, for all real \(u,v\),
\[
|u^\top A v|
\le K_GB
\left(\sum_i\mu_i u_i^2\right)^{1/2}
\left(\sum_j\nu_j v_j^2\right)^{1/2}.
\tag{2}
\]
For \(0<\varepsilon\le1\), delete
\[
T=\{i:\mu_i>2/(\varepsilon n)\}
\cup
\{i:\nu_i>2/(\varepsilon n)\}.
\]
Since both weights have total mass one,
\[
|T|<\varepsilon n.
\]
On \(R=[n]\setminus T\), (2) implies
\[
\|A[R]\|_{\rm op}
\le\frac{2K_GB}{\varepsilon n}
\le\frac{4K_GQ(A)}{\varepsilon n}.
\tag{3}
\]

Equivalently, for every \(L>0\) one may take
\[
\varepsilon=\frac{4K_GQ(A)}{nL}
\]
whenever this is at most one, and obtain
\[
\boxed{
|T|<\frac{4K_GQ(A)}L,\qquad
\|A[[n]\setminus T]\|_{\rm op}\le L.}
\tag{4}
\]
Thus if \(Q(A)=O(n^{3/2})\), every threshold \(L/\sqrt n\to\infty\)
admits a common deletion set \(T=o(n)\) which removes **all**
positive and negative spectral modes above \(L\).

This proves the common-support theorem sought in the spectral-localization
route.  It is stronger and cleaner than unioning localization sets of
individual eigenvectors.

## 3. What this does and does not yet settle

At fixed deletion density \(\varepsilon\), (3) gives
\[
\|A[R]\|_{\rm op}=O_\varepsilon(\sqrt n)
\]
for every near minimizer.  If \(\varepsilon\to0\), however, the coefficient
is \(O(1/\varepsilon)\).  In a Hanson--Wright selector bound this produces
an exponent of order \(\varepsilon^2n\), which does not uniformly pay for
the \(2^n\) possible Boolean witnesses.  Therefore (4) solves the
positive/negative cancellation and common-support problem, but does not by
itself prove proportional restriction or convergence of
\(M_n/n^{3/2}\).

