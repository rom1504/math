# Grothendieck--Pietsch spectral regularization

## 1. Norm convention and preliminary inequality

Let \(A\) be symmetric with zero diagonal and define

\[
Q(A)=\max_{x\in\{\pm1\}^n}|x^\top Ax|.
\]

Thus \(Q(A)=2M(A)\) when
\(M(A)=\max_x|\sum_{i<j}a_{ij}x_ix_j|\).

Zero diagonal makes \(q(z)=z^\top Az\) affine in each coordinate separately,
so

\[
\sup_{\|z\|_\infty\le1}|z^\top Az|=Q(A).
\tag{1}
\]

For sign vectors \(u,v\), put

\[
a=\frac{u+v}{2},\qquad b=\frac{u-v}{2}.
\]

The vectors \(a,b\) lie in the cube and have disjoint supports. Symmetry gives

\[
u^\top Av=a^\top Aa-b^\top Ab.
\]

Consequently,

\[
\boxed{\|A\|_{\infty\to1}\le2Q(A).}
\tag{2}
\]

Also \(\|A\|_{1\to\infty}=\max_{ij}|a_{ij}|\le1\). Interpolation at
\(\theta=1/2\) therefore gives the independently audited global estimate

\[
\boxed{\|A\|_{\mathrm{op}}\le\sqrt{2Q(A)}.}
\tag{3}
\]

If \(Q(A)=O(n^{3/2})\), this is \(O(n^{3/4})\), improving the earlier
\(O(n^{5/6})\) bootstrap.

## 2. Large-core spectral removal theorem

Grothendieck--Pietsch factorization for the bilinear form
\((u,v)\mapsto u^\top Av\) supplies probability weights
\(\mu,\nu\) on \([n]\) such that

\[
|u^\top Av|
\le
K_G\|A\|_{\infty\to1}
\left(\sum_i\mu_i u_i^2\right)^{1/2}
\left(\sum_j\nu_j v_j^2\right)^{1/2}
\tag{4}
\]

for all real \(u,v\). Here \(K_G\) is the real Grothendieck constant.

For \(\epsilon\in(0,1)\), remove every index satisfying

\[
\mu_i>\frac{2}{\epsilon n}
\quad\text{or}\quad
\nu_i>\frac{2}{\epsilon n}.
\]

Each exceptional set has fewer than \(\epsilon n/2\) members. The remaining
set \(R\) therefore has

\[
|R|\ge(1-\epsilon)n.
\]

For vectors supported on \(R\), (4) gives

\[
|u^\top A[R]v|
\le
\frac{2K_G}{\epsilon n}
\|A\|_{\infty\to1}\|u\|_2\|v\|_2.
\]

Hence

\[
\boxed{
\|A[R]\|_{\mathrm{op}}
\le
\frac{2K_G}{\epsilon n}\|A\|_{\infty\to1}
\le
\frac{4K_G}{\epsilon n}Q(A).}
\tag{5}
\]

In particular:

\[
\boxed{
Q(A)\le Cn^{3/2}
\quad\Longrightarrow\quad
\exists\,|R|\ge(1-\epsilon)n:
\ \|A[R]\|_{\mathrm{op}}
\le\frac{4K_GC}{\epsilon}\sqrt n.}
\tag{6}
\]

The theorem does not require the entries to equal \(\pm1\); symmetry, zero
diagonal, and the bilinear norm bound suffice.

## 3. Audit against selector/union-bound constants

Write

\[
M(A)=cn^{3/2},\qquad Q(A)=2cn^{3/2}.
\]

Then (6) becomes

\[
\|A[R]\|_{\mathrm{op}}
\le\frac{8K_Gc}{\epsilon}\sqrt n.
\tag{7}
\]

Take the retained set itself to have target size \(m=\alpha n\), so
\(\epsilon=1-\alpha\), and seek

\[
M(A[R])\le L:=c\alpha^{3/2}n^{3/2}.
\]

For a uniform spin on \(R\), Hanson--Wright at quadratic threshold \(2L\)
has exponent, up to its universal constant \(c_{\rm HW}\),

\[
c_{\rm HW}
\min\left\{
\frac{(2L)^2}{\|A[R]\|_F^2},
\frac{2L}{\|A[R]\|_{\mathrm{op}}}
\right\}
\ge
c_{\rm HW}
\min\left\{
4c^2\alpha,\,
\frac{\alpha^{3/2}(1-\alpha)}{4K_G}
\right\}n.
\tag{8}
\]

To union-bound over \(2^m\) spins, this coefficient would have to exceed
\(\alpha\log2\). The spectral branch would require

\[
\frac{c_{\rm HW}}{4K_G}\sqrt\alpha(1-\alpha)>\log2.
\tag{9}
\]

But

\[
\max_{0<\alpha<1}\sqrt\alpha(1-\alpha)
=\frac{2}{3\sqrt3}=0.384900\ldots .
\]

Thus (9) fails even if one unrealistically replaces both \(c_{\rm HW}\) and
\(K_G\) by \(1\). This is a structural constant loss, not a matter of
optimizing the published Hanson--Wright constant.

Small multiplicative steps do not help. For \(\alpha=1-d\), the usable tail
rate from (8) is \(O(d)n\) (and in variance-dominated formulations
\(O(d^2)n\)), whereas the cube entropy remains
\((1-d)(\log2)n\).

## 4. Retain first, then select

Suppose instead that the regularized core has size
\(\rho n\), with \(\rho=1-\epsilon\), and a random \(m=\alpha n\) subset is
selected inside it. The conditional mean of a restricted energy can be as
large as

\[
\left(\frac{\alpha}{\rho}\right)^2M(A).
\]

For the target \(L=\alpha^{3/2}M(A)\), a positive selector gap requires

\[
\rho>\alpha^{1/4}.
\tag{10}
\]

With (7), the normalized gap is

\[
d(\alpha,\rho)
=\alpha^{3/2}-\frac{\alpha^2}{\rho^2},
\]

and the variance-dominated selector exponent is at best a constant multiple
of

\[
\epsilon^2d(\alpha,\rho)^2n.
\tag{11}
\]

Under \(\epsilon<1-\alpha^{1/4}\), this is far below the
\(\alpha(\log2)n\) entropy that a black-box union or raw layer count must
overcome. The multiplicity-weighted layer lemma therefore also does not close
from (6) alone.

## 5. Verdict

The spectral removal statement requested in this route is true, with the
explicit bound

\[
K(C,\epsilon)=\frac{4K_GC}{\epsilon}.
\]

It cleanly removes the \(\sqrt n\)-hub obstruction and converts the
operator-norm exponent to the correct \(\sqrt n\) scale on a
\((1-\epsilon)n\)-vertex core. Nevertheless, black-box spectral tails still
lose too much constant against the Boolean-cube entropy, and this remains
true under iteration.

The missing ingredient after factorization must use the geometry or entropy
of the indexed spin process, rather than union-bounding its \(2^m\) values.

