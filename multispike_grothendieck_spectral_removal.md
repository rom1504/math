# Multi-spike Boolean rounding and Grothendieck--Pietsch spectral removal

## Normalization

Let \(A=A^\top\) have zero diagonal and \(|a_{ij}|\le1\).  Write

\[
Q(A)=\max_{x\in\{\pm1\}^n}|x^\top A x|,
\qquad
B(A)=\|A\|_{\infty\to1}
=\max_{x,y\in\{\pm1\}^n}|y^\top A x|.
\]

Coordinatewise affinity on the box and polarization give

\[
\max_{z\in[-1,1]^n}|z^\top Az|=Q(A),
\qquad
B(A)\le2Q(A).
\tag{1}
\]

Also \(\|A\|_{1\to\infty}=\max_{ij}|a_{ij}|\le1\), so interpolation
gives the useful exact estimate

\[
\boxed{\|A\|_{\mathrm{op}}^2\le B(A)\le2Q(A).}
\tag{2}
\]

Thus \(Q(A)=O(n^{3/2})\) already implies
\(\|A\|_{\mathrm{op}}=O(n^{3/4})\), improving the earlier
\(O(n^{5/6})\) bootstrap.

## 1. Weighted multi-spike inequality

Let \((\lambda_r,u_r)\), \(r\in\mathcal R\), be orthonormal eigenpairs
whose eigenvalues all have the same sign.  For arbitrary weights \(w_r\ge0\),
put

\[
L(w)=\max_i\sum_{r\in\mathcal R}w_r u_r(i)^2.
\]

Then

\[
\boxed{
Q(A)\ge
\frac1{2K_G}\,
\frac{\sum_{r\in\mathcal R}w_r|\lambda_r|}
     {\max_i\sum_{r\in\mathcal R}w_r u_r(i)^2}.
}
\tag{3}
\]

Here \(K_G\) is the real Grothendieck constant.

Indeed, the vectors

\[
p_i=\frac{(\sqrt{w_r}\,u_r(i))_{r\in\mathcal R}}{\sqrt{L(w)}}
\]

have norm at most one, and their vector-valued bilinear objective is

\[
\sum_{ij}a_{ij}\langle p_i,p_j\rangle
=\frac1{L(w)}\sum_r w_r\lambda_r.
\]

Grothendieck's inequality bounds this by \(K_GB(A)\), and (1) completes
the proof.

Two useful choices are

\[
w_r=\|u_r\|_\infty^{-2}
\quad\text{and}\quad
w_r=\|u_r\|_1^2.
\tag{4}
\]

For disjoint coherent blocks of sizes \(s_r\), the first choice has
\(L(w)=1\) and (3) gives, up to the Grothendieck constant,
\(\sum_r|\lambda_r|s_r\).  This is the correct scale: concatenating the
block sign vectors produces exactly that Boolean energy.  Thus (3)
genuinely aggregates disjoint localized spikes, unlike the one-eigenvector
bound

\[
Q(A)\ge |\lambda_r|/\|u_r\|_\infty^2.
\]

For a uniformly diffuse high-eigenvalue subspace,
\(\|u_r\|_1^2\asymp n\) and
\(\sum_r u_r(i)^2\asymp r/n\).  Equation (3) then forces
\(|\lambda_r|=O(Q/n)=O(\sqrt n)\).  Eigenvalues much larger than
\(\sqrt n\) therefore cannot occupy a uniformly diffuse subspace.

## 2. Exact Grothendieck--Pietsch removal theorem

The factorization form of Grothendieck's theorem gives probability
vectors \(\mu,\nu\) such that

\[
|x^\top Ay|
\le K_G B(A)
\left(\sum_i\mu_i x_i^2\right)^{1/2}
\left(\sum_j\nu_j y_j^2\right)^{1/2}
\quad\text{for all }x,y\in\mathbb R^n.
\tag{5}
\]

Fix \(0<\varepsilon<1\), and delete

\[
T=
\left\{i:\mu_i>\frac2{\varepsilon n}\right\}
\cup
\left\{i:\nu_i>\frac2{\varepsilon n}\right\}.
\]

Each set in the union has fewer than \(\varepsilon n/2\) members, so
\(|T|<\varepsilon n\).  On \(R=[n]\setminus T\), (5) gives

\[
\boxed{
\|A_{R,R}\|_{\mathrm{op}}
\le \frac{2K_GB(A)}{\varepsilon n}
\le \frac{4K_GQ(A)}{\varepsilon n}.
}
\tag{6}
\]

In particular, if \(Q(A)=O(n^{3/2})\), then deleting \(O(n/t)\) vertices
leaves operator norm \(O(t\sqrt n)\).  Taking \(t\to\infty\) proves a
precise spectral-localization statement:

> all spectral mass above a diverging multiple of \(\sqrt n\) can be
> removed on \(o(n)\) exceptional vertices.

This statement holds for every dense signing and does not require a
conference or pseudorandom hypothesis.

An equivalent symmetric proof uses the SDP dual.  The positive and
negative vector SDPs yield diagonal matrices \(D_+,D_-\succeq0\) with

\[
D_+-A\succeq0,\quad D_-+A\succeq0,
\qquad
\operatorname{tr}(D_++D_-)\le4K_GQ(A).
\]

Thus

\[
-(D_++D_-)\preceq A\preceq D_++D_-,
\tag{7}
\]

which is a one-measure quadratic version of (5), with a comparable
constant.

## 3. Selector consequence and the remaining obstruction

For a Bernoulli selector \(\delta_i\) of density \(\alpha\), a fixed
Boolean vector \(x\), and \(B=A_{R,R}\), write \(\xi_i=\delta_i-\alpha\).
The exact centered decomposition is

\[
x^\top D_\delta B D_\delta x-\alpha^2x^\top Bx
=
2\alpha\,\xi^\top D_xBx
+\xi^\top D_xBD_x\xi.
\tag{8}
\]

Hanson--Wright for a fixed \(x\) at deviation \(s\asymp n^{3/2}\)
has exponent

\[
\min\left\{
\Theta\!\left(\frac{s^2}{\|B\|_F^2}\right),
\Theta\!\left(\frac{s}{\|B\|_{\mathrm{op}}}\right)
\right\}.
\tag{9}
\]

After deleting \(\varepsilon n\) vertices, (6) makes the second exponent
only \(\Theta(\varepsilon n)\).  A direct union bound over the
\(2^{\alpha n}\) Boolean vectors therefore requires \(\varepsilon\)
bounded away from zero.  It cannot simultaneously use the
\(o(n)\)-vertex deletion needed for scale transfer.

Consequently (6) closes the spectral-localization half of the proposed
restriction argument, but not the selector half.  A successful completion
needs a chaining or entropy bound exploiting the correlation among the
\(2^{\alpha n}\) selector processes; treating them as independent events
loses exactly the gain supplied by deleting only \(o(n)\) vertices.
