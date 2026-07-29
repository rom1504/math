# Regular microblock lifts: exact covariance and the independent-absorption wall

## 1. Setup

Let \(A=(a_{ij})_{i,j\le n}\) be symmetric, have zero diagonal, and
\(a_{ij}\in\{\pm1\}\) for \(i\ne j\). Let \(H\in\{\pm1\}^{r\times r}\)
be a regular Hadamard matrix:

\[
HH^\top=rI,\qquad H\mathbf 1=\sqrt r\,\mathbf 1,\qquad
H^\top\mathbf 1=\sqrt r\,\mathbf 1.
\]

For each unordered macro edge \(i<j\), choose independent uniform
permutation matrices \(P_{ij},Q_{ij}\), and put

\[
L_{ij}=a_{ij}P_{ij}HQ_{ij},\qquad L_{ji}=L_{ij}^{\top}.
\]

The diagonal \(r\times r\) blocks are omitted here; they may be filled
separately and contribute at most \(nQ(D_r)\) to the full quadratic
norm.

For a fiber spin \(x_i\in\{\pm1\}^r\), write

\[
m_i=\frac1r\mathbf 1^\top x_i,\qquad
u_i=x_i-m_i\mathbf 1,\qquad u_i\perp\mathbf 1.
\]

The half-sum Hamiltonian of the off-diagonal lift is

\[
\mathcal H_L(x)=\sum_{i<j}x_i^\top L_{ij}x_j.
\]

## 2. Exact macro/residual decomposition

Because permutations preserve \(\mathbf 1\) and a regular Hadamard
maps \(\mathbf 1^\perp\) into itself, the two cross terms vanish
identically. Hence, edge by edge,

\[
x_i^\top L_{ij}x_j
=a_{ij}r^{3/2}m_im_j+Z_{ij},
\]

where

\[
Z_{ij}=a_{ij}u_i^\top P_{ij}HQ_{ij}u_j.
\]

Consequently,

\[
\boxed{
\mathcal H_L(x)=r^{3/2}H_A(m)+R(x)
}
\]

with

\[
H_A(m)=\sum_{i<j}a_{ij}m_im_j,\qquad
R(x)=\sum_{i<j}Z_{ij}.
\]

The macro channel is therefore exactly the soft-spin quadratic form
of the seed. Random orientations affect only the centered residual.

## 3. Exact finite-\(r\) covariance theorem

For two fiber configurations \(x^a,x^b\), define their centered
fiber overlaps

\[
c_i^{ab}
=\frac1r\langle u_i^a,u_i^b\rangle
=\frac1r\langle x_i^a,x_i^b\rangle-m_i^am_i^b.
\]

Let

\[
P_0=I-\frac1rJ
\]

be the orthogonal projection onto \(\mathbf 1^\perp\). For
\(u,v\perp\mathbf 1\), uniform permutation averaging gives

\[
\mathbb E_P[P^\top uv^\top P]
=\frac{\langle u,v\rangle}{r-1}P_0.
\]

Regularity of \(H\) gives

\[
H^\top P_0H=rP_0,\qquad
\operatorname{Tr}(H^\top P_0HP_0)=r(r-1).
\]

Applying the permutation identity first to \(P_{ij}\), then to
\(Q_{ij}\), yields

\[
\mathbb E Z_{ij}^a=0
\]

and the exact covariance formula

\[
\boxed{
\mathbb E[Z_{ij}^aZ_{ij}^b]
=\frac{r^3}{r-1}c_i^{ab}c_j^{ab}.
}
\tag{3.1}
\]

In particular,

\[
\boxed{
\operatorname{Var}(Z_{ij})
=\frac{r^3}{r-1}(1-m_i^2)(1-m_j^2).
}
\tag{3.2}
\]

Independence across macro edges therefore gives

\[
\boxed{
\operatorname{Cov}(R^a,R^b)
=\frac{r^3}{r-1}\sum_{i<j}c_i^{ab}c_j^{ab}.
}
\tag{3.3}
\]

The covariance is independent of the seed signs \(a_{ij}\) and,
beyond regularity, independent of the chosen Hadamard matrix.

For fixed \(r\) and fixed \(k\), the vector
\((R(x^1),\ldots,R(x^k))\), after normalization by its standard
deviation, satisfies the usual triangular-array central limit theorem
as \(n\to\infty\), provided no single macro edge carries a nonvanishing
fraction of the total covariance. If \(r\) grows as well, the same
conclusion requires the corresponding Lindeberg condition.

Thus every independent regular-Hadamard microblock lift appends the
same Gaussian \(k\)-profile, determined only by the centered overlap
array \((c_i^{ab})\). It cannot encode higher seed structure in the
residual channel.

## 4. The independent-absorption obstruction

Write

\[
Q(A)=\max_{x\in\{\pm1\}^n}|x^\top Ax|
=2\max_x|H_A(x)|.
\]

An independent-noise absorption proof with limiting half-sum constant
\(c<1/2\) would need the macro soft-spin energy to contract in
proportion to its retained Boolean mass. The natural necessary
inequality is

\[
|H_A(m)|
\le
\max_{x\in\{\pm1\}^n}|H_A(x)|
\frac{\|m\|_2^2}{n},
\qquad m\in[-1,1]^n.
\tag{4.1}
\]

In full quadratic normalization this says

\[
|m^\top Am|
\le Q(A)\frac{\|m\|_2^2}{n}.
\tag{4.2}
\]

Take a top absolute eigenvector \(z\), scaled so that
\(\|z\|_\infty\le1\). Then (4.2) forces

\[
\|A\|_{\mathrm{op}}\|z\|_2^2
\le Q(A)\frac{\|z\|_2^2}{n},
\]

or

\[
\boxed{n\|A\|_{\mathrm{op}}\le Q(A).}
\tag{4.3}
\]

The reverse inequality \(Q(A)\le n\|A\|_{\mathrm{op}}\) always holds.
Therefore (4.1) can hold in every direction only at exact spectral
saturation, when a Boolean vector is an extremal eigenvector.

For conference-like seeds,

\[
\|A\|_{\mathrm{op}}\sim\sqrt n,
\]

so (4.3) requires \(Q(A)\sim n^{3/2}\), corresponding to the
half-sum spectral constant \(1/2\). Hence:

> **Independent-absorption wall.** A microblock proof that keeps the
> macro channel unchanged and appends a seed-independent centered
> residual profile cannot prove an amplification theorem below the
> spectral constant \(1/2\). The obstruction already appears in an
> infinitesimal soft-spin eigenvector direction.

This is a methodological obstruction, not a lower bound on the true
minimization problem. To beat it, the microblock residual must be
dependently coupled to the macro local fields so that its conditional
mean or covariance supplies a compensating Onsager/cavity term.

## 5. Target for a dependent correction

A useful dependent lift would have to replace the exact decomposition

\[
r^{3/2}H_A(m)+R
\]

by a uniform finite-\(r\) inequality of the form

\[
\left|r^{3/2}H_A(m)+R_A(m,u)\right|
\le
r^{3/2}\frac{Q(A)}2
+o\!\left((nr)^{3/2}\right)
\]

simultaneously for every Boolean fiber configuration. Since the
failure of (4.2) is controlled by

\[
m^\top\!\left(A-\frac{Q(A)}n I\right)m,
\]

the first candidate correction must be correlated with the seed local
field \(Am\), rather than with centered overlaps alone.

## 6. General block lifts and the exact microcanonical ANOVA identity

The preceding obstruction can be sharpened without assuming
independence, a Hadamard template, or regular blocks.

Let \(B_{ij}\in\{\pm1\}^{r\times r}\) be arbitrary deterministic
off-diagonal blocks, with \(B_{ji}=B_{ij}^{\top}\). Put

\[
e_0=\frac{\mathbf 1}{\sqrt r},\qquad P_0=I-e_0e_0^\top,
\]

and decompose each block orthogonally as

\[
B_{ij}
=s_{ij}e_0e_0^\top
+q_{ij}e_0^\top
+e_0c_{ij}^\top
+R_{ij},
\tag{6.1}
\]

where

\[
\begin{aligned}
s_{ij}&=e_0^\top B_{ij}e_0,\\
q_{ij}&=P_0B_{ij}e_0,\\
c_{ij}&=P_0B_{ij}^{\top}e_0,\\
R_{ij}&=P_0B_{ij}P_0.
\end{aligned}
\]

The four terms in (6.1) are Frobenius-orthogonal, so the sign
constraint gives the exact budget

\[
\boxed{
r^2=s_{ij}^2+\|q_{ij}\|_2^2+\|c_{ij}\|_2^2
+\|R_{ij}\|_F^2.
}
\tag{6.2}
\]

Fix a grid magnetization

\[
m_i\in\mathcal G_r
=\left\{-1,-1+\frac2r,\ldots,1-\frac2r,1\right\}.
\]

Independently for each \(i\), let \(X_i\) be uniform on the Hamming
slice

\[
\{x\in\{\pm1\}^r:\mathbf 1^\top x=rm_i\},
\]

and set \(U_i=X_i-m_i\mathbf 1\). Then

\[
\mathbb E U_i=0,\qquad
\mathbb E U_iU_i^\top
=\alpha_iP_0,\qquad
\alpha_i=\frac r{r-1}(1-m_i^2).
\tag{6.3}
\]

For the multipartite Hamiltonian

\[
\mathcal H_B(X)=\sum_{i<j}X_i^\top B_{ij}X_j,
\]

define the constant channel

\[
\mu_B(m)=r\sum_{i<j}s_{ij}m_im_j
\tag{6.4}
\]

and the centered one-fiber field

\[
g_i(m)
=\sqrt r\sum_{j\ne i}m_jq_{ij}^{(i)}.
\tag{6.5}
\]

Here \(q_{ij}^{(i)}=P_0B_{ij}e_0\) when \(i<j\), and
\(q_{ij}^{(i)}=P_0B_{ji}^{\top}e_0\) when \(j<i\).

Expanding \(X_i=\sqrt r\,m_ie_0+U_i\) gives the exact Hoeffding
decomposition

\[
\mathcal H_B(X)
=\mu_B(m)
+\sum_i U_i^\top g_i(m)
+\sum_{i<j}U_i^\top R_{ij}U_j.
\tag{6.6}
\]

All nonconstant terms displayed in (6.6) are mutually orthogonal in
\(L^2\) of the product of Hamming slices. Therefore:

\[
\boxed{
\begin{aligned}
\mathbb E\mathcal H_B(X)^2
={}&\mu_B(m)^2\\
&+\sum_i\alpha_i\|g_i(m)\|_2^2\\
&+\sum_{i<j}\alpha_i\alpha_j\|R_{ij}\|_F^2.
\end{aligned}
}
\tag{6.7}
\]

Since the maximum of \(|\mathcal H_B|\) dominates its root mean
square, (6.7) yields the finite-\(r\) energy inequality

\[
\boxed{
\begin{aligned}
\max_X|\mathcal H_B(X)|
\ge
\sup_{m\in\mathcal G_r^n}
\Bigg[
&\left(r\sum_{i<j}s_{ij}m_im_j\right)^2\\
&+\frac r{r-1}\sum_i(1-m_i^2)\|g_i(m)\|_2^2\\
&+\left(\frac r{r-1}\right)^2
\sum_{i<j}(1-m_i^2)(1-m_j^2)\|R_{ij}\|_F^2
\Bigg]^{1/2}.
\end{aligned}
}
\tag{6.8}
\]

If diagonal fiber blocks \(D_i\) are added, the norm of their sum is
at most \(\sum_i M(D_i)\). Thus the right side of (6.8), minus
\(\sum_iM(D_i)\), remains a lower bound for the full Hamiltonian norm.
For \(M(D_i)=O(r^{3/2})\), this is an
\(O(n^{-1/2})\) normalized error on the \((nr)^{3/2}\) scale.

For regular Hadamard blocks,

\[
s_{ij}=a_{ij}\sqrt r,\qquad q_{ij}=c_{ij}=0,\qquad
\|R_{ij}\|_F^2=r^2-r,
\]

and (6.7) reduces exactly to the covariance theorem in Section 3.

## 7. Why a static Onsager field has the wrong sign

The most direct field-correlated ansatz is to choose row-sum
deviations of the form

\[
q_{ij}^{(i)}=\gamma a_{ij}v_i
\]

for centered microvectors \(v_i\). Equation (6.5) then gives

\[
g_i(m)=\gamma\sqrt r\,(Am)_i v_i.
\]

But its complete contribution to the exact conditional second moment
is

\[
\boxed{
\frac{\gamma^2r^2}{r-1}
\sum_i(1-m_i^2)(Am)_i^2\|v_i\|_2^2\ge0.
}
\tag{7.1}
\]

It has no negative covariance with the macro term: the macro term is
the degree-zero Hoeffding component, while the local-field correction
is degree one. They are orthogonal exactly, not just asymptotically.

More generally, correlations among the block orientations cannot
alter (6.7), because (6.7) averages over fiber spins after the blocks
have been fixed. Incident-edge row effects may cancel inside
\(g_i(m)\), but the best possible cancellation merely removes the
positive middle line of (6.7). It can never compensate the macro
energy.

The Frobenius budget (6.2) makes the remaining tradeoff explicit:

- regularity sets the field channels \(q_{ij},c_{ij}\) to zero, but
  leaves essentially all microscopic mass in the irreducible
  two-fiber channel \(R_{ij}\);
- moving mass out of \(R_{ij}\) creates row or column fields, which
  enter (6.7) as positive squares;
- making \(g_i(m)=0\) for every \(m\) forces every incident
  coefficient \(q_{ij}^{(i)}=0\), returning to the regular case.

Consequently:

> **One-step Onsager no-go.** A static cavity correction implemented
> by centered row/column biases of sign blocks cannot yield the desired
> amplification inequality. At fixed macro magnetization its exact
> effect is a nonnegative variance term. Any successful dependent
> construction must instead alter the degree-zero coarse block means,
> or use a genuinely multistep/nonlinear mechanism not representable
> as a one-step local-field orientation bias.

This no-go does not exclude an \(A\)-dependent choice of coarse means
\(s_{ij}\), for example an odd spectral filter involving \(A^3\).
Such a filter changes the quotient itself, however, and must pay for
the microscopic Frobenius mass displaced by (6.2). That is the next
remaining finite-type possibility.
