# Boolean spectral-profile compactification

## Purpose

For a symmetric signing \(A\) of order \(n\), put
\[
 U_A=\frac{A}{\sqrt n}.
\]
After the purification theorem, one may restrict at arbitrarily small
variational cost to sequences satisfying
\[
 \|U_A\|_{\mathrm{op}}\le C. \tag{1}
\]
Ordinary graphons lose the \(n^{-1/2}\)-scale matrix, and ordinary
probability/traffic limits can lose exponentially rare maximizing spin
vectors.  This note gives a compactification which retains the exact
Boolean quadratic supremum by definition.

It does **not** prove all-order realization.  It isolates that missing
statement without any continuity or uniform-integrability ambiguity.

## 1. One-spin Boolean spectral profile

Let \(E_A\) be the projection-valued spectral measure of \(U_A\).  For a
sign vector \(x\in\{\pm1\}^n\), define a probability measure on
\([-C,C]\) by
\[
 \mu_{A,x}(B)
 =\frac1n\langle x,E_A(B)x\rangle. \tag{2}
\]
Indeed, \(\mu_{A,x}([-C,C])=\|x\|_2^2/n=1\).

Define the compact finite set
\[
 \mathcal K_1(A)
 =\{\mu_{A,x}:x\in\{\pm1\}^n\}
 \subset\mathcal P([-C,C]). \tag{3}
\]
Then
\[
 \frac{x^\top A x}{n^{3/2}}
 =\frac1n x^\top U_Ax
 =\int t\,d\mu_{A,x}(t), \tag{4}
\]
and hence
\[
\boxed{
 \frac{1}{n^{3/2}}\max_{x\in\{\pm1\}^n}|x^\top Ax|
 =\sup_{\mu\in\mathcal K_1(A)}
 \left|\int t\,d\mu(t)\right|.} \tag{5}
\]

The hyperspace of nonempty compact subsets of
\(\mathcal P([-C,C])\), equipped with Hausdorff distance induced by
\(W_1\), is compact.  Therefore every sequence satisfying (1) has a
subsequence for which
\[
 \mathcal K_1(A_n)\longrightarrow\mathcal K_1^\infty. \tag{6}
\]

Since \(\mu\mapsto\int t\,d\mu\) is \(1\)-Lipschitz for \(W_1\), (5)
implies the exact continuity statement
\[
\frac{\max_x|x^\top A_nx|}{n^{3/2}}
\longrightarrow
\sup_{\mu\in\mathcal K_1^\infty}
\left|\int t\,d\mu(t)\right|. \tag{7}
\]
Unlike weak empirical-energy convergence, (6) retains isolated or
subexponential families of maximizing configurations.

## 2. Finite-tuple enrichment

The one-spin set is enough for continuity of the objective but not for
gluing or realization.  For a \(k\)-tuple
\(\mathbf x=(x^1,\ldots,x^k)\) of sign vectors, define the
matrix-valued spectral measure
\[
 \mathbf M_{A,\mathbf x}(B)_{ab}
 =\frac1n\langle x^a,E_A(B)x^b\rangle,
 \qquad 1\le a,b\le k. \tag{8}
\]
It is a positive-semidefinite \(k\times k\) matrix measure.  Its total
matrix records the overlaps:
\[
 \mathbf M_{A,\mathbf x}([-C,C])_{ab}
 =\frac1n\langle x^a,x^b\rangle. \tag{9}
\]
Its first moment records every bilinear action:
\[
 \int t\,d\mathbf M_{A,\mathbf x}(t)_{ab}
 =\frac1n\langle x^a,U_Ax^b\rangle. \tag{10}
\]

Also record the empirical joint spin law
\[
 \nu_{\mathbf x}(\sigma)
 =\frac1n\#\{i:(x_i^1,\ldots,x_i^k)=\sigma\},
 \qquad \sigma\in\{\pm1\}^k. \tag{11}
\]
Let
\[
 \mathcal K_k(A)
 =\{(\nu_{\mathbf x},\mathbf M_{A,\mathbf x}):
 \mathbf x\in(\{\pm1\}^n)^k\}. \tag{12}
\]
For fixed \(k,C\), the ambient space in (12) is compact: the scalar
total variations of all matrix entries are bounded by Cauchy--Schwarz,
and the measures have common compact support.  Thus the Hausdorff
hyperspace is compact.

By a diagonal subsequence, every sequence obeying (1) has a full Boolean
spectral profile limit
\[
 \mathfrak K_\infty=(\mathcal K_k^\infty)_{k\ge1},
 \qquad
 \mathcal K_k(A_{n_j})\to\mathcal K_k^\infty
 \quad\text{for every fixed }k. \tag{13}
\]

This profile retains:

- every finite overlap pattern of selected Boolean configurations;
- every finite matrix of their normalized bilinear energies;
- the spectral distribution seen from each selected configuration;
- isolated ground states, because the profile is a set of attainable
  laws rather than a law under a uniformly random spin.

It is an exact zero-temperature counterpart of action profiles.

## 3. Why this compactification is sufficient for subsequential limits

Let \(\mathscr S_n(C)\) be the set of full profiles of symmetric,
zero-diagonal sign matrices of order \(n\) satisfying (1).  Let
\[
 \overline{\mathscr S}(C)
 =\overline{\bigcup_{n\ge1}\mathscr S_n(C)}
 \tag{14}
\]
in the product Hausdorff topology.

Define
\[
 \Phi(\mathfrak K)
 =\sup_{\mu\in\mathcal K_1}
 \left|\int t\,d\mu(t)\right|. \tag{15}
\]
Then \(\Phi\) is continuous and every purified liminf sequence converges,
after passage to a subsequence, to some
\(\mathfrak K_*\in\overline{\mathscr S}(C)\) with
\[
 \Phi(\mathfrak K_*)
 =\liminf_n\frac{2F_n}{n^{3/2}} \tag{16}
\]
in the \(x^\top Ax\) normalization.

Thus compactness and objective continuity are completely resolved on the
purified class.  No endpoint-action theorem is needed.

## 4. Exact recovery theorem still required

To deduce convergence, one needs the following property at an extremal
limit:

> **All-order Boolean spectral-profile recovery.**  If
> \(\mathfrak K_*\) is an extremal limit of purified sign matrices, then
> for every sufficiently large \(m\) there is a symmetric zero-diagonal
> sign matrix \(B_m\) with
> \[
> d(\mathfrak K(B_m),\mathfrak K_*)\to0
> \quad\text{and}\quad
> \|B_m\|_{\rm op}=O(\sqrt m). \tag{17}
> \]

Continuity of \(\Phi\) would then give
\[
 \limsup_m\frac{2F_m}{m^{3/2}}
 \le\Phi(\mathfrak K_*)
 =\liminf_m\frac{2F_m}{m^{3/2}},
 \]
proving existence of the limit.

The tuple enrichment makes (17) strong enough to support every finite
Boolean test; it does not make (17) automatic.  Compactness yields
microstates only along the original subsequence.

## 5. Exact obstruction to graphon-style blow-up recovery

Suppose a base entry is to be replaced by an \(r\times r\) sign block
\(L\), while preserving the normalized action on the constant fiber
mode.  The required row sum is
\[
 L\mathbf1=(1+o(1))a\sqrt r\,\mathbf1. \tag{18}
\]
Writing
\[
 L=\frac{a}{\sqrt r}J_r+R,\qquad R\mathbf1=o(\sqrt r)\mathbf1,
\]
Frobenius orthogonality gives, in the exact-row-sum case,
\[
\boxed{\|R\|_F^2=r^2-r=(1-o(1))r^2.} \tag{19}
\]
Hence virtually all of the microscopic sign energy is forced into
fiber-orthogonal modes.  A graphon rational blow-up cannot converge to
the original full Boolean spectral profile: it necessarily adds a
leading spectral/action component.

The possibilities already audited are:

- independent residual blocks: a fresh Wigner/Gaussian component;
- Hadamard residual blocks: a persistent tensor-fiber component;
- constant blocks: amplification of the base action by \(\sqrt r\).

The deterministic balanced-channel theorem in
`independent_regular_hadamard_lift.md` shows that independently oriented
regular-Hadamard blocks carry a Boolean quadratic contribution of order
\((nr)^{3/2}\); it is not an invisible profile error.

Therefore standard graphon sampling, rational blow-up, and independent
microscopic completion do not prove (17).

## 6. Compulsory-noise interpretation

Equation (19) suggests that the correct limit object, if one admits a
sampling theorem, cannot consist only of a macroscopic operator.  Exact
sign realization at a new scale injects a full amount of microscopic
variance.  A recoverable extremal object would have to be **absorbing**
under that compulsory noise:

\[
\text{macro profile}
\quad+\quad
\text{fresh fiber profile}
\quad\longmapsto\quad
\text{the same extremal profile}. \tag{20}
\]

Proving such an absorption/fixed-point theorem would be sufficient and is
strictly sharper than generic action convergence.  At present no
contraction or uniqueness theorem for the transformation (20) is known.

## 7. Verdict

The Boolean spectral profile is a compact, exact limit object for the
normalized Boolean quadratic supremum on the purified
\(O(\sqrt n)\)-operator class.  It removes the two weaknesses of ordinary
limits:

1. it retains exponentially rare maximizing configurations;
2. the objective is exactly continuous.

The existence problem is therefore reduced without loss to all-order
recovery (17), or equivalently to an absorption theorem for the compulsory
microscopic sign variance (19).  The variance identity sharply refutes all
standard profile-preserving blow-ups, but does not refute the possibility
that extremal profiles themselves are absorbing.
