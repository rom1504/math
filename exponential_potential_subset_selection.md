# Exponential-potential subset selection: exact contraction and the
# changing-temperature obstruction

## Setup

For \(S\subset[n]\), define

\[
H_S(x)=\sum_{\{i,j\}\subset S}a_{ij}x_ix_j
\]

and the two-sided soft maximum

\[
Z_S(\lambda)
=\sum_{x\in\{\pm1\}^S}2\cosh(\lambda H_S(x)).
\]

Then

\[
\max_x|H_S(x)|\le\lambda^{-1}\log Z_S(\lambda).
\tag{1}
\]

## 1. Exact conditional-Jensen inequality

Fix \(x_S\), and fill the coordinates outside \(S\) with independent
uniform signs \(X_{S^c}\).  Every edge touching \(S^c\) has conditional
mean zero, so

\[
\mathbb E[H_{[n]}(X)\mid X_S=x_S]=H_S(x_S).
\]

Convexity of \(\cosh\) therefore gives

\[
2\cosh(\lambda H_S(x_S))
\le
2^{-|S^c|}
\sum_{\substack{y\in\{\pm1\}^n\\y_S=x_S}}
2\cosh(\lambda H_{[n]}(y)).
\]

Summing over \(x_S\) proves, for every individual subset,

\[
\boxed{
Z_S(\lambda)\le2^{|S|-n}Z_{[n]}(\lambda).
}
\tag{2}
\]

Consequently averaging over \(m\)-subsets, taking the minimum, or using
the elementary-symmetric generating function gives no improvement over
(2) without additional information about the full energy distribution.

The exact subset generating function is

\[
\sum_{S\subset[n]}u^{|S|}Z_S(\lambda)
=
\sum_{z\in\{0,\pm1\}^n}
u^{|\operatorname{supp}z|}
2\cosh\!\left(
\lambda\sum_{i<j}a_{ij}z_iz_j
\right),
\tag{3}
\]

but its coefficient positivity alone reproduces only conditional
contraction, not the desired \(3/2\)-homogeneous scaling.

## 2. Normalized temperature necessarily changes

Let \(m=\alpha n\).  The natural zero-temperature scaling for the subset
is

\[
\lambda=\frac{\beta}{\sqrt m}
=\frac{\beta/\sqrt\alpha}{\sqrt n}.
\]

Thus (2), evaluated at the subset temperature \(\beta\), invokes the full
system at temperature \(\beta/\sqrt\alpha\).  Restriction heats the
normalized full pressure by the fixed factor \(\alpha^{-1/2}\).

Using only \(|H_{[n]}|\le M(A)\) in (2) yields

\[
M(A_S)
\le
M(A)+\frac{(m+1)\log2}{\lambda}.
\tag{4}
\]

Optimizing (4) sends \(\lambda\to\infty\) and returns only
\[
M(A_S)\le M(A).
\]
It cannot produce
\(\alpha^{3/2}M(A)+o(n^{3/2})\).

## 3. Why universal conditional-variance corrections are insufficient

One might seek a reverse Jensen gain from the variables outside \(S\).
Their quadratic energy has variance \(\Theta(|S^c|^2)\), but variance
alone forces only

\[
\mathbb E\cosh(\lambda Y)
\ge
\cosh\!\left(\lambda\sqrt{\mathbb EY^2}\right).
\tag{5}
\]

At \(\lambda=\Theta(n^{-1/2})\), the logarithm of (5) is only
\(O(\sqrt n)\), whereas removing the spin-entropy term in (1) requires
an \(O(n)\) gain.

The same limitation appears algebraically.  Edge characters
\(x_ix_j\) span the even-weight character space, of rank only \(n-1\);
one can extract at most \(O(n)\) independent quadratic characters.
A matching or a star therefore supplies only
\((\cosh\lambda)^{O(n)}=\exp(O(1))\) at
\(\lambda=\Theta(n^{-1/2})\), again far short of \(\exp(\Theta(n))\).

An exponential gain is possible only from a genuine large-deviation
statement about the full energy distribution.  That statement is
precisely the missing entropic/ROM-versus-SK decision, not a consequence
of the subset potential formalism.

## Verdict

The exponential-potential method gives a clean exact contraction, but at
the wrong normalized temperature.  Elementary-symmetric averaging and
interlacing based solely on positive partition-function coefficients do
not remove this mismatch.  Progress by this route would require a new
uniform free-energy large-deviation theorem strong enough to compare
temperatures \(\beta\) and \(\beta/\sqrt\alpha\); absent that input, the
method collapses to the trivial restriction inequality.
