# A one-sided energy-product theorem for near-minimal signings

## Normalization

Let \(A\) be a symmetric \(n\times n\) matrix with zero diagonal and
off-diagonal entries in \(\{\pm1\}\).  Write
\[
P(A)=\max_{x\in\{\pm1\}^n}x^\top Ax,\qquad
N(A)=-\min_{x\in\{\pm1\}^n}x^\top Ax,
\]
\[
Q(A)=\max\{P(A),N(A)\},\qquad R(A)=P(A)+N(A).
\]
Thus \(Q(A)=2M(A)\) in the half-energy normalization of the original
problem.

## External input

Bollobás and Scott proved the following theorem.  If \(G\) is an
\(n\)-vertex graph of density \(p\), and
\[
\operatorname{disc}^{+}(G)
=\max_U\left(e_G(U)-p\binom{|U|}{2}\right),\qquad
\operatorname{disc}^{-}(G)
=\max_U\left(p\binom{|U|}{2}-e_G(U)\right),
\]
then, provided \(p(1-p)\ge 1/n\),
\[
\operatorname{disc}^{+}(G)\operatorname{disc}^{-}(G)
\ge \frac{p(1-p)n^3}{6400}. \tag{1}
\]

Reference: B. Bollobás and A. D. Scott, *Discrepancy in graphs and
hypergraphs*, More Sets, Graphs and Numbers (2006), 33--56.

## Translation to Boolean quadratic energies

Switch \(A\) by a positive ground state, so that
\(\mathbf1^\top A\mathbf1=P(A)\).  Let \(G\) contain precisely the
edges \(ij\) for which \(a_{ij}=1\).  Its density is
\[
p=\frac12+\frac{P(A)}{2n(n-1)}
=\frac{1+r_P}{2},\qquad
r_P=\frac{P(A)}{n(n-1)}. \tag{2}
\]

For \(U\subset[n]\), put
\[
S_A(U)=\sum_{\{i,j\}\subset U}a_{ij}.
\]
Then
\[
D_G(U):=e_G(U)-p\binom{|U|}{2}
=\frac12\left(
S_A(U)-r_P\binom{|U|}{2}
\right). \tag{3}
\]

Fix all spins in \(U\) to \(1\), and choose the spins outside \(U\)
independently and uniformly.  The expected doubled energy is
\[
\mathbb E\,X^\top AX=2S_A(U).
\]
Every realization lies in \([-N(A),P(A)]\), hence
\[
-\frac{N(A)}2\le S_A(U)\le\frac{P(A)}2. \tag{4}
\]
Since \(r_P\ge0\), (3)--(4) imply the exact bounds
\[
\operatorname{disc}^{+}(G)\le\frac{P(A)}4,\qquad
\operatorname{disc}^{-}(G)\le\frac{R(A)}4. \tag{5}
\]

Combining (1), (2), and (5) gives
\[
\boxed{
P(A)R(A)
\ge
\frac{(1-r_P^2)n^3}{1600}
} \tag{6}
\]
whenever \((1-r_P^2)/4\ge1/n\).

Apply the same argument to \(-A\), switched by a positive ground state
of \(-A\).  This exchanges \(P\) and \(N\), and yields
\[
\boxed{
N(A)R(A)
\ge
\frac{(1-r_N^2)n^3}{1600},
\qquad r_N=\frac{N(A)}{n(n-1)}.
} \tag{7}
\]

## Near-minimizer consequence

For every sequence \(A_n\) satisfying \(Q(A_n)=O(n^{3/2})\), both
\(r_P\) and \(r_N\) tend to zero.  Equations (6)--(7) therefore give
\[
P(A_n)R(A_n),\,N(A_n)R(A_n)
\ge (1-o(1))\frac{n^3}{1600}. \tag{8}
\]
As \(R\le2Q\), if \(m=\min(P,N)\), then
\[
mQ\ge (1-o(1))\frac{n^3}{3200}.
\]
Since \(Q=\max(P,N)\), this is the product bound
\[
\boxed{
P(A_n)N(A_n)
\ge (1-o(1))\frac{n^3}{3200}.
} \tag{9}
\]

In particular, for every fixed \(C\), if \(Q(A)\le Cn^{3/2}\), then
both one-sided extrema are at least
\[
(1-o(1))\frac{n^{3/2}}{3200C}. \tag{10}
\]
Thus a signing competitive at the \(n^{3/2}\) scale cannot hide almost
all of its absolute norm in only one orientation.

## What this does and does not settle

The result is stronger than the usual absolute discrepancy lower
bound because it controls the two one-sided extrema simultaneously.
It is directly usable in recursive block arguments whose obstruction
was orientation imbalance.

It does not by itself prove convergence: its constant is small and it
does not control how the one-sided energy is distributed among
principal blocks.  The next tests are:

1. combine (8) with macroscopic ground-closure splits;
2. combine (9) with affine-ground type bounds;
3. test whether regularity or the local-field constraints improve
   \(\operatorname{disc}^{-}\le R/4\) to a bound involving \(N\)
   alone, which would sharpen the constant and the recursion.
