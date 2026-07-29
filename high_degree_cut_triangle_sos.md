# Growing-degree cut-triangle/SOS audit: induced hafnians and a sharp cancellation barrier

## 1. Setup

Let
\[
q_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad a_{ij}\in\{\pm1\},
\qquad x\in\{\pm1\}^n.
\]
Writing \(v_{ij}=x_ix_j\), the identities
\[
v_{ij}v_{jk}v_{ki}=1
\]
are already built into the Walsh algebra on the vertex cube.  The aim
of this audit was to test a genuinely growing-degree Fourier/SOS
certificate, rather than another fixed-moment relaxation.

The outcome is an exact hierarchy whose top layer is governed by
induced signed hafnians.  It gives a universal degree-\(\Theta(n)\)
certificate, but that certificate is only \(O(n)\).  Moreover, an
explicit signing has complete hafnian cancellation at every even
order, proving that random-scale anti-cancellation is false by a
superexponential factor at the endpoint.

## 2. Exact top-Fourier identity

For an even set \(S\), let
\[
\operatorname{haf}(A_S)
=\sum_{\mathcal M\in\operatorname{PM}(S)}
\prod_{\{i,j\}\in\mathcal M}a_{ij}
\]
be the signed hafnian of the principal submatrix on \(S\).

### Proposition 2.1

For every \(1\le k\le n/2\) and every \(S\subset[n]\) with
\(|S|=2k\),
\[
\boxed{
\widehat{q_A^k}(S)=k!\operatorname{haf}(A_S).}
\tag{2.1}
\]
Consequently,
\[
\boxed{
\mathbb E_x q_A(x)^{2k}
\ge
(k!)^2
\sum_{\substack{S\subset[n]\\|S|=2k}}
\operatorname{haf}(A_S)^2.}
\tag{2.2}
\]

### Proof

Expand \(q_A^k\) as a sum over ordered \(k\)-tuples of edges.
The Walsh support of a product is the set of vertices of odd degree
in that edge multigraph.  To obtain support of size \(2k\), all \(k\)
edges must be pairwise disjoint and have union \(S\).  Thus they form
a perfect matching of \(S\), and each matching occurs in all \(k!\)
orders.  This proves (2.1).  Parseval applied to \(q_A^k\) gives
\[
\mathbb E q_A^{2k}
=\sum_T\widehat{q_A^k}(T)^2;
\]
retaining only the level-\(2k\) terms gives (2.2).

Every \(\operatorname{haf}(A_S)\) is a sum of
\((2k-1)!!\) signs.  Since this number is odd,
\[
\operatorname{haf}(A_S)\in2\mathbb Z+1,
\qquad
|\operatorname{haf}(A_S)|\ge1.
\tag{2.3}
\]
Hence
\[
\boxed{
\|q_A\|_\infty
\ge
\left[
\binom n{2k}(k!)^2
\right]^{1/(2k)}.}
\tag{2.4}
\]

This is a valid certificate at degree \(2k=\Theta(n)\), but it has
the wrong scale.  If \(k=\alpha n\), \(0<\alpha\le1/2\), then
\[
\left[
\binom n{2k}(k!)^2
\right]^{1/(2k)}
=
n\,\frac{\alpha}{e}
\exp\left\{\frac{H(2\alpha)}{2\alpha}\right\}
(1+o(1)),
\tag{2.5}
\]
where the product in (2.5) is intended as
\[
n\left[
\frac{\alpha}{e}
\exp\left\{\frac{H(2\alpha)}{2\alpha}\right\}
\right](1+o(1)).
\]
Equivalently, and without typographical ambiguity,
\[
\frac1n
\left[
\binom n{2k}(k!)^2
\right]^{1/(2k)}
\longrightarrow
\frac{\alpha}{e}
\exp\left\{\frac{H(2\alpha)}{2\alpha}\right\}.
\tag{2.6}
\]
Thus it is \(\Theta(n)\), a factor \(\sqrt n\) below the required
\(n^{3/2}\) scale.

## 3. Random-scale anti-cancellation is false

A tempting strengthening of (2.3) is
\[
\sum_{|S|=2k}\operatorname{haf}(A_S)^2
\stackrel{?}{\gtrsim}
\binom n{2k}(2k-1)!!\,e^{-O(k)}.
\tag{3.1}
\]
The factor \((2k-1)!!\) is the exact mean square of a random signed
hafnian: only equal matching pairs survive after averaging the edge
signs.  The following construction disproves (3.1) at
\(n=2k\), even with an arbitrary \(e^{-O(k)}\) loss.

### Theorem 3.1

For every \(k\ge1\), there is a symmetric zero-diagonal
\(\{\pm1\}\) matrix \(B_k\) of order \(2k\) such that
\[
\boxed{\operatorname{haf}(B_k)=1.}
\tag{3.2}
\]

Label the vertices
\[
L_1,R_1,L_2,R_2,\ldots,L_k,R_k.
\]
Put every edge sign equal to \(+1\), except
\[
b_{R_r,L_s}=-1\qquad(r<s).
\tag{3.3}
\]

### Proof

On polynomials in commuting variables \(u,v\), define
\[
\mathsf L=u+\partial_u-\partial_v,\qquad
\mathsf R=v+\partial_u+\partial_v.
\tag{3.4}
\]
When the word
\[
\mathsf L,\mathsf R,\mathsf L,\mathsf R,\ldots,
\mathsf L,\mathsf R
\]
is Wick-expanded and evaluated at \(u=v=0\), the pair contraction
between two operators is:

\[
\begin{array}{c|cc}
\text{earlier}\backslash\text{later}&\mathsf L&\mathsf R\\ \hline
\mathsf L&1&1\\
\mathsf R&-1&1
\end{array}.
\tag{3.5}
\]

These are exactly the edge signs in (3.3).  Therefore
\[
\operatorname{haf}(B_k)
=
\left[(\mathsf R\mathsf L)^k1\right]_{u=v=0}.
\tag{3.6}
\]

For completeness, this vacuum expectation can be evaluated exactly.
Let
\[
\mathcal F(t,u,v)
=e^{t\mathsf R\mathsf L}1.
\]
Direct multiplication gives
\[
\mathsf R\mathsf L
=1+uv+(u+v)\partial_u+(u-v)\partial_v
+\partial_u^2-\partial_v^2.
\tag{3.7}
\]
Use the Gaussian ansatz
\[
\mathcal F
=C(t)\exp\left\{
\frac{a(t)}2(u^2+v^2)+b(t)uv
\right\}.
\tag{3.8}
\]
The Riccati equations preserve
\[
a+a^2-b^2=0,
\tag{3.9}
\]
because
\[
a'=2b,\qquad b'=1+2a.
\tag{3.10}
\]
At \(t=0\), \(a=b=0\), so (3.9) holds for all \(t\).  The constant
equation is therefore
\[
\frac{C'}C=1+a-a=1,\qquad C(0)=1.
\]
Hence
\[
\mathcal F(t,0,0)=C(t)=e^t.
\tag{3.11}
\]
Comparing coefficients of \(t^k/k!\) in (3.6)--(3.11) proves
\(\operatorname{haf}(B_k)=1\).

At \(n=2k\), the sum in (3.1) has only the full set \(S=[n]\), so its
left side is \(1\), whereas its proposed random benchmark is
\[
(2k-1)!!
=\exp\{k\log(2k)-k+O(\log k)\}.
\]
Thus the gap is superexponential in \(k\).  Parity/nonvanishing in
(2.3) is exactly sharp at the linear-degree endpoint.

## 4. A fixed-degree contrast: the \(k=2\) average is random-scale

The endpoint counterexample does not mean that cancellation is free
at every fixed degree.  For \(|S|=4\),
\[
\operatorname{haf}(A_S)
=a_{12}a_{34}+a_{13}a_{24}+a_{14}a_{23}
\in\{\pm1,\pm3\}.
\]
Let \(B_4(A)\) be the number of four-sets on which its magnitude is
\(3\).  Then
\[
\sum_{|S|=4}\operatorname{haf}(A_S)^2
=\binom n4+8B_4(A).
\tag{4.1}
\]

Fix a root and switch all incident edges to \(+1\).  A four-set
containing the root is counted by \(B_4(A)\) exactly when the
remaining triangle is monochromatic in the induced two-coloring.
Goodman's monochromatic-triangle bound therefore gives
\[
B_4(A)\ge\left(\frac14-o(1)\right)\binom n4,
\]
and hence
\[
\boxed{
\sum_{|S|=4}\operatorname{haf}(A_S)^2
\ge(3-o(1))\binom n4.}
\tag{4.2}
\]
The constant \(3\) is the random-sign benchmark
\((2k-1)!!\) at \(k=2\).

So anti-cancellation is real at fixed degree and can be proved from a
local Ramsey theorem, but it collapses completely at the
degree-\(n\) endpoint.  There is no uniform fixed-to-growing-degree
transfer.

## 5. Verdict for high-degree triangle/SOS methods

The top Fourier layer of \(q_A^k\) is a clean, exact way to use the
triangle/cut algebra at degree \(\Theta(n)\).  It yields the following
rigorous barrier:

> Any high-degree certificate that uses only the existence, parity,
> or random-scale lower bounds of the top Fourier coefficients
> \(k!\operatorname{haf}(A_S)\) cannot prove the
> \(n^{3/2}\) lower bound.  Parity gives only \(O(n)\), and the
> random-scale strengthening is false by a superexponential factor.

This does **not** rule out the full degree-\(\Theta(n)\) SOS cone.  The
explicit \(B_k\) can have large information in lower Fourier levels,
and a full certificate is allowed to combine all of them.  A viable
successor would need a multilevel conservation law: when the
level-\(2k\) hafnians cancel as in Theorem 3.1, mass must be forced
into lower Fourier levels in a quantitatively usable way.  Parseval
alone is circular because the sum over all levels is exactly
\(\mathbb E q_A^{2k}\).

The next precise target is therefore a signing-uniform inequality of
the form
\[
\sum_{j\le k}w_{k,j}
\sum_{|S|=2j}\widehat{q_A^k}(S)^2
\ge
\exp\{2k(\tfrac32\log n+O(1))\},
\]
with explicit nonnegative weights \(w_{k,j}\) that are compatible
with restriction or scale composition.  Theorem 3.1 shows that no
single top level can carry such an estimate.
