# Spectral localization versus proportional restriction

## 1. Main correction

Let
\[
Q(A)=\max_{x\in\{\pm1\}^n}|x^\top Ax|
\]
for a symmetric zero-diagonal matrix with \(|a_{ij}|\le1\).
The previously used bootstrap
\[
\|A\|_{\rm op}^3\le nQ(A)
\]
is not the strongest elementary relation.

### Theorem 1.1

\[
\boxed{\|A\|_{\rm op}^2\le2Q(A).}
\tag{1.1}
\]

### Proof

Because the diagonal is zero, \(z^\top Az\) is affine in each
coordinate separately.  Therefore
\[
\max_{z\in[-1,1]^n}|z^\top Az|=Q(A).
\tag{1.2}
\]
For sign vectors \(x,y\), put
\[
u=\frac{x+y}{2},\qquad v=\frac{x-y}{2}.
\]
Symmetry gives the exact polarization identity
\[
x^\top Ay=u^\top Au-v^\top Av.
\]
Since \(u,v\in[-1,1]^n\),
\[
\|A\|_{\infty\to1}
=\max_{x,y\in\{\pm1\}^n}|x^\top Ay|
\le2Q(A).
\tag{1.3}
\]
Also
\[
\|A\|_{1\to\infty}=\max_{i,j}|a_{ij}|\le1.
\]
Riesz--Thorin interpolation between
\(\ell_1\to\ell_\infty\) and
\(\ell_\infty\to\ell_1\) at parameter \(1/2\) gives
\[
\|A\|_{2\to2}
\le
\sqrt{\|A\|_{1\to\infty}\|A\|_{\infty\to1}},
\]
which proves (1.1).

For a near-minimizing signing, \(Q(A)=O(n^{3/2})\), so
\[
\boxed{\|A\|_{\rm op}=O(n^{3/4}),}
\tag{1.4}
\]
improving the previous \(O(n^{5/6})\) bootstrap.

## 2. Correct selector exponent

For an independent Bernoulli-\(p\) selector
\(\delta=p{\bf1}+\xi\), the restricted energy around a switched spin
has the decomposition
\[
\frac12\delta^\top B\delta
=p^2H_A(x)+p\,\xi^\top B{\bf1}
+\frac12\xi^\top B\xi.
\]
The linear field has variance proxy
\[
\|B{\bf1}\|_2^2\le n\|A\|_{\rm op}^2.
\]
At a fixed proportional restriction, the relevant gap is
\(u=\Theta(n^{3/2})\).  With the old
\(\|A\|_{\rm op}=O(n^{5/6})\), the two operator-norm exponents were
\[
\frac{u^2}{n\|A\|_{\rm op}^2}=n^{1/3},
\qquad
\frac{u}{\|A\|_{\rm op}}=n^{2/3}.
\]
The minimum is \(n^{1/3}\); any earlier \(n^{2/3}\) report omitted the
dominant linear selector term.

Using (1.4), the same branches become
\[
\boxed{
\frac{u^2}{n\|A\|_{\rm op}^2}=\Theta(n^{1/2}),
\qquad
\frac{u}{\|A\|_{\rm op}}=\Theta(n^{3/4}).}
\tag{2.1}
\]
Thus the corrected norm-only selector exponent is \(n^{1/2}\).  This
is a material improvement, but still subexponential in \(n\), so it
does not by itself compensate a positive energy-layer entropy rate.

## 3. Exact clipped-tail localization of one eigenvector

There is a useful quantitative localization lemma which improves the
single-coordinate spike statement.

### Lemma 3.1

Let
\[
Au=\lambda u,\qquad \|u\|_2=1,\qquad\lambda>0.
\]
Order the coordinates so that
\[
|u_1|\ge|u_2|\ge\cdots\ge|u_n|.
\]
For \(T=\{1,\ldots,t\}\), put
\[
a_t=\|u_T\|_2^2.
\]
Then
\[
\boxed{
a_t\ge
\frac{\bigl[\lambda-Q(A)/(t+1)\bigr]_+}
{2\lambda+t-1}.}
\tag{3.1}
\]

### Proof

Put \(w=u_{T^c}\).  The eigenvector equation gives
\[
\begin{aligned}
w^\top Aw
&=(u-u_T)^\top A(u-u_T)\\
&=\lambda-2\lambda a_t+u_T^\top Au_T.
\end{aligned}
\]
The principal \(t\times t\) signing has operator norm at most \(t-1\),
so
\[
w^\top Aw\ge\lambda-a_t(2\lambda+t-1).
\tag{3.2}
\]
On the other hand,
\[
\|w\|_\infty^2\le\frac1{t+1}.
\]
Scaling \(w\) into the cube and using (1.2) gives
\[
|w^\top Aw|\le\frac{Q(A)}{t+1}.
\tag{3.3}
\]
Combining (3.2) and (3.3) proves (3.1).

If
\[
t+1=s\frac{Q(A)}{\lambda},\qquad s>1,
\]
then
\[
\boxed{
a_t\ge
\frac{1-1/s}{2+sQ(A)/\lambda^2}.}
\tag{3.4}
\]
In particular, an eigenvalue
\(\lambda\asymp\sqrt{Q(A)}\) has a constant fraction of its
eigenvector mass on
\[
O\left(\frac{Q(A)}{\lambda}\right)
=O(\sqrt{Q(A)})
\]
coordinates.  For \(Q(A)=O(n^{3/2})\), this is
\(O(n^{3/4})=o(n)\).

For smaller eigenvalues, however, (3.4) captures only a fraction
\(\Theta(\lambda^2/Q)\).  Iterating the estimate separately over all
eigenvectors down to \(C\sqrt n\) does not currently fit in an
\(o(n)\) deletion budget.

## 4. Sharp obstruction to a norm-only selector improvement

The \(n^{1/2}\) exponent in (2.1) is attained by the natural coherent
rectangular spike.

Take a set \(T\) of size \(s=\lfloor\sqrt n\rfloor\), put all cross
edges between \(T\) and \(T^c\) equal to \(+1\), and fill the two
diagonal blocks by any signings with \(O(n^{3/2})\) Boolean norm.
The complete bipartite part has
\[
\|A\|_{\rm op}=\Theta(\sqrt{s(n-s)})
=\Theta(n^{3/4})
\]
and Boolean quadratic norm \(\Theta(s(n-s))=\Theta(n^{3/2})\).
In the coherent switching, the \(s\) rows on the small side have
field size \(\Theta(n)\), while the other rows have field
\(\Theta(\sqrt n)\).  Hence
\[
\|A{\bf1}\|_2^2=\Theta(n^{5/2}),
\]
and a deviation \(u=\Theta(n^{3/2})\) has linear-tail exponent
\[
\frac{u^2}{\|A{\bf1}\|_2^2}=\Theta(n^{1/2}).
\]

This example is also instructive rather than purely negative:
deleting the coherent side \(T\), only \(O(\sqrt n)=o(n)\) vertices,
removes the entire spike.  It suggests the exact missing theorem:

> If a near-minimizing signing has a switched field with
> \(\ell_2^2\) large enough to reduce the selector exponent below
> \(\Theta(n)\), then all such bad fields share a common
> \(o(n)\)-vertex coherent support whose deletion leaves
> \(O(\sqrt n)\) operator norm.

Lemma 3.1 proves this only for a single eigenvector near the
\(\sqrt Q\) endpoint.  It does not yet globalize the support across
all high eigenvectors or all high-energy spin switchings.

## 5. Current verdict

The spectral route has made two rigorous advances:

1. the universal near-minimizer bootstrap is \(n^{3/4}\), not
   \(n^{5/6}\);
2. eigenvalues near that cap force constant mass onto an
   \(o(n)\)-coordinate set.

But no proof was found that deleting \(o(n)\) vertices reduces every
near-minimizer to \(O(\sqrt n)\) operator norm.  The quantitative wall
is the band
\[
\sqrt n\ll\lambda\ll\sqrt{Q(A)},
\]
where (3.4) captures only \(\lambda^2/Q(A)\) of each eigenvector.
Naively unioning the individual localization sets can cost order
\(n\) or more.

A successful continuation must use common structure across this whole
spectral band, not one-eigenvector truncation.  The coherent
complete-bipartite spike shows what the extremal obstruction should
look like and why a common-support theorem, if true, would be exactly
the right input to the multiplicity-refined restriction lemma.
