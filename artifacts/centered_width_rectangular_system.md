# Centered width and its rectangular two-parameter system

For a symmetric zero-diagonal signing \(A\), let
\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\quad
P(A)=\max_xH_A(x),
\quad
Q(A)=-\min_xH_A(x),
\]
\[
W(A)=\frac{P(A)+Q(A)}2,
\qquad
G_n=\min_AW(A).
\]
For an \(m\times n\) rectangular signing \(B\), put
\[
R_{m,n}
=\min_B\|B\|_{\infty\to1}.
\]

## 1. Exact cut identity

For two Boolean vectors \(x,y\), let
\[
S=\{i:x_i\ne y_i\}.
\]
Edges internal to \(S\) or \(S^c\) contribute equally to \(H_A(x)\) and
\(H_A(y)\), while crossing edges change sign.  Hence
\[
H_A(x)-H_A(y)
=2\sum_{i\in S,j\notin S}a_{ij}x_ix_j.
\]
Maximizing first over \(x,y\), and then over \(S\), gives
\[
\boxed{
W(A)
=\max_{S\subset[n]}
\|A_{S,S^c}\|_{\infty\to1}.}
\tag{1}
\]
Thus \(G_n\) is the minimum possible maximum signed rectangular cut norm.

## 2. Superadditivity

For every vertex partition \(S_1,\ldots,S_k\),
\[
\sum_rP(A[S_r])\le P(A),
\qquad
\sum_rQ(A[S_r])\le Q(A).
\]
Indeed, choose an endpoint witness in each block and independently flip the
whole block; cross contributions average to zero.  Therefore
\[
\sum_rW(A[S_r])\le W(A).
\]
Apply this to an optimizer of order \(n_1+\cdots+n_k\):
\[
\boxed{
G_{n_1+\cdots+n_k}\ge\sum_rG_{n_r}.}
\tag{2}
\]
In particular, \(G_n\) is superadditive.

## 3. Sharp rectangular gluing inequality

Take signings \(A_n,D_m\) and an \(n\times m\) cross signing \(B\), and form
\[
A=\begin{pmatrix}A_n&B\\B^\top&D_m\end{pmatrix}.
\]
For a vertex cut \(I\cup J\), where \(I\subset[n]\) and \(J\subset[m]\),
the cross-block edges in the cut form the two disjoint rectangles
\[
B_{I,J^c},\qquad B_{I^c,J}.
\]
Their signed norms add.  This sum is at most
\(\|B\|_{\infty\to1}\): choose signs optimizing both rectangles, combine
them into full row/column sign vectors, and compare the two full bilinear
values obtained by simultaneously reversing the signs on \(I\) and \(J\).
The unwanted diagonal rectangles are unchanged and the two target
rectangles reverse, so one of the two values has absolute value at least
their norm sum.

Using (1) and the triangle inequality for the two internal cut blocks,
\[
W(A)\le W(A_n)+W(D_m)+\|B\|_{\infty\to1}.
\]
Optimizing gives
\[
\boxed{
G_{n+m}\le G_n+G_m+R_{n,m}.}
\tag{3}
\]

Conversely, the cut separating the two vertex blocks shows
\[
\boxed{R_{n,m}\le G_{n+m}.}
\tag{4}
\]

Rectangular concatenation gives
\[
\boxed{
R_{m_1+m_2,n}\le R_{m_1,n}+R_{m_2,n},\qquad
R_{m,n_1+n_2}\le R_{m,n_1}+R_{m,n_2}.}
\tag{5}
\]
Thus \(R\) is symmetric and separately subadditive.

## 4. These exact axioms do not force a \(3/2\)-scale limit

Let, for sufficiently small fixed \(\varepsilon>0\),
\[
g(x)=x^{3/2}
\left[1+\varepsilon\sin\log\log(x+e^e)\right],
\tag{6}
\]
with an immaterial smooth adjustment near zero, and let
\[
r(x,y)=K\sqrt{xy(x+y)}.
\tag{7}
\]

For small enough \(\varepsilon\), \(g\) is increasing and convex, with
\[
c_-x^{3/2}\le g(x)\le c_+x^{3/2},
\qquad
\sup_x\frac{g'(x)}{\sqrt x}<2c_-.
\]
Choose
\[
\sup_x\frac{g'(x)}{\sqrt x}\le K\le2c_-.
\tag{8}
\]

Convexity and \(g(0)=0\) imply superadditivity:
\[
g(x+y)\ge g(x)+g(y).
\]
If \(y\le x\), then
\[
\begin{aligned}
g(x+y)-g(x)-g(y)
&\le g(x+y)-g(x)\\
&\le K y\sqrt{x+y}\\
&\le K\sqrt{xy(x+y)}
=r(x,y).
\end{aligned}
\]
Thus the analogues of (2)--(3) hold.  Moreover
\[
r(x,y)\le\frac K2(x+y)^{3/2}\le g(x+y),
\]
which is the analogue of (4).

Finally, \(r\) is separately subadditive.  For example, after squaring,
\[
\sqrt{(a+b)n(a+b+n)}
\le
\sqrt{an(a+n)}+\sqrt{bn(b+n)}
\]
reduces to
\[
ab\le\sqrt{ab(a+n)(b+n)}.
\]

Rounding \(g(n),r(m,n)\) by bounded amounts preserves all leading-scale
statements and can be arranged to preserve the inequalities exactly.
Nevertheless,
\[
\frac{g(n)}{n^{3/2}}
=1+\varepsilon\sin\log\log(n+e^e)
\]
does not converge.

Therefore the exact scalar system (2)--(5) admits slow oscillation.  No
Fekete theorem based only on \(G_n\) and \(R_{m,n}\) can prove convergence
at the \(n^{3/2}\) scale.

## 5. Transfer back to the original norm

Let
\[
F_n=\min_A\max\{P(A),Q(A)\}.
\]
For every signing,
\[
W(A)\le\max\{P(A),Q(A)\},
\]
so \(G_n\le F_n\).

If, for every \(n\), there is a **centered-width minimizer** \(A_n\) with
\[
W(A_n)=G_n,\qquad |P(A_n)-Q(A_n)|=O(1),
\]
then
\[
\max\{P(A_n),Q(A_n)\}
=W(A_n)+\frac12|P(A_n)-Q(A_n)|
=G_n+O(1).
\]
Consequently
\[
\boxed{0\le F_n-G_n=O(1),}
\tag{9}
\]
and convergence (and the limiting value) for \(G_n/n^{3/2}\) transfers to
\(F_n/n^{3/2}\).

Centering only a minimizer of \(F_n\) does **not** by itself prove the
reverse comparison \(F_n\le G_n+O(1)\); the required hypothesis must apply
to a \(G_n\)-minimizer, or come with a centering operation that preserves
width.

