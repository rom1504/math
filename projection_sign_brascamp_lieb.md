# Conference projection stable-tail audit: Brascamp--Lieb and equitable sets

## Status

This note records two rigorous reductions for the proposed stable-tail
statement

\[
\Sigma_n(\varepsilon)
=\frac1n\log\#\left\{
y\in\{\pm1\}^n:
y=\operatorname{sign}(P_+y),\
\|P_-y\|_2^2\le\frac{\varepsilon n}{2}
\right\},
\]

where

\[
U=\frac{C}{\sqrt{n-1}},\qquad
P_\pm=\frac{I\pm U}{2},
\qquad C^2=(n-1)I.
\]

The desired conclusion is

\[
\lim_{\varepsilon\downarrow0}\limsup_n\Sigma_n(\varepsilon)=0.
\tag{T}
\]

The results below do not prove (T).  They identify an exact universal
Laplace bound and its sharp \(2^{n/2}\) wall, and reduce the zero-defect
case to counting equitable switch sets in a conference strongly regular
graph.

## 1. Exact local-field identities

For \(y\in\{\pm1\}^n\), put

\[
r(y)=\frac{y^\top Uy}{n},
\qquad
h_i=y_i(Uy)_i.
\]

Then

\[
\|P_-y\|_2^2=\frac{n-y^\top Uy}{2}
=\frac{1-r(y)}2n,
\tag{1}
\]

and orthogonality of \(U\) gives

\[
\boxed{\|h-\mathbf1\|_2^2=2n(1-r(y))
=4\|P_-y\|_2^2.}
\tag{2}
\]

Indeed, \(\sum_i h_i=y^\top Uy\) and
\(\sum_i h_i^2=\|Uy\|_2^2=n\).  Stability is exactly

\[
y=\operatorname{sign}(P_+y)
\quad\Longleftrightarrow\quad
h_i\ge-1\quad\text{for every }i,
\tag{3}
\]

with ties resolved in favor of the current sign.

## 2. A finite-\(n\) Gaussian Brascamp--Lieb bound

The following theorem uses only that \(P=P_+\) is a rank-\(n/2\)
orthogonal projection with constant diagonal \(1/2\).

### Theorem 2.1

For every \(t\ge0\),

\[
\boxed{
\mathbb E_{y\in\{\pm1\}^n}
\exp\!\left(t\,y^\top Uy\right)
\le
e^{-tn}
\left(\frac{1+e^{4t}}2\right)^{n/2}.}
\tag{4}
\]

Consequently, for \(0<\varepsilon\le1\),

\[
\boxed{
\frac1n\log
\#\{y:y^\top Uy\ge(1-\varepsilon)n\}
\le B_{\rm BL}(\varepsilon),}
\tag{5}
\]

where

\[
B_{\rm BL}(\varepsilon)
=\frac{\log2}{2}
+\frac{\varepsilon}{4}
\log\frac{2-\varepsilon}{\varepsilon}
-\frac12\log\!\left(1-\frac{\varepsilon}{2}\right).
\tag{6}
\]

In particular,

\[
B_{\rm BL}(\varepsilon)
=\frac{\log2}{2}
+\frac{\varepsilon}{4}\log\frac2\varepsilon
+O(\varepsilon),
\tag{7}
\]

so constant leverage and geometric Brascamp--Lieb alone stop at the
universal \(2^{n/2}\) bound and cannot prove (T).

### Proof

Let \(g\) be a standard Gaussian in \(\operatorname{Range}(P)\).
Since \(U=2P-I\), the Hubbard--Stratonovich identity gives

\[
\begin{aligned}
\mathbb E_y e^{t y^\top Uy}
&=e^{-tn}\,
\mathbb E_g\mathbb E_y
\exp(2\sqrt t\,g^\top y)\\
&=e^{-tn}\,
\mathbb E_g\prod_{i=1}^n\cosh(2\sqrt t\,g_i).
\end{aligned}
\tag{8}
\]

On \(\operatorname{Range}(P)\), define

\[
v_i=\sqrt2\,Pe_i.
\]

Then \(\|v_i\|_2=1\), \(g_i=\langle v_i,g\rangle/\sqrt2\), and

\[
\sum_{i=1}^n\frac12\,v_iv_i^\top=I_{\operatorname{Range}(P)}.
\tag{9}
\]

The geometric Gaussian Brascamp--Lieb inequality, with weights \(1/2\)
and one-dimensional functions

\[
f_i(s)=\cosh^2(\sqrt{2t}\,s),
\]

therefore yields

\[
\mathbb E_g\prod_i
\cosh(\sqrt{2t}\langle v_i,g\rangle)
\le
\prod_i
\left(\mathbb E_{G\sim N(0,1)}
\cosh^2(\sqrt{2t}\,G)\right)^{1/2}.
\tag{10}
\]

The scalar expectation is

\[
\mathbb E\cosh^2(\sqrt{2t}\,G)
=\frac{1+e^{4t}}2.
\tag{11}
\]

Substitution into (8) proves (4).  Chernoff's inequality gives

\[
\frac1n\log\#\{y:y^\top Uy\ge(1-\varepsilon)n\}
\le
\frac{\log2}{2}
+\varepsilon t
+\frac12\log(1+e^{-4t}).
\tag{12}
\]

The minimizer satisfies

\[
e^{4t}=\frac2\varepsilon-1.
\]

Substitution gives (5)--(7).

## 3. Exact cap states are equitable switch sets

Assume an exact positive Boolean eigenvector exists.  Necessarily
\[
k:=\sqrt{n-1}\in\mathbb Z,
\tag{13}
\]
because \(Cy=ky\) has an integral left side and a nonzero integral
right vector.

Switch by one such eigenvector so that

\[
C\mathbf1=k\mathbf1.
\]

Let

\[
A=\frac{J-I-C}{2}.
\]

Then \(A\) is the adjacency matrix of a regular graph of degree

\[
d=\frac{n-1-k}{2}.
\tag{14}
\]

For a second positive Boolean eigenvector \(w\), put

\[
S=\{i:w_i=-1\},\qquad s=|S|.
\]

Since \(w=\mathbf1-2\mathbf1_S\),

\[
Cw=kw
\quad\Longleftrightarrow\quad
C\mathbf1_S=k\mathbf1_S.
\tag{15}
\]

Equivalently,

\[
\boxed{
\begin{aligned}
\deg_S(v)&=\frac{s}{2},
&&v\notin S,\\
\deg_{G[S]}(v)&=\frac{s-1-k}{2},
&&v\in S.
\end{aligned}}
\tag{16}
\]

Thus every external vertex sees exactly half of \(S\), while the
induced graph on \(S\) is regular.  In particular,

\[
s\ge k+1=\sqrt{n-1}+1.
\tag{17}
\]

The centered indicator

\[
z_S=\mathbf1_S-\frac{s}{n}\mathbf1
\]

lies in the restricted adjacency eigenspace

\[
Az_S=-\frac{k+1}{2}z_S.
\tag{18}
\]

Conversely, every subset satisfying (16) gives a second positive
Boolean eigenvector.  Hence the zero-defect counting problem is exactly
the problem of counting these equitable two-cell partitions in a
conference strongly regular graph.

Equation (18) by itself again gives only the dimension bound
\(2^{n/2}\).  Expander mixing is also an equality here rather than an
entropy loss: the discrepancy of \(S\) is precisely supported on the
restricted eigenspace.  Any successful count must use the simultaneous
coordinate constraints in (16), not just the graph spectrum.

## 4. Hadamard/bent obstruction audit

For a symmetric normalized Hadamard involution, Boolean top
eigenvectors are self-dual bent sequences.  These examples show that a
coherence-only theorem cannot be stronger than subexponential:

- a Bush-type Hadamard matrix of order \(N=4u^2\) has at least
  \(2^{2u}=2^{\sqrt N}\) Boolean top eigenvectors;
- for the Sylvester Hadamard matrix at truth-table length \(N=64\),
  the published exact count is \(42{,}896\) self-dual bent functions.

Neither is a counterexample to (T).  The normalized Hadamard
involutions have diagonal entries \(\pm N^{-1/2}\), while an exact
conference involution has

\[
U_{ii}=0,\qquad |U_{ij}|=(n-1)^{-1/2}\quad(i\ne j).
\tag{19}
\]

Deleting the Hadamard diagonal destroys the involution identity, so
the large self-dual-bent families do not transfer directly.  No
positive-rate family of exact Boolean eigenvectors, or stable
vanishing-defect vectors, was found for exact conference matrices.

## 5. Verdict and exact missing lemma

The generic projection route is exhausted at (6).  The remaining
conference-specific target can be stated finitely:

> For every conference strongly regular graph \(G\) of order \(n\),
> the number of subsets satisfying (16), and the number of approximate
> versions obtained from stable defect-\(\varepsilon\) switchings, is
> at most
> \[
> \exp\{n\,\sigma(\varepsilon)+o(n)\},
> \qquad \sigma(\varepsilon)\to0.
> \]

At zero defect this is an equitable-partition enumeration theorem for
conference regular two-graphs.  At positive defect, (2) says that all
but \(O(\varepsilon n/\tau^2)\) switched row sums lie within
\(\tau\sqrt n\) of \(+\sqrt n\).  A proof must convert those
simultaneous near-degree constraints into an encoding; neither
spectral mixing nor constant-leverage Brascamp--Lieb supplies that
conversion.
