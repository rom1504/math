# High-order flat quadratic chaos: exact reductions and barriers

## Scope and normalization

For a symmetric zero-diagonal signing \(A=(a_{ij})\), \(a_{ij}\in\{\pm1\}\), write

\[
H_A(x)=\sum_{1\le i<j\le n}a_{ij}x_ix_j,\qquad
M(A)=\|H_A\|_{L_\infty(\{\pm1\}^n)}.
\]

Thus the original quantity is \(F(n)=\min_A M(A)\).  The doubled matrix
normalization used elsewhere is

\[
Q(A)=\max_x|x^\top A x|=2M(A).
\]

This note audits whether moments of order \(k=\Theta(n)\), Boolean
Bohnenblust--Hille/Sidon inequalities, or positivity-preserving Eulerian
subfamilies can improve the current universal lower bound or prove convergence.
Claims below are marked as proved, literature, or open.

## 1. Exact convolution identity (proved)

Let \(G=\mathbb F_2^n\), and identify an edge \(\{i,j\}\) with
\(e_i+e_j\in G\).  Define

\[
f_A(e_i+e_j)=a_{ij},\qquad f_A(g)=0
\quad\text{off the weight-two sphere}.
\]

Use the unnormalized convolution

\[
(f*g)(z)=\sum_{y\in G}f(y)g(z-y)
\]

and Fourier transform

\[
\widehat f(u)=\sum_{z\in G}f(z)(-1)^{u\cdot z}.
\]

If \(x_i=(-1)^{u_i}\), then \(\widehat f_A(u)=H_A(x)\).  Therefore

\[
\boxed{\quad
\mathbb E_x H_A(x)^{2k}
=\sum_{z\in G}\left(f_A^{*k}(z)\right)^2.
\quad} \tag{1}
\]

Equivalently, if

\[
H_A(x)^k=\sum_{T\subset[n]}c_{k,T}\chi_T(x),
\]

then

\[
c_{k,T}
=\sum_{\substack{e_1,\ldots,e_k\\
e_1\triangle\cdots\triangle e_k=T}}
\prod_{\ell=1}^k a_{e_\ell},
\qquad
\mathbb E H_A^{2k}=\sum_Tc_{k,T}^2. \tag{2}
\]

This is the useful positivity: the *final Fourier coefficients* are squared.
It does not make the signed sums defining those coefficients positive.

### The fourth moment

Put \(N=\binom n2\), and let

\[
\operatorname{haf}(A_T)
=\sum_{\text{perfect matchings }P\text{ of }T}\prod_{e\in P}a_e
\]

for a four-set \(T\).  Splitting the boundary of an ordered pair of edges
into weights \(0,2,4\) gives the exact identity

\[
\boxed{\quad
\mathbb E H_A^4
=N^2
+4\sum_{i<j}(A^2_{ij})^2
+4\sum_{\substack{T\subset[n]\\|T|=4}}
\operatorname{haf}(A_T)^2.
\quad} \tag{3}
\]

Every even-order principal hafnian of a sign matrix is odd: it is a sum of
an odd number \((2k-1)!!\) of signs.  Hence it never vanishes.  At fixed
order, however, (3) and this parity fact only yield an \(O(n)\) lower bound
for the \(L_4\) norm, far below the \(n^{3/2}\) scale.

### Top Fourier layer at arbitrary order

If \(|T|=2k\), all \(k\) edges contributing to \(c_{k,T}\) must be disjoint.
Consequently

\[
\boxed{\quad c_{k,T}=k!\operatorname{haf}(A_T),\qquad |T|=2k.\quad} \tag{4}
\]

Thus

\[
\mathbb E H_A^{2k}
\ge (k!)^2
\sum_{\substack{T\subset[n]\\|T|=2k}}
\operatorname{haf}(A_T)^2
\ge (k!)^2\binom n{2k}. \tag{5}
\]

The last inequality is only the odd-parity bound.  For
\(k=\alpha n\), it has \(2k\)-th root of order \(n\), not
\(n^{3/2}\).

A natural finer object is therefore

\[
\mathcal H_{n,k}(A)
=\sum_{|T|=2k}\operatorname{haf}(A_T)^2. \tag{6}
\]

The exact recursion

\[
\operatorname{haf}(A_T)
=\frac1k\sum_{\{i,j\}\subset T}
a_{ij}\operatorname{haf}(A_{T\setminus\{i,j\}}) \tag{7}
\]

views the hafnian vectors as an orbit of a signed up-operator on the even
subset lattice.  A lower singular-growth theorem for this particular orbit
would be genuine progress.  Ordinary Johnson-scheme diagonalization does
not apply: the edge signs destroy permutation equivariance, and the
off-diagonal terms of the Gram operator contain truncated row correlations
and signed four-vertex terms.

Even the random-sign benchmark for (6) is not enough by itself.  If one
could prove uniformly that

\[
\mathcal H_{n,k}(A)\ge
e^{-\gamma k}\binom n{2k}(2k-1)!!,\qquad k=\alpha n,
\]

then the contribution (4) would give the asymptotic constant

\[
c_{\rm top}(\alpha,\gamma)
=\alpha\sqrt{2\alpha}\,
\exp\!\left(
\frac{h(2\alpha)}{2\alpha}-\frac32-\frac\gamma2
\right), \tag{8}
\]

where \(h(p)=-p\log p-(1-p)\log(1-p)\).  With \(\gamma=0\), the
maximum of (8) is only about \(0.15925\), attained near
\(\alpha=0.2914\).  Thus a random-scale hafnian-minor theorem restricted to
the top layer would still not beat the present \(0.33649\) lower constant;
one needs many boundary layers or a much stronger correlated lower bound.

## 2. A rigorous no-go for retaining “positive diagrams” (proved)

In the raw expansion of \(\mathbb E H_A^{2k}\), consider terms in which
exactly \(k\) distinct edges each occur twice.  Every such term has product
\(+1\), independent of \(A\), and their total is

\[
D_{N,k}=\binom Nk\frac{(2k)!}{2^k}. \tag{9}
\]

For \(N=\binom n2\) and \(k=\lfloor\alpha n\rfloor\), Stirling's formula
gives

\[
D_{N,k}^{1/(2k)}
=\left(\sqrt{\frac{\alpha}{e}}+o(1)\right)n^{3/2}. \tag{10}
\]

Along the infinite Paley conference orders there are signings \(C\) with
\(C^2=(n-1)I\), and hence

\[
\|H_C\|_\infty\le \frac12n\sqrt{n-1}. \tag{11}
\]

It follows from (10)--(11) that for every fixed
\(\alpha>e/4\),

\[
D_{N,k}>\mathbb E H_C^{2k}
\]

for all sufficiently large conference orders.  Therefore all the remaining
Eulerian diagrams have a *negative net contribution* which cancels the
manifestly positive exact-pairing family at exponential-in-\(n\) scale.

This rules out the tempting proof step

> expand the \(2k\)-th moment at \(k=\Theta(n)\), retain the
> coefficient-independent edge-even/pairing diagrams, and discard the rest.

The same warning applies to the full edge-even subfamily, which is the
\(2k\)-th moment of a sum of \(N\) independent signs.  Positivity appears
only after collecting all diagrams with the same boundary, as in (1), not
diagram by diagram.

## 3. Linear-order moments are already close to the original extremum (proved)

There are only \(2^{n-1}\) distinct values after quotienting by the global
flip \(x\mapsto-x\).  Hence

\[
\boxed{\quad
\|H_A\|_{2k}\le M(A)
\le 2^{(n-1)/(2k)}\|H_A\|_{2k}.
\quad} \tag{12}
\]

For \(k=\alpha n\), moment control and supremum control differ by the fixed
factor \(2^{1/(2\alpha)}+o(1)\).  For \(k/n\to\infty\), the factor is
\(1+o(1)\).  Consequently:

* fixed moments cannot see the target scale sharply;
* moments of order \(\Theta(n)\) can give constants, but require controlling
  the same signed additive energy that drives the original problem;
* superlinear moments are asymptotically equivalent to the original
  \(L_\infty\) problem rather than an evident shortcut around it.

The convolution formula (1) is therefore best viewed as a precise
reformulation and a source of structured subproblems, not by itself a
convergence mechanism.

## 4. Boolean Sidon and Bohnenblust--Hille audit

### 4.1 The exact flat-Sidon reformulation (proved)

For the full level-two support

\[
\mathcal S_2=\{\{i,j\}:1\le i<j\le n\},
\]

the unrestricted Sidon constant is

\[
\operatorname{Sid}(\mathcal B_{=2}^n)
=\sup_{c\ne0}
\frac{\sum_{i<j}|c_{ij}|}
{\left\|\sum_{i<j}c_{ij}x_ix_j\right\|_\infty}. \tag{13}
\]

Define its flat restriction by requiring \(|c_{ij}|=1\).  Then exactly

\[
\boxed{\quad
\operatorname{Sid}_{\rm flat}(\mathcal B_{=2}^n)
=\frac{\binom n2}{F(n)}.
\quad} \tag{14}
\]

Thus the original limit question is equivalent to convergence of the
normalized *flat* Sidon constant, with the reciprocal normalization

\[
\frac{\operatorname{Sid}_{\rm flat}(\mathcal B_{=2}^n)}{\sqrt n}
=\frac{1+o(1)}{2\,F(n)/n^{3/2}}. \tag{15}
\]

### 4.2 What the 2024 Boolean Sidon paper does and does not imply

Defant--Galicer--Mansilla--Mastyło--Muro study the unrestricted quantity
(13).  Their Corollaries 5.11--5.12 imply the correct order

\[
\operatorname{Sid}(\mathcal B_{=2}^n)=\Theta(\sqrt n), \tag{16}
\]

with constants depending on the fixed degree.  Their Theorem 5.1 compares
this Sidon constant, again up to a degree-dependent factor, with the
projection constant of level one.  These statements do **not** identify
the leading constant for (13), do not prove that flat coefficients are
extremal, and therefore do not determine or improve (14).

This is a clean scope distinction:

\[
\operatorname{Sid}_{\rm flat}(\mathcal B_{=2}^n)
\le \operatorname{Sid}(\mathcal B_{=2}^n),
\]

but lower bounds for the unrestricted Sidon constant need not be witnessed
by flat coefficients.  An explicit sharp upper bound

\[
\operatorname{Sid}(\mathcal B_{=2}^n)
\le (K+o(1))\sqrt n
\]

would imply \(F(n)/n^{3/2}\ge 1/(2K)+o(1)\), but to improve
\(0.3364933644\) one needs

\[
K<1.485876\ldots. \tag{17}
\]

No such explicit constant follows from that paper; its degree-dependent
comparison constants are far too large for (17).

The exact projection formula in the same paper is also a different
quantity:

\[
\lambda(\mathcal B_{=2}^n)
=\mathbb E\left|\sum_{i<j}x_ix_j\right|,
\]

which is order \(n\), whereas the unrestricted level-two Sidon constant is
order \(\sqrt n\) and the flat minimax norm is order \(n^{3/2}\).

### 4.3 Bohnenblust--Hille constants: exact threshold and a scope correction

Let \(B^{\rm tet}_{2,\mathbb R}\) be the dimension-free degree-two
homogeneous Boolean/tetrahedral Bohnenblust--Hille constant:

\[
\left(\sum_{i<j}|c_{ij}|^{4/3}\right)^{3/4}
\le B^{\rm tet}_{2,\mathbb R}
\left\|\sum_{i<j}c_{ij}x_ix_j\right\|_\infty. \tag{18}
\]

For flat coefficients, (18) would give

\[
\liminf_n\frac{F(n)}{n^{3/2}}
\ge \frac{2^{-3/4}}{B^{\rm tet}_{2,\mathbb R}}. \tag{19}
\]

It beats \(0.3364933644\) exactly if

\[
B^{\rm tet}_{2,\mathbb R}<1.767059\ldots. \tag{20}
\]

The often quoted \(1.83737\ldots\) is **not** this constant.  It is the
sharp real polynomial Bohnenblust--Hille constant for general
two-homogeneous polynomials in **two variables**, as computed by
Jiménez-Rodríguez et al.  That class permits square monomials, the result is
dimension-specific, and it supplies neither a dimension-free upper bound
nor a sharp no-diagonal/tetrahedral constant.  Even if it could be inserted
in (19), it would give only \(0.32362\ldots\), below the current bound.

Defant--Mastyło--Pérez prove that the Boolean/tetrahedral constants are
dimension-free and subexponential in the degree, but do not compute the
sharp degree-two constant.  Astashkin--Lykov likewise prove order-sharp
equivalences for second-order Rademacher chaos, with non-sharp universal
constants; their estimates do not approach (20).

For comparison, a balanced bipartition followed by sharp real Littlewood
\(4/3\) gives only

\[
M(A)\ge \frac14n^{3/2}+o(n^{3/2}) \tag{21}
\]

for flat complete-graph coefficients.  Polarization of a symmetric
quadratic form and then Littlewood is weaker still (constant
\(2^{-5/2}\approx0.17678\)).  Thus the classical bilinear inequality does
not imply a competitive tetrahedral constant.

The promising analytic target is now precise:

> Prove the weighted no-diagonal inequality (18) with constant below
> \(1.767059\), or prove that its extremizers may be taken flat and determine
> its sharp asymptotic constant.

This is broader than the original flat problem in the first formulation,
and essentially contains it in the second.

## 5. Conditional linearization does not bootstrap the greedy bound

For an \(r\times s\) cross block \(B\), condition on \(y\) and write
\(c_i=(By)_i\).  Exact positive pairings in the linear Rademacher moment give

\[
\mathbb E_x\left(\sum_i c_ix_i\right)^{2k}
\ge \frac{(2k)!}{2^k}e_k(c_1^2,\ldots,c_r^2). \tag{22}
\]

If \(0\le c_i^2\le L\) and \(\sum_i c_i^2=S\), then

\[
k!\,e_k(c_1^2,\ldots,c_r^2)
\ge \prod_{j=0}^{k-1}(S-jL), \tag{23}
\]

by expanding ordered distinct tuples and exposing one coordinate at a time.
This is rigorous but loses badly when a few row fields are large.

A Laplace/product variant leads to

\[
\max_{x,y}x^\top By
\ge r\sqrt{s}\,
\frac{\mathbb E\log\cosh(tG)}{t},
\qquad G\sim N(0,1), \tag{24}
\]

after the usual asymptotic replacement of a length-\(s\) Rademacher sum.
Optimizing \(t\), the right-hand side increases to
\(r\sqrt{s}\,\mathbb E|G|\), exactly the ordinary
greedy/Khintchine value.  Thus independent conditional linear moments do
not by themselves improve the known one-round extraction mechanism.

## 6. Surviving research targets

1. **All-boundary convolution inequality.**  Lower-bound
   \(\sum_Tc_{k,T}^2\) for \(k=\alpha n\) using simultaneous constraints
   among boundary layers.  A top-layer hafnian theorem alone is too weak.
2. **Signed up-operator orbit.**  Control the orbit (7) from
   \(h_\varnothing=1\) despite the lack of Johnson-scheme equivariance.
   A theorem only about the minimum singular value of the full operator is
   likely false or irrelevant; it must exploit the special orbit.
3. **Sharp weighted tetrahedral BH.**  Establish (18) with the explicit
   threshold (20).  Existing Boolean BH and Sidon literature supplies only
   order estimates here.
4. **Flat versus unrestricted Sidon.**  Determine whether
   \(\operatorname{Sid}_{\rm flat}/\operatorname{Sid}\) has a nonzero
   limit, or whether an unrestricted extremizer can be flattened without
   increasing its Boolean norm by more than \(1+o(1)\).
5. **Positivity after grouping, not before.**  Any Eulerian/free-energy
   approach must group diagrams into the squared boundary coefficients of
   (1), or find another manifestly positive identity.  The conference
   counterexample in Section 2 forbids termwise deletion of signed
   diagrams at the required moment order.

## 7. Bottom line

No convergence proof or improved numerical lower constant came from this
route.  It did produce two reusable rigorous conclusions:

* high-order moments are exactly a signed weight-two convolution energy,
  with hafnian principal-minor layers as concrete subobjects;
* at \(k=\Theta(n)\), positive pairing diagrams can exceed the *entire*
  moment on conference matrices, so any proof that keeps those diagrams
  and discards signed Eulerian terms is invalid at exponential scale.

The Boolean Sidon paper does not subsume the flat problem: it controls the
unrestricted level-two Sidon constant only up to fixed-degree factors,
whereas \(F(n)\) is exactly the reciprocal of its flat-coefficient
restriction.

## Primary sources

1. A. Defant, D. Galicer, M. Mansilla, M. Mastyło, S. Muro,
   “Asymptotic insights for projection, Gordon--Lewis and Sidon constants
   in Boolean cube function spaces,” IMRN (2024),
   <https://arxiv.org/abs/2302.00233>.
2. A. Defant, M. Mastyło, A. Pérez,
   “On the Fourier spectrum of functions on Boolean cubes,”
   *Mathematische Annalen* 374 (2019), 653--680,
   <https://arxiv.org/abs/1706.03670>.
3. P. Jiménez-Rodríguez, G. A. Muñoz-Fernández,
   M. Murillo-Arcila, J. B. Seoane-Sepúlveda,
   “Sharp values for the constants in the polynomial
   Bohnenblust--Hille inequality,”
   <https://arxiv.org/abs/1502.02173>.
4. S. V. Astashkin, K. V. Lykov,
   “Random unconditional convergence of Rademacher chaos in
   \(L_\infty\) and sharp estimates for discrepancy of weighted graphs and
   hypergraphs,” <https://arxiv.org/abs/2412.20107>.
5. J. E. Littlewood, “On bounded bilinear forms in an infinite number of
   variables,” *Quarterly Journal of Mathematics* 1 (1930), 164--174.
