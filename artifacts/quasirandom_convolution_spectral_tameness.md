# Quasirandom convolution removes the logarithmic spectral loss

## Status

**Proved.**  The inverse-orbit bounded-dependency construction from
`compatible_graphing_hyperfiniteness.md` is uniformly spectrally tame on
any quasirandom group family whose smallest nontrivial irreducible
dimension grows faster than \(\log |G|\).  In particular this applies to
\(\mathrm{PSL}_2(q)\).

Consequently there are symmetric flat signings whose normalized-square
correlation graph contains a bounded-degree Cayley expander **and** whose
normalized squares have uniformly bounded operator norm.  Thus adding a
uniform spectral-tameness hypothesis does not restore hyperfiniteness.

The logarithm in the earlier matrix-Bernstein audit was an artifact of
using an operator-variance bound.  Scalar Schur orthogonality gains a
factor \(1/d_\rho\), and an \(\varepsilon\)-net costs only
\(\exp(O(d_\rho))\), so the two effects balance at scale \(\sqrt{|G|}\).

## 1. Setup

Let \(G\) be a finite group of order \(N\).  Let
\[
\mathcal O=G/(g\sim g^{-1})
\]
be the inverse-orbit set.  For \(v\in\mathcal O\), define
\[
M_v(\rho)=\sum_{g\in v}\rho(g),
\tag{1}
\]
where \(\rho:G\to U(d_\rho)\) is an irreducible unitary
representation.  Every \(M_v(\rho)\) is Hermitian.

Suppose \((\eta_v)_{v\in\mathcal O}\) are centered signs with a
dependency graph of maximum degree at most \(\Gamma\).  Thus the graph
can be properly colored with
\[
\chi\le \Gamma+1
\tag{2}
\]
colors, and the signs in each color class are mutually independent.
Lift them to the inversion-symmetric group function
\[
a(g)=\eta_{[g]}.
\tag{3}
\]
Changing \(a(e)\) to zero at the end changes every Fourier block by an
operator of norm at most one.

Before that harmless change, the \(\rho\)-Fourier block is
\[
\widehat a(\rho)
=\sum_{v\in\mathcal O}\eta_v M_v(\rho).
\tag{4}
\]

## 2. The Schur-orthogonality variance gain

For unit vectors \(u,w\in\mathbb C^{d_\rho}\), Schur orthogonality gives
\[
\sum_{g\in G}|u^*\rho(g)w|^2=\frac{N}{d_\rho}.
\tag{5}
\]
Since every inverse orbit has size at most two, Cauchy--Schwarz gives
\[
\begin{aligned}
\sum_{v\in\mathcal O}|u^*M_v(\rho)w|^2
&\le
\sum_{v\in\mathcal O}|v|
\sum_{g\in v}|u^*\rho(g)w|^2\\
&\le \frac{2N}{d_\rho}.
\end{aligned}
\tag{6}
\]
The same bound holds after restricting the sum to any color class.

This is the key estimate missed by matrix Bernstein: a fixed scalar
matrix coefficient has variance \(O(N/d_\rho)\), not \(O(N)\).

## 3. A dimension-sensitive Fourier-block tail

Fix a color class \(C\) and put
\[
F_C(\rho)=\sum_{v\in C}\eta_vM_v(\rho).
\tag{7}
\]
This is Hermitian.  For each unit vector \(u\), the summands in
\(u^*F_C(\rho)u\) are independent centered real random variables.
Hoeffding's inequality and (6) imply
\[
\Pr\{|u^*F_C(\rho)u|>t\}
\le
2\exp\left(-\frac{d_\rho t^2}{4N}\right).
\tag{8}
\]

Let \(\mathcal N\) be a \(1/4\)-net of the unit sphere in
\(\mathbb C^{d_\rho}\), viewed as a real sphere of dimension
\(2d_\rho\).  It may be chosen with
\[
|\mathcal N|\le 9^{2d_\rho}=81^{d_\rho}.
\tag{9}
\]
For Hermitian \(F\),
\[
\|F\|_{\rm op}\le
2\max_{u\in\mathcal N}|u^*Fu|.
\tag{10}
\]
Combining (8)--(10), for every \(K>0\),
\[
\boxed{
\Pr\{\|F_C(\rho)\|_{\rm op}>K\sqrt N\}
\le
2\exp\left[
d_\rho\left(\log81-\frac{K^2}{16}\right)
\right].
}
\tag{11}
\]

By summing the \(\chi\) color classes,
\[
\|\widehat a(\rho)\|_{\rm op}
\le\sum_{C=1}^{\chi}\|F_C(\rho)\|_{\rm op}.
\tag{12}
\]
Therefore
\[
\boxed{
\Pr\{\|\widehat a(\rho)\|_{\rm op}>
\chi K\sqrt N\}
\le
2\chi\exp\left[
d_\rho\left(\log81-\frac{K^2}{16}\right)
\right].
}
\tag{13}
\]

No independence between different color classes is needed.

## 4. Uniformity over all nontrivial irreducibles

Let
\[
d_{\min}(G)=
\min_{\rho\ne\mathbf 1}d_\rho.
\tag{14}
\]
A finite group has at most \(N\) inequivalent irreducible
representations.  Choose any fixed \(K>4\sqrt{\log81}\), and set
\[
c_K=\frac{K^2}{16}-\log81>0.
\tag{15}
\]
Taking a union bound in (13) gives
\[
\Pr\left\{
\max_{\rho\ne\mathbf1}\|\widehat a(\rho)\|_{\rm op}
>\chi K\sqrt N
\right\}
\le
2\chi N e^{-c_Kd_{\min}(G)}.
\tag{16}
\]
Hence if
\[
\frac{d_{\min}(G_N)}{\log |G_N|}\longrightarrow\infty,
\tag{17}
\]
then with probability \(1-o(1)\),
\[
\max_{\rho\ne\mathbf1}
\|\widehat a(\rho)\|_{\rm op}=O_\Gamma(\sqrt N).
\tag{18}
\]
More generally it is enough that
\(d_{\min}(G_N)\ge C\log N\), with \(K\) chosen large enough in terms
of \(C\) and \(\Gamma\).

For \(G=\mathrm{PSL}_2(q)\),
\[
N\asymp q^3,\qquad d_{\min}(G)\ge (q-1)/2,
\tag{19}
\]
so (17) holds.

## 5. The trivial representation

For the trivial representation, (4) is the scalar
\[
\widehat a(\mathbf1)=\sum_{g\in G}a(g).
\tag{20}
\]
Coloring the dependency graph and applying scalar Hoeffding to each
class shows that for every fixed \(\varepsilon>0\), there is a constant
\(K_0=K_0(\Gamma,\varepsilon)\) such that
\[
\Pr\{|\widehat a(\mathbf1)|>K_0\sqrt N\}<\varepsilon
\tag{21}
\]
uniformly in \(N\).  Indeed, the weight of each orbit is at most two
and the sum of squared weights is at most \(2N\).

Combining (16) and (21), there is a fixed constant
\(C_\Gamma<\infty\) for which
\[
\Pr\left\{
\max_\rho\|\widehat a(\rho)\|_{\rm op}
\le C_\Gamma\sqrt N
\right\}
\ge1-\varepsilon-o(1).
\tag{22}
\]

Since convolution diagonalizes into the irreducible Fourier blocks,
(22) is exactly
\[
\boxed{\|A\|_{\rm op}\le C_\Gamma\sqrt N.}
\tag{23}
\]

## 6. Simultaneous nonamenability and spectral tameness

In the inverse-orbit Gaussian factor construction, each sign is the
sign of a centered Gaussian field depending only on its vertex noise
and incident edge noises.  Nonadjacent orbit vertices depend on
disjoint noises.  Hence the orbit signs have a bounded-degree
dependency graph and satisfy the hypotheses above.

The previous correlation calculation and bounded-differences estimate
show, with probability \(1-o(1)\), that for every Cayley generator
\(s\),
\[
q(s)\ge \xi>0,
\tag{24}
\]
while the number of \(h\) with \(|q(h)|>\xi\) stays bounded.  Here this
is exactly the normalized-square entry
\[
q(x^{-1}y)
=\frac{(A^2)_{x,y}}{N-1}
=C_{x,y}.
\tag{25}
\]
Thus the fixed-threshold graph of the very same normalized square
contains the prescribed bounded-degree Cayley expander.

The nontrivial Fourier-block event in (18) also has probability
\(1-o(1)\), and the trivial-block event (21) has probability at least
\(1-\varepsilon\).  Taking \(\varepsilon<1/2\), all events have
joint probability at least \(1-\varepsilon-o(1)>0\) for large \(N\).
Therefore one can choose one deterministic realization satisfying
both the expander-correlation event and the spectral event
simultaneously; no conditioning or second construction is being used.

Normalize
\[
B=\frac{A}{\sqrt{N-1}},\qquad C=B^2.
\tag{26}
\]
Then
\[
\boxed{
\|C\|_{\rm op}
=\frac{\|A\|_{\rm op}^2}{N-1}
\le O_\Gamma(1).
}
\tag{27}
\]
At the same time, the fixed-threshold graph of off-diagonal entries of
\(C\) contains a bounded-degree expander and is therefore not
hyperfinite.

## 7. Consequence for the convergence campaign

The proposed classification

> compatible flat-square-root graphings become hyperfinite once their
> normalized square roots are uniformly spectrally tame

is false.

Thus the bounded-degree compatible-graphing route cannot prove the
existence of the min--max limit by spectral tameness plus
hyperfiniteness.  Any positive theorem must use a property stronger
than

1. flat symmetric square root,
2. positive-semidefinite correlation kernel,
3. bounded fixed-threshold degree, and
4. uniformly bounded operator norm.

The construction also lies at the competitive scale, because
\(\|A\|_{\rm op}=O(\sqrt N)\) implies the spectral certificate
\[
Q(A)=\max_{x\in\{\pm1\}^N}|x^\top Ax|
\le N\|A\|_{\rm op}=O(N^{3/2}).
\tag{28}
\]
What it need not satisfy is near-optimality for the unknown sharp
constant.
