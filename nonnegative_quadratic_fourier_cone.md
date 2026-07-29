# Nonnegative quadratic Fourier cone and signed-Johnson audit

## 1. Normalization and exact graph-theoretic identification

Let \(A=(a_{ij})\) be symmetric, zero diagonal, and
\(a_{ij}\in\{\pm1\}\).  Use the doubled energy

\[
q_A(x)=x^\top A x=2\sum_{i<j}a_{ij}x_ix_j,
\qquad
R=\max_x|q_A(x)|.
\]

After changing \(A\) to \(-A\), if necessary, and switching vertices,
assume

\[
q_A(\mathbf 1)=R.
\]

Then the two degree-two Walsh polynomials

\[
f_\pm(x)=\frac{R\pm q_A(x)}4
\tag{1.1}
\]

are nonnegative on the whole cube.  Their constant Fourier coefficient is
\(R/4\), and every nonconstant Fourier coefficient has magnitude \(1/2\).
Also

\[
f_-(x)=\sum_{\substack{i\in S\\j\notin S}}a_{ij},
\qquad S=\{i:x_i=-1\}.
\tag{1.2}
\]

Thus this Fourier cone is exactly the cone of complete signed graphs whose
every signed cut is nonnegative, with the additional upper bound

\[
0\le f_-(x)\le R/2.
\tag{1.3}
\]

There is a useful standard graph translation.  Let \(H\) have edge \(ij\)
when \(a_{ij}=+1\), and put \(L=\binom n2\).  Then

\[
|E(H)|=\frac L2+\frac R4.
\tag{1.4}
\]

Moreover, (1.2) is equivalent to

\[
e_H(S,S^c)\ge\frac12|S||S^c|
\qquad\text{for every }S.
\tag{1.5}
\]

Hence \(H\) is a graph of maximum size in its Seidel-switching class
(``s-maximal'' in the graph-theory literature).  The second half of the
absolute bound gives the exact two-sided strengthening

\[
\boxed{
0\le
e_H(S,S^c)-\frac12|S||S^c|
\le \frac R4
\qquad\text{for every }S.
}
\tag{1.6}
\]

Equivalently, the original problem asks how small the maximum deviation
from \(L/2\) can be over all edge counts in one Seidel-switching class.
This is not merely analogous to switching-minimal graphs; it is exactly
that problem with both endpoints controlled.

## 2. What low-degree positivity gives, and why it stops at order \(n\)

Pointwise nonnegativity of (1.1) says that the Fourier convolution matrix

\[
\bigl(\widehat f_\pm(\alpha+\beta)\bigr)_{\alpha,\beta\in\mathbb F_2^n}
\]

is positive semidefinite.  Restricting to characters indexed by
\(k\)-subsets gives the signed Johnson operator

\[
(T_k)_{S,S-\{i\}+\{j\}}=a_{ij}
\]

and the endpoint hierarchy

\[
-\frac R2I\preceq T_k\preceq\frac R2I.
\tag{2.1}
\]

Entrywise flatness fixes

\[
\frac{\operatorname{tr}T_k^2}{\binom nk}=k(n-k).
\tag{2.2}
\]

Equations (2.1)--(2.2) imply only

\[
R\ge 2\sqrt{k(n-k)}=O(n).
\]

Thus degree-one SOS, principal Fourier-convolution PSD, and
Frobenius/second-trace arguments miss a factor \(\sqrt n\).  This is the
same signed-Johnson obstruction already recorded in the main research
ledger; equal coefficient magnitudes do not remove it.

The first new correlated trace is

\[
\boxed{
\frac{\operatorname{tr}T_k^3}{\binom nk}
=
6\,\frac{k(n-k)}{n(n-1)}
\sum_{i<j<\ell}a_{ij}a_{j\ell}a_{\ell i}.
}
\tag{2.3}
\]

Indeed, every length-three closed walk in \(J(n,k)\) is supported on a
base triangle, and each base triangle occurs in
\(\binom{n-2}{k-1}\) states and six orientations.  Even the extreme
bound obtained from (2.3),

\[
|\operatorname{tr}T_k^3|
\le \frac R2\operatorname{tr}T_k^2,
\]

still gives only \(R=O(n)\) when the signed triangle sum is
\(\Theta(n^3)\).  Therefore a fixed number of correlated traces cannot
reach the \(n^{3/2}\) scale.  A successful trace/SOS route must use order
growing with \(n\), or an inverse theorem converting many coherent cycles
into a direct Boolean witness.

## 3. Exact structural facts from the switching-minimal cone

The complement of an s-maximal graph is s-minimal.  Known elementary
structure theorems translate as follows.

* Every vertex of \(H\) has degree at least \((n-1)/2\); this is just the
  one-spin local-field inequality.
* If \(\omega=\omega(H)\), then

  \[
  |E(H)|
  \ge
  \frac{n(n-1)+\omega(\omega-1)}4.
  \]

  Combining this with (1.4) gives the exact restriction

  \[
  \boxed{R\ge\omega(H)(\omega(H)-1).}
  \tag{3.1}
  \]

* A triangle-free s-maximal graph is necessarily
  \(K_{\lfloor n/2\rfloor,\lceil n/2\rceil}\).  Its switching class also
  contains the empty graph, so its doubled absolute quadratic norm is
  \(n(n-1)\).  Thus the low-positive-energy extreme of the one-sided cone
  is maximally bad for the two-sided problem.

These facts provide a rigorous structural dichotomy but not yet a sharper
constant: a competitive signing must avoid the near-bipartite extreme,
whereas a cut-pseudorandom s-maximal graph can satisfy (1.6) with
\(R=\Theta(n^{3/2})\).  In fact (1.6) with
\(R=o(n^2)\) already forces the ordinary graphon limit of \(H\) to be the
constant \(1/2\) graphon.  Consequently every fixed subgraph density and
every fixed-order SOS/cycle statistic is asymptotically blind to the
desired \(n^{3/2}\)-order fluctuation.

The useful missing statement would have to be genuinely second-order:
either

1. triangle/cycle fluctuations of an s-maximal graph force a Boolean
   witness above the current rounding constant, or
2. small such fluctuations force \(H\) quantitatively close to a
   bipartite switching extreme, which then creates a much larger witness
   of the opposite orientation.

No such quantitative theorem was obtained here.

## 4. Hard-core-boson versus fermionic lift

The operator \(T_k\) is the \(k\)-particle hard-core-boson hopping lift of
\(A\).  The corresponding fermionic lift is the additive exterior power
\[
\mathrm d\Gamma_k(A),
\]
whose eigenvalues are

\[
\lambda_{i_1}(A)+\cdots+\lambda_{i_k}(A).
\tag{4.1}
\]

This makes the fermionic norm naturally of order \(n^{3/2}\) for a flat
spectrum.  It is tempting to compare it to \(\|T_k\|\) by a diamagnetic or
random-gauge inequality.

There is a precise obstruction.  The two lifts have the same absolute
matrix entries, but differ by the exchange signs around cycles of the
Johnson graph.  Frustrated hopping has no universal boson-below-fermion
ordering; rigorous counterexamples to such an ordering are known for
hard-core particles.  Thus the Perron--Frobenius/diamagnetic comparison is
available only in the unfrustrated case, which is exactly the rank-one-like
case where the original Boolean norm is already quadratic.

There is also a direct benchmark against the hoped-for constant.  For
symmetric Paley conference matrices, exact sparse diagonalization of
\(T_{\lfloor n/2\rfloor}\) gives:

| \(n\) | \(\|T_{n/2}\|/n^{3/2}\) | \(2\|T_{n/2}\|/n^{3/2}\) | \(\|T_{n/2}\|/\|\mathrm d\Gamma_{n/2}(A)\|\) |
|---:|---:|---:|---:|
| 6 | 0.2882279184 | 0.5764558367 | 0.6314757303 |
| 14 | 0.3020338379 | 0.6040676759 | 0.6268706527 |
| 18 | 0.3058612597 | 0.6117225195 | 0.6294572795 |

For \(n=14\), checking every \(k\le n/2\) shows that \(k=n/2\) is the
best sector, with normalized doubled bound \(0.6040676759\).
These are certified numerical eigenvalues of explicitly constructed
Paley matrices, not an asymptotic theorem.  Nevertheless they decisively
show that a direct many-particle norm bound, or a constant-fraction
fermionic comparison near the observed ratio \(0.63\), cannot improve the
already proved doubled constant

\[
c_*=0.6729867\ldots.
\]

Any useful bosonic hierarchy would have to combine several correlated
states/traces with the classical endpoint information; the bare operator
norm is weaker on the canonical conference benchmark.

## 5. Status and surviving target

No improvement of \(c_*\) was proved by this route.  The exact reductions
above rule out the following as standalone arguments:

1. degree-one SOS or PSD plus Frobenius norm;
2. any fixed signed-Johnson trace;
3. a plain diamagnetic comparison with the fermionic exterior-power lift;
4. the bare half-filled hard-core-boson operator norm.

The strongest surviving formulation is a **switching-cone fluctuation
dichotomy**.  Starting from the exact band (1.6), one must prove that an
s-maximal graph with all switching-class edge counts in a band of width
\(R/2\) either has a coherent growing-order cycle profile, detectable by a
correlated SOS certificate, or is close enough to a switching-extreme
join/bipartite structure to force an opposite-orientation energy larger
than \(R\).  Fixed graphon, fixed-cycle, and uncorrelated PSD information
cannot establish this dichotomy.

## Reference pointers

* S. Kozerenko, *On graphs with maximum size in their switching classes*
  (2015), for s-maximal graph structure.
* W. Nie, H. Katsura, and M. Oshikawa, *Ground-State Energies of Spinless
  Free Fermions and Hard-core Bosons* (2013), for the limitation of
  boson--fermion energy ordering under frustrated hopping.

