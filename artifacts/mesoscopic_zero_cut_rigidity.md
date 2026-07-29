# Mesoscopic zero-cut rigidity: a linear rank wall and its sharp reuse obstruction

## Status

This note attacks the prime-order exact endpoint-cut branch of the
min--max Boolean quadratic form problem.  Its main new result is that
the previously recorded \(\Omega(\sqrt n)\) ground-span cost of a
balanced exact zero cut can be upgraded to a **linear** cost.

The upgrade uses the rank of a flat cross block through
\(\gamma _2\)-duality, rather than its stable rank through the global
operator norm.  If
\[
A=\begin{pmatrix}B&C\\ C^{\mathsf T}&D\end{pmatrix}
\]
is split across \(r+s=n\) vertices, then
\[
\|C\|_{\infty\to1}\le W(A)
\quad\hbox{and}\quad
\|C\|_{\infty\to1}
\ge \frac{rs}{K_G\sqrt{\operatorname{rank}C}}.
\]
Consequently
\[
\operatorname{rank}C
\ge
\left(\frac{rs}{K_GW(A)}\right)^2.
\]
For a balanced split in a signing with \(W(A)=O(n^{3/2})\), this is
\(\Omega(n)\).

For an exact zero cut, Cartesian ground factorization gives
\[
\operatorname{rank}C
\le
\operatorname{codim}L_+(B)+\operatorname{codim}L_+(D).
\]
Thus every balanced exact zero cut consumes a linear amount of
component positive-ground-span codimension.

This is a genuine strengthening, but it does **not** prove
convergence.  A square-order Hadamard residual construction already
present in `affine_type_closure_recursion.md` has exponentially many
zero cuts, no large atom, centered width exactly
\(\frac12n^{3/2}\), and total atom ground-span codimension
\(n-\sqrt n\).  The same linear codimension budget is reused across
all its zero cuts.  It is therefore a sharp obstruction to any
argument which only sums the new costs.  Global minimality below the
spectral \(1/2\) branch, or a new exchange theorem for the residual
linear-dimensional channel, is still required.

Throughout,
\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
p(A)=\max_xH_A(x),\qquad
\nu(A)=-\min_xH_A(x),
\]
and
\[
W(A)=\frac{p(A)+\nu(A)}2.
\]
The real Grothendieck constant is denoted by \(K_G\).

## 1. A rank lower bound for every flat cross block

For an \(r\times s\) matrix \(C\), use
\[
\|C\|_{\infty\to1}
=
\max_{\substack{u\in\{\pm1\}^r\\v\in\{\pm1\}^s}}
|u^{\mathsf T}Cv|.
\]

### Lemma 1.1 (range controls every cross block)

For every signing \(A\) and every vertex partition \(U\sqcup V\),
with cross block \(C=A[U,V]\),
\[
\boxed{\|C\|_{\infty\to1}\le W(A).}
\tag{1.1}
\]

#### Proof

Fix \(u,v\) attaining the bilinear norm.  The two full Boolean states
\((u,v)\) and \((u,-v)\) have the same internal energy and their
energies differ by
\[
2u^{\mathsf T}Cv.
\]
The full energy range is \(p(A)+\nu(A)=2W(A)\), proving (1.1).
\(\square\)

The point is that (1.1) is true for every partition; exact zero-cut
structure is needed only later, when matrix rank is converted into
ground-span codimension.

### Lemma 1.2 (flat rank versus Boolean bilinear norm)

Let \(C\in\{\pm1\}^{r\times s}\) have real rank \(k\).  Then
\[
\boxed{
\|C\|_{\infty\to1}
\ge \frac{rs}{K_G\sqrt{k}}.
}
\tag{1.2}
\]

#### Proof

We use the factorization norm
\[
\gamma _2(C)
=
\inf_{c_{ij}=\langle u_i,v_j\rangle}
\left(\max_i\|u_i\|_2\right)
\left(\max_j\|v_j\|_2\right).
\]
For a rank-\(k\) matrix with \(\|C\|_{\max}\le1\),
\[
\gamma _2(C)\le\sqrt k.
\tag{1.3}
\]
For completeness, take any full rank-\(k\) factorization
\(c_{ij}=\langle u_i,v_j\rangle\) and form the symmetric convex body
\[
K=\{x:|\langle u_i,x\rangle|\le1\ \text{for every }i\}.
\]
Every \(v_j\) lies in \(K\).  After sending the John ellipsoid of
\(K\) to the Euclidean unit ball, the symmetric form of John's
theorem gives
\[
B_2^k\subseteq K\subseteq\sqrt k\,B_2^k.
\]
In the transformed factorization,
\(\max_i\|u_i\|_2\le1\) and
\(\max_j\|v_j\|_2\le\sqrt k\), proving (1.3).

Let \(\gamma _2^*\) be the dual norm.  Real Grothendieck gives
\[
\gamma _2^*(C)
\le K_G\|C\|_{\infty\to1}.
\tag{1.4}
\]
Since every entry of \(C\) has magnitude one,
\[
rs=\langle C,C\rangle
\le\gamma _2(C)\gamma _2^*(C)
\le\sqrt k\,K_G\|C\|_{\infty\to1}.
\]
This is (1.2). \(\square\)

Combining the two lemmas gives the basic rank theorem.

### Theorem 1.3 (linear flat-cross rank)

For every vertex partition \(U\sqcup V\), with
\(|U|=r\), \(|V|=s\),
\[
\boxed{
\operatorname{rank}A[U,V]
\ge
\left(\frac{rs}{K_GW(A)}\right)^2.
}
\tag{1.5}
\]

In particular, if
\[
W(A)\le C_0n^{3/2},\qquad
r=\alpha n+o(n),
\]
then
\[
\boxed{
\operatorname{rank}A[U,V]
\ge
\left(
\frac{\alpha(1-\alpha)}{K_GC_0}
\right)^2n-o(n).
}
\tag{1.6}
\]

At the spectral upper scale \(C_0=1/2\), a half--half block has rank
at least
\[
\left(\frac1{4K_G^2}-o(1)\right)n
\tag{1.7}
\]
Using the standard rigorous upper bound \(K_G<1.783\), this is at
least \(0.0786n-o(n)\).

This immediately rules out the scalable rank-one four-region
obstruction under an \(O(n^{3/2})\) cap.  More generally, every
macroscopic flat cross block must have linear rank.

## 2. Exact zero cuts convert rank into ground-span codimension

Suppose a positive ground has been switched to \(\mathbf1\), and let
\(S\) be an exact zero cut:
\[
\sum_{i\in S,j\notin S}a_{ij}=0.
\tag{2.1}
\]
Put
\[
B=A[S],\qquad D=A[S^c],\qquad C=A[S,S^c].
\]
Let
\[
d_+(E)
=
\dim\operatorname{span}
\{x\in\{\pm1\}^{|E|}:H_E(x)=p(E)\},
\]
and define
\[
\Phi(E)=|E|-d_+(E).
\tag{2.2}
\]
Switching a principal signing only diagonally changes its ground
vectors, so \(d_+\) and \(\Phi\) are switching invariant.

Exact zero-cut factorization says
\[
x^{\mathsf T}Cy=0
\quad
\text{for every }x\in\mathcal G_+(B),\
y\in\mathcal G_+(D).
\]
Therefore
\[
\operatorname{rank}C\le\Phi(B)+\Phi(D).
\tag{2.3}
\]

### Corollary 2.1 (linear ground-span payment)

Every exact zero-cut split satisfies
\[
\boxed{
\Phi(A[S])+\Phi(A[S^c])
\ge
\left(
\frac{|S|(n-|S|)}{K_GW(A)}
\right)^2.
}
\tag{2.4}
\]

Thus a balanced exact zero cut in a competitive signing forces a
linear combined component-ground-span codimension.  This upgrades
the earlier Frobenius/operator-norm estimate
\[
\Phi(A[S])+\Phi(A[S^c])
\ge
\frac{|S|(n-|S|)}{\|A\|_{\mathrm{op}}^2},
\]
which gave only \(\Omega(\sqrt n)\) under the general
\(\|A\|_{\mathrm{op}}=O(n^{3/4})\) bootstrap.

For an exact conference matrix, the older estimate is already
linear because \(\|A\|_{\mathrm{op}}^2=n-1\):
\[
\Phi(A[S])+\Phi(A[S^c])
\ge\frac{|S|(n-|S|)}{n-1}.
\tag{2.5}
\]
Hence the new theorem is consistent with the conference benchmark
and extends its linear conclusion to every \(O(n^{3/2})\)-competitive
signing.

## 3. The three Klein-four payments of a crossing pair

Let \(x_0,x_1,x_2\) be three positive grounds.  Partition the
vertices according to the two relative signs
\[
(x_0x_1,\ x_0x_2)\in\{\pm1\}^2
\]
and denote the four cell sizes by
\[
n_{00},n_{10},n_{01},n_{11}.
\]
The three pairwise ground differences give the three nontrivial
two--two pairings of these four cells.  Each is an exact zero cut
after switching by one of the two relevant grounds.

Let their shore sizes be \(s_1,s_2,s_3\), and put
\[
\mathcal P(x_0,x_1,x_2)
=
\sum_{h=1}^3
\bigl[
\Phi(A[S_h])+\Phi(A[S_h^c])
\bigr].
\tag{3.1}
\]

### Proposition 3.1 (Klein-four codimension inequality)

\[
\boxed{
\mathcal P(x_0,x_1,x_2)
\ge
\frac{\sum_{h=1}^3s_h^2(n-s_h)^2}
{K_G^2W(A)^2}.
}
\tag{3.2}
\]
Moreover,
\[
\sum_{h=1}^3s_h(n-s_h)
=
n^2-\sum_{\omega\in\mathbb F_2^2}n_\omega^2.
\tag{3.3}
\]
Consequently, by Cauchy,
\[
\boxed{
\mathcal P(x_0,x_1,x_2)
\ge
\frac{
\left(n^2-\sum_\omega n_\omega^2\right)^2
}{3K_G^2W(A)^2}.
}
\tag{3.4}
\]

#### Proof

Apply (2.4) to the three pairwise difference cuts.  Every pair of
distinct Venn cells is separated by exactly two of the three
nonzero characters of \(\mathbb F_2^2\), giving (3.3).
Cauchy's inequality gives (3.4). \(\square\)

This gives a precise local dichotomy.  If one Venn cell has size
\(n-h\), then
\[
n^2-\sum_\omega n_\omega^2\ge2h(n-h),
\]
so
\[
\mathcal P(x_0,x_1,x_2)
\ge
\frac{4h^2(n-h)^2}{3K_G^2W(A)^2}.
\tag{3.5}
\]
If \(W(A)=O(n^{3/2})\) and
\(\mathcal P=o(n)\), then \(h=o(n)\).  All three ground differences
then have an \(o(n)\) shore and give \(n-o(n)\) principal
restrictions.

At the other extreme, if all four cells have size at least
\(\varepsilon n\), then
\[
n^2-\sum_\omega n_\omega^2
\ge6\varepsilon(1-2\varepsilon)n^2,
\]
and (3.4) is a linear lower bound.  For four equal cells and
\(W(A)\le(1/2+o(1))n^{3/2}\),
\[
\mathcal P(x_0,x_1,x_2)
\ge
\left(\frac{3}{4K_G^2}-o(1)\right)n
\tag{3.6}
\]
With \(K_G<1.783\), the coefficient is at least \(0.2359\).

The inequality is insensitive to the scalar four-region uncrossing
payment.  In particular, it remains linear when that payment is the
smallest parity-allowed constant.  This is exactly where the matrix
rank sees information which signed block totals miss.

## 4. Iterated zero-payment closure

Suppose a family of exact zero cuts closes under all zero-payment
unions and intersections and ends in a Boolean algebra of zero cuts.
Let its atoms be
\[
P_1,\ldots,P_\ell,\qquad m_i=|P_i|.
\]
Every union of atoms is a zero cut, every interatom signed block total
is zero, and positive grounds concatenate Cartesianly over the atoms.
Put
\[
\psi_i=m_i-d_+(A[P_i]),\qquad
\Psi=\sum_{i=1}^{\ell}\psi_i.
\tag{4.1}
\]

For a union \(U=\bigcup_{i\in I}P_i\), iterative Cartesian
factorization gives
\[
d_+(A[U])\ge\sum_{i\in I}d_+(A[P_i]).
\]
The same holds for \(U^c\), hence
\[
\Phi(A[U])+\Phi(A[U^c])\le\Psi.
\tag{4.2}
\]
Combining this upper bound with (2.4) yields:

### Theorem 4.1 (zero-algebra atom inequality)

For every union \(U\) of zero-algebra atoms,
\[
\boxed{
|U|(n-|U|)
\le K_GW(A)\sqrt{\Psi}.
}
\tag{4.3}
\]

If no atom has more than \(2n/3\) vertices, there is an atom union
with size between \(n/3\) and \(2n/3\).  Therefore
\[
\boxed{
\Psi
\ge
\frac{4n^4}{81K_G^2W(A)^2}.
}
\tag{4.4}
\]
For \(W(A)\le C_0n^{3/2}\), this is
\[
\Psi\ge
\frac{4}{81K_G^2C_0^2}n.
\tag{4.5}
\]

If the largest atom has size \(n-h>2n/3\), applying (4.3) to that
atom gives
\[
\boxed{
h(n-h)\le K_GW(A)\sqrt{\Psi}.
}
\tag{4.6}
\]
Consequently,
\[
\boxed{
W(A)=O(n^{3/2}),\ \Psi=o(n)
\quad\Longrightarrow\quad
\max_i|P_i|=n-o(n).
}
\tag{4.7}
\]

Thus iterative zero-payment refinement now has a sharp quantitative
outcome:

* sublinear total atom ground-span codimension forces an
  \(n-o(n)\) principal atom;
* absence of such an atom forces a linear residual dimension.

If the complementary shore has \(O(1)\) vertices, exact positive
maximum additivity also gives an \(O(1)\)-defect principal
puncture.  For a merely \(o(n)\) shore the statement is a
scale-preserving structural descent, not yet the summable
one-vertex recurrence needed for convergence.

## 5. Why the linear costs cannot simply be summed

### 5.1 The dimension potential is not submodular

The function
\[
U\longmapsto\Phi(A[U])
\]
is neither submodular nor supermodular.

For any two-vertex signing, each singleton has \(\Phi=0\), while the
two-vertex block has only one antipodal positive-ground pair and
\(\Phi=1\).  Thus disjoint singletons violate submodularity.

For the all-negative triangle, a two-vertex principal edge has
\(\Phi=1\), but the full triangle has six maximizing spins spanning
\(\mathbb R^3\), so \(\Phi=0\).  Taking a singleton \(U\) and the
opposite edge \(V\) gives
\[
\Phi(U)+\Phi(V)=1>
\Phi(U\cap V)+\Phi(U\cup V)=0,
\]
violating supermodularity.

Therefore ordinary Shearer or polymatroid inequalities do not
prevent the same codimension directions from paying for many
crossing cuts.

### 5.2 A sharp square-order reuse model

Let \(k\) be an even Hadamard order, \(n=k^2\), and partition the
vertices into \(k\) classes \(V_1,\ldots,V_k\), each of size \(k\).
Put \(+1\) on every internal edge.  In each \(V_a\), choose balanced
orthogonal Hadamard vectors
\[
\{v_{a,b}:b\ne a\}.
\]
For \(a\ne b\), set
\[
A[V_a,V_b]=v_{a,b}v_{b,a}^{\mathsf T}.
\tag{5.1}
\]

The exact calculation in `affine_type_closure_recursion.md` gives
\[
p(A)=\frac{n(k-1)}2,\qquad
\nu(A)=\frac{n(k+1)}2,
\]
and hence
\[
\boxed{
W(A)=\frac12n^{3/2},\qquad
M(A)=\frac12n^{3/2}+\frac n2.
}
\tag{5.2}
\]
Every type-constant spin is a positive ground.  Therefore every
union of the \(k\) types is an exact zero cut: there are exponentially
many zero cuts and no atom larger than \(\sqrt n\).

Each atom is a positive clique, whose positive-ground span has
dimension one.  Hence
\[
\boxed{\Psi=k(k-1)=n-\sqrt n.}
\tag{5.3}
\]
For a union of \(a\) types, the cross block has exact rank
\[
\boxed{\operatorname{rank}A[U,U^c]=a(k-a).}
\tag{5.4}
\]
In particular a balanced type union has cross rank \(n/4\).
Indeed, on the normalized residual vectors
\[
e_{a,b}=v_{a,b}/\sqrt k,
\]
the cross operator sends
\[
e_{b,a}\longmapsto k e_{a,b}
\]
for each ordered type pair with one endpoint in \(U\) and one in
\(U^c\), and vanishes on the remaining residual directions.  These
channels are mutually orthogonal, so the nonzero singular values of
the cross block are exactly \(k\), with multiplicity \(a(k-a)\).
This also checks
\[
\|A[U,U^c]\|_F^2
=a(k-a)k^2
=|U||U^c|.
\]

Thus the linear ground-span budget is reused across all
\(2^{\sqrt n}\) zero cuts.  The construction has the optimal
\(\sqrt n\) spectral scale and centered-width constant exactly
\(1/2\); it is not a rank-one pathology.  It proves:

> Flatness, exact Cartesian factorization, linear cross rank,
> spectral \(O(\sqrt n)\), and even the complete zero-cut Boolean
> algebra do not force a Boolean witness beyond the cap.

Any theorem eliminating this branch must use either global
minimizer exchange information not present in (5.1), or a strict
asymptotic gap below the centered spectral value \(1/2\).

## 6. Finite audit at the globally optimal order-nine example

For the globally optimal order-nine signing in
`prime_face_cover_quotient.md`, work in the bottom gauge based at
mask \(50\).  The two exact bottom-face difference cuts with masks
\[
11,\qquad99
\]
have four-region sizes
\[
(2,1,2,4).
\]
Their three pair-difference bipartitions have shore sizes
\[
3+6,\qquad4+5,\qquad3+6.
\]
Exact enumeration gives:

\[
\begin{array}{c|c|c|c}
\text{split}&\operatorname{rank}C&
\bigl(d_+(\text{left}),d_+(\text{right})\bigr)&
\Phi(\text{left})+\Phi(\text{right})\\ \hline
3+6&3&(1,1)&7\\
4+5&4&(1,1)&7\\
3+6&3&(1,1)&7.
\end{array}
\tag{6.1}
\]

Thus the three Klein-four payments are real, but all six component
ground spans are minimal antipodal lines.  The same seven-dimensional
codimension scale pays for every pairing.  This is a finite,
globally optimal example of maximal local reuse.  It also has the
smallest nontrivial scalar crossing payment allowed by parity, so
scalar uncrossing and dimension uncrossing fail for complementary
reasons.

Here is a direct reproducibility recipe.  Number vertices
\(0,\ldots,8\), let \(A\) be matrix (5.1) of
`prime_face_cover_quotient.md`, let \(z_m\) be \(-1\) on the set bits
of \(m\), and form
\[
\widetilde A=-\operatorname{diag}(z_{50})
A\operatorname{diag}(z_{50}).
\tag{6.2}
\]
For the three rows of (6.1), use respectively
\[
(\text{gauge mask},\text{shore mask})
=(0,11),\ (0,99),\ (11,104).
\tag{6.3}
\]
The resulting cross blocks, with shore vertices first in increasing
order, are
\[
\begin{aligned}
C_{11}&=
\begin{pmatrix}
-1&1&1&1&-1&-1\\
1&-1&-1&1&-1&1\\
1&1&-1&-1&1&-1
\end{pmatrix},\\[2mm]
C_{99}&=
\begin{pmatrix}
-1&1&1&-1&-1\\
1&1&-1&-1&1\\
-1&-1&1&1&1\\
-1&-1&1&-1&1
\end{pmatrix},\\[2mm]
C_{104}^{(11)}&=
\begin{pmatrix}
1&1&-1&-1&-1&1\\
-1&1&-1&1&1&1\\
-1&-1&-1&1&-1&1
\end{pmatrix}.
\end{aligned}
\tag{6.4}
\]
Exact row reduction gives ranks \(3,4,3\).  Exhaustive evaluation of
the two principal cubes gives positive maxima
\[
(3,9),\qquad(4,8),\qquad(3,9),
\tag{6.5}
\]
and exactly one antipodal maximizing pair in every component.
Equations (6.4)--(6.5) reproduce every entry of (6.1) without a
floating-point computation.

The order-six conference optimizer from
`capped_bilinear_inverse.md` is also consistent.  Its displayed
\(2+4\) zero cut has cross rank \(2\), while the two component
positive-ground spans both have dimension one, giving total
codimension \(4\).

## 7. Application to the puncture tangent torsor

The puncture construction produces full Boolean configurations
\(x_i\), one for each selected child puncture, with endpoint deficits
\(a_i\le d_i\).  Their pairwise difference shores obey the exact
torsor identity
\[
D_{ij}\triangle D_{jk}=D_{ik},
\tag{7.1}
\]
and in the gauge switched by \(x_i\), the signed traffic of \(D_{ij}\)
is
\[
\frac{a_j-a_i}{2}.
\tag{7.2}
\]
Within a common endpoint-orientation class, if all relevant
\(d_i=0\), the \(x_i\)'s are exact same-face grounds and every
\(D_{ij}\) is an exact pair-difference zero cut.

For such an exact class, (2.4) gives, for every pair,
\[
\boxed{
\Phi(A[D_{ij}])+\Phi(A[D_{ij}^c])
\ge
\left(
\frac{|D_{ij}|(n-|D_{ij}|)}{K_GW(A)}
\right)^2.
}
\tag{7.3}
\]
Hence an \(o(n)\)-shore pair gives an \(n-o(n)\) principal descent,
whereas every \(\varepsilon\)-balanced pair has a linear dimension
payment.

There is one exact limit on reuse inside the torsor.  For a balanced
pair \(i,j\), put
\[
\mathcal K_{ij}
=
\{k:x_ix_jx_k\text{ is again a positive ground}\}.
\tag{7.4}
\]
Switch by \(x_i\).  The cuts \(D_{ij}\), \(D_{ik}\), and
\(D_{ij}\triangle D_{ik}\) are all zero exactly when
\(k\in\mathcal K_{ij}\).  Zero-payment four-region refinement then
shows that the restrictions of all \(x_k\),
\(k\in\mathcal K_{ij}\), are component positive grounds on both
shores.  Therefore
\[
\boxed{
\operatorname{rank}\{x_k:k\in\mathcal K_{ij}\}
\le
n-\Phi(A[D_{ij}])-\Phi(A[D_{ij}^c]).
}
\tag{7.5}
\]

Let the whole tangent family have real rank \(r\).  Since adjoining
one excluded column raises rank by at most one,
\[
\boxed{
\#\{k:x_ix_jx_k\text{ is not a positive ground}\}
\ge
\left[
r+\Phi(A[D_{ij}])+\Phi(A[D_{ij}^c])-n
\right]_+.
}
\tag{7.6}
\]
Thus a full-rank tangent family and a balanced pair force linearly
many nonclosed cubic products.  At orders \(n\equiv1\pmod4\), every
such product has energy deficit at least four.

This still does not close the asymptotic argument.  With \(O(n^2)\)
balanced pairs, (7.6) can force \(O(n^3)\) nonclosed triples, but the
parity charge per triple may remain \(O(1)\).  After averaging over
all triples this is only a constant cubic defect, not an
\(n^{3/2}\)-scale Boolean witness.  The order-nine example realizes
constant crossing payments, so this loss is genuine at finite order.

There is also a stability wall.  If the \(a_i\)'s are merely equal
or close but positive, (7.2) gives a zero or small signed traffic
between **near** grounds, not an exact top-ground zero cut.
The range estimate (1.1) and flat rank theorem remain valid, but the
rank-to-ground-codimension implication (2.3) becomes only scalar
near-annihilation.  Turning it into approximate rank requires frame
lower bounds for the principal near-ground clouds, exactly the
stability input missing elsewhere in the campaign.

## 8. Verdict and sharpened next lemma

The proved chain is now:
\[
\begin{gathered}
\text{flat macroscopic cross block}
\Longrightarrow
\text{linear matrix rank},\\
\text{exact zero cut}
\Longrightarrow
\text{linear component-ground-span codimension},\\
\text{three crossing grounds}
\Longrightarrow
\text{three simultaneous Klein-four dimension payments},\\
\text{zero-payment Boolean closure}
\Longrightarrow
\bigl[\text{an }n-o(n)\text{ atom}\bigr]
\ \text{or}\ 
\bigl[\text{a linear residual dimension}\bigr].
\end{gathered}
\]

The second branch is sharp at centered constant \(1/2\), by the
Hadamard residual construction.  Therefore the exact missing
statement cannot be a generic rank or Shearer inequality.  A useful
completion would have to be one of the following genuinely
global-minimality statements:

1. **Sublinear residual theorem.**  In the maximal zero-payment
   closure of a prime-order global minimizer, the atom potential
   \(\Psi\) is \(o(n)\).  Theorem 4.1 would then force an
   \(n-o(n)\) principal atom.
2. **Strict residual gap.**  A linear-dimensional mesoscopic
   residual compatible with the prime exact edge-replacement cover
   has centered width at least
   \((1/2-o(1))n^{3/2}\), with equality only in a replaceable
   Hadamard-type switching class.
3. **Stable tangent factorization.**  Near-ground zero traffic in
   the puncture torsor implies approximate Cartesian annihilation
   with a frame constant strong enough to preserve the linear rank
   payment.

Any one of these would materially advance convergence.  The present
note proves the linear mesoscopic rigidity needed by all three and
identifies, by an exact infinite construction, why that rigidity
alone cannot finish the problem.
