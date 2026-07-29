# Prime-order exact face covers: quotient consequences and the laminar wall

## Status

This note uses the parity sharpening from `near_cap_insertion.md` at
orders
\[
n\equiv1\pmod4.
\]
Its main outputs are:

1. a careful version of the exact edge-replacement cover, including
   the case in which only one absolute orientation is active;
2. an exact finite quotient obtained from affine bases of the two
   endpoint faces;
3. a quantitative theorem showing that the one-sided quotient must
   have at least \(\Omega(\sqrt n)\) atoms (and hence face dimension at
   least \(\frac12\log_2 n-O(1)\)) for an \(O(n^{3/2})\)-competitive
   signing;
4. a two-sided "almost endpoint-determining" alternative in the
   balanced opposite-ground branch; and
5. an explicit globally optimal order-nine counterexample showing
   that the exact face cover need not possess any laminar subcover.

Thus the parity improvement is real, but it does not by itself turn
the favorable-edge cover into a laminar principal recursion.

Throughout,
\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
p(A)=\max_xH_A(x),\qquad
\nu(A)=-\min_xH_A(x),
\]
and \(M(A)=\max\{p(A),\nu(A)\}\).

## 1. Exact replacement witnesses at \(1\bmod4\)

Let \(A\) be a globally optimal signing of order
\(n\equiv1\pmod4\), and switch and globally negate if necessary so
that
\[
H_A(\mathbf1)=p(A)=M(A)=M_n=:M.
\tag{1.1}
\]
For \(T\subseteq[n]\), put
\[
c(T)=\sum_{ij\in\delta(T)}a_{ij}.
\]
Then
\[
H_A(\mathbf1^T)=M-2c(T),\qquad 0\le c(T)\le M.
\tag{1.2}
\]
Since \(n\) is odd, \(|T|(n-|T|)\) is even.  Therefore every signed
cut sum \(c(T)\) is even.  Also \(M\) is even because
\(\binom n2\) is even.  Consequently all energies of \(A\) are
congruent to \(M\pmod4\).

Fix an edge \(e=ij\) with \(a_e=+1\), and flip only this coefficient,
obtaining \(A'\).  The known top energy drops from \(M\) to \(M-2\).
Global optimality gives \(M(A')\ge M\).

If \(H_{A'}(w)\ge M\), then
\[
H_A(w)-2w_iw_j\ge M.
\]
The upper bound \(H_A(w)\le M\), together with the mod-four
congruence, forces
\[
H_A(w)=M,\qquad w_iw_j=-1.
\tag{1.3}
\]
If \(H_{A'}(w)\le-M\), the corresponding argument gives
\[
H_A(w)=-M,\qquad w_iw_j=+1.
\tag{1.4}
\]
Indeed the alternative signs would contradict \(|H_A(w)|\le M\),
and the apparent two-unit gap is excluded modulo four.

Hence:

> **Exact replacement law.** Every \(+1\) coefficient in the active
> top gauge is either separated by another exact top ground, or is
> unseparated by an exact bottom ground.

If \(\nu(A)<M\), the second alternative is absent.  Thus every
positive edge is variable on the top face.

When a bottom ground \(z\) exists, only the edges with
\(z_iz_j=-1\) automatically become same-face difference cuts in
both alternatives: a top witness differs from \(\mathbf1\) on that
edge, while a bottom witness differs from \(z\) on that edge.
Writing
\[
S=\{i:z_i=-1\},
\]
the covered favorable graph is
\[
F=\{ij\in S\times S^c:a_{ij}=+1\},
\tag{1.5}
\]
and
\[
|F|=\frac{|S||S^c|+M}{2}.
\tag{1.6}
\]
If \(S\) is unbalanced, its larger shore is already an
\(n-o(n)\)-vertex principal descent.  If both shores are linear,
\(F\) has \(\Theta(n^2)\) edges and is covered by exact top- and
bottom-face difference cuts.

The qualification in the preceding paragraph matters.  A bottom
witness with \(w_iw_j=+1\) only gives a bottom *difference* cut
crossing \(ij\) when the chosen reference bottom ground has
\(z_iz_j=-1\).

## 2. The canonical face-signature quotient

Let \(\mathcal V_+\) be the set of top edge-correlation vectors
\((x_ix_j)_{i<j}\), and define \(\mathcal V_-\) analogously.  Write
\[
d_\pm=\dim\operatorname{aff}\mathcal V_\pm.
\tag{2.1}
\]
Choose a reference in each nonempty face and choose
\(d_\pm\) face differences spanning its affine tangent space.
Their supports are vertex cuts
\[
\delta(D_1),\ldots,\delta(D_d),\qquad d=d_++d_-.
\tag{2.2}
\]
Every edge variable in either face lies in at least one of these
cuts.

Give vertex \(i\) its signature
\[
\sigma(i)=
\bigl(1_{i\in D_1},\ldots,1_{i\in D_d}\bigr)\in\mathbb F_2^d
\tag{2.3}
\]
and partition the vertices into the nonempty signature atoms
\[
V_\sigma=\{i:\sigma(i)=\sigma\}.
\tag{2.4}
\]
There are at most \(2^d\) atoms.

### Proposition 2.1 (sign-exclusion quotient)

In the one-sided case \(\nu(A)<M\), every atom \(V_\sigma\) induces
an all-\((-1)\) principal clique.

In the two-sided branch (1.5), if
\[
a_\sigma=|V_\sigma\cap S|,\qquad
b_\sigma=|V_\sigma\cap S^c|,
\]
then every edge between \(V_\sigma\cap S\) and
\(V_\sigma\cap S^c\) has sign \(-1\).

#### Proof

Two vertices in the same atom are not separated by any basis
difference cut.  The corresponding edge coordinate is therefore
constant on both affine face hulls.  In the one-sided case every
positive edge is top-variable, so a same-atom edge cannot be
positive.  In the two-sided case every edge of \(F\) is variable in
one of the two faces, so a same-atom cross edge cannot belong to
\(F\). \(\square\)

This is an actual finite quotient, not a generic low-rank cloud:
the atoms are simultaneous coordinate types of exact endpoint
faces, and their forbidden blocks have a fixed sign.

## 3. Quantitative consequences

We use the elementary principal restriction inequality
\[
\boxed{M(A[U])\le M(A)\quad(U\subseteq[n]).}
\tag{3.1}
\]
For a fixed spin on \(U\), extend it by independent random signs on
\(U^c\).  The expected full energy is exactly the principal energy,
so the maximum absolute full energy dominates it.

### Theorem 3.1 (one-sided face dimension must grow)

Assume \(\nu(A)<M\).  If \(m_\sigma=|V_\sigma|\), then
\[
\binom{m_\sigma}{2}\le M
\tag{3.2}
\]
for every signature atom.  Consequently
\[
\#\{\text{nonempty atoms}\}
\ge
\frac{n}{\sqrt{2M}+1}
\tag{3.3}
\]
and
\[
\boxed{
d_+\ge
\log_2\frac{n}{\sqrt{2M}+1}.
}
\tag{3.4}
\]
For \(M=O(n^{3/2})\),
\[
d_+\ge\frac14\log_2n-O(1).
\tag{3.5}
\]

There is a stronger centered-energy form.  Let \(r\) be the number
of atoms.  Multiply independently a spanning collection of the
top face-difference spins.  Averaging uniformly over the generated
Boolean group kills every inter-atom correlation and leaves every
intra-atom correlation equal to one.  Proposition 2.1 therefore
gives
\[
\boxed{
\nu(A)\ge\sum_{\sigma}\binom{m_\sigma}{2}.
}
\tag{3.6}
\]
By Cauchy,
\[
\nu(A)\ge\frac12\left(\frac{n^2}{r}-n\right),
\]
so
\[
\boxed{
r\ge\frac{n^2}{2\nu(A)+n},\qquad
d_+\ge
\log_2\frac{n^2}{2\nu(A)+n}.
}
\tag{3.7}
\]
For any one-sided \(O(n^{3/2})\)-competitive signing this improves
(3.5) to
\[
\boxed{d_+\ge\frac12\log_2n-O(1).}
\tag{3.8}
\]

#### Proof

An all-negative clique of order \(m\) has norm
\(\binom m2\), proving (3.2) by (3.1).  This gives
(3.3)--(3.5).

For (3.6), let \(g_1,\ldots,g_k\) be face-difference spins whose
coordinate signatures distinguish precisely the \(r\) atoms.
For uniform \(\varepsilon\in\mathbb F_2^k\), the product
\(\prod_jg_j^{\varepsilon_j}\) has expected pair correlation one
inside an atom and zero between distinct atoms.  Its expected energy
is therefore
\(-\sum_\sigma\binom{m_\sigma}{2}\).  Some spin has at most this
energy.  Cauchy and \(r\le2^{d_+}\) finish the proof. \(\square\)

Thus global edge-flip optimality does **not** presently force a
second absolute orientation, but it does force a large exact top
face whenever that orientation is absent.

### Proposition 3.2 (two-sided mixed-atom bound)

In the balanced two-sided branch,
\[
\boxed{a_\sigma b_\sigma\le M}
\tag{3.9}
\]
for every atom.  Hence the total minority mass obeys
\[
\boxed{
\sum_\sigma\min(a_\sigma,b_\sigma)
\le2^d\sqrt M.
}
\tag{3.10}
\]
In particular, if
\[
d\le\left(\frac14-\varepsilon\right)\log_2n
\quad\text{and}\quad M=O(n^{3/2}),
\tag{3.11}
\]
then, after changing \(o(n)\) exceptional vertices, the opposite
ground shore \(S\) is a union of signature atoms.

#### Proof

On \(U=V_\sigma\cap S\) and \(V=V_\sigma\cap S^c\), the cross block
is the all-\((-1)\) matrix.  For arbitrary fixed internal spins,
reversing all signs on one shore reverses the cross energy and
leaves both internal energies unchanged.  One of the two
orientations has absolute energy at least \(a_\sigma b_\sigma\).
Use (3.1), then
\(\min(a_\sigma,b_\sigma)\le\sqrt{a_\sigma b_\sigma}\), and sum over
at most \(2^d\) atoms. \(\square\)

This is a precise low-dimensional quotient alternative.  It still
does not say that the union of atoms approximating \(S\) is itself a
zero cut: arbitrary Boolean combinations of face-difference cuts
need not remain on either endpoint face.

## 4. Exact crossing-cut disintegration

There is a sharp local statement for two crossing zero cuts, but its
payment need not have the natural scale.

Let \(R,Q\) be top-face difference cuts, so
\[
c(R)=c(Q)=0.
\]
Write their four regions as
\[
X=R\cap Q,\quad Y=R\setminus Q,\quad
Z=Q\setminus R,\quad T=[n]\setminus(R\cup Q),
\]
and write \(w_{UV}\) for the signed sum of the block between regions
\(U,V\).  Put
\[
k=\frac12c(R\triangle Q)
=-\bigl(w_{XT}+w_{YZ}\bigr)\ge0.
\tag{4.1}
\]

### Proposition 4.1 (four-region payment)

The four adjacent block sums are nonnegative and satisfy
\[
\boxed{
w_{XY}+w_{ZT}
=w_{XZ}+w_{YT}
=k.
}
\tag{4.2}
\]
More precisely,
\[
\boxed{
\begin{aligned}
2w_{XY}&=c(X)+c(Y),&
2w_{ZT}&=c(Z)+c(T),\\
2w_{XZ}&=c(X)+c(Z),&
2w_{YT}&=c(Y)+c(T).
\end{aligned}}
\tag{4.3}
\]

If \(k=0\), then all six inter-region block sums vanish and all four
corner cuts are zero.  Thus \(R,Q\) refine to an exact four-way
principal Cartesian closure.

For general \(k\), regard \(R\) as the first zero-cut principal
split.  The restrictions of the \(Q\)-ground to its two children
have respective positive-ground deficits
\[
2w_{XY},\qquad 2w_{ZT},
\]
whose sum is \(2k\), and their cross equality payment is also
\(2k\).  Hence a crossing pair has the exact dichotomy

* \(k=0\): exact uncrossing and four-way closure;
* \(0<k<L\): an \(L\)-near Cartesian refinement in both children;
* \(k\ge L\): a capped-bilinear cross payment at least \(2L\), while
  the product of the two grounds has full-energy deficit at least
  \(4L\).

#### Proof

Set
\[
x=w_{XY},\ y=w_{XZ},\ z=w_{XT},\
u=w_{YZ},\ v=w_{YT},\ t=w_{ZT}.
\]
The two zero-cut equations are
\[
y+z+u+v=0,\qquad x+z+u+t=0.
\]
Since \(z+u=-k\), these give
\[
x+t=y+v=k.
\]
Directly adding the four nonnegative corner cut values gives
\[
c(X)+c(Y)=2x,\quad c(Z)+c(T)=2t,
\]
and the other two identities in (4.3).  This proves (4.2)--(4.3).

If \(k=0\), the four adjacent sums vanish.  Then
\(z+u=0\), while
\[
c(X)=c(T)=z,\qquad c(Y)=c(Z)=u.
\]
Nonnegativity forces \(z=u=0\).  The deficit and equality-payment
statements follow by evaluating the \(Q\)-ground in the principal
split \(R\sqcup R^c\). \(\square\)

At odd order every corner cut value is even.  Thus if none of the
four corners is a new zero cut, (4.3) gives
\[
k\ge4.
\tag{4.4}
\]
This quantization is only constant, and cannot be upgraded to a
bound in terms of the four region sizes from zero-cut structure
alone.

### A scalable irreducible-crossing obstruction

Let \(n\) be odd, choose \(t\in\{\pm1\}^n\) with
\(\sum_it_i=1\), and put
\[
a_{ij}=-t_it_j.
\tag{4.5}
\]
For \(U\subseteq[n]\), writing \(T_U=\sum_{i\in U}t_i\), one has
\[
\boxed{c(U)=T_U(T_U-1)\ge0.}
\tag{4.6}
\]
Partition the vertices into four linear-sized regions
\((X,Y,Z,T)\) whose \(t\)-sums are respectively
\[
-1,\quad2,\quad2,\quad-2.
\tag{4.7}
\]
Such partitions exist for all sufficiently large odd \(n\), after
choosing the region-size parities appropriately.  Then
\[
R=X\cup Y,\qquad Q=X\cup Z
\]
are crossing zero cuts, all four corner cuts are strictly positive,
and
\[
c(R\triangle Q)=12,\qquad k=6,
\tag{4.8}
\]
independently of \(n\), although all four regions are macroscopic.

This signing is deliberately not competitive:
\[
\min_xH_A(x)=-\binom n2.
\]
It proves that any scale-relevant strengthening of Proposition 4.1
must use the \(O(n^{3/2})\) absolute-cap hypothesis or global
near-minimality, not only exact zero-cut factorization.

## 5. Laminarization is false at the first relevant computed order

The following order-nine signing is written in a positive-ground
gauge:
\[
A=
\begin{pmatrix}
0&1&1&-1&1&1&-1&1&1\\
1&0&1&1&1&1&1&-1&1\\
1&1&0&-1&1&-1&1&-1&-1\\
-1&1&-1&0&1&-1&1&-1&1\\
1&1&1&1&0&-1&1&1&1\\
1&1&-1&-1&-1&0&1&1&1\\
-1&1&1&1&1&1&0&1&-1\\
1&-1&-1&-1&1&1&1&0&-1\\
1&1&-1&1&1&1&-1&-1&0
\end{pmatrix}.
\tag{5.1}
\]
Exhaustive evaluation of the \(2^8\) switching-reduced spin states
gives
\[
p(A)=\nu(A)=12,
\tag{5.2}
\]
with nine positive and twelve negative representatives.  Since the
certified exact value is \(M_9=12\), (5.1) is a global minimizer.

Represent a cut by the bit mask of its shore not containing vertex
8.  The positive ground-difference cuts are
\[
\mathcal D_+
=
\{0,4,8,94,128,160,161,224,247\}.
\tag{5.3}
\]
Taking the bottom state of mask \(50\) as reference, the negative
ground-difference cuts are
\[
\mathcal D_-
=
\{0,11,14,15,113,99,97,67,75,156,142,140\}.
\tag{5.4}
\]
Their union contains nineteen nontrivial cuts.  Together they cover
all twenty-four positive coefficients of (5.1).  (For a generic
two-sided signing the replacement law only guarantees this statement
on the favorable block (1.5); the complete cover is an additional
feature of this example.)

However:
\[
\boxed{
\text{no laminar subfamily of }
\mathcal D_+\cup\mathcal D_-
\text{ covers all positive edges.}
}
\tag{5.5}
\]
This was checked in exact integer arithmetic by enumerating all
\(2^{19}\) subfamilies.  Exactly 792 subfamilies are pairwise
laminar after orienting every shore away from vertex 8, and their
maximum positive-edge coverage is \(23\) out of \(24\).

The verification uses only:

* the exact energy evaluation of (5.1);
* the two explicit mask lists (5.3)--(5.4);
* the four-region crossing test
  \[
  S\cap T,\quad S\setminus T,\quad T\setminus S,\quad
  [n]\setminus(S\cup T)
  \]
  (two rooted shores are laminar iff at least one region is empty);
  and
* bitwise union of their cut-edge masks.

Thus even global optimality, the relevant prime parity, two exact
opposite faces, and complete favorable-edge coverage do not force a
laminar exact-face recursion.

The obstruction already contains irreducible crossing pairs.  In
the bottom gauge, the pair of difference masks \(11,99\) has region
sizes
\[
(2,1,2,4),
\]
all four corner cut values equal to \(2\), and \(k=4\).  Thus even a
global optimum at the relevant parity can realize the smallest
nontrivial constant-order crossing payment allowed by (4.4).

## 6. What survives

The exact parity cover leaves the following rigorous trichotomy.

1. If only one absolute orientation is active, its exact face has
   affine dimension at least
   \(\frac12\log_2n-O(1)\), and its generated quotient has
   \(\Omega(\sqrt n)\) all-negative atoms.
2. If two orientations are active and their difference shore is
   \(o(n)\), there is an \(n-o(n)\) principal descent.
3. If the shore is balanced, the dense favorable block is covered
   by exact difference cuts.  Low total face dimension makes the
   shore almost a union of simultaneous face-signature atoms, while
   high dimension supplies a genuinely growing endpoint face.

The missing bridge is quantitative: logarithmic face dimension is
too small for the existing correlation-rank route, while the
order-nine example rules out replacing the exact cover by a laminar
one for free.  Any continuation must either extract energy from the
crossing pairs of zero cuts, or prove that a low-dimensional
nonlaminar signature quotient itself admits a scale-preserving
exchange.
