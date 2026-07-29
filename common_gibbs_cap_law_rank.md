# A common Gibbs cap law and the rank-growth obstruction

## Status

This note addresses the consistency gap in the blockwise dual measures.
The separate laws \(\mu_T\) can be replaced by one law which is
simultaneously balanced on every sparse edge block.  The price is that
the signing is an \(o(n^{3/2})\)-near-minimizer rather than necessarily
an exact minimizer.

The common-law theorem is exact and useful.  It does not by itself
force \(\Theta(n)\) tangent rank.  An explicit Boolean-cloud model shows
that exact cap constancy plus the resulting coordinatewise balance can
have rank only \(O(\lambda^{-1})\).  The cut triangle identities exclude
the literal model; a Welch bound quantifies part of the extra rigidity.

Throughout,
\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
W(A)=\frac{\max H_A-\min H_A}{2}.
\]
For a pair \(p=(x,y)\), put
\[
s_e(p)=\frac{x_ix_j-y_iy_j}{2},\qquad
S_A(p)=\sum_ea_es_e(p).
\]
Then \(\max_pS_A(p)=W(A)\).

## 1. A common finite-temperature minimizer

For \(\lambda>0\), define
\[
\mathcal Z_A(\lambda)
=
\sum_{x,y}\exp\{\lambda S_A(x,y)\}.
\tag{1.1}
\]
Choose a signing \(A_\lambda\) minimizing \(\mathcal Z_A(\lambda)\)
over all order-\(n\) signings, and let
\[
\mu_\lambda(p)
=
\frac{\exp\{\lambda S_{A_\lambda}(p)\}}
{\mathcal Z_{A_\lambda}(\lambda)}.
\tag{1.2}
\]

### Theorem 1.1: simultaneous edge balance

For every edge \(e\),
\[
\boxed{
\left(a_e\mathbb E_{\mu_\lambda}s_e\right)_+
\le \tanh\lambda\le\lambda .
}
\tag{1.3}
\]
Consequently the same single law satisfies, for every edge block \(T\),
\[
\boxed{
\sum_{e\in T}
\left(a_e\mathbb E_{\mu_\lambda}s_e\right)_+
\le \lambda |T|.
}
\tag{1.4}
\]

#### Proof

Flip only the coefficient \(a_e\).  If
\(u=a_es_e\in\{-1,0,1\}\), minimality of (1.1) gives
\[
1\le
\frac{\mathcal Z_{A_\lambda^{(e)}}(\lambda)}
{\mathcal Z_{A_\lambda}(\lambda)}
=
\mathbb E_{\mu_\lambda}e^{-2\lambda u}.
\tag{1.5}
\]
Writing \(p_+=\Pr(u=1)\) and \(p_-=\Pr(u=-1)\), (1.5) is equivalent to
\[
p_-\ge e^{-2\lambda}p_+.
\tag{1.6}
\]
Therefore
\[
\mathbb Eu=p_+-p_-
\le
\frac{1-e^{-2\lambda}}{1+e^{-2\lambda}}(p_++p_-)
\le\tanh\lambda.
\]
This proves (1.3), and summation proves (1.4). \(\square\)

This is stronger than merely coupling a prescribed finite collection
of sparse graphings: one law works for all \(2^{\binom n2}\) edge
blocks simultaneously.

## 2. The law lives on a common near-cap shell

Let \(W_n=\min_AW(A)\).  Since there are \(4^n\) ordered spin pairs,
\[
e^{\lambda W(A_\lambda)}
\le\mathcal Z_{A_\lambda}(\lambda)
\le4^ne^{\lambda W_n}.
\]
Hence
\[
\boxed{
W(A_\lambda)
\le W_n+\frac{2n\log2}{\lambda}.
}
\tag{2.1}
\]

Put
\[
\delta(p)=W(A_\lambda)-S_{A_\lambda}(p).
\]
After extracting the factor \(e^{\lambda W(A_\lambda)}\) from
(1.2), the ground term in the denominator is at least one, so
\[
\boxed{
\mu_\lambda\{\delta\ge t\}
\le4^ne^{-\lambda t}.
}
\tag{2.2}
\]

Take
\[
\lambda=\frac{b_n}{\sqrt n},\qquad
b_n\longrightarrow\infty,\qquad b_n=o(\sqrt n),
\tag{2.3}
\]
and
\[
t_n=\frac{3n\log2}{\lambda}
=\frac{3\log2}{b_n}n^{3/2}.
\tag{2.4}
\]
Then \(A_\lambda\) is an \(o(n^{3/2})\)-near-minimizer, \(t_n\) is
\(o(n^{3/2})\), and the bad-shell probability in (2.2) is at most
\(2^{-n}\).

Condition \(\mu_\lambda\) on \(\{\delta<t_n\}\), and call the resulting
law \(\nu_n\).  Since \(|a_es_e|\le1\), conditioning changes every
coordinate mean by \(O(2^{-n})\).  Thus
\[
\boxed{
\begin{aligned}
\operatorname{supp}\nu_n
&\subseteq
\{p:S_{A_\lambda}(p)\ge W(A_\lambda)-t_n\},\\
\sum_{e\in T}
\left(a_e\mathbb E_{\nu_n}s_e\right)_+
&\le
\left(\frac{b_n}{\sqrt n}+O(2^{-n})\right)|T|
\quad\text{for every }T.
\end{aligned}
}
\tag{2.5}
\]

In particular, every \(O(n)\)-edge compatible graphing has total
aligned gradient \(o(n)\), under one and the same
\(o(n^{3/2})\)-cap law.  This removes the mutually singular
\(\mu_T\) obstruction from the blockwise minimax theorem.

## 3. Why coordinate balance alone does not imply tangent rank

The following finite model obeys the two scalar conclusions one would
like to exploit: every state has exactly the same positive cap energy,
and every coordinate of the common mean is small and aligned.  Its
affine rank is nevertheless small.

Let \(k\ge4\) be even, and let the columns of a \(k\times m\) sign
matrix consist, with equal multiplicity, of every vector
\(\omega\in\{\pm1\}^k\) having exactly \(k/2+1\) positive entries.
Let \(w^{(1)},\ldots,w^{(k)}\) be its rows and let \(\nu\) be uniform
on these rows.  Then, for every state \(r\) and every coordinate \(e\),
\[
\boxed{
\sum_{e=1}^m w^{(r)}_e=\frac{2m}{k},
\qquad
\mathbb E_\nu w_e=\frac2k.
}
\tag{3.1}
\]
Yet
\[
\boxed{
\dim\operatorname{aff}\{w^{(1)},\ldots,w^{(k)}\}\le k-1.
}
\tag{3.2}
\]

Thus at balance scale \(\lambda\), taking \(k\asymp\lambda^{-1}\)
gives an exact positive cap, coordinate means \(O(\lambda)\), and
rank only \(O(\lambda^{-1})\).  For
\(\lambda=b_n/\sqrt n\), this is \(O(\sqrt n/b_n)\), not
\(\Theta(n)\).

This model is not a cut-feature cloud: its columns were prescribed
independently and need not satisfy the triangle identities.  It is an
exact convex obstruction to every proof which uses only cap constancy
and the common first moments (1.3).

### 3.1 A cut-compatible \(\sqrt n\)-rank obstruction

The triangle identities alone still do not force linear rank.  Let
\(n=q^2\), partition the vertices into \(q\) groups
\(G_1,\ldots,G_q\) of size \(q\), and put
\[
x^{(j)}=\mathbf1,\qquad
y^{(j)}_i=
\begin{cases}
-1,&i\in G_j,\\
+1,&i\notin G_j.
\end{cases}
\]
The corresponding pair features are
\[
s^{(j)}=\mathbf1_{\delta(G_j)},\qquad
r^{(j)}=\mathbf1-\mathbf1_{\delta(G_j)}.
\tag{3.3}
\]
They obey all exact endpoint identities
\[
s_er_e=0,\qquad s_e^2+r_e^2=1,
\]
and are genuine cut features.  Moreover,
\[
\|s^{(j)}\|_1=q(n-q)=n^{3/2}-n,
\tag{3.4}
\]
while, under the uniform law on \(j\),
\[
\mathbb Es^{(j)}_e=
\begin{cases}
2/q=2/\sqrt n,&
e\text{ joins two different groups},\\
0,&e\text{ is internal to one group}.
\end{cases}
\tag{3.5}
\]
Nevertheless,
\[
\boxed{
\dim\operatorname{span}\{s^{(j)},r^{(j)}:1\le j\le q\}
\le q+1=\sqrt n+1.
}
\tag{3.6}
\]

Each \(s^{(j)}\) is also a nonnegative point of the dual cut cone and
has the desired \(n^{3/2}\) total mass.  Thus even the following data
do not imply \(\Theta(n)\) rank:

* exact cut-triangle compatibility;
* exact complementary midpoint features;
* equal \(n^{3/2}\)-scale profile mass;
* coordinatewise common-law balance \(O(n^{-1/2})\);
* membership of every gradient in the dual cut cone.

For the all-positive fixed signing these profiles are not global cap
pairs: a balanced cut has quadratic-scale score.  Hence global
near-activity for one competitive signing is precisely the hypothesis
which must eliminate this obstruction.

### 3.2 Global activity kills the partition version

There is an exact lemma explaining why the preceding grouped cloud
cannot be the cap cloud of one signing when its pairs share a top
endpoint.

Switch a top state \(x\) of \(A\) to \(\mathbf1\), and call the
switched signing \(B\).  Let \(D_1,\ldots,D_k\) be disjoint vertex
sets, and suppose
\[
B\cdot\delta(D_j)\ge L\qquad(1\le j\le k).
\tag{3.7}
\]
If the \(D_j\)'s partition all vertices, then
\[
\sum_jB\cdot\delta(D_j)
=2\sum_{\substack{i<j\\i,j\text{ in different }D\text{-blocks}}}B_{ij}.
\]
A uniformly random union of the blocks therefore has expected signed
cut weight one quarter of the left side.  Since every signed cut
weight is at most \(W(A)\),
\[
\boxed{kL\le4W(A).}
\tag{3.8}
\]

The same conclusion holds when the \(D_j\)'s do not cover all
vertices.  Indeed, the extra block \(R=[n]\setminus\bigcup_jD_j\)
contributes a boundary term
\(B\cdot\delta(\bigcup_jD_j)\ge0\), because \(\mathbf1\) is a top
state.  The random-union expectation is then at least one quarter of
\(\sum_jB\cdot\delta(D_j)\).

Thus only \(O(W/L)\) disjoint difference cuts sharing one top endpoint
can have pair score at least \(L\).  In particular, only \(O(1)\) such
cuts can be near-width profiles.  Any genuine low-rank obstruction
must rotate its endpoint gauge; it cannot be the grouped-star model in
one common top gauge.

The same proof gives the useful near-top and packing forms.  Suppose
\[
H_A(x)\ge \max H_A-\alpha
\]
and \(D_1,\ldots,D_k\) are disjoint.  In the gauge switched by \(x\),
the boundary of their union is at least \(-\alpha/2\), while every cut
weight is at most \(W(A)\).  Therefore
\[
\boxed{
\sum_{j=1}^k
\frac{H_A(x)-H_A(x^{D_j})}{2}
\le4W(A)+\frac{\alpha}{2}.
}
\tag{3.9}
\]
In particular, if every displayed pair score is at least \(L>0\),
\[
\boxed{
\nu(\mathcal D)
\le\frac{4W(A)+\alpha/2}{L},
}
\tag{3.10}
\]
where \(\nu(\mathcal D)\) is the maximum number of pairwise disjoint
members of the difference-set family \(\mathcal D\).

For a genuine near-width family \(L=W-o(W)\) and
\(\alpha=o(W)\), its matching number is at most \(4+o(1)\).
If every difference set has size at most \(s\), a maximal matching
therefore supplies a vertex hitting set of size
\[
\boxed{O(sW/L).}
\tag{3.11}
\]
Consequently a common-endpoint cap family made of \(o(n)\)-vertex
differences has an \(o(n)\)-vertex transversal: every difference set
meets one fixed exceptional set of size \(o(n)\).  This is a genuine
localization handle, but it does **not** say that the differences are
contained in that exceptional set.

This common-endpoint lemma actually applies to an entire one-sided cap
family without any entropy loss.  Fix one exact top state \(x^+\).
Every negative cap state \(y\) gives
\[
\frac{H_A(x^+)-H_A(y)}2=W(A)-o(n^{3/2}),
\]
so all of its difference cuts from \(x^+\) form a family of matching
number at most \(4+o(1)\).  Symmetrically, all positive cap states can
be gauged from one exact bottom state.

The obstruction to combining this with Theorem 1.1 is more precise
than endpoint rotation.  The pair Gibbs law factorizes into a positive
temperature law and a negative temperature law, and
\[
\mathbb E s_{ij}
=\frac12\left(
\mathbb E_+x_ix_j-\mathbb E_-y_iy_j
\right).
\tag{3.12}
\]
The coordinate balance (1.3) controls only this **difference** of the
two one-sided correlation clouds.  Either one-sided cloud can have
large coordinate biases, cancelled by the other.  Fixing \(x^+\) makes
the negative-cap matching theorem available but destroys the common
first-moment balance.  Thus the exact surviving obstruction is
top--bottom marginal cancellation, not merely a changing endpoint
gauge.

## 4. What a stronger cut hypothesis adds: a Welch bound

There is a clean partial rank theorem when the common law is uniform
on \(k\) ordinary (same-orientation) cut features
\[
v^{(r)}_{ij}=x_i^{(r)}x_j^{(r)}.
\]
Put
\[
\tau_i=(x_i^{(1)},\ldots,x_i^{(k)})\in\{\pm1\}^k.
\]
Then
\[
m_{ij}:=\frac1k\sum_{r=1}^k v^{(r)}_{ij}
=\frac1k\langle\tau_i,\tau_j\rangle.
\tag{4.1}
\]
If
\[
|m_{ij}|\le\eta\qquad(i\ne j),
\tag{4.2}
\]
the Gram matrix \(G=(k^{-1}\langle\tau_i,\tau_j\rangle)\) has
\(\operatorname{rank}G\le k\), diagonal one, and
\[
\operatorname{tr}G^2
\le n+n(n-1)\eta^2.
\]
On the other hand,
\(\operatorname{tr}G^2\ge(\operatorname{tr}G)^2/k=n^2/k\).
Therefore
\[
\boxed{
k\ge\frac{n}{1+(n-1)\eta^2}.
}
\tag{4.3}
\]

At \(\eta=b/\sqrt n\), this forces \(k\ge(1+o(1))n/(1+b^2)\).
For fixed \(b\) this is linear rank.  For the
\(b_n\to\infty\) regime needed in (2.1), it gives only
\(n/b_n^2\).

There are two further gaps between (4.3) and the actual common Gibbs
law:

1. (1.3) controls only the aligned positive part, not
   \(|\mathbb E v_{ij}|\);
2. top and bottom orientations produce an indefinite signed Gram
   mixture rather than the positive semidefinite Gram matrix in
   (4.1).

The unrestricted model in Section 3 shows that these are substantive,
not cosmetic, unless the cut-triangle and competitive-scale
hypotheses are used.

### 4.1 Asymmetric temperature exposes the cancellation set

The indefinite top--bottom cancellation can be made quantitative.
For \(0<\theta<1\), put
\[
J_\theta(A)=\theta P(A)+(1-\theta)Q(A),
\]
where \(P=\max H_A\) and \(Q=-\min H_A\), and define
\[
\mathcal Z_{A,\theta}(\lambda)
=
\sum_{x,y}
\exp\{\lambda[\theta H_A(x)-(1-\theta)H_A(y)]\}.
\tag{4.4}
\]
Choose \(A\) minimizing (4.4), and use its Gibbs law.  Let
\[
C_{ij}
=
\theta\mathbb E(x_ix_j)
-(1-\theta)\mathbb E(y_iy_j).
\tag{4.5}
\]
Then \(C\) is symmetric,
\[
C_{ii}=2\theta-1,
\tag{4.6}
\]
and single-edge minimality plus Hoeffding's lemma gives
\[
\boxed{a_{ij}C_{ij}\le\lambda\quad(i\ne j).}
\tag{4.7}
\]
Indeed, for
\[
u=a_{ij}[\theta x_ix_j-(1-\theta)y_iy_j]\in[-1,1],
\]
one has
\[
1\le\mathbb Ee^{-2\lambda u},
\qquad
\log\mathbb Ee^{-2\lambda u}
\le-2\lambda\mathbb Eu+2\lambda^2.
\]

This minimizer remains competitive:
\[
\boxed{
J_\theta(A)
\le M_n+\frac{2n\log2}{\lambda},
}
\tag{4.8}
\]
because an absolute minimizer is a competitor and
\(J_\theta\le\max(P,Q)\).  The Gibbs expected weighted energy
\[
L=\sum_{i<j}a_{ij}C_{ij}
\]
satisfies
\[
L\ge J_\theta(A)-\frac{2n\log2}{\lambda}.
\tag{4.9}
\]

Put \(g_{ij}=a_{ij}C_{ij}\), and assume \(L\ge0\).  From (4.7),
\[
\sum(g_{ij})_+\le\lambda\binom n2,
\qquad
\sum(-g_{ij})_+
\le\lambda\binom n2-L.
\tag{4.10}
\]
For \(0<\lambda\le1\),
\[
\sum_{i<j}C_{ij}^2
=\sum_{i<j}g_{ij}^2
\le2\lambda\binom n2.
\tag{4.11}
\]
The trace--Frobenius rank bound therefore yields the exact estimate
\[
\boxed{
\operatorname{rank}C
\ge
\frac{n^2(2\theta-1)^2}
{n(2\theta-1)^2+4\lambda\binom n2}.
}
\tag{4.12}
\]
At \(\lambda=b/\sqrt n\) and fixed
\(\theta\ne1/2\), this is \(\Omega_\theta(\sqrt n/b)\), but not yet
\(\Theta(n)\).

Equations (4.10)--(4.12) identify the precise loss.  Away from a set
carrying total negative aligned mass \(O(\lambda n^2)\), all
correlations are bounded above by \(\lambda\).  The possible large
negative entries contribute \(O(\lambda n^2)\) to the Frobenius
square and collapse the rank bound from linear to square-root scale.
A theorem showing that global near-minimality forbids this sparse
orientation-cancellation set, or distributes it at amplitude
\(O(\lambda)\), would upgrade (4.12) to \(\Theta(n)\).

### 4.2 The dual cut cone permits the full cancellation loss

Switching a coefficient signing by a vertex cut preserves every
partition function.  Thus for the common Gibbs mean \(g_{ij}=a_{ij}C_{ij}\),
Jensen's inequality gives
\[
\sum_{ij\in\delta(S)}g_{ij}\ge0
\quad\text{for every }S,
\tag{4.13}
\]
so \(g\) lies in the dual cut cone.  This extra fact still does not
improve (4.11).

Indeed, let \(G\) be any \(d\)-regular graph and define
\[
g_{ij}=
\begin{cases}
-1,&ij\in E(G),\\
\lambda,&ij\notin E(G).
\end{cases}
\tag{4.14}
\]
If
\[
d\le\frac{\lambda n}{2(1+\lambda)},
\tag{4.15}
\]
then \(g\in\operatorname{CUT}_n^*\).  To see this, take
\(|S|=k\le n/2\).  Since
\[
e_G(S,S^c)\le dk
\le\frac{\lambda}{1+\lambda}k(n-k),
\]
one has
\[
g\cdot\delta(S)
=\lambda k(n-k)-(1+\lambda)e_G(S,S^c)\ge0.
\tag{4.16}
\]

At \(\lambda=b/\sqrt n\), take
\(d\asymp b\sqrt n\).  Then
\[
\sum_e g_e=\Theta(bn^{3/2}),
\qquad
\sum_e(g_e)_-=\Theta(bn^{3/2}),
\qquad
\sum_e(g_e)_-^2=\Theta(bn^{3/2}).
\tag{4.17}
\]
Thus coordinatewise upper balance, competitive total mass, and every
dual-cut inequality coexist with the full sparse-negative Frobenius
loss.  The missing input cannot be first-moment cut-cone geometry; it
must use that \(C\) is a difference of two actual cap correlation
matrices for one near-minimizing signing.

### 4.3 Even Boolean correlation realizability permits square-root rank

The preceding obstruction can be realized as the difference of two
genuine Boolean correlation matrices.  Fix
\(\theta\in(1/2,1)\) and put
\[
d_0=2\theta-1.
\]
Partition the vertices into \(q=n/s\) groups of size \(s\).  For
\(0\le\rho\le1\), let \(R(\rho)\) be the correlation matrix with
\[
R(\rho)_{ij}=
\begin{cases}
1,&i,j\text{ lie in the same group},\\
\rho,&i,j\text{ lie in different groups}.
\end{cases}
\tag{4.18}
\]
This is a Boolean correlation matrix: with probability \(\rho\), give
all groups one common random sign; with probability \(1-\rho\), give
the groups independent random signs; all vertices in a group inherit
its sign.

Set
\[
C=\theta R(\rho)-(1-\theta)R(0),
\qquad
\lambda_0=\theta\rho.
\tag{4.19}
\]
Then
\[
C_{ij}=
\begin{cases}
d_0,&i,j\text{ lie in the same group},\\
\lambda_0,&i,j\text{ lie in different groups},
\end{cases}
\qquad
\operatorname{rank}C\le q.
\tag{4.20}
\]
Choose \(a_{ij}=-1\) within groups and \(a_{ij}=+1\) across groups.
The aligned matrix \(g_{ij}=a_{ij}C_{ij}\) is
\(-d_0\) within groups and \(+\lambda_0\) across them.  If
\[
s\le
\frac{\lambda_0 n}{2(d_0+\lambda_0)},
\tag{4.21}
\]
then \(g\in\operatorname{CUT}_n^*\).  Indeed, for
\(|S|=k\le n/2\), the number of within-group edges crossing \(S\) is
at most \(sk\), and hence
\[
g\cdot\delta(S)
\ge
\lambda_0 k(n-k)-(d_0+\lambda_0)sk
\ge0.
\tag{4.22}
\]

Take \(s\) of the order on the right side of (4.21).  Then
\[
\sum_eg_e=\Theta(\lambda_0n^2),
\qquad
\operatorname{rank}C
\le\frac ns
=O\!\left(\frac{d_0+\lambda_0}{\lambda_0}\right).
\tag{4.23}
\]
At \(\lambda_0=b/\sqrt n\), this gives competitive
\(\Theta(bn^{3/2})\) total mass but rank only
\(O_\theta(\sqrt n/b)\), matching (4.12) in order.

Thus all of the following still fail to force linear rank:

* the exact diagonal \(2\theta-1\);
* coordinate upper balance at scale \(n^{-1/2}\);
* the full dual cut cone;
* representation as a difference of two Boolean correlation
  matrices.

What the example does not provide is the decisive joint property:
the two correlation laws are not proved to be the positive and
negative near-cap Gibbs laws of the **same** competitive signing.
That one-signing cap compatibility is now the irreducible missing
hypothesis.

### 4.4 Cap-face compatibility kills this model, but not yet the class

The grouped construction (4.18)--(4.23) cannot itself be upgraded to
a same-signing cap example.  Let
\[
B_{ab}=\sum_{i\in G_a,\ j\in G_b}a_{ij}
\]
be the intergroup block sums of a hypothetical signing \(A\).  On
group-constant spin states, the nonconstant part of the energy is the
quotient quadratic form
\[
\sum_{a<b}B_{ab}z_az_b.
\tag{4.24}
\]
The law realizing \(R(0)\) is uniform on all group-sign vectors
\(z\in\{\pm1\}^q\).  Its correlation point is the barycenter, hence a
relative-interior point, of the quotient cut-correlation polytope.

If \(R(0)\) were a bottom-face correlation matrix for \(A\), every
state in this uniform representation would have to be a bottom state:
their average energy equals the minimum only when every summand equals
the minimum.  Therefore (4.24) would be constant on the Boolean cube.
Fourier independence of the quadratic characters forces
\[
\boxed{B_{ab}=0\quad(a\ne b).}
\tag{4.25}
\]
But then every group-constant correlation law, including
\(R(\rho)\), has the same expected energy.  It cannot form the
opposite cap face with a positive \(n^{3/2}\)-scale gap.

Thus same-signing cap-face compatibility genuinely removes the sharp
model above.  The proof uses more than PSD or a four-cycle inequality:
it uses that one marginal is relative interior in the relevant
quotient cut polytope.  A boundary correlation law can avoid (4.25),
and cut polytopes have many low-dimensional exposed faces.  No general
argument was found which forces one of two opposite competitive cap
faces to contain such an interior quotient.

The exact remaining face-rigidity statement is:

> If \(R_+\) and \(R_-\) lie in opposite exposed faces of the same
> signing \(A\), their weighted difference has diffuse
> \(O(n^{-1/2})\) aligned mean, and \(A\) is globally near-minimal,
> then the joint edge-feature affine span is \(\Theta(n)\), or a
> macroscopic principal block violates global minimality.

The grouped model proves that PSD, Boolean correlation realizability,
and dual-cut first moments are insufficient without this opposite-face
condition.

## 5. Precise surviving lemma

The common-law consistency problem is solved by Theorems 1.1--2.1.
The remaining task is an inverse theorem:

> Let a competitive signing admit a common
> \(o(n^{3/2})\)-cap pair law satisfying (2.5).  Use the cut triangle
> identities, two-orientation structure, and global near-minimality to
> prove that its midpoint/tangent directions have effective rank
> \(\Theta(n)\), rather than the \(O(\lambda^{-1})\) behavior of
> Section 3.

The affine positive-ground inverse results rule out the most direct
bounded-type implementation of Section 3 at competitive scale.  A
general proof still needs to control indefinite orientation
cancellation, or to show that such cancellation localizes a
macroscopic cross-energy witness which violates principal-block
minimality.
