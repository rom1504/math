# Margin-preserving refills in the affine mesoscopic branch

Checkpoint: 2026-07-26.

## 1. Scope and status

Let \(A\) be a symmetric signing of \(K_n\), with one-copy energy

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j.
\]

Suppose that an exact affine family of positive ground states has
already produced type classes \(V_\alpha\), of sizes \(k_\alpha\).
The exact closure theorem gives

\[
\sum_{i\in V_\alpha,\ j\in V_\beta}a_{ij}=0
\qquad(\alpha\ne\beta)
\tag{1.1}
\]

and, writing

\[
r_{i,\beta}=\sum_{j\in V_\beta}a_{ij},
\]

the vertexwise domination

\[
r_{i,\alpha}\ge
\sum_{\beta\ne\alpha}|r_{i,\beta}|
\qquad(i\in V_\alpha).
\tag{1.2}
\]

This note studies refills of every intertype block which preserve
**all row and column sums**, rather than only the block total.

The main proved conclusions are:

1. In the nonrecursive branch
   \(\max_\alpha k_\alpha\le(1-\delta)n\), the simultaneous
   fixed-margin fibre always contains a Boolean cube of
   dimension \(\Omega_\delta(n^2)\).  This remains true despite
   examples in which an individual block has a singleton fibre.
2. Uniform doubly balanced refills have an exact two-replica
   covariance kernel.  Its natural Gaussian/checkerboard limit has
   overlap correlation \(\rho^2\).
3. The exponential second-moment calculation is stable precisely
   up to the constant \(1/2\).  Thus a generic isotropic
   margin-preserving refill already has Boolean witnesses of size
   \((1/2-o(1))n^{3/2}\); it cannot construct a sub-\(1/2\)
   signing.
4. If an exceptional refill is spectrally flat, exact affine
   closure forces a sharp alternative.  Zero cross margins force
   type sizes \((1+o(1))\sqrt n\) and positive cap
   \((1/2-o(1))n^{3/2}\).  Arbitrary compensated margins give the
   weaker but exact floor \(1/(2\sqrt2)\), and approaching that
   floor forces nearly extreme margin channels.

No convergence theorem is claimed.  The unresolved inverse statement
is now: a sub-\(1/2\) exceptional member of the quadratic refill
fibre must either have a nonisotropic trade covariance or use
substantial, highly biased margin channels; turn either failure into
a principal quotient/descent.

## 2. What fixed-margin replacement preserves

For an \(a\times b\) intertype block \(B\), prescribe

\[
B\mathbf1_b=r,\qquad
\mathbf1_a^\mathsf TB=c^\mathsf T.
\tag{2.1}
\]

Replacing \(B\) by any other sign matrix with the same \(r,c\)
preserves:

* the total of the block, hence the energy of every type-constant
  state;
* every intertype contribution to every vertex row sum;
* the row-domination inequalities (1.2), since the internal blocks
  are held fixed;
* every singleton-flip margin at every type-constant state.

It does **not** automatically preserve global ground maximality away
from the type cube.  This distinction is important below.

## 3. The global \(L^1\) margin budget

For a block \(B_{\alpha\beta}\), put

\[
L_{\alpha\beta}
=
\sum_{i\in V_\alpha}|r_{i,\beta}|
+
\sum_{j\in V_\beta}|r_{j,\alpha}|.
\tag{3.1}
\]

Summing (1.2) over all vertices gives the exact global budget

\[
\boxed{
L:=\sum_{\alpha<\beta}L_{\alpha\beta}
\le
\sum_\alpha\sum_{i\in V_\alpha}r_{i,\alpha}
=2p(A).
}
\tag{3.2}
\]

Thus, for a competitive signing with
\(M(A)\le Cn^{3/2}\),

\[
\boxed{L\le2Cn^{3/2}.} \tag{3.3}
\]

This is much smaller than the \(\Theta(n^2)\) number of cross edges
in the nonrecursive branch.

## 4. Alternating rectangles and quadratic fibre rank

Consider a zero-total \(a\times b\) sign matrix \(B\).  Let its row
and column sums be \(r_i,c_j\), and put

\[
L_B=\sum_i|r_i|+\sum_j|c_j|.
\]

An alternating rectangle is a pair of rows \(i<i'\) and columns
\(j<j'\) on which the signs are

\[
\begin{pmatrix}+&-\\-&+\end{pmatrix}
\quad\hbox{or}\quad
\begin{pmatrix}-&+\\+&-\end{pmatrix}.
\]

Flipping all four signs preserves every row and column sum.

### Lemma 4.1 (exact rectangle formula)

For a pair of rows \(i<i'\), let \(h_{ii'}\) be their Hamming
distance.  The number \(T(B)\) of alternating rectangles satisfies

\[
\boxed{
T(B)
=\frac14\sum_{i<i'}h_{ii'}^2
-\frac1{16}\sum_{i<i'}(r_i-r_{i'})^2.
}
\tag{4.1}
\]

Moreover,

\[
\sum_{i<i'}h_{ii'}
=\frac14\left(a^2b-\sum_jc_j^2\right),
\qquad
\sum_{i<i'}(r_i-r_{i'})^2
=a\sum_i r_i^2.
\tag{4.2}
\]

#### Proof

For two rows, let \(u\) be the number of columns of pattern
\((+,-)\) and \(v\) the number of pattern \((-,+)\).  Their
alternating rectangles number \(uv\), while

\[
u+v=h_{ii'},\qquad
u-v=\frac{r_i-r_{i'}}2.
\]

This gives (4.1).  Summing columnwise gives the first identity in
(4.2); the second is the standard pairwise-square identity, using
\(\sum_i r_i=0\). \(\square\)

### Lemma 4.2 (dense alternating trades)

If

\[
L_B\le\frac18ab,
\tag{4.3}
\]

then

\[
\boxed{
T(B)\ge\frac{33}{2048}a^2b^2.
}
\tag{4.4}
\]

Consequently \(B\) contains at least

\[
\boxed{
\frac{33}{8192}ab
}
\tag{4.5}
\]

cell-disjoint alternating rectangles, up to an additive rounding
error.

#### Proof

Since \(|c_j|\le a\) and \(|r_i|\le b\),

\[
\sum_jc_j^2\le aL_B,\qquad
\sum_i r_i^2\le bL_B.
\]

Writing \(D=\sum_{i<i'}h_{ii'}\), (4.2) and (4.3) give

\[
D\ge\frac{7}{32}a^2b.
\]

Cauchy--Schwarz in (4.1) now yields

\[
\begin{aligned}
T(B)
&\ge
\frac{D^2}{4\binom a2}
-\frac a{16}\sum_i r_i^2\\
&\ge
\left(\frac{(7/8)^2}{32}-\frac1{128}\right)a^2b^2
=\frac{33}{2048}a^2b^2.
\end{aligned}
\]

A fixed cell belongs to at most
\((a-1)(b-1)\le ab\) rectangles.  Greedily choosing a rectangle and
discarding all rectangles meeting one of its four cells loses at
most \(4ab\) candidates per choice.  This proves (4.5). \(\square\)

### Theorem 4.3 (quadratic fixed-margin cube)

Let

\[
E_{\rm cross}
=\sum_{\alpha<\beta}k_\alpha k_\beta
=\frac12\left(n^2-\sum_\alpha k_\alpha^2\right).
\tag{4.6}
\]

The simultaneous row-and-column-margin fibre contains a Boolean cube
of dimension at least

\[
\boxed{
\frac{33}{8192}\bigl(E_{\rm cross}-8L\bigr)-O(q^2).
}
\tag{4.7}
\]

In particular, if
\(\max_\alpha k_\alpha\le(1-\delta)n\) and
\(M(A)\le Cn^{3/2}\), its dimension is

\[
\boxed{
\left(\frac{33\delta}{16384}-o_{C,\delta}(1)\right)n^2.
}
\tag{4.8}
\]

Hence the fibre has at least

\[
\exp\!\left[
\left(\frac{33\delta\log2}{16384}-o(1)\right)n^2
\right]
\tag{4.9}
\]

members.

#### Proof

Call a block bad if
\(L_{\alpha\beta}>k_\alpha k_\beta/8\).  The total number of cells
in bad blocks is at most \(8L\).  Apply Lemma 4.2 to every other
block and choose the alternating rectangles cell-disjoint inside
each block.  Rectangles in different blocks are already
cell-disjoint.  Flipping any subcollection independently preserves
all margins, and distinct subcollections give distinct signings.

Finally,
\[
E_{\rm cross}\ge\frac\delta2n^2
\]
under the hypothesis on the largest type, while (3.3) is
\(o(n^2)\).  The affine closure theorem also gives
\(q=O_C(n^{5/6})\), so the rounding loss \(O(q^2)\) is
\(o_C(n^2)\). \(\square\)

### 4.4 A genuine per-block obstruction

The global budget is essential.  Individual fixed-margin fibres can
be singletons.  If \(a<b\), \(b\) is even, and half of the columns
are constantly \(+1\) while half are constantly \(-1\), then all
row sums are zero and the prescribed column sums are \(\pm a\).
These margins determine the block uniquely.

They are compatible with row domination whenever the target type has
internal row sum at least \(a\), for example when \(b-1\ge a\).
What (3.2) proves is that such rigid or almost rigid blocks can cover
only \(O(n^{3/2})\) cross cells in a competitive nonrecursive
instance.

## 5. Exact fixed-spin laws on a rectangle cube

Let \(\mathcal R\) be a cell-disjoint family of alternating
rectangles.  On each rectangle choose independently between its two
orientations.  For spins \(u\in\{\pm1\}^a\),
\(v\in\{\pm1\}^b\), write the original orientation on a rectangle
\(\rho=(i,i';j,j')\) as

\[
a_\rho
\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.
\]

Then the random block energy is exactly

\[
\boxed{
u^\mathsf TBv
=m(u,v)+
\sum_{\rho\in\mathcal R}
\varepsilon_\rho
a_\rho(u_i-u_{i'})(v_j-v_{j'}),
}
\tag{5.1}
\]

where the \(\varepsilon_\rho\) are independent signs and \(m(u,v)\)
is the contribution of cells outside the randomized rectangles.

If \(N_{\rm act}(u,v)\) is the number of rectangles on which both
spin differences are nonzero, then

\[
\boxed{
\operatorname{Var}(u^\mathsf TBv)
=16N_{\rm act}(u,v),
}
\tag{5.2}
\]

\[
\boxed{
\mathbb E e^{t(u^\mathsf TBv-\mathbb Eu^\mathsf TBv)}
=\cosh(4t)^{N_{\rm act}(u,v)},
}
\tag{5.3}
\]

and

\[
\boxed{
\Pr\{|u^\mathsf TBv-\mathbb Eu^\mathsf TBv|\ge s\}
\le
2\exp\!\left(-\frac{s^2}{32N_{\rm act}(u,v)}\right).
}
\tag{5.4}
\]

Thus the quadratic fibre rank from Theorem 4.3 is real, not merely
enumerative.  It need not be isotropic: an adversarial spin can be
constant on many selected row or column pairs.  This is the first
place where a fusion-frame/inverse theorem is needed.

## 6. Uniform doubly balanced blocks: exact two-replica covariance

Assume \(a,b\) are even, and let \(B\) be uniform on

\[
\Omega^0_{a,b}
=\{B\in\{\pm1\}^{a\times b}:B\mathbf1=0,\
\mathbf1^\mathsf TB=0\}.
\tag{6.1}
\]

Row and column permutation symmetry gives, for two cells,

\[
\mathbb E B_{ij}B_{i\ell}=-\frac1{b-1},
\quad
\mathbb E B_{ij}B_{h j}=-\frac1{a-1},
\quad
\mathbb E B_{ij}B_{h\ell}
=\frac1{(a-1)(b-1)}
\tag{6.2}
\]

when the indicated indices are distinct.  Equivalently,

\[
\boxed{
\mathbb E[\operatorname{vec}B\,\operatorname{vec}B^\mathsf T]
=
\frac{ab}{(a-1)(b-1)}
(P_a\otimes P_b),
}
\tag{6.3}
\]

where \(P_a=I-a^{-1}J\).

For four spin vectors \(u,u'\in\{\pm1\}^a\) and
\(v,v'\in\{\pm1\}^b\), this gives the exact two-replica law

\[
\boxed{
\begin{aligned}
&\mathbb E
[(u^\mathsf TBv)(u'^\mathsf TBv')]\\
&\quad=
\frac{a\langle u,u'\rangle-s_us_{u'}}
{a-1}
\cdot
\frac{b\langle v,v'\rangle-s_vs_{v'}}
{b-1},
\end{aligned}
}
\tag{6.4}
\]

where \(s_u=\sum_i u_i\), and similarly for the other spins.
In particular,

\[
\boxed{
\operatorname{Var}(u^\mathsf TBv)
=
\frac{a^2-s_u^2}{a-1}
\frac{b^2-s_v^2}{b-1}.
}
\tag{6.5}
\]

The product in (6.5), rather than a sum, is intended:

\[
\operatorname{Var}(u^\mathsf TBv)
=
\left(\frac{a^2-s_u^2}{a-1}\right)
\left(\frac{b^2-s_v^2}{b-1}\right).
\]

Equation (6.4) is the exact gain from preserving both margins.  It
annihilates every fluctuation for which either endpoint spin is
constant on its type.  For many independent type blocks, put

\[
g_\alpha(x,y)
=
\frac{
k_\alpha\langle x_\alpha,y_\alpha\rangle
-s_\alpha(x)s_\alpha(y)
}{k_\alpha-1}.
\tag{6.6}
\]

Then the cross-energy covariance is

\[
\boxed{
K(x,y)
=\sum_{\alpha<\beta}g_\alpha(x,y)g_\beta(x,y)
=\frac12\left[
\left(\sum_\alpha g_\alpha\right)^2
-\sum_\alpha g_\alpha^2
\right].
}
\tag{6.7}
\]

This is the exact centered-overlap kernel.

## 7. A fully discrete tractable refill

There is a convenient margin-preserving ensemble with an exact
Rademacher representation.

For each even \(a\times b\) block:

1. choose a perfect matching of its \(a\) rows and one of its \(b\)
   columns;
2. tile the block by the resulting \(2\times2\) rectangles;
3. put an independently oriented checkerboard on every tile.

Every row and column sum is exactly zero.  If \(R_u\) is the number
of bichromatic row pairs for a spin \(u\), and \(R_v\) the analogous
column count, then, conditional on the matchings,

\[
\boxed{
u^\mathsf TBv
=4\sum_{\ell=1}^{R_uR_v}\varepsilon_\ell.
}
\tag{7.1}
\]

For a uniform perfect matching,

\[
\mathbb E R_u
=\frac{a^2-s_u^2}{4(a-1)}.
\tag{7.2}
\]

Consequently the averaged variance in (7.1) is exactly (6.5).
For two replicas, if

\[
\alpha_e(u)=u_i-u_{i'}
\]

on a matched row edge \(e=\{i,i'\}\), then

\[
\mathbb E_{\rm matching}
\sum_e\alpha_e(u)\alpha_e(u')
=
\frac{a\langle u,u'\rangle-s_us_{u'}}{a-1}.
\tag{7.3}
\]

Thus this discrete checkerboard ensemble has exactly the same
averaged two-replica covariance as the uniform doubly balanced
fibre.

## 8. The \(1/2\) two-replica barrier

Take \(q\) equal types of even size \(k\), with

\[
n=qk,\qquad q,k\to\infty,\qquad q\asymp k.
\]

Restrict to spins balanced inside every type.  There are

\[
\binom{k}{k/2}^q
=\exp(n\log2-o(n))
\tag{8.1}
\]

such states.

For the checkerboard ensemble, the exact Rademacher moment generating
function (7.1), at tilt \(t/\sqrt n\), has the Gaussian limit

\[
\frac1n\log
\mathbb E\exp\!\left(\frac t{\sqrt n}X_x\right)
\longrightarrow\frac{t^2}{4}.
\tag{8.2}
\]

The matching fluctuations contribute only \(o(n)\) to the logarithm.
This follows uniformly on every macroscopic overlap profile by the
switching bounded-difference inequality for a uniform perfect
matching: every required bichromatic or signed-overlap matching
count differs from its mean by \(O(\sqrt{k\log k})\) outside
\(\exp(-\Omega(\log k))\), while its contribution to a block
log-mgf is \(O(1/k)\) times that deviation.  Summed over
\(\Theta(q^2)=\Theta(n)\) blocks the error is \(o(n)\).
For two states whose within-type overlaps are
\(\rho_1,\ldots,\rho_q\), their limiting correlation is

\[
R(\rho)
=
\frac{(\sum_\alpha\rho_\alpha)^2-\sum_\alpha\rho_\alpha^2}
{q(q-1)}
=\bar\rho^2+o(1).
\tag{8.3}
\]

The joint upper-tail rate at level \(cn^{3/2}\) is therefore

\[
\boxed{
I_2(c,\rho)
=\frac{2c^2}{1+\bar\rho^2}+o(1).
}
\tag{8.4}
\]

The number of ordered pairs with a common overlap \(\rho\) has
entropy

\[
n\left[
\log2+h\!\left(\frac{1-\rho}{2}\right)
\right]+o(n).
\tag{8.5}
\]

After division by the square of the first moment, the two-replica
exponent is

\[
\boxed{
\Phi_c(\rho)
=
h\!\left(\frac{1-\rho}{2}\right)-\log2
+
\frac{2c^2\rho^2}{1+\rho^2}.
}
\tag{8.6}
\]

Pinsker's inequality gives

\[
\log2-h\!\left(\frac{1-\rho}{2}\right)
\ge\frac{\rho^2}{2}.
\tag{8.7}
\]

Hence

\[
\boxed{
\Phi_c(\rho)\le0
\quad\hbox{for every }\rho\in[-1,1]
\quad\Longleftrightarrow\quad c\le\frac12
}
\tag{8.8}
\]

at the local threshold.  More explicitly, near zero,

\[
\Phi_c(\rho)
=
\left(2c^2-\frac12\right)\rho^2+O(\rho^4).
\tag{8.9}
\]

For varying type overlaps, concavity of binary entropy and
\(R(\rho)\le\bar\rho^2+o(1)\) reduce the exponent to (8.6).
Thus the exponential-scale second moment has no clustering
obstruction below \(1/2\), and overlap zero becomes unstable
immediately above \(1/2\).

The same calculation for the centered Gaussian block process is
fully rigorous and gives

\[
\max_x X_x\ge(1/2-o(1))n^{3/2}
\]

with high probability: the second moment is
\(\exp(o(n))\) times the square of the first, and Gaussian
concentration upgrades the resulting
\(\exp(-o(n))\) probability.  The checkerboard process has the same
fixed-replica logarithmic moment generating functions; its uniform
second-moment transfer follows by perfect-matching concentration.

This is a statement about the generic isotropic refill ensemble, not
about every member of the fixed-margin fibre.

### 8.1 Why the first moment already fails

At the fully balanced microprofile, the cross variance is

\[
\left(\frac12+o(1)\right)n^2.
\]

The annealed threshold is therefore

\[
\sqrt{\log2}\,n^{3/2}
=0.83255\ldots\,n^{3/2},
\tag{8.10}
\]

well above \(n^{3/2}/2\).  The two-replica calculation shows more:
the natural refill itself already contains witnesses up to the
spectral half-scale before overlap clustering begins.  A first
moment, ordinary union bound, or unconditioned random refill cannot
produce the desired sub-\(1/2\) signing.

## 9. Spectral flatness plus affine closure

Assume now that every internal type block is all positive.  Then

\[
p(A)=\sum_\alpha\binom{k_\alpha}{2}.
\tag{9.1}
\]

Put \(N=n-1\), \(q=\#\{\alpha\}\), and
\[
u_\alpha=k_\alpha^{-1/2}\mathbf1_{V_\alpha}.
\]

Define the margin compensation

\[
\eta_\alpha
=
\sum_{\beta\ne\alpha}
\frac{\|B_{\beta\alpha}\mathbf1_{V_\alpha}\|_2^2}
{k_\alpha}.
\tag{9.2}
\]

The exact type-indicator identity is

\[
\boxed{
u_\alpha^\mathsf TA^2u_\alpha
=(k_\alpha-1)^2+\eta_\alpha.
}
\tag{9.3}
\]

Since every individual margin has magnitude at most its source type
size,

\[
\boxed{
\sum_\alpha\eta_\alpha
\le L\le2p(A).
}
\tag{9.4}
\]

### Theorem 9.1 (exact flat compensated floor)

If

\[
A^2=(n-1)I,
\tag{9.5}
\]

then

\[
\boxed{
2\sum_\alpha k_\alpha^2
\ge q(n-2)+3n.
}
\tag{9.6}
\]

Combining this with
\(\sum k_\alpha^2\ge n^2/q\) gives

\[
\boxed{
p(A)\ge
\left(\frac1{2\sqrt2}-o(1)\right)n^{3/2}.
}
\tag{9.7}
\]

#### Proof

Summing (9.3) under (9.5) gives

\[
\sum_\alpha\eta_\alpha
=q(n-1)-\sum_\alpha(k_\alpha-1)^2.
\]

Use (9.4),
\(2p=\sum k_\alpha^2-n\), and
\(\sum(k_\alpha-1)^2=\sum k_\alpha^2-2n+q\).
This is exactly (9.6).  If \(q=y\sqrt n\), the two lower bounds on
\(\sum k_\alpha^2\), followed by (9.1), give

\[
\frac{p(A)}{n^{3/2}}
\ge
\max\left\{\frac1{2y},\frac y4\right\}-o(1).
\]

The minimum occurs at \(y=\sqrt2\), proving (9.7). \(\square\)

### Corollary 9.2 (zero margins force the half ceiling)

If, in addition,

\[
B_{\alpha\beta}\mathbf1=0,\qquad
\mathbf1^\mathsf TB_{\alpha\beta}=0
\tag{9.8}
\]

for every cross block, then \(\eta_\alpha=0\), and (9.3)--(9.5)
force

\[
\boxed{
k_\alpha=1+\sqrt{n-1}
\quad\hbox{for every occupied type}.
}
\tag{9.9}
\]

Therefore

\[
\boxed{
p(A)=\frac12n\sqrt{n-1}
=\left(\frac12-o(1)\right)n^{3/2}.
}
\tag{9.10}
\]

Thus a conference-flat doubly centered affine residual cannot lie
strictly below the half ceiling: the type-constant Boolean
eigenvectors themselves attain it.

### 9.3 Approximate form

Assume \(q=O(\sqrt n)\) and

\[
\Delta:=\|A^2-(n-1)I\|_F^2=o(n^{5/2}).
\tag{9.11}
\]

Since the \(u_\alpha\) are orthonormal,

\[
\sum_\alpha
\left|
u_\alpha^\mathsf T(A^2-(n-1)I)u_\alpha
\right|
\le\sqrt q\,\sqrt\Delta=o(n^{3/2}).
\tag{9.12}
\]

The proof of Theorem 9.1 therefore gives

\[
p(A)\ge
\left(\frac1{2\sqrt2}-o(1)\right)n^{3/2}.
\tag{9.13}
\]

If also

\[
\sum_\alpha\eta_\alpha=o(n^{3/2}),
\tag{9.14}
\]

then

\[
\sum_\alpha
\left|(k_\alpha-1)^2-(n-1)\right|
=o(n^{3/2}).
\]

It follows that

\[
q=(1+o(1))\sqrt n,\qquad
p(A)=\left(\frac12-o(1)\right)n^{3/2}.
\tag{9.15}
\]

### 9.4 The precise projected fourth-moment alternative

The full Frobenius hypothesis in (9.11) can be replaced by exactly
the part seen by the type indicators.  Put

\[
d_\alpha
=u_\alpha^\mathsf T(A^2-(n-1)I)u_\alpha,
\qquad
\Delta_{\rm type}=\sum_\alpha d_\alpha^2.
\tag{9.16}
\]

Repeating the proof without discarding the error gives the exact
inequality

\[
\boxed{
2\sum_\alpha k_\alpha^2
\ge
q(n-2)+3n-\sqrt{q\Delta_{\rm type}}.
}
\tag{9.17}
\]

Consequently

\[
\boxed{
p(A)\ge
\frac12\left[
\max\left\{
\frac{n^2}{q},
\frac{q(n-2)+3n-\sqrt{q\Delta_{\rm type}}}{2}
\right\}
-n
\right].
}
\tag{9.18}
\]

For \(q=y\sqrt n+o(\sqrt n)\) and
\(\Delta_{\rm type}\le\delta n^{5/2}\), this becomes

\[
\boxed{
\frac{p(A)}{n^{3/2}}
\ge
\max\left\{
\frac1{2y},
\frac y4-\frac{\sqrt{y\delta}}4
\right\}
-o(1).
}
\tag{9.19}
\]

Thus there is a clean fourth-moment alternative:

* projected defect \(o(n^{5/2})\) gives the compensated
  \(1/(2\sqrt2)\) floor;
* with negligible margin compensation it gives the \(1/2\) floor;
* otherwise
  \(\Delta_{\rm type}=\Omega(n^{5/2})\).

The last branch cannot presently be converted into a Boolean witness.
The known Gaussian-rounding stability inequality controls

\[
\|A^2-(n-1)I\|_F^2
\]

only after division by \(n^{5/2}\) in an **additive** lower bound for
\(M(A)\).  A cap gap of order \(n^{3/2}\) therefore permits a defect
as large as order \(n^4\), whereas (9.19) needs control at order
\(n^{5/2}\).  This is a factor \(n^{3/2}\) gap in the defect and is
not a bookkeeping loss.  Fixed trace moments see typical row
correlations but do not see the exponentially sparse Boolean
witnesses which determine the ground energy.

Accordingly, this note does not infer spectral flatness merely from
\(M(A)<(1/2-\varepsilon)n^{3/2}\).  Proving that implication, or
replacing it by a direct high-degree/fusion-frame inverse theorem, is
the precise remaining obstruction.

## 10. What a sub-half flat exception must look like

Under exact flatness, write

\[
S_2=\sum_\alpha k_\alpha^2,\qquad
p(A)=c\,n^{3/2}+o(n^{3/2}),\qquad
q=y\sqrt n+o(\sqrt n).
\]

Then

\[
\sum_\alpha\eta_\alpha
=(y-2c)n^{3/2}+o(n^{3/2}),
\tag{10.1}
\]

while Cauchy gives \(y\ge1/(2c)+o(1)\).  Hence

\[
\boxed{
\frac{\sum_\alpha\eta_\alpha}{L}
\ge
\frac1{4c^2}-1-o(1).
}
\tag{10.2}
\]

Indeed \(L\le2p\).  More explicitly,

\[
\sum_\alpha\eta_\alpha
=
\sum_{\alpha\ne\beta}
\sum_{j\in V_\beta}
|r_{j,\alpha}|
\frac{|r_{j,\alpha}|}{k_\alpha}.
\tag{10.3}
\]

Thus (10.2) says that the \(L^1\)-margin-weighted mean of the
normalized bias
\(|r_{j,\alpha}|/k_\alpha\) is bounded away from zero.

Put

\[
\gamma(c)=\frac1{4c^2}-1.
\]

At least a fraction

\[
\boxed{
\frac{\gamma(c)}{2-\gamma(c)}
}
\tag{10.4}
\]

of the absolute margin mass lies on entries satisfying

\[
|r_{j,\alpha}|\ge\frac{\gamma(c)}2k_\alpha.
\tag{10.5}
\]

At the ROM value \(c=\sqrt{15}/8\),

\[
\gamma(c)=\frac1{15},
\]

so at least \(1/29-o(1)\) of the absolute margin mass lies on
columns with bias at least \(k_\alpha/30\).  At the compensated
floor \(c=1/(2\sqrt2)\), \(\gamma(c)=1\): essentially all margin
mass must be carried by nearly constant columns.

This is the promised exceptional-set dichotomy:

* negligible or isotropic margin channels force the type sizes to
  \(\sqrt n\) and the cap to \(1/2\);
* a genuinely sub-half flat exception must use a macroscopic
  \(n^{3/2}\)-scale family of biased columns;
* near the \(1/(2\sqrt2)\) endpoint those columns are almost
  constant and are natural candidates for a Ferrers/paired quotient
  and principal descent.

The last implication has not yet been proved for moderate bias such
as \(k_\alpha/30\).

## 11. Deterministic refill audit

The obvious exact tight-frame constructions do not close the gap.

1. A square zero-row/zero-column sign block cannot be an exact tight
   frame on \(\mathbf1^\perp\) for \(k>2\).  Such a frame would have
   Gram matrix
   \[
   \frac{k^2}{k-1}P_{\mathbf1^\perp},
   \]
   whose off-diagonal entries are \(-k/(k-1)\), incompatible with
   integer row inner products.
2. Removing the first row from a normalized Hadamard matrix of
   order \(k\) gives a \((k-1)\times k\) sign block with orthogonal
   rows and total zero.  Its column margins are
   \(k-1,-1,\ldots,-1\).  The one extreme column consumes an entire
   target vertex's row-domination budget, while all remaining
   columns still cost one unit.  Such blocks cannot populate a dense
   \(\Theta(\sqrt n)\)-type quotient.
3. Tensoring a Hadamard core with a \(2\times2\) checkerboard cancels
   the margins, but introduces repeated/opposite rows and exactly the
   paired quotient/descent already identified in the principal
   restriction audit.

Random dense regular bipartite sign blocks are approximately flat,
but Section 8 shows that their isotropic Boolean process already
reaches the half scale.  A useful deterministic construction would
therefore have to be simultaneously:

* margin preserving;
* spectrally/fusion-frame flat;
* arithmetically nonresonant on Boolean vectors; and
* free of repeated/opposite fibre quotients.

No known Hadamard or regular-bipartite gadget satisfies all four.

## 12. Surviving inverse target

The results above reduce the affine refill branch to the following
specific theorem target.

> Let a competitive exact affine-ground signing have
> \(\Omega(n^2)\) switchable fixed-margin directions.  If no
> member of its fixed-margin fibre has a Boolean witness at the
> half scale under an isotropic refill, then either:
> 
> 1. the rectangle-trade covariance has an
>    \(\Omega(n)\)-dimensional anisotropic kernel which yields an
>    affine/Ferrers quotient; or
> 2. the margin compensation in (9.2) has the biased-column mass
>    from (10.4), and that biased incidence structure yields a
>    smaller principal signing with no worse normalized objective.

Theorem 4.3 supplies the quadratic rank.  Section 8 identifies the
precise overlap obstruction that must fail.  Sections 9--10 show
what spectral flatness would then force.  The missing step is the
inverse passage from nonisotropic rectangle features or biased
columns to a bounded-depth principal descent.
