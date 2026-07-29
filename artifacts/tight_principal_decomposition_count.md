# Tight principal decompositions and ground-layer counting

## 1. Setup

For a symmetric zero-diagonal signing \(A\) of order \(n\), use the
one-copy energy

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
P(A)=\max_x H_A(x),
\qquad
Q(A)=\max_x|H_A(x)|.
\]

Let

\[
\mathcal G^+(A)=
\{x\in\{\pm1\}^n/\{\pm\mathbf1\}:H_A(x)=P(A)\}.
\]

For a nonempty proper vertex set \(S\), put

\[
\operatorname{Tight}_A(S)
\iff
P(A_S)+P(A_{S^c})=P(A).
\]

The aim of this note is to count such tight bipartitions under the
global assumption \(Q(A)=O(n^{3/2})\).

## 2. Exact difference--decomposition equivalence

Represent projective spins by fixing one vertex and identify them with
\(\mathbb F_2^{n-1}\).  Let \(Z\) be the resulting positive ground
family and

\[
r_Z(d)=|Z\cap(Z+d)|.
\]

Then, for every nonzero \(d\), whose two shores are \(S,S^c\),

\[
\boxed{
r_Z(d)>0
\quad\Longleftrightarrow\quad
\operatorname{Tight}_A(S).
}
\tag{2.1}
\]

To prove the forward implication, switch by one representation of
\(d=z+(z+d)\), and write

\[
A=
\begin{pmatrix}D&B\\B^\top&E\end{pmatrix}.
\]

The two ground energies are

\[
H_D(u)+H_E(v)\pm u^\top Bv=P(A),
\]

so the cross term is zero and the internal sum is \(P(A)\).  Hence
\(P(D)+P(E)\ge P(A)\), while choosing independent block maximizers and
their better relative sign gives

\[
P(D)+P(E)+|s^\top Bt|\le P(A).
\]

Thus equality holds.  Conversely, if \(P(D)+P(E)=P(A)\), then every
pair of positive block maximizers has zero cross term by the same
inequality, and their two relative signs give two full positive ground
states differing on \(S\).

More precisely, if \(\gamma_S,\gamma_{S^c}\) are the projective
positive-ground degeneracies of the two blocks, then

\[
\boxed{r_Z(d)=2\gamma_S\gamma_{S^c}.}
\tag{2.2}
\]

Therefore

\[
\boxed{
|\{S:\operatorname{Tight}_A(S)\}|
=|\operatorname{supp}(Z+Z)|-1,
}
\tag{2.3}
\]

with complements identified if desired.  A tight-decomposition count
is exactly a difference-set count for the positive ground face.

The additive energy consequently satisfies

\[
E(Z)=\sum_d r_Z(d)^2\ge3|Z|^2-2|Z|,
\tag{2.4}
\]

with equality exactly when \(Z\) is Sidon, equivalently when both
blocks of every realized tight decomposition have unique projective
positive grounds.

## 3. Determinant and rank information

For a tight partition, let \(U_S\) and \(U_{S^c}\) be the real spans
of the positive block-ground vectors.  Cartesian annihilation says

\[
\boxed{U_S^\top A_{S,S^c}U_{S^c}=0.}
\tag{3.1}
\]

Hence

\[
\operatorname{rank}(A_{S,S^c})
\le |S|+|S^c|-\dim U_S-\dim U_{S^c}.
\tag{3.2}
\]

This yields a determinant/minor obstruction only when both block
ground layers have large real span.  It gives no information in the
Sidon regime, where both spans may be one-dimensional.  In particular,
the rank route and the additive-energy route are complementary rather
than cumulative: high difference multiplicity forces rank deficiency,
while minimum difference multiplicity evades it completely.

There is a parity obstruction that is sometimes useful.  Since
\(s^\top A_{S,S^c}t\) is a sum of \(|S||S^c|\) signs, a tight
partition requires

\[
\boxed{|S||S^c|\ \text{even}.}
\tag{3.3}
\]

For even \(n\), every difference between two same-orientation ground
states therefore has even Hamming weight.

## 4. A scalable counterfamily at the \(n^{3/2}\) scale

A general bound \(\exp(o(\sqrt n))\) on either positive ground states
or tight decompositions is false, even under the sharp spectral-scale
condition \(Q(A)=O(n^{3/2})\).

Let \(r\ge1\), put

\[
N=2^{2r},\qquad b=2^r=\sqrt N,
\]

and index rows and columns by pairs
\((u,v)\in\mathbb F_2^r\times\mathbb F_2^r\).  Define the symmetric
Hadamard matrix

\[
K_{(u,v),(x,y)}
=(-1)^{v\cdot x+u\cdot y}.
\tag{4.1}
\]

Then

\[
K=K^\top,\qquad K^2=NI,\qquad K_{zz}=1.
\]

Consequently

\[
A=K-I
\tag{4.2}
\]

is a symmetric zero-diagonal \(\{\pm1\}\)-matrix, and

\[
\boxed{Q(A)\le \frac12N(\sqrt N+1)=O(N^{3/2}).}
\tag{4.3}
\]

For every Boolean function

\[
g:\mathbb F_2^r\to\mathbb F_2,
\]

define the Boolean vector

\[
X_g(x,y)=(-1)^{x\cdot y+g(x)}.
\tag{4.4}
\]

A direct Walsh summation gives

\[
\begin{aligned}
(KX_g)(u,v)
&=\sum_{x,y}
(-1)^{v\cdot x+u\cdot y+x\cdot y+g(x)}\\
&=2^r(-1)^{u\cdot v+g(u)}
=\sqrt N\,X_g(u,v).
\end{aligned}
\tag{4.5}
\]

Therefore every \(X_g\) is a positive Boolean top eigenvector.  Since
\(x^\top Kx\le N\sqrt N\) for every Boolean \(x\), and the deleted
diagonal contributes the constant \(N\),

\[
\boxed{
P(A)=\frac12(N\sqrt N-N)
}
\tag{4.6}
\]

in the one-copy normalization, and

\[
\boxed{
|\mathcal G^+(A)|\ge2^{\sqrt N-1}.
}
\tag{4.7}
\]

The products \(X_gX_h=(-1)^{g(x)+h(x)}\) form a projective binary
subspace of dimension \(\sqrt N-1\).  Hence this one matrix has at
least

\[
\boxed{2^{\sqrt N-1}-1}
\tag{4.8}
\]

nontrivial tight principal bipartitions.  The vertices fall into
\(\sqrt N\) evaluation types of size \(\sqrt N\), exactly realizing
the frozen-type structure forced by Fourier analysis when an entire
zero-cut subspace is present.

This construction is a coordinate-explicit version of the
Bush-type-Hadamard phenomenon: a Bush-type Hadamard matrix of order
\(4u^2\) has at least \(2^{2u}\) Boolean top eigenvectors.  Thus the
\(\exp(\Theta(\sqrt n))\) barrier is structural, not a small-order
artifact.

## 5. Consequences for the counting strategy

The following possible statements are now separated.

1. **False:** uniformly \(O(1)\), polynomially many, or
   \(\exp(o(\sqrt n))\) tight decompositions whenever
   \(Q(A)=O(n^{3/2})\).
2. **Still possible:** \(\exp(O(\sqrt n\,\mathrm{polylog}\,n))\)
   decompositions under only \(Q(A)=O(n^{3/2})\).
3. **What scale transfer actually needs:** after choosing a deletion
   scale \(h\gg\sqrt n\,\mathrm{polylog}\,n\), the forced exact layer
   has \(\exp(\Omega(h))\) members.  Therefore even a uniform
   \(\exp(O(\sqrt n\,\mathrm{polylog}\,n))\) bound would still close
   that version of the repair argument.

The explicit family also shows exactly why a bare Littlewood--Offord
argument cannot work.  For its tight partitions, the cross block has
a common Boolean nullspace generated by the repeated evaluation
types.  Inverse Littlewood--Offord can only rediscover this low-rank
block structure; it does not contradict the \(n^{3/2}\) norm bound.

There is nevertheless an exact rigidity theorem at the opposite,
small-block endpoint.  It helps explain why the construction above
uses blocks of order \(\sqrt n\), rather than bounded blocks.

### Proposition 5.1 (paired-coordinate cube forces quadratic norm)

Partition \(2d\) vertices into pairs \(V_1,\ldots,V_d\).  Suppose all
\(2^d\) configurations that are constant on each pair are positive
ground states of a signing \(A\).  Then

\[
\boxed{Q(A)\ge 2d^2-d=\frac{n^2-n}{2}.}
\tag{5.1}
\]

Here \(Q\) is in the one-copy normalization.

To prove this, write \(B_{ij}\) for the \(2\times2\) block between
two pairs.  Constancy of the ground energy under independent pair
signs first gives

\[
\mathbf1^\top B_{ij}\mathbf1=0.
\tag{5.2}
\]

Changing one pair from its uniform mode \(\mathbf1\) to its
antiuniform mode \(v=(1,-1)\), while choosing all other pair signs
adversarially, shows that the internal edge of every pair is \(+1\)
and

\[
\sum_{j\ne i}|v^\top B_{ij}\mathbf1|\le2.
\tag{5.3}
\]

Under (5.2), every summand in (5.3) is either \(0\) or \(4\), hence
all are zero.  Applying the same argument at the other endpoint gives
zero column sums as well.  Therefore

\[
B_{ij}=c_{ij}vv^\top,\qquad c_{ij}\in\{\pm1\}.
\tag{5.4}
\]

Change any three pairs to their antiuniform modes.  Relative to the
ground energy, the change is

\[
-6+4(c_{ij}t_it_j+c_{ik}t_it_k+c_{jk}t_jt_k).
\tag{5.5}
\]

It must be nonpositive for every \(t_i,t_j,t_k\).  A signed triangle
has maximum one-copy energy \(3\) when its sign product is \(+1\) and
\(1\) when its sign product is \(-1\).  Thus

\[
c_{ij}c_{ik}c_{jk}=-1
\]

for every triangle.  Hence \(c_{ij}=-\sigma_i\sigma_j\) for suitable
signs \(\sigma_i\).  Putting every pair in antiuniform mode with
\(t_i=\sigma_i\) yields energy

\[
d-2d-4\binom d2=-2d^2+d,
\]

which proves (5.1).

Thus a linear-dimensional product cube of bounded two-vertex
directions is incompatible with the \(O(n^{3/2})\) regime.  The
natural generalization suggested by (5.1) and the
\(\exp(\Theta(\sqrt n))\) Hadamard example is:

\[
\dim W=O(\sqrt n\,\mathrm{polylog}\,n)
\]

for every *coordinate-block* ground subspace \(W\) of a signing with
\(Q(A)=O(n^{3/2})\).  Extending the triangle argument from
\(2\times2\) balanced blocks to larger blocks is the most concrete
remaining determinant/inverse-Littlewood--Offord target.

The viable inverse theorem would need to say:

> If a signing with \(Q(A)=O(n^{3/2})\) has more than
> \(\exp(C\sqrt n\,\mathrm{polylog}\,n)\) tight decompositions, then
> the corresponding ground differences contain a large additive
> subspace, and the frozen-type quotient generated by that subspace
> creates a Boolean quadratic witness larger than \(O(n^{3/2})\).

The first conclusion is not supplied by ordinary
Balog--Szemerédi--Gowers because exact minimizers through order seven
already have Sidon positive-ground layers.  Proving it would require a
new inverse theorem using the fact that every realized difference is a
tight *principal* decomposition, not merely an additive difference.
