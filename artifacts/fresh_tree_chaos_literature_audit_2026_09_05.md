# Elementary odd-tree chaos CLT: independent audit and joint-root extension

Date: 2026-09-05. Auditor: fresh-limit literature subagent.

## Verdict

**Passed.** The proposed argument gives a direct Rademacher moment proof for
every fixed tree with external root of degree one and all free vertices of
odd degree. No Gaussian invariance or Wiener-chaos limit theorem is needed.
The same argument gives joint moments at finitely many external roots, with
the expected row-Gram covariance.

This is a theorem about the explicitly defined injective tree polynomials.
It does not by itself justify removing injectivity, replacing nonlinear
responses by tree expansions, controlling their tails, or taking tree size
or response depth to infinity with dimension. Those remain separate proof
obligations in any application.

## 1. Exact statement and normalization

Let `A=A_n` be hollow symmetric with off-diagonal entries in `{+1,-1}`, and
write

\[
\beta(A)=\max_{x,y\in\{\pm1\}^n}|x^\top Ay|.
\]

The application assumes `beta(A)=O(n^(3/2))`. The proof actually only needs
`beta(A)/n²->0`, with the explicit moment errors below.

Let `T` be a fixed finite tree, with distinguished external root `o` of
degree one. Its `d` remaining vertices are called free vertices. Assume
every free vertex has odd degree. For `i in [n]`, define

\[
X_i(T)=n^{-d/2}\sum_{\substack{\varphi:V(T)\hookrightarrow[n]\\
                               \varphi(o)=i}}
 \left(\prod_{uv\in E(T)}A_{\varphi(u),\varphi(v)}\right)
 \left(\prod_{v\ne o}S_{\varphi(v)}\right),
\]

where `S` consists of independent Rademacher signs. The sum is over
injective embeddings; abstract vertices of the tree are individually
distinguished in this sum. Let `a(T)` be the number of automorphisms of `T`
fixing its external root.

Then uniformly over permitted matrices and selected roots,

\[
X_i(T)\Rightarrow N(0,a(T)).
\]

For a finite collection of nonisomorphic externally rooted trees, the
same-root limits are independent. Rooted-isomorphic tree polynomials are
in fact identical after renaming their abstract vertices, so such copies
should not be treated as independent coordinates.

The handshake lemma shows that `d` is odd: all `d+1` vertices have odd
degree. Thus every such polynomial is an odd homogeneous multilinear
polynomial in the probe signs, and its odd moments vanish exactly.

## 2. Moment partitions and the negligible terms

First consider a product of `l` copies at the same external root, allowing
different fixed trees satisfying the hypotheses. Let `D` be the total number
of their free-vertex occurrences; there are also exactly `D` edge
occurrences, since a tree with `d` free vertices has `d` edges.

Expand the product and partition the free occurrences according to their
image labels. Within one copy, no two occurrences can belong to the same
class, by injectivity. A nonzero spin expectation requires every class to
have even size. If `D` is odd, every expectation is zero. Otherwise there
are at most `V=D/2` distinct free labels.

For every partition with fewer than `V` classes, the number of embeddings
is `O(n^(V-1))`. All edge products have absolute value at most one, whereas
the normalization is `n^-V`. Thus the total contribution of these terms is
`O(1/n)`. There are only finitely many partitions, with their number
depending on the fixed moment and trees, not on the matrix or external root.

It remains to consider partitions in which every class has size exactly
two. In particular, each free label belongs to exactly two different
copies. This last fact is crucial in the propagation argument below.

## 3. Cancellation when an odd-multiplicity edge remains

For a leading partition form the quotient multigraph of all tree copies,
with the external roots identified. Reduce every edge multiplicity modulo
two, obtaining a simple graph `H`. There are no loops in an admissible
partition: the endpoints of every original edge belong to one injectively
embedded copy.

Every free quotient vertex has even degree in `H`, because it joins two
odd-degree free occurrences. The external root also has even degree: the
number of copies is even whenever `D` is even, since each tree has odd `d`.
Thus `H` is Eulerian. If it is nonempty, it must contain an edge between
two free vertices. Indeed a graph containing only external-root-to-free
edges would give degree one to every incident free vertex.

Choose such an edge `uv` and fix all other labels. Every remaining factor
depending on the label of `u` or `v` is a bounded unary weight: the only
factor jointly depending on the two labels is `Auv`. Therefore the sum over
these two labels is bounded in absolute value by `beta(A)`. Forbidden
previous labels can be imposed by setting the corresponding unary weights
to zero. The constraint `u!=v` is automatic for this sum because `Auu=0`.

There are at most `n^(V-2)` choices of the remaining labels. Consequently
the normalized contribution of this partition is at most

\[
\frac{n^{V-2}\beta(A)}{n^V}=\frac{\beta(A)}{n^2}.
\]

This estimate is uniform in the external root. It controls the entire
signed diagram and does not require estimating the signs of its individual
embeddings. In the low-cap application it is `O(n^-1/2)`.

## 4. Empty mod-two graph forces complete-copy pairing

Suppose `H` is empty. The full quotient graph is connected, has `V+1`
vertices including its external root, and has `2V` edge occurrences.
Every edge multiplicity is even. If `E` is the number of distinct edges,
connectivity gives `E>=V`, while even multiplicities give `E<=V`.
Therefore `E=V`, the underlying graph is a tree, and every distinct edge
has multiplicity exactly two.

The first edge of any copy must be paired with the first edge of another
copy, because each copy has only one edge incident to the external root.
Let these copies be `P` and `R`. Their first-child label occurs exactly
twice, hence it occurs only in these two copies. Every edge incident to that
label appears at most once in each copy and must appear twice in total.
Its incident edges therefore match between `P` and `R`.

Every adjacent free label now occurs in both `P` and `R`; its two occurrences
are exhausted too. Repeating the same argument propagates through the
connected image of `P`, forcing a root-preserving graph isomorphism between
the whole of `P` and the whole of `R`. It cannot transfer to a third copy,
because every encountered free label already has its two occurrences.

Different pairs have disjoint free-label sets. Thus the surviving partitions
are exactly pairings of whole rooted-isomorphic copies, with one rooted
isomorphism per pair. There is no additional family of partial pairings.

For `l` copies of a single tree, their number is

\[
(l-1)!!\,a(T)^{l/2}.
\]

Each contributes `(n-1)_V/n^V=1+O(1/n)`. These are precisely the centered
Gaussian Wick moments of variance `a(T)`. For mixed trees, only pairings
between rooted-isomorphic copies survive.

## 5. Joint external roots and the row-Gram covariance

The same counting works for finitely many copies placed at finitely many
possibly distinct fixed external roots. First discard terms in which a free
occurrence uses the label of any fixed external root. A nonzero spin
expectation requires at least two free occurrences of any such fixed label,
so at most `D/2-1` genuinely variable labels remain. Their total normalized
contribution is `O(1/n)`. Thus in leading terms all free labels may be assumed
distinct from every fixed external-root label.

As before, nonleading partitions cost `O(1/n)`. For a leading pair
partition, if the mod-two graph contains a free-to-free edge, the same
bilinear argument bounds its contribution by `beta(A)/n²`.

Consider the remaining case, where all free-to-free edges have even
multiplicity. A free label belongs to exactly two copies and can be adjacent
to the external root at most once in each copy. Thus it has at most two
external-root incidences, counted with multiplicity. Its total degree is
even, and all its free-edge multiplicities are even, so the number of its
external-root incidences is either zero or two.

In particular, a first-child occurrence is paired with another first-child
occurrence. The two corresponding copies may have different external roots.
Every free-to-free edge has multiplicity at most two, because either
endpoint occurs in only these two copies; its positive even multiplicity
is therefore exactly two. The same propagation argument now pairs the
entire free trees, preserving their first-child attachment. Hence the whole
copies are again externally rooted-isomorphic.

For a paired copy rooted at `i` and a copy rooted at `j`, all internal edge
weights cancel, while the first-child label `u` leaves the weight
`Aiu Aju`. After temporarily removing the distinctness and fixed-root
exclusions, its normalized sum is

\[
\frac1n\sum_u A_{iu}A_{ju}=\frac{(A^2)_{ij}}n.
\]

Restoring those exclusions changes each fixed diagram by `O(1/n)`, since
all weights are bounded and only `O(n^(V-1))` assignments are excluded.
The remaining free variables have no surviving weights. Therefore the
full joint moment expansion is the Gaussian Wick sum with covariance

\[
\boxed{\quad
\operatorname{Cov}_{\rm limiting}(X_i(T),X_j(U))
=\begin{cases}
 a(T)\displaystyle\lim\frac{(A^2)_{ij}}n,
       &T\cong U\text{ as externally rooted trees},\\
 0,&\text{otherwise},
\end{cases}\quad}
\]

whenever the finitely many Gram entries have limits. More precisely, the
joint moments differ from the Wick sums formed with the **finite-order**
entries `(A²)ij/n` by
`O(1/n+beta(A)/n²)`, with constants depending only on the fixed trees and
moment orders.

For `B=A/sqrt(n-1)` this covariance is asymptotically `a(T)(B²)ij`.
Consequently the familiar condition `(B²)ij->0` gives independence of
different external-root blocks. No independent Gaussian-chaos theorem is
needed to justify this two-coordinate extension.

## 6. Moment convergence really gives the claimed distributional limit

The estimates hold for every fixed moment, uniformly over the matrices and
roots. Second moments give tightness; bounded moments of higher order give
uniform integrability needed to pass any lower moment to subsequential
limits. The limiting moments are Gaussian.

Gaussian moment determinacy can be seen directly here: if a random variable
has the centered Gaussian moments of variance `sigma²`, Tonelli applied to
the nonnegative power series for `cosh(tX)` gives
`E cosh(tX)=exp(sigma²t²/2)` for every real `t`. This supplies exponential
integrability, and expansion of its moment-generating function identifies
the Gaussian law. Apply this argument to every fixed linear combination of
the finitely many tree variables; their mixed moment limits give the
Gaussian moments with the corresponding covariance quadratic form.

Thus the joint distributional conclusion follows, including singular
limiting covariance matrices. Uniformity follows by applying the same
argument to any contrary sequence of matrices and selected roots.

## 7. Hypothesis checks and scoped counterexamples

The external degree-one condition is essential. If the external root has
two leaf children, all free vertices still have odd degree, but the
normalized polynomial is a centered square of its linear row field and
converges to `Z²-1`, not a Gaussian. Partial pairings can then exchange
branches at the external root, so whole-copy propagation fails there.

Odd free degrees are also doing actual work. A length-two path has one
even-degree free vertex. At the first row of a hollow Sylvester matrix,
its quadratic rooted field has variance tending to two, although the
rooted path has only one automorphism. This is the row-sum obstruction
identified in the general-mask audit. It corresponds precisely to a
nonempty mod-two diagram supported only on external-root edges, for which
there is no free-to-free edge to expose to the bilinear cap.

The original low-cap assumption supplies the needed discrepancy bound by
ordinary polarization, `beta(A)<=4q(A)=O(n^(3/2))`. No spectral regularity,
conference structure, gauge, or hidden Gaussian replacement is required
for the odd-tree theorem itself.
