# Speed-\(n^2\) Laplace principle: the necessary state space and a no-go for standard limits

## 1. Setup and outcome

Write
\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
M(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,\qquad
G_n(A)=\frac{M(A)}{n^{3/2}},
\]
and let \(U_n\) denote the uniform law on the \(2^{N_n}\) edge
signings, where \(N_n=\binom n2\).  The canonical pressure is
\[
\Phi_n(\beta)
=\frac{N_n\log2}{n^2}
 +\frac1{n^2}\log\mathbb E_{U_n}
  \exp\{-\beta n^2G_n(A)\}.
\tag{1.1}
\]

This investigation did **not** prove convergence of (1.1).  It did
produce:

1. an explicit compact projective state space which retains all
   hypercube rank-one projectors and on which \(G_n\) is continuous at
   every bounded-energy limit;
2. automatic exponential tightness at speed \(n^2\);
3. a proof that a full speed-\(n^2\) LDP on this state space would
   already force convergence of \(F(n)/n^{3/2}\);
4. a deterministic pair of signing sequences which are
   asymptotically identical to graphons, all fixed-replica laws, and
   normalized spectral noncommutative distributions, but whose
   \(G_n\)'s differ by a fixed positive constant;
5. a support-set coupling under one-vertex gauge deletion, retaining
   all replica levels simultaneously.

Thus ordinary graphon, spectral-microstate, and fixed-order traffic
LDPs cannot be contracted to the desired observable.  The minimal
state that can see the observable contains its extremal support
already at level one; proving uniqueness of its LDP is essentially the
open problem in a correct topological form.

## 2. Exact convex-body, counting, and entropy formulas

For \(x\in\{\pm1\}^n\), put
\[
v_x=(x_ix_j)_{i<j}\in\{\pm1\}^{N_n},\qquad
K_n=\operatorname{conv}\{\pm v_x:x\in\{\pm1\}^n\}.
\]
Then
\[
M(A)=h_{K_n}(A),
\tag{2.1}
\]
so the sublevel set is the cube-lattice section
\[
\{A:G_n(A)\le c\}
=\{\pm1\}^{N_n}\cap c n^{3/2}K_n^\circ.
\tag{2.2}
\]
If
\[
Z_n(c)=\#\{A:G_n(A)\le c\},
\]
then layer-cake integration gives the exact identity
\[
\boxed{
\sum_Ae^{-\beta n^2G_n(A)}
=\beta n^2\int_0^\infty
 e^{-\beta n^2c}Z_n(c)\,dc.}
\tag{2.3}
\]
Thus a convex-body support-function LDP here is precisely an LDP for
the cube-lattice counts in (2.2); replacing lattice counts by
Euclidean volume is not justified.

The Gibbs variational formula is
\[
\boxed{
\Phi_n(\beta)
=\sup_{\nu\in\mathcal P(\{\pm1\}^{N_n})}
\left\{
\frac{H(\nu)}{n^2}
-\beta\mathbb E_\nu G_n(A)
\right\}.}
\tag{2.4}
\]
Equivalently, for the uniform-expectation term in (1.1),
\[
\frac1{n^2}\log\mathbb E_{U_n}e^{-\beta n^2G_n}
=\sup_\nu
\left\{
-\frac{D(\nu\Vert U_n)}{n^2}
-\beta\mathbb E_\nu G_n
\right\}.
\tag{2.5}
\]
Since \(G_n\) is invariant under vertex permutations, vertex
switchings, and global negation, averaging \(\nu\) over this group
preserves the energy term and increases entropy.  Hence the supremum
in (2.4) may be restricted exactly to invariant laws on Seidel
switching classes.  This removes labels but does not remove the
second-order extremal landscape.

## 3. A compact projective state retaining every Boolean projector

Let
\[
P_x=\frac{xx^\top}{n},\qquad
e_A(x)=\frac{H_A(x)}{n^{3/2}}
=\frac12\operatorname{Tr}\left(\frac A{\sqrt n}P_x\right).
\]
For every \(k\ge1\), define the finite support set
\[
\mathcal T_{n,k}(A)
=
\left\{
\left(
  \bigl(\operatorname{Tr}(P_{x^a}P_{x^b})\bigr)_{a<b},
  \bigl(\tanh e_A(x^a)\bigr)_{a\le k}
\right):
x^1,\ldots,x^k\in\{\pm1\}^n
\right\}.
\tag{3.1}
\]
Here
\[
\operatorname{Tr}(P_xP_y)
=\left(\frac{x\cdot y}{n}\right)^2,
\]
so (3.1) records the finite Gram geometry of the actual rank-one
hypercube projectors, not merely sampled energy moments.

For
\[
X_k=[0,1]^{\binom k2}\times[-1,1]^k,
\]
let \(\mathcal K(X_k)\) be its hyperspace of nonempty compact subsets
with the Hausdorff metric.  Each \(\mathcal K(X_k)\) is compact.  The
family
\[
\mathcal T_n(A)=(\mathcal T_{n,k}(A))_{k\ge1}
\]
satisfies all deletion, repetition, and permutation consistencies.
Those consistencies define a closed subset
\[
\mathfrak T\subseteq\prod_{k\ge1}\mathcal K(X_k).
\tag{3.2}
\]
Consequently \(\mathfrak T\) is compact and metrizable.

At level one, put
\[
\rho(T)=\max\{|y|:y\in T_1\}.
\]
Then
\[
\boxed{G_n(A)=\operatorname{arctanh}\rho(\mathcal T_n(A)).}
\tag{3.3}
\]
The right side is understood as \(+\infty\) at \(\rho=1\).  In
particular, \(G\) is continuous at every state with bounded energy.

Let \(\mu_n\) be the pushforward of \(U_n\) to \(\mathfrak T\).
Because \(\mathfrak T\) is compact, \((\mu_n)\) is exponentially tight
at **every** speed, in particular at speed \(n^2\).

There is a useful general subsequence fact.  From any subsequence one
can extract a further subsequence on which
\[
\frac1{n^2}\log\int e^{n^2f(T)}\,d\mu_n(T)
\]
converges for every \(f\) in a countable uniformly dense subset of
\(C(\mathfrak T)\).  The \(1\)-Lipschitz dependence on \(f\) extends
this to all continuous \(f\), and the inverse Varadhan lemma yields a
good subsequential LDP.  Therefore the obstruction is not existence
of subsequential rate functions; it is uniqueness of the rate
function (or merely uniqueness of its variational value for
\(-\beta G\)).

### 3.1 A full LDP would already solve the original problem

Suppose the full sequence \((\mu_n)\) obeys a speed-\(n^2\) LDP on
\(\mathfrak T\), with good rate \(I\).  Then \(F(n)/n^{3/2}\)
converges.

To prove this, suppose instead that its liminf is \(a\) and its limsup
is \(b>a\).  Choose minimizers \(A_{n_j}\) along a subsequence tending
to \(a\).  Compactness gives, after extraction,
\[
\mathcal T_{n_j}(A_{n_j})\longrightarrow T,\qquad G(T)=a.
\]
Every neighborhood of \(T\) has \(\mu_{n_j}\)-mass at least one cube
atom:
\[
2^{-N_{n_j}}
=\exp\left\{-\left(\frac{\log2}{2}+o(1)\right)n_j^2\right\}.
\]
The LDP upper bound on shrinking closed neighborhoods therefore gives
\[
I(T)\le\frac{\log2}{2}<\infty.
\tag{3.4}
\]
Choose \(c\in(a,b)\) and an open neighborhood \(O\ni T\) on which
\(G<c\).  The LDP lower bound gives
\[
\liminf_n\frac1{n^2}\log\mu_n(O)
\ge-\inf_O I>-\infty.
\]
Thus \(O\) is nonempty at every sufficiently large order, forcing
\(F(n)/n^{3/2}<c\), contrary to the limsup \(b\).

This proof also isolates the sparse-phase issue: a single signing, or
a single switching orbit, has finite maximal rate
\(\log2/2\), not rate \(+\infty\).  Such a phase cannot be discarded
from a speed-\(n^2\) LDP when \(\beta\) is large.

If the full LDP were available, truncating \(G\) and applying Varadhan
would also give
\[
\lim_n\frac1{n^2}\log
\mathbb E e^{-\beta n^2G_n}
=-\inf_{T\in\mathfrak T}\{I(T)+\beta G(T)\}.
\tag{3.5}
\]
The truncation is harmless: deterministic signings with
\(G_n=O(1)\) give an exponential lower bound, while the contribution
of \(G>L\) is at most \(e^{-\beta L n^2}\).

The price of making the topology correct is visible in (3.3):
the one-replica **support**, rather than its sampled law, already
contains the original maximum.  Hence a full LDP on \(\mathfrak T\)
is a correct reformulation, but not a black-box shortcut.

## 4. A fixed-gap pair invisible to graphon and spectral states

The following construction gives a direct no-go theorem.

### Proposition 4.1

There are two sequences of zero-diagonal symmetric sign matrices
\(A_n,B_n\), for \(n=2^d\), such that
\[
\left\|\frac{A_n}{\sqrt n}\right\|_{\rm op}
+\left\|\frac{B_n}{\sqrt n}\right\|_{\rm op}=O(1),
\tag{4.1}
\]
\[
\left\|\frac{A_n-B_n}{\sqrt n}\right\|_{2,\tau}
\longrightarrow0,
\qquad
\|X\|_{2,\tau}^2=\frac1n\operatorname{Tr}(X^\top X),
\tag{4.2}
\]
but
\[
\limsup_nG_n(A_n)\le\frac12,\qquad
\liminf_nG_n(B_n)\ge\frac23.
\tag{4.3}
\]
Moreover every fixed signed graph density and every fixed-replica
normalized energy law has the same limit for the two sequences.

### Proof

Index rows and columns by \(\mathbb F_2^d\), and let
\[
W_{uv}=(-1)^{u\cdot v}.
\]
Then \(W=W^\top\) and \(W^2=nI\).  Let \(A\) be \(W\) with its diagonal
replaced by zero.  Hence
\[
\|A\|_{\rm op}\le\sqrt n+1
\]
and
\[
G_n(A)\le\frac{n\|A\|_{\rm op}/2}{n^{3/2}}
\le\frac12+o(1).
\tag{4.4}
\]
Also
\[
\sum_{u<v}a_{uv}=\frac n2.
\tag{4.5}
\]
Indeed \(W{\bf1}=n e_0\), while
\(\operatorname{Tr}W=\sum_u(-1)^{u\cdot u}=0\).

Let \(\Gamma\) be the graph of negative off-diagonal entries of \(A\).
By (4.5),
\[
|E(\Gamma)|=\frac{n(n-2)}4.
\]
Independently retain every edge of \(\Gamma\) with probability
\[
p_n=\frac{4}{3\sqrt n}.
\]
With probability tending to one, the retained graph \(F\) has
\[
|E(F)|=\left(\frac13+o(1)\right)n^{3/2}.
\tag{4.6}
\]
Matrix Bernstein gives
\[
\|\operatorname{Adj}(F)-p_n\operatorname{Adj}(\Gamma)\|_{\rm op}
=O(n^{1/4}\sqrt{\log n}),
\]
while
\[
\|p_n\operatorname{Adj}(\Gamma)\|_{\rm op}\le p_nn=O(\sqrt n).
\]
Thus there is a deterministic realization satisfying (4.6) and
\[
\|\operatorname{Adj}(F)\|_{\rm op}=O(\sqrt n).
\tag{4.7}
\]

Obtain \(B\) by flipping the entries on \(F\) from \(-1\) to \(+1\).
Then (4.1) follows from
\[
B-A=2\operatorname{Adj}(F).
\]
At the all-ones spin,
\[
H_B({\bf1})=H_A({\bf1})+2|E(F)|
=\left(\frac23+o(1)\right)n^{3/2},
\]
which proves (4.3).

Since each flipped undirected edge changes two matrix entries by
magnitude \(2\),
\[
\left\|\frac{A-B}{\sqrt n}\right\|_{2,\tau}^2
=\frac{8|E(F)|}{n^2}
=O(n^{-1/2}),
\]
proving (4.2).  On the operator-norm bounded set (4.1), telescoping
and Cauchy--Schwarz give, for every fixed \(r\),
\[
\left|
\frac1n\operatorname{Tr}(A/\sqrt n)^r
-\frac1n\operatorname{Tr}(B/\sqrt n)^r
\right|
\le
rC^{r-1}
\left\|\frac{A-B}{\sqrt n}\right\|_{2,\tau}
=o(1).
\tag{4.8}
\]
The same proof applies to every fixed trace polynomial.

For a fixed signed test graph with \(e\) edges, a union bound over
which test edge lands in \(F\) gives a difference in homomorphism
density at most
\[
O_e\left(\frac{|E(F)|}{n^2}\right)=O_e(n^{-1/2}).
\tag{4.9}
\]
Finally, for a uniform spin \(X\),
\[
H_B(X)-H_A(X)=2\sum_{\{i,j\}\in E(F)}X_iX_j.
\]
Orthogonality of distinct edge characters gives
\[
\mathbb E_X
\left|
\frac{H_B(X)-H_A(X)}{n^{3/2}}
\right|^2
=\frac{4|E(F)|}{n^3}
=O(n^{-3/2}).
\tag{4.10}
\]
For every fixed number of replicas, coupling them with the same spins
therefore makes the joint normalized energy laws converge in
Wasserstein distance.  This completes the proof.

### Consequences

1. \(G\) is discontinuous in the normalized spectral-microstate
   topology, even on a common operator-norm bounded set.
2. Every LDP state consisting of empirical spectral measure, finitely
   many trace moments, fixed signed-subgraph densities, or fixed
   replica distributions identifies the two sequences but cannot
   contract to \(G\).
3. Standard Wigner-traffic observables in the double-tree sector are
   also unchanged: for a fixed connected test graph in which each
   underlying simple edge occurs at least twice, the changed-edge
   bound after Wigner normalization is
   \(O(|E(F)|/n^2)=o(1)\).
4. The missing information is an exponentially sparse rank-one
   direction.  In normalized noncommutative trace,
   every \(P_x\) has trace \(1/n\) and disappears.  Retaining finitely
   many such infinitesimal projections is still insufficient; the
   supremum ranges over \(2^{n-1}\) of them.

Thus a second-order graph limit which retains only fixed traffic
coordinates is not enough.  It must retain the entire extremal support
of the Boolean projectors, as in Section 3, or an equivalent
exponential-complexity object.

## 5. Exact all-replica coupling under one-vertex gauge deletion

The support state from Section 3 does have a genuine projective
coupling which holds at every replica level simultaneously.

Given a uniform order-\(n\) signing \(A\), switch it so that all edges
incident to vertex \(n\) are \(+1\).  If \(B\) is the remaining
order-\((n-1)\) core, then the map
\[
A\longleftrightarrow
\left((a_{in})_{i<n},B\right),\qquad
b_{ij}=a_{ij}a_{in}a_{jn},
\tag{5.1}
\]
is a bijection.  Consequently \(B\) is itself a uniform signing.
Switching only relabels the Boolean projectors, so it leaves
\(\mathcal T_n(A)\) unchanged.

For a core spin \(u\in\{\pm1\}^{n-1}\) and last spin
\(\varepsilon\in\{\pm1\}\),
\[
H_A(u,\varepsilon)
=H_B(u)+\varepsilon\sum_{i<n}u_i.
\tag{5.2}
\]
Writing \(e_n=H_A/n^{3/2}\) and
\(e_{n-1}=H_B/(n-1)^{3/2}\), this is
\[
e_n(u,\varepsilon)
=\left(\frac{n-1}{n}\right)^{3/2}e_{n-1}(u)
+\frac{\varepsilon\sum_i u_i}{n^{3/2}}.
\tag{5.3}
\]
Uniformly over every signing and every spin,
\[
|e_n(u,\varepsilon)-e_{n-1}(u)|=O(n^{-1/2}).
\tag{5.4}
\]
The first term costs \(O(n^{-1})|e_{n-1}|\), and the universal
bound \(|e_{n-1}|\le O(\sqrt n)\) makes this \(O(n^{-1/2})\);
the affine field has the same bound.

For two extended spins,
\[
q_n((u,\varepsilon),(v,\delta))
=\frac{n-1}{n}q_{n-1}(u,v)+\frac{\varepsilon\delta}{n},
\]
so
\[
|q_n^2-q_{n-1}^2|\le\frac4n.
\tag{5.5}
\]
Equip every \(X_k\) with the sup metric and \(\mathfrak T\) with
\[
d_{\mathfrak T}(S,T)
=\sum_{k\ge1}2^{-k}
\min\{1,d_H(S_k,T_k)\}.
\tag{5.6}
\]
Every \(k\)-tuple in the order-\(n\) landscape restricts to a core
\(k\)-tuple, and every core tuple can be extended, for example with
all last spins \(+1\).  Equations (5.4)--(5.5) therefore prove the
deterministic, all-level estimate
\[
\boxed{
d_{\mathfrak T}(\mathcal T_n(A),\mathcal T_{n-1}(B))
\le\frac C{\sqrt n}.}
\tag{5.7}
\]
In particular, the laws \(\mu_n\) and \(\mu_{n-1}\), coupled by
(5.1), are exponentially equivalent at speed \(n^2\): for every
fixed \(\eta>0\), the event that their distance exceeds \(\eta\) is
empty for all sufficiently large \(n\).

This is stronger than fixed-replica convergence, but it still does not
force a unique LDP.  Adjacent exponential equivalence only says that
two neighboring terms have the same local asymptotics.  Point masses
at
\[
z_n=e^{i\log n}
\]
on the unit circle satisfy \(d(z_n,z_{n-1})=O(1/n)\), hence the
analogous adjacent exponential equivalence at every speed, while
having a continuum of subsequential limits and no full LDP.

Iterating (5.7) through \(h\) deletions gives only
\[
d_{\mathfrak T}(\mathcal T_n,\mathcal T_{n-h})
=O\left(\frac h{\sqrt n}\right)
\tag{5.8}
\]
when \(h=o(n)\).  Thus the coupling is effective only for
\(h=o(\sqrt n)\).  At proportional scales its accumulated affine
fields are of leading order.  This is the precise point at which the
all-replica projective coupling stops short of convergence.

## 6. Restriction inequality: equivalent to the existing Shearer bound

The natural Finner calculation gives a useful normalization check, but
it is **exactly equivalent** to the previously derived
Shearer/relative-entropy restriction inequality; it is not an
additional result.

Let
\[
L_n(\lambda)=\mathbb E_{U_n}e^{-\lambda M(A)},\qquad
\Psi_n(\beta)=\frac1{n(n-1)}
\log L_n(\beta\sqrt n).
\]
For \(m<n\), define
\[
\beta_{m,n}
=\beta\sqrt{\frac mn}\frac{m-1}{n-1}.
\tag{6.1}
\]
Then
\[
\boxed{\Psi_n(\beta)\le\Psi_m(\beta_{m,n}).}
\tag{6.2}
\]

To prove it, let \(\mathcal S\) be all \(m\)-subsets of \([n]\),
\[
K=\binom nm,\qquad c=\binom{n-2}{m-2}.
\]
For every fixed signing,
\[
M(A_S)\le M(A).
\tag{6.3}
\]
Indeed, extend a maximizing spin on \(S\) by independent random signs
off \(S\); the conditional mean of the full Hamiltonian is the
restricted Hamiltonian.  Therefore
\[
e^{-\lambda M(A)}
\le
\prod_{S\in\mathcal S}
e^{-\lambda M(A_S)/K}.
\]
Every edge coordinate occurs in exactly \(c\) factors.  Generalized
Hölder (Finner's inequality) gives
\[
L_n(\lambda)
\le
\prod_{S\in\mathcal S}
\left(
\mathbb E e^{-\lambda c M(A_S)/K}
\right)^{1/c}
=
L_m\left(\lambda\frac cK\right)^{K/c}.
\tag{6.4}
\]
Since
\[
\frac Kc=\frac{n(n-1)}{m(m-1)},
\]
substituting \(\lambda=\beta\sqrt n\) gives (6.2).

For \(m/n\to\alpha\), the temperature changes as
\[
\beta_{m,n}=\beta\alpha^{3/2}+o(1).
\tag{6.5}
\]
This is the same scale drift that obstructs ordinary free-energy
subadditivity.  There is no reverse Finner inequality: even for two
independent signs \(z_1,z_2\), the centrally symmetric functions
\[
\mathbf1_{\{z_1=z_2\}},\qquad
\mathbf1_{\{z_1=-z_2\}}
\]
have product expectation \(0\), while the product of expectations is
\(1/4\).  Switching symmetry therefore supplies no positive
association.

At zero entropy, (6.2) reduces only to the restriction monotonicity
\[
F(n)\ge F(m),
\qquad
\frac{F(n)}{n^{3/2}}
\ge
\left(\frac mn\right)^{3/2}
\frac{F(m)}{m^{3/2}},
\]
which permits slowly oscillating normalized sequences.  Hence (6.2)
is exact but not a convergence proof.

## 7. What remains

The full speed-\(n^2\) route is now sharply formulated.

* A bulk graphon or noncommutative microstate LDP is too coarse by
  Proposition 4.1.
* The compact support-projector state \(\mathfrak T\) is fine enough
  and automatically exponentially tight; (5.7) couples all its
  replica levels under adjacent deletion.
* A full LDP on \(\mathfrak T\), or merely uniqueness of
  \(\inf(I+\beta G)\) across all subsequential rate functions for an
  unbounded set of \(\beta\)'s, proves the desired limit.
* Proving that uniqueness requires an order-composition theorem for
  extremal support sets.  The restriction inequality below loses the factor
  \(\alpha^{3/2}\) in temperature and has no reverse inequality.

An adequate next state should probably combine:

1. the support sets \(\mathcal T_{n,k}\), to retain zero-entropy
   resonant projectors;
2. logarithmic multiplicities of neighborhoods of those support sets,
   to retain the \(e^{\Theta(n)}\) spin complexity;
3. switching-class entropy at speed \(n^2\).

The missing lemma is then not exponential tightness, which is free,
but a balanced-block transport or cavity identity proving that the
resulting subsequential rate functional is unique.  Without such a
lemma, a “traffic LDP” simply assigns rate \(+\infty\) to the sparse
phase or merges it with the bulk, and either choice loses the
observable that controls the minimum.
