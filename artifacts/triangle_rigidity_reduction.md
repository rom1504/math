# Triangle Rigidity in the Deep-Hole Certificate

## 1. Setup

Let

\[
E=\binom{[n]}2,\qquad N=|E|,
\]

and encode a sign vector by a binary edge cochain:

\[
s_e=(-1)^{\alpha_e},\qquad \alpha\in\mathbb F_2^E.
\]

The augmented cut code is

\[
\mathcal C_n
=
\left\{
c(t,z)_{ij}=t+z_i+z_j:
t,z_i\in\mathbb F_2
\right\}.
\]

In sign notation these are precisely

\[
v_{ij}=(-1)^t x_i x_j.
\]

For a fixed signing \(a\), let

\[
h_v=a\cdot v,\qquad
M=\max_{v\in\mathcal C_n}h_v,
\qquad
g_v=\frac{M-h_v}{2}.
\]

If

\[
N_v=\{e:a_ev_e=-1\},
\]

then

\[
|N_v|=d+g_v,
\qquad
d=\frac{N-M}{2}.
\]

If \(a\) is a global minimizer and \(S\subseteq E\) is any edge-flip
set, the exact deep-hole certificate is

\[
\exists v:\qquad
2|S\cap N_v|\ge |S|+g_v. \tag{1}
\]

Equivalently,

\[
|S\mathbin\triangle N_v|\le d. \tag{2}
\]

The purpose of this note is to determine exactly what additional force
comes from the triangle identities obeyed by \(v\).

## 2. The exact local rank of the triangle identities

Let \(H=(V,S)\) be a graph with no isolated vertices, with \(c(H)\)
connected components. Write

\[
\beta(H)=|S|-|V|+c(H)
\]

for its cycle rank.

### Theorem 2.1

The restriction of the augmented cut code to \(S\) has binary dimension

\[
\dim(\mathcal C_n|_S)
=
\begin{cases}
|V|-c(H),&H\text{ is bipartite},\\[2mm]
|V|-c(H)+1,&H\text{ is nonbipartite}.
\end{cases}
\]

Consequently its codimension in \(\mathbb F_2^S\) is

\[
q(H)
=
\begin{cases}
\beta(H),&H\text{ is bipartite},\\[2mm]
\beta(H)-1,&H\text{ is nonbipartite}.
\end{cases} \tag{3}
\]

#### Proof

The ordinary cut restrictions form the image of the binary
vertex-edge incidence map

\[
z\longmapsto (z_i+z_j)_{ij\in S}.
\]

Its rank is \(|V|-c(H)\). Adding the global bit \(t\) adds the all-one
edge vector \(\mathbf1_S\). This vector already lies in the ordinary cut
space exactly when the equations

\[
z_i+z_j=1\qquad(ij\in S)
\]

are solvable, which is equivalent to every component of \(H\) being
bipartite. This proves the dimension formula and hence (3). \(\square\)

### Corollary 2.2

The dual local constraints are exactly the edge sets \(F\subseteq S\)
such that

1. every vertex has even degree in \(F\), and
2. \(|F|\) is even.

Indeed, orthogonality to all stars gives the Eulerian condition, while
orthogonality to the global bit gives even cardinality.

Thus:

- every forest has \(q(H)=0\);
- every matching and every star sees the full local cube;
- a single triangle, and more generally a connected odd unicyclic
  graph, also sees the full local cube;
- the first genuine constraint is a \(4\)-cycle;
- two independent odd cycles give one constraint: their odd-cycle
  products must agree.

This is the precise local content of

\[
v_{ij}v_{jk}v_{ki}=\text{one global constant}.
\]

A single triangle does not constrain an augmented cut, because the
global sign supplies its missing parity bit.

### Corollary 2.3: quantitative repair

Every binary pattern on \(S\) is at Hamming distance at most \(q(H)\)
from \(\mathcal C_n|_S\).

To see this, take a full-rank \(q(H)\times |S|\) parity-check matrix and
choose \(q(H)\) pivot columns. Any syndrome can be canceled by changing
only a subset of those pivot coordinates.

Therefore the triangle identities alone can change a majority score on
\(S\) by at most \(2q(H)\).

For a subcritical random sparse graph with \(|S|=\alpha n\),
\(\alpha<1/2\), the cycle excess \(q(H)\) is tight in probability.
Hence such a test support contains only \(O_{\mathbb P}(1)\) independent
triangle-consistency bits. To obtain \(\Theta(n)\) independent local
constraints one needs \(\Theta(n)\) cycle excess.

## 3. Stars give a tautological certificate

If \(S\) is a cut, in particular a full star, then flipping \(S\) is a
Seidel switching and leaves \(M(a)\) unchanged for every signing \(a\),
whether or not \(a\) is globally minimizing.

More explicitly, take a top state \(v_0\), so \(N_{v_0}\) has size
\(d\). Multiplying \(v_0\) by the cut codeword supported on \(S\)
produces a state \(v\) with

\[
N_v=N_{v_0}\mathbin\triangle S
\]

and

\[
g_v=|S|-2|S\cap N_{v_0}|.
\]

Substitution in (1) gives equality. Thus star tests contain no
deep-hole information.

## 4. Exact two-graph/cofilling reduction

Define the simplicial coboundary

\[
(\delta\alpha)_{ijk}
=
\alpha_{ij}+\alpha_{jk}+\alpha_{ki}
\quad\text{in }\mathbb F_2.
\]

Then

\[
\mathcal C_n=\delta^{-1}(\langle\mathbf1\rangle):
\]

an edge vector is an augmented cut if and only if its triangle
coboundary is constant.

The image of \(\delta\) consists of the triangle vectors \(\tau\) such
that every tetrahedron contains an even number of selected triangles:

\[
\tau_{ijk}+\tau_{ij\ell}+\tau_{ik\ell}+\tau_{jk\ell}=0. \tag{4}
\]

These are precisely two-graphs. Conversely, (4) gives an edge filling
explicitly: fix vertex \(1\), set \(\alpha_{1i}=0\) and
\(\alpha_{ij}=\tau_{1ij}\) for \(i,j\ne1\), then use (4).

The augmented covering radius is therefore

\[
\rho(\mathcal C_n)
=
\max_{\tau\ {\rm satisfying}\ (4)\ {\rm mod}\ \mathbf1}
\min\left\{
|\alpha|:\delta\alpha\in\{\tau,\tau+\mathbf1\}
\right\}. \tag{5}
\]

This is the maximum edge-cofilling norm of a two-graph, with a
two-graph identified with its complement. In this notation,

\[
M_n=N-2\rho(\mathcal C_n). \tag{6}
\]

In signed-graph terminology, the unaugmented version is the maximum
**line index of imbalance** (or maximum frustration index) over
switching classes, equivalently the covering radius of the cocycle
code. This correspondence is treated by Solé and Zaslavsky, *A Coding
Approach to Signed Graphs*, SIAM J. Discrete Math. 7 (1994),
DOI 10.1137/S0895480189174374. The present augmented problem further
identifies every switching class with its edge complement.

Moreover, (1) for every \(S\) says exactly that every translated edge
cochain \(\alpha+S\) has an augmented-cut representative within distance
\(\rho\). Hence:

> The full majority-with-gap certificate plus the triangle identities
> is exactly the original covering-radius/maximal-cofilling problem.

The triangle identities are not an independent extra hypothesis after
the code is identified; they are the defining equations of the quotient
in (5).

### Rooted self-similarity

The tetrahedron identities do not make the switching invariant
low-dimensional in the asymptotic sense relevant here. Switch a signing
so that

\[
a_{1i}=1\qquad(i>1).
\]

Then its rooted triangle signs are

\[
\tau_{1ij}=a_{1i}a_{ij}a_{j1}=a_{ij},
\]

and these \(\binom{n-1}{2}\) signs are completely arbitrary. Every other
triangle sign is recovered from the tetrahedron relation. After fixing
the global spin by \(x_1=1\), the energy becomes

\[
\sum_{i=2}^n x_i+
\sum_{2\le i<j\le n}\tau_{1ij}x_ix_j. \tag{7}
\]

Thus rooted triangle coordinates reparameterize an arbitrary signing on
\(K_{n-1}\), with only an added all-\(+1\) linear field. This is the
exact affine recursion already implicit in the original problem.
Triangle consistency therefore does not create a pseudorandomness or
regularity assumption; a global argument in these coordinates must still
handle an arbitrary smaller signing.

## 5. The precise missing nonlocal statistic

For a support graph \(H=(V,S)\), let

\[
\Gamma_H(y)
=
\min\{g_v:N_v|_S=y\},
\qquad
y\in \alpha|_S+\mathcal C_n|_S.
\]

Then the restriction of the full certificate to flip sets
\(T\subseteq S\) is exactly

\[
\forall T\subseteq S,\qquad
\max_y
\left(
2|T\cap y|-|T|-\Gamma_H(y)
\right)
\ge0. \tag{8}
\]

The triangle identities determine the affine domain of \(y\), whose
codimension is (3). They give no control on the conditional extension
gap \(\Gamma_H(y)\).

Controlling \(\Gamma_H\) is the missing inequality. On a spanning sparse
support it already records essentially the complete quadratic energy
landscape, so (8) is not a scalar local consequence of triangle parity.

## 6. A necessary active-layer counting inequality

Let

\[
A_j=\#\{v:g_v=j\}.
\]

For a fixed state in layer \(j\), its mismatch set has size \(d+j\).
Among all \(k\)-edge flip sets it certifies exactly

\[
C_{N,d}(k,j)
=
\sum_{\ell\ge\lceil(k+j)/2\rceil}
\binom{d+j}{\ell}
\binom{N-d-j}{k-\ell}. \tag{9}
\]

Double-counting all certified \(k\)-sets yields the necessary condition

\[
\boxed{
\binom Nk
\le
\sum_{j=0}^{k}
A_j C_{N,d}(k,j).
} \tag{10}
\]

This is the strongest conclusion available from majority coverage and
the global layer weights without using overlaps between different cut
states.

If \(M=cn^{3/2}\), \(k=\alpha n\), and \(j=O(\sqrt n)\), the probabilities
\(C_{N,d}(k,j)/\binom Nk\) are constant-order hypergeometric tails.
Thus (10) can force only a constant number of states in the genuinely
near-active window \(g=O(\sqrt n)\). It cannot by itself improve the
\(n^{3/2}\) constant. Taking \(k\gg n\) eventually gives exponential
requirements, but then it uses broad energy layers and becomes the
ordinary Hamming sphere-covering count.

Any improvement must use high-order overlap information between the
certifying events. For augmented cuts those overlaps are governed by the
full cut association scheme, so this is again a high-order
covering-radius problem rather than a local triangle estimate.

## 7. Exact order-\(10\) diagnostic

For an exact optimizer with

\[
N=45,\qquad M=13,\qquad d=16,
\]

the positive augmented energy layers have sizes

\[
A_0=40,\qquad A_1=80,\qquad A_4=192,\qquad A_5=200.
\]

Direct enumeration gives:

| Test family | Number | Top-layer result | Next layer |
|---|---:|---|---|
| triangles | 120 | all certified; margins \(1\) or \(3\) | unnecessary |
| \(4\)-cycles | 630 | all certified; margins \(0,2,4\) | unnecessary |
| full stars | 10 | every star misses by exactly \(1\) | \(g=1\) certifies at equality |
| perfect matchings | 945 | 18 fail at \(g=0\) | all 18 rescued by \(g=1\) |

For one fixed triangle, all eight local patterns occur with conditional
minimum gaps

\[
7\text{ patterns at }g=0,\qquad1\text{ pattern at }g=1.
\]

For a representative perfect matching, its \(32\) patterns have
conditional minimum gaps

\[
15\text{ at }g=0,\qquad16\text{ at }g=1,\qquad1\text{ at }g=4.
\]

For a full star, all \(512\) patterns occur, but their conditional
minimum gaps are

\[
\begin{array}{c|rrrrrr}
g&0&1&4&5&8&9\\ \hline
\#&40&80&169&153&47&23.
\end{array}
\]

The first top-layer failures occur on a forest (perfect matchings), where
triangle consistency imposes no local restriction. Cyclic tests are
easier in this model, not harder. This supports the structural conclusion
that the useful missing input is the conditional gap profile
\(\Gamma_H\), not another finite triangle-parity consequence.

## 8. Verdict

The triangle route does produce a rigorous structural theorem and a
precise stopping point:

1. the proposed matching, star, and single-triangle tests cannot exploit
   triangle rigidity at all;
2. bounded-cycle-excess supports contain only boundedly many independent
   constraints;
3. the complete triangle system plus the all-flip certificate is exactly
   the two-graph cofilling/augmented-cut covering-radius problem;
4. progress requires a genuinely nonlocal estimate on conditional
   extension gaps or on high-order overlaps of cut-code balls.

No lower bound beyond the verified field-plus-spin constant follows from
these local tests alone.
