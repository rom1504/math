# Independent audit: proportional random conference thinning

Date: 2026-09-06. This reconstructs the theorem in
`continued_convergence_proportional_thinning_2026_09_06.md` from its finite
graph and Gaussian-rounding arguments. No spectral-limit or eigenvector
theorem is imported as a substitute for the required diagonal control.

## 1. Verdict and exact scope

For any unbounded family of symmetric hollow conference signings C_N,
and a uniform subset S of size `m=floor(pN)`, the proof supports

```math
\liminf_{p\downarrow0}\liminf_{N\to\infty}
\mathbb E_S\frac{Q(C_N[S])}{m^{3/2}}\ge\frac2\pi. \tag{1}
```

Already one fixed degree-three polynomial supplies
`2cos(pi/5)/pi>.5` in the iterated limit. Therefore for some fixed positive
retention proportion the average normalized child cap exceeds one-half,
whereas every conference parent has normalized cap at most one-half.

This disproves the universal AVERAGE proportional-transfer inequality at
the original `N^(3/2)` scale. It does not prove that every subset is bad,
that the parents are exact minimizers, or that the original minimum sequence
does not converge. Exceptional selection and minimizer-specific transfer
remain open.

## 2. The rooted signed-walk estimate really uses conference orthogonality

For a connected loopless Eulerian multigraph with v vertices, e edge
occurrences, and one fixed root label, let T_G be its injective signed sum.
The required bound is `|T_G|<=C_G N^(e/2)`.

Absolute counting proves this whenever `v<=e/2+1`. Above that threshold,
there is a nonroot degree-two vertex: otherwise the sum of degrees is at
least `4(v-1)+2>2e`. This also proves existence of such a vertex in the
equality case, except for the trivial single root.

A doubled leaf can be summed with the EXACT multiplier `N-v+1`, leaving
v-1 vertices and e-2 edge occurrences. A degree-two vertex u with distinct
neighbors a,b instead satisfies

```math
\sum_{z\notin\phi(V\setminus\{u\})}C_{\phi(a),z}C_{z,\phi(b)}
=-\sum_{w\in V\setminus\{u\}}
 C_{\phi(a),\phi(w)}C_{\phi(w),\phi(b)}. \tag{2}
```

This uses `(C²)_ab=0` for a!=b, not an approximate discrepancy estimate.
The a,b terms vanish by hollowness. Every surviving term has one fewer
vertex and the same e edge occurrences, with a-u-b replaced by a-w-b.
The graph remains connected and Eulerian. A degree-two vertex in a
connected Eulerian graph lies on a cycle, so its removal leaves its two
neighbors connected by the rest of that cycle. Parallel edges retain their
full occurrence count for the normalization.

Induction on the number of vertices proves the bound. At
`v=e/2+1`, doubled-leaf removal preserves equality. Unless the graph is
a doubled tree, a distinct-neighbor removal eventually occurs and reduces
the free-label count by one without reducing e. Absolute counting then
gives `O(N^(e/2-1))`. For a doubled tree the exact sum is
`(N-1)_(v-1)`. This verifies the essential one-power improvement.

## 3. Fixed-cardinality rooted moments and products

Normalize the compression by `B=C_N[S]/sqrt(m)`. Conditional on its root
being selected, the exact inclusion weight for an injective v-vertex graph
is

```math
\frac{(m-1)_{v-1}}{(N-1)_{v-1}}.
```

For fixed p, its rooted moment contribution is that weight times
`m^(-e/2)T_G`. Graphs below the equality threshold vanish with N. Equality
leaves doubled trees and therefore semicircle Catalan moments. Graphs above
the threshold have bound `O(p^(v-1-e/2))`, which vanishes as p decreases.
There are finitely many graphs for every fixed moment order.

For the second moment of a diagonal polynomial, wedge its two closed walks
at the root. On a tree, an individual closed walk uses every edge an even
number of times. Hence two walks cannot share an edge when the COMBINED
graph has every edge doubled. Their edge-disjoint connected subtrees can
meet only at the root. The leading count therefore factors into the two
Catalan counts. This supplies the actual rooted diagonal concentration

```math
\lim_{p\downarrow0}\limsup_{N\to\infty}\max_i
\mathbb E\left[|(f(B))_{ii}-\tau_{\rm sc}(f)|^2
\mid i\in S\right]=0 \tag{3}
```

for each fixed polynomial f. Empirical spectral convergence alone would
not have supplied (3).

## 4. Gaussian covariance trimming and all constants

Fix a polynomial g, put `K=g(B)²`, `nu=tau_sc(g²)>0`, and
`lambda=tau_sc(t g²)`. For fixed delta>0 let P retain those coordinates
with `K_ii<=nu+delta`. The matrix `PKP/(nu+delta)` is positive semidefinite
with diagonal at most one; independent diagonal Gaussian noise fills the
deficit. The resulting unit-variance Gaussian has off-diagonal correlations
`R_ij=(PKP)_ij/(nu+delta)`.

Its sign vector has EXPECTED ORIGINAL half-energy

```math
\frac{\sqrt m}{\pi}\sum_{ij}B_{ij}\arcsin(R_{ij}). \tag{4}
```

The factor is correct: half-energy contributes 1/2 and Gaussian sign
covariance contributes 2/pi. The uniform inequality
`|arcsin(r)-r|<=(pi/2-1)r²` follows by summing the positive power-series
coefficients and using `|r|^(2k+1)<=r²` for k>=1. Because
`sqrt(m)|B_ij|=1` off the diagonal, the unnormalized error in (4) is at
most `C_(nu,delta) Tr(K²)`. For fixed p, bounded polynomial moments make
its normalized expectation `O_p(m^-1/2)`.

The trimming trace error also reconstructs exactly. Split
`K-PKP=(I-P)K+PK(I-P)`. Applying row Cauchy--Schwarz separately to the
two traces gives

```math
|\operatorname{Tr}(BK)-\operatorname{Tr}(BPKP)|
\le2\sqrt{(m-1)/m}\sum_{i\ {m discarded}}\sqrt{(K²)_{ii}}. \tag{5}
```

Averaged root expectations are valid since
`E[(1/m)sum_(i in S) a_i]=(1/N)sum_i E[a_i|i in S]`.
Equation (3) bounds the expected discarded proportion by a quantity tending
to zero as N then p tend to their limits. Another Cauchy--Schwarz inequality
and the g^4 rooted moment bound make (5), divided by m, negligible. Taking
N, then p, then delta to their limits gives `lambda/(pi nu)`.

Thus neither a small maximum covariance entry nor an operator-norm loss
from normalizing each row is assumed. Both common gaps are avoided by
the explicit trace-square bound and diagonal filling.

## 5. Polynomial optimization and limit order

In the semicircle orthonormal basis `U_j(t/2)`, multiplication by t has
the finite path adjacency matrix on degrees 0,...,d. Its top eigenvalue is
`2cos(pi/(d+2))`, with the displayed sine coefficients in the source
artifact. Each FIXED polynomial gives that quotient divided by pi.
Only after the N,p,delta limits does d increase, yielding 2/pi.

The resulting fixed positive p obstruction requires only d=3 and a small
enough fixed delta, followed by a small enough fixed p. No polynomial
degree grows with matrix order. No extension from special conference
orders to all original minimizer orders is asserted or required.

The fixed-cardinality theorem is therefore verified. The Bernoulli version
uses the same graph proof with inclusion weight p^(v-1); exponential
binomial concentration and the conference norm bound handle normalization
by the random selected size. It is not necessary to couple adjacent-sized
matrices with an unproved negligible energy difference.
