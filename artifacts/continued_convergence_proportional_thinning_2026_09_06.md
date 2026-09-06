# Continued convergence campaign: proportional thinning (2026-09-06)

Notation throughout this note is

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,\qquad M_n=\min_AQ(A).
```

## 1. Independent assessment frozen before reading historical assessments

The original convergence problem led to three concrete architectures.

1. **Sharp insertion at minimizers.** A theorem
   `M_{n+1} <= (1+1/n)^{3/2} M_n + O(n^{1/2-delta})`, for some
   fixed `delta>0`, would give convergence by summing normalized errors.
   A more structural first theorem would find, for at least one exact
   minimizer `A` of each large order, a signed row `b` satisfying
   `max_x (|H_A(x)|+|b.x|) <= (1+1/n)^{3/2}M_n+O(n^{1/2-delta})`.
   The actual finite falsifier of the zero-defect, all-order version is
   the order-six minimizer: its best extension is 9, whereas `M_6=5`.
   This does not falsify an asymptotic minimizer theorem.
2. **Proportional extraction.** An appropriately uniform theorem selecting
   near-minimizers and principal subsets with
   `Q(A[S]) <= (|S|/n)^{3/2} Q(A) + o(n^{3/2})` would transfer good
   orders downward. A universal average-over-subsets version is a
   concrete, stronger test. The theorem below asymptotically falsifies
   that stronger version, at a fixed positive retained proportion.
3. **Algebraic amplification.** For each of two multiplicatively
   independent integer multipliers, a minimizer-specific lift obeying
   `M_{rn} <= r^{3/2}M_n + O(n^{3/2-delta})` would imply convergence.
   The first explicit candidate was
   `L(A)=[[A,A+I],[A+I,-A]]`. After freezing it, archive inspection found
   that this exact candidate, its finite tests, and the two-multiplier
   convergence criterion had already been treated in
   `fresh_limit_algebra_2026_09_05.md`. They are not new progress here.

The exact split identity remains

```math
Q\!\begin{pmatrix}A&C\\C^T&D\end{pmatrix}
=\max_{x,y}\bigl(|H_A(x)+H_D(y)|+|x^TCy|\bigr).
```

In particular, global spin reversal cannot cancel the bridge.

## 2. Result and exact scope

**Verified by an independent reconstruction, pending final root audit.**
Let `C_N` be any family of symmetric conference matrices with unbounded
orders:

```math
(C_N)_{ii}=0,\qquad (C_N)_{ij}\in\{\pm1\}\ (i\ne j),\qquad
C_N^2=(N-1)I.
```

Take `m=floor(pN)`, and choose `S` uniformly among the `m`-vertex
subsets. The same conclusion holds for independent Bernoulli-`p`
selection, with the ratio set to zero on the empty subset. Then

```math
\boxed{
\liminf_{p\downarrow0}\ \liminf_{N\to\infty}
\mathbb E_S\frac{Q(C_N[S])}{|S|^{3/2}}
\ \ge\ \frac2\pi.
} \tag{T}
```

Both lower limits may be restricted to the available conference orders.
The estimates are uniform over the choice of conference matrix at each
order. In particular, there exists a fixed `p>0` for which

```math
\liminf_{N\to\infty}
\mathbb E_S\frac{Q(C_N[S])}{|S|^{3/2}}>\frac12,
```

whereas `Q(C_N)/N^{3/2} <= sqrt(1-1/N)/2`. Thus universal proportional
**average** contraction at exponent `3/2` is false even inside the exact
conference class, at the original `N^{3/2}` scale.

This does not prove nonconvergence of `M_n/n^{3/2}`. It does not rule out
an exceptional selected subset, or a transfer theorem that uses global
minimality in an essential way. The distinction is material: a positive
average lower bound does not imply that every principal subset is bad.

Already a degree-three polynomial gives the explicit iterated-limit
coefficient

```math
\frac{2\cos(\pi/5)}\pi=0.515036\ldots>\frac12.
```

## 3. Rooted conference-walk estimate

Let `G` be a finite connected Eulerian multigraph without loops, with
`v` vertices, `e` edge occurrences, and one distinguished root. Parallel
edges are retained as occurrences. For a fixed root label `i`, let

```math
T_{G,N}(i)=\sum_{\substack{\phi:V(G)\hookrightarrow[N]\\
                         \phi(\mathrm{root})=i}}
\prod_{ab\in E(G)}(C_N)_{\phi(a),\phi(b)}.
```

For each fixed `G`, uniformly in `N`, `i`, and the conference matrix,

```math
|T_{G,N}(i)|=O_G(N^{e/2}). \tag{3.1}
```

At `v=e/2+1`, the sharper statement is

```math
T_{G,N}(i)=
\begin{cases}
(N-1)_{v-1},&G\text{ is a tree with every edge doubled},\\
O_G(N^{e/2-1}),&\text{otherwise}.
\end{cases} \tag{3.2}
```

Here the first line is exact because every matrix entry is squared.

**Proof.** If `v<=e/2+1`, absolute counting gives (3.1). Otherwise there
is a degree-two vertex other than the root. More generally, such a
vertex exists in the equality case as well, unless the graph is the
single root: if every nonroot degree were at least four, Eulerian
connectivity would give `2e>=4(v-1)+2`, a contradiction.

If the two neighbors coincide, the vertex is a doubled leaf. Summing
its label removes it and its two edge occurrences, with the exact
factor `N-v+1`. If its neighbors `a,b` are distinct, conference
orthogonality gives

```math
\sum_{z\notin\phi(V(G)\setminus\{u\})}
 C_{\phi(a),z}C_{z,\phi(b)}
=-\sum_{w\in V(G)\setminus\{u\}}
 C_{\phi(a),\phi(w)}C_{\phi(w),\phi(b)}.
```

Terms with `w=a` or `w=b` vanish. Each other term is the injective sum
of a connected Eulerian graph with one fewer vertex and the same
number of edge occurrences: replace `a-u-b` by `a-w-b`.
Induction on `v` proves (3.1). Connectivity is preserved because a
degree-two vertex with distinct neighbors in an Eulerian connected
graph lies on a cycle, so deleting that vertex leaves its remaining
vertices connected.

For (3.2), doubled-leaf removals preserve `v=e/2+1`. Unless the graph
is a doubled tree, a distinct-neighbor removal eventually occurs.
It reduces the vertex count without reducing `e`, after which absolute
counting saves a full factor `N`. Restoring any previously removed
doubled leaves proves the stated error. This proof also covers repeated
edges; loops created during substitution contribute zero.

## 4. Rooted semicircle moments after small proportional thinning

Put

```math
B=B_{N,p}=\frac{C_N[S]}{\sqrt m}.
```

Condition on the root `i` belonging to `S`. For each fixed polynomial
`f`, with `tau_sc` integration against the variance-one semicircle law,

```math
\lim_{p\downarrow0}\limsup_{N\to\infty}\max_i
\mathbb E\left[
\left|(f(B))_{ii}-\tau_{\rm sc}(f)\right|^2
\,\middle|\,i\in S\right]=0. \tag{4.1}
```

To verify this directly, expand a diagonal power into closed walks,
then partition positions according to equality of labels. For a
rooted walk graph with `v` distinct vertices and `e` edge occurrences,
the term is

```math
\frac{(m-1)_{v-1}}{(N-1)_{v-1}}m^{-e/2}T_{G,N}(i). \tag{4.2}
```

At fixed `p`, the inclusion factor is `p^{v-1}(1+O_G(1/N))`, and
`m/(pN)->1`. Terms with `v<e/2+1` vanish by absolute counting. At
`v=e/2+1`, (3.2) leaves precisely doubled trees. Terms with
`v>e/2+1` are `O_G(p^{v-1-e/2})` and tend to zero as `p` tends to zero.
The doubled-tree counts for a single closed walk are the Catalan
moments of the semicircle law.

For a product of two diagonal powers, wedge their walks at the root.
The same argument applies. A closed walk on a tree traverses every
edge an even number of times. Therefore, when the combined graph is
a doubled tree, the two walks cannot share an edge; their two subtrees
meet only at the root. The leading count is the product of the two
Catalan counts. This proves the second-moment assertion (4.1), not
merely convergence of the empirical spectral distribution.

All arguments involve finitely many graphs for each fixed polynomial,
so the iterated limits and their uniformity are legitimate. No growing
polynomial degree is used before taking `N` and `p` limits.

## 5. Gaussian covariance trimming without an operator-norm loss

Fix a real polynomial `g`, and write

```math
K=g(B)^2,\qquad \nu=\tau_{\rm sc}(g^2)>0,\qquad
\lambda=\tau_{\rm sc}(t g(t)^2).
```

Fix `delta>0`. Retain coordinates for which `K_ii<=nu+delta`; let `P`
be their diagonal projection. There is a centered Gaussian vector
with unit variances and off-diagonal covariance

```math
R_{ij}=\frac{(PKP)_{ij}}{\nu+\delta}\qquad(i\ne j). \tag{5.1}
```

Indeed use covariance `PKP/(nu+delta)`, then add independent Gaussian
noise to fill every diagonal to one. The added noise also gives unit
independent coordinates at the discarded locations.

The expected energy of the signs of this Gaussian vector is

```math
\frac{\sqrt m}\pi\sum_{i,j}B_{ij}\arcsin R_{ij}.
```

The bound

```math
|\arcsin r-r|\le(\pi/2-1)r^2\qquad(|r|\le1)
```

shows that its difference from
`sqrt(m) Tr(BPKP)/(pi(nu+delta))` is at most a constant times
`Tr(K^2)`. After dividing by `m^{3/2}`, its expectation tends to zero
for fixed `p`; bounded fixed-polynomial moments suffice. There is no
unproved maximum-entry delocalization requirement here.

Each row of `B` has squared norm `(m-1)/m`. Hence

```math
\left|\operatorname{Tr}(BK)-\operatorname{Tr}(BPKP)\right|
\le2\sqrt{\frac{m-1}{m}}
\sum_{i\ {m discarded}}\sqrt{(K^2)_{ii}}. \tag{5.2}
```

Apply Cauchy--Schwarz, (4.1) for `g^2`, and the bounded moment for
`g^4`. The expected right side divided by `m` tends to zero as
`N->infinity`, then `p->0`, for each fixed `delta`.

Using (4.1) again for `t g(t)^2`, we obtain

```math
\liminf_{p\downarrow0}\liminf_{N\to\infty}
\mathbb E\frac{Q(C_N[S])}{m^{3/2}}
\ge\frac{\lambda}{\pi(\nu+\delta)}.
```

Let `delta` decrease to zero. Thus every fixed polynomial gives the
lower bound `lambda/(pi nu)`.

For Bernoulli selection the proof is identical with inclusion weight
`p^{v-1}` and normalization `sqrt(pN)`. Binomial concentration and the
conference bound `Q(C_N[S]) <= |S| sqrt(N-1)/2` then permit the final
normalization by `|S|^{3/2}`. This is not a coupling argument claiming
that vertex-by-vertex energy differences are negligible; the fixed-size
statement was proved directly by its own inclusion probabilities.

## 6. Polynomial choice

Let `P_j(t)=U_j(t/2)` be the semicircle orthonormal polynomials. In
their basis, multiplication by `t` has ones on both adjacent diagonals.
On polynomials of degree at most `d`, its maximal Rayleigh quotient is

```math
\max_{\deg g\le d}\frac{\tau_{\rm sc}(t g^2)}{\tau_{\rm sc}(g^2)}
=2\cos\frac\pi{d+2}.
```

A maximizing polynomial is

```math
g_d(t)=\sqrt{\frac2{d+2}}
\sum_{j=0}^d\sin\!\left(\frac{(j+1)\pi}{d+2}\right)U_j(t/2).
```

For each fixed `d`, the preceding proof yields
`2 cos(pi/(d+2))/pi`. Only after those limits do we let `d` grow,
giving (T). For the fixed proportional-contraction falsifier, `d=3`
already suffices.

## 7. Literature and archive scope

The spectral analogue was proved by M. Magsino, D. G. Mixon, and
H. Parshall, *Kesten--McKay law for random subensembles*,
[author paper](https://par.nsf.gov/servlets/purl/10303784).
Its Theorem 1 treats Bernoulli principal submatrices of arbitrary
symmetric conference matrices at fixed retained proportion, along
lacunary order sequences. The argument here does not import an
eigenvector conclusion from that spectral theorem: (4.1) proves the
additional rooted diagonal concentration needed by Gaussian rounding.

Archive comparison found only the fixed-sample-then-growing-sample
obstruction in `mesoscopic_induced_sampling_no_go.md`, and an explicitly
unproved proportional extraction problem in
`exact_conference_proportional_extraction.md`. The new scope here is a
fixed positive proportion, arbitrary conference input family, and the
same-spin absolute quadratic objective. It remains an obstruction to
one convergence architecture, not a solution of the original problem.
