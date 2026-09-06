# Uniform sparse Hadamard compression: local moments and a cap above one half

Date: 2026-09-06. Seed-transfer track. The proof below is elementary and
uniform over the symmetric Hadamard parent. It covers simultaneous
`n->infinity`, `n/D->0`, not only an iterated fixed-retention limit.

## Theorem

Let `H_D` be ANY real symmetric Hadamard matrix of order `D`, with full
sign entries and `H_D^2=D I`. Choose a uniform subset `T` of exact size
`n<=D`, put `A=hollow((H_D)_T)`, and let `L=A/sqrt(n)`. Assume
`n->infinity` and `p=n/D->0`.

For every fixed nonnegative integer `k`, put
`m_(2j)=Catalan(j)` and `m_(2j+1)=0`. Then

```math
\mathbb E\frac1n\sum_i\big((L^k)_{ii}-m_k\big)^2
=O_k(p+1/n).                                               (1)
```

The constants depend only on `k`, not on `H_D,n,D`. In particular all
fixed global moments converge and the needed local polynomial
diagonals concentrate. The finite cubic Gaussian witness from
`transfer_seed_finite_spectral_gaussian_cap_2026_09_06.md` therefore gives

```math
\frac{Q(A)}{n^{3/2}}\ge\frac{576}{361\pi}-o_P(1)
>\frac{2016}{3971}-o_P(1)>\frac12-o_P(1).                    (2)
```

More precisely, the left side is at least any fixed number below
`576/(361pi)>0.50768` with probability tending to one. This concerns
UNIFORM random principal selectors. Carefully chosen exceptional
selectors and the original minima are not bounded below by (2).

## 1. A uniform graph contraction estimate

Allow finite undirected graphs with loops and count a loop as one edge
and twice toward its vertex degree. Parallel edges will always be
cancelled in pairs because Hadamard entries are signs. For a graph `F`
on `v` vertices, including any isolated vertices, define

```math
t_F(H)=D^{-v}\sum_{z_1,...,z_v}\prod_{uv\in E(F)}H_{z_u,z_v},
\qquad b(F)=v-c(F)-|E(F)|/2,
```

where `c(F)` counts components including isolated vertices.

Suppose every vertex degree is even. A vertex `w` of degree two with
two distinct neighbors `u,v` has no loop. Averaging its color gives

```math
D^{-1}\sum_z H_{z_u,z}H_{z,z_v}=1_{z_u=z_v}.
```

Identify `u,v`, remove `w`, and cancel all resulting parallel edges
in pairs, including pairs of loops. Denote the new graph by `F'`.
Then exactly

```math
t_F(H)=D^{-1}t_{F'}(H).                                    (3)
```

The graph retains all remaining isolated vertices. If `h` pairs of
edges were cancelled, then

```math
v'=v-2,\qquad e'=e-2-2h,\qquad
b(F)-b(F')=1+c(F')-c(F)-h\le1.                             (4)
```

Before cancellations the contraction does not increase the component
count. Removing each cancelled edge pair can increase it by at most
one, so `c(F')-c(F)<=h`; this proves the last inequality. Degree parity
remains even under all these operations.

Continue until no such vertex exists. A residual nontrivial component
other than a single loop has all degrees at least four, so its `b`
value is at most `-1`. A single-loop component has `b=-1/2`; an isolated
vertex contributes zero. Thus the final `b` is nonpositive.
If `q` eliminations were made, (4) gives `q>=b(F)`, while (3) gives

```math
|t_F(H)|\le D^{-q}\le D^{-\max(b(F),0)}.                    (5)
```

Moreover, if `F` is nonempty and `b(F)>=0`, at least one elimination
must be possible: a graph with no eliminable vertex has negative `b`
unless it is empty. In that case

```math
|t_F(H)|\le D^{-\max(b(F),1)}.                             (6)
```

The use of (6) when `b(F)=0` is important; (5) alone would not remove
all unwanted leading diagrams.

## 2. Closed-walk diagrams and the correct power balance

Let a connected closed-walk multigraph have `k` edge occurrences and
`v` vertices. Loops are allowed in this auxiliary statement. Reduce
each edge multiplicity modulo two to form `F`. It is Eulerian.
Put `a=v-1-k/2`.

Every connected component of `F`, including isolated vertices, must
be joined to the others by edges of positive even multiplicity. If
there are `c(F)` components, at least `c(F)-1` such edges are required.
Therefore

```math
k\ge |E(F)|+2(c(F)-1),\qquad a\le b(F).                    (7)
```

Combining (5) and (7) shows, uniformly for `1<=n<=D`,

```math
|n^a t_F(H)|\le1.                                         (8)
```

Indeed `a<=0` is immediate, while `a>0` gives
`n^a|t_F|<=(n/D)^a<=1`.

There is a sharper classification for a fixed original walk:

- If `F` is empty and `a=0`, the underlying graph is a tree and every
  edge is traversed exactly twice. There are no loops.
- If `a<0`, its contribution is at most `n^a`.
- If `F` is nonempty and `a=0`, (6) bounds it by `1/D`.
- If `a>0`, its contribution is at most `p^a`.

For even `k`, nonleading terms are consequently `O_k(p+1/n)`.
For odd `k`, there is no doubled-tree term and the bound is
`O_k(sqrt(p)+1/sqrt(n))`. There are only finitely many abstract walk
patterns at every fixed length.

## 3. Exact-size sampling and why no collision-probability approximation is used

An abstract walk with `v` distinct SAMPLE vertices assigns them to
`v` distinct ambient Hadamard colors. Let `t_F^inj(H)` be its average
over these injective ambient assignments. Partition-lattice
inclusion-exclusion gives exactly

```math
t_F^{\rm inj}(H)
=\frac{D^v}{(D)_v}
  \sum_{\pi}\mu(\pi)D^{-c_\pi}t_{F/\pi}(H),
\quad c_\pi=v-|\pi|,
\quad \mu(\pi)=\prod_{B\in\pi}(-1)^{|B|-1}(|B|-1)! .         (9)
```

The quotient may create loops, which is why Section 1 allowed them.
It remains the parity graph of the quotient closed walk of the SAME
length `k` but with `v-c_pi` vertices. Its exponent in (8) is
`a-c_pi`. Hence each term with a nontrivial partition satisfies

```math
|n^a D^{-c_\pi}t_{F/\pi}(H)|
=p^{c_\pi}|n^{a-c_\pi}t_{F/\pi}(H)|
\le p^{c_\pi}\le p.                                      (10)
```

The fixed-length prefactor `D^v/(D)_v` is `1+O_k(1/D)`.
Also `(n)_v/n^v=1+O_k(1/n)` when counting the sample labels.
Thus all injectivity corrections are `O_k(p+1/n)`. In particular,
the proof does NOT require `n^2/D->0`; the total probability of a
collision in an independently drawn entire sample may be large.

The doubled-tree walks yield the Catalan moments by the usual rooted
plane-tree/contour-walk bijection: a length-`2j` walk which traverses
each tree edge exactly twice is the contour of a rooted plane tree
with `j` edges. There are `Catalan(j)` such abstract patterns.
Sections 2--3 therefore prove the stated global moment expectations.

## 4. Rooted double walks give the needed local diagonal concentration

For `(L^k)_(ii)(L^ell)_(ii)`, expand TWO closed walks with a common
sample root. Their union is connected and their combined length is
`k+ell`, so the entire preceding diagram estimate applies unchanged.
The sample-label count is now `(n-1)_(v-1)` rather than `(n)_v/n`;
these are identical.

A leading doubled-tree union has an additional rigidity. Any closed
walk on a tree uses every edge an even number of times. Since the two
walks jointly use each edge exactly twice, every edge belongs to
exactly one of the walks. Their connected edge sets can meet only at
the common root: two connected subtrees of a tree sharing another
vertex would share the intervening path as well.

Consequently the leading patterns are exactly independent pairs of
rooted plane trees when both lengths are even. Their number is
`m_k m_ell`. If either length is odd, there is no such leading pattern.
When `k+ell` is even,

```math
\mathbb E (L^k)_{ii}(L^\ell)_{ii}
=m_km_\ell+O_{k,\ell}(p+1/n).                              (11)
```

For even `k`, combine the single-walk estimate with (11) at `ell=k`.
For odd `k`, `m_k=0` and (11) already bounds its square. This proves
(1), after averaging over the exchangeable sample root.

## 5. The actual Boolean-cap consequence

Let `P(x)=x^3+(8/5)x^2-(2/5)x-3/5` and
`R=I/10+P(L)^2`. Equation (1), by finite linear combinations and
Cauchy--Schwarz, gives convergence of the global traces through degree
7 and mean-square concentration of the degree-6 diagonal of `P(L)^2`.
Global moments through degree 12 are bounded (indeed converge), so
`tr(R^2)/n=O_P(1)`.

These are exactly the hypotheses of the previously proved finite
Gaussian covariance inequality. Its error terms vanish in probability,
and its deterministic moment values give (2).

This is a random-selector obstruction, not a universal lower bound on
all principal restrictions. It also does not assert the exact limiting
cap, prove an SK limit, or treat a bounded-operator non-Hadamard parent.
The identity `H^2=D I` was used exactly in (3), and cannot silently be
replaced by an operator inequality.

Exact replay source:
`computations/transfer_seed_hadamard_graph_verify_2026_09_06.py`.
It checks all 2198 Eulerian loop graphs on at most five labeled vertices
for the elimination/component inequalities, and checks 150 exact
homomorphism identities and 150 exact partition-Mobius identities on
the order-four Sylvester Hadamard.
