# Exact exponential pressure for nearly homogeneous complete-graph row profiles

Date: 2026-09-06. The theorem below is a combinatorial counting result.
Its lower bound is for the actual edge contraction, not the product-norm
upper bound supplied by Finner. No second-moment or typical-instance
conclusion follows from it alone.

## 1. Statement, including heterogeneous rows

Fix a finite alphabet `A` of size `q`, a probability vector `nu` on it,
and a symmetric kernel `K:A times A -> (0,infinity)`. The alphabet and
kernel do not depend on the graph order. For each `m`, let
`d=m-1`, `E=m(m-1)/2`, and prescribe integer endpoint-color counts
`c_(i,a)>=0`, `sum_a c_(i,a)=d`, at each vertex `i` of `K_m`. Assume

```math
\sum_{i=1}^m\sum_{a\in A}|c_{i,a}-d\nu_a|=o(m^2).
\tag{1}
```

Independently at each vertex, uniformly assign its multiset of colors
to its incident edges. Write `a_ij` for the color at the endpoint `i`
of edge `{i,j}`, and define

```math
Z_m=\mathbb E\prod_{i<j}K(a_{ij},a_{ji}).
```

Then

```math
\boxed{\quad
\lim_{m\to\infty}\frac1{m^2}\log Z_m
=\frac12\sup_{\substack{\gamma=\gamma^T\ge0\\\gamma\mathbf1=\nu}}
\left\{H(\gamma)-2H(\nu)+\sum_{a,b}\gamma_{ab}\log K(a,b)\right\}.
\quad}
\tag{2}
```

The convention is `0 log0=0`. In particular the theorem permits zero
cells of `gamma`; it even permits zero coordinates of `nu`. For positive
`nu` and `K` the entropy optimizer is in fact strictly positive, but this
extra fact is not needed.

The homogeneous case has `c_(i,a)=c_a` and `c_a/d ->nu_a`. The requested
extension with every row within `o(m)` of a common count vector follows
from (1). Condition (1) also allows an `o(m)` exceptional set of arbitrary
rows. Positivity of `K`, rather than positive semidefiniteness, is the
relevant hypothesis for the lower bound.

Equivalently the right side is `-F_K(nu)/2`, where

```math
F_K(\nu)=\inf_{\gamma\in\Pi(\nu,\nu)}
\left\{D(\gamma\Vert\nu\otimes\nu)-\sum_{a,b}\gamma_{ab}\log K(a,b)\right\}.
```

Transposition averaging shows that the unrestricted coupling infimum has
a symmetric optimizer.

## 2. Normalization and the upper bound

An endpoint configuration is an array `(a_ij)_(i!=j)`. The number meeting
all prescribed row counts is exactly

```math
N_m=\prod_{i=1}^m\frac{d!}{\prod_a c_{i,a}!}.
\tag{3}
```

Every such configuration has probability `1/N_m`. Put
`nu_i=c_i/d` and `bar nu_m=m^-1 sum_i nu_i`. Uniform Stirling bounds give

```math
\log N_m=d\sum_iH(\nu_i)+O_q(m\log(m+1))
=m^2H(\nu)+o(m^2).
\tag{4}
```

For the last equality, (1) implies that the average `l^1` distance of
`nu_i` from `nu` tends to zero. Uniform continuity of entropy on the
finite probability simplex then gives average entropy convergence.

Orient edges by `i<j`. If `r_ab` is the empirical law of their ordered
color pairs, its averaged first and second marginals equal `bar nu_m`.
There are at most `(E+1)^(q^2)` possible tables, and the number of arbitrary
edge assignments with a given table is at most `exp(E H(r))`. Ignoring
the individual row constraints is an upper bound. Symmetry of `K` gives

```math
H(r)+\sum r_{ab}\log K(a,b)
\le H((r+r^T)/2)+\sum\frac{r_{ab}+r_{ba}}2\log K(a,b),
```

and the symmetrized table has marginal `bar nu_m`. Divide the resulting
weighted count by (3). Compactness of the finite table simplex and
continuity of entropy show that the limsup is at most (2): any sequence
of maximizing tables has a subsequential limit with marginal `nu`.
This proves the upper bound without any PSD or graph-contraction theorem.

## 3. A count-preserving endpoint repair

Fix any symmetric coupling `gamma` with marginal `nu`. Independently on
each of the `E` oriented edges, sample its ordered endpoint pair from
`gamma`. Let `b_(i,a)` be the resulting row counts. At a fixed vertex,
the `d` incident endpoint colors are independent with law `nu`, because
the edge variables are independent and both marginals of `gamma` equal
`nu`. Therefore

```math
\mathbb E\sum_{i,a}|b_{i,a}-d\nu_a|\le q m\sqrt d.
```

Markov gives, with probability tending to one,

```math
\sum_{i,a}|b_{i,a}-d\nu_a|\le q m^{7/4}.
\tag{5}
```

For each row, change one endpoint of a surplus color to a deficit color
until its exact target counts are reached. Fix an ordering of vertices,
colors and neighbors so that this is a deterministic map. It changes
exactly `sum_a |b_(i,a)-c_(i,a)|/2` entries of that row and does not change
any other row's counts: the two endpoints of an edge are separate color
variables. The total number of changes is at most

```math
R_m=\frac12\left[qm^{7/4}+\sum_{i,a}|c_{i,a}-d\nu_a|\right]+1=o(m^2).
\tag{6}
```

No graph-degree realization, graphicality, parity, or exact pair-table
rounding theorem is required. The output pair histogram need only remain
within `o(1)` of `gamma`, not equal it exactly.

## 4. Typical source count, preimages, and weights

The empirical edge-pair table converges to `gamma`. Restrict to samples
for which its `l^1` error is at most `E^-1/4` and (5) holds; this event
has probability tending to one. Only cells with `gamma_ab>0` occur in
these raw samples. Since the positive support is a fixed finite set,

```math
\log\mathbb P\{\text{a particular raw sample}\}
=-E H(\gamma)+o(m^2),
\qquad
\log\prod_{i<j}K(a_{ij},a_{ji})
=E\sum\gamma_{ab}\log K(a,b)+o(m^2),
\tag{7}
```

uniformly on this event. Its number of raw configurations is consequently
at least `exp(E H(gamma)-o(m^2))`.

For any repaired output, every preimage is within `R_m` changed positions
among `2E` labeled endpoints. If `q>=2`, its preimage count is at most

```math
\sum_{r\le R_m}\binom{2E}{r}(q-1)^r
\le\exp(o(m^2)).
\tag{8}
```

The usual binomial-entropy bound proves the last statement since
`R_m/(2E)->0`. For `q=1` it is immediate. Hence the set of repaired
outputs has size at least `exp(E H(gamma)-o(m^2))`.

Let `osc(log K)=max log K-min log K`, which is finite by strict
positivity. At most `R_m` edges are touched by the repair, and therefore
its logarithmic weight decreases by at most `R_m osc(log K)=o(m^2)`.
This remains valid if a repaired pair lies in a cell where `gamma_ab=0`:
only `K`, not `gamma`, weights the repaired output. Zero cells therefore
cause no singular likelihood comparison.

Combining this weight bound with (7)--(8), the exact-count numerator is
at least

```math
\exp\left\{E\left[H(\gamma)+\sum\gamma_{ab}\log K(a,b)\right]-o(m^2)\right\}.
```

Divide by (3), use (4), and optimize `gamma`. The coupling simplex is
compact and its entropy-plus-energy objective is continuous, so an
optimizer exists. This proves the lower bound and completes (2).

## 5. Relation to Finner and exact scope

For a PSD kernel, the Finner/Gram contraction upper bound is the product
of the square-root self-permanents of the row multisets. Their finite-type
exponents are `-F_K(nu_i)/2`; under (1) their total exponent is exactly
the right side of (2). Thus Finner is exponentially sharp on nearly
homogeneous finite row profiles on `K_m`, although it need not be an
equality at finite order.

The lower bound applies to the actual positive edge expectation. It says
that merely retaining the joint contraction, while leaving all typical
rows asymptotically homogeneous and retaining the same independent edge
law, cannot produce an exponentially smaller value than its Finner rate.
It does not assert concentration of the partition sum or existence of
high-energy spins in typical realizations.

The assumptions that the alphabet and positive kernel are fixed matter.
For a growing alphabet or a vanishing minimum kernel entry, (8) and the
weight-repair payment need new uniform bounds. Such cases can instead be
handled by first fixing a finite truncation/quantization, taking `m` to
infinity, and only afterward refining it, provided the application proves
the required approximation inequalities separately.

## 6. A signed hard-tail counting version

In addition to the endpoint assignments above, independently give each
edge a fair sign `s_ij`. Let `psi(a,b,s)` be any fixed bounded score
symmetric in `a,b`. Let `theta(a,b,s)` be a symmetric probability table
whose two endpoint marginals are `nu`. For every fixed
`u<sum theta psi`,

```math
\liminf_{m\to\infty}\frac1{m^2}
\log\mathbb P\left\{\frac1E\sum_{i<j}
\psi(a_{ij},a_{ji},s_{ij})\ge u\right\}
\ge-\frac12D(\theta\Vert\nu\otimes\nu\otimes\mathrm{fair}).
\tag{9}
```

To prove this, sample independent edge triples from `theta` and apply
the same endpoint repair, leaving every edge sign unchanged. The typical
raw triple count is `exp(E H(theta)-o(m^2))`; the preimage bound remains
(8), since no signs are edited. The score changes by at most
`R_m osc(psi)=o(m^2)`, so the strict threshold is eventually met. Under
the actual uniform-row/fair-sign model, every repaired triple array has
probability `1/(N_m 2^E)`. Dividing its count by this denominator gives
(9), since both endpoint marginals of `theta` are `nu`. This proof
retains the fair-sign entropy cost `E log2` explicitly.

## 7. Independent exact finite checks

Run

```sh
.venv/bin/python computations/transfer_reconstruction_homogeneous_edge_count_exact_checks_2026_09_06.py
```

The checker computes the exact weighted endpoint partition by edgewise
dynamic programming for orders two through seven, for both the PSD kernel
`[[2,1],[1,2]]` and the positive non-PSD kernel `[[1,2],[2,1]]`. A
constant-kernel count independently verifies the product of row
multinomials. It checks the PSD Finner inequalities in exact rational
arithmetic, and all 64 endpoint configurations of a three-vertex repair,
including exact counts, Hamming changes, preimages and logarithmic-weight
loss after exponentiation. All checks passed. The limiting pressure of
both sample kernels is `(1/2)log(3/2)` by (2), not inferred from these
finite values.
