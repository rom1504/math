# Exact positive edge-transport limit for finitely many vertex profiles

Date: 2026-09-06. Status: proved finite-alphabet theorem; application to the
full signing construction is a separate obligation. This does not claim a
limit for M_n. It avoids the leading complex-frame cancellation loss by
working with nonnegative physical edge weights throughout.

## Theorem

Fix finite vertex classes a=1,...,r, finite colors c=1,...,D, a strictly
positive symmetric matrix K(c,d), and class proportions pi_a>0 summing to one.
At each vertex of class a, independently choose a uniform arrangement of
n-1 outgoing colors with prescribed counts k_(a,c)^(n), where
k_(a,c)^(n)/(n-1)->nu_(a,c). Let n_a/n->pi_a. Define

```math
Z_n=\mathbb E\prod_{i<j}K(X_{ij},X_{ji}).
```

Then

```math
\lim_n n^{-2}\log Z_n
=\frac12\max_\gamma\sum_{a,b}\pi_a\pi_b
 \left[H(\gamma_{ab})+\sum_{c,d}\gamma_{ab}(c,d)\log K(c,d)\right]
-\sum_a\pi_a H(\nu_a).                                      \tag{1}
```

The maximum is over probability matrices gamma_ab with

```math
\gamma_{ba}(d,c)=\gamma_{ab}(c,d),\qquad
\sum_b\pi_b\sum_d\gamma_{ab}(c,d)=\nu_{a,c}.                 \tag{2}
```

No constraint says that the marginal toward EACH other class is nu_a:
only its average over neighboring classes is prescribed. Omitting that
distinction would incorrectly forbid allocation of colors by future neighbor.

## Upper bound

Write Z as the weighted sum over the legal directed color arrays, divided
by the product of their row multinomial counts. The latter has logarithm
n^2 sum_a pi_a H(nu_a)+o(n^2), by Stirling, with fixed D,r.

For any probability distribution P on legal arrays, entropy subadditivity
over UNDIRECTED edges gives H(P)<=sum_(i<j) H(P_(Xij,Xji)). Gibbs' finite
variational principle bounds its energy-plus-entropy by the same sum of
edge entropy plus expected log K. Average those pair laws within each
ordered class pair. Entropy concavity increases the bound, and exact row
counts impose (2) up to o(1); edges inside a class omit only n_a diagonal
slots. Compactness of the finite collection of probability simplices and
continuity of entropy give the upper bound in (1).

## Lower bound: repair costs only o(n^2)

Fix any admissible gamma. Independently for every undirected edge choose
its pair of colors with distribution gamma_ab. The expected outgoing count
at a vertex of class a differs from its prescribed count by o(n), uniformly
over vertices: proportions and prescribed profiles converge, and deleting
the self slot contributes at most one. Hoeffding and a union bound show
that, with probability tending to one, every outgoing count differs by at
most o(n)+O(sqrt(n log n)).

For each such array repair each ROW separately: replace surplus colors by
deficit colors. An endpoint color can be changed without changing its
opposite endpoint, so this is always feasible. The total number of changed
incidences is s_n=o(n^2). There are at most

```math
\sum_{j\le s_n}\binom{n(n-1)}j D^j=\exp(o(n^2))
```

arrays mapping to any given repaired array, for a deterministic repair rule.
Since K has positive fixed minimum and finite maximum, the repair changes
the log edge weight by at most s_n times the oscillation of log K, or o(n^2).

For completeness, under the independent-edge law the information density
-log P(array) and log weight concentrate within o(n^2) of their expectations.
All nonzero gamma entries are fixed positive numbers; zero entries are never
drawn. Bounded independent summands and Chebyshev (or Hoeffding) suffice.
Intersect this event with the count event. Its probability tends to one,
so it contains at least exp(sum_edges H(gamma_ab)-o(n^2)) arrays, each with
weight at least exp(sum_edges <gamma_ab,log K>-o(n^2)). Repair, divide by
the subexponential preimage bound, and then divide by the exact row-count
denominator. This gives the lower bound for gamma. Optimize to obtain (1).

## Homogeneous case and normalization

For r=1 the constraint is gamma=gamma^T and gamma 1=nu. Formula (1) becomes

```math
\frac12\max_{\gamma=\gamma^T,\,\gamma1=\nu}
\left[\langle\gamma,\log K\rangle-D(\gamma\Vert\nu\otimes\nu)\right].
```

Thus positive permutation interaction has an exact self-transport limit,
not merely a Hilbert-norm upper bound. Heterogeneous profiles require (2),
not separately fixed pairwise marginals. The proof uses ordinary finite
entropy and row repair; its usefulness is the exact nonlocal compatibility
constraint and its uniform all-order realization, not a new name for entropy.

## Application boundary

The grouped weave has positive physical kernel
exp[-t(||x||^2+||y||^2)] cosh(2t x^T R_A y). The companion bridge proof
uses a deterministic energy cutoff, then finite physical quantization, to
retain R_A with arbitrarily small pressure error. Applying (1) also requires
tracking the distribution/count of the joint vertex profiles after summing
over spins and randomized bases. This theorem alone does not estimate that
profile rate, does not identify the maximum with the scalar Bellman envelope,
and does not propagate an arbitrary seed's Boolean cap.
