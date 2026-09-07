# Symmetric type pressure is seed-universal

Date: 2026-09-07. Status: **Pending independent audit**. This note contains a
complete proof attempt with an explicit subleading error. It concerns an
annealed kernel partition, not the maximum energy of an actual parent.

## 1. Statement

Let `E` be a fixed finite alphabet equipped with a fixed-point-allowed
involution `a -> -a`. Let `nu` be a positive probability law on `E` satisfying
`nu(a)=nu(-a)`. Fix a strictly positive kernel `K:E x E -> (0,infinity)` with

```math
K(a,b)=K(b,a)=K(-a,-b).
```

For each integer `m` for which `d nu(a)` is integral, where `d=m-1`, let
`u_i=(u_ij)_{j != i}` independently be a uniform permutation of the type
having exactly `d nu(a)` copies of every `a`. Let `S` be ANY symmetric
off-diagonal signing of order `m`. Define

```math
Z_S=\mathbb E\prod_{i<j}K(u_{ij},S_{ij}u_{ji}).
```

Write `e=m(m-1)/2` and

```math
\mathcal V_K(\nu)=
\max_{\pi\in\Pi(\nu,\nu)}
\left\{\sum_{a,b}\pi(a,b)\log K(a,b)
-D(\pi\Vert\nu\otimes\nu)\right\}.
```

Then, with constants depending only on `(E,nu,K)`, uniformly over EVERY `S`,

```math
\log Z_S=e\,\mathcal V_K(\nu)+O(m^{3/2}\log m).
\tag{1}
```

In particular the difference between the log partitions for any two seeds
is `O(m^(3/2) log m)=o(m^2)`. No cap, spectral, pseudorandomness, or optimality
condition on the seed is assumed.

The Gaussian application is `K(a,b)=exp[-t(a-b)^2]`, for which the hypotheses
hold at every fixed `t>0` and fixed finite alphabet. The theorem does not
require explicit independent sign phases at the output columns.

## 2. Symmetric Sinkhorn optimizer

The variational objective is strictly concave in `pi` and the transport
polytope is compact. Its unique maximizer `pi_*` is strictly positive: if a
zero cell were present, mixing slightly toward `nu tensor nu` gains a
positive `epsilon log(1/epsilon)` entropy term, dominating the bounded
linear and other entropy changes.

Lagrange multipliers therefore give positive scalings `r(a),c(b)` such that

```math
\pi_*(a,b)=\nu(a)\nu(b)K(a,b)r(a)c(b).
```

Both transpose and simultaneous involution preserve feasibility and the
objective. Uniqueness makes `pi_*` invariant under both. Symmetry permits a
common scaling: the identity `r(a)c(b)=r(b)c(a)` makes `r/c` constant, so
absorb its square root and write

```math
\pi_*(a,b)=\nu(a)\nu(b)K(a,b)q(a)q(b).
\tag{2}
```

Since all entries are positive, involution symmetry gives
`q(-a)q(-b)=q(a)q(b)` for all a,b; taking a=b yields `q(-a)=q(a)`.
Substituting (2) into the objective gives

```math
\mathcal V_K(\nu)=-2\sum_a\nu(a)\log q(a).
\tag{3}
```

## 3. Exact partition identity

Let `P_0` be the law under which every directed incidence `u_ij` is
independent with law `nu`. Let `T` be the event that all m row types are
exactly `d nu`. Conditioning `P_0` on `T` produces precisely the independent
uniform row permutations in the theorem.

Let `P_S` be the edge-independent law whose pair at edge i<j is

```math
P_S(u_{ij}=a,u_{ji}=b)=\pi_*(a,S_{ij}b).
```

Every directed incidence has marginal nu. Equations (2) and the evenness of
q yield the exact density ratio

```math
\frac{dP_S}{dP_0}(u)
=\prod_{i<j}K(u_{ij},S_{ij}u_{ji})
 \prod_i\prod_{j\ne i}q(u_{ij}).
```

On the event T, the second product is the constant

```math
\exp\!\left(md\sum_a\nu(a)\log q(a)\right)
=\exp[-e\,\mathcal V_K(\nu)].
```

Consequently

```math
Z_S=\exp[e\,\mathcal V_K(\nu)]\frac{P_S(T)}{P_0(T)}.
\tag{4}
```

This identity isolates every dependence on the signing S in the exact-type
probability `P_S(T)`. In particular the leading variational exponent is
already seed-independent, provided exact types are not exponentially
expensive at speed m².

## 4. Uniform exact-type repair

Let `L=|E|`. Under `P_S`, the incidences in any ONE row lie on different
edges, so they are independent with common law nu. Different rows need not
be independent. Let `N_i(a)` denote the row counts. The minimum number of
directed entries requiring replacement to repair row i is

```math
r_i=\tfrac12\sum_a|N_i(a)-d\nu(a)|.
```

Cauchy--Schwarz, followed by the binomial variance bound, gives

```math
\mathbb E\sum_i r_i
\le \tfrac m2\sum_a\sqrt{d\nu(a)(1-\nu(a))}
\le \tfrac L2m\sqrt d.
```

Hence the set G of arrays with `sum_i r_i <= L m sqrt(d)` has `P_S(G)>=1/2`.
Use a deterministic rowwise rule to replace surplus symbols by deficient
symbols. This maps G into T and alters at most
`R=ceil(L m sqrt(d))` directed incidences, thus at most R undirected edges.

Let

```math
\rho=\frac{\min_{a,b}\pi_*(a,b)}{\max_{a,b}\pi_*(a,b)}>0.
```

If an array u maps to v, edge independence gives
`P_S(v)>=rho^R P_S(u)`. For each fixed v the number of possible preimages is
at most

```math
B_R=\sum_{r=0}^{R}\binom{md}{r}(L-1)^r
\le\exp[O_{L}(m^{3/2}\log m)].
```

(For the finitely many orders with R>=md, use the trivial full-space count;
the asserted asymptotic bound concerns sufficiently large m.) Summing over
G and its image proves

```math
P_S(T)\ge\frac{\rho^R}{2B_R}
\ge\exp[-O_{E,\nu,K}(m^{3/2}\log m)],
\tag{5}
```

uniformly over S. Also `P_S(T)<=1`. The multinomial type probability for one
row under `P_0` is at least `(d+1)^(-L)`; multiplying over independent rows
gives `-log P_0(T)=O_L(m log m)`. Inserting these two estimates into (4)
proves (1).

## 5. Scope and the remaining escape routes

This proves that signed cycle holonomies can affect finite partitions but
cannot change their leading m² rate for fixed symmetric row types and
strictly positive fixed kernels. It does NOT show that every actual weave
has the same maximum; converting a first-moment partition to an actual cap
can be loose in a seed-dependent way.

The proof does not cover:

- asymmetric row types with an order-one signed imbalance;
- kernels whose minimum tends to zero too quickly with m (the repair cost
  explicitly depends on `log(1/rho)`);
- alphabets or recursive depths growing so quickly that the repair error is
  leading order;
- retained cross-row dependencies beyond independent row permutations; or
- a quenched optimization over row bases/signings that changes the relation
  between partition size and the extremal parent.

Thus the positive reusable principle is: **a dense edge-product model with
finite symmetric local types admits a uniformly subextensive exact-type
compiler**. The negative implication applies only to extracting a leading
seed-frustration gain from this particular annealed kernel pressure.
