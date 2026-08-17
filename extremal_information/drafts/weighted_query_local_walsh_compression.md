# Weighted query-local Walsh compression

Status: task-local theorem draft.  It combines the exact unrooted Walsh
orbit collapse with a scale-sensitive deletion principle.  The result is a
response upper theorem, not a lower bound for every possible carrier.

## 1. Setup

Let `q=2^m`, `n=q^2`, and let `W` be the order-`n` Walsh bridge.  A marked
linear-label Walsh graph has one Boolean block `x_v in {+-1}^n` at every
active vertex, child energy `H_a(x)=x^TC_ax/2`, arbitrary real onsite
weights `h_v`, and weighted bridge energy

```math
E_(G,a)(x)=sum_v h_vH_(a_v)(x_v)
 +sum_(uv in E(G)) w_(uv)x_u^TWx_v.                    \tag{WQ.1}
```

The Walsh normalization gives

```math
|x^TWy|<=n^(3/2)                                       \tag{WQ.2}
```

for all Boolean `x,y`, and equality is attainable.  For a vertex set `C`,
write

```math
I_a(C)=(G(a|_C),R(a|_C))                               \tag{WQ.3}
```

for its unrooted Gram/relation state.  The ambient-Witt theorem says that
`I_a(C)` determines the complete Boolean landscape on every weighted graph
supported on `C`, up to a common coordinate permutation.

## 2. A mass-truncated response carrier

For a partition `P` of the vertices of `G`, let

```math
d_G(P)=sum_(uv in E(G):[u]_P ne [v]_P)|w_(uv)|         \tag{WQ.4}
```

be its deleted interaction mass, and let `G[P]` be the graph obtained by
deleting all cross-part edges.

### Theorem WQ.1 (interaction-mass/local-memory law)

For every public partition `P`, fixed from the graph query and error budget
rather than from hidden labels, the decoder that stores the states

```math
(I_a(C))_(C in P)                                      \tag{WQ.5}
```

and returns the sum of the exact component upper optima on `G[P]`
approximates the upper optimum in (WQ.1) with error at most

```math
d_G(P)n^(3/2).                                         \tag{WQ.6}
```

For a finite declared family `Theta` of weighted graph queries, one may
publish a different partition `P_theta` for each query and store the local
states on all inclusion-maximal parts that occur.  The same code answers
each `theta` with its own error
`d_(G_theta)(P_theta)n^(3/2)`.

If

```math
L(P_Theta)=sum_(C maximal among all P_theta)|C|^2,      \tag{WQ.7}
```

then a direct binary presentation uses `O(L(P_Theta))` bits.  Thus the
scale-sensitive response complexity is bounded by the graph functional

```math
C_Theta(delta)=
 min\left\{L(P_Theta):
       d_(G_theta)(P_theta)<=delta_theta
       \text{ for every }theta\right\}.                \tag{WQ.8}
```

The statement remains valid for independent rooted component queries after
replacing `I_a(C)` by the rooted relation-form state on `C`.  A shared
auxiliary variable joins the corresponding supports and cannot be decoded
componentwise.

#### Proof

Deleting one edge changes the energy of every Boolean configuration by at
most `|w_e|n^(3/2)` by (WQ.2).  Summing before taking a maximum proves

```math
|max E_(G,a)-max E_(G[P],a)|<=d_G(P)n^(3/2).            \tag{WQ.9}
```

The truncated graph is a disjoint union.  Its upper optimum is the sum of the
component upper optima, and each component landscape is determined by (WQ.3).
This proves the decoder and (WQ.6).  Taking the union of the local states
needed by all declared queries proves the simultaneous claim.  If `C` is
contained in a stored maximal part `C'`, its Gram form is a principal
restriction and its relation kernel is `R(C') cap F_2^C`.  A Gram
matrix and a row-reduced basis of a relation kernel on `k` labels use
`O(k^2)` bits, proving (WQ.7)--(WQ.8).  The rooted version uses the same
argument and one additional affine root fibre per part. `square`

For the absolute objective, the same pointwise perturbation bound holds, but
the exact truncated decoder is

```math
max\left\{sum_C max E_C,-sum_C min E_C\right\},         \tag{WQ.9a}
```

not the sum of component absolute maxima.  The local orbit state preserves
the entire component landscape and therefore supplies both extrema.

The error in (WQ.6) is a one-copy decoder error.  Consequently two systems
with the same carrier (WQ.5) are at response distance at most
`2d_G(P)n^(3/2)` for this query.  No factor two is needed when estimating
either system from its decoded value.

## 3. Regime separation

The graph functional in (WQ.8) distinguishes three elementary regimes
without appealing to algebraic matrix rank.

### Corollary WQ.2 (bounded clusters, paths, and dense graphs)

1. If every query is already a union of components of size at most `w` and
   the total vertex incidence is `L`, then `delta=0` and the exact carrier
   costs `O(wL)` bits.

2. For the unit-weight path on `t` vertices, allowing error
   `eta t n^(3/2)` permits components of size at most
   `ceil(1/eta)+1` and hence a carrier of `O(t/eta)` bits for fixed
   `0<eta<1`.

3. If `G=K_t` and every edge has absolute weight at least `c>0`, every
   partition with deleted mass at most `eta t^2` satisfies

```math
sum_(C in P)|C|^2
 >=\left(1-{2eta\over c}\right)t^2.                    \tag{WQ.10}
```

   Hence this component-deletion architecture remains quadratic whenever
   `eta<c/2` is fixed.

#### Proof

The first claim is immediate.  For the path, cut at most `eta t` edges at
approximately equal spacings `1/eta`; the sum of squared component sizes is
`O(t/eta)`.  For the dense claim, the number of cross-part edges is

```math
{1\over2}\left(t^2-sum_C|C|^2\right),                 \tag{WQ.11}
```

exactly.  Multiplication by `c` and the deleted-mass assumption give
(WQ.10). `square`

The dense conclusion is only a lower bound on this *specific upper
architecture*.  It is not an information lower bound on arbitrary semantic
codes.  Conversely the path result is approximate query-local compression,
not a claim that treewidth alone gives a bounded exact Walsh state: nonzero
edges synchronize coordinate frames along the whole connected path, and
global label relations can remain response-visible.

## 4. Interpretation

The theorem gives a quantitative form of

```text
reusable compatibility memory
    <= local orbit information remaining after subscale interactions
       are deleted.
```

It is not ordinary graph sparsification.  The retained object is an exact
algebraic response quotient on each surviving component, and the discarded
quantity is the total interaction mass measured directly at the Boolean
extremal scale.  The two inputs play different roles: the orbit theorem
compresses strong local interactions, while (WQ.2) certifies that weak
cross-component interactions cannot restore the omitted compatibility at
more than the declared error.

The next lower question is whether a dense family can have small semantic
response complexity despite the quadratic value of (WQ.8), or whether a
positive fraction of its Gram/relation information is exposed by ordinary
unrooted graph maxima.

The exact verifier checks the pointwise deletion inequality for signed
four-block queries and all set partitions at `m=1`, together with the dense
cross-edge identity.  It is a regression test for the Lipschitz step; the
ambient orbit theorem supplies the semantic component equivalence.
