# Exact homogeneous positive permutation pressure

Date: 2026-09-06. Status: proved finite-alphabet asymptotic theorem, with
uniform bounded-domain quantization. This is an ensemble statement with
FIXED outgoing empirical profiles, not an identification of the original
minimum cap or of its optimized recursive ensemble pressure.

## Theorem

Fix a finite alphabet A of size s and a symmetric matrix K_ab>0. For each
n let c_a be nonnegative integers summing to d=n-1, with p_n=c/d -> p.
At every vertex i independently assign its n-1 outgoing directed slots a
uniform random permutation of a multiset containing c_a copies of a.
Let X_ij be the color at endpoint i of edge ij. Then

`lim n^-2 log E product_(i<j) K_(X_ij,X_ji)
 = (1/2) max_(gamma=gamma^T, gamma 1=p)
       [sum_ab gamma_ab log K_ab - D(gamma || p tensor p)]`.

Zero coordinates of p can be removed by o(n²) directed-slot repairs; thus
it suffices to prove the theorem with all p_a>0. The maximum exists by
compactness. The functional is strictly concave on feasible distributions,
so the maximizer is unique (up to coordinates forced zero).

## Upper bound, including the finite-n correction

Let Omega be the set of exact-count directed configurations. Its cardinality
is `multinomial(d;c)^n`. Apply the finite Gibbs variational formula to
`sum_(omega in Omega) product_edges K`. For any law rho on Omega, let
gamma_ij be the joint endpoint-color distribution on each undirected edge,
with the endpoint order chosen i<j. Symmetrize their average to gamma_bar.
The exact row counts imply gamma_bar 1=p_n. Entropy subadditivity over
the edge-pair variables and concavity of Shannon entropy give

`log Z_n <= binom(n,2) max_(gamma symmetric, gamma1=p_n)
                   [H(gamma)+sum gamma logK]
               - n log multinomial(n-1;c)`.

Indeed the expected edge energy depends only on the symmetrized average
because K is symmetric. Stirling with fixed alphabet gives
`log multinomial(d;c)=d H(p_n)+O_s(log(n+1))`, including zero counts.
This yields the claimed limsup. No independence of edges under rho is used.

## Lower bound by repairing independently colored edges

First fix a strictly positive symmetric gamma_n with marginal p_n, tending
to a strictly positive feasible gamma. Such approximations suffice because
mixing any target gamma with p tensor p makes it positive at arbitrarily
small objective cost. Draw edge-pair colors independently with law gamma_n;
orientation is harmless by symmetry. With probability tending to one:

1. every vertex/color count differs from its target c_a by
   O_s(sqrt(n log n));
2. every global edge-pair frequency differs from gamma_n by o(1).

These follow directly from Hoeffding and a union bound. Since gamma_n has
uniformly positive entries, each configuration in this typical set has
probability `exp[-binom(n,2)H(gamma_n)+o(n²)]`. Hence the number of typical
configurations is at least `exp[binom(n,2)H(gamma_n)-o(n²)]`.

Repair each vertex's outgoing endpoint colors independently: replace a
surplus color by a deficit color until its counts are exactly c. This alters
at most r=O_s(n^(3/2)sqrt(log n)) directed endpoint colors in total. The
other endpoint's count is not changed, so no iterative graph-degree repair
is necessary. Make the map deterministic, for example by lexicographic
choice of slots and colors. Each repaired image has at most

`sum_(j<=r) binom(n(n-1),j) s^j = exp(o(n²))`

preimages. Consequently the image still contains
`exp[binom(n,2)H(gamma_n)-o(n²)]` configurations. Bounded logK implies each
has total log weight `binom(n,2)sum gamma_n logK-o(n²)`. Divide their total
weight by `multinomial(d;c)^n` to obtain the lower bound. Continuity and
positive approximation finish the proof.

If a limiting color has zero mass, at most o(n²) slots carry that color.
Replace all of them by one fixed positive-mass color in both directions of
comparison. The weight error is o(n²) because logK is bounded. The number
of choices of replacement positions and original labels is exp(o(n²)).
This proves the reduction asserted above, including normalization entropy.

## Uniform physical quantization for the seed kernel

Take `K_R(x,y)=exp[-t(||x||²+||y||²)]cosh(2t x^T R y)` for orthogonal R.
On the Euclidean ball of radius B,

`||gradient_x logK_R|| <= 2t||x||+2t||y|| <=4tB`,

and the same bound holds for y. Choose a finite eta-net INSIDE that ball
and map each point to a nearest representative. Then

`|logK_R(x,y)-logK_R(Qx,Qy)| <= 8tB eta`.

The segment remains in the ball, so this is a direct mean-value estimate.
Every graph product is therefore compared multiplicatively within
`exp[8tB eta binom(n,2)]`, uniformly over every array of source values.
This uses physical dimension q=k², not an infinite signed feature space;
all quantized kernels remain strictly positive. Exact row permutations
descend to exact color-count row permutations after quantization.

For unbounded arrays obeying the deterministic pervertex energy bound,
combine this bounded-domain step with the bad-incidence deletion argument
in decisive_bridge_global_coherent_truncation_2026_09_06.md. That deletion
is an UPPER-bound reduction to graphs with missing edges; it does not prove
an unbounded exact pressure limit. The director's heterogeneous/defective
graph formulation is the relevant next extension. In particular, do not
silently apply the homogeneous complete-graph theorem after deleting edges.

## Scope for actual sign constructions

In the grouped Hadamard weave, independent uniform block-column permutations
give the row-permutation mechanism conditionally on the spectra before those
permutations. The outgoing joint block empirical profiles generally depend
on the full spin tuples and on each group's basis randomness. They need not
be identical across vertices. Thus this theorem gives an exact homogeneous
sector and explains its self-coupling entropy; a full spin union bound needs
the heterogeneous-profile source entropy and cannot assume homogeneous
profiles are optimal without proving the needed concavity/compatibility.

There is a harmless slot-count detail: a block permutation acts on n blocks,
whereas the offdiagonal graph has n-1 slots. Conditional on the omitted self
block, the remaining slots are a uniform permutation of the remaining
multiset. Altering one omitted symbol per vertex changes only O(n) directed
colors. For a bounded positive alphabet kernel its energy and counting
effects are O(n log n), negligible at speed n². This does NOT justify
discarding the energetic within-group edges without checking the intended
direction of the defect inequality; discarding nonnegative defect is an
upper-bound step only.

## Exact source-enumerator interface

For a fixed deterministic group basis H and finite quantization, let N_H(p)
be the weighted number of retained k-spin tuples whose outgoing quantized
block multiset has empirical profile p. Precisely, sum over spin tuples the
probability, over the uniformly omitted self-block, of the resulting profile
p. This can be fractional. If the n-1 slot multiset is deterministic from
the outset, N_H is the ordinary integer count. If groups use this same basis
and independent uniform block permutations, the spin-summed positive defect
partition decomposes exactly as

`sum_(p_1,...,p_n) [product_i N_H(p_i)] Z_K(p_1,...,p_n)`.

Here Z_K is the independent row-permutation partition with the indicated
possibly heterogeneous profiles. This identity is for the selected positive
kernel factors; if within-group factors were discarded it represents that
upper-bound partition. With a fixed alphabet there are at most (n+1)^s
profiles, so the number of profile tuples is exp(O_s(n log n)). Thus its
leading log is the maximum of `sum_i log N_H(p_i)+log Z_K(profiles)` up to
O_s(n log n). No asymptotic limit of N_H is asserted. Independent group
basis randomness replaces each enumerator by its expectation, with the same
conditioning variables retained. This is an exact location of the missing
joint source-entropy input; inserting the scalar E certificate there without
an actual comparison theorem would be unjustified.
