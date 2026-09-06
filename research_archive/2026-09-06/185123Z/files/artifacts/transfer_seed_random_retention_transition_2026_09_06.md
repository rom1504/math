# One parent family has good high-retention and bad sparse random restrictions

Date: 2026-09-06. Seed-transfer track. This is an actual random-restriction
transition, not a claim about nonconvergence of the original minima.

## Theorem

There is a sequence of deterministic symmetric sign Hadamard matrices
`W_D`, of orders `D=m^2` tending to infinity, with the following properties.
Let `A_D` be their hollowings and put `p=31/32`.

1. `Q(A_D)/D^(3/2)->1/2`.
2. If `T` is a uniform subset of EXACT size `N=pD`, then

```math
\Pr\{Q((A_D)_T)/N^{3/2}\le0.4995\}\longrightarrow1.        (1)
```

3. On the SAME deterministic parent sequence, for every size sequence
   `n->infinity`, `n/D->0`, and every fixed `eta>0`, a uniform exact
   `n`-subset `U` satisfies

```math
\Pr\{Q((A_D)_U)/n^{3/2}\ge2/\pi-\eta\}\longrightarrow1.    (2)
```

The number `0.4995` in (1) can be replaced by any fixed
`c>c_*`, where `c_*=1/2-a/(8 sqrt(31/32))<0.499432220485404`
is the pre-cluster strict upper constant and
`a=91470529542342299/20460000000000000000`.
Different choices of `c` may initially choose different parent families.
A final diagonal selection can make the high-retention cap have
probabilistic limsup at most `c_*` on a single sequence; no uniform
growing-depth estimate is used.

This disproves a possible extension of the sparse-random lower theorem
to all fixed retentions below one. It shows an actual high-retention
improvement on the parent cap, followed by substantial degradation at
vanishing retention. It does NOT imply that `M_n/n^(3/2)` oscillates.

## 1. What the existing strict-upper ensemble actually supplies

Fix a constant `c_b` strictly between `c_*` and the desired `c<1/2`.
The standalone reconstruction
`transfer_reconstruction_standalone_2026_09_06.md`, Section 6,
supplies a fixed-depth ensemble of full symmetric Hadamard weaves

```math
W_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i),\quad i,j,a,b\in[m]. (3)
```

Take `m` along arbitrarily large admissible Hadamard orders divisible
by 32 and put `k=pm`, so `N=mk=pD`. For every FIXED balanced selector
`T'=union_i({i} times T_i')` with `|T_i'|=k`, the ensemble has

```math
\Pr_W\{Q((\operatorname{hollow}W)_{T'})>c_bN^{3/2}\}
\le e^{-\kappa m^2}                                      (4)
```

for all sufficiently large admissible `m`, with some `kappa>0`.
This is a high-probability statement, not merely existence of one
successful realization.

Here is the exact source of the margin in (4). Choose `a'<a` and
`eta<a'/4` so that `1/2-eta/(2 sqrt(p))<c_b`. At a sufficiently deep
FIXED recursion, the row estimate and all-spin Markov bound in that
proof give failure probability
`exp[(-a'+4 eta+o(1))m^2]`. On its complement the normalized hollow
cap is at most `1/2-eta/(2 sqrt(p))+1/(2 sqrt(N))`.
The strict margins absorb the diagonal cost and all error terms.
The bound is uniform in the fixed balanced selectors because each
row input has exactly the same empirical source and the ensemble has
fresh input signed-permutation invariance. No selector is chosen
after viewing a favorable realization in this step.

Average (4) over the product-uniform balanced selector. The mean,
over the parent ensemble, of its balanced bad-selector fraction is at
most `exp(-kappa m^2)`. Hence a deterministic weave `W_D` exists for
which

```math
\Pr_{T'\ \mathrm{balanced}}
\{Q((A_D)_{T'})>c_bN^{3/2}\}\le e^{-\kappa m^2}.           (5)
```

The choice is of ONE parent for each order. It does not assert the
same property for every Hadamard parent.

## 2. Exact coupling from a uniform selector to a balanced selector

Partition `[D]` into `m` fibres of size `m`. Choose `T` uniformly among
all `N=mk` subsets and write `K_i=|T intersect fibre i|`.
Conditional on the count vector, the sets in the fibres are independent
and uniform at their prescribed sizes.

Repair each fibre independently as follows. If `K_i>k`, retain a
uniform `k`-subset of its selected coordinates. If `K_i<k`, add a
uniform `(k-K_i)`-subset of its unselected coordinates. If `K_i=k`,
make no change. Call the resulting set `T'`.

Conditional on EVERY count vector, the repaired `T_i'` are independent
uniform `k`-subsets, by permutation symmetry within each fibre.
Their conditional law does not depend on the counts. Thus the marginal
law of `T'` is exactly the product-uniform balanced law used in (5).
The symmetric difference is exactly

```math
h:=|T\mathbin{\triangle}T'|=\sum_{i=1}^m|K_i-k|.
```

The hypergeometric mean of `K_i` is `k` and its variance is

```math
\operatorname{Var}K_i
=k(1-1/m)\frac{D-N}{D-1}\le k.
```

Therefore

```math
\mathbb Eh\le
B_m:=m\sqrt{k(1-1/m)(D-N)/(D-1)}
\le m\sqrt k\le m^{3/2}.                                (6)
```

In particular `h/D->0` in probability. This coupling does not use
approximation by independent vertex sampling, a condition such as
`N^2/D->0`, or independence of the fibre counts.

## 3. Exact cap continuity for nearby equal-size selectors

For any real symmetric hollow parent `A` and equal-size selectors
`T,T'` of size `N`,

```math
|Q(A_T)-Q(A_{T'})|\le\|A\|_{op}\sqrt{N|T\triangle T'|}.    (7)
```

Proof. Take any spin on `T`; put zeros outside `T`. Make a spin on
`T'` agreeing on the intersection and arbitrary on the new coordinates,
again putting zeros outside. For their ambient vectors `x,y`,
`||x-y||=sqrt(h)` and `||x+y||<=2 sqrt(N)`. Symmetry gives

```math
\tfrac12|x^TAx-y^TAy|
=\tfrac12|(x-y)^TA(x+y)|\le\|A\|_{op}\sqrt{Nh}.
```

Take absolute values inside the two caps, maximize, and reverse the
roles of the selectors. This proves (7). No maximizer information is
required.

For the Hadamard parent, `||A_D||op<=m+1`. Combining (5)--(7) gives
the explicit finite bound, for every `epsilon>0`,

```math
\Pr_T\left\{\frac{Q((A_D)_T)}{N^{3/2}}>c_b+\epsilon\right\}
\le e^{-\kappa m^2}
 +\frac{(m+1)^2B_m}{\epsilon^2N^2}.                       (8)
```

The second term is `O_p(epsilon^(-2)m^(-1/2))` at fixed `p`, so (8)
proves (1) by taking `epsilon=c-c_b`.

## 4. Both ends belong to the same parent sequence

Every full weave (3) is a symmetric Hadamard: column orthogonality in
the inner fibre sum, followed by row orthogonality, gives `W_D^2=D I`.
For an even `m`, pair its fibre indices by a fixed-point-free involution
`pi`. The Boolean vector

```math
x_{i,a}=u_iH_i(a,\pi(i)),\qquad
u_{\pi(i)}=\sigma S_{i,\pi(i)}u_i
```

satisfies `W_D x=sigma m x`, for either `sigma=+1` or `-1`.
Spectral upper and these actual Boolean eigenvectors imply that
hollowing has normalized cap `1/2+O(1/m)`. In particular the selected
parents satisfy a uniform bounded-cap hypothesis.

Now apply
`transfer_seed_sparse_bounded_cap_parent_theorem_2026_09_06.md`
to these SAME deterministic parents. It gives (2) at every simultaneous
vanishing-retention scale. No relation between the two independently
sampled selectors in (1) and (2) is asserted or required.

For completeness, to realize the endpoint `c_*` on one parent sequence,
choose successively smaller fixed margins, their fixed recursion depths,
and admissible orders so large that both terms in (8) are less than the
corresponding error tolerance. This is an ordinary diagonal subsequence
argument. It neither proves a growing-depth type bound nor changes
the original all-order limit proof.

## Interpretation for seed transfer

The random restriction operation is not uniformly bad or uniformly
good. Its performance depends essentially on retention and parent
structure. The high-retention phenomenon supplies a concrete family of
good selectors occupying probability tending to one, while the sparse
theorem proves that the same behavior cannot persist down to vanishing
retention. This does not yet propagate an arbitrary liminf seed, nor
identify a critical retention or the exact cap curve.

The finite coupling law and nearby-selector cap inequality are checked
by `computations/transfer_seed_retention_transition_verify_2026_09_06.py`.
The inherited strict-upper certificate remains in its canonical source;
the finite checks here do not replace its analytic proof.
