# Strategic steering

Evidence cutoff: Wave 32, §10.85 (2026-07-30). Regenerate by Wave 37, or earlier after a decisive proof, counterexample, or change of leading route.

## User-stated research objective

Determine whether $`\lim_{n\to\infty}M_n/n^{3/2}`$ exists. The conjectural value $`1/2`$ is not proved and is not an additional user objective.

## Present judgment

No route proves convergence. The rigorous interval remains $`0.336493364431\ldots\le\liminf\le\limsup\le1/2`$.

The leading route remains **adaptive optimized principal restriction**, implemented by a **fractional row-good block-coset cover**. Its decisive bare statistic is the partial-collision power sum (10.907). Its most concrete structured successor is now the low-row center-star moment (10.926), while the exact forest geometry is the coherent-orientation terminal min-cut (10.923)--(10.924).

The new complete-signing wall makes the source of difficulty unusually clear: the known `Q`, operator-norm, favorable-tolerance, child-optimality, and global low-row-abundance facts do not force favorable fibers near low-row centers. A successful structured proof must use exact **global signing minimality** in a way not yet present in the ledger.

The constant-shortfall pressure criterion (10.617) is still **not** the strongest route. No new evidence produces its square-root reward or asymptotically falsifies it.

## Leading route and exact sufficient lemma

Put $`q_n=2M_n`$, choose fixed $`c_0\in(0,1/4)`$, and set

```math
L_0=n^{3/4-c_0},\qquad
k_0=\Theta(L_0/\log n),\qquad
r=\left\lceil\frac{n\log(2k_0)}{L_0}\right\rceil.
```

For `T<=n^eta`, `eta<c_0`, let `s=ceil(r/T)` and let `H_a` be the hit sets of all eligible row-good cosets. The sharp bare sufficient lemma is, uniformly over every selector law `w`,

```math
\boxed{
\sum_a w(H_a)^s\ge e^{-O(rL_0)}.
}
```

This is equivalent at the required exponent to an iid batch having probability `e^{-O(rL_0)}` that one coset hits `s` positions. The partial-collision theorem gives `log tau_*<=O(TL_0)=O(n^(3/4-(c_0-eta)))`; rounding and pigeonhole then produce an arbitrary cut with a power-saving restriction edge. Geometric-window summability and exact landing prove convergence.

The strongest currently testable center-based sufficient statement is the total-cost moment. Let `nu_2` be uniform on `C_2={z:R_2(z)<=2n(n-1)}` and let `a_z(S)` be the minimum projective distance from `z` to a favorable completion for `S`. Prove, uniformly over `w`,

```math
\boxed{
\mathbb E_{z\sim\nu_2}
\Pr_{S_1,\ldots,S_s\sim w}
\left\{\sum_{j=1}^sa_z(S_j)\le Ck_0\right\}
\ge e^{-O(rL_0)}.
}
```

On this event the completions form a low-row star of total length `O(k_0)`; the verified forest bound and balanced refinement produce the required row-good coset. A more general forest may succeed even if this star moment fails.

The exact geometric alternative is to partition the batch into `T<=n^eta` groups and, in every group, choose favorable partial labels, a low-row root, a tree, and shared node orientations such that the total coordinate terminal min-cut is `O(k_0)`. Proving this event with probability `e^{-O(rL_0)}` is also sufficient.

## Obstructions and falsification

- A normalized-degree/Jensen proof of the center-star moment would require essentially `2^{m-o(n)}` favorable labels per selector. Low-row centers occupying half the cube does not help a fixed conditional fiber.
- Independent random selector labels give an abstract system with every center degree `e^{-Omega(n)}` and center moment `e^{-omega(rL_0)}`. Cross-selector correlation is indispensable.
- A planted clique of size `Kn^(3/4)` plus a low-norm random signing remainder is a complete signing with `Q=O(n^(3/2))` and `||A||op=O(n^(3/4))`; for every fixed-density child `S` containing the planted block, its full project-favorable fiber remains `Omega(n^(3/4))` from every `O(n^2)`-row center. It is not a global minimizer, so it isolates—rather than falsifies—the missing use of exact global minimality.
- Completion forests are exactly coordinate terminal separators after one shared orientation per node. Pairwise overlap consistency and signed exchange cycles do not bound the resulting conflict-path packing.
- Degree-two Hamming noise contracts row-square excess only by the two-coordinate correlation. At every relevant sublinear radius the multiplier is `1-o(1)`; a conference-block signing has a whole high-row basin.
- Uniform-selector averages, fixed-order agreement, KKT incidence, pointwise planting, and generic codebook entropy remain blocked by the recorded walls.

For any fixed sufficient choice `(c_0,eta,T(n))`, failure is witnessed along infinitely many target pairs by laws `w_n` with `sum_a w_n(H_a)^s=e^{-omega(rL_0)}`. Falsifying the entire fractional implementation requires such a sequence for every fixed sufficient parameter choice; the witnessing laws and subsequences may depend on that choice. Failure of the center-star or forest event alone does not falsify the bare coset power sum or arbitrary-cut coverage. A direct restriction falsifier still needs minimum excess `Omega(n^(3/2))` along infinitely many orders at one fixed ratio.

## Harmonic route

Matched parent-Gibbs interpolation remains the strongest independent route. Selector constants are pure gauge. The exact parent target is the centered Doob energy (10.920); the heat-bath screened-response bound (10.940) is an exact sufficient local criterion which retains the combined selector-common drift and posterior-KL cancellation. Posterior variance alone misses common response. The adjacent-selector Johnson-Hellinger estimate is a separate required input.

Centered flatness (10.921) is demoted as a general proxy: on the exact minimizer `A_6,m=3`, its oscillation grows like `4 beta` while the actual endpoint entropy decays like `e^(-4 beta)`. This is a finite-parameter mechanism wall, not a project-scale asymptotic falsifier.

## Constant-shortfall audit

The exact sufficient result remains fixed `rho,beta>0` and `K<infinity` such that every target pair has an exact root and deletion path landing at `m`, with

```math
\kappa_{\beta,i}(B_r)\ge
\frac{q_n}{n^{3/2}}[r^{3/2}-(r-1)^{3/2}]-K
```

at every step. Pressure telescoping would prove convergence. A literal fixed-temperature falsifier is an unbounded exact-minimizer family with `max_i kappa_(beta,i)=o(sqrt r)`; falsifying the path version requires this at every eligible root or an unavoidable deletion cutset. None is known.

## Ranked alternatives

1. **Global-minimality cut compression / conditional fibers.** Derive an edge-replacement or signing-variation inequality that rules out the planted localized quadratic mode in exact minimizers and forces (10.926) or the terminal-min-cut event.
2. **Bare minimizer-specific power sum.** Prove (10.907) directly, without coherent completions or low-row centers, by excluding the symmetric and independent-label incidence walls.
3. **Matched parent screened response.** Bound (10.940) for the endpoint likelihood, retain drift--KL cancellation, and separately prove selector Hellinger.
4. **Quadratic boundary susceptibility.** Use discrete global minimality to control (10.914); ordinary exchange-cycle signs are insufficient.
5. **Conditional row-Laplace alignment.** Combine a boundary scalar reverse tail with (10.878)--(10.880); both minimizer-specific inputs remain open.
6. **Constant shortfall or global replacement/compactness.** Exact formulations remain falsifiable, but their missing inputs are less local and testable.

## Decision rule

Wave 33 should test a global-minimality variation against localized quadratic modes, independently attack cross-selector correlation in the center-star moment, and retain one combined screened-response route. Do not restart normalized-degree counting, small-radius noise smoothing, coordinatewise projective orientations, posterior-variance-only screening, uncentered reverse-KL accumulation, pointwise planting, or generic incidence arguments as if their recorded walls were absent. Do not infer flatness from selector centering or use it as an entropy proxy without a new project-scale minimizer mechanism.
