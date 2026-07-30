# Strategic steering

Evidence cutoff: Wave 31, §10.84 (2026-07-30). Regenerate by Wave 36, or earlier after a decisive proof, counterexample, or change of leading route.

## User-stated research objective

Determine whether $`\lim_{n\to\infty}M_n/n^{3/2}`$ exists. The conjectural value $`1/2`$ is not proved and is not an additional user objective.

## Present judgment

No route proves convergence. The rigorous interval remains $`0.336493364431\ldots\le\liminf\le\limsup\le1/2`$.

The leading route is **adaptive optimized principal restriction**, implemented by a **fractional row-good block-coset cover**. Wave 31 weakens its mesoscopic endpoint: one coset need only hit a fraction `1/T` of the batch, and `T<=n^eta` still proves convergence when `eta<c_0`. Its sharpest structured successor is now low-row forest clustering of favorable completions.

The constant-shortfall pressure criterion (10.617) is **not** the strongest route. Existing signed-cycle, plaquette, cavity, and deficit audits neither produce nor asymptotically falsify its required square-root reward.

## Constant-shortfall audit

The exact sufficient result is fixed `rho,beta>0` and `K<infinity` such that every target pair has an exact root and deletion path landing at `m`, with

```math
\kappa_{\beta,i}(B_r)\ge
\frac{q_n}{n^{3/2}}[r^{3/2}-(r-1)^{3/2}]-K
```

at every step. Pressure telescoping then proves convergence. A literal fixed-temperature falsifier is an unbounded exact-minimizer family with `max_i kappa_(beta,i)=o(sqrt(r))`; falsifying the path version requires this for every eligible root or at an unavoidable deletion cutset. None is known.

## Leading route and exact sufficient lemma

Put $`q_n=2M_n`$, $`a_n=q_n/n^{3/2}`$, and $`p=m/n`$. Convergence follows if fixed $`\rho\in[1/2,1)`$, `c>0`, and `K` exist such that every large `n` and $`m\in[\rho n,n)`$ admit an exact minimizer `A` and `m`-set `S` with

```math
Q(A[S])\le p^{3/2}q_n+Kn^{3/2-c}.
```

Then $`a_m\le a_n+O_\rho(n^{-c})`$; geometric-window errors are summable and exact landing proves convergence. The arbitrary-cut lemma (10.795) supplies this edge from one cut having row cost `O(n^(9/4-c))` and selector coverage `e^{-O(n^(3/4-c))}` at tolerance `O(n^(3/2-c))`.

For the structured implementation, choose `c_0 in (0,1/4)` and put

```math
L_0=n^{3/4-c_0},\quad k_0=\Theta(L_0/\log n),\quad
r=\left\lceil\frac{n\log(2k_0)}{L_0}\right\rceil.
```

Let `H_a` be the hit sets of all eligible row-good cosets. For `T<=n^eta`, `eta<c_0`, set `s=ceil(r/T)`. The sharp bare sufficient lemma is, uniformly for every selector law `w`,

```math
\sum_aw(H_a)^s\ge e^{-O(rL_0)}.
```

Equivalently at this exponent, an iid batch has probability `e^{-O(rL_0)}` that one coset hits `s` positions. The partial-collision theorem gives

```math
\log\tau_*\le O(TL_0)=O(n^{3/4-(c_0-\eta)}).
```

Greedy rounding and pigeonhole yield the arbitrary cut with degraded saving `c'=c_0-eta>0`, hence convergence.

The sharp completion-based sufficient event is stronger but concrete: with probability `e^{-O(rL_0)}` against every `w`, partition the batch into `T<=n^eta` groups; in each group choose favorable completions joined to an auxiliary `R_2=O(n^2)` center by a projective-Hamming tree of total length `D=O(k_0)`. The verified forest theorem gives

```math
J\le1+D,\qquad
C_{\rm sig}\le2R_2(x^0)+8\|A\|_{\rm op}^2D
=O(n^{9/4-c_0}).
```

Wave 30 then supplies one eligible row-good coset per group. Proving this adversarial-law forest event would therefore prove convergence.

## Obstructions and falsification

- Pointwise planting gives `T=r`, costing `Theta(n log n)`. Local completion cylinders and hashing alone can remain at `e^{-Theta(n)}` coverage.
- A transitive abstract system whose candidates are all fixed-size subsets has exact KKT, regular degrees, and every fixed-order common cover, yet its partial-collision probability is `e^{-omega(rL_0)}` throughout the sufficient `T` range. Generic incidence or saddle-point optimality is insufficient.
- Low-row centers are abundant by `E_U R_2=n(n-1)`, but no theorem places favorable principal completions within `O(k_0)` total forest distance of fewer than `n^{c_0}` such centers.
- Principal-ground exchange cycles control signed first-order shore totals. `C_sig` depends on a quadratic boundary-susceptibility vector; a checked weighted ground-state example separates the two. Discrete signing minimality must add genuinely new information.
- Uniform-selector means remain circular. Pair/triple agreement and low-information pointwise selection do not give the mesoscopic statistic.
- All bounds must hold for every target order in one fixed ratio window, uniformly against every selector law.

For fixed constants, the fractional implementation is falsified by infinitely many target pairs for which every exact minimizer has `tau_*>e^{K n^(3/4-c)}`. At the batch level, an adversarial law with

```math
\sum_aw(H_a)^s=e^{-\omega(rL_0)}
```

for every sufficient parameter choice is equivalent evidence. This does not by itself falsify arbitrary fixed-cut coverage. A direct restriction falsifier still requires minimum excess `Omega(n^(3/2))` along infinitely many orders at one fixed ratio.

## Ranked alternatives

1. **Matched parent-Gibbs interpolation.** Control the centered harmonic Doob energy (10.920) or centered flatness (10.921), plus selector Hellinger, at fixed density. Reverse-KL first moments and naive block summation lose common modes; `A_6` gives a finite mechanism wall.
2. **Bare minimizer-specific power sum.** Prove (10.907) directly without coherent completions, by showing exact-signing hit incidence cannot resemble the symmetric fixed-subset wall.
3. **Quadratic boundary susceptibility.** Derive a discrete-minimality inequality controlling (10.914) directly; ordinary exchange-cycle signs are insufficient.
4. **Conditional row-Laplace alignment.** Combine a boundary scalar reverse tail with (10.878)--(10.880); two minimizer-specific inputs remain.
5. **Constant shortfall or global replacement/compactness.** Exact formulations remain falsifiable, but their missing inputs are currently less local and testable.

## Decision rule

Wave 32 should attack the adversarial-law low-row forest event, independently attack the power sum using exact signing structure, and retain one centered harmonic-screening route. Do not restart all-`r` coherence, pointwise planting, fixed-order gluing, every-refinement control, uncentered harmonic accumulation, uniform-mean truncation, generic parent LSI, or unsigned-cycle bounds as if their recorded walls were absent.
