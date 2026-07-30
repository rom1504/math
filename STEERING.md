# Strategic steering

Evidence cutoff: Wave 29, §10.82 (2026-07-30). Regenerate by the Wave 35 boundary, or earlier after a decisive proof, counterexample, or change of leading route.

## User-stated research objective

Determine whether $`\lim_{n\to\infty}M_n/n^{3/2}`$ exists. The conjectural value $`1/2`$ is not proved and is not an additional user objective.

## Present judgment

No route proves convergence. The rigorous interval remains $`0.336493364431\ldots\le\liminf\le\limsup\le1/2`$.

The leading route is **adaptive optimized principal restriction**. Its weakest endpoint is one arbitrary full cut with low row square and stretched-exponential centered-selector coverage, (10.795). Its leading constructive implementation is now the **fractional row-good block-coset cover** (10.873)--(10.875), not one globally coherent coset. Wave 29 showed that fractional mixing can stay strong when every deterministic coset misses selectors, while an adjacent-triangle wall defeats pairwise gluing at a fixed cap.

The constant-shortfall pressure criterion (10.617) is **not** the strongest route. It is sufficient but stronger than convergence requires: it demands a fixed-temperature `Theta(sqrt(r))` deletion reward, up to a constant, along a complete landing path. Signed-cycle, plaquette, cavity, and finite-deficit audits have not supplied that reward. No asymptotic family falsifies it, so it remains a sharply testable secondary route.

## Leading route and exact sufficient lemma

Put $`q_n=2M_n`$, $`a_n=q_n/n^{3/2}`$, $`p=m/n`$, and $`p_2=(m)_2/(n)_2`$. Convergence follows if fixed $`\rho\in[1/2,1)`$, `c>0`, and `K` exist such that every large `n` and every $`m\in[\rho n,n)`$ admit an exact order-`n` minimizer `A` and `m`-set `S` with

```math
Q(A[S])\le p^{3/2}q_n+Kn^{3/2-c}.
```

This gives $`a_m\le a_n+O_\rho(n^{-c})`$; geometric-window errors are summable and exact target landing through (10.529) proves convergence.

For an oriented cut `d`, define

```math
\widehat\ell(S,d)=Q(A[S])-c_A(S,d)-p_2\Delta_A(d)-B_{n,m},
\qquad B_{n,m}=(p^{3/2}-p_2)q_n.
```

The arbitrary-cut lemma (10.795) asks, uniformly in one active ratio window, for

```math
R_2(d)=O(n^{9/4-c}),\qquad
U_m\{\widehat\ell(S,d)\le O(n^{3/2-c})\}\ge e^{-O(n^{3/4-c})}.
```

Conditioning on this event in (10.794), with $`\lambda\asymp n^{-3/4}`$, proves the power-saving restriction above.

The exact structured lemma now sought is (10.875). Let $`L=n^{3/4-c}`$, let $`\mathscr A_C`$ contain all eligible balanced block cosets whose whole support has $`R_2\le C=K_Rn^{9/4-c}`$, and put

```math
H_a=\{S:Q(P_a^{\mathsf T}D_aH_SD_aP_a)
\ge Y_A(S)-K_tn^{3/2-c}\}.
```

For fixed `rho,c` and uniform finite constants, prove every large target pair admits an exact minimizer for which

```math
\tau_*=\min\left\{\sum_au_a:
u_a\ge0,\ \sum_{a:S\in H_a}u_a\ge1\ \forall S\right\}
\le e^{K L}.
```

By minimax, every selector law then gives mass at least `e^{-KL}` to one row-good coset hit set. Greedy rounding produces an `e^{O(L)}` low-row codebook covering every selector; pigeonhole gives the fixed cut above. Thus this lemma proves convergence.

## Obstructions and falsification

- Uniform-selector means are circular: every fixed cut has mean effective loss equal to the unknown restriction excess. Linear row prices, mean child entropy, and pointwise low-information completion do not escape this.
- Uniform random diagonals have hit probability `e^{-omega(L)}` above the target boundary. At the boundary their Hanson--Wright exponent only scale-matches, and row-good conditioning can erase a rare hit set.
- Pointwise planting gives only `tau_*<=binom(n,m)`, with logarithm `Theta(n)`. The missing theorem is overlap against every adversarial selector law, not self-incidence.
- Deterministic coherence needs few coordinate signatures and low-row closure under every balanced refinement. The exact `A_9` triangle is pairwise but not triple compatible at cap `80`; signature count plus low-row ternary span is insufficient.
- A scalar aligned reverse tail need not survive the row cap. Conditional row-Laplace control (10.878)--(10.880) is an exact repair, but is a second missing input beyond the scalar tail.
- Every statement must hold for all target orders in one fixed ratio window, with uniform constants and summable landing costs.

For fixed constants, the coset implementation is falsified by infinitely many target pairs for which every exact minimizer has $`\tau_*>e^{KL}`$. Failure at `e^{omega(L)}` for every fixed constant and every candidate `rho,c` defeats the implementation, but not arbitrary fixed-cut coverage. The latter is falsified by maximal eligible-cut coverage `e^{-omega(L)}` for every fixed row and tolerance constant. A direct restriction falsifier is a fixed ratio with infinitely many orders where every exact minimizer has minimum restriction excess `Omega(n^{3/2})`.

## Audit of constant shortfall (10.617)

A convergence-sufficient version asks for fixed `rho,beta>0` and finite `K` such that every target pair admits an exact root and deletion path landing at `m`, with every order-`r` step selecting `i` so that

```math
\kappa_{\beta,i}(B_r)
\ge \frac{q_n}{n^{3/2}}[r^{3/2}-(r-1)^{3/2}]-K.
```

Pressure telescoping gives (10.618), an `O(n^{-1/2})` cost per fixed-ratio window, hence convergence. The signed-cycle form is (10.780), and the cavity sufficient condition is $`v_i\ge(9/16)\alpha^2r`$ in (10.819).

Known walls are precise: competitor floors point the wrong way, absolute cycle coefficients erase cancellation, sharp plaquette curvature yields only constant reward, and audited minimizers fail the pointwise cavity premise at several temperatures. A literal fixed-temperature falsifier is an unbounded exact-minimizer family with $`\max_i\kappa_{\beta,i}=o(\sqrt r)`$. Falsifying the target-specific path version requires this for all eligible roots or at an unavoidable deletion cutset. Deficit profiles (10.820) and (10.847) are the cleanest tests; no such family is known.

## Ranked alternatives

1. **Matched parent-Gibbs interpolation.** Prove both estimates (10.800) at `beta=gamma=Theta(n^{-1/2+c})`. Matching removes the old linear wall on audited minimizers, but parent entropy contains a common centered-row-square mode invisible to selector Hellinger variation.
2. **Higher-order deterministic coset agreement.** Construct globally low-signature witnesses and control every balanced refinement with slack. Pairwise Johnson compatibility is explicitly insufficient.
3. **Conditional row-Laplace alignment.** Combine a boundary-scale scalar reverse tail with (10.878)--(10.880). It is exact but currently requires two independent minimizer-specific inputs.
4. **Constant-shortfall signed-cycle/cavity route.** Direct and sharply falsifiable, but no complete-signing mechanism supplies its square-root deletion reward.
5. **External-surplus replacement or global compactness.** Exact Hall and replacement formulations survive, as do universal sharp lower-bound and thermodynamic possibilities, but their missing inputs are less local and testable.

## Decision rule

Prioritize the adversarial-selector dual of the fractional cover, with one independent attempt on higher-order/refinement structure and one on the matched parent common mode. Use literature on agreement, list recovery, fractional covers, or discrepancy only after matching hypotheses and exponents. Do not restart uniform-coset tails above the boundary, independent completion, pairwise gluing, scalar mean truncation, linear row pricing, generic parent LSI, unsigned-cycle bounds, capped boosting, or two-ground mixing as if their recorded walls were absent.
