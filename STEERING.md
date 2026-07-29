# Strategic steering

Evidence cutoff: Wave 24, §10.77 (2026-07-29). Regenerate by Wave 29, or earlier after a decisive result or route change.

## User-stated research objective

Determine whether $`\lim_{n\to\infty}M_n/n^{3/2}`$ exists. The conjectural value
$`1/2`$ is not proved and is not an additional user objective.

## Present judgment

No route proves convergence. The rigorous interval remains
$`0.336493364431\ldots\le\liminf\le\limsup\le1/2`$.

The leading route remains **adaptive optimized principal restriction**, but Wave 24 strictly weakened its best concrete implementation: seek global selector--parent-ground overlap with a captured row-square budget, not a capped worst-selector learner or large codebook. Its exact optimizer uses at most two positive-coverage grounds.

The constant-shortfall criterion (10.617) is **not** the strongest route; it is stronger than convergence requires. The complete-signing competitor family now has an exact signed-cycle formulation, but every proved generic upper-affinity mechanism gives only constant reward instead of `Theta(sqrt(n))`. It remains the cleanest sharply falsifiable local alternative.

## Leading route and exact sufficient lemma

Put $`q_n=2M_n`$ and $`a_n=q_n/n^{3/2}`$. Convergence follows if there are
fixed $`\rho\in[1/2,1)`$, $`c>0`$, and `C` such that, for all large `n` and
every $`m\in[\rho n,n)`$, some exact order-`n` minimizer `A` and `m`-set `S`
satisfy

```math
Q(A[S])\le(m/n)^{3/2}q_n+Cn^{3/2-c}.
```

Different targets may use different minimizers and subsets. This gives $`a_m\le a_n+O_\rho(n^{-c})`$; fixed-ratio iterations have a geometric error sum, and exact target landing gives the adaptive tail (10.529), hence convergence.

The strongest implementation to prove is (10.764)--(10.766). For a parent
ground `g`, set

```math
u_g=\Pr_{S\sim U_m}\{Q(A[S])-c_A(S,g)\le t\},
\qquad c_g=R_2(g)=\sum_i r_i(g)^2.
```

For fixed $`0<c<1/4`$, seek a target-specific exact minimizer, threshold
$`t=O(n^{3/2-c})`$, and a law `nu` on at most two parent grounds such that

```math
Z=\sum_g\nu_gu_g\ge\exp\{-O(n^{3/4-c})\},
\qquad
\frac{\sum_g\nu_gu_gc_g}{Z}=O(n^{9/4-c}).
```

Conditioning `U_m x nu` on the good event costs exactly `-log Z`; the selector KL plus mutual information is no larger. The domain-free slice mgf (10.759), with $`\lambda\asymp n^{-3/4}`$, then proves the displayed power-saving restriction lemma through (10.761).

Equivalently, for row budget `C`, prove

```math
Z_*(C)=\min_{\theta\ge0}\max_g
u_g[1+\theta(C-c_g)]
\ge\exp\{-O(n^{3/4-c})\},
\qquad C=O(n^{9/4-c}).
```

The fallback within the leading route is the parent-Gibbs global tail (10.771), requiring overlap `exp{-O(n^(1/2-2c))}` without row improvement. Its exact full-interpolation alternative is endpoint transport (10.774)--(10.775).

## Obstructions and falsification

Known obstructions are:

- switching rewrites the loss as replenishment $`\ell(S,g)=R_T(g)+b_T(g)-d_T`$, but its mean is exactly the unknown optimized-restriction excess and the scalar constraints force no lower tail;
- edge-flip and block-replacement witnesses can have positive parent deficit,
  so they need not supply columns of the exact-ground overlap LP;
- the concentration converse (10.769) matches the `n^(3/4)` exponent: a genuine leading gap forces overlap too small for the sufficient lemma, but current constants do so only below density `0.452911...`;
- the finite `A_9` hard-cap separation disappears in the soft LP, but this is
  not asymptotic evidence; captured compatibility can still favor high `R_2`;
- parent interpolation needs exponential rare-tail transport; base variance, unweighted selector gradients, mean truncation, and generic temperature monotonicity do not provide it;
- every comparison must cover all requested target orders with summable costs;
  one uncontrolled landing order or one nonuniform window is insufficient.

For proposed constants, the implementation is falsified by infinitely many active-ratio pairs where every exact minimizer has $`Z_*(O(n^{9/4-c}))=\exp\{-\omega(n^{3/4-c})\}`$ for every allowed threshold. A fixed positive leading restriction gap implies this through (10.769), but does not falsify all optimized restriction. A direct power-saving falsifier is a fixed ratio and infinitely many `n` where every exact minimizer has restriction excess at least $`c_0n^{3/2}`$. Falsifying the adaptive route requires nonvanishing cost on every available directed comparison path.

## Audit of constant shortfall (10.617)

The exact result needed is fixed $`\rho,\beta>0`$, finite `K`, and, for every
large exact minimizer and every $`m\in[\rho n,n)`$, a path landing at `m` whose
order-`r` step selects `i` with

```math
\kappa_{\beta,i}(B_r)\ge
\frac{q_n}{n^{3/2}}[r^{3/2}-(r-1)^{3/2}]-K.
```

For convergence, “every minimizer” may be weakened to a target-specific exact root. Temperature, shortfall, ratio window, and landing must remain uniform; (10.618) then costs `O(n^(-1/2))` per geometric window.

Wave 24 identifies the exact complete-signing object. With
$`\rho_\beta=\tanh(2\beta)`$ and signed even-cycle polynomial
$`\mathcal P_B`$,

```math
\operatorname{BC}_i=
\frac{\mathcal P_{B[-i]}(\rho_\beta)}
{\cosh(2\beta)^{r-1}\mathcal P_B(\rho_\beta)}.
```

Thus (10.780) is the exact deletion-ratio theorem needed. Exact minimality gives only lower floors on edge-character twists, while absolute cycle bounds destroy enormous signed cancellation. Sharp two-spin plaquette curvature gives `max kappa >= beta^(-1) log(2/C_beta)`, only a constant and therefore `o(sqrt(r))`.

For the literal criterion, an asymptotic falsifier is fixed `beta` and an unbounded sequence containing an exact minimizer with $`\max_i\kappa_{\beta,i}=o(\sqrt n)`$. The convergence-weakened version needs all exact minimizers at those orders to be bad. A mandatory deletion cutset with no qualifying edge also falsifies a proposed `(beta,K)` path statement.

## Ranked alternatives

1. **Parent-Gibbs full interpolation.** Part of optimized restriction but independent of exact-ground row control; it needs the exponential endpoint transport bound (10.774)--(10.775).
2. **Constant-shortfall signed-cycle route.** Direct and sharply falsifiable, but requires minimizer-specific cancellation stability under deletion.
3. **Terminal external surplus with spatial/temporal service.** Exact Hall and replacement formulations exist; the missing input is a localized exact-minimizer surplus upper tail.
4. **Universal sharp lower bound or global thermodynamic compactness.** Either could settle the problem, but heavy-tail, scalar-pressure, fixed-replica, and finite-profile methods remain below the required resolution.

## Decision rule

Prioritize direct lower bounds or structural dichotomies for the two-ground
LP `Z_*(C)`. Keep one attack on the full exponential endpoint transport.
Test constant shortfall only through signed-cycle deletion/cancellation or
another genuinely complete-signing upper mechanism. Do not restart capped
worst-selector boosting, hard row regularity, base-variance, mean-truncation,
star-cover, or coefficientwise-cycle arguments as if their recorded walls
were absent.
