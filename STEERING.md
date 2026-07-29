# Strategic steering

Evidence cutoff: Wave 23, §10.76 (2026-07-29). Regenerate by Wave 28, or earlier after a decisive result or route change.

## User-stated research objective

Determine whether $`\lim_{n\to\infty}M_n/n^{3/2}`$ exists. The conjectural value
$`1/2`$ is not proved and is not an additional user objective.

## Present judgment

No route proves convergence. The rigorous interval remains
$`0.336493364431\ldots\le\liminf\le\limsup\le1/2`$.

The leading route remains **adaptive optimized principal restriction**. Its best
engine is now the regular capped weak-learner/codebook criterion
(10.741)--(10.742), tolerating `O(n^{3/4-c})` selector complexity rather than
`O(n^{1/2-2c})`. This strongest proved interface still has an open premise.

The constant-shortfall criterion (10.617) is **not** the strongest route. It is
stronger than convergence requires, and Waves 19--23 show generic Gibbs, cover,
and edge-flip arguments pointing the wrong way. It remains the cleanest sharply
falsifiable local alternative if a proof uses complete-signing/minimizer structure.

## Leading route and exact sufficient lemma

Put $`q_n=2M_n`$ and $`a_n=q_n/n^{3/2}`$. The route would close if there are
fixed $`\rho\in[1/2,1)`$, $`c>0`$, and $`C<\infty`$ such that, for all large
`n` and every $`m\in[\rho n,n)`$, some exact order-`n` minimizer `A` and some
`m`-set `S` satisfy

```math
Q(A[S])\le (m/n)^{3/2}q_n+Cn^{3/2-c}.
```

Different targets may use different minimizers and subsets. This gives
$`a_m\le a_n+O_\rho(n^{-c})`$. Fixed-ratio iterations have a geometrically
summable total error, and target-specific landing gives the adaptive tail
(10.529), hence convergence.

The strongest exact implementation is this joint statement: for fixed
$`0<c<1/4`$ and ratio window, choose an exact minimizer `A` whose row-regular
parent-ground game satisfies, with constants independent of `n,m`,

```math
\begin{aligned}
B&=O(n^{3/4-c}), & t&=O(n^{3/2-c}),\\
\varepsilon&\asymp n^{-1/2-c}, &
\alpha_{t,\varepsilon}^{\mathrm{gr,cap}}(B)
&\ge \exp\{-O(n^{3/4-c})\}.
\end{aligned}
```

The capped minimax theorem then produces a parent-ground codebook of log-size
`O(n^{3/4-c})` and average distortion `O(n^{3/2-c})`. Substitution in (10.742)
with $`\lambda\asymp n^{-3/4}`$ proves the displayed restriction lemma.

The general fallback is parent tilted-tail overlap (10.727): find a selector law
whose KL cost plus average $`-\log\Omega_{\beta,S}(t)`$ is `O(n^{1/2-2c})`, with
$`t=O(n^{3/2-c})`. It avoids a row cap but demands much stronger overlap.

## Obstructions and falsification

Known obstructions to the leading engine are:

- the expected-loss weak learner is exactly the unknown optimized-restriction
  excess, so ordinary boosting is circular;
- fixed-slice concentration is a converse to coverage: a pre-existing leading
  gap forces exponentially small coverage. Present constants make this
  unconditional only below density `0.452911...`, outside the active window;
- overlap for one deletion controls only the deleted row. Heavy rows can
  migrate, and simple trimming loses the full critical `n^{3/2}` scale;
- the finite `A_9` table shows row regularity and selector coverage cannot be
  proved separately and then combined;
- conditional free energy is smoother by a square root, but under the wrong
  marginal and at an insufficient scale. The exact `A_9` limit has vanishing
  base variance with a positive tilted Jensen gap, so only full interpolation
  or tilted-tail control can work;
- a comparison must cover every target order with summable costs; one good
  subset, root, deletion chain, or uncontrolled landing order is insufficient.

For fixed proposed constants, this engine is falsified by infinitely many
active-ratio pairs where every exact minimizer has no admissible regular grounds
or has
$`\alpha_{t,\varepsilon}^{\mathrm{gr,cap}}(B)
=\exp\{-\omega(n^{3/4-c})\}`$ throughout the permitted parameter scales.
That would not falsify all optimized restriction. A direct falsifier of its
power-saving form is a fixed ratio and infinitely many `n` for which every
exact minimizer has minimum restriction excess at least $`c_0n^{3/2}`$.
Falsifying the full adaptive route requires a nonvanishing cost on every
available directed comparison path, not merely failure of this engine.

## Audit of constant shortfall (10.617)

The exact result needed is fixed $`\rho,\beta>0`$, finite `K`, and, for every
large exact minimizer and every $`m\in[\rho n,n)`$, a deletion path landing at
`m` on which each order-`r` step has

```math
\kappa_{\beta,i}(B_r)\ge
\frac{q_n}{n^{3/2}}[r^{3/2}-(r-1)^{3/2}]-K.
```

For convergence, “every minimizer” may be weakened to a target-specific exact
root. The temperature, shortfall, fixed-ratio window, and landing statement
must remain uniform. Equation (10.618) then costs only `O(n^{-1/2})` per
window, a summable geometric tail. Equivalently, the selected coordinate needs
Gibbs Hellinger affinity $`e^{-\Theta(\sqrt r)}`$.

Evidence now weighs against proofs based only on generic local geometry:
pressure stability, Hölder aggregation, witness incidence, simultaneous-flip
covers, nonuniform nonnegative weights, and the exact weighted-cover dual all
give the wrong inequality direction. All nonnegative star moments and pair
detailed balance permit only a constant parity reward. This is not a complete
signing counterexample, so (10.617) remains open.

For the literal criterion, an asymptotic falsifier is a fixed `beta` and an
unbounded sequence containing an exact minimizer with
$`\max_i\kappa_{\beta,i}=o(\sqrt n)`$. For the convergence-weakened version,
all exact minimizers at those orders must be bad. For a proposed `(beta,K)`, a
mandatory deletion cutset with no qualifying outgoing edge also falsifies the
required path.

## Ranked alternatives

1. **Constant-shortfall Hellinger/cavity route.** Direct, tail-summable, and
   sharply falsifiable, but it now requires a new complete-signing congestion
   or competitor-partition upper bound rather than another star aggregation.
2. **Terminal external surplus with spatial/temporal service.** Exact
   replacement and Hall formulations exist; the missing input is an
   exact-minimizer upper tail that localizes surplus and makes service causal.
3. **Universal sharp lower bound.** A proof
   $`q_n\ge(1-o(1))n^{3/2}`$ would combine with conference upper bounds and
   settle the limit. Existing field, spectral, and capped gains remain below
   the positive-heavy-tail barrier.
4. **Global thermodynamic or signing-space compactness.** A planted-sensitive
   pressure/LDP or amplification theorem could settle convergence, but scalar
   pressure, fixed replicas, finite boundary states, and ordinary graphon or
   traffic summaries are already known to be too coarse.

## Decision rule

Attack row-regular parent-ground existence and capped coverage jointly, not as
independent lemmas. Keep a parent tilted-tail attempt active as the general
fallback. Test constant shortfall only through complete-signing structure that
could reverse the established wrong-way inequalities. Do not recycle generic
entropy, scalar pressure, random restriction, base variance, star-cover, or
purely spatial arguments without an explicit mechanism overcoming their
recorded obstruction.
