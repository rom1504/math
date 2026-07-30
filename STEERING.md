# Strategic steering

Evidence cutoff: Wave 30, §10.83 (2026-07-30). Regenerate by Wave 35, or earlier after a decisive proof, counterexample, or change of leading route.

## User-stated research objective

Determine whether $`\lim_{n\to\infty}M_n/n^{3/2}`$ exists. The conjectural value $`1/2`$ is not proved and is not an additional user objective.

## Present judgment

No route proves convergence. The rigorous interval remains $`0.336493364431\ldots\le\liminf\le\limsup\le1/2`$.

The leading route is **adaptive optimized principal restriction**, whose
sharpest constructive implementation is the **fractional row-good
block-coset cover**. Wave 30 makes this an exact mesoscopic collision problem
and proves that one controlled refinement follows from two batch statistics.

The constant-shortfall pressure criterion (10.617) is **not** the strongest
route. It demands a fixed-temperature `Theta(sqrt(r))` deletion reward along
a complete landing path. Existing audits neither supply that reward nor give
an asymptotic falsifier, so it remains a sharply testable secondary route.

## Constant-shortfall audit

The exact sufficient result is fixed `rho,beta>0` and `K<infinity` such that
every target pair has an exact root and a deletion path landing at `m`, with

```math
\kappa_{\beta,i}(B_r)\ge
\frac{q_n}{n^{3/2}}[r^{3/2}-(r-1)^{3/2}]-K
```

at every step. Pressure telescoping then proves convergence. A literal
fixed-temperature falsifier is an unbounded exact-minimizer family with
`max_i kappa_(beta,i)=o(sqrt(r))`; falsifying the path version requires this
for every eligible root or at an unavoidable deletion cutset. None is known.

## Leading route and exact sufficient lemma

Put $`q_n=2M_n`$, $`a_n=q_n/n^{3/2}`$, $`p=m/n`$, and
$`p_2=(m)_2/(n)_2`$. Convergence follows if fixed
$`\rho\in[1/2,1)`$, `c>0`, and `K` exist such that every large `n` and every
$`m\in[\rho n,n)`$ admit an exact order-`n` minimizer `A` and `m`-set `S`
with

```math
Q(A[S])\le p^{3/2}q_n+Kn^{3/2-c}.
```

This gives $`a_m\le a_n+O_\rho(n^{-c})`$; geometric-window errors are
summable and exact target landing proves convergence.

The weakest current endpoint remains the arbitrary-cut lemma (10.795): find
one cut `d` with

```math
R_2(d)=O(n^{9/4-c}),\qquad
U_m\{\widehat\ell(S,d)\le O(n^{3/2-c})\}
\ge e^{-O(n^{3/4-c})}.
```

For its structured implementation, put `L=n^{3/4-c}`, use
`k=Theta(L/log n)` blocks, and let `H_a` be the hit set of each eligible
row-good coset. The fractional target is

```math
\tau_*=\min\left\{\sum_au_a:
u_a\ge0,\ \sum_{a:S\in H_a}u_a\ge1\ \forall S\right\}
\le e^{O(L)}.
```

Wave 30 gives an exponent-equivalent exact target. Set

```math
r=\left\lceil\frac{n\log(2k)}L\right\rceil
=\Theta(n^{1/4+c}\log n)
```

and let `J_r` be the event that `r` iid selectors share a hit coset. Prove,
uniformly for every selector law `w`,

```math
\Pr_{w^{\otimes r}}(J_r)\ge e^{-O(rL)}.
```

The `(2k)^n` candidate bound and collision sandwich yield
`tau_*=e^{O(L)}`. Greedy rounding and pigeonhole then supply the arbitrary
fixed cut, hence the restriction edge and convergence.

The sharp completion-based sufficient lemma is now:

> Against every `w`, an `r`-sample has probability `e^{-O(rL)}` of admitting
> favorable completions whose coordinate-signature count is `J=O(L/log n)`
> and whose coarsest signature coset has
> `C_sig=O(n^{9/4-c})`.

The constrained-refinement theorem turns those statistics into one eligible
row-good common coset. This lemma is stronger than bare `J_r`, not equivalent.

## Obstructions and falsification

- Uniform-selector means are circular: every fixed cut has mean effective
  loss equal to the unknown restriction excess. Linear row pricing and mean
  child entropy do not escape it.
- Pointwise planting gives only `tau_*<=binom(n,m)`, with logarithm
  `Theta(n)`. The abstract independent-label wall (10.897)--(10.898) proves
  that local completion cylinders and hashing alone can genuinely remain at
  `e^{-Theta(n)}` coverage. Exact-minimizer correlations are indispensable.
- Pairwise and triple common cover are insufficient. Both exact `A_9`
  systems have `J_2=J_3` surely for every law, no globally covering coset,
  and fractional values below one. The required order is mesoscopic.
- Generic agreement/list-recovery theorems begin from a prescribed overlap
  premise and do not provide adversarial-law incidence, signature entropy,
  or the coarsest quadratic cap. Perfect global list consistency can still
  have `J=n` coordinate signatures.
- The coarsest cap `C_sig` is necessary. Wave 30 proves that, with the target
  signature count, it is also sufficient up to allowed slack for one
  size-capped refinement; do not reinstate “every refinement” as an open
  requirement.
- Scalar aligned reverse tails need not survive a row cap. The exact
  row-Laplace repair still requires both a scalar tail and a minimizer-specific
  conditional covariance estimate.
- Every asymptotic statement must hold for every target order in one fixed
  ratio window, with uniform constants and summable landing costs.

For fixed constants, the fractional implementation is falsified by infinitely
many target pairs for which every exact minimizer has `tau_*>e^{KL}`, or
equivalently at exponent scale by a law with `Pr(J_r)=e^{-omega(rL)}`. A
direct restriction falsifier requires minimum excess `Omega(n^{3/2})` along infinitely many orders at one fixed ratio.

## Ranked alternatives

1. **Matched parent-Gibbs interpolation.** Control the exact harmonic-cost variance (10.900) and selector Hellinger term at fixed density; one deletion works for `c<=1/12`, while `A_6` defeats centered-`R_2` reduction.
2. **Bare mesoscopic common-coset incidence.** Prove `J_r` directly from exact minimization without the stronger coherent-completion hypothesis.
3. **Conditional row-Laplace alignment.** Combine a boundary scalar reverse tail with (10.878)--(10.880); two minimizer-specific inputs remain.
4. **Constant-shortfall signed-cycle/cavity route.** Direct and falsifiable, but no mechanism supplies its square-root reward.
5. **External-surplus replacement or global compactness.** Exact formulations survive, but their missing inputs are less local and testable.

## Decision rule

Wave 31 should put one attack on minimizer-specific control of `J` and `C_sig`,
one independent attack on bare mesoscopic collision probability,
and one on fixed-density harmonic parent cost. Check the collision law
against adversarial selector distributions, not only the uniform slice. Do
not restart self-incidence, independent completion, fixed-order gluing,
every-refinement control, scalar mean truncation, generic parent LSI,
unsigned-cycle bounds, or two-ground mixing as if their recorded walls were
absent.
