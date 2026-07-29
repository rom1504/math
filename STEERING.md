# Strategic steering

Evidence cutoff: `ledger.md` through Wave 18, §10.71 (2026-07-29). Next
mandatory regeneration: no later than the boundary after Wave 23, and earlier
after a decisive proof, counterexample, or change of leading route.

## User-stated research objective

Determine whether $`\lim_{n\to\infty}M_n/n^{3/2}`$ exists. The conjectural
value $`1/2`$ is not presently proved and is not being treated as an additional
user objective.

## Present assessment

No route currently proves convergence. The rigorous interval remains
$`0.336493364431\ldots\le\liminf\le\limsup\le1/2`$.

The leading route is **adaptive optimized principal restriction**, expressed
through the directed comparison tail (10.529). It is the weakest concrete
scale-transfer framework that directly matches the objective: it permits a
different exact minimizer and endpoint for each requested order, tolerates bad
local deletions, and asks only that the accumulated normalized cost vanish.

The constant-shortfall pressure criterion (10.617) is **not** the leading
route. It is the cleanest and most falsifiable local sufficient condition now
available, but it is substantially stronger than the endpoint statement
needed for convergence and currently has a severe generic-geometry barrier.
It remains the first ranked alternative and a useful source of tests.

## Leading route: adaptive optimized restriction

Write $`q_n=\min_A Q(A)=2M_n`$ and $`a_n=q_n/n^{3/2}`$. A concrete sufficient
lemma to seek is:

> **Power-saving restriction lemma.** There are constants
> $`0<\rho<1`$, $`c>0`$, $`C<\infty`$, and $`n_0`$ such that, for every
> $`n\ge n_0`$ and every integer $`m\in[\rho n,n)`$, some exact order-$`n`$
> minimizer $`A=A_{n,m}`$ and some $`m`$-vertex set $`S`$ satisfy

```math
Q(A[S])\le\left(\frac mn\right)^{3/2}q_n+C n^{3/2-c}.
```

Different targets may use different minimizers and subsets. Requiring one root
or one nested deletion chain for all targets would be unnecessarily stronger.

This lemma gives $`a_m\le a_n+O_\rho(n^{-c})`$. Repeated fixed-ratio steps,
with the last step landing at the requested order, have a geometric total
$`O(N^{-c})`$. Thus the adaptive distance in (10.529) satisfies
$`\Omega_{\rm ad}(N)\to0`$, which compares every tail order to a liminf
subsequence and proves convergence.

Candidate engines are the selected terminal-excess potential
(10.530)--(10.531), mean puncture stability (10.304)--(10.306), and a
matrix-specific upper bound for optimized restricted pressure (10.628)--(10.630).

Known obstructions:

- random restriction does not control the ground state selected after seeing
  the subset, and raw near-ground entropy is too large;
- the order-nine minimizer has no improving child, and excess can remain flat
  through two deletion levels;
- $`S\mapsto Q(A[S])`$ is neither submodular nor supermodular;
- fixed-temperature cavity rewards are an exact subset-lattice gradient, so
  reordering or banking rewards cannot beat endpoint optimization;
- Finner gives a lower bound on average child pressure, the wrong direction;
- a uniform $`o(n^{3/2})`$ error without a rate need not be summable over an
  arbitrarily long comparison path.

Falsification criterion for a proposed $`\rho`$: for some fixed
$`\rho_0\in[\rho,1)`$ and infinitely many $`n`$, prove that at
$`m\sim\rho_0n`$ every exact order-$`n`$ minimizer obeys

```math
\min_{|S|=m}\left[Q(A[S])-\left(\frac mn\right)^{3/2}q_n\right]
\ge c_0 n^{3/2}
```

for a constant $`c_0>0`$. Ruling out every possible $`\rho`$ requires such
barriers at ratios arbitrarily close to one. Even then, a sparser comparison
graph could satisfy (10.529); fully falsifying the adaptive route requires a
nonvanishing lower bound on every available directed comparison path.

## Audit of the constant-shortfall criterion (10.617)

To prove (10.617) literally, one must find fixed
$`\rho\in(0,1)`$, $`\beta>0`$, $`K<\infty`$, and $`n_0`$ such that every exact
order-$`n`$ minimizer, for every target $`m\in[\rho n,n)`$, has a deletion path
$`B_n,\ldots,B_m`$ on which at each order $`r`$ some deleted coordinate obeys

```math
\kappa_{\beta,i}(B_r)
\ge\frac{q_n}{n^{3/2}}\left[r^{3/2}-(r-1)^{3/2}\right]-K.
```

For convergence alone, “every minimizer” may be weakened to a target-specific
choice of exact root minimizer. The temperature and shortfall must remain
uniform, and the path must land at the requested order.

By (10.623)--(10.624), this is exactly a coordinate Hellinger near-isolation
bound with affinity $`e^{-\Theta(\sqrt r)}`$. Telescoping gives (10.618), an
$`O(n^{-1/2})`$ normalized cost per fixed-ratio window, hence convergence.

Generic cube log-Sobolev, support entropy, antipodality, and layer projection
cannot reach that scale. The $`A_9`$ wall only forces $`K>0.886740446\ldots`$;
it is finite and does not falsify a universal constant.

For the literal every-minimizer criterion, a decisive root falsifier is, for
each fixed $`\beta>0`$, an unbounded sequence containing an exact minimizer
with $`\max_i\kappa_{\beta,i}=o(\sqrt n)`$. Falsifying the convergence-weakened
version requires all exact minimizers at those orders to be bad. Equivalently,
their coordinate affinities are $`e^{-o(\sqrt n)}`$. For a proposed
$`(\beta,K)`$, a mandatory cutset with no outgoing edge meeting (10.617) also
falsifies the required path statement.

## Ranked alternatives

1. **Constant-shortfall Hellinger/cavity route.** Stronger than needed but
   exact, scalar, tail-summable, and sharply falsifiable. It needs genuinely
   quadratic, exact-minimizer ancestry beyond generic Gibbs geometry.
2. **Spatial partition plus temporal service.** The target (10.600) is
   algebraically feasible at critical scale, but the symplectic witness is
   nonminimal and has zero raw resource. It still needs two theorems: force a
   useful partition from global minimality, then convert correlated response
   into causal allocation/Hall service.
3. **Universal sharp lower bound.** Proving
   $`Q(A)\ge(1-o(1))n^{3/2}`$ for every signing would combine with conference
   upper bounds to give the limit and its value. The verified field-plus-spin,
   capped-field, and $`A^2`$ gains are real, but one-probe rounding ceilings,
   manufactured heavy rows, and unresolved positive-heavy tails block the
   current machinery.
4. **Global signing-space entropy or thermodynamic compactness.** A full
   planted-direction-sensitive pressure/LDP or uniform amplification theorem
   would settle convergence. Fixed replicas, spectral moments, scalar pressure
   axioms, finite boundary states, and ordinary graphon/traffic summaries have
   all proved too coarse, so this remains a long-range route.

Mixed conference products, finite flat-child examples, and further small-order
searches should be used primarily to falsify proposed lemmas. They are not
currently positive routes unless they reveal an asymptotic minimizer-specific
mechanism.

## Decision rule for the next wave

Prioritize attempts that prove, weaken, or asymptotically falsify the
power-saving restriction lemma. Keep one independent probe of (10.617) or the
spatial route when it tests a genuinely minimizer-specific mechanism. Do not
reopen generic entropy, scalar telescoping, random-restriction, or purely
spatial arguments without a stated way around the obstructions above.
