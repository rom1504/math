# Strategic steering

Evidence cutoff: Wave 33, §10.86 (2026-07-30). Regenerate by Wave 38, or earlier after a decisive change.

## User-stated research objective

Determine whether $`\lim_{n\to\infty}M_n/n^{3/2}`$ exists. The conjectural value $`1/2`$ is not an additional user objective.

## Present judgment

No route proves convergence. The rigorous interval remains $`0.336493364431\ldots\le\liminf\le\limsup\le1/2`$.

The leading route is still **adaptive optimized principal restriction**, via
the **fractional row-good block-coset cover**. Its bare sufficient statistic is
the partial-collision power sum (10.907). Its strongest convex, dualizable
center certificate is the soft capacity (10.947)--(10.950), a one-kernel
sufficient surrogate for the logically weaker hard total-cost event (10.926).

Wave 33 changes one diagnosis: the specific sufficiently-large-`K` planted
construction cannot be an exact minimizer, by (10.943)--(10.944). Missing is a conditional theorem
transferring this parent constraint to prescribed favorable child fibers while
controlling orientation and witness migration.

The constant-shortfall criterion (10.617) is **not** the strongest route. No
evidence supplies its square-root reward along a deletion path, while the
fractional and harmonic routes have sharper exact reductions.

## Leading route and exact sufficient lemma

Put $`q_n=2M_n`$, choose fixed $`c_0\in(0,1/4)`$, and set

```math
L_0=n^{3/4-c_0},\qquad
k_0=\Theta(L_0/\log n),\qquad
r=\left\lceil\frac{n\log(2k_0)}{L_0}\right\rceil.
```

For `T<=n^eta`, `eta<c_0`, let `s=ceil(r/T)` and `H_a` be the eligible
row-good-coset hit sets. Uniformly over every selector law `w`, it suffices that

```math
\boxed{\sum_a w(H_a)^s\ge e^{-O(rL_0)}.}
```

This gives an iid batch where one coset hits `s` positions. Partial collision,
rounding, and pigeonhole give a power-saving restriction edge; geometric-window
summability and exact landing then prove convergence.

The strongest concrete structured sufficient lemma uses
$`K_\lambda(z,S)=e^{-\lambda a_z(S)}`$ and

```math
V_{s,\lambda}
=\inf_{u\in\Delta(\Omega)}\|K_\lambda u\|_{L^s(\nu_2)}
=\max_{h\ge0,\ \|h\|_{L^q(\nu_2)}\le1}
  \min_S\mathbb E_{\nu_2}[h(z)e^{-\lambda a_z(S)}].
```

Here `q=s/(s-1)`. With $`D=\Theta(k_0)`$ and
$`\lambda=\Lambda rL_0/(D+1)`$, it is enough to prove, for uniform constants
`0<A<Lambda` and every relevant exact-minimizer target pair, the
existence of one common weight `h` (independent of `S` and `w`) such that

```math
\|h\|_{L^q(\nu_2)}\le1,\qquad
\min_S\mathbb E_{\nu_2}[h(z)e^{-\lambda a_z(S)}]\ge e^{-ATL_0}.
```

Then $`V_{s,\lambda}^s\ge e^{-AsTL_0}`$; the exact margin
`Lambda r>A sT` makes it beat $`e^{-\Lambda rL_0}`$ and proves (10.926), the
bare power sum, and convergence. A common hard core shows this certificate is
weaker than a uniform normalized-degree theorem.

## Known obstructions and falsification criteria

- The common parent law (10.945) admits a negative-orientation shield, and its
  witnesses need not lie in prescribed child fibers; exact `A_9` exhibits
  witness migration.
- Soft capacity is nearly hard: $`\lambda(D+1)=\Omega(rL_0)`$. One abstract
  assignment makes it $`e^{-\Omega(n)}`$ simultaneously for every
  `lambda>=lambda_min`; any proof must use signing-minimizer structure.
- Uniform degrees, first moments, fixed-order agreement, generic KKT incidence,
  and pointwise planting do not supply shared overlap.
- A tree certificate that pays Johnson coordinate entry is exponentially too
  rare under the uniform selector law by (10.954)--(10.955). Only exact
  terminal-min-cut screening, which may exploit free coordinates, survives.
- Generic harmonic likelihood algebra does not bound screened response. The
  abstract square in (10.962) makes it unbounded; finite exact examples defeat
  `C_scr<=1`, pointwise restoring drift, and omission of the orientation edge.
- Degree-two Hamming noise contracts row-square excess only by `1-o(1)` at
  relevant radii, and conference blocks contain a high-row basin.

For fixed sufficient parameters, fractional coverage is falsified by laws with
$`\sum_a w(H_a)^s=e^{-\omega(rL_0)}`$. Falsifying the whole implementation
requires witnesses for every fixed choice. Failure of the center or forest
event does not falsify the bare power sum. Direct falsification still requires
$`\Omega(n^{3/2})`$ along infinitely many orders at one fixed ratio.

The adaptive soft-kernel certificate is falsified if, along actual-minimizer
targets, $`\sup_{\lambda\in\mathcal A_n}V_{s,\lambda}^s
\le e^{-\omega(rL_0)}`$ for each fixed sufficient parameter choice. This does
not falsify (10.926) or the bare power sum; abstract labels are insufficient.

## Harmonic alternative

Matched parent-Gibbs interpolation is the strongest independent route. Its
centered Doob target (10.920) is implied by a uniform restoring-defect estimate
(10.961), `kappa_t>=kappa_0>0`, plus
$`\int_0^1t\mathcal E_t(g)\,dt=O(n^{1/2-2c})`$. Then (10.900) gives the parent
entropy target. Selector Johnson-Hellinger control remains separate.

## Constant-shortfall audit

It suffices to find fixed `rho,beta>0`, `K<infinity`, and for every target pair
an exact root and deletion path landing at `m`, with at every step

```math
\kappa_{\beta,i}(B_r)\ge
\frac{q_n}{n^{3/2}}[r^{3/2}-(r-1)^{3/2}]-K
```

Pressure telescoping then proves convergence. A falsifier needs
$`\max_i\kappa_{\beta,i}=o(\sqrt r)`$ on an unbounded exact-minimizer family;
the path version requires this at every root or an unavoidable deletion
cutset. Neither is known.

## Ranked alternatives

1. **Conditional replacement/common weighted centers.** Transfer (10.945) to
   prescribed favorable fibers and prove (10.950).
2. **Bare minimizer-specific power sum.** Prove (10.907) directly without
   committing to centers or coherent completions.
3. **Matched parent resonant response.** Prove (10.961), resonant control, and
   the separate selector-Hellinger bound.
4. **Exact terminal-min-cut screening.** Compress coordinate conflict paths
   without paying Johnson entry for every selector move.
5. **Quadratic boundary susceptibility / conditional row-Laplace.** Extract a
   usable child-fiber inequality from discrete global minimality.
6. **Constant shortfall or global replacement/compactness.** Keep as exact
   fallbacks, but current missing inputs are less local and less testable.

## Next-wave decision

Wave 34 should independently test conditional witness transfer, a direct
common-weight capacity theorem, and minimizer-specific resonant/restoring
control. Do not restart normalized-degree counting, Johnson-entry trees,
small-radius noise smoothing, posterior-variance-only screening, uncentered
reverse-KL accumulation, or generic incidence arguments without a mechanism
that explicitly escapes their recorded obstruction.
