# Strategic steering

Evidence cutoff: Wave 51, ledger §10.104 (2026-07-30). This early refresh
follows the spectral-excess theorem (10.1261) and the moment-support wall
(10.1255)--(10.1260). The next mandatory refresh is Wave 56, or earlier after
another decisive result, with a fresh blank-slate audit.

## User-stated research objective

Determine whether $\lim_{n\to\infty}M_n/n^{3/2}$ exists. The conjectural value
$1/2$ is not an additional user objective.

## Present judgment

No route proves convergence; the rigorous interval remains
$0.336493364431\ldots\le\liminf\le\limsup\le1/2$.

Adaptive optimized principal restriction remains the overarching framework,
and the bare arbitrary-cut tail (10.795) remains its sharp general sufficient
lemma. Project-row spectral excess is now its leading concrete implementation.
Wave 51 showed that the former non-strict coarea/FKN target was unnecessarily
strong and that the two- and three-moment completion mechanisms cannot bypass
the local principal recurrence.

## Leading sufficient package and convergence chain

Fix a compact window $[p_0,p_1]\subset(1/2,1)$, $m/n$ in that window, and
$0<c<1/4$. For an exact minimizer $A$, first prove a child-ground box witness

~~~math
\min_{S,\ y\text{ child ground}}
\min_{w\in\{\pm1\}^{S^c}}
\lVert A[:,S]y+A[:,S^c]w\rVert_2^2
\le R_0=O(n^{9/4-c}).
~~~

For the project-row class $C_R$, zero-deficit families $f_z$, and a
common-core kernel $K_\ell$, let $D_R=\mathbb E[\mathbf1_{C_R}a_z^2]$ and
let $P_\ell(R)$ be the aggregate triple retention (10.1248). The exact
missing overlap lemma is that some

~~~math
R_0\le R\le C(R_0+n^2),\qquad 2\le\ell<m
~~~

satisfies

~~~math
P_\ell(R)-\lambda_2(\ell)
\ge\exp\{-O(n^{3/4-c})\}.
~~~

The box witness gives $D_R>0$. Equation (10.1261) then gives a center with

~~~math
\max_{z\in C_R}a_z\ge
\frac{P_\ell(R)-\lambda_2(\ell)}
{(1-\lambda_1)+n(\lambda_1-\lambda_2)}
=\exp\{-O(n^{3/4-c})\}.
~~~

That center has the required project row and saved selector degree, hence
(10.795). The established inverse tail, geometric-window summability, and
exact landing force $q_n/n^{3/2}$ and therefore $M_n/n^{3/2}$ to converge.

## Evidence, obstructions, and falsification criteria

- The spectral-excess theorem is exact. Retention above $\lambda_2$ by the
  saved amount suffices; $P\ge\lambda_1$, constant harmonic gap, and slice
  FKN are not needed.
- The Wave 50 row-optimal cap-ten example already passes the sharper target
  at scales one and two. It falsifies only the overstrong first-eigenvalue
  package, not the convergence-scale overlap route.
- Scale one is circular up to density constants:
  $r_C\le P_1\le r_C/p$. Independent-resampling admixture is exactly
  redundant: its extraction bound never beats the pure-kernel certificate or
  the tautology $\max a_z\ge r_C$. The live scale must have
  $2\le\ell<m$.
- Cap inflation creates more completions of a fixed child ground but does not
  by itself create overlap among distinct favorable selectors. Spectral
  excess changes by a signed row-shell sum with no law-free monotonicity.
- An abstract singleton-family model defeats box-plus-cap-only reasoning, but
  is not an exact-minimizer counterexample. The leading route would be
  falsified by an unbounded exact-minimizer family that has a target box
  witness yet fails saved $P_\ell-\lambda_2$ at every allowed cap and every
  $2\le\ell<m$. The box component is separately falsified by uniform failure of
  the target row bound.
- On compact fixed-density windows with $t\ge0$, the entire Wave 50
  two-moment envelope is supported on $g\ge-O(n)$. Bonami and
  constant-probability linear--quadratic pairing reach only
  $g\ge-O(n^{5/4})$. Saved local-state mass at the target entropy and
  tolerance scales in either band already implies the principal recurrence
  by Hanson--Wright inversion.
- Cubic endpoint localization is support-dominated by the old reverse
  numerator. Moments through order three therefore do not supply a far-tail
  mechanism. These facts retire the implementations, not the full annealed
  incidence $Z_t$.
- For the same $c$, a live annealed alternative must obtain saved contribution
  not confined to $r=O(n^{3/2-c})$, with loss at most
  $\exp\{O(n^{3/4-c})\}$, or avoid reducing to saved local-state mass.
  Generic upper concentration and constant-probability reverse bounds are
  insufficient.
- Complement incidence, robust cylinders, positive-core first-eigenvalue
  mixtures, canonical Bellman pressure, strict pressure, high replicas,
  local harmonic Poincaré, monotone deletion, nested chains, greedoids, and
  constant shortfall without a square-root cavity reward remain retired.

## Most recent blank-slate abstraction audit

The Wave 45 audit began from the original problem before route comparison.
Its candidates and judgments are agent-authored hypotheses, not user
directives. This early decisive-result refresh does not regenerate it; Wave
56 will.

1. **Finite-temperature zero-sum interpolation.** Balanced $o(n)$ soft-minimax
   additivity would imply convergence, but temperature mismatch and overlap
   geometry still lack a signing-specific bridge.
2. **Summable cavity derivative.** Extension error
   $O(n^{1/2-\delta})$ would make normalized positive variation summable.
   This remains dormant until a square-root deletion reward is proved.
3. **Second-order signed-kernel compactness.** Continuity plus all-order
   recovery would prove convergence; existing compactness models lack the
   recovery operation.

## Ranked routes and Wave 52 direction

1. **Project-row box plus spectral excess.** Prove the box witness and saved
   $P_\ell-\lambda_2$ for some controlled cap and $2\le\ell<m$. Attack the
   overlap lemma through exact-minimizer selector intersections or a positive
   signed-shell budget, not through $K_0$ or scale one.
2. **Far-negative annealed completion tail.** Use regular linear coefficients,
   a multiscale signed spectrum, or an exact minimizer dichotomy to obtain
   saved contribution beyond the $O(n^{3/2-c})$ target local-margin band, or
   a lower bound that does not collapse to saved local-state mass.
3. **Box discrepancy itself.** Prove (10.1237) using cancellation between
   selected and outside columns; generic uniform completion remains too weak.
4. **Tight decomposition, dormant.** Reopen only with a non-switching global
   certificate absent from (10.1251)--(10.1252).
5. **Other global routes, dormant.** Reopen only with a testable composition,
   recovery, or square-root cavity theorem at the required scale.

Wave 52 must independently rank ten ideas and assign three falsifiable attacks.
