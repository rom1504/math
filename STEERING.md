# Strategic steering

Evidence cutoff: Wave 52, ledger §10.105 (2026-07-30), after (10.1266)--(10.1270)
and (10.1275)--(10.1278). Next mandatory refresh: Wave 57, with a blank-slate audit.

## User-stated research objective

Determine whether $\lim_{n\to\infty}M_n/n^{3/2}$ exists. The conjectural value
$1/2$ is not an additional user objective.

## Present judgment

No route proves convergence; the rigorous interval remains $0.336493364431\ldots\le\liminf\le\limsup\le1/2$.

Adaptive optimized principal restriction remains the overarching framework;
the bare arbitrary-cut tail (10.795) is its sharp general sufficient lemma.
Large-core spectral excess leads structurally, while the dependence-safe
K-functional tail is a precise direct alternative. Both lack abundance theorems.

## Leading structured package and convergence chain

Fix a compact window $[p_0,p_1]\subset(1/2,1)$, $m/n$ in that window, and
$0<c<1/4$. Put $L_n=n^{3/4-c}$. First prove a child-ground box witness

~~~math
\min_{S,\ y\text{ child ground},\ w\in\{\pm1\}^{S^c}}
\lVert A[:,S]y+A[:,S^c]w\rVert_2^2
\le R_0=O(n^{9/4-c}).
~~~

For the project-row class $C_R$, zero- or admissible positive-deficit
families $f_z$, and common-core kernel $K_\ell$, let $D_R>0$ be the
double-incidence mass, $P_\ell(R)$ the aggregate retention, and
$\rho_\ell=\binom m\ell/\binom n\ell$. The genuinely spectral missing lemma
is that some scale and cap satisfy

~~~math
R_0\le R\le C(R_0+n^2),\qquad -\log\rho_\ell=\omega(L_n),\qquad
P_\ell(R)-\lambda_2(\ell)\ge\exp\{-O(L_n)\}.
~~~

Equation (10.1261) then gives a center of degree $\exp\{-O(L_n)\}$ and row
$O(n^{9/4-c})$, hence (10.795). The inverse tail, geometric-window
summability, and exact landing force $q_n/n^{3/2}$ and therefore
$M_n/n^{3/2}$ to converge.

## Exact direct alternative

For a local state, write $Z=L+Q_T$, $L=\beta^{\mathsf T}w$, and
$r=-g/p_2>0$. Put

~~~math
H=\lceil L_n\rceil,\qquad T_n=n^{3/2-c},\qquad
\mathcal K_H(\beta)=K_{1,2}(\beta,\sqrt H),\qquad
b_H=C_{\rm HW}\{\lVert A[S^c]\rVert_F\sqrt H+
\lVert A[S^c]\rVert_{\rm op}H\}=O(T_n).
~~~

Montgomery-Smith plus a dependence-safe Hanson-Wright union bound gives

~~~math
r+b_H\le\tfrac12\mathcal K_H(\beta)
\quad\Longrightarrow\quad
\Pr\{Z\le-r\}\ge\tfrac12\,12^{-H}.
~~~

For a deterministic $\omega_n\uparrow\infty$, let $G_n(\omega)$ be the local
states satisfying this condition and $r\ge\omega_nT_n$. An open sufficient
direct lemma is

~~~math
\text{For some }\omega_n\uparrow\infty\text{ and }C_0<\infty,\qquad
\nu_m(G_n(\omega))\ge e^{-C_0H}
~~~

uniformly over the required exact minimizers, window, and tolerances. It yields
$Z_t\ge\exp\{-O(H)\}$ and convergence through (10.795); Holmstedt makes the
profile explicit up to universal constants.

## Evidence, obstructions, and falsification criteria

- Core-load counting gives
  $M\ge\max\{\rho_\ell P_\ell,[P_\ell-\lambda_2]_+/D_\ell\}$, where
  $D_\ell=(1-\lambda_1)+n(\lambda_1-\lambda_2)$. Every $\ell=O(L_n)$ is direct degree
  counting at the target exponent, not spectral amplification. A distinct
  spectral proof must use $-\log\rho_\ell=\omega(L_n)$; a linear core suffices.
- At level two, $P_2\ge\lambda_2(2)$ directly forces
  $M=\Omega(n^{-2})$ on the fixed-density window. All stored target-side
  tests of this nonnegative criterion are inconclusive because
  $h_2\ge\lambda_2$ makes the criterion automatic; strict excess can contain
  off-diagonal collisions. At fixed density the self-loop is exponentially small.
- Maximal-selector port inequalities give a real positive-threshold
  certificate for only polynomially many one-exchange labels; it need not
  apply to an arbitrary box-witness selector. Every one-seed sublinear-radius
  ball has $e^{o(n)}$ labels. A10 has no adjacent row-ten selector pair and no
  port preserving $q_*$ by saturation, although 3,200/3,840 maximizing-sign
  exchanges have zero deficit at their smaller selector cap.
- The structured route would be falsified by an unbounded exact-minimizer
  family that has target box witnesses but fails saved large-core excess at
  every controlled cap. The box component is separately falsified by uniform
  failure of the target row bound. No such exact family is known.
- The K-functional theorem is exact and makes no independence error. Its
  success requires cross energy $\gg n^{9/4-c}$, but high cross energy alone
  does not force the required head/diffuse coefficient profile or saved
  far-margin population.
- The direct implementation would be falsified if, for every diverging
  $\omega_n$ and every usable $c$, tolerance, and compact window, uniformly
  $\log\nu_m(G_n(\omega))/H\to-\infty$. That would not falsify the full
  annealed incidence $Z_t$.
- Exact parent-cap pairing removes $Q_T$ pointwise but is recurrence-circular:
  its K-functional condition forces $e\ge r$. Two- and three-moment bounds,
  Bonami, and constant-probability pairing remain trapped in the local band.
- Complement incidence, robust cylinders, $K_0$ mixing, small-core vocabulary
  changes, canonical Bellman pressure, strict pressure, high replicas, local
  harmonic Poincaré, nested chains, greedoids, and constant shortfall without
  a square-root cavity reward remain retired.

## Most recent blank-slate abstraction audit

The Wave 45 audit began from the original problem before route comparison.
Its candidates and judgments are agent-authored hypotheses, not user
directives. This early decisive-result refresh does not regenerate it; Wave
57 will.

1. **Finite-temperature zero-sum interpolation.** Balanced $o(n)$ soft-minimax
   additivity would work, but temperature mismatch and overlap lack a bridge.
2. **Summable cavity derivative.** Error $O(n^{1/2-\delta})$ would make positive
   variation summable; this awaits a square-root deletion reward.
3. **Second-order signed-kernel compactness.** Continuity plus all-order recovery
   would work; existing compactness models lack that recovery operation.

## Ranked routes and Wave 53 direction

1. **Box plus genuinely large-core spectral excess.** Prove saved excess at
   $-\log\rho_\ell=\omega(L_n)$ using global selector geometry or a signed-shell
   budget; small cores and one-seed sublinear exchange balls are not substitutes.
2. **Far K-functional completion abundance.** Prove the exact-minimizer state
   bound for $G_n(\omega)$, or derive a checkable coefficient-profile dichotomy that
   forces it on the high-cross branch.
3. **High-cross versus box discrepancy.** Use (10.1278) to separate states
   where the K route can operate from low-cross child grounds, then control
   the internal Gram needed for the box witness in the latter branch.
4. **Tight decomposition, dormant.** Reopen only with a non-switching global
   certificate absent from (10.1251)--(10.1252).
5. **Other global routes, dormant.** Reopen only with a testable composition,
   recovery, or square-root cavity theorem at the required scale.

Wave 53 must independently rank ten ideas and assign three falsifiable attacks.
