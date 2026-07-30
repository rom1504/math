# Strategic steering

Evidence cutoff: Wave 50, ledger §10.103 (2026-07-30). This early refresh
follows (10.1251)--(10.1252). The next mandatory refresh is Wave 55, or
earlier after another decisive result, with a fresh blank-slate audit.

## User-stated research objective

Determine whether $\lim_{n\to\infty}M_n/n^{3/2}$ exists. The conjectural value
$1/2$ is not an additional user objective.

## Present judgment

No route proves convergence; the rigorous interval remains
$0.336493364431\ldots\le\liminf\le\limsup\le1/2$.

Adaptive optimized principal restriction remains the leading framework, and
the bare arbitrary-cut tail (10.795) remains its sharp sufficient lemma.
Annealed completion incidence is its cleanest direct scalar implementation;
project-row non-strict coarea is the strongest structured alternative. The
canonical Bellman-pressure implementation of tight decomposition is retired:
its full pressure is a switching-symmetry average, not a descent direction.

## Leading sufficient lemma and convergence chain

Fix a compact window $[p_0,p_1]\subset(1/2,1)$, $m/n$ in that window, and
$c\in(0,1/4)$. For an exact minimizer $A$, let $\Pi_n$ be uniform on oriented
projective cuts and define

~~~math
Z_t=(U_m\otimes\Pi_n)
\{(S,d):\widehat\ell(S,d)\le t\}.
~~~

The leading direct sufficient lemma is

~~~math
t=O(n^{3/2-c}),\qquad -\log Z_t=O(n^{3/4-c}).
~~~

Conditioning the product law on this event costs $-\log Z_t$ in relative
entropy. Equations (10.1229)--(10.1230) then give one cut with

~~~math
R_2(d)=O(n^{9/4-c}),\qquad
U_m\{S:\widehat\ell(S,d)\le t\}\ge Z_t/2.
~~~

This is (10.795). The inverse tail gives

~~~math
q_m\le(m/n)^{3/2}q_n+O(n^{3/2-c}),
~~~

and geometric-window summability plus exact landing forces
$q_n/n^{3/2}$, hence $M_n/n^{3/2}$, to converge.

The leading structured sufficient package is cap-inflated project-row coarea.
First prove a child-ground box witness

~~~math
\min_{S,\ y\text{ child ground}}
\min_{w\in\{\pm1\}^{S^c}}
\lVert A[:,S]y+A[:,S^c]w\rVert_2^2
\le R_0=O(n^{9/4-c}).
~~~

Then prove that some cap $R_0\le R\le C(R_0+n^2)$ and a pure or mixed
common-core kernel with
$\kappa=(\lambda_1-\lambda_2)/(1-\lambda_1)\ge\kappa_0>0$ satisfies

~~~math
\Pr\{f_z(T)=1\mid C,f_z(S_0)=f_z(S_1)=1\}\ge\lambda_1.
~~~

The box witness remains inside the inflated cap, so the conditional
denominator is positive. The displayed inequality is exactly non-strict
coarea by (10.1235); slice FKN then gives a constant-degree project-row
center and hence (10.795).

## Evidence, obstructions, and falsification criteria

- The exact completion CDF has the two-moment floor $Z_t\ge\mathcal H_t$
  from (10.1243)--(10.1246). Its positive branch relaxes local deficit; its
  negative branch requires near-saturation of conditional Parseval. Either
  saved population proves convergence.
- The two-moment implementation would be falsified only by a relevant
  unbounded exact-minimizer/window family failing the saved bound for every
  admissible tolerance and usable $c$; only analogous uniform failure for
  $Z_t$ falsifies the full route. At $A_9,m=8$, $88.66\%$ of favorable mass
  has negative margin while two moments capture only $6.97\%$ of $Z_0$.
- Generic concentration supplies upper tails, not the needed lower bound.
  A robust project-codimension cylinder already implies the desired
  recurrence by (10.1233), so that shortcut is circular.
- Every fixed-density project-row complement column pays
  $e^{-\Omega(n^{3/4})}$ by (10.1216). Complement mass, high slack, scalar
  optimization, arithmetic retention, and omitted-block excess remain retired.
- Coarea mass is partial-completion box discrepancy. Principal shortfall is
  rigid; finite examples instead show strong outside-column cancellation.
- Coarea is exactly triple retention. An exact $n=10,m=6$ minimizer has the
  globally optimal box row $10$ but fails every positive-core scale at cap
  $10$, including positive-core-only mixtures. Cap inflation repairs two
  scales, while $\frac34K_0+\frac14K_1$ passes at cap $10$. In general its
  $K_0$ contribution is $r_C\le\max a_z$, so constant-$\kappa$ repair may
  be degree-circular. An unbounded exact family failing box mass and every repaired kernel/cap would falsify coarea.
- The canonical Bellman point $p_e=(1-w_e)/4$ is exactly a mixture of
  switchings $A\mapsto D_UAD_U$. Equation (10.1240) is expected slack along
  that orbit, and a ground migrates explicitly with every outcome while its
  bad child restriction stays fixed. Canonical rounding, complete-family
  anti-migration, and state-dependent pressure based only on this point are retired. Tight
  decomposition remains open only through a non-switching certificate or
  higher-order incompatibility.
- Strict pressure, high replicas, local harmonic Poincaré, monotone deletion,
  nested chains, greedoids, ground-face-only pressure, and constant shortfall
  without a square-root cavity reward remain retired.

## Most recent blank-slate abstraction audit

The Wave 45 audit began from the original problem before route comparison.
Its candidates and judgments are agent-authored hypotheses, not user
directives. This early decisive-result refresh does not regenerate it;
Wave 55 will.

1. **Finite-temperature zero-sum interpolation.** Balanced $o(n)$ soft-minimax
   additivity would imply convergence, but temperature mismatch and overlap
   geometry still lack a signing-specific bridge.
2. **Summable cavity derivative.** Extension error
   $O(n^{1/2-\delta})$ would make normalized positive variation summable.
   This remains dormant until a square-root deletion reward is proved.
3. **Second-order signed-kernel compactness.** Continuity plus all-order
   recovery would prove convergence; existing compactness models lack the
   recovery operation.

## Ranked routes and Wave 51 direction

1. **Annealed bare completion tail.** Prove saved positive-margin or
   near-Parseval mass in (10.1246), or a higher-moment negative-tail bound.
2. **Box/coarea with controlled repair.** Prove retention after
   $R_0\le R\le C(R_0+n^2)$; separately test whether $K_0$ mixing avoids
   the degree-circularity tradeoff. Pure and positive-core-only same-cap kernels are finitely false.
3. **Annealed profile/spectral dichotomy.** Separate high-operator-norm
   profiles, where reverse chaos may supply the tail, from flat profiles,
   where a minimizer-specific restriction theorem is required.
4. **Tight decomposition, dormant.** Reopen only with a non-switching global
   certificate not encoded by (10.1240); tropical affine endpoints alone are
   insufficient.
5. **Other global routes, dormant.** Reopen only with a testable composition,
   recovery, or square-root cavity theorem at the required scale.

Wave 51 must independently rank ten ideas and assign three falsifiable attacks.
