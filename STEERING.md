# Strategic steering

Evidence cutoff: Wave 55, ledger §10.108 (2026-07-30), through (10.1314)--(10.1338).
Early decisive-result refresh; next mandatory refresh: Wave 60, with a blank-slate abstraction audit.

## User-stated research objective

Determine whether the limit of `M_n/n^(3/2)` exists. The conjectural value `1/2` is not an additional user objective.

## Present judgment

No route proves convergence; the rigorous interval remains

~~~math
0.336493364431\ldots\le\liminf\frac{M_n}{n^{3/2}}
\le\limsup\frac{M_n}{n^{3/2}}\le\frac12.
~~~

The selected-cut project-row route is now the clear leader. Its strongest
implementation is a low-row box witness followed by full weighted spectral
excess. The uniform arbitrary-cut annealed completion route is no longer an
independent alternative: joint Hanson--Wright concentration proves that any
saved incidence under that product law already contains recurrence-strong
local-band mass. The K-profile route is decisively falsified.

## Leading route: selected low-row cut via box and full histogram

Fix a compact density window `p=m/n in (1/2,1)`, choose `0<c<1/4`, and put

~~~math
H=\lceil n^{3/4-c}\rceil,qquad T_n=n^{3/2-c}.
~~~

For an exact order-`n` minimizer, a selector `S`, and a child ground `y`, let

~~~math
\mathcal V(S,y)=\min_{w\in\{\pm1\}^{S^c}}
\lVert A[:,S]y+A[:,S^c]w\rVert_2^2.
~~~

The exact sufficient package is:

1. **Box:** for some child-ground incidence,

   ~~~math
   \mathcal V(S,y)\le R_0=O(n^{9/4-c}).
   ~~~

2. **Full weighted excess:** for some controlled
   `R_0<=R<=C(R_0+n^2)` and some nontrivial core `2<=ell<m`, the associated
   low-row center class satisfies

   ~~~math
   P_\ell(R)-\lambda_2(\ell)\ge e^{-O(H)}.
   ~~~

The extraction theorem (10.1261) then gives one row-`R` global cut favorable
on `e^{-O(H)}` selector mass. This is the bare fixed-cut tail (10.795), which
gives

~~~math
q_m\le p^{3/2}q_n+O(T_n).
~~~

Geometric-window summability and exact landing then force
`q_n/n^(3/2)`, hence `M_n/n^(3/2)`, to converge.

The current box attack uses exact minimality. If the internal field square
`I` is much larger than `n^(9/4-c)`, the Wave 55 attenuated-block theorem
extracts a field-excess edge set and forces a parent state of deficit
`o(T_n)` that reverses a fixed fraction of it. The exact missing bridge is an
aggregate child-fibre/parent-witness compatibility theorem: charge those
migrated states back to low-cross favorable incidences or use their
multiplicity to produce a box-canceling center.

The overlap attack must retain the complete self-loop-free intersection
histogram (10.1292). At a fixed linear core it must generate
`exp(n v_*(p,alpha)+o(n))` partners for a positive fraction of biased bases,
spread across the `Theta(sqrt(n))` transition-typical intersection levels.

## Decisive retirements and falsification criteria

- For the full uniform local-state/completion incidence,

  ~~~math
  Z_t\le\nu_m\{g\ge-p_2DT_n\}+2e^{-cDH}.
  ~~~

  Thus every `Z_t>=e^{-O(H)}` lower bound already proves the local
  recurrence after fixed-depth truncation. This does not apply to a selected
  global cut.

- K success at depth `D T_n` has mass at most `2e^{-cD^2H}`. Hence the
  diverging-depth lemmas (10.1277), (10.1285), and the K branch of (10.1297)
  are impossible. Fixed-depth K abundance is recurrence-local. Scalar
  delocalization, canonical heads, and head exchange cannot repair this.

- Generic cap, operator, trace, child optimality, and block algebra do not
  force the box: scalable nonminimal signings hide `Theta(n^(5/2))` internal
  curvature with zero cross and correction terms. Exact minimality now forces
  witness migration, but finite exact examples falsify pointwise common-
  witness compatibility.

- Every linear-core one-threshold histogram compression is impossible. Away
  from its tuned surface the shortfall is exponential; on the tuned surface
  it is `Theta(sqrt(n))`. Nonnegative core mixtures are convex repackagings
  and cannot beat their best component.

- The leading package would be falsified as an implementation by an
  unbounded exact-minimizer family on which either every box value exceeds
  `n^(9/4-c)` for every usable `c`, or every controlled cap and nontrivial
  core has spectral excess below `e^{-omega(H)}`. That would not disprove
  convergence or the bare fixed-cut lemma.

## Ranked alternatives

1. **Aggregate migration compatibility.** Bound how many excess-edge choices
   can share one escaping parent state and transfer the mass to a low-row child fibre.

2. **Full weighted-histogram inequality.** Use signing compatibility and row
   geometry to beat the degree-two baseline without threshold compression.

3. **Direct selected-cut construction.** Prove (10.795) through a nonuniform
   low-row cut law, near-ground face, or deterministic selection principle.

4. **Tight decomposition, dormant.** Reopen only with a non-switching global
   certificate that rounds all-bad pressure with state-dependent budgets.

5. **Other global routes, dormant.** Reopen only with a testable composition,
   recovery, approximate-subadditivity, compactness, or square-root cavity
   theorem at the required scale.

Uniform annealed completion, orientation/K abundance, scalar K profiles,
canonical-head rescue, one-threshold clusters, nonnegative core mixtures,
constant shortfall, complement incidence, robust cylinders, harmonic
pressure, nested chains, greedoids, and unweighted edge-cube covers are
retired in their recorded forms.

## Blank-slate audit status

The most recent blank-slate audit was Wave 45. Its finite-temperature
interpolation, summable cavity derivative, and signed-kernel compactness
candidates supplied no recovery theorem. This Wave 55 refresh is early and
decisive-result-driven, not a scheduled mandatory audit. Wave 60 must begin
from the original convergence problem and perform the next blank-slate audit
before route comparison.

All rankings and hypotheses above are the main agent's assessment, not user
directives. Wave 56 must independently rank ten ideas against this steering
document and the Wave 55 Updated frontier.
