# Independent algebra audit: uniform scalar-hierarchy escape

Date: 2026-09-05. Read the complete proof in
`fresh_uniform_scalar_hierarchy_escape_2026_09_05.md`, with particular
attention to Sections 4--6. The stated strict escape from the scalar
central-mask supremum passes this audit. It is not original convergence.

1. For every mass in the compact interval `[9/50,5/8]`, the normalized
   central indicator has Hermite derivative energy infinity. Thus its
   nonnegative-coefficient noise kernel has `K(q)<q` for some `q<1`.
   Continuity at each such fixed `q`, the proved height barrier, and a
   finite mass subcover give a single positive residual lower bound.
   The resulting lower bound on `Var(V|W)` is uniform in the entire
   first-chaos unit sphere, not merely finite-height vectors.

2. Enlarging the feasible covariance triples to the displayed compact
   positive-semidefinite set is legitimate. The lower bound on the
   edge-linear coefficient uses only that covariance matrix and the
   central-mask geometry. In particular it remains valid on every
   possibly non-realizable added triple. The numerical endpoint check
   is strict: `sqrt(35/88)-5/8>1/200`, and concavity checks the other
   endpoint. Hence no added parameter can have all higher odd Hermite
   coefficients zero: its bounded conditional response would otherwise
   equal a nonzero unbounded linear function.

3. The coefficient maps are continuous even when the full triple is
   singular. A common Gaussian realization through a positive matrix
   square root gives almost-sure convergence; the only response
   discontinuities have zero probability under the nondegenerate
   `(V,W)` marginal. Fixed Gaussian Hermite moments provide domination.
   The finite open-cover argument therefore fixes both a finite degree
   bound and a strictly positive coefficient threshold.

4. The compact parameter range and the uniform determinant lower bound
   for `(V,W)` bound its density below on the fixed rectangle
   `1<=V<=2`, `|W|<=t/4`, and its reflected copy. The central thresholds
   are all below `.9`, so this is a uniformly feasible local-update
   region. The smoothing loss, local mass, and selected coefficient
   bounds consequently have constants independent of the target.

5. All-odd local normalization leaves variance at most one, and the
   odd Schur standard-deviation theorem supplies a uniform positive
   mean coefficient. The finite degree bound is fixed before `L`, the
   finite smooth approximant, and the matrix-size limit. Dependence of
   finite-matrix error constants on that degree or on the approximant
   is therefore harmless. The clipping radius proportional to
   `sqrt(log L)` produces the same uniform `c/log L` improvement.

6. Selecting a target within a fixed fraction of that improvement of
   the scalar supremum, then a finite smooth approximation, and only
   then taking the matrix limit does not interchange the tree-depth
   and matrix limits. Principal deletion loses `O(epsilon)` while the
   fixed-operator gain is order `1/log(1/epsilon)`. For some fixed small
   `epsilon` the latter dominates, yielding a strict gap above the
   entire scalar central-mask supremum.

This audit depends on the already independently checked height barrier,
all-odd weighted projection, and normalized localized-gain theorems.
It does not claim that the enlarged response family is complete.
