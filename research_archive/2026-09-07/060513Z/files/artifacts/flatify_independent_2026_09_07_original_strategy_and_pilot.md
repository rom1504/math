# Independent original-problem route freeze and ternary construction pilot

Date: 2026-09-07. Status: research record, not a convergence proof.

The original problem was supplied without historical route vocabulary. Before
consulting ACTIVE_STATE, STEERING, or the ledger, I froze the augmented
cut-code covering-radius mapping and sought a code-composition theorem.
Primary searches found Sole--Zaslavsky's cocycle-code paper
(https://people.math.binghamton.edu/zaslav/Tpapers/cas.sidma1994.pdf), but the
ordinary cocycle code needs the all-one augmentation for this problem.
The archived audit already contains the relevant normal-code/shortening
obstructions. No new applicable code theorem was found. The proposed
composition inequality itself was equivalent to the missing all-order
recovery obligation, and was abandoned as a new route.

A second independent derivation tried powered block aggregation. Exact
endpoint-width superadditivity does not supply the required 2/3-power
superadditivity. The archive contains scalable random-restriction and
near-minimizer planting obstructions; no such recurrence is claimed.

The concrete surviving bounded task is the already-realized all-order
ternary ensemble, whose rigorous scalar bound is

```math
\limsup M_n/n^{3/2}
 \le {t+p\log 2+E_t(\nu_p)\over 2t\sqrt p},
\qquad
\nu_p=(1-p)\delta_0+{p\over2}(\delta_{-1/\sqrt p}+\delta_{1/\sqrt p}).
```

I reconstructed the archived ternary reproduction-support reduction before
using its exact two-variable variational formula. This is an optimization
of an actual ensemble, not an identification of its optimum with M_n.

## Numerical pilot and a detected optimizer failure

At p=24/25,t=24/5 the pilot indicates Gaussian phase and the scalar value
0.4936690966983393. At p=24/25,t=97/20 an informative branch exceeds the
Gaussian branch; the actual pilot scalar value is about 0.493607276804.
The previous t=4,p=31/32 scalar point is not stationary in t within its
Gaussian phase.

**Important numerical failure.** A single differential-evolution run at
p=.96,t=4.85 returned the near-endpoint informative branch, with value
-0.824293760983 and a spuriously low scalar bound 0.493594019441.
A separate local optimizer found the other branch at lambda about
3.1885677,z about .95999445, with larger value -0.824167762739.
Therefore the smaller pilot number is NOT an upper bound. The rigorous
verifier covers every point of the precision/reproduction rectangle and
does not trust numerical critical-point completeness.

The numerical best neighborhood has three almost tied branches: Gaussian,
interior ternary reproduction, and near-endpoint ternary reproduction.
The pilot triple contact is p about .95977855,t about 4.86035103, giving
about .493597061588. This is not a certified optimum, and the task does not
depend on identifying it.

The clean fixed rational certificate target is p=24/25,t=97/20 and
E_t(nu_p)<=-10302/12500=-.82416. If verified, this gives a bound below
.49361. Its verification and independent audit are separate obligations.

## Artifacts

- `tmp/flatify_independent_2026_09_07_ternary_pilot.py`: numerical pilot,
  with fixed seed 73207 and explicit record of all coarse grid points.
- `computations/results/flatify_independent_2026_09_07_ternary_pilot.json`:
  numerical coarse-grid output; some individual maxima may be missed.
- `computations/flatify_independent_2026_09_07_ternary_interval.py`:
  directed interval verifier, with exact rational subdivision endpoints.

The interval script's new mean-value bound differentiates the exact
optimized mixture-weight envelope. It encloses the clipped optimum weight
throughout each rectangle, then pays the supremum of both coordinate
derivatives times the corresponding half-width. Every accepted rectangle
uses directed intervals; floating point never prunes a rectangle.
