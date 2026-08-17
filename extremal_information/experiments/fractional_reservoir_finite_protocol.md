# Frozen protocol: finite fractional-reservoir audit

**Status:** protocol frozen before solving any of the LPs described below.
This is a finite experimental audit, not a theorem or asymptotic claim.

## Question

For a hollow sign matrix `A`, let `Q=Q(A)` and orient every projective cut
word `z` so that `<a,z>>=0`.  At shell threshold `m`, solve

```math
\min_{0\le w_e\le1}\sum_e w_e
\quad\text{subject to}\quad
\sum_e w_e a_ez_e\ge m
\quad\text{for every oriented shell word }z.       \tag{P}
```

The two predeclared shells are:

- `active`: `m=Q` and `|<a,z>|=Q`;
- `deficit_2`: `m=Q-2` and `|<a,z>|>=Q-2`.

The primary statistic is `C_inst=(sum_e w_e)/m`.  We also record the support,
number of fractional coordinates, maximum normalized atom, inverse-Herfindahl
effective support, and constraint/bound residuals.  These observables assess
how cheaply the *whole finite shell* can be supported simultaneously.  They
do not certify the asymptotic constants in the fractional-reservoir theorem.

## Frozen corpus rule

Input is
`extremal_information/experiments/nearmin_blind_structural_results.json`.
Matrices are selected without looking at the output of (P):

1. **exact:** every authoritative orbit representative at orders 3--8, and
   the first two distinct repository exact representatives by matrix SHA-256
   at each order 9--14;
2. **one_step_near:** every distinct repository one-step-near representative;
3. **heuristic_low_cap:** the first distinct independently generated greedy
   witness by SHA-256 in each `(n,cap_delta)` stratum;
4. **random_low_cap:** the first distinct uniform-random low-cap draw by
   SHA-256 in each `(n,cap_delta)` stratum;
5. **structured_control:** the first distinct cyclic-distance low-cap control
   by SHA-256 in each `(n,cap_delta)` stratum, followed by every distinct
   `control_extremes` matrix.

Deduplication is global in the stated priority order.  Thus a matrix already
selected as exact is not counted again as a control.  Within a source and
stratum, SHA-256 order—not a structural statistic—chooses the representative.
No matrix will be added or removed after inspecting its LP result.

## Numerical validation

Each LP is solved independently with HiGHS dual simplex and HiGHS interior
point through SciPy.  A record is accepted only if:

- both solvers report success;
- their objectives agree to `1e-7 * max(1, objective)`;
- maximum constraint and box violation is at most `1e-7`;
- the dual objective reconstructed from HiGHS marginals agrees with the
  primal objective to `1e-7 * max(1, objective)`.

Entries below `1e-9` are treated as zero only for reporting support.  Raw
minimum/maximum coordinates and residuals remain in the JSON.  Aggregate
summaries are computed by provenance, order, shell type, and cap delta.

## Stopping rule

Stop rather than enlarge the corpus if either shell enumeration exceeds the
available projective state space through order 14 or the frozen batch becomes
computationally unwieldy.  A surprising result may motivate a later protocol,
but it will not alter this one.
