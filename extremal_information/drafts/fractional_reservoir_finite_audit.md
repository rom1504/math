# Finite audit of shell-wide fractional reservoirs

**Status:** residual-checked finite computation; not an asymptotic theorem.

This note audits the fractional-reservoir observable from
[`fractional_reservoir_localized_flip.md`](fractional_reservoir_localized_flip.md)
on the frozen low-cap corpus.  The selection and LP protocol were written
before solving any instance:

- [`fractional_reservoir_finite_protocol.md`](../experiments/fractional_reservoir_finite_protocol.md),
- [`fractional_reservoir_finite_audit.py`](../experiments/fractional_reservoir_finite_audit.py),
- [`fractional_reservoir_finite_results.json`](../experiments/fractional_reservoir_finite_results.json).

Run with

```bash
source .venv/bin/activate
python extremal_information/experiments/fractional_reservoir_finite_audit.py
```

The input, protocol, and script SHA-256 hashes are embedded in the result.
The frozen run selected 64 distinct matrix hashes and solved 128 LPs: one
active-shell and one deficit-2-shell problem for each matrix.

## 1. Observable

For each projective cut word, orient its edge word `z` so that
`<a,z>>=0`.  If `S_m` is the shell with oriented response at least `m`, solve

```math
W_*(A,S_m)=min\left\{
\sum_e w_e:
0\le w_e\le1,\quad
\sum_e w_ea_ez_e\ge m\ \text{for every }z\in S_m
\right\}.
```

The reported ratio is

```math
C_{\rm inst}=W_*(A,S_m)/m.
```

The active shell uses `m=Q(A)`; the deficit-2 shell uses `m=Q(A)-2`.
Full edge weight is always feasible, so this asks whether the complete finite
shell admits a substantially cheaper simultaneous fractional reservoir.

## 2. Corpus and validation

The 64 matrices consist of:

| stratum | matrices |
|---|---:|
| certified/exhaustive exact representatives | 21 |
| repository one-step-near representatives | 9 |
| independently generated heuristic low-cap | 6 |
| uniform-random draws that happened to be low-cap | 12 |
| cyclic/structured controls | 16 |

Orders range from 3 through 14.  Distinct matrix hashes at the larger orders
are not asserted to be distinct switching-permutation orbits.

Every LP was solved independently by HiGHS dual simplex and interior point.
Across all 128 records, the maximum cross-solver objective discrepancy was
`4.15e-12`, maximum constraint violation `3.07e-11`, maximum box violation
`5.77e-12`, and maximum reconstructed primal-dual gap `4.55e-13`.

## 3. Main finite result: broad shells often force full edge mass

Across all strata:

| shell | `C_inst` range | median | mean | numerical full-edge optima |
|---|---:|---:|---:|---:|
| active | `1`--`4.333333` | `2.333333` | `2.323287` | 36/64 |
| deficit 2 | `1`--`5` | `3` | `2.942564` | 31/64 |

The effect is stronger on the exact-minimizer stratum:

| shell | `C_inst` range | median | mean | numerical full-edge optima | no literal common-correct edge |
|---|---:|---:|---:|---:|---:|
| active | `1`--`4.333333` | `3` | `2.814544` | 13/21 | 17/21 |
| deficit 2 | `2`--`5` | `3.666667` | `3.568767` | 15/21 | 21/21 |

Thus disappearance of the literal common-correct reservoir does not make the
fractional LP infeasible, but the whole-shell replacement is frequently not
cheap: it often retains every physical edge.

For the two selected order-14 exact representatives, the active shell has 156
oriented words and

```math
W_*=91={14\choose2},
\qquad C_{\rm inst}={91\over21}={13\over3}.
```

The deficit-2 shell has 520 words and again `W_*=91`, now giving
`C_inst=91/19=4.789473...`.  At order 13, the active shell forces all 78
edges, whereas the deficit-2 shell is supported by the diffuse uniform weight
`0.9` on all 78 edges, giving `W_*=70.2` and `C_inst=3.9`.

The complete exact-representative table is:

| `n` | reps | active shell sizes | active `C_inst` | active full | deficit-2 shell sizes | deficit-2 `C_inst` | deficit-2 full |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 1 | 1 | 1.000 | 1/1 | 4 | 3.000 | 1/1 |
| 4 | 1 | 2 | 1.000 | 0/1 | 6 | 2.000 | 0/1 |
| 5 | 1 | 10 | 2.500 | 1/1 | 10 | 2.500 | 0/1 |
| 6 | 1 | 12 | 3.000 | 1/1 | 32 | 5.000 | 1/1 |
| 7 | 3 | 3--7 | 1.667--2.333 | 1/3 | 13--15 | 3.000 | 3/3 |
| 8 | 2 | 8 | 2.000 | 0/2 | 28 | 3.500 | 2/2 |
| 9 | 2 | 21--25 | 3.000 | 2/2 | 21--25 | 2.600--2.700 | 0/2 |
| 10 | 2 | 40 | 3.462 | 2/2 | 120 | 4.091 | 2/2 |
| 11 | 2 | 5--24 | 2.647--3.235 | 1/2 | 41--53 | 3.667 | 2/2 |
| 12 | 2 | 20 | 3.333 | 0/2 | 92 | 4.125 | 2/2 |
| 13 | 2 | 78 | 3.900 | 2/2 | 78 | 3.900 | 0/2 |
| 14 | 2 | 156 | 4.333 | 2/2 | 520 | 4.789 | 2/2 |

No monotone or asymptotic claim is made from this short table.  Descriptively,
however, `C_inst` is strongly associated with shell size in the frozen corpus:
Spearman correlations are `0.935` on the active shell and `0.858` on the
deficit-2 shell.  These are exploratory correlations on selected finite data,
not inferential statistics for an ensemble.

### Does it detect near-minimality?

There is a finite within-order signal, but not a clean classifier.  At each of
the ten orders 4--13 containing both cap-delta-0 and cap-delta-2 records, the
mean active-shell `C_inst` was larger for delta 0 in 9/10 orders and tied in
one; the mean paired difference was `0.719`.  For the deficit-2 observable it
was larger in all 10/10 orders, with mean paired difference `1.125`.  Mean
support fraction moved in the same direction in only 6/10 active comparisons
and 8/10 deficit-2 comparisons.  Thus `C_inst` is the clearer of the two
finite signals.

This signal should not be overread.  Exact representatives and non-exact-
provenance controls having the same cap delta had essentially identical
within-order values (many are structurally equivalent copies), ranges overlap
substantially, and order and shell size are strong confounders.  The observable
appears to measure breadth/compatibility of the exposed shell more directly
than it measures the provenance label “exact minimizer.”

## 4. An exact finite full-mass certificate

Some numerical full-mass conclusions have a simple exact certificate.  Let
`R` be the shell-by-edge response matrix, let `K` be the shell size, and let
`E` be the edge count.  If

```math
{1\over K}\sum_{z\in S_m}R_{z,e}\le {m\over E}
\qquad\text{for every edge }e,                 \tag{FA.1}
```

then averaging the constraints of any feasible solution gives

```math
m\le {1\over K}\sum_z\sum_eR_{z,e}w_e
\le {m\over E}\sum_ew_e.
```

Hence `W>=E`; since `w_e<=1`, necessarily `W=E` and every weight is one.
Condition (FA.1) was checked with exact integer arithmetic through
`E sum_z R_{z,e} <= Km`.

It certifies 18 active-shell records, including 7 of the 21 exact records:
orders 3, 5, 6, both selected order-13 representatives, and both selected
order-14 representatives.  The other numerical full-mass cases may have
nonuniform dual certificates; this audit does not claim exact rational
certificates for them.

## 5. Diffuseness and integrality

The LP was designed to allow diffuse fractional weights, but a HiGHS extreme
optimum was integral in 62/64 active-shell instances and 47/64 deficit-2
instances.  This does **not** prove integrality of the polytope or uniqueness
of a sparse optimum.  The full-mass instances are unique because the box
constraint and total weight `E` force `w_e=1` for every edge.

At the opposite extreme, `C_inst=1` occurred in 10 active and 4 deficit-2
records.  In every such record the literal common-correct edge set already
had at least `m` edges, so fractionalization supplied no new reduction.  The
interesting middle cases exist—for example, the order-13 deficit-2 uniform
`0.9` solution—but they are not the dominant exact-minimizer pattern here.

## 6. Research judgment

This experiment gives a useful **ceiling**, not a route falsifier.

- The finite-phase theorem FR.1 fixes an anchor count `K` first.  Its constant
  is allowed to depend on `K`.  The present LP instead inserts the *entire*
  active or deficit-2 shell, with `K` as large as 520.  Full-mass solutions do
  not contradict FR.1 or the growing-packing argument.
- The data argue against upgrading FR.1 to a shell-uniform `O(M)` reservoir
  using only simultaneous positive response constraints.  On the selected
  order-14 exact matrices, the required mass is all `E=91` edges, already on
  the active shell.
- Consequently the promising use of fractional reservoirs remains
  **sequential finite-anchor extension**.  Treating the complete near-ground
  shell as one jointly supported phase loses the desired localization and can
  cost the full landscape interface.
- The next discriminating mathematical question is whether the sequential
  construction can choose a growing anchor set with extra geometry that keeps
  its optimal fractional reservoir at `O(M)`, rather than whether arbitrary
  broad shells enjoy that bound.  The present computation supplies explicit
  finite targets and dual certificates for attempts at such a strengthening.

All conclusions above remain finite evidence unless an exact certificate is
explicitly stated.

The strongest rigorous finite falsifier is therefore narrow: **fractional
weights do not automatically turn the full active shell of an exact minimizer
into a proper sub-interface.**  The exact uniform-average certificate forces
`w_e=1` on both selected order-13 and both selected order-14 exact records.
It does not by itself falsify an order-uniform asymptotic theorem; that would
require a scalable family carrying the same certificate.
