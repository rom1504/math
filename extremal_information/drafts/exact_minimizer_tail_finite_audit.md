# Finite audit of the exact-minimizer tail proposal

**Status:** FINITE EXPERIMENTAL EVIDENCE ONLY / NO ASYMPTOTIC THEOREM

**Verdict:** `L_tail` and its stronger spectral sufficient condition survive
this finite counterexample search, but neither is materially established.
The exact corpus through order fourteen has a positive measured upper-tail
rate at every predeclared threshold.  At the principal value `d_0=1/8`, the
smallest rate among the 61 available repository exact representatives is

```math
I_{\rm tail}=0.2143980062\ldots\quad\hbox{nats per vertex},
```

attained at order ten by a tail of 60 out of 512 projective spins.  The
largest observed exact-minimizer spectral ratio is

```math
\|A\|_{2\to2}/\sqrt n=1.3744068037\ldots
```

at order eleven.

These are finite minima and maxima, not candidate universal constants.
Exact representatives are not exhaustive above order eight, and orders
12--14 are strongly construction-biased.  The audit therefore produces no
valid lower bound on the asymptotic `kappa` in SB.8 and no valid uniform
constant in SB.9.

## 1. Reproducibility and frozen input

The preregistered protocol is

[`../experiments/exact_minimizer_tail_finite_protocol.md`](../experiments/exact_minimizer_tail_finite_protocol.md),

and the deterministic program and result are

[`../experiments/audit_exact_minimizer_tail.py`](../experiments/audit_exact_minimizer_tail.py)
and
[`../experiments/exact_minimizer_tail_finite_results.json`](../experiments/exact_minimizer_tail_finite_results.json).

Their frozen hashes are

```text
protocol
  bb032b565257b245060540a5c3c94ab5ca0b92295e3e2c78c7607d8de050b27f

script
  3d1725eef58f5b46e98788014903d0d4a151334b80d2e4c0ac032fa044a8a435

result
  e3ecc65727f73e8d56beceddd6b976c3879d2c762cb58d69e9cf189bd054b1b1
```

The low-cap corpus was not regenerated or enlarged.  It is the already
frozen

```text
nearmin_blind_structural_results.json
sha256 2c086cf7523ead804942948e800c6231eac33d954e049b5aa113c9fb0cca47a5
```

The protocol froze the mathematical target at source SHA
`87f517bdf945f71e2c45a82bcdb75bee4c55d5498b067a86845339b0b7a5c5ea`.
The source was subsequently scope-repaired to SHA
`37cc9807b847ab1b2935dd5f5d1f9c2dfd9a0630e536c37d761154e49c3bcd18`;
SB.8 and SB.9 themselves are unchanged.  The repair correctly limits the
consequence of `L_tail` to conditional boundary-response profiles until a
further low-cap selector is supplied.
After independent proof audit, a status-only update produced final source
SHA `5278e6cb96a3a554141fe52cfd31dbb1ca38cf7b1260a33a554c116bf6074e8f`.

Reproduce the audit with

```bash
.venv/bin/python \
  extremal_information/experiments/audit_exact_minimizer_tail.py
```

The run takes about nine seconds in the current environment.  It reproduced
every frozen uniform-random and cyclic-control cap histogram and the complete
root-gauged minimizer counts `2,6,12,12,3240` at orders `3,...,7`.

## 2. Quantity audited and absolute-cap orientation

For

```math
H_A(x)=\sum_{i<j}A_{ij}x_ix_j,
\qquad Q(A)=\max_x|H_A(x)|,
```

put

```math
P_+=\max_xH_A(x),\qquad P_-=-\min_xH_A(x).
```

The `L_tail` statement first chooses a global matrix sign for which the
positive maximum equals `Q(A)`.  The audit does the same.  If only one of
`P_+,P_-` equals `Q(A)`, that orientation is mandatory.  If both do, the
audit deliberately chooses, separately at each threshold, the orientation
with the **larger** upper tail.  This is the adverse choice for SB.8 and is
also the correct uniform interpretation because both `A` and `-A` belong to
the exact-minimizer class.

For the selected orientation define

```math
p_A(d_0)=2^{-(n-1)}
 #\{x\text{ projective}:Q(A)-H_{\pm A}(x)<d_0n^{3/2}\}
```

and

```math
I_A(d_0)=-{1\over n}\log p_A(d_0).               \tag{FT.1}
```

Because quadratic energies are antipodally even, `p_A` is also the full
Boolean-cube density.  Hence the full upper-level count is exactly

```math
\exp\{(\log2-I_A(d_0))n\}.                       \tag{FT.2}
```

Thus `I_A` is the finite empirical counterpart of the `kappa` in SB.8,
with natural logarithms.  The predeclared grid was

```text
d_0 in {1/64, 1/32, 1/16, 1/8}.
```

The strict shell inequality was tested by exact integer squaring, not a
floating comparison.  This matters at the small orders because energy
deficits lie on a coarse parity lattice.

Orientation is not cosmetic.  For example, the worst exact spectral witness
at order eleven has positive maximum 15 but negative maximum 17.  Testing
only the original positive tail would not audit SB.8 at all.  Moreover many
two-sided-cap matrices have unequal positive and negative shell counts; the
adverse tie rule is active rather than formal.

## 3. Frozen strata

After byte-level matrix deduplication, the nonexhaustive strata contain:

| stratum | matrices | role |
|---|---:|---|
| repository exact representatives | 61 | exact `Q=M_n`; incomplete orbit coverage above order 8 |
| repository one-step-near representatives | 9 | exact recomputed `Q=M_n+2` |
| cap-constrained adversarial samples | 88 | distinct frozen samples; search never left `Q<=M_n+2` |
| independent greedy low-cap samples | 28 | heuristic discovery, exact recomputed cap |
| unconditioned uniform controls | 521 distinct of 576 frozen draws | random baseline |
| uniform draws conditioned to be low-cap | 128 | selection-bias control |
| cyclic-distance controls | 378 | deterministic structured baseline |

The full 48-draw-per-order random controls were reconstructed from the
frozen seed only because their matrices were not individually stored in the
original JSON.  Their cap histogram at every order exactly matches the
frozen summary.  No descent, rejection, or minimizer search was performed.

## 4. Exhaustive root-gauged orders three through seven

Every root-gauged signing was enumerated through order seven.  The table
reports the range of `I_A(1/8)` and of the spectral ratio
`S(A)=||A||_(2 to2)/sqrt n` in the exact and one-step cap strata.

| `n` | class | count | min/median/max `I_A(1/8)` | min/median/max `S(A)` |
|---:|---|---:|---:|---:|
| 3 | exact | 2 | .4621 / .4621 / .4621 | 1.1547 / 1.1547 / 1.1547 |
| 3 | one-step | 0 | -- | -- |
| 4 | exact | 6 | .5199 / .5199 / .5199 | 1.1180 / 1.1180 / 1.1180 |
| 4 | one-step | 2 | .5199 / .5199 / .5199 | 1.5000 / 1.5000 / 1.5000 |
| 5 | exact | 12 | .2326 / .2326 / .2326 | 1.0000 / 1.0000 / 1.0000 |
| 5 | one-step | 30 | .4159 / .4159 / .4159 | 1.3416 / 1.3416 / 1.3416 |
| 6 | exact | 12 | .2790 / .2790 / .2790 | .9129 / .9129 / .9129 |
| 6 | one-step | 180 | .4621 / .4621 / .4621 | 1.2247 / 1.2247 / 1.2247 |
| 7 | exact | 3240 | .3161 / .3961 / .4372 | 1.1339 / 1.3461 / 1.3461 |
| 7 | one-step | 16884 | .4372 / .4951 / .5941 | 1.3324 / 1.4470 / 1.5496 |

At these orders all four predeclared thresholds usually select the same
energy lattice layer.  The table is therefore a correctness and
counterexample check, not evidence of a rate curve.  In particular, the
order-five exact tail density `5/16` gives the pooled exact minimum
`I=.23263...` at all four thresholds, yet says nothing about large-order
behavior.

## 5. Available exact representatives at orders eight through fourteen

For each threshold the next table displays the largest selected projective
tail found at that order and its corresponding rate in parentheses.  The
number of available byte-distinct exact representatives and the largest
spectral ratio are also shown.

| `n` | exact reps | `S_max` | `d_0=1/64` | `1/32` | `1/16` | `1/8` |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 8 | 1.2748 | 4/128 (.4332) | 4/128 (.4332) | 4/128 (.4332) | 14/128 (.2766) |
| 9 | 6 | 1.3725 | 15/256 (.3152) | 15/256 (.3152) | 15/256 (.3152) | 15/256 (.3152) |
| 10 | 9 | 1.3038 | 20/512 (.3243) | 20/512 (.3243) | 20/512 (.3243) | 60/512 (.2144) |
| 11 | 9 | 1.3744 | 24/1024 (.3412) | 24/1024 (.3412) | 24/1024 (.3412) | 88/1024 (.2231) |
| 12 | 8 | 1.1355 | 10/2048 (.4435) | 10/2048 (.4435) | 46/2048 (.3163) | 112/2048 (.2422) |
| 13 | 2 | 1.0000 | 39/4096 (.3580) | 39/4096 (.3580) | 39/4096 (.3580) | 195/4096 (.2342) |
| 14 | 3 | .9636 | 78/8192 (.3324) | 78/8192 (.3324) | 260/8192 (.2464) | 260/8192 (.2464) |

No rate in this table is zero or visibly decaying over the short range.
That is the strongest finite fact favorable to `L_tail`.  It remains weak
evidence because there are only one or two known structural families behind
many of the last rows, and no exhaustive search above order eight.

The worst recorded exact `d_0=1/8` witness has hash

```text
8a7654f72476d29aa4e24d8131078a41a8b8cefc884e8b23e5e8a89b7bfb22e1
```

at order ten.  It has two-sided cap 13 and 60 selected projective states,
giving density `0.1171875` and rate `0.2143980062817407` nats per vertex.
The full matrix is retained in the machine result.

## 6. Counterexample-first comparison at `d_0=1/8`

Each entry below is

```text
number of matrices : minimum rate / median rate ; maximum spectral ratio.
```

`E`, `N`, `A`, `R`, and `C` denote repository exact, repository one-step,
adversarial low-cap, unconditioned random, and cyclic structured strata.

| `n` | `E` | `N` | `A` | `R` | `C` |
|---:|---:|---:|---:|---:|---:|
| 8 | 8:.277/.277;1.275 | 1:.363/.363;1.458 | 40:.277/.363;1.454 | 48:.332/.469;1.832 | 16:.296/.451;2.475 |
| 9 | 6:.315/.323;1.373 | -- | 12:.308/.417;1.488 | 48:.385/.494;1.863 | 16:.372/.555;2.667 |
| 10 | 9:.214/.214;1.304 | 2:.347/.396;1.244 | 6:.214/.324;1.419 | 48:.384/.514;1.962 | 32:.347/.624;2.846 |
| 11 | 9:.223/.223;1.374 | 1:.280/.280;1.518 | 13:.223/.302;1.483 | 48:.304/.421;1.905 | 32:.349/.349;3.015 |
| 12 | 8:.242/.242;1.135 | 4:.255/.309;1.336 | 7:.242/.300;1.352 | 48:.358/.486;1.915 | 64:.395/.635;3.175 |
| 13 | 2:.234/.234;1.000 | 1:.265/.265;1.374 | 5:.234/.285;1.238 | 48:.392/.463;2.128 | 64:.234/.389;3.328 |
| 14 | 3:.246/.246;.964 | -- | 5:.246/.246;1.193 | 48:.351/.446;1.901 | 128:.340/.644;3.474 |

The adverse finding is not a zero exponent but a systematic direction:
available exact minimizers generally have **smaller** tail exponents, hence
denser near-top sets, than generic random controls.  This matches the prior
soft-shell observations and warns against importing a random-energy tail
constant.  The desired statement only asks that the exponent remain
positive, so the comparison does not falsify it.

The adversarial-walk minima frequently coincide with an exact seed.  They
do not produce a larger finite upper tail than the worst exact examples at
orders 10--14.  This is a meaningful failed falsification attempt, but the
walk explores only frozen local low-cap components.

## 7. Spectral sufficient condition

The observed pooled maxima of `||A||op/sqrt(n)` are:

| stratum | maximum |
|---|---:|
| repository exact | 1.374407 |
| repository one-step near | 1.518486 |
| cap-constrained adversarial | 1.488034 |
| greedy low-cap | 1.522902 |
| unconditioned random | 2.127620 |
| cyclic structured | 3.474396 |

When the complete exhaustive one-step populations through order seven are
combined with the repository list, the true finite one-step maximum is
`1.5496363663...` at order seven, hash
`8d590a27dc380730f8a0d4d037e1d345e006f346a44fe4774183ebd7d340a620`.
This explicit combined witness block is retained in the machine result; the
`1.518486` entry above is intentionally scoped to the nine repository
one-step representatives.

The worst exact witness is the order-eleven matrix with hash

```text
eca03a38db2d890177f0682ef72312bc886c40a45b776a653297c4d218cd55f0.
```

It has `Q=17`, positive maximum 15, negative maximum 17, and spectral ratio
`1.3744068037446782`.  Correct orientation selects the negative landscape;
at `d_0=1/8` it has 88 of 1024 projective states and rate
`0.2231031810110224`.

Nothing here proves SB.9.  Any finite corpus has some finite maximum, and
the exact-minimizer inventory is sparse at the largest orders.  Conversely,
the cyclic family includes the all-positive matrix, whose spectral ratio is
`(n-1)/sqrt(n)` and grows, while its extreme upper tail is very sparse.
This gives a scalable calibration that SB.9 is sufficient rather than
necessary for a positive tail rate.  It is not a counterexample to SB.9 on
exact minimizers because the all-positive matrix is far from minimizing.

## 8. What was and was not falsified

### Survived this finite audit

- Every available exact representative has positive `I_A(d_0)` on the
  entire predeclared grid.
- Exhaustive exact populations through order seven have the same property.
- Neither the available one-step corpus nor the cap-constrained adversarial
  corpus exhibits a zero-rate finite tail.
- Available exact representatives have spectral ratios below `1.375`.

### Not established

- A finite positive minimum is not a uniform asymptotic `kappa`.
- No order above eight has exhaustive orbit coverage in this audit.
- The spectral ratios do not establish a universal operator bound.
- A positive `L_tail` result would yield the conditional response packing of
  Theorem 21.8, not by itself a scalar all-spins-free parent-cap packing.
- Neither candidate compares different orders or advances convergence of
  `M_n/n^(3/2)`.

### Strongest live falsifier

The most direct falsifier remains a certified sequence of exact minimizers
`A_n` for which, for every fixed `d_0>0`,

```math
-{1\over n}\log
\Pr\{Q(A_n)-H_{\pm A_n}(X)<d_0n^{3/2}\}\longrightarrow0
```

in at least one admissible global orientation.  A spectral falsifier would
be exact minimizers with `||A_n||op/sqrt(n) -> infinity`; this would kill the
convenient sufficient condition but would not by itself kill `L_tail`.

## 9. Scoped research judgment

The data justify keeping `L_tail` as a precise, falsifiable exact-minimizer
target.  They do not justify promoting it above the other live route merely
because the observed rates look flat.  Its advantages are mathematical:

1. it is strictly weaker than reconstructing the Boolean landscape;
2. Theorem 21.8 consumes exactly this one upper-level count; and
3. SB.9 gives a recognizable sufficient condition.

Its disadvantages are equally clear:

1. exact minimizers have unusually dense near-top shells, so random-control
   heuristics point in the wrong quantitative direction;
2. the largest-order inventory is thin and structured; and
3. even a proof leaves a separate low-cap selector between boundary-profile
   separation and scalar parent-cap separation.

The best next finite discriminator would be a certified, symmetry-diverse
set of exact representatives beyond order fourteen, evaluated with this
unchanged threshold grid and orientation rule.  The best next theorem is
not a fit to the displayed constants: it is either a uniform
`||A||op=O(sqrt n)` consequence of exact minimality, or a direct
minimality-to-upper-tail inequality which bypasses spectral flatness.
