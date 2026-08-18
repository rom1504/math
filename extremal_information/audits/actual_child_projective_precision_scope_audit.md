# Precision-scope audit of the projective synchronization data

## Verdict

The committed file
`computations/results/actual_child_projective_synchronization.json` contains
correct bridge calculations for the signings it enumerates, but **66 of its
102 records were not induced by optimizing children**.  The cause is a
precision-initialization bug in the driver, not in the analytic projective,
curvature, influence, or KL identities.

Accordingly, every statement which calls all 102 records "actual minimizing
child laws" is withdrawn.  Theorems 37.24, 37.27, and their analytic proofs
do not depend on this enumeration and remain valid.  The finite numerical
reports can be repaired either by regenerating them after the accompanying
driver patch or, for diagnostic purposes, by filtering the old records by
the high-precision winner SHA values.

## 1. Exact cause

The child selector uses

```python
tie_tolerance = mp.mpf(10) ** (-(mp.mp.dps - 20))
```

to reserve twenty guard digits.  The projective synchronization driver did
not initialize `mp.mp.dps`.  At mpmath's default `dps=15`, this becomes

```text
tie_tolerance = 10^5,
```

so every finite pressure difference in the enumerated range is declared a
tie.  For example, at child order four, parent order eight, and `beta=1`,
the all-positive signing and the true minimizing histogram have log
pressures separated by

```text
0.052469120187727813...
```

but both were labelled minimizing in the old run.

The patched driver now:

1. sets `mp.mp.dps=80` before the first child selection;
2. records this precision in the output parameters; and
3. rejects `--mp-dps<=20`, for which the selector's guard-digit convention
   would again be meaningless.

A corrected smoke run at `N=4,beta=1` returns exactly the two orientations
of one left/right minimizer pair.  Across `N=4,...,9` and `beta=1,2,4`, the
80-digit selector has one minimizing class for each child order and a
minimum gap to the next absolute-energy histogram of

```text
0.0372541530230792926864101620013847331494...
```

This is more than fifty-eight decimal orders above the `10^(-60)` tie
tolerance used by the corrected run.

## 2. Which committed claims are affected

The two directly affected reports are:

- `experiments/actual_child_projective_synchronization_report.md`;
- `experiments/actual_child_tilted_average_influence_report.md`.

The following formulations in those reports are invalid as written:

- "every contracted-temperature minimizer class" and "all 102
  actual-child laws";
- every displayed range whose extrema were taken over all 102 records;
- the claimed 102-case actual-minimizer falsification of worst-context
  projective synchronization;
- the claimed 102-case optimizer evidence for the relative sharpness of
  tilted/conditional influence and for the total-correlation versus
  marginal-drift split;
- the phrase "all 102 balanced actual minimizer-class/orientation cases" in
  `PROGRESS_FRONTIER.md` checkpoint 19.

The reports are numerical evidence rather than inputs to a proof.  No
displayed theorem inequality, normalization, or KL chain rule was inferred
from the 102-case extrema.  In particular, the following remain intact:

- the exact identities IC.7, ES.21--ES.28 and Theorems 37.24 and 37.27;
- the actual `N=8` target-reaching reverse-product certificate, whose driver
  separately sets `mp.mp.dps=80`;
- the order-eight radial/sector actual-minimizer falsifier, which imports
  only the `projective_record` calculation and supplies independently
  certified order-eight children;
- the low-transport feature ceiling and the sector--Gram/cluster promotion
  theorem, neither of which uses the finite projective JSON.

## 3. What survives after filtering the frozen data

The old JSON records child SHA values, so one can compare them against the
80-digit winners without rerunning a bridge calculation.  Exactly `36` of
the `102` records survive: six parent orders, three beta values, and two
orientations of the unique minimizing pair.  On this valid subset the main
qualitative diagnostics happen to survive, but the published ranges and
case count do not.

```text
beta   IC.15/J on valid records     IC.23/J on valid records
 1     [1.990832, 2.064579]         [1.989529, 2.062491]
 2     [1.697290, 2.137839]         [1.658494, 2.082232]
 4     [1.734799, 2.491237]         [1.486425, 2.126007]
```

At `N=9`, the valid worst-projective-bound ratios are respectively
`97.13`, `43.94`, and `13.52` for `beta=1,2,4`.  Thus worst-context diameter
is still loose on the actual finite winners, but the old statement "at
least 74.2 over all 102 actual cases" used a false population description.

Likewise, on the valid subset the integrated row-total-correlation shares
range from approximately `.9721` to `1` at `beta=1`, `.6121` to `1` at
`beta=2`, and `.3462` to `1` at `beta=4`.  Hence both dependence and
marginal drift are still visible on genuine finite winners.  This is a
post-hoc filtered diagnostic, not a replacement frozen experiment and not
an asymptotic theorem.

## 4. Required repository handling

The frozen JSON and the two reports should not be cited with the old
102-case optimizer interpretation.  This audit deliberately does not
overwrite them: regeneration changes a large committed result file and its
summary tables and should be done as a separate, explicit evidence-repair
campaign.  The present patch only prevents recurrence of the selector bug.
