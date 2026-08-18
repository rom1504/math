# Finite audit of the raw actual-child negative overlap

Classification: **complete finite child and bridge enumeration through
balanced `N=9`, with numerical transcendental evaluation and quadrature; no
asymptotic claim**.  The protocol and reproduction commands are in
[`actual_child_negative_overlap_protocol.md`](actual_child_negative_overlap_protocol.md).

The implementation evaluates the exact target

```math
\widehat\rho_N^-(\lambda)
={1\over\lambda mn}\int_{-\lambda}^0
 E_{\widehat\Pi_s}\sum_e r_e(B_{-e})^2\,ds
```

on contracted-temperature pressure-minimizing children.  It enumerates all
bridges, not a conference/Paley surrogate.  The pressure-flip formula agrees
with an independent direct Gibbs cavity calculation on the deterministic
check masks to maximum absolute error `2.30e-14`.

## Complete-cube values at `lambda=1`

The table gives the range over all minimizing child classes and both relative
orientations.  Singleton ranges are displayed as one value.

| `N` | split | `beta=1` | `beta=2` | `beta=4` |
|---:|---:|---:|---:|---:|
| 4 | `2+2` | `.017438` | `.193866` | `.359952` |
| 5 | `2+3` | `.025222` | `.248944` | `.783333` |
| 6 | `3+3` | `[.027896,.030501]` | `[.219979,.235506]` | `[.498687,.606257]` |
| 7 | `3+4` | `.027588` | `.218858` | `.486662` |
| 8 | `4+4` | `.027121` | `.216336` | `.492500` |
| 9 | `4+5` | `.023450` | `.202726` | `.556071` |

The separately enumerated comparable `3+7` split at `N=10` gives, over its
two orientations,

| `beta` | path-average range at `lambda=1` |
|---:|---:|
| 1 | `[.020203,.020300]` |
| 2 | `[.188342,.188759]` |
| 4 | `[.464091,.464537]` |

At the recorded `N=8,beta=4` target value
`lambda=5.382104195764755`, the complete-cube path average is `.525128`,
the midpoint value is `.525278`, and the actual negative-disorder endpoint
is `.575897`.  Thus the large value is not an artifact of including only the
fair endpoint `s=0`.

## Held-out balanced sampling

The separate reproducible sampled-bridge audit exhaustively selects the
actual children, evaluates the full finite Gibbs response on every sampled
bridge, and checks its `N=8` estimates against the complete cube.  At
`lambda=1` it gives:

| `N` | samples | `beta=1` | `beta=2` | `beta=4` |
|---:|---:|---:|---:|---:|
| 8 | `100000` | `.027064(79)` | `.216292(330)` | `.492755(404)` |
| 10 | `100000` | `.021491(58)` | `.194996(310)` | `.469165(360)` |
| 12 | `50000` | `.016301(64)` | `.176755(427)` | `.446383(488)` |
| 14 | `10000` | `.021445(184)` | `.201932(1078)` | `.461008(1286)` |

An additional held-out `N=16` run exhaustively selected both minimizing
signed-permutation classes and sampled every ordered pair of classes.  With
`2000` independent bridges per pair, the ranges over the four pairs are:

| `N` | child-class pairs | `beta=1` | `beta=2` | `beta=4` |
|---:|---:|---:|---:|---:|
| 16 | `4` | `[.018077,.018793]` | `[.187500,.191936]` | `[.451017,.456478]` |

Individual one-standard-error estimates are retained in the JSON record;
the largest standard errors are `.000379`, `.002354`, and `.003247`,
respectively.  This extension is still finite numerical evidence, but it
checks that the observed strong-channel floor is not an artifact of choosing
one optimizer class.

Parentheses contain one self-normalized importance-sampling standard error
in the last displayed digits.  All three `N=8` estimates are within one
standard error of the complete-cube values.

## Falsification judgment

At `beta=2,4`, neither the complete values through `N=10` nor the held-out
balanced values through `N=14` support a naive decay theorem; the observed
path averages remain of constant order.  At `beta=1`, the values are much
smaller but the range is too short and nonmonotone to certify a decay rate.

This is the strongest finite falsifier presently available, but it is **not**
a scalable obstruction: the sequence of child minimizers is certified only
at the displayed finite orders, and no theorem proves a positive limiting
overlap.  Consequently these data alone are a STRIKE rather than a RESET
under the campaign rule.

Complete records:

- [`../../computations/results/actual_child_negative_overlap_exact.json`](../../computations/results/actual_child_negative_overlap_exact.json)
- [`../../computations/results/actual_child_negative_overlap_exact_n10_3x7.json`](../../computations/results/actual_child_negative_overlap_exact_n10_3x7.json)
- [`../../computations/results/actual_child_negative_overlap_sample.json`](../../computations/results/actual_child_negative_overlap_sample.json)
- [`../../computations/results/actual_child_negative_overlap_sample_n16.json`](../../computations/results/actual_child_negative_overlap_sample_n16.json)
