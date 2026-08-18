# Exact finite audit of the actual parent response norm

Status: **complete finite computation, asymptotic status open**.

This experiment tests the intrinsic response statistic in (SH.0c)--(SH.0h):

```math
R(B)=\left\|\mathbb E_{\nu_{\epsilon,B}}
             [\tau XY^{\mathsf T}]\right\|_F.
```

Run

```bash
.venv/bin/python \
  extremal_information/experiments/actual_child_response_norm_exact.py
```

to regenerate
[`../../computations/results/actual_child_response_norm_exact.json`](../../computations/results/actual_child_response_norm_exact.json).
At every total order `4<=N<=8`, the program:

1. enumerates all child signings and selects every contracted-temperature
   pressure-minimizing signed-permutation class;
2. takes the balanced split and both relative orientations;
3. enumerates the complete bridge cube;
4. evaluates the exact finite Gibbs sum for every entry of
   `E[tau XY^T]` in floating point.

No conference/Paley surrogate and no bridge sampling is used.  A direct
central-difference check of one response coordinate at `N=6,beta=1` agreed
to `2.5e-11`.

## Main finite observation

The statistic sharply separates the displayed inverse temperatures already
at these small orders.  The table gives the mean of `R(B)/N` for one
orientation (the other agrees except where inequivalent child orientations
exist):

| `N` | `beta=0.5` | `beta=1` | `beta=2` | `beta=4` |
|---:|---:|---:|---:|---:|
| 4 | 0.1225 | 0.2299 | 0.3528 | 0.4159 |
| 5 | 0.1081 | 0.2078 | 0.3253 | 0.3599 |
| 6 | 0.1012 | 0.1996 | 0.3338 | 0.4239 |
| 7 | 0.0930 | 0.1859 | 0.3183 | 0.4107 |
| 8 | 0.0882 | 0.1788 | 0.3145 | 0.4087 |

For `beta=0.5`, multiplying the mean ratio by `sqrt(N)` gives approximately
`0.25` throughout, consistent with the tight `R(B)=Theta(sqrt(N))` response
certified by strict-high-temperature covariance.  At `beta=4`, the ratio is
already roughly constant and close to the maximal balanced value `1/2`,
indicating macroscopic response rather than tight response.

At `N=8`:

- for `beta=0.5`, the median, 99th percentile, and maximum of `R/N` are
  `0.0871, 0.1072, 0.1457`;
- for `beta=2`, they are `0.2975, 0.4861, 0.4990`, and `98.5%` of bridges
  have `R/N>0.2`;
- for `beta=4`, they are `0.4255, 0.4998, 0.5000`, and `60.4%` have
  `R/N>0.4`.

## Evidentiary judgment

This is a finite falsifier of the blanket hypothesis that actual optimized
children always induce submacroscopic parent response at all temperatures.
It does **not** prove an asymptotic lower bound: orders through eight cannot
distinguish a positive limiting ratio from slow decay.  Nor does a uniform-
bridge quantile itself construct the convex typical carrier required in
Lemma SH.0.

The experiment supports a narrower campaign split:

- at sufficiently small `beta`, prove a convex-carrier
  `R(B)=O(sqrt(N))` theorem (SH.1 gives one route under a child spectral
  margin);
- at larger `beta`, test whether the apparent macroscopic response is
  asymptotically unavoidable and whether it lies in the row-product gain or
  directed-dependence side of (AC.24).

