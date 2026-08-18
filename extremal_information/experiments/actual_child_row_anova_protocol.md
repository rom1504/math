# Actual-child row-ANOVA verifier protocol

This is the finite normalization check for
[`../drafts/actual_child_row_anova_infinitesimal.md`](../drafts/actual_child_row_anova_infinitesimal.md).

Run:

```bash
.venv/bin/python \
  extremal_information/experiments/verify_actual_child_row_anova.py
```

The default case is the completely enumerated `N=6`, split `3+3`,
`beta=1`, orientation `epsilon=+1` instance.  Each child is selected by
complete signing enumeration and high-precision comparison of its contracted-
temperature pressure.  The bridge cube is enumerated completely.  The
row-product objective at the displayed small `lambda` values is solved by
coordinate Gibbs updates from the uniform law; this is a numerical
critical-point verification, not an independent global certificate (global
uniqueness for sufficiently small `lambda` follows from Theorem RA.1).

The generated result is
[`../../computations/results/actual_child_row_anova_verify.json`](../../computations/results/actual_child_row_anova_verify.json).
It records:

- physical row ANOVA variances and the exact mixed-row response check;
- the child overlap coefficient in RA.17;
- convergence of `sigma_cross^2(L_u)/u^4` and
  `sigma_add^2(L_u)/u^4` to RA.17--RA.18;
- convergence of the four small-`lambda` ratios in RA.5--RA.9 to one.

For the default case, the overlap limits are

```text
K_cross = 0.7951276391239326
K_add   = 0.2650425463746443
```

At bridge amplitudes `0.1,0.05,0.025,0.0125`, the measured cross ratios are

```text
0.7666277212, 0.7877780018, 0.7932758992, 0.7946638038.
```

At disorder temperatures `0.1,0.05,0.025,0.0125`, the ratios of the exact
escort gain, evaluated row-product gain, reverse projection, and total
correlation to their respective leading terms all approach one.  At the
smallest value they are respectively

```text
0.9984256109, 0.9993091417, 0.9982035246, 0.9972846690.
```

These are finite numerical checks of coefficients and conventions only.
They are not evidence for uniformity in `N`.

