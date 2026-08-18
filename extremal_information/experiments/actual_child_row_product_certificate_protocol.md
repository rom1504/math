# Actual-child row-product global-certificate protocol

Run from the repository root with the project environment:

```bash
.venv/bin/python \
  extremal_information/experiments/actual_child_row_product_certificate.py
```

The output is
[`../../computations/results/actual_child_row_product_certificate.json`](../../computations/results/actual_child_row_product_certificate.json).

## Scope frozen by the script

The three audited laws are:

| total order | split | `beta` | `lambda` | relative orientation | row-Walsh masks |
|---:|:---:|---:|---:|---:|:---:|
| 8 | `4+4` | 4 | `5.382104195764755` | -1 | `(2,1,4,8)` |
| 9 | `4+5` | 2 | 1 | -1 | `(1,1,2,2)` |
| 9 | `4+5` | 4 | 1 | -1 | `(1,1,4,4)` |

For every child order and temperature, the child signing is selected by the
existing complete signing/absolute-energy-histogram enumeration.  The program
then independently encloses every distinct histogram pressure with interval
arithmetic and records a positive interval gap proving that the selected
histogram is strictly minimal.  Its signed-permutation/global-sign
classification is exact.  The same-temperature target and the negative-moment
soft pressure are likewise recomputed as outward intervals.

## Evidentiary classes

1. **Exact finite combinatorics.**  Every bridge is enumerated.  Its pressure
   is represented by the integer histogram of
   `(abs(internal child energy), abs(rank-one bridge energy))`.  Equality of
   signatures and equality of signature-count columns are integer facts.
2. **Outward interval evaluation.**  Every distinct integer signature is
   evaluated with `mpmath.iv`; coarse escort probabilities and Walsh
   coefficients are saved as outward intervals.
3. **Computer-assisted inequality.**  The order-eight coarse reverse-product
   KL lower bound uses two-dimensional interval subdivision after the exact
   twin reduction.  Entropy endpoints are evaluated by outward interval
   arithmetic from their exact dyadic inputs.  Every remaining binary64 box
   lower bound includes an additional `1e-12` downward safety subtraction;
   the fixed arithmetic has fewer than 64 operations, all intermediates have
   magnitude below 16, and the corresponding conservative IEEE-754 error
   budget is below `2.3e-13`.  The proof target `1.075` is separated from the
   feasible value by more than `6e-4`.
4. **Rigorous analytic certificate.**  At order nine, an outward upper bound
   on the coarse rectangle-matrix row sum is below four.  The rectangle-
   Hessian theorem therefore makes the complement-symmetric uniform product
   the exact unique coarse minimizer.
5. **Numerical falsification diagnostic only.**  The complete-row rectangle
   spectral radii are computed from the complete pressure cube in
   `numpy.longdouble/float64`.  They exceed the required threshold by factors
   between `3.75` and `47.99`; they are not used for a positive proof.
6. **Feasible upper bounds.**  Full row-product upper endpoints are imported
   from the preregistered coordinate-Gibbs files and padded upward by `1e-6`.
   They are evaluated feasible-law bounds, not claims of global optimality.
7. **Threshold audit.**  The target-excess lower bound compares the outward
   lower endpoint for `V_row` directly with the outward upper endpoint for the
   same-temperature child target.  It does not assume that the binary64
   threshold root is exact.

The theorem and proof are in
[`../drafts/actual_child_row_product_global_certificate.md`](../drafts/actual_child_row_product_global_certificate.md).
