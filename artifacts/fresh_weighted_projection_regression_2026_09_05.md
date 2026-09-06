# A nonconference regression for the weighted projection theorem

Date: 2026-09-05. The finite variance identity is exact algebra; the sampled
energy measurements below are statistical diagnostics, not proof certificates.

The reproducible program
`computations/fresh_weighted_projection_tensor_regression.py` uses a fixed
symmetric full signing C of order five and the regular Hadamard matrices
`H_s=(J_4-2I_4)^(tensor k)`. It hollows `C tensor H_s` before constructing
`B=A/sqrt(N-1)`. Matrix products use the exact Kronecker representation;
the smallest instance is checked against a dense implementation.

For the full normalized sequence, the five local cubic variances are

```
diag[B (B^2)^(circ3) B] = (733,733,733,733,533)/625.
```

Thus a pointwise lower bound of one is false, on a scalable bounded-operator
family. The proved mean-standard-deviation inequality is the needed statement.
Hollowing has a vanishing effect on these fixed-degree variances.

The regression tests the actual cubic field `Z=B h_3(BS)`, the normalized
weights `d_i=1/sqrt(max(v_i,1/4))`, and the even gate
`M_i=1{|X_(three-vertex tree),i|<1}`. It checks both a predicted zero
projection against the old tree and a nonzero projection against the cubic
channel. The exact injective old field includes its own-root corrections;
it is not replaced by a naive recursion.

With seed 2026090523 and 4096 independent spin inputs at each order
20, 80, 320, 1280, 5120, 20480, the largest instance gives:

| Quantity | Measured value | Standard error or prediction |
|---|---:|---:|
| Old-tree weighted projection | -0.000109153 | SE 0.000127775 |
| Cubic weighted projection | 0.717745918 | SE 0.000740150 |
| Cubic limiting-form prediction at this order | 0.717533727 | analytic finite variance input |

These are consistency checks on a family with nonconstant local variances,
not confidence bounds uniform over signings or an asymptotic extrapolation.
The proof of the universal weighted identity is separate.

The program emits per-order JSON to standard output and saves the complete
run to `computations/results/fresh_weighted_projection_tensor_regression.json`.
Reproduce with `.venv/bin/python -B computations/fresh_weighted_projection_tensor_regression.py`;
the seed, sample count, batch size, and order list are fixed in the program.
