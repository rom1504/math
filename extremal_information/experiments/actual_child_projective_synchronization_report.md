# Finite actual-child test of projective cavity synchronization

Classification: **complete finite bridge enumeration; numerical evidence,
not an interval certificate and not an asymptotic claim**.

The implementation
[`actual_child_projective_synchronization.py`](actual_child_projective_synchronization.py)
tests CC.14 on every contracted-temperature minimizer class for balanced
orders `N=4,...,9`, both orientations, `beta=1,2,4`, and `lambda=1`.
For every bridge context it computes exactly up to floating transcendental
evaluation

```math
 \delta_i(B_{-i})
 =\operatorname{osc}_{b_i}
   \{\log p(b_i,B_{-i})-\log p_i(b_i)\},
```

then records the uniform `Delta^2`, its canonical-product average, the
Efron--Stein terms in CC.17b, `Var_r(h)`, and the exact-cube canonical
cumulant `J`.

## Result

The uniform sufficient criterion is already very loose on the actual
optimized children.  At the largest enumerated order:

| `N=9` | range of `Delta^2/N` | range of `bar Delta^2/N` | range of `Var_r(h)/N` | range of `J/N` |
|---|---:|---:|---:|---:|
| `beta=1` | `[.6731,1.2761]` | `[.07039,.18468]` | `[.001840,.004514]` | `[.000866,.002151]` |
| `beta=2` | `[8.9259,13.9253]` | `[2.3697,4.4881]` | `[.05442,.09156]` | `[.02539,.04583]` |
| `beta=4` | `[46.6241,76.9124]` | `[15.0962,39.4735]` | `[.29434,1.77863]` | `[.12640,.43119]` |

For `beta=1`, where the canonical interaction itself is tiny, the best
CC.2 bound-to-`J` ratio rises from `2.75` at `N=4` to at least `74.2` at
`N=9`.  Thus worst-context projective synchronization is not tracking the
observed interaction scale.  At `beta=2,4`, both the average curvature and
the canonical cumulant remain visibly extensive over this finite range.

This is a finite **falsifier of CC.14 as an efficient next target**, not a
proof that `Delta^2=o(N)` or `J=o(N)` fails asymptotically.  The more credible
observable is an averaged cavity sensitivity propagated along the hybrid
interaction path, because CC.17b controls only its product endpoint.

## Internal checks

- All row marginal log-likelihoods agree to at most
  `3.56e-15`, independently confirming CR.0 on all 102 cases.
- Product probabilities normalize with log error at most `3.56e-15`.
- Every case satisfies the rigorous finite inequality
  `J<=Delta^2/8`; the smallest numerical bound-to-`J` ratio is `2.747`.
- Every case satisfies the independently computed Efron--Stein inequality
  in CC.17b.

The machine-readable records are in
[`../../computations/results/actual_child_projective_synchronization.json`](../../computations/results/actual_child_projective_synchronization.json).
