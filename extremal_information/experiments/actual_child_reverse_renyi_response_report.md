# Exact finite audit of the reverse-Renyi response identity

Status: **complete finite enumeration with numerical transcendental
evaluation**.  This is a theorem normalization check and a finite structural
diagnostic, not asymptotic evidence by itself.

The executable is
[`actual_child_reverse_renyi_response.py`](actual_child_reverse_renyi_response.py)
and the frozen output is
[`actual_child_reverse_renyi_response.json`](../../computations/results/actual_child_reverse_renyi_response.json).
It uses complete-enumeration thermal child minimizers, then enumerates all
`2^16` bridges for the balanced order-eight split at `beta=4`.

## Verification

For both relative orientations, `lambda=1` and the independently located
target threshold `lambda=5.382104`, and four bridge amplitudes, the exact
analytic derivative in (RR.7) agrees with a centered finite difference of
the reverse Renyi divergence to at worst `6.1e-10`.  The conditional-bias
identity (RR.20) agrees directly to at worst `4.1e-15`.

The two orientations are numerically identical for this equal-child case.
At the physical amplitude `u=4/sqrt(8)`:

| `lambda` | `R_lambda` | `D(U||Pi)` | centered gain | `S/(mn)` | leave-one-out MI sum |
|---:|---:|---:|---:|---:|---:|
| `1` | `5.86045` | `4.39026` | `1.47019` | `.46336` | `3.30071` |
| `5.382104` | `7.93297` | `4.39026` | `3.54271` | `.57590` | `11.02956` |

Thus the target-reaching finite escort has a dense planted-response phase:
more than half of the maximum possible squared extrinsic-coordinate response
is present on average at the physical amplitude, and the summed
coordinate-to-rest mutual information already exceeds the total order.
This is fully consistent with Corollary RR.2.

The same enumeration also verifies the canonical row-erased decomposition
(CR.5)--(CR.6) to absolute residual below `9e-16`.  At the physical
amplitude its values are:

| `lambda` | sum of erased-row works | `E_r h` | canonical inverse work | centered interaction cumulant `J` | `J/lambda` |
|---:|---:|---:|---:|---:|---:|
| `1` | `5.11633` | `2.17356` | `2.94277` | `2.91768` | `2.91768` |
| `5.382104` | `8.21719` | `4.72380` | `3.49339` | `23.89429` | `4.43958` |

Thus the explicit iid-row certificate obtained by erasing all other rows is
not close to the target-reaching order-eight escort.  This is a finite
falsifier of *canonical row-erased synchronization*, not an asymptotic
lower bound: a different row-product law can do better, and the separately
certified optimum has reverse gap only known to lie between `1.075` and
`4.506450` at this order.

## Scope

The physical raw amplitude at order eight is `1.414`, whereas for fixed
`beta` it tends to zero asymptotically.  The numbers therefore do not prove a
positive limiting response density.  They do decisively test the algebra and
show that the response statistic is nonvacuous on the actual target-reaching
child instance.  The canonical decomposition isolates a large finite
cross-row interaction cumulant, but it does not prove that the *optimal*
directed row-product gap is extensive.  Its asymptotic scale on actual
optimizing children is the remaining theorem target.
