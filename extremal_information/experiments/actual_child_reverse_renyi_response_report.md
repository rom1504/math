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

The collision--cavity projective ranges in CC.14 were also enumerated.  At
the physical amplitude every row has supremum range `10.10199`, so
`Delta^2=408.20046`.  The resulting CC.15 Hoeffding upper bounds are
`51.0251` at `lambda=1` and `1478.0451` at `lambda=5.382104`, compared with
the exact cumulants `2.91768` and `23.89429`.  Even replacing each supremum
by its squared product-law mean leaves sums `205.7731` and `217.7466`.
Thus worst-case projective synchronization is a rigorous sufficient
criterion but is very loose on this finite target-reaching law.  Any useful
asymptotic proof will probably need a tilted/typical-range inequality rather
than a raw supremum, unless optimizing-child rigidity improves the ranges
dramatically with order.

The hybrid interaction path IC.3 is much sharper on the same data.  A
17-point trapezoidal evaluation gives

| `lambda` | exact `J` | curvature integral | optimized conditional-entropic integral | unoptimized tilted-influence integral | projective-sup bound |
|---:|---:|---:|---:|---:|---:|
| `1` | `2.91768` | `2.92321` | `4.60221` | `5.27313` | `51.0251` |
| `5.382104` | `23.89429` | `24.60695` | `48.46409` | `96.60183` | `1478.0451` |

The curvature column checks the exact identity IC.7 (the displayed residual
is numerical quadrature error).  Both influence columns are also quadrature,
not rigorous finite upper certificates.  Optimizing the row-deleted
comparison as in IC.4 retains another substantial amount of cancellation:
at the target threshold it halves the unoptimized bound.  The exact slack
identity ES.22--ES.23 shows that the remaining gap above `J` is integrated
dual total correlation.  Thus conditional-entropic influence is a sharp
localization diagnostic, but not by itself a strict reduction of the
canonical error.

## Scope

The physical raw amplitude at order eight is `1.414`, whereas for fixed
`beta` it tends to zero asymptotically.  The numbers therefore do not prove a
positive limiting response density.  They do decisively test the algebra and
show that the response statistic is nonvacuous on the actual target-reaching
child instance.  The canonical decomposition isolates a large finite
cross-row interaction cumulant, but it does not prove that the *optimal*
directed row-product gap is extensive.  Its asymptotic scale on actual
optimizing children is the remaining theorem target.
