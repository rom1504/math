# Actual-child PT.1 transport-gap candidate diagnostic

## Question and scope

For the centrally symmetric rank-one child prior `mu`, put `d=mn` and define
the normalized PT.1 gap of a nonzero candidate
`A in conv(supp(mu))` by

```math
c_\mu(A)
={1\over2}
-{\sqrt d\over\|A\|_F^2}
 \log E_\mu\exp\left\{{\langle A,Q\rangle\over2\sqrt d}\right\}. \tag{TG.1}
```

The transport hypothesis PT.1 with constant `c>0` is exactly the assertion
that `c_mu(A)>=c` for **every** nonzero `A` in the child-prior convex hull.
This experiment asks only whether two principled small candidate families
already make (TG.1) negative.

It uses the certified actual contracted-temperature minimizing children at
`beta=4`, splits `3+3`, `4+4`, and `3+7`, and both relative orientations.
Child signing and energy-histogram enumeration is exhaustive, every competing
histogram is interval separated at 80 digits, and every bridge cube is
complete.  The LP, covariance, posterior, and MGF calculations are numerical
finite evaluations.  Most importantly, the search over `A` is **not
exhaustive**.

## Candidate families

### Covariance-top radial family

Let

```math
\Gamma_\mu=E_\mu[\operatorname{vec}(Q)\operatorname{vec}(Q)^T]
```

and let `v` be its unit top eigenvector.  Central symmetry makes the mean
zero.  The infinitesimal normalized gap along this ray is

```math
c_{\rm tan}={1\over2}-{\lambda_{\max}(\Gamma_\mu)\over8\sqrt d}. \tag{TG.2}
```

The program computes the atomic norm of `v` by linear programming against
all projective rank-one prior atoms.  This gives the exact radial boundary
`a_max v` of `conv(+-supp(mu))` up to LP precision.  It evaluates (TG.1) at
fractions `.01,.05,.10,.25,.50,.90` of that radius.  Thus the finite ray
candidates have explicit convex-hull witnesses; the tangent value is only
their small-scale limit.

### Collision-selected posterior barycenters

For every complete bridge, the full posterior barycenter

```math
M(B)=E_\mu[Q\mid B]
```

belongs automatically to the required convex hull.  After computing the
exact doubled-channel collision factor `K_0(B)`, the program retains:

- the 24 antipodal bridge pairs with largest `K_0`;
- the 24 antipodal bridge pairs with largest individual `q(B)K_0(B)`
  contribution, where `q` is used only for candidate ranking and has
  `lambda=1`.

The posterior is reconstructed from pressure flips and independently from a
direct latent Gibbs sum.  The maximum coordinate discrepancy in the output
is below `6e-14`.

## Results

No tested candidate has negative normalized gap.

| `N`, split, orientation | covariance tangent | minimum finite covariance ray | minimum posterior barycenter | `||A||_F^2/d` at posterior minimum |
|---|---:|---:|---:|---:|
| `6`, `3+3`, `-` | .416989 | .416989 | .455948 | .403978 |
| `6`, `3+3`, `+` | .132611 | .132624 | .440314 | .250296 |
| `8`, `4+4`, `-` | .290295 | .290298 | .467784 | .940957 |
| `8`, `4+4`, `+` | .290295 | .290298 | .467784 | .940957 |
| `10`, `3+7`, `-` | .206034 | .206041 | .496339 | .708071 |
| `10`, `3+7`, `+` | .393349 | .393349 | .482178 | .983524 |

The smallest observed finite-candidate gap is `.132624`, on the
`N=6,+` covariance ray.  The corresponding tangent limit is `.132611`.
The posterior candidates remain farther from violation:

| `N`, orientation | minimum among top-`K_0` pairs | minimum among top-`qK_0` pairs |
|---|---:|---:|
| `6`, `-` | .455948 | .455948 |
| `6`, `+` | .451978 | .440314 |
| `8`, either | .467784 | .473642 |
| `10`, `-` | .496339 | .496646 |
| `10`, `+` | .482178 | .491180 |

This is not an artifact of testing only tiny barycenters.  At `N=10`, the
posterior minimizers in the displayed family have normalized squared norms
`.7081` and `.9835` in the two orientations.

The maximum LP reconstruction residual for covariance-ray membership is
below `1e-15`.  The complete machine-readable output records every candidate,
its MGF, norm, bridge mask and collision severity, plus the child selection
certificates.

## Finite research judgment

These natural child-only candidates do **not** falsify PT.1.  At the tested
orders they are consistent with a fixed positive transport gap, with an
observed candidate-family margin above `.13`.  This is useful directional
evidence because the tested posterior barycenters include the strongest
collision bridges and are often macroscopically normed.

It is not a certificate for PT.1.  The convex hull contains many directions
outside the covariance-top ray and the 24+24 selected posterior pairs.  A
negative direction could occur elsewhere, and six finite records cannot
establish an all-order constant.  The experiment therefore supports the
transport branch without changing its smallest missing lemma: one still
needs either an analytic uniform inequality over the entire actual-child
convex hull or a genuine violating direction.

## Reproduction

```bash
.venv/bin/python \
  extremal_information/experiments/actual_child_prior_transport_gap_diag.py
```

Machine-readable output:
[`../../computations/results/actual_child_prior_transport_gap_diag.json`](../../computations/results/actual_child_prior_transport_gap_diag.json).
