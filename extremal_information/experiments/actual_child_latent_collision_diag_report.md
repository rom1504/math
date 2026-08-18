# Actual-child latent posterior collision diagnostic

## Question and finite scope

This computation evaluates the posterior collision factors in the global
latent-coreset theorem on the **actual contracted-temperature minimizing
children**, rather than on a conference or generic surrogate.  It uses
`beta=4`, `lambda=1`, the splits `3+3`, `4+4`, and `3+7`, both relative child
orientations, and every bridge in each complete cube.

The children are selected by exhaustive signing and exact energy-histogram
enumeration.  Every competing histogram is separated by interval arithmetic
at 80 decimal digits.  The hyperbolic evaluations and the collision summaries
are reproducible floating computations.  All conclusions below are therefore
**finite evidence only**, not asymptotic lower bounds.

## Exact doubled-channel computation

Put `t=beta/sqrt(N)`.  For fixed children and orientation `epsilon`, let

```math
z_u(B)=E_{x,y}
 \left[
 \cosh\{t(E_L(x)+\epsilon E_R(y))\}
 \cosh\{u x^TBy\}
 \right],
\qquad z_0=E_{x,y}\cosh\{t(E_L+\epsilon E_R)\}.
```

The child prior remains at temperature `t` when the bridge channel is
doubled.  The exact posterior-collision identity becomes

```math
K_0(B)=1+\chi^2(\mu_B\Vert\mu)
      ={z_{2t}(B)z_0\over z_t(B)^2}.             \tag{LC.1}
```

For a deleted edge,

```math
z_{u,e}(B_{-e})
 ={z_u(B_e=+,B_{-e})+z_u(B_e=-,B_{-e})\over2\cosh u},
\qquad
K_e={z_{2t,e}z_0\over z_{t,e}^2}.               \tag{LC.2}
```

The actual inverse escort is `q(B) proportional to z_t(B)^(-lambda)`.  The
program evaluates

```math
\overline K_{\rm full}=E_qK_0,
\qquad
\overline K_{\rm del}={1\over mn}\sum_eE_qK_e. \tag{LC.3}
```

by XOR--Walsh convolution.  Three direct Gibbs sums per record independently
check both the `t` and `2t` arrays; the largest log-pressure error is below
`3e-8`.

## Collision growth on the complete cubes

| `N`, split, orientation | `log Kbar_full/N` | `log Kbar_del/N` | `E_q log K_0/N` | `E_(q,e) log K_e/N` |
|---|---:|---:|---:|---:|
| `6`, `3+3`, `-` | 1.006858 | 1.010416 | .850074 | .730024 |
| `6`, `3+3`, `+` | .995969 | 1.026246 | .896136 | .879468 |
| `8`, `4+4`, `-` | 1.045065 | 1.028473 | .815143 | .773770 |
| `8`, `4+4`, `+` | 1.045065 | 1.028473 | .815143 | .773770 |
| `10`, `3+7`, `-` | .905747 | .895942 | .690290 | .670028 |
| `10`, `3+7`, `+` | .816051 | .808743 | .630309 | .605768 |

The two `N=8` orientations are numerically identical, as expected from their
signed equivalence.  Orientation is not silently optimized or averaged in
any row of the table.

At these orders the collision is not produced only by an isolated maximum.
The following are the median, 90th percentile, and 99th percentile of
`log K/N` under the declared base law (`q` for `K_0`, and uniform edge times
the deleted-edge marginal of `q` for `K_e`):

| `N`, orientation | full `q50/q90/q99` | deleted `q50/q90/q99` |
|---|---:|---:|
| `6`, `-` | `.954/1.132/1.132` | `.611/1.110/1.362` |
| `6`, `+` | `.904/1.170/1.170` | `.952/1.137/1.205` |
| `8`, either | `.779/1.137/1.333` | `.800/1.171/1.318` |
| `10`, `-` | `.646/.962/1.153` | `.645/.942/1.154` |
| `10`, `+` | `.667/.828/1.113` | `.629/.832/1.058` |

Here the first full entry in the `N=6,-` row is `.954`; the leading period in
the compact triple is decimal notation, not an omitted digit.

Fixed-threshold tails make the finite distinction from a sublinear-size
posterior especially visible.  Each entry below is `full/deleted` base mass:

| `N`, orientation | mass with `log K/N>.50` | mass with `log K/N>.75` | mass with `log K/N>1` |
|---|---:|---:|---:|
| `6`, `-` | `.774/.654` | `.731/.430` | `.299/.430` |
| `6`, `+` | `.888/.783` | `.888/.783` | `.209/.412` |
| `8`, either | `.958/.885` | `.561/.548` | `.248/.161` |
| `10`, `-` | `.934/.912` | `.276/.243` | `.077/.071` |
| `10`, `+` | `.711/.625` | `.226/.258` | `.029/.019` |

In particular, more than 62% of both `N=10` base laws have
`log K>.5N`; the large mean is not explained solely by a negligible tail.

## How rare is the mass carrying the collision mean?

For a collision variable `K` under base law `p`, define its size-biased law

```math
d\widehat p={K\over E_pK},dp.
```

The table reports (i) the smallest base mass encountered by sorting `K`
downward until 50% of `E_pK` has been accumulated, and (ii)
`exp{-D(widehat p||p)}`, an intrinsic effective base mass.

| `N`, orientation | full mass for 50% | deleted mass for 50% | full `exp(-D)` | deleted `exp(-D)` |
|---|---:|---:|---:|---:|
| `6`, `-` | .21761 | .14284 | .62506 | .39028 |
| `6`, `+` | .17720 | .24384 | .66697 | .67325 |
| `8`, either | .03241 | .03406 | .19363 | .16994 |
| `10`, `-` | .01952 | .01553 | .10128 | .08118 |
| `10`, `+` | .01534 | .00968 | .11695 | .09935 |

Thus the mean becomes increasingly tail-amplified in this short sequence:
at `N=10`, roughly 1--2% of the relevant base mass supplies half of the
collision mean.  Nevertheless the median normalized log collision remains
between `.629` and `.667` in the two `N=10` records.  The finite obstruction
is therefore not merely a vanishing collection of singular bridges at the
tested orders.

## Research judgment

The complete-cube data are compatible with **extensive latent posterior
collision**, not with an already visible `log Kbar=o(N)` regime:
`log Kbar/N` remains between `.8087` and `1.0451` from `N=6` through `10`.
The decline between the balanced `N=8` and the comparable but unequal
`N=10` split prevents any asymptotic lower-bound inference.  In particular,
these six points do not prove `Kbar>=exp(cN)`.

What they do falsify is the finite expectation that actual child optimality
has already made the prior-sampled global coreset cheap.  Both the typical
posterior collision and its tail-amplified mean are exponential-sized on all
tested actual laws.  Any proof of `log Kbar_del=o(N)` must therefore reveal a
new asymptotic mechanism not apparent in these complete cubes; alternatively,
an actual-minimizer lower bound would turn this pattern into a rigorous
obstruction to prior-sampled frame synchronization.

## Reproduction

```bash
.venv/bin/python \
  extremal_information/experiments/actual_child_latent_collision_diag.py
```

Machine-readable output:
[`../../computations/results/actual_child_latent_collision_diag.json`](../../computations/results/actual_child_latent_collision_diag.json).
