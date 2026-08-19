# Outward director review after the cross-order strikes

**Status.**  Rigorous identity and campaign-selection audit.  This note does
not improve the permanent actual-child defect exponent.  It records the
outward review required after two substantive strikes and selects one
different, directly falsifiable architecture for a later campaign.

## 1. Verdict on the completed branch

The best unconditional comparable-order estimate for the actual
finite-temperature optimizing children remains

```math
P_{m+n}(\beta)-P_m(\beta)-P_n(\beta)=O_\beta(m+n).
```

The Gaussian replacement error is `O_beta(sqrt(m+n))`, the fixed-outer
completion error in the replica construction is sublinear, and the
inhomogeneous discrepancy target has an exact power-saving implication.
None controls its remaining endpoint/core/sign-realization term.  Those
terms can still be linear.  Hence none changes the total exponent, and the
permanent SML is unchanged.

The director ranking is:

1. canonical-disorder entropy--restriction compensation;
2. scale-growing algebraic replicas;
3. joint-sector Gaussian interpolation;
4. an unspecified rare-event smoothing transform;
5. pointwise proportional restriction or near-order insertion.

Only the first item is selected below.  It leaves the struck actual-child
implementation rather than pretending to solve it.

## 2. A canonical-disorder cross-order identity

Let `A_n` be the set of hollow symmetric signings of order `n`, put

```math
L_n={n\choose2},\qquad
Q(A)=\max_x|H_A(x)|,
```

and, for fixed `beta>0`, define

```math
Z_n^{\rm dis}(\beta)
=\sum_{A\in\mathcal A_n}e^{-\beta\sqrt n\,Q(A)},
\qquad
\psi_n(\beta)={1\over L_n}\log Z_n^{\rm dis}(\beta).       \tag{DR.1}
```

Write `mu_(n,beta)` for the exact canonical law proportional to the
summand.  Fix `2<=m<N`, let `S` be a uniformly random `m`-subset of the
vertices, let `mu_S` be the marginal law of `A[S]` when
`A~mu_(N,beta)`, and put

```math
q={L_m\over L_N}.
```

The canonical law is permutation invariant, but averaging over `S` keeps
the formulas label-free.  Define the Shearer slack and restriction-energy
excess

```math
\begin{aligned}
\mathcal S_{N,m}
&=q^{-1}\mathbb E_S H(\mu_S)-H(\mu_{N,\beta}),\\
\mathcal D_{N,m}
&={\sqrt m\over q}\,
  \mathbb E_{A,S}Q(A[S])
  -\sqrt N\,\mathbb E_AQ(A).
\end{aligned}                                             \tag{DR.2}
```

Here entropy is in nats.  Shearer's inequality, applied to the family of
induced edge sets, gives

```math
\mathcal S_{N,m}\ge0.                                    \tag{DR.3}
```

Equivalently, if `U_k` is uniform disorder at order `k`, cancellation of
the ambient edge entropies gives the exact information-loss form

```math
\boxed{
\mathcal S_{N,m}
=D(\mu_{N,\beta}\Vert U_N)
-q^{-1}\mathbb E_SD(\mu_S\Vert U_m).}                   \tag{DR.3a}
```

Thus ECR asks the loss of canonical-disorder information under proportional
restriction to pay the rescaled cap-energy excess.  This is more specific
than merely asserting `S>=0`, and it identifies the only independent
inequality the proposed campaign is allowed to seek.

### Proposition DR.1 (exact entropy-compensated restriction identity)

For every finite `N,m,beta`,

```math
\boxed{
\log Z_N^{\rm dis}-q^{-1}\log Z_m^{\rm dis}
=\beta\mathcal D_{N,m}-\mathcal S_{N,m}
-q^{-1}\mathbb E_S
D(\mu_S\Vert\mu_{m,\beta}).}                            \tag{DR.4}
```

**Proof.**  The canonical entropy identity is

```math
H(\mu_{N,\beta})
=\log Z_N^{\rm dis}+\beta\sqrt N\,\mathbb E_AQ(A).       \tag{DR.5}
```

For each restriction marginal,

```math
D(\mu_S\Vert\mu_{m,\beta})
=-H(\mu_S)+\beta\sqrt m\,\mathbb E_{\mu_S}Q
+\log Z_m^{\rm dis}.                                     \tag{DR.6}
```

Average (DR.6), divide by `q`, and substitute (DR.2) and
(DR.5).  All energy and entropy terms cancel to give (DR.4). `square`

Equation (DR.4) is useful discipline: any argument that proves the proposed
bound below merely by substituting (DR.5) back into it is circular.

There is an equivalent variational formulation which makes the strength of
the missing statement transparent.  For a law `nu` on order-`k` signings,
put

```math
\mathcal G_{k,\beta}(\nu)
=\beta\sqrt k\,\mathbb E_\nu Q+D(\nu\Vert U_k).          \tag{DR.6a}
```

The canonical law uniquely minimizes this functional.  Using (DR.3a), ECR
is exactly

```math
\boxed{
q^{-1}\mathbb E_S\mathcal G_{m,\beta}(\mu_S)
\le \mathcal G_{N,\beta}(\mu_{N,\beta})
   +C_\beta N^{2-\delta}.}                              \tag{DR.6b}
```

Thus it asserts that proportional restriction of the *large-system Gibbs
optimizer* is an approximate scaled competitor at the same temperature.
Ordinary Shearer plus `Q(A[S])<=Q(A)` proves (DR.6b) only at

```math
\beta' =\beta q\sqrt{N/m}
        =\beta(m/N)^{3/2}+o(1),                         \tag{DR.6c}
```

which is the archived temperature drift.  The entire new content of ECR is
an `o(N^2)` payment for raising `beta'` back to `beta`; generic entropy
contraction leaves that payment of order `N^2`.

## 3. The selected sufficient lemma and its exact consequence

The proposed **entropy-compensated restriction** statement is: for every
fixed `beta>0`, there are `delta>0` and `C_beta<infinity`, uniform for

```math
{m\over N}\in[1/3,2/3],
```

such that

```math
\boxed{
\beta\mathcal D_{N,m}
\le\mathcal S_{N,m}+C_\beta N^{2-\delta}.}             \tag{ECR}
```

This statement has the mandatory immediate quantitative arrow.  From
(DR.4) and nonnegativity of relative entropy,

```math
\boxed{
\psi_N(\beta)
\le\psi_m(\beta)+O_\beta(N^{-\delta})}                 \tag{DR.7}
```

uniformly on the same ratios.  Indeed `q^{-1}L_m=L_N`, so division of
(DR.4) by `L_N` converts its second partition function exactly to
`psi_m`.

The uniformity in the ratio interval is essential.  Given a selected order
`r`, repeatedly replace a much larger `N` by `floor(N/2)` until the current
order lies in `[3r/2,3r]`, and then restrict directly to `r`.  Every ratio
lies in `[1/3,2/3]`, while the sum of normalized errors is `O(r^{-delta})`.
Taking `r` along a liminf subsequence proves

```math
\limsup_N\psi_N(\beta)\le\liminf_N\psi_N(\beta),         \tag{DR.8}
```

so `psi_n(beta)` converges.

Finally, if `c_n=M_n/n^(3/2)`, the exact disorder soft-minimum squeeze is

```math
-{n-1\over2\beta n}\psi_n(\beta)
\le c_n\le
{n-1\over2\beta n}\{\log2-\psi_n(\beta)\}.             \tag{DR.9}
```

Thus (ECR) for every member of an unbounded set of fixed temperatures gives

```math
\limsup_nc_n-\liminf_nc_n\le{\log2\over2\beta},          \tag{DR.10}
```

and then convergence after `beta->infinity`.

## 4. Why this is a real but high-risk architecture change

The desired recurrence controls

```math
\beta\mathcal D_{N,m}-\mathcal S_{N,m}
-q^{-1}\mathbb E_SD(\mu_S\Vert\mu_{m,\beta}),            \tag{DR.11}
```

whereas (ECR) controls the first two terms before the favorable relative
entropy is subtracted.  Hence (ECR) is a **strictly stronger sufficient
condition than the recurrence**, not a reformulation of it.  A large
marginal relative entropy can make (DR.11) small while (ECR) fails.

On the other hand, (ECR) is strictly weaker than the archived pointwise
`3/2` restriction inequality.  The latter would force

```math
{\sqrt m\over q}Q(A[S])\le\sqrt NQ(A)
```

for every relevant `A,S`, hence `\mathcal D_{N,m}\le0` and (ECR)
immediately by
(DR.3).  ECR asks only for a canonical average and permits bad restrictions
when their excess is paid by lost disorder entropy.

There is also an unavoidable scope warning.  At finite `beta`,
`mu_(N,beta)` has full support: it is the exact variational optimizer of a
new canonical-disorder free energy, not a law supported on actual cap
minimizers.  Consequently (DR.4) is not a RESET of the permanent
actual-child SML.  It is admissible only as the explicitly outward campaign
chosen after that implementation received two strikes.

Archive comparison imposes a second warning.  The partition function
(DR.1), its soft-minimum squeeze, the changing-temperature Shearer
inequality, and the support-sensitive speed-`N^2` obstruction were already
isolated in `good_signing_entropy_threshold.md`,
`traffic_laplace_principle.md`, and
`microcanonical_disorder_counting_composition.md`.  Proposition DR.1 and
the specific ECR falsifier sharpen the balance, but do not turn the old
canonical pressure into a strict reduction.  The selected campaign is
therefore a **single high-risk falsification test**, not a promoted theory:
if no independent cap-layer/incidence proof of ECR appears, this scalar
pressure architecture must be frozen immediately.

A focused primary-literature check found no theorem supplying that missing
step.  Madiman--Tetali fractional-cover entropy gives exactly the archived
temperature-rescaled Shearer inequality, not ECR.  Recent Gaussian Ising
work proves the upper ground-state tail at speed `N` and explicitly leaves
the expected lower tail at speed `N^2`; the available speed-`N^2` lower-tail
bound is spherical, Gaussian, and one-sided.  Fixed-replica negative-moment
limits for Gaussian SK are not uniform in the replica order `Theta(N)`
needed here.  Dense-graph LDPs see the zero graphon but not its
`N^{-1/2}` extremal tangent.  Relevant sources are
[Madiman--Tetali](https://arxiv.org/abs/0901.0044),
[Chen et al. (2026)](https://arxiv.org/abs/2603.06368),
[Huang--Sellke](https://arxiv.org/abs/2311.15495), and
[Chen](https://arxiv.org/abs/2311.08351).  Thus ECR would be a genuinely
new entropy--energy theorem, not an application of a known LDP.

## 5. Exact small-order falsification audit

Exhaustive enumeration through order six verifies (DR.4) to absolute error
at most `5.5e-12`.  At `beta=1`, representative tuples

```math
(\mathcal S_{N,m},\mathcal D_{N,m},
 \beta\mathcal D_{N,m}-\mathcal S_{N,m},
 q^{-1}D(\mu_S\Vert\mu_{m,\beta}))
```

are

| restriction | tuple |
|---|---|
| `4 -> 2` | `(0.257323, 0.461009, 0.203686, 0)` |
| `5 -> 3` | `(1.519574, 8.250203, 6.730628, 0)` |
| `6 -> 3` | `(3.822531, 13.217732, 9.395201, 0)` |
| `6 -> 4` | `(3.248521, 7.376322, 4.127801, 0.009374)` |
| `6 -> 5` | `(1.788389, 1.120771, -0.667618, 0.033124)` |

Exact switching-quotient enumeration extends the audit to order eight.  At
`beta=1` it gives

| restriction | tuple |
|---|---|
| `7 -> 4` | `(2.063025, 6.253040, 4.190014, 1.304102)` |
| `8 -> 4` | `(4.577970, 10.668867, 6.090897, 1.691843)` |
| `8 -> 5` | `(4.186990, 5.135347, 0.948357, 5.545085)` |

The last row is especially informative: the residual tested by ECR is
positive while the exact partition defect is `-4.596728`, because the
scaled marginal KL contributes `5.545085`.  Discarding KL is therefore a
substantial strengthening even at order eight.

The sharper warning is frozen.  For **every** `2<=m<N<=6`, restriction of
the uniform order-`N` minimizing fibre is exactly the uniform order-`m`
minimizing fibre.  Hence the marginal KL term tends to zero as
`beta->infinity`.  If

```math
c_k={M_k\over k^{3/2}},
```

then, **when the restriction marginal remains on the smaller minimizing
fibre**, the frozen energy excess is exactly

```math
\mathcal D_{N,m}^{\infty}
=N(N-1){m\over m-1}c_m-N^2c_N
=N^2(c_m-c_N)+O(N)                                    \tag{DR.12}
```

on proportional restrictions.  Thus a temperature-uniform version of ECR
is already false whenever this excess is positive, and the frozen form is
quantitatively the desired normalized-cap comparison rather than a simpler
entropy correction.  No scalable counterexample to **fixed-temperature**
ECR appears through order six, but the data disprove the hope that marginal
KL automatically supplies the missing cancellation.

The accidental support property stops at the next order.  Exact orbit
enumeration gives

| restriction | nonminimal marginal mass | expected restricted cap | `D_infty` |
|---|---:|---:|---:|
| `7 -> 4` | `0.155556` | `4.311111` | `6.366016` |
| `8 -> 4` | `0.142857` | `4.285714` | `11.715729` |
| `8 -> 5` | `0.657143` | `5.428571` | `5.703962` |

In these cases the frozen marginal KL against the smaller canonical law
itself has a term linear in `beta`, because the restriction gives positive
mass to nonminimal children.  That term may cancel part of the energy
excess in the exact recurrence, but ECR deliberately discards it.  This is
further evidence that ECR is a strong new entropy--energy assertion, not a
free corollary of canonical marginalization.

The computation is reproduced by

```bash
.venv/bin/python computations/audit_canonical_disorder_restriction.py \
  --max-n 6 --betas 0.5 1 2 \
  --output computations/results/canonical_disorder_restriction_audit.json

.venv/bin/python computations/audit_frozen_minimizer_restrictions.py \
  --output computations/results/frozen_minimizer_restrictions_audit.json

.venv/bin/python computations/audit_canonical_disorder_root_gauge.py \
  --max-n 8 --min-parent 7 --betas 0.5 1 2 \
  --output computations/results/canonical_disorder_root_gauge_audit.json
```

## 6. Falsifier and next campaign boundary

The decisive falsifier is a sequence of comparable pairs and fixed
`beta>0` for which

```math
\beta\mathcal D_{N,m}-\mathcal S_{N,m}=\Omega(N^2).       \tag{DR.13}
```

This kills (ECR), even if the relative-entropy term in (DR.4) happens to
save the true recurrence.  A proof of (ECR) must instead come independently
from restriction incidence, cap layers, or an entropy inequality.  It may
not use the unknown limiting partition functions.

The selected one-shot campaign was therefore only:

> **Canonical-disorder entropy--restriction compensation.**  Audit
> (DR.4)--(DR.10), enumerate (DR.13) at larger orders and in the
> zero-temperature uniform-minimizer limit, and then either prove (ECR)
> from an independent cap-layer/entropy argument or construct the scalable
> falsifier (DR.13).

## 7. Outcome of the one-shot scalar audit

The campaign did not prove or falsify ECR for the actual quadratic cap.
It did prove a sharp proof-class obstruction, recorded in
`cross_order_scalar_entropy_restriction_no_go.md`.  There are hereditary,
parity-correct, switching-symmetric scalar cap systems with exact extension
fibres and edge Lipschitzness for which

```math
\beta\mathcal D_{N_j,\lfloor N_j/2\rfloor}
-\mathcal S_{N_j,\lfloor N_j/2\rfloor}
=(0.02\beta+o(1))N_j^2.                               \tag{DR.14}
```

The obstruction survives arbitrary uniformly `o(n^(3/2))` scalar
perturbations, including a parity-correct `O(n)` function of the genuine
cut cap.  It also survives an exact affine representation by the `2^n`
multiplicatively closed cut directions when their coefficient is `1/n`.
On the other hand, a flat positively homogeneous response with only a
bounded scalar range requires `exp(Theta(n^2))` queries.  Therefore the
remaining possible positive input is sharply localized: it must use the
zero-offset, unit-leading-coefficient, joint restriction incidence of the
actual cut characters.  Generic cap layers, extension counts, canonical
variationality, and Shearer contraction cannot supply a sublinear defect.

This is a criterion-four architecture boundary, not Level 6 and not an
actual-child recurrence.  The struck finite-temperature child branch stays
frozen.  Canonical-disorder ECR deserves at most one cut-incidence-specific
campaign; if that campaign cannot produce an inequality unavailable to the
countermodels, the scalar-pressure architecture must also be frozen.
