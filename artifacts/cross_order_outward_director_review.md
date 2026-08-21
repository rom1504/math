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

This is a sharp architecture boundary, not Level 6 and not an actual-child
recurrence.  Under the campaign's strict accounting it is nevertheless a
**STRIKE**, not a criterion-four reset.  The archived theorem in
`cross_order_fixed_temperature_centering_no_go.md` had already proved that
fixed-temperature scalar pressure plus universal scalar regularity cannot
yield an `o(N)` own-scale defect.  The present theorem gives a stronger and
more targeted falsifier for ECR, but it still leaves actual unit-leading cut
incidence untouched and therefore does not force a new recurrence
architecture beyond the signing-specific input already known to be needed.
The actual-child and scalar-pressure branches are frozen.

## 8. Sole post-strike campaign: rare-event restriction shadow

The outward alternative is not another scalar ECR inequality.  For the
actual quadratic cap define the cumulative low-cap layer and its entropy
density by

```math
\mathcal K_n(c)=\{A:Q(A)\le c n^{3/2}\},
\qquad
s_n(c)=L_n^{-1}\log|\mathcal K_n(c)|.                 \tag{DR.15}
```

Set `s_n(c)=-infinity` when the layer is empty.

For fixed `beta`, put `a_n=beta n^2/L_n` and choose an **exposed cap
level**

```math
c_{n,\beta}^*\in
\mathop{\rm argmax}_c\{s_n(c)-a_nc\}.                 \tag{DR.15a}
```

When there are ties, the target asks for the existence of one selection of
maximizers across all orders for which the estimates below hold uniformly;
it does not require the estimates for every maximizer.

The distance from the ground-state layer is quantitatively bounded.  If
`c_n=M_n/n^(3/2)`, optimality in (DR.15a) gives

```math
\boxed{
0\le c_{n,\beta}^*-c_n
\le {\log2-s_n(c_n)\over a_n},
\qquad a_n={2\beta n\over n-1}.}                      \tag{DR.15b}
```

Indeed,
`s_n(c*)-a_nc*>=s_n(c_n)-a_nc_n` and `s_n(c*)<=log2`.
Thus fixed temperature permits a genuine entropy buffer of at most
`(log2)/(2beta)+o(1)`, while no uniform ground-state restriction is assumed.

### Proposition DR.2 (an exposed cumulative layer is one shell)

Let `q_0<...<q_J` be the possible cap values at order `n`, put

```math
K_j=|\{A:Q(A)\le q_j\}|,
```

and let `q_(j*)/n^(3/2)=c_(n,beta)^*` be a selected maximizer in
(DR.15a).  Then

```math
\boxed{
{K_{j^*-1}\over K_{j^*}}
\le \exp\{-\beta\sqrt n(q_{j^*}-q_{j^*-1})\}
\le e^{-2\beta\sqrt n},}                              \tag{DR.15c}
```

with the left side interpreted as zero when `j*=0`.  Consequently the total
variation distance between the uniform law on the cumulative exposed layer
and the uniform law on the exact shell
`{A:Q(A)=q_(j*)}` is at most `e^(-2 beta sqrt(n))`.

**Proof.**  Maximization of
`log K_j-beta sqrt(n)q_j` at `j*` gives

```math
\log{K_{j^*}\over K_{j^*-1}}
\ge\beta\sqrt n(q_{j^*}-q_{j^*-1}).
```

Every quadratic energy, and hence every cap, is congruent to `L_n` modulo
two, so consecutive distinct cap values differ by at least two.  Finally,
the cumulative uniform law puts mass `K_(j*-1)/K_(j*)` below its top shell;
conditioning it on that shell gives the claimed total-variation identity.
`square`

Thus any restriction-tail theorem proved for the uniform exact exposed
shell implies (ERSR) with the explicit transfer

```math
\eta_{N,m,\beta}^{\rm cumulative}
\le\eta_{N,m,\beta}^{\rm shell}+e^{-2\beta\sqrt N}.   \tag{DR.15d}
```

This is the smallest currently identified parent ensemble: no full cap
histogram or mixture of low-cap layers is required.

The weakest unproved target is the following **exposed-layer
restriction-shadow statement**.  Uniformly over every integer pair
`m/N in [1/3,2/3]`, let `A` be uniform on
`K_N(c_(N,beta)^*)` and `S` a uniform `m`-set.  For each fixed `beta`,
prove some `delta_beta>0` and nonnegative errors such that

```math
\Pr\left\{Q(A[S])>
 (c_{N,\beta}^*+\epsilon_{N,m,\beta})m^{3/2}\right\}
\le\eta_{N,m,\beta},
\qquad
\epsilon_{N,m,\beta}+\eta_{N,m,\beta}
=O_\beta(N^{-\delta_\beta}).                         \tag{ERSR}
```

The finite data give a real warning rather than support.  For a uniform
exact order-eight minimizer, every four-vertex restriction has cap `4` or
`6`, so all restrictions have normalized cap at least `0.5`, while the
parent has normalized cap `10/8^(3/2)=0.4419417...`.  Thus the bad mass in
(ERSR), whenever the exposed level is the minimum, is exactly one at this
order for every
`epsilon<0.0580582...`.  This is not an asymptotic falsifier, but it makes
vanishing escaping mass the first required experiment rather than an
assumption.  The exact histogram is reproduced by
`audit_frozen_minimizer_restrictions.py`.

The finite-temperature exposed-layer audit is also adverse.  At order
eight and `beta=0.25` or `0.5`, the exposed normalized level is
`12/8^(3/2)=0.530330...`; its bad restriction masses at child orders
`3,4,5` are respectively `1`, `0.193526...`, and `0.722379...`.  At
`beta=1,2` the exposed layer is already the minimizing layer and the masses
are `1`, `1`, and `0.657143...`.  These exact values are reproduced by
`audit_exposed_layer_restrictions.py`; they are a finite warning, not a
scalable falsifier.

Passing to the exact top shell does not remove the warning.  At order eight,
`beta=0.25` or `0.5`, its bad masses for child orders four and five are
`0.194157...` and `0.723192...`, versus `0.193526...` and `0.722379...` for
the cumulative exposed layer.  This agrees with the exponentially small
comparison in Proposition DR.2.

This is an ensemble theorem about the actual cut landscape, not a scalar
axiom or a surrogate cap.  It has an immediate quantitative cross-order
consequence.  Let `Y` be the canonically relabelled restriction and let
`p<=eta_(N,m,beta)` be its actual bad probability.  The parent layer is
permutation invariant, so edge-Shearer and conditioning on the event in
(ERSR), with `c=c_(N,beta)^*`, give

```math
\begin{aligned}
\log|\mathcal K_N(c)|
&\le q^{-1}H(Y),\\
H(Y)
&\le h_2(p)
 +(1-p)\log|\mathcal K_m(c+\epsilon_{N,m,\beta})|
 +pL_m\log2.
\end{aligned}                                         \tag{DR.16}
```

For all large `N`, `eta_(N,m,beta)<=1/2`, and the right side is bounded by
the expression with `p` replaced by `eta_(N,m,beta)`.  After division by
`L_N=q^{-1}L_m`, this is

```math
s_N(c)\le s_m(c+\epsilon_{N,m,\beta})
 +\eta_{N,m,\beta}\log2
 +{h_2(\eta_{N,m,\beta})\over L_m}.                 \tag{DR.17}
```

Moreover, elementary grouping by the at most
`L_n+1` cap values yields the exact Laplace sandwich

```math
\sup_c\{s_n(c)-a_nc\}
\le\psi_n(\beta)
\le\sup_c\{s_n(c)-a_nc\}
 +{\log(L_n+1)\over L_n}.                             \tag{DR.18}
```

For fixed `beta`, the exposed levels lie in a compact `c` interval:
`s_n(c)<=log2`, while any uniform all-order bound
`M_n<=Cn^(3/2)` (the standard random-union bound suffices) gives a finite
lower competitor.  Enlarge the compact interval slightly to absorb the
shift by `epsilon_(N,m,beta)`.  Since
`a_n=2 beta n/(n-1)` and `a_N-a_m=O_beta(N^{-1})` on comparable orders,
(DR.17)--(DR.18), evaluated only at the exposed parent level, prove with
`delta'_beta=min(delta_beta,1)`,

```math
\boxed{
\psi_N(\beta)\le\psi_m(\beta)
 +O_\beta(N^{-\delta'_\beta})}                        \tag{DR.19}
```

after decreasing `delta_beta` if necessary.  The geometric restriction argument
from (DR.8), followed by the soft-minimum squeeze (DR.9), would then prove
convergence.  Thus (ERSR) has the required direct quantitative arrow.

Its scope warning remains severe.  If the exposed level happens to be
`c_(N,beta)^*=M_N/N^(3/2)`, then `eta_(N,m,beta)<1` already guarantees one
restriction with

```math
{M_m\over m^{3/2}}
\le {M_N\over N^{3/2}}+\epsilon_{N,m,\beta}.          \tag{DR.20}
```

More generally, before any tail analysis, the mere existence of one good
restriction gives the necessary condition

```math
\boxed{
{M_m\over m^{3/2}}
\le c_{N,\beta}^*+\epsilon_{N,m,\beta}.}              \tag{DR.20a}
```

Thus a fixed positive gap between the exact child optimum and the exposed
parent level on any comparable subsequence is a decisive scalable
falsifier.  This scalar check must precede attempts to estimate the full
restriction distribution.

However, ERSR is required only at one finite-temperature exposed layer for
each `beta`; it ignores every unexposed layer, its constants may deteriorate
with `beta`, and no uniform minimizer statement is assumed before taking
`N->infinity`.  It therefore has strictly fewer layer requirements than
all-layer RSR.  It does not directly yield (DR.20) with the ERSR error, nor
exhibit a good minimizer restriction, unless the exposed layer actually freezes.
Whether that distinction creates leverage is precisely the one-shot test.
The campaign must be counterexample-first: prove ERSR from genuinely new
cut-incidence input or produce a scalable exposed-layer falsifier, then
stop.  Another conditional scalar identity is not progress.

The per-signing projection boundary is now rigorous; see
[`cross_order_exposed_shell_projection_barriers.md`](cross_order_exposed_shell_projection_barriers.md).
Every scalar parent-plus-residual coordinate projection has an
`Omega(N^(3/2))` residual on comparable restrictions, and low cap—even
`o(N^(3/2))` additive near-minimality—does not improve the forced
`O(N^(3/4))` operator scale.  Consequently Hanson--Wright,
hypercontractive, and bounded-difference arguments applied one signing at
a time cannot prove ERSR.  These are proof-class barriers, not an ERSR
falsifier: joint incidence across the uniform exposed shell remains the
sole live mechanism.

One completely quantitative form of that incidence input is available.
Let `E_N(q)={A:Q(A)=q}` and `L_N(q)={A:Q(A)<=q-2}`, where `q` is the exact
exposed cap selected in Proposition DR.2.  Fix
`epsilon_N=O(N^(-delta))`, and call `(A,S)` bad when `A in E_N(q)`,
`|S|=m`, and

```math
Q(A[S])>(q/N^{3/2}+\epsilon_N)m^{3/2}.               \tag{DR.21}
```

Suppose, uniformly on comparable pairs, that every bad incidence has a
repair `R(A,S) in L_N(q)` and a label `kappa(A,S)` in a set of size
`R_N=exp(o(sqrt N))`, such that

```math
(A,S)\longmapsto(R(A,S),S,\kappa(A,S))               \tag{DR.22}
```

is injective.  If `K_j=|{A:Q(A)<=q}|` and
`r=|L_N(q)|/K_j`, double counting and (DR.15c) give

```math
\boxed{
\Pr_{A\in E_N(q),S}\{(A,S)\text{ is bad}\}
\le {R_Nr\over1-r}
\le \exp\{-2\beta\sqrt N+o(\sqrt N)\}.}             \tag{DR.23}
```

Adding the exact-shell-to-cumulative total-variation error from (DR.15d)
preserves the same exponential bound.  Thus (DR.22) implies ERSR with the
stated `epsilon_N` and exponentially small `eta_N`, and hence by (DR.19)

```math
\boxed{
\psi_N(\beta)\le\psi_m(\beta)
 +O_\beta(N^{-\min(\delta,1)}).}                    \tag{DR.24}
```

This is a conditional arrow, not a new recurrence.  It makes the next
test concrete: construct or falsify a low-multiplicity cap-lowering repair.
Another per-signing norm estimate without (DR.22) cannot use the
exponential exposed-shell gap.

In particular, (DR.22) needs no abstract transport construction.  If every
bad `(A,S)` can be repaired to cap at most `q-2` by flipping at most

```math
k_N=o(\sqrt N/\log N)                                 \tag{DR.25}
```

edges, choose one repair deterministically and use its flipped-edge set as
the label.  Given the repaired signing and that set, the parent is recovered
exactly, while

```math
R_N\le\sum_{i\le k_N}{L_N\choose i}
     =\exp\{O(k_N\log N)\}=\exp(o(\sqrt N)).          \tag{DR.26}
```

Hence this sparse cap descent implies (DR.23)--(DR.24).  It is a strictly
concrete sufficient lemma: no child histogram, Gibbs table, or target-order
optimizer is part of its statement.  It is not the archived sparse
cap-shaving target, which used `Theta(N^(3/2))` flips to seek a
`Theta(N^(3/2))` cap decrease.  Here only bad exposed-shell incidences are
repaired, the desired cap decrease is the parity step two, and the allowable
radius is `o(sqrt(N)/log N)` because shell entropy—not uniform energy
approximation—pays for the repair.

The first, cheapest repair rule is already false at finite order.  Exact
switching-quotient enumeration at `N=8` shows that, on the minimum cap-10
shell, none of the bad incidences for `m=3,4,5` can be moved into cap at
most eight by flipping one edge inside the restricted set.  On the cap-12
shell, the corresponding repairable fractions are only
`0.037406...`, `0.103486...`, and `0.122167...`.  This is not a scalable
falsifier of (DR.22), but it rules out beginning with an internal one-edge
descent map.  Even allowing the flipped edge anywhere in the parent repairs
only `0.289276...`, `0.263486...`, and `0.281773...` of the cap-12 bad
incidences.  The exact audit is reproduced by

```bash
.venv/bin/python computations/audit_exposed_shell_edge_repairs.py \
  --output computations/results/exposed_shell_edge_repair_audit.json
```

### DR.27 Sparse repair is exceptional, and the uniform version is false

The follow-up sparse-repair campaign is recorded in
[`cross_order_exposed_shell_sparse_repair_no_go.md`](cross_order_exposed_shell_sparse_repair_no_go.md).
It proves three boundaries.

First, write `E={A:Q(A)=q}` and `L={A:Q(A)<=q-2}` at a fixed-`beta`
exposed cap, and let `rho=|L|/(|E|+|L|)`.  For every radius `k`,

```math
\Pr_{A\sim U(E)}\{d_H(A,L)\le k\}
\le \left(\sum_{i\le k}{\binom{N}{2}\choose i}\right)
   {\rho\over1-\rho}.                              \tag{DR.27}
```

Since `rho<=exp(-2 beta sqrt(N))`, at
`k=gamma sqrt(N)/log N` this is at most

```math
\exp\left[-\left(2\beta-{3\gamma\over2}+o(1)\right)\sqrt N\right].
                                                               \tag{DR.28}
```

Thus almost every actual exposed-shell root is farther than every radius
allowed by (DR.25).  Sparse repair could still prove ERSR only by showing
that bad incidences concentrate on this exponentially exceptional subset.

Second, if the repair in (DR.22) has no checkable locality, algebraic, or
congestion constraint, then its existence is equivalent, fibre by fibre in
`S`, to the cardinality inequality

```math
p_{\rm bad}\le R_N{|L|\over|E|}.                    \tag{DR.29}
```

It is therefore the desired rare-tail theorem itself, not a strict
reduction.

Third, the natural extension of the sparse condition from exposed shells to
arbitrary exact quadratic-cap shells is scalably false.  For the
square-field Paley conference signing of order
`N=r^2+1`, an explicit edge-balanced law on Boolean ground states gives,
for every `k`-edge flip set `F`,

```math
Q(C^F)\ge Q(C)-{2k\over r}.                          \tag{DR.30}
```

Hence lowering the cap by two needs at least `r` flips, where
`r=sqrt(N-1)`.
At the same time, the restriction to an even union of `h~r/2` parallel
affine subfield fibres has the exact cap `rm/2`, so its normalized cap tends
to `1/sqrt(2)`, against the parent's limit `1/2`.  This is a bad comparable
restriction whose repair distance is at least project scale.  It does not
falsify DR.25 in its original exposed-shell scope or falsify ERSR: the
conference shell is not known to be fixed-temperature exposed with positive
bad-incidence mass.

An exact order-eight breadth-first search complements the asymptotic
obstruction.  Every cap-12 root is within three arbitrary edge flips of the
cap-at-most-10 layer, whereas the cap-10 layer is bottom and cannot reach
cap at most eight.  The calculation is reproduced by

```bash
.venv/bin/python computations/audit_exposed_shell_repair_distances.py \
  --output computations/results/exposed_shell_repair_distance_audit.json
```

These results supersede any all-shell version of (DR.25) as a viable
**uniform** target.  A
mass-sensitive repair theorem may discard a vanishing exceptional set, but
it must control the bad-incidence mass directly on the actual exposed
shell.  No actual-child cross-order defect improves: the bound remains
`O_beta(N)`, so this implementation receives a strike and is frozen.
