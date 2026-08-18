# Actual-child projective atom and collision audit

Status: **exhaustive finite child selection through order eight; numerical
Gibbs summaries**.  This experiment tests the common-sign/shared-latent
mechanism on actual contracted-temperature pressure minimizers.  It uses no
conference, Paley, or other surrogate signing.

Reproducible sources:

- [`actual_child_projective_atom_entropy.py`](actual_child_projective_atom_entropy.py)
- [`../../computations/results/actual_child_projective_atom_entropy.json`](../../computations/results/actual_child_projective_atom_entropy.json)
- [`../../computations/logs/actual_child_projective_atom_entropy.log`](../../computations/logs/actual_child_projective_atom_entropy.log)

From the repository root, the committed result is reproduced by

```bash
.venv/bin/python \
  extremal_information/experiments/actual_child_projective_atom_entropy.py \
  --orders 2 3 4 5 6 7 8 --betas 1 2 4 --mp-dps 80 \
  --output computations/results/actual_child_projective_atom_entropy.json \
  2>&1 | tee computations/logs/actual_child_projective_atom_entropy.log
```

The recorded run used Python 3.9.2 and NumPy 2.0.2 and took 47.1 seconds on
the current host.  The JSON records those versions, all selector
certificates, representative hashes, pressure gaps, class multiplicities,
and every class-pair/orientation result.

All entropies below use natural logarithms.  The computation uses balanced
parent order `N=2m` and the contracted raw temperature
`t=beta/sqrt(N)` exactly as requested.

## 1. Exact finite identities

Let `A,D` be two children, fix the relative orientation `epsilon`, and use
one representative from each projective spin class.  Put

```math
E_epsilon(x,y)=H_A(x)+epsilon H_D(y),
qquad
mathcal Z_epsilon
=sum_{[x],[y]}cosh(tE_epsilon(x,y)).             \tag{PA.1}
```

Conditioning the two exact augmented child Gibbs laws on `epsilon` and
summing the two sign representatives in each projective class gives

```math
\boxed{
P_epsilon([Q]=[xy^T])
={\cosh(tE_epsilon(x,y))\over\mathcal Z_epsilon}.} \tag{PA.2}
```

If the augmented sector `s=tau_1` is retained, then

```math
\boxed{
P_epsilon(s,[x],[y])
={e^{tsE_epsilon(x,y)}\over2\mathcal Z_epsilon}.} \tag{PA.3}
```

Conditional on `[Q]`, the two signed words `Q,-Q` have equal probability.
Consequently, if `w` denotes the projective law in (PA.2),

```math
\boxed{
\max_QP(Q)={1\over2}\max_{[Q]}w([Q]),
\quad
\sum_QP(Q)^2={1\over2}\sum_{[Q]}w([Q])^2,
\quad
H(Q)=H([Q])+\log2.}                               \tag{PA.4}
```

Thus the mass of the best literal common-sign two-word component is exactly
`max_[Q] w([Q])`.

Writing `K=max_(x,y)|E_epsilon(x,y)|`, (PA.2)--(PA.3) also give the exact
maximum-atom identities

```math
\boxed{
w_{max}={\cosh(tK)\over\mathcal Z_epsilon},
\qquad
a_{max}={e^{tK}\over2\mathcal Z_epsilon}
={w_{max}\over1+e^{-2tK}}.}                      \tag{PA.5}
```

Their collision probabilities are

```math
\boxed{
C_{\mathrm{proj}}
={\sum_{[x],[y]}\cosh^2(tE_epsilon(x,y))
  \over\mathcal Z_epsilon^2},
\qquad
C_{\mathrm{aug}}
={\sum_{[x],[y]}\cosh(2tE_epsilon(x,y))
  \over2\mathcal Z_epsilon^2}.}                 \tag{PA.6}
```

The same formulas with `E_epsilon` replaced by `H_A` give the one-child
projective law on `[x]` and augmented-projective law on `(s,[x])`.

There is also an exact shell certificate.  If `g` of the `K_0` projective
atoms have absolute energy `K`, and every other atom has absolute energy at
most `K_2<K`, then

```math
\boxed{
{1\over g+(K_0-g){\cosh(tK_2)\over\cosh(tK)}}
\le w_{max}\le {1\over g}.}                     \tag{PA.7}
```

The program verifies (PA.5) numerically for every selected class pair and
records (PA.7) from exact integer energy histograms.

## 2. Exhaustive scope

For every `m=2,...,8` and `beta in {1,2,4}`, the program:

1. enumerates all `2^((m-1)(m-2)/2)` root-gauged child signings;
2. computes exact projective integer energy histograms;
3. selects every histogram minimizing
   `log E_x cosh(beta H_A(x)/sqrt(2m))` using 80-decimal `mpmath`;
4. classifies all winners modulo switching, permutation, and coefficient
   negation;
5. evaluates every minimizing-class pair and both values of `epsilon`.

At `m=8`, this checks `2,097,152` root-gauged signings, `96` distinct
absolute-energy histograms, and `4,200` minimizing signings in two classes.
The winning histogram is unique at each tested beta.  Order nine would
require `268,435,456` root-gauged signings and a much larger exact energy
array, so order eight is the largest feasible order for the present complete
enumerator.

The signing selection and energies are exact.  Ordering the transcendental
histogram pressures is high-precision numerical, with the gap to the next
histogram recorded in the JSON.  Gibbs probabilities and information
quantities are double-precision numerical evaluations.

## 3. Pair projective concentration

For each `(m,beta)`, the table reports the largest pair-projective atom over
all minimizing class pairs and orientations, the corresponding largest
augmented-projective atom, and the smallest collision-effective projective
support

```math
N_{\mathrm{eff},2}={1\over\sum_{[Q]}P([Q])^2}.
```

The full projective rank-one alphabet has `2^(2m-2)` atoms.

| `m` | `beta` | minimizing classes | max projective two-word mass | max augmented-projective atom | min `N_eff,2` | full atoms |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 0.303388 | 0.267223 | 3.826 | 4 |
| 2 | 2 | 1 | 0.395006 | 0.387902 | 2.993 | 4 |
| 2 | 4 | 1 | 0.482337 | 0.482176 | 2.146 | 4 |
| 3 | 1 | 1 | 0.223382 | 0.221729 | 11.098 | 16 |
| 3 | 2 | 1 | 0.627275 | 0.627240 | 2.483 | 16 |
| 3 | 4 | 1 | **0.978594** | **0.978594** | **1.044** | 16 |
| 4 | 1 | 1 | 0.065334 | 0.065107 | 40.642 | 64 |
| 4 | 2 | 1 | 0.186659 | 0.186657 | 11.320 | 64 |
| 4 | 4 | 1 | 0.394624 | 0.394624 | 3.166 | 64 |
| 5 | 1 | 1 | 0.010004 | 0.009941 | 158.190 | 256 |
| 5 | 2 | 1 | 0.016471 | 0.016470 | 72.581 | 256 |
| 5 | 4 | 1 | 0.019697 | 0.019697 | 51.544 | 256 |
| 6 | 1 | 1 | 0.002995 | 0.002985 | 629.694 | 1,024 |
| 6 | 2 | 1 | 0.005820 | 0.005820 | 301.688 | 1,024 |
| 6 | 4 | 1 | 0.010221 | 0.010221 | 128.684 | 1,024 |
| 7 | 1 | 1 | 0.003849 | 0.003849 | 1,377.956 | 4,096 |
| 7 | 2 | 1 | 0.020638 | 0.020638 | 189.585 | 4,096 |
| 7 | 4 | 1 | **0.080753** | **0.080753** | **16.953** | 4,096 |
| 8 | 1 | 2 | 0.000959 | 0.000959 | 4,871.339 | 16,384 |
| 8 | 2 | 2 | 0.004507 | 0.004507 | 779.737 | 16,384 |
| 8 | 4 | 2 | **0.015837** | **0.015837** | **113.697** | 16,384 |

At `m=8,beta=4`, all eight class-pair/orientation records agree up to
floating-point error.  More detail for that largest case is

```math
\begin{aligned}
\max_{[Q]}P([Q])&=0.0158372167138032,\\
\max_QP(Q)&=0.00791860835690158,\\
C_{\mathrm{proj}}&=0.00879533023790067,\\
H_2([Q])&=4.73353435312907,\\
N_{\mathrm{eff},2}&=113.696697332730,\\
H([Q])&=5.52418379323915.
\end{aligned}                                     \tag{PA.8}
```

For comparison, the full projective entropy is
`log(16384)=9.70406052783923`.  Thus the collision-effective support is only
about `0.694%` of the full alphabet, even though no single two-word component
dominates.

## 4. One-child concentration

The one-child data explain the remaining low effective pair support.  At the
largest orders:

| `m` | `beta` | max child projective atom | max child augmented-projective atom | min child `N_eff,2` | full child projective atoms |
|---:|---:|---:|---:|---:|---:|
| 5 | 1 | 0.076120 | 0.070502 | 14.827 | 16 |
| 5 | 2 | 0.091324 | 0.090748 | 11.813 | 16 |
| 5 | 4 | 0.099244 | 0.099240 | 10.152 | 16 |
| 6 | 1 | 0.040789 | 0.038635 | 30.306 | 32 |
| 6 | 2 | 0.054112 | 0.053944 | 24.222 | 32 |
| 6 | 4 | 0.071488 | 0.071487 | 16.042 | 32 |
| 7 | 1 | 0.044225 | 0.043868 | 47.758 | 64 |
| 7 | 2 | 0.101647 | 0.101640 | 20.784 | 64 |
| 7 | 4 | 0.215107 | 0.215107 | 6.818 | 64 |
| 8 | 1 | 0.022041 | 0.021893 | 91.540 | 128 |
| 8 | 2 | 0.047472 | 0.047470 | 39.449 | 128 |
| 8 | 4 | 0.088987 | 0.088987 | 15.080 | 128 |

At `m=8,beta=4`, each child therefore has only about `15.1`
collision-effective projective states out of `128`.  The pair law is not a
product of the child cosh marginals--the common augmented sector couples
their energy signs--but the child concentration supplies a concrete finite
source for its small effective support.

The program also reconstructs every pair law independently from the exact
sector factorization

```math
\pi_s^{(\epsilon)}\,
\bar\mu_{A,s}\otimes\bar\mu_{D,\epsilon s}.
```

Across all orders, temperatures, class pairs, and orientations, the largest
coordinate residual between this reconstruction and (PA.2)--(PA.3) was
`1.666e-16`.  This is a direct numerical guard against a normalization or
orientation error in the reported projective laws.

## 5. Comparison with the rigorous actual-child entropy floor

The subsequently proved Hamming-sphere theorem in
[`../audits/actual_child_sector_min_entropy.md`](../audits/actual_child_sector_min_entropy.md)
now gives, for comparable actual minimizing children,

```math
\begin{aligned}
-{1\over N}\log\max_QP(Q)&\ge\eta_\beta-o(1),\\
{1\over N}H_2(Q)&\ge\eta_\beta-o(1),\\
\eta_\beta
&=\sup_{0<q<1/2}
 \{h(q)-4(\log2+\beta^2/4)q(1-q)\}>0.          \tag{PA.9}
\end{aligned}
```

Here `Q` is signed, whereas the table uses `[Q]`; exactly
`max_Q P(Q)=w_max/2` and `H_2(Q)=H_2([Q])+log 2`.  The order-eight observed
rates and the explicit asymptotic constants in (PA.9) are:

| `beta` | explicit `e^(-beta^2)/16` | optimized `eta_beta` | observed signed min-entropy rate | observed signed `H_2` rate |
|---:|---:|---:|---:|---:|
| 1 | 0.0229925 | 0.0250484 | 0.477693 | 0.574017 |
| 2 | 0.00114473 | 0.00115307 | 0.380958 | 0.459506 |
| 4 | 7.03345e-9 | 7.03345e-9 | 0.302409 | 0.339168 |

The finite rates are far above this deliberately uniform bound, especially
at large `beta`.  This comparison is not an estimate of their limit: the
theorem is asymptotic and the exhaustive sequence is short and nonmonotone.
Its rigorous content is that actual-child projective atom and collision
diffuseness *does* hold exponentially at every fixed `beta`; the computation
only measures how conservative the proved exponent is at the accessible
orders.

## 6. Research judgment

The finite evidence makes two different statements.

1. **The literal two-word mechanism is not algebraically excluded by actual
   minimizing children.**  At `m=3,beta=4`, one projective rank-one word and
   its two signed representatives carry `97.86%` of the law.  Thus actual
   thermal minimization plus the exact augmented Gibbs construction alone do
   not forbid the generic common-sign geometry at finite order.
2. **It does not persist literally at the largest exhaustive orders.**  At
   `m=8,beta=4`, the best two-word mass is only `1.58%`; the exact top
   absolute shell contains `32` projective pair atoms, already forcing the
   upper bound `w_max<=1/32` from (PA.7).
3. **A softer shared-latent concentration remains substantial at these
   finite orders.**  The same
   order-eight law uses only about `114` collision-effective projective words
   out of `16,384`, and order seven at `beta=4` uses about `17` out of
   `4,096`.  Bounded conditional row complexity therefore remains compatible
   with a strongly compressed, common-sector multiword law.

The sequence is highly nonmonotone (`m=7` is much more concentrated than
`m=8`), so the data alone imply nothing asymptotic.  The Hamming-sphere
theorem now rigorously excludes persistent atoms and subexponential latent
catalogues.  What remains compatible with both theorem and experiment is
**diffuse exponential-rate retuning**: a posterior can reorganize mass across
exponentially many individually tiny rank-one words.  Controlling that
distributed retuning, rather than proving atom or collision diffuseness, is
the remaining optimizer-specific obligation.  The present experiment
supplies finite targets for that question and no asymptotic claim.
