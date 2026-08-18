# Actual-minimizer radial data do not control the ES path

Status: **rigorous low-channel theorem and exact actual-minimizer
falsifier**, plus a complete finite numerical check at the physical channel
amplitude.  This note combines the optimizer-specific order-eight witness
from FC.5 with the actual-child overlap tangent RA.2 and the canonical
interaction path ES.21--ES.29.

The result is deliberately a no-go for a strict coarsening of FC.8.  It does
not use conference matrices or a generic row law.  It shows that all radial
consequences of actual minimization, even the complete pressure/entropy
curve, fail to determine the first nonzero row-dependence coefficient of the
actual child-induced law.

## 1. The low-channel canonical interaction

Fix finite children `A,D`, their internal raw temperature `t`, an
orientation `epsilon`, and inverse-disorder parameter `lambda`.  Let
`mu_epsilon` be the actual augmented child law in RA.14 and put

```math
Q_{ij}=\tau X_iY_j,
\qquad
\Gamma_{ik;j\ell}
=E_{\mu_\epsilon}X_iX_kY_jY_\ell.                  \tag{RP.1}
```

For a separate bridge amplitude `u`, write `rho=tanh u`, let `p_u` be the
forward bridge likelihood, let `r_u` be the product of its negative-escort
row marginals as in CC.5, and let `q_(u,s)` be the canonical interaction
path

```math
{dq_{u,s}\over dr_u}\propto e^{-s h_u},
\qquad 0\le s\le\lambda.
```

Define the two exact ES.28 contributions

```math
\begin{aligned}
T_u&=\lambda\int_0^\lambda
 {\operatorname {TC}_{\rm row}(q_{u,s})\over s^2}\,ds,\\
R_u&=\lambda\int_0^\lambda
 {\sum_iD(q_{u,s,i}\Vert r_{u,i})\over s^2}\,ds.
\end{aligned}                                      \tag{RP.2}
```

Thus the canonical error is exactly `J_u=T_u+R_u`.

**Theorem RP.1 (actual overlap curvature of the canonical path).**  At
fixed finite child orders and fixed `t,lambda`, as `u` tends to zero,

```math
\boxed{
\begin{aligned}
\mathcal J_u
 &= {\lambda^2\rho^4\over2}
    \sum_{i<k}\sum_{j,\ell}\Gamma_{ik;j\ell}^2
    +O(\rho^6),\\
T_u
 &= {\lambda^2\rho^4\over2}
    \sum_{i<k}\sum_{j,\ell}\Gamma_{ik;j\ell}^2
    +O(\rho^6),\\
R_u&=O(\rho^8).
\end{aligned}}                                      \tag{RP.3}
```

In particular, the first nonzero canonical interaction is genuine row total
correlation.  Moving the one-row factors cannot absorb it to this order.

*Proof.*  Relative to the fair bridge law, the binary channel has the exact
Fourier form

```math
p_u(B)=E_{\mu_\epsilon}
       \prod_{i,j}(1+\rho B_{ij}Q_{ij}).             \tag{RP.4}
```

Global flip symmetry of either child kills every odd Fourier level.  At
level two, subtracting the logarithms of the exact row marginals deletes
precisely the pairs contained in one row.  Hence, uniformly on the fixed
finite bridge cube,

```math
h_u(B)=\rho^2H_2(B)+O(\rho^4),
\qquad
H_2(B)=\sum_{i<k}\sum_{j,\ell}
 \Gamma_{ik;j\ell}B_{ij}B_{k\ell}.                 \tag{RP.5}
```

Also `r_u=U+O(rho^2)`.  The displayed Walsh characters are distinct and
orthonormal, so

```math
\operatorname {Var}_{r_u}(h_u)
=\rho^4\sum_{i<k,j,\ell}\Gamma_{ik;j\ell}^2
 +O(\rho^6).                                       \tag{RP.6}
```

The centered cumulant identity IC.7 now gives the first line of RP.3.

For the marginal statement, every term of `H_2` contains exactly one bit
from each of two different rows.  Every canonical row factor is centrally
symmetric.  Thus integrating `H_2` over all rows but one is identically
zero.  Exponential tilting consequently changes each one-row marginal only
at order `rho^4`, uniformly for `0<=s<=lambda`; its KL from `r_(u,i)` is
`O(s^2rho^8)`.  This proves the third line after division by `s^2` and
integration.  ES.28 then gives the second line. `square`

The theorem is compatible with RA.2: its coefficient is exactly the
cross-row ANOVA coefficient `K_cross`.  RP.3 adds that this same coefficient
is the leading curvature of the finite-`lambda` canonical path and locates
it on the total-correlation, rather than marginal-retuning, side of ES.29.

## 2. An exact actual-minimizer falsifier

Let `A_0,A_1` be the two certified order-eight minimizer classes displayed
in FC.5.  They have the identical signed energy histogram

```math
\#\{H=-10,-8,\ldots,10\}
=(4,10,12,16,16,12,16,16,12,10,4),                 \tag{RP.7}
```

with the entries read in increasing energy order.  Consequently they have
identical pressure and entropy at every `t`, and identical values of every
homogeneous or fixed-size flip average FC.10--FC.21.  The exhaustive
order-eight classification and FC.22 prove that both are actual pressure
minimizers for every `t>=3`.

Use the unique order-two signing as the left child.  Its sector correlation
is

```math
E[X_1X_2\mid\tau=s]=s\tanh t.
```

For an order-eight signing `D`, define its signed tangent matrix

```math
S_D(t)
={E_Y[YY^{\mathsf T}\sinh(tH_D(Y))]
  \over E_Y\cosh(tH_D(Y))}.                         \tag{RP.8}
```

Directly from the sector decomposition RA.20,

```math
\Gamma_{12;j\ell}^{(\epsilon)}
=\epsilon\tanh(t)S_D(t)_{j\ell}.                   \tag{RP.9}
```

At zero temperature only the eight projective states with `|H|=10`
survive.  Exact integer summation gives

```math
\boxed{
\lim_{t\to\infty}\operatorname {Tr}S_{A_0}(t)^2=14,
\qquad
\lim_{t\to\infty}\operatorname {Tr}S_{A_1}(t)^2=10.}          \tag{RP.10}
```

The two integer numerators, before division by the eight ground states, are
stored in the reproducible result file.  Their squared Frobenius norms are
`896` and `640`; division by `8^2` gives RP.10.

Combining RP.3, RP.9, and RP.10 proves the no-go:

```math
\boxed{
\text{complete radial FC data do not determine even the leading}
\quad \rho^{-4}\mathcal J_u
\quad\text{or weighted row-TC mass.}}               \tag{RP.11}
```

Indeed, along actual pressure minimizers their limiting `rho^4`
coefficients are `7lambda^2` and `5lambda^2`.  This is stronger than a
failure of a few pressure moments: the two children have the same entire
radial pressure function.  It does not obstruct a genuinely nonradial
coarsening of FC.8.

There is also an exact separation on the **physical diagonal** `u=t`, not
only at the low-channel tangent.  Fix orientation `epsilon=-1`.  For a
bridge `B`, define its tropical
parent rate in the projective gauge by

```math
\kappa_D(B)=\max_{x_1=y_1=1}
 \left\{|H_A(x)-H_D(y)|+|x^{\mathsf T}By|\right\}.  \tag{RP.12}
```

The sum rather than `|H_A-H_D+x^TBy|` is forced by the exact projective
identity `cosh(t(H_A-H_D))cosh(tx^TBy)`.  Put

```math
\kappa_i(b_i)=\max_{B_{-i}}\kappa_D(B).
```

For `lambda=1`, exact tropical enumeration gives

```math
\boxed{
\lim_{t\to\infty}{\mathcal J_t(A_0)\over t}={25\over12},
\qquad
\lim_{t\to\infty}{\mathcal J_t(A_1)\over t}={9\over5}.}        \tag{RP.13}
```

*Proof.*  Each finite bridge pressure has an expansion
`L_t(B)=t kappa_D(B)+log c_D(B)+o(1)` with a positive integer-multiple
leading coefficient `c_D(B)`.  A forward row marginal is dominated by
`kappa_i(b_i)`, so its negative escort concentrates on the row words
minimizing `kappa_i`; within that support its weight is proportional to the
inverse leading coefficient.  The endpoint escort concentrates on bridges
minimizing `kappa_D`.  Therefore

```math
\lim_{t\to\infty}{\mathcal J_t\over t}
=E_{r_\infty}\kappa_D-\min_B\kappa_D(B).             \tag{RP.14}
```

All quantities in RP.14 are finite integer tables.  For `A_0`, both row
minimax rates are `21`, each tropical row support has size `24`, the global
minimum is `15`, and the inverse-coefficient product law gives
`E_(r_infty)kappa=205/12`.  For `A_1`, the corresponding support sizes are
`8`, the same rates are `21` and `15`, and the expectation is `84/5`.
Subtraction proves RP.13.  The exact active multiplicities and rational
calculation are recorded by the experiment. `square`

Thus identical complete radial optimizer data fail to determine both the
infinitesimal ES curvature and the leading canonical error along an actual
physical-temperature ray.  RP.13 does not by itself split that physical
leading error between the two ES.42 alternatives.

## 3. Physical-channel finite check

At `t=u=3`, equivalently `beta=3sqrt(10)` for the `2+8` split, complete
enumeration of all `2^16` bridges gives:

| right child | `J/N` | integrated TC share | marginal-retuning share |
|---|---:|---:|---:|
| `A_0` | `.443573173` | `.816813379` | `.183186621` |
| `A_1` | `.523848938` | `.474757829` | `.525242171` |

Both orientations agree.  These values are numerical rather than interval
certificates, but the separation is large and the ES.28 identity is checked
to floating accuracy.  They show that the exact tangent separation persists
at the physical channel and, more sharply, that identical radial optimizer
data do not even predict which of TC or row retuning carries most of the
canonical mismatch.

The computation is reproduced by
[`../experiments/actual_child_radial_path_curvature_falsifier.py`](../experiments/actual_child_radial_path_curvature_falsifier.py),
with output in
[`../../computations/results/actual_child_radial_path_curvature_falsifier.json`](../../computations/results/actual_child_radial_path_curvature_falsifier.json).

## 4. Consequence for the optimizer-specific attack

This is a substantive strike against every proposed consequence of FC.8
which retains only scalar pressure/entropy information, homogeneous flip
noise, fixed-size flip noise, or their full radial transforms.  Those data
can be identical on actual minimizing children while the ES path has
different leading dependence curvature and a different physical TC/retuning
split.

The exact surviving obligation is nonradial: find a compressed statistic of
the sector overlap tensor in RP.1 whose control at zero bridge field extends
uniformly to `u=t=beta/sqrt(N)` and fixed `lambda`, or prove that such a
uniform transport requires the higher connected overlap hierarchy.  The
full labelled tensor is not automatically admissible: its exact field-grid
oracle can reconstruct the child response landscape by CC.12.
