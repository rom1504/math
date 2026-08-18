# Exact smoothing derivative for the capped actual-child overlap

Status: **rigorous scalar interpolation identity and route-specific
minimality blockade**.  The capped moderate-tilt overlap from MT.3--MT.4 is
exactly a first derivative of one Gaussian smoothing curve based at the
actual spiked bridge law.  This gives a lower-information formulation of the
remaining lemma.  It also shows why the currently available consequences of
contracted-temperature child minimality do not bound it: at replica number
one the overlap cancels identically, while away from one the needed result is
an `o(N)` derivative estimate and minimality supplies only `O(N)` endpoint
bounds for a different, zero-cross-edge variational problem.

The Gaussian variable below is an infinitesimal interpolation device around
the actual sign bridge.  The base law, pressure, children, cap, and tilted
expectations are the actual ones; no Gaussian or conference child replaces
them.

## 1. Setup and retained-tail estimate

Let `A,D` be actual contracted-temperature minimizing children of orders
`m,n`, put `N=m+n` and `t=beta/sqrt(N)`, and extend the bridge pressure to
real bridge fields `h` by

```math
 L(h)=\log E_{x,z}\cosh\left(t\{H_A(x)+\epsilon H_D(z)
                                  +x^{\mathsf T}hz\}\right).    \tag{CI.1}
```

Let `B` have the actual spiked marginal law `mu_y^(otimes m)`, let `G` have
iid standard Gaussian entries, and let `G` be independent of `B`.  Fix a
finite tilt window `S`, choose an integer `k>S`, and choose a cap

```math
 T=CN                                                        \tag{CI.2}
```

with `C` sufficiently large in terms of `beta,k,S`.  Write
`L_T=L wedge T` and

```math
 {d\Pi_{s,T}\over d\mu_y^{\otimes m}}(B)
 ={e^{sL_T(B)}\over E e^{sL_T(B)}}.               \tag{CI.3}
```

The fixed-replica proof of MT.1, with
`||e_y||_2^2<=3/2` in place of `||q_(v,y)||_2^2<=5/2`, gives

```math
 \log E_{\mu_y^{\otimes m}}e^{kL}
 \le {m\over2}\log(3/2)+mnk^2t^2
      +k\{p_A(t)+p_D(t)+\log2\}
 \le C_{\beta,k}N.                               \tag{CI.4}
```

Consequently, if `C>C_(beta,k)/(k-S)`, then uniformly for `|s|<=S`,

```math
 \boxed{
 \Pi_{s,T}\{L\ge T\}\le e^{-c_{\beta,k,S,C}N}.} \tag{CI.5}
```

Indeed, for `0<=s<=S` the numerator of the capped tilted mass is at most
`e^(sT)P(L>=T)`, the denominator is at least one, and Markov's inequality
with (CI.4) applies.  For `-S<=s<0`, the denominator is at least `e^(sT)`,
so the tilted tail is at most the untilted tail.  This is the compact-tilt
version of response-specific truncation.

## 2. Exact Gaussian smoothing derivative

For `s ne 0`, define the scalar curve

```math
 \mathscr F_{s,T}(u)
 ={1\over s}\log E_{B,G}
  \exp\left\{s\min\bigl(L(B+\sqrt uG),T\bigr)\right\},
 \qquad u\ge0.                                    \tag{CI.6}
```

At `s=0`, use its continuous extension
`F_(0,T)(u)=E_(B,G)min(L(B+sqrt(u)G),T)`.

Choose `T` outside the finite set `{L(B):B in {+-1}^(m times n)}`.  This
does not affect CI.5 or the response truncation theorem and avoids a cap
kink at `u=0`.

Let `nu_B` be the Gibbs law associated with (CI.1), and put

```math
 m_{ij}(B)=E_{\nu_B}[\tau X_iZ_j].                \tag{CI.7}
```

**Theorem CI.1 (capped overlap derivative).**  For every fixed
`0<|s|<=S`,

```math
 \boxed{
 \mathscr F_{s,T}'(0+)
 ={t^2\over2}E_{\Pi_{s,T}}
 \left[
  1_{\{L<T\}}
  \left\{mn+(s-1)\sum_{i,j}m_{ij}(B)^2\right\}
 \right].}                                       \tag{CI.8}
```

Therefore, by (CI.5),

```math
 \boxed{
 {2\mathscr F_{s,T}'(0+)\over t^2mn}
 =1+(s-1)
 E_{\Pi_{s,T},\nu_B^{\otimes2}}
   [\tau^1\tau^2R_XR_Z]
 +O_{\beta,S}(e^{-cN}).}                         \tag{CI.9}
```

*Proof.*  The heat-semigroup derivative at zero is

```math
 {d\over du}E_G f(B+\sqrt uG)\big|_{u=0+}
 ={1\over2}\Delta_hf(B).                         \tag{CI.10}
```

At a bridge with `L(B)<T`, elementary Gibbs differentiation gives

```math
 \partial_{h_{ij}}L=t m_{ij},
 \qquad
 \partial_{h_{ij}}^2L=t^2(1-m_{ij}^2).           \tag{CI.11}
```

Hence

```math
 \Delta_he^{sL}
 =st^2e^{sL}
 \left\{mn+(s-1)\sum_{i,j}m_{ij}^2\right\}.     \tag{CI.12}
```

At a bridge with `L(B)>T`, the capped function is locally constant and its
Laplacian is zero.  Average (CI.12), divide by the value and by `s`, and use
(CI.10), proving (CI.8).  Finally

```math
 {1\over mn}\sum_{i,j}m_{ij}^2
 =E_{\nu_B^{\otimes2}}[\tau^1\tau^2R_XR_Z],      \tag{CI.13}
```

and CI.5 removes the retained-tail indicator at exponentially small error.
This proves (CI.9). `square`

The edge-cavity squares used in `rho_N(S)` differ from the full-Gibbs
squares in CI.13 by `O(t)=O(N^(-1/2))`, uniformly in the capped tilt, by
FI.17.  Thus, for every fixed `s ne 1`, the desired statement

```math
 \rho_N(s)=o(1)                                   \tag{CI.14}
```

is equivalent to the scalar derivative rigidity

```math
 \boxed{
 \mathscr F_{s,T}'(0+)={t^2mn\over2}+o(N).}       \tag{CI.15}
```

This is formally one scalar derivative rather than the complete bridge
pressure or Gibbs table.  No finite-precision generation theorem is proved,
so this analytic reduction is not yet an operational information bound.

At `s=1`, the overlap is recovered by a mixed derivative rather than by the
value of the smoothing derivative:

```math
 \boxed{
 E_{\Pi_{1,T},\nu_B^{\otimes2}}[\tau^1\tau^2R_XR_Z]
 ={2\over t^2mn}
  \partial_s\mathscr F_{s,T}'(0+)\big|_{s=1}
 +O(e^{-cN}).}                                    \tag{CI.15a}
```

This follows by differentiating CI.9; the fixed-replica tail estimate
justifies differentiation uniformly near one.  Uniform overlap control on
a window containing `s=1` therefore requires secant rigidity at the
corresponding precision, not merely the unscaled estimate CI.15.

## 3. The replica-one cancellation and the minimality sign mismatch

Taking `s=1` in the limiting form of CI.8 gives

```math
 \boxed{
 \mathscr F_{1,T}'(0+)={t^2mn\over2}+O(Ne^{-cN}),}              \tag{CI.16}
```

independently of the children and independently of their overlap.  Thus the
one-replica annealed pressure estimate used in every current consequence of
child minimality is exactly blind to the new observable.  One must control
a genuine replicated derivative away from `s=1`.

Contracted-temperature minimality currently gives

```math
 p_A(t)\le {m\choose2}\log\cosh t,
 \qquad p_D(t)\le {n\choose2}\log\cosh t,         \tag{CI.17}
```

and the discrete internal-edge flip inequalities at zero external field.
Together with CI.4 these imply only

```math
 |\mathscr F_{s,T}(u)|=O_{\beta,S,u}(N)           \tag{CI.18}
```

on fixed smoothing intervals.  The baseline derivative in CI.15 is itself
`Theta(N)`, and a fixed positive overlap changes it by another `Theta(N)`.
Therefore an `O(N)` endpoint or secant estimate has exactly the wrong
precision to imply CI.15.

There is also no admissible variational comparison hidden in (CI.6):
Gaussian smoothing changes the **cross edges** of the parent field.
Child minimality compares discrete internal signings at orders `m` and `n`;
it neither declares `B+sqrt(u)G` an admissible child perturbation nor
compares its replicated pressure with the zero-cross-edge child values.
If one instead smooths internal child edges, the analogue of CI.9 contains
the internal replica overlap, and converting it to the bridge product
overlap is precisely an unproved synchronization statement.

This is a route-specific blockade to the presently available direct
interpolation, not a theorem that actual minimizers cannot satisfy CI.15 or
that a new synchronization argument cannot derive it.  It proves that the
existing zero-bridge minimality value, annealed bound, and internal flip
comparisons do not enter the exact interpolation directly at the precision
or on the edge set required.  The actual
order-eight pressure-minimizer collision in RP.7--RP.11 gives a finite sharp
warning: two minimizer classes have identical complete radial pressure and
fixed-size flip data but different first bridge-interaction curvature.
That collision rules out repairing the mismatch by retaining more scalar
zero-bridge derivatives, while not ruling out a new asymptotic
synchronization theorem.

## 4. Revised missing lemma

For the spiked response, the remaining optimizer-specific statement can now
be written without a pressure table or all-tilt MGF:

> **Capped replicated smoothing rigidity.**  Uniformly in the declared
> spike direction and on the moderate tilt window needed by MT.4,
>
> ```math
> \sup_{s\ne1}
> {\left|2\mathscr F_{s,T}'(0+)/(t^2mn)-1\right|\over|s-1|}
> =o(1),
> \qquad
> {2\over t^2mn}\partial_s\mathscr F_{s,T}'(0+)\big|_{s=1}=o(1),
> ```
>
> with a quantitative rate on the corresponding growing window.

By CI.9 this is exactly the capped moderate-overlap lemma, not another
reformulation of the full bridge optimization.  No present consequence of
contracted-temperature child minimality proves it.  A valid next move would
need a replicated external-field comparison or a synchronization theorem
between internal and cross-edge overlaps; repeating annealed or zero-field
minimality interpolation cannot change the scale.
