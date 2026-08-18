# The rank-one-support overlap floor need not be directional

Status: **rigorous scalable generic sharpness example**.  This note does not
use actual minimizing children.  It proves that the support-driven overlap
floor in the actual-child theorem can, at the level of rank-one channel
algebra alone, be carried entirely inside independent rows.  Thus positive
raw negative-tilt overlap cannot imply reverse-product dependence without
an additional child-specific row-irreducibility theorem.

## 1. A fixed-projective rank-one prior

Let `m,n` tend to infinity with

```math
{n\over N}\longrightarrow\alpha\in(0,1),
\qquad N=m+n,
\qquad t={\beta\over\sqrt N},                       \tag{RNS.1}
```

where `beta>0` is fixed.  Fix any `z in {+-1}^n`, draw
`X_1,...,X_m` independently and fairly, and put

```math
Q_{ij}=X_i z_j.                                      \tag{RNS.2}
```

For `B in {+-1}^{m times n}`, write

```math
S_i(B)=\sum_{j=1}^n B_{ij}z_j.
```

The exact channel pressure is

```math
L(B)=c+\log E_Xe^{t\langle B,Q\rangle}
    =c+\sum_{i=1}^m\log\cosh(tS_i(B)).               \tag{RNS.3}
```

Consequently, for every real `s`, the disorder tilt

```math
dq_s={e^{sL}\over E_Ue^{sL}}dU                      \tag{RNS.4}
```

is exactly a product over bridge rows:

```math
q_s=\bigotimes_{i=1}^m q_{s,n},
\qquad
{dq_{s,n}\over dU_n}(b)
={\cosh(t\langle b,z\rangle)^s
  \over E_{U_n}\cosh(t\langle B,z\rangle)^s}.       \tag{RNS.5}
```

In particular, the negative endpoint `q_(-lambda)` is itself an admissible
row product, so its reverse row-product projection is exactly

```math
\boxed{\mathcal I_\lambda^{\leftarrow}=0.}           \tag{RNS.6}
```

## 2. Exact cavity response

Delete edge `(i,j)` and put

```math
S_{i,-j}=\sum_{k\ne j}B_{ik}z_k.
```

All other rows are irrelevant to inference of `Q_(ij)`.  Direct Bayesian
calculation gives

```math
\boxed{
r_{ij}(B_{-(ij)})
=E[Q_{ij}\mid B_{-(ij)}]
=z_j\tanh(tS_{i,-j}).}                              \tag{RNS.7}
```

Thus every bit of the raw cavity response is row-explainable:

```math
E[Q_{ij}\mid B_{-(ij)}]
=E[Q_{ij}\mid B_{i,-j}].                            \tag{RNS.8}
```

The full overlap density at tilt `s` is therefore

```math
{1\over mn}E_{q_s}\sum_{i,j}r_{ij}^2
=E_{q_{s,n}}\tanh^2(tS_{1,-1}).                     \tag{RNS.9}
```

## 3. Uniform Gaussian limit on the negative path

Let `Z` be standard Gaussian and put `b=beta sqrt(alpha)`.  For
`s in [-lambda,0]`, define

```math
R_{b}(s)
={E[\cosh(bZ)^s\tanh^2(bZ)]
  \over E\cosh(bZ)^s}.                              \tag{RNS.10}
```

### Theorem RNS.1 (positive row-factorized overlap limit)

For every fixed `beta,lambda>0`, uniformly in `s in [-lambda,0]`,

```math
\boxed{
{1\over mn}E_{q_s}\sum_{i,j}r_{ij}^2
\longrightarrow R_b(s).}                           \tag{RNS.11}
```

Consequently

```math
\boxed{
\widehat\rho_N^-(\lambda)
\longrightarrow
\mathfrak r_{\alpha,\beta,\lambda}
:={1\over\lambda}\int_{-\lambda}^0
 {E[\cosh(\beta\sqrt\alpha Z)^s
             \tanh^2(\beta\sqrt\alpha Z)]
  \over E\cosh(\beta\sqrt\alpha Z)^s}\,ds>0.}    \tag{RNS.12}
```

This holds simultaneously with (RNS.6).

*Proof.*  Gauge `z` to the all-one vector.  Under `U_n`,

```math
{S_{1,-1}\over\sqrt n}\Rightarrow Z,
\qquad
tS_1-tS_{1,-1}=O(N^{-1/2}),                         \tag{RNS.13}
```

and hence `(tS_1,tS_(1,-1))` converges jointly to `(bZ,bZ)`.
For fixed `s<=0`, the two functions

```math
\cosh(x)^s,
\qquad
\cosh(x)^s\tanh^2(y)
```

are bounded by one.  The central limit theorem and (RNS.13) therefore give
pointwise convergence of the numerator and denominator in (RNS.9) to those
in (RNS.10).  The limiting denominator is strictly positive.

It remains only to make convergence uniform in `s`.  For `s<=0`,

```math
\left|{\partial\over\partial s}\cosh(x)^s\right|
\le\log\cosh x\le |x|.                              \tag{RNS.14}
```

The same bound holds after multiplying by `tanh^2(y)`.  The variables
`tS_1` are uniformly subgaussian, so their first absolute moments are
uniformly bounded.  Thus both finite numerator/denominator families and
their Gaussian limits are uniformly Lipschitz in `s`.  Pointwise convergence
on a finite mesh, followed by mesh refinement, is uniform on
`[-lambda,0]`.  Since the Gaussian denominator has a positive minimum on
this compact interval, the ratios converge uniformly.  This proves
(RNS.11).  Dominated integration gives (RNS.12).

Finally, `b>0`; the Gaussian integrand in the numerator is nonnegative and
strictly positive off the null event `Z=0`.  Its denominator is finite and
positive, so `R_b(s)>0` throughout the interval and the final integral is
strictly positive. `square`

## 4. Exact lesson for the actual-child theorem

The example has precisely the exact rank-one support and weak-channel
normalization used by the actual-child obstruction, yet

```text
raw negative-tilt overlap density   -> a positive constant,
reverse row-product projection      = 0.
```

The missing information is not the size of the latent support.  It is
whether the recovered planted coordinate needs other bridge rows.  At a
minimum, any directional use of the support floor must remove the
row-explainable response

```math
r_{ij}^{\rm row}=E[Q_{ij}\mid B_{i,-j}]              \tag{RNS.15}
```

and prove that a positive fraction of the floor remains in a genuinely
cross-row residual.  Equivalently, an actual-child hypothesis must exclude
asymptotic projective freezing of the opposite child sector, of which
`Y=+-z` is the exact endpoint.

That exclusion alone is necessary but not sufficient for an extensive
reverse projection.  Generic weak common-latent examples already show that
even extensive cross-row erasure information can coexist with
`I^leftarrow=O(1)`.  A sufficient directional theorem must additionally
provide an actual-channel reverse tensorization inequality, schematically

```math
\mathcal I_\lambda^{\leftarrow}
\ge c\,t^2\int_{-\lambda}^0
 \sum_{i,j}E_{q_s}
  |r_{ij}-r_{ij}^{\rm row}|^2\,ds-o(N),             \tag{RNS.16}
```

or a weaker certified coarse-row version of it.  The conjunction of

1. a positive cross-row residual after (RNS.15), and
2. a reverse tensorization/no-product-background property such as
   (RNS.16)

is the minimal evidence-backed extra structure needed to turn the raw
support floor into a directional reverse-product certificate.  Neither
property follows from support size or one-bit regularity alone.
