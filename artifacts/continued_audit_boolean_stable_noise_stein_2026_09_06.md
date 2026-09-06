# Stable Gaussian noise on the finite Boolean cube

Date: 2026-09-06. Independent finite-cube proof requested by the director.
This module retains the literal coherent Boolean law and does not require
moment determinacy of that law. It is separate from the special marked
field's stronger subexponential-tail argument.

The underlying Rademacher integration-by-parts and Stein framework is
classical; see Nourdin, Peccati and Reinert,
[*Stein's method and stochastic analysis of Rademacher functionals*](https://arxiv.org/abs/0810.2890).
The calculation below is self-contained and verifies the particular
stable comparison needed here. No novelty claim is made for Stein's
method or the discrete integration-by-parts method itself.

## 1. Statement

Let S have independent unbiased Boolean coordinates. Let Z_n be a centered
homogeneous Boolean Walsh polynomial of fixed degree q>=2, with variance
sigma_n² uniformly bounded. Normalize its symmetric diagonal-free tensor
K so that ||K||_F²=E Z_n². Assume every proper fixed-root flattening of K
has operator norm at most epsilon_n, where epsilon_n->0.

Let W_n=(W_n^1,...,W_n^k) be a fixed finite family of Boolean polynomials
of degrees at most d<q, whose variances are uniformly bounded. High
individual influences in W are allowed. Then for fixed real t and u,

```math
\boxed{
\mathbb E e^{itZ_n+iu\cdot W_n}
=e^{-t^2\sigma_n^2/2}\mathbb E e^{iu\cdot W_n}
+O_{q,d,k,t,u}(\epsilon_n).
} \tag{1}
```

Constants also depend on the stated variance bounds, but not on means
of W. The error is uniform on compact sets of t,u. The same conclusion
holds for d=q if every same-degree covariance E[Z_n W_n^j] is O(epsilon_n).

If the W means are bounded as well, their laws are tight. Equation (1)
then implies comparison of the actual joint law with the product of its
actual W law and an independent N(0,sigma_n²), for fixed bounded continuous
tests, in the usual ordered asymptotic sense. No Gaussian law is assigned
to W, and no assumption that its moments determine its law is used.

## 2. Exact cube integration by parts

Set

```math
\Delta_aF=\frac{F(S)-F(S^{(a)})}{2S_a}.
```

For homogeneous Z of degree q, exact Walsh orthogonality gives

```math
\mathbb E[ZF]=\frac1q\sum_a
\mathbb E[\Delta_a Z\,\Delta_a F]. \tag{2}
```

It holds for every function F on the finite cube, with no smoothness
assumption. Also sum_a E(Delta_a Z)²=q E Z². For a polynomial W of
degree at most d, its total derivative energy is at most d Var(W).

Define

```math
\Gamma_{ZZ}=\frac1q\sum_a(\Delta_a Z)^2,
\qquad
\Gamma_{ZW}=\frac1q\sum_a\Delta_aZ\,\Delta_aW.
```

Then

```math
\|\Gamma_{ZZ}-\sigma^2\|_2=O(\epsilon_n),\qquad
\|\Gamma_{ZW}\|_2=O(\epsilon_n) \tag{3}
```

under the statement's degree and covariance conditions.

## 3. Why the mixed derivative bounds follow from proper cuts

Write Z=sum_{|I|=q} z_I S_I and a homogeneous degree-e part of W as
sum_{|J|=e} w_J S_J. The exact product identity is

```math
\sum_a\Delta_aZ\Delta_aW
=\sum_{I,J}|I\cap J|\,z_Iw_J S_{I\triangle J}. \tag{4}
```

At fixed intersection size r>=1, its coefficient tensor is a fixed
degree-dependent multiple of a symmetrized r-slot contraction of the
two diagonal-free tensors, followed by projection onto distinct remaining
marked labels. The base diagonal-free tensors automatically exclude any
summed common label equaling a remaining label. Distinctness between the
two unmatched marked sets is precisely the final Hilbert projection.

When r<q,

```math
\|K_Z\mathbin{\star_r}K_W\|_F
\le\|K_Z^{(r\mid q-r)}\|_{op}\|K_W\|_F
\le\epsilon_n\|K_W\|_F. \tag{5}
```

The symmetrization and diagonal-free projection do not increase Hilbert
norm apart from fixed normalization constants. If e<q, all possible
r<=e are proper. If e=q, the sole exception r=q is the constant q E[ZW],
which is small under the stated extra covariance condition. Summing the
fixed finite homogeneous parts proves the mixed estimate in (3).
For Z against itself, r=q is exactly q sigma² and all other contractions
are proper; this proves the first estimate. Thus (3) does not rely on a
coherent-network factorization or on moment-determinacy arguments.

The one-slot cut also gives

```math
\max_a\mathbb E(\Delta_aZ)^2=O_q(\epsilon_n^2). \tag{6}
```

Indeed its tensor slice has squared norm equal to the Boolean influence
divided by q. This is where small noise influences enter the next step;
no analogous bound is imposed on W.

## 4. Characteristic functions and the discrete chain remainder

Take F=exp(itZ+iu dot W). Since
`Z(S^a)=Z-2S_a Delta_aZ` and similarly for W, the elementary exponential
Taylor bound gives

```math
\Delta_aF
=iF\left(t\Delta_aZ+\sum_j u_j\Delta_aW^j\right)+R_a,
\quad
|R_a|\le C_{t,u}
\left(|\Delta_aZ|+\sum_j|\Delta_aW^j|\right)^2. \tag{7}
```

Fixed-degree Boolean hypercontractivity bounds the L3 and L4 norms of
each derivative by a degree-dependent multiple of its L2 norm. Together
with (6) and the total derivative energies, it gives

```math
\sum_a\mathbb E|\Delta_aZ|^3=O(\epsilon_n),
\qquad
\sum_a\mathbb E|\Delta_aZ|\,|\Delta_aW|^2=O(\epsilon_n),
\qquad
\sum_a\mathbb E|\Delta_aZ|^2|\Delta_aW|=O(\epsilon_n). \tag{8}
```

For example, the middle sum is at most a constant times
`max_a sqrt(Inf_a Z) sum_a Inf_a W`. The last is bounded by
`max_a sqrt(Inf_a Z) sqrt(sum_a Inf_a Z sum_a Inf_a W)`.
Finite sums over coherent components cause only fixed constants.
Therefore the contribution of R_a in (2) is O(epsilon_n), even if a
coherent coordinate has an order-one own-spin influence.

Substitution of (3) and (7) into (2) yields

```math
\mathbb E[ZF]=it\sigma_n^2\mathbb E F+O_{t,u}(\epsilon_n).
```

For phi_n(t,u)=E F this is the differential equation

```math
\partial_t\phi_n(t,u)
=-t\sigma_n^2\phi_n(t,u)+O_{t,u}(\epsilon_n).
```

Integrating with initial condition phi_n(0,u)=E exp(iu dot W_n) proves
(1). There is no division by sigma_n, so vanishing variances are harmless.

## 5. Several fixed noise degrees and exact scope

The same proof handles a fixed finite sum Z=sum_q Z_q, with all noise
degrees at least d and the necessary equal-degree coherent covariances
small. In (2) use the inverse-generator field U=sum_q Z_q/q, so
E[ZF]=sum_a E[Delta_aU Delta_aF]. Cross-noise derivative products are
again finite contraction sums. Different noise degrees have zero total
covariance; their full contraction of the smaller tensor is proper for
the larger tensor, so (5) bounds it. Equal-degree whole contractions
give precisely Var(Z). The same maximum-influence and Taylor estimates
apply to the finite sum. A fixed finite-dimensional noise vector can be
treated by this argument for each linear combination.

This is a stable LOCAL law module. It does not prove a normalized-nuclear
cross-root covariance theorem, mixed-star ENERGY classification, a
countable-response approximation, or a new feedback iteration by itself.
In particular q<d requires additional full-noise contraction bounds and
is not covered just by the noise's proper cuts. In the marked application
the degree-three h3(G) exception needs its separately verified small
same-degree covariance with D,Y,QD; higher noise degrees satisfy q>d.

## 6. Reproducible high-influence finite-cube check

The script `computations/continued_audit_boolean_stable_noise_stein_2026_09_06.py`
uses Z=m^-1/2 sum_j product_{a=1}^5 S_{j,a} and the fixed coherent variable
W=sum_{a=1}^5 2^-a S_{1,a}. The latter has maximum influence 1/4 and a
fixed 32-atom non-Gaussian law. The noise has proper-cut norm m^-1/2,
variance one, and Gamma_ZZ=1 exactly. Its mixed Gamma L2 norm is exactly
`|| (2^-1,...,2^-5) ||_2/(5sqrt(m))`.

All characteristic functions are computed by the 32 first-block patterns
and an independent-block cosine product, not Monte Carlo. The comparison
error tends to zero at the predicted scale while retaining the literal
coherent characteristic function. This is a test of the general finite-cube
module, not an extremal-signing certificate.
