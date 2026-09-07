# Heterogeneous fibre retention cannot improve the scalar E certificate

Date: 2026-09-07. Status: analytic derivation, pending independent audit.
Scope: the scalar direct-E weave certificate, not arbitrary actual signings.

## 1. Source convexity of E

The exact rate-distortion dual is

```math
J_\lambda(\nu)
 =-\sup_\pi\int\log\left(\int e^{-\lambda(x-y)^2}\,d\pi(y)\right)d\nu(x).
```

Consequently

```math
E_t(\nu)=\sup_{0<\lambda\le t,\pi}
\left[c_t(\lambda)+
 \int\log\left(\int e^{-\lambda(x-y)^2}\,d\pi(y)\right)d\nu(x)\right].   (1)
```

For each fixed lambda,pi the bracket is affine in the source law. Thus
E_t is convex in the source law. For the finite source alphabets below all
integrals are finite; no interchange involving divergent entropy is needed.
This is convexity in the SOURCE, not a claim that an infimum over optimizing
signings preserves convexity of a coupling or temperature parameter.

## 2. Actual unequal-size fibres and their certified exponent

Use m fibres with a common Hadamard column count m, but retain k_i rows in
fibre i. Let p_i=k_i/m and pbar=m^(-1)sum_i p_i, so the total output order
is N=m^2 pbar. For a spin vector x_i on its retained rows put
h_i=H_i[T_i,:]^T x_i. Orthogonality gives

```math
\sum_i\|h_i\|^2=m\sum_i k_i=mN.
```

The rank-one weave's exact defect is still
D_sigma=2(mN-sigma x^T W x). Normalize every row by sqrt(m pbar), not by
its individual sqrt(k_i). Its input law for the row Bellman estimate is

```math
\mu_i=(1-p_i)\delta_0+{p_i\over2}
 (\delta_{-1/\sqrt{\bar p}}+\delta_{1/\sqrt{\bar p}}).       (2)
```

It has second moment p_i/pbar. The entropy cost for its row-spin sum is
p_i log2 per Hadamard column. Averaging over independent fibres therefore
produces the scalar sufficient upper bound

```math
C_{\rm het}(t;\{p_i\})
 ={t+\bar p\log2+m^{-1}\sum_i E_t(\mu_i)\over2t\sqrt{\bar p}}.   (3)
```

For a fixed finite list of retention types, (3) follows from the same
direct-E realization proof: finite alphabets and recursion depth remain
fixed, all row moments are at most 1/pbar, independent row factors multiply,
and the diagonal deletion costs O(N). Type multiplicities and retained
row counts may be rounded by O(1); their total is o(m^2). The same dense
Hadamard order sequence fills all large output orders. There is no claim
here about a growing number of types without uniform estimates.

## 3. Homogeneous retention is optimal for this certificate

The average of (2) is EXACTLY the normalized ternary law nu_pbar. Convexity
(1) gives

```math
{1\over m}\sum_i E_t(\mu_i)\ge E_t(\nu_{\bar p}),
\qquad
C_{\rm het}(t;\{p_i\})\ge C(\bar p,t).                   (4)
```

The same conclusion holds for any limiting finite probability mixture of
retention types. Varying their row dimensions therefore cannot improve this
scalar certificate at fixed average retention. This answers one natural
multi-type closure proposal without introducing a new unknown variational
quantity.

One must not normalize each fibre separately and then average
E_t(nu_pi): that changes the defect kernel. The shared normalization in
(2) is essential to both the actual weave identity and Jensen's argument.

Equation (4) concerns the bound produced by the mean-E argument. It does
not assert that every actual heterogeneous weave has cap at least C(pbar,t),
and does not exclude an improved analysis retaining joint fibre correlations.
In particular, it is not a proof that the scalar construction exhausts
original near-minimizers or that M_n/n^(3/2) converges.
