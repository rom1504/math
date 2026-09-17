# First frozen mapping: biased GS rounding of scalar responses

2026-09-17, before consulting old discrepancy artifacts.

Let `E=binom([n],2)`, `m=|E|`, `q_x(e)=x_i x_j`, and
`R(b)=max_x |<b,q_x>|` for a fractional full edge vector `b in [-1,1]^E`.
This is the edge-sum convention; the hollow symmetric matrix convention is
`x^T B x=2<b,q_x>`.

For any positive definite matrix `Q` on edge space with `Q_ee<=1`, run GS
on the columns of `Q^(1/2)` from initial coloring `b`. Its output `s` is
full sign, `E s=b`, and the sharpened Harshaw--Savje--Spielman--Zhang
analysis gives

```math
\mathbb E\exp\langle t,s-b\rangle
\le \exp\bigl(\tfrac12t^TQ^{-1}t\bigr).
```

For a desired cap `H>=R(b)`, set `d_(x,sigma)=H-sigma<b,q_x>` and
`v_x=q_x^T Q^(-1)q_x`. If

```math
\sum_{x:x_1=1}\sum_{\sigma=\pm1}
\exp\left(-\frac{d_{x,\sigma}^2}{2v_x}\right)<1,
```

there exists an output with `R(s)<=H`. Every Boolean query is retained;
the positive and negative cap have their own slack. This is a concrete
exact transport criterion, not yet a novel theorem beyond GS plus a union
bound. The point to investigate is whether near-optimal response slacks
admit a cheap ellipsoid or a random pivot-phase compensator.

The 2024 refinement actually establishes

```math
\mathbb E\exp\left(\langle\theta,V(s-b)\rangle
-\tfrac12\sum_p\|P_pv_p\|^2\|P_p\theta\|^2\right)\le1,
```

where the random mutually orthogonal pivot-phase projectors `P_p` sum to
the projector onto the input column space. Crucially, one may not replace
the random compensator by its expectation inside this exponential.

Two immediate barriers need to be included in any claimed bias refinement.

1. For a single biased sign with mean `b`, the variance is `1-b^2`, whereas
   its optimal all-lambda subGaussian proxy is
   `kappa(b)=2b/log((1+b)/(1-b))` (continuous value `1` at zero).
   The ratio diverges near an endpoint. Thus a universal subGaussian proxy
   comparable to the true covariance is impossible already in dimension one.
2. Uniform Boolean averaging gives `E_x q_x q_x^T=I_m`. Consequently
   `E_x v_x=tr(Q^-1)>=m` under `Q_ee<=1`; ellipsoid preconditioning alone
   cannot lower every all-query variance below the independent scale.

Primary inputs:

- Bansal--Dadush--Garg--Lovett, *The Gram--Schmidt Walk: A Cure for the
  Banaszczyk Blues*, Theory of Computing 15(21), 2019,
  https://theoryofcomputing.org/articles/v015a021/v015a021.pdf.
- Harshaw--Savje--Spielman--Zhang, *Balancing Covariates in Randomized
  Experiments with the Gram--Schmidt Walk Design*, JASA 2024,
  https://arxiv.org/abs/1911.03071, supplement S3.5, S8.1. The proof is for
  arbitrary starting bias and improves variance proxy 40 to exactly 1.

Status: primary GS mechanisms reconstructed; new convergence input not yet
obtained. Need compare against archive only after this frozen mapping.
