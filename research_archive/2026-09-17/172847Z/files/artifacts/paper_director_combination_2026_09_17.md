# Combining fixed-law comparison, localization, and exact sign rounding

2026-09-17. LIVE DERIVATION. Unless individually marked reconstructed or
audited, statements below are candidates under independent verification.
This is not a convergence proof or a new numerical signing bound.

## Primary inputs selected by mechanism

1. El Alaoui--Montanari, *An Information-Theoretic View of Stochastic
   Localization*, [arXiv2109.00709v2](https://arxiv.org/html/2109.00709v2):
   exact posterior covariance dissipation and information cost.
2. Bansal--Dadush--Garg--Lovett, *The Gram--Schmidt Walk*, with the sharper
   biased-start proof of Harshaw--Savje--Spielman--Zhang,
   [arXiv1911.03071](https://arxiv.org/abs/1911.03071): actual vertex-valued
   rounding and subGaussian imbalance. The latter constant is being audited.
3. Bednorz--Latala, *On the boundedness of Bernoulli processes* (2014),
   followed by Liu--Zadik,
   [*A Bayesian Proof of the Bernoulli Theorem*,2608.11031v1](https://arxiv.org/html/2608.11031v1).
   The new prescribed-law theorem is materially stronger than a set-width
   theorem for retaining arbitrary energy deficits.
4. Supporting bridge discovered during reconstruction: van Handel,
   [*On the subgaussian comparison theorem*,2512.18588v2](https://arxiv.org/html/2512.18588v2).
   Its convex-order/fixed-law strengthening preserves affine offsets.

The fourth paper is included because it provides a precise missing link,
not to expand a literature catalogue. Main-agent reconstruction has read
the complete Liu--Zadik proof, including the scalar harmonic kernel,
product data processing, multiscale clipping and exact-type lifting.
Decisive steps of van Handel's short proof are also reconstructed below.

## 1. An offset-preserving decomposition candidate

For a finite index set T with vectors t in R^d and offsets h_t, put

```math
\Psi_Z(T,h)=\mathbb E\max_{t\in T}\{\langle Z,t\rangle+h_t\}.
```

For a prescribed law mu on T, define B(mu) as the maximum expected
inner product of a mu-distributed index and an independently-marginal
Rademacher vector; the coupling itself is unrestricted. Define G(mu)
analogously for a standard Gaussian marginal. The cited fixed-law result
says B(mu) is comparable, by absolute constants, to

```math
\delta(\mu)=\inf_{a:T\to\mathbb R^d}
 \{\mathbb E_\mu\|a(t)\|_1+G((t-a(t))_\#\mu)\}.
```

The main-agent candidate is the exact minimax identity

```math
\inf_a\mathbb E\max_{t\in T}
 \{\|a(t)\|_1+\langle G,t-a(t)\rangle+h_t\}
=\sup_{\mu\in\mathcal P(T)}\{\delta(\mu)+\mathbb E_\mu h_t\}.
```

For fixed a, the left expectation is exactly the supremum over index
laws of its fixed-law coupling value. For fixed mu, the objective is
convex continuous in a; for fixed a it is concave upper-semicontinuous
in mu. Compactness of the finite probability simplex supplies the compact
side of Sion minimax. Thus there is no exchange of min and max over
nonconvex signings here.

Consequently this quantity is sandwiched between
Psi_(c epsilon)(T,h) and Psi_(C epsilon)(T,h), with the SAME offsets h.
This does not retain a sharp leading noise constant, but it avoids
discarding deficits or paying separate layer union bounds.
Root proof written; independent audit and novelty check pending.

The Bernoulli specialist is independently checking the related consequence
that a centered L-subGaussian vector with coordinate bound B is below
C(B+L) times independent signs in convex order. This uses the prescribed
law, not merely setwise expected maxima.

## 2. Localization-to-independent-rounding obstruction

The localization specialist has derived an actual-class candidate. For
ANY fixed full signing A, random switching X_ij=a_ij xi_i xi_j preserves
its scalar cap exactly. Given any observation Y, let b_ij=E[X_ij|Y] and
I=I(X;Y). Then

```math
\mathbb E\sum_{i<j}b_{ij}^2\le nI.
```

Root independently reconstructed this inequality: lift the conditional
law to globally symmetric xi. Each coordinate remains uniform. At a
fixed vertex i, conditional entropy subadditivity gives
sum_(j!=i) I(xi_j;xi_i|Y=y)<=D(P_(xi|y)||U_n). Binary Pinsker gives
correlation squared <=2 mutual information. Sum over i, divide by two,
and average. The signing signs cancel on squaring. This applies to
exact minimizing A, not just conference or Hadamard examples.

Independent edge rounding with means b then has expected cap at least

```math
\frac23\sqrt{\frac2\pi}\,n^{3/2}
-\sqrt{\frac2\pi}\,n\sqrt I-O(n^{4/3}).
```

Root independently derived the greedy lower bound using smooth absolute
values and coordinate Lindeberg replacement. This step keeps a rigorous
o(n^(3/2)) remainder without importing a sharper Wasserstein CLT.
Full proof and constants are being written by the localization specialist.

## 3. Additional resource tradeoff: candidate under audit

Let Z be ANY full-sign output law conditional on Y, with conditional
means r(Y). Set

```math
R=\mathbb E\|r(Y)-b(Y)\|_2^2,\qquad
D=\mathbb E D_{\rm KL}(P_{Z|Y}\|\bigotimes_e P_{Z_e|Y}).
```

Conditional product transport and the fact that Q is 2-Lipschitz in
edge Hamming distance suggest the fully quantified bound

```math
\mathbb E Q(Z)\ge
\frac23\sqrt{\frac2\pi}\,n^{3/2}
-\sqrt{\frac2\pi}\,(n\sqrt I+\sqrt{nR})
-\sqrt{n(n-1)D}-O(n^{4/3}).
```

This follows if the product-transport coupling bound
E d_H <=sqrt(binomial(n,2) D/2) is applied conditionally and averaged.
It separates information acquisition, coherent mean retuning, and
remaining coordinate dependence. The result would be about a specified
response-preserving architecture, NOT an impossibility of global rewrites:
returning a deterministic copy of A changes all posterior means and
escapes through the retuning term R.

The intended next task is to prove the transport step self-containedly,
test constants on exact examples, and see whether the positive comparison
theorems can construct a low-cost non-product replacement instead.
