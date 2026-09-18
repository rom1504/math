# Initial mechanism frozen before archive inspection

Time: 2026-09-17, first localization reading. Status: independently derived
mapping; not a convergence claim. No previous synthesis or ledger was read.

## Primary theorem selected

El Alaoui--Montanari, *An Information-Theoretic View of Stochastic Localization*,
arXiv:2109.00709v2, Theorems 1--2. The Gaussian-channel construction is selected
as the main paper, with Eldan arXiv:1811.11530 as the original result.

For finite-support X, covariance Sigma, precision L positive definite, observe
Y=sqrt(tau) X + L^{-1/2} Z, tau uniform on [1,2]. The conditional laws mu_theta
are supported on the exact original support. If C_theta is posterior covariance,
then

```math
\mathbb E C_\theta\preceq L^{-1},\qquad
I(X;\theta)\leq\tfrac12\log\det(I+2L^{1/2}\Sigma L^{1/2}),\qquad
\mathbb E[C_\theta L C_\theta]\preceq\Sigma.
```

The decisive proof identity is d E C_t/dt = -E[C_t L C_t], obtained either by
Gaussian integration by parts or the posterior mean martingale
dm_t=C_t L^{1/2} dW_t and total covariance. This finite-support formulation
avoids singular inverse and integrability ambiguities.

## Exact spin mapping and quantitative barrier

Write H_A(x)=x^T A x/2 for full hollow symmetric sign A. For any law mu on
spin signs and posterior mean m_theta,

```math
\mathbb E_\mu H_A(X)
=\mathbb E_\theta H_A(m_\theta)
 +\tfrac12\mathbb E_\theta\operatorname{Tr}(A C_\theta).
```

Independent sign rounding of m_theta preserves H_A(m_theta) in conditional
expectation but need not preserve the original support or near-extremality.
The localization theorem and weighted Hilbert--Schmidt Cauchy--Schwarz give

```math
\mathbb E|\operatorname{Tr}(A C_\theta)|
\leq \sqrt{\operatorname{Tr}(A L^{-1}A)\operatorname{Tr}\Sigma}.
```

Thus L=r I yields energy error at most n sqrt(n-1)/(2 sqrt(r)). The entropy
certificate is at most n log(1+2r)/2, using Tr Sigma<=n and concavity. Taking
r small to certify o(n) information cannot certify o(n^{3/2}) energy error.
This is a limitation of the generic bound, not yet an impossibility theorem.
Anisotropy replaces n(n-1)/r by Tr(A L^{-1}A); a useful extension must exploit
the interaction between Sigma, A, and signed trace cancellations, not merely
the ambient dimension.

## Exact gradient complexity check

For f_beta(x)=beta H_A(x)/sqrt(n), the Gaussian width of the continuous
gradient image of the cube is exactly

```math
\mathbb E\sup_{x\in\{-1,1\}^n}\langle g,\nabla f_\beta(x)\rangle
=\frac{\beta}{\sqrt n}\mathbb E\|Ag\|_1
=\beta\sqrt{\frac2\pi}\,n\sqrt{1-\frac1n}.
```

Every row of A has squared Euclidean norm n-1. Hence at fixed positive beta
the width is linear in n for every full sign A, not just a worst-case example.
No low-gradient-complexity premise may be silently imported at the spin-glass
normalization. The scale beta/n instead has width O(sqrt(n)), but energy itself
then lives at a different scale.

## Alternative edge-sign mapping

One can take X to be the N=binom(n,2) coefficient signs instead, and the scalar
objective becomes the supremum of a linear Bernoulli process indexed by spin
signs. Posterior laws retain full sign support; posterior means do not.
Small posterior coordinate covariance alone is insufficient to transfer this
supremum through independent rounding. The needed quantity is a process
metric/complexity bound for the posterior residual, with loss o(n^{3/2}).

## First extension candidates

1. Adaptive anisotropic observation aimed only at the signed scalar
   Tr(A C_t), charged by dI/dt=Tr(L_t C_t)/2.
2. Near-extremal spin prior plus entropy-to-mean estimates: decide whether
   scalar energy-preserving productization requires Omega(n) information.
3. Exact edge-sign coupling: use the decomposition only to condition a good
   coefficient law, and retain a correlated posterior sample rather than
   replacing it by independent coordinates.

These are hypotheses to test. None presently yields convergence.
