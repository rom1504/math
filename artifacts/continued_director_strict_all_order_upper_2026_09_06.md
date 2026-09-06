# A strict all-order asymptotic upper bound for quadratic signings

Date: 2026-09-06. Director reconstruction, independently audited. This is
an unconditional improvement for the ORIGINAL minimax sequence, not merely
a special-class lower bound or a finite-order experiment. It does not prove
convergence or nonconvergence.

## Theorem

Let

```math
a=\frac{91470529542342299}{20460000000000000000},\qquad
c_* = \frac12-\frac{a}{8\sqrt{31/32}}.
```

Then

```math
\limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
\le c_* < 0.499432220485404 < \frac12.
```

The approximate value of c_* is .49943222048540312077. The displayed
decimal upper bound is outward checked in rational arithmetic: replacing
sqrt(31/32) by the existing exact square-root upper endpoint gives
`80459630021641337701/161102201102367360000`, strictly below
`499432220485404/10^15`.

In particular a limiting value of 1/2 is ruled out. Existence of any limit
remains open. Combined with the preserved universal lower theorem, the
current rigorous interval is

```math
0.4333221116640807\le\liminf_n M_n/n^{3/2}
\le\limsup_n M_n/n^{3/2}\le c_*.
```

## 1. The new factorization, not an assumed alignment of temperatures

Let Phi_t(nu)=-F_t(nu)/2 be the entropic Gaussian self-transport potential,
g_t(v)=Phi_t(N(0,v)), and B the exact signed-pair Bellman operator of the
recursive Hadamard ensemble. Define

```math
T_t(\nu)=\sup_L\left\{
 \mathbb E_L g_t(\operatorname{Var}(X\mid L))-I(X;L)\right\}.
```

For finite symmetric sources the [conditional-variance supersolution
theorem](continued_convergence_conditional_variance_supersolution_2026_09_06.md)
proves `B T_t<=T_t` and `G_t<=T_t<=Phi_t`, where G_t(nu)=g_t(m_2(nu)).
This is NOT the still-open inequality for the smaller envelope
g_t(E Var(X|L))-I(X;L).

The essential new operation is sequential labeling: for child labels M,N
and rotated inputs A,B, label A by (M,N,B), but B only by (M,N). The sum
of information costs is at most the two child costs plus I(A;B).
Conditional linear regression produces a Schur-complement variance.
Concavity of g_t as a function of log variance and ordinary convexity of
g_t compare this Schur pair through the covariance eigenvalues to the two
rotated diagonal variances. This preserves both child rewards without
assuming equal temperatures, Gaussian posteriors, or independent residuals.
Every source in a finite-depth tree is finite, so the sequential labels
are finite as well.

## 2. Removing the terminal-potential gap

The [unbounded self-transport proof](continued_convergence_unbounded_sinkhorn_rigidity_2026_09_06.md)
establishes strict source concavity on every finite-second-moment law.
It proves that zero one-step Bellman drift occurs exactly for centered
Gaussian laws. The [Gaussian-boundary replacement theorem](continued_convergence_terminal_gap_reduction_2026_09_06.md)
then gives

```math
\lim_{r\to\infty}\mathcal B^r\Phi_t(\nu)
=\lim_{r\to\infty}\mathcal B^rG_t(\nu)\le T_t(\nu).
```

The first limit decreases and the second increases. The proof controls
moment escape using the logarithmic Gaussian self-cost and a bounded-time
stopped-tree estimate, not a uniform central limit theorem for controlled
policies. The supersolution is needed only on finite reachable states.

## 3. An exact negative root certificate

For p=31/32, t=4 and
`nu_p=(1-p)delta_0+(p/2)delta_(+1/sqrt(p))+(p/2)delta_(-1/sqrt(p))`,
the [exact rational certificate](continued_feedback_conditional_variance_exact_certificate_2026_09_06.md)
proves

```math
p\log2+T_4(\nu_p)+4(1-\sqrt p)\le-a<0.
```

The posterior optimization is reduced exactly to a one-moment concave
hull over two scalar coordinates. Its verifier checks every point of a
2501-by-2501 grid with exact integer/rational upper values; an explicit
mean-preserving rounding modulus pays for all unsampled posterior laws.
No numerical optimizer or machine-learning prediction enters the bound.
The director and independent auditor both reran the complete verifier and
obtained the identical fraction.

Therefore, for every a' with 0<a'<a, there exists a FINITE depth r such that

`p log2+(B^r Phi_4)(nu_p)+4(1-sqrt(p))<=-a'`.

No effective value of r is asserted or needed for the asymptotic existence
theorem. This is the only passage from an infinite recursion to a finite
construction.

## 4. From the exact ensemble to all orders

The [uniform orbital/type theorem](continued_convergence_recursive_orbit_bound_2026_09_06.md)
applies at that fixed r to exact randomized Hadamard bases. Its terminal
loss is exp(O(sqrt(m))) and its finite-alphabet type errors are o(m).
The [restricted-weave identity and counting inequality](continued_convergence_restricted_weave_2026_09_06.md)
use fresh independent bases in different fibres: the averaged bound is
(E Z)^m, not E[Z^m]. The exact squared defect is
`D_sigma=2(m^2 k-sigma x^T K x)`, and the exponential weight is
`exp(-t D_sigma/(2k))`. The PSD-kernel graph contraction and spin union
bound convert the negative row exponent into existence of an actual
symmetric full sign matrix. Removing its diagonal costs at most N/2.

For every fixed 0<eta<a'/4 this gives normalized hollow cap at most
`1/2-eta/(2sqrt(p))+o(1)` along the constructed orders. The
[all-order realization proof](continued_director_recursive_weave_all_order_implication_2026_09_06.md)
uses arbitrary terminal Hadamards, a relatively dense supply of their
orders from primes 3 modulo4, and principal restriction. It does not
assume that near-minimizers are conference matrices or that conference
extremal responses describe this ensemble.

Thus every such eta bounds the ALL-ORDER limsup. First take the order
limit at fixed r,a',eta; then let eta approach a'/4 and a' approach a.
This gives exactly c_*. There is no growing-depth type estimate, exchange
of uncontrolled limits, or accumulated insertion error.

## Verification and remaining question

The director reconstructed the factorization, numerical certificate, full
weave/Fock/type chain, normalization and all-order extension. A separate
reviewer freshly reconstructed the same chain and replayed the exact
certificate; another researcher checked the new algebra and all imported
interval primitives. The detailed audit is
`continued_audit_exact_conditional_envelope_to_all_order_cap_2026_09_06.md`;
the supporting algebra audit is
`continued_audit_temperature_alignment_boundary_2026_09_06.md`.

The result provides a uniform upper construction, not a matching lower
bound and not a transfer theorem for arbitrary optimizing seeds. The
remaining original obligation is still
`limsup_n M_n/n^(3/2) <= liminf_n M_n/n^(3/2)`
or a strict separation of two infinite subsequences. No equality with
c_* or with the optimum of the weave ensemble is conjectured as a proof
step. The new theorem invalidates convergence to1/2, not convergence itself.
