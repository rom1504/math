# Actual designer Gibbs laws are macroscopically nonprojective

Date: 2026-09-06. This is an original-signing consequence of the proved
bounded-cap principal-restriction theorem. It obstructs a particular
same-temperature pressure interpolation, not convergence of the optimized
pressure or of the original extremal sequence.

Write `e_n=n(n-1)/2`, `Q(A)=max_x |x^T A x|/2`, and
`m_n=M_n/n^(3/2)`. On ALL hollow symmetric signings of order `n`, define

```math
 Z_n(\beta)=\sum_A e^{-\beta\sqrt n Q(A)},\qquad
 \Phi_n(\beta)=n^{-2}\log Z_n(\beta),\qquad
 \mu_n^\beta(A)=Z_n(\beta)^{-1}e^{-\beta\sqrt n Q(A)}.
```

This partition function is over the designer's EDGE signs. It is not
the ordinary positive-temperature spin partition function. Entropies
and relative entropies below use natural logarithms.

## 1. Statement with flexible constants

Assume the established all-order bound `limsup m_n<=U<2/pi`.
Fix `beta>0` and `c<2/pi` such that

```math
             g:=\beta(c-U)-\frac{\log2}{2}>0.             (1)
```

Then there exists `rho>0` such that, whenever `k=k(n)->infinity` and
`k<=rho n`, the principal `k`-coordinate marginal `nu_(n,k)` of
`mu_n^beta` satisfies

```math
 \nu_{n,k}\{Q\ge c k^{3/2}\}=1-O(1/k)-e^{-\Omega(n^2)},  (2)
```

whereas the standalone law at that same temperature satisfies

```math
 \mu_k^\beta\{Q\ge c k^{3/2}\}
       \le \exp\{-(g-o(1))k^2\}.                         (3)
```

Consequently

```math
 \|\nu_{n,k}-\mu_k^\beta\|_{\rm TV}\longrightarrow1,
 \qquad
 \liminf\frac{D(\nu_{n,k}\Vert\mu_k^\beta)}{k^2}\ge g.  (4)
```

The marginal can be taken on the first `k` coordinates, with no random
selector retained in the state. This uses exact vertex-permutation
invariance of the designer Gibbs law.

For a concrete choice requiring no new numerical upper certificate,
take `U=1/2`, `beta=4`, and `c=5/8<2/pi`. Then

```math
                  g=\frac{1-\log2}{2}>0.1534.             (5)
```

In particular (4) holds along every sufficiently small FIXED positive
retention `k=floor(pn)`, not just a vanishing-retention limit.

## 2. Proof, including all probabilistic quantifiers

For any `C0>U+(log2)/(2 beta)`, the single minimizing signing in the
partition function gives

```math
 \mu_n^\beta\{Q>C0 n^{3/2}\}
 \le \exp\{e_n\log2-\beta(C0-m_n)n^2\}
 =e^{-\Omega(n^2)}.                                     (6)
```

Choose such a FIXED `C0`. The actual bounded-cap theorem proved in
`transfer_adversary_fixed_retention_random_loss_2026_09_06.md`, with its
spectral-core and fixed-degree Gaussian-witness proof in
`transfer_seed_sparse_bounded_cap_parent_theorem_2026_09_06.md`, says:
for this `C0,c` there exist `rho,K,k0` such that EVERY parent obeying
`Q(A)<=C0 n^(3/2)` and EVERY `k0<=k<=rho n` has

```math
 \Pr_T\{Q(A_T)<c k^{3/2}\}\le K/k.                       (7)
```

The constants are independent of the particular parent. Thus (7) can
be integrated over the Gibbs law; no replacement of a random parent
by a presumed typical deterministic optimizer is involved. Add (6).
Since `mu_n^beta` is permutation invariant, the resulting averaged
restriction law equals its first-`k` principal marginal. This proves
(2). The same one-atom partition lower bound at order `k` gives

```math
 \mu_k^\beta\{Q\ge c k^{3/2}\}
 \le\exp\{e_k\log2-\beta(c-m_k)k^2\},                   (8)
```

which is (3). If this event is empty then (2) already rules that out
for all sufficiently large relevant orders.

Let `E_k={Q>=c k^(3/2)}`, `a=nu_(n,k)(E_k)`, and
`b=mu_k^beta(E_k)`. The event itself proves total-variation separation.
Binary data processing gives

```math
 D(\nu_{n,k}\Vert\mu_k^\beta)
 \ge a\log(1/b)-h(a),                                  (9)
```

because the omitted term `-(1-a)log(1-b)` is nonnegative.
Equations (2)--(3) prove (4).

For (5), one may use `C0=1` in (6), since
`4(1-1/2)>(log2)/2`. The elementary inequality `pi<16/5` verifies
`5/8<2/pi`.

## 3. Exact variational interpretation

For a law `nu` on order-`k` signings put

```math
 {cal V}_{k,\beta}(\nu)
       =H(\nu)/k^2-\beta\,\mathbb E_\nu Q/k^{3/2}.
```

There is the exact identity

```math
 \Phi_k(\beta)-{\cal V}_{k,\beta}(\nu)
       =D(\nu\Vert\mu_k^\beta)/k^2.                     (10)
```

Hence the ACTUAL Gibbs marginal loses at least the fixed amount `g`
from the optimum child variational value. Also, directly using
`H(nu)<=e_k log2`, (2), and `Phi_n>=-beta m_n`,

```math
 \liminf\{\Phi_n(\beta)
              -{\cal V}_{k,\beta}(\nu_{n,k})\}\ge g.    (11)
```

This diagnoses the obstruction to a naive same-temperature principal
restriction/Gibbs--Shearer proof: the actual induced child law is not
even an asymptotically optimal competitor on its own scale. No entropy
estimate that merely reveals its exact entropy can remove (10).

It does NOT refute a comparison of the optimized numbers `Phi_n,Phi_k`
using a different child law, a designed rare selector, a change of
temperature, or a genuinely non-principal construction. In particular,
nonprojectivity of finite-volume Gibbs laws is not nonconvergence of
their normalized free energies.

## 4. Dependency audit and finite replay

For (7), the main dependency was reread in full, including its uniform
connected-graph covariance bound, fixed polynomial degree before order
limits, simultaneous two-sided spectral majorant, and conditional
hypergeometric core intersection. Only its polynomial `K/k` bound is
needed here; the later convex-distance exponential selector theorem
is not imported. The proof above adds no Gaussian or Hadamard assumption
on the designer Gibbs samples.

`computations/decisive_audit_designer_gibbs_checks_2026_09_06.py`
exhaustively forms the finite designer laws through order six and checks
the exact marginal variational identity and finite one-atom tail bounds.
Those checks do not replace the asymptotic restriction theorem.
