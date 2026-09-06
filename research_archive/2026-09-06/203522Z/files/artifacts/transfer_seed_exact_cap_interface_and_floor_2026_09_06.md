# An exact arbitrary-seed cap criterion, its Gaussian floor, and the residual

Date: 2026-09-06. Seed-transfer track. This document connects the actual
sign realization to a finite variational criterion BEFORE discussing
kernel regularization. The criterion has a rigorous universal floor
`sqrt(15)/8`; the regularized bound additionally retains a full quadratic
Gram optimization. Neither statement is a lower bound on the actual
ensemble's caps or on the original minimax sequence.

## 1. Put the arbitrary seed gate before the orthogonal children

Let `B` be ANY full symmetric sign seed of order `d`, and put `M=B/sqrt(d)`.
For terminal Hadamard order `q`, use the forward gate

```math
V=\operatorname{diag}(U_q^{(1)},\ldots,U_q^{(d)})
                          (M\otimes I_q)g,\qquad H=\sqrt{dq}\,V^T.
```

Every entry of `H` is a sign, exactly as for the reversed gate. Here the
input source is mixed by the finite seed BEFORE it reaches the independent
orthogonal child bases. This order permits an exact `d`-tuple method of
types. No orthogonality of `M` is needed.

For a symmetric finite law `nu`, define

```math
(B_M f)(\nu)=\sup_\pi\left\{
\frac1d\sum_i f(\mu_i)-\frac1dD(\pi\Vert\nu^{\otimes d})\right\},
```

where the input law `pi` has average absolute marginal `|nu|`, and `mu_i`
is the symmetrized law of `(MX)_i`. Global reversal lets one assume that
`pi` is centrally symmetric; then its actual signed marginals average to
`nu`. The type cost is exactly the displayed relative entropy divided
by `d`.

Apply depth-`r` recursive Hadamard children, use concatenation of the
tilted orbital norms, and apply their uniform terminal bounds. The
one-row pressure is at most

```math
p\log2+(B_M\Psi_{t,r})(\nu_p),\qquad
\Psi_{t,r}(\mu)=t\,m_2(\mu)+(\mathcal B^r\Phi_t)(\mu).
```

The alphabets after the seed gate are finite for each FIXED `d` and `p`.
The reachable output-law sets are compact finite simplices. The functions
`B^r Phi_t` are continuous there and decrease pointwise to a function at
most the continuous finite-source envelope `T_t`. Consequently, for
every positive `eta`, one common finite depth works on the entire reachable
set with `B^r Phi_t<=T_t+eta`. To check this uniform statement directly,
the compact sets `{mu:B^r Phi_t(mu)-T_t(mu)>=eta}` are nested, and their
intersection is empty; hence one of them is empty. No continuity of the
infinite-depth limit is asserted or needed.

The exact tilted weave theorem therefore implies

```math
\boxed{\quad
\limsup_{N\to\infty}\frac{M_N}{N^{3/2}}
\le C_B(p,t):=
\frac{p\log2+B_M(t m_2+T_t)(\nu_p)}{2t\sqrt p}.
\quad}                                                     (1)
```

As usual, prove every strict larger cap using a fixed depth, then let
the numerical slack go to zero. Relatively dense terminal Hadamard orders
and principal restriction give all target orders. This is a genuine
finite-seed sufficient inequality for the ORIGINAL cap, not an assumed
equality of the minimax value with an ensemble value.

For a hollow near-optimal seed `A` completed by a sign diagonal to `B`,
the desired transfer would require choices with

```math
C_B(p,t)\le Q(A)/d^{3/2}+o_d(1).
```

This has not been established.

## 2. A universal sqrt(15)/8 floor for this entire criterion

Take the admissible independent input law `pi=nu_p^d`. Every row of `M`
has squared Euclidean norm one, so every output marginal has second
moment one. The trivial label gives `T_t(mu)>=g_t(1)`. Therefore

```math
C_B(p,t)\ge\frac{p\log2+t+g_t(1)}{2t\sqrt p}.               (2)
```

This is exact for EVERY seed order and seed matrix; no central limit
theorem, spectral cutoff, or near-optimality assumption is used.

Write `2t=rho/(1-rho^2)`, `0<rho<1`. Since
`t+g_t(1)=t rho+(1/4)log(1-rho^2)`, differentiation in `t` shows that
the minimum of the right side of (2), at fixed `p`, occurs at

```math
1-\rho^2=2^{-4p},\qquad
\inf_{t>0}\frac{p\log2+t+g_t(1)}{2t\sqrt p}
=\frac{\sqrt{1-2^{-4p}}}{2\sqrt p}.
```

The square of the last expression decreases with `p`: for `a=4log2`,
the derivative of `(1-exp(-ap))/p` has numerator
`exp(-ap)(ap+1)-1<0`. Thus

```math
\boxed{\quad C_B(p,t)\ge\frac{\sqrt{15}}8
                           =0.484122918275927\ldots .\quad} (3)
```

The infimum of the lower bound is approached as `p` tends to one; its
temperature tends to `2sqrt(15)`. The same floor applies before replacing
finite-depth terminal potentials by `T_t`, since the independent child
policy and `Phi_t(mu)>=g_t(m_2(mu))` give the identical lower witness.

Hence this whole certificate architecture cannot transfer seeds whose
normalized caps approach a constant strictly below `sqrt(15)/8`.
The original liminf may or may not be below that number: the current
rigorous interval does not decide this. This floor concerns a proposed
UPPER certificate, not the actual matrix cap. Its equality with the
numerical value of an older GFOM ceiling does not identify the two
architectures.

## 3. A defect-insensitive regularization of the untilted kernel

This section records exactly what softening the precision factorization
does and why it does not itself solve the cap inequality (1).

Let `G=M^TM`, `epsilon>0`, `Lambda=diag(lambda_i)` with `0<lambda_i<=t`,
and define

```math
P_\epsilon=M^T\Lambda M+\epsilon tI,
\qquad t'=(1+\epsilon)t.
```

Adding `epsilon t||x||^2` permits completion of the square even for
singular `M`; the unused constant residual is nonnegative. The pivots
`delta_i` of `P_epsilon` obey
`delta_i<=bar(lambda)+epsilon t<=t'`. Exact determinant monotonicity gives

```math
\frac{\det P_\epsilon}{\det\Lambda}
=\det(MM^T+\epsilon t\Lambda^{-1})\ge\det(G+\epsilon I).
```

For `c_t(lambda)=(1/4)log[(lambda/t)(2-lambda/t)]`, the kernel prefactor
is consequently at most

```math
\sum_i c_t(\lambda_i)-\sum_i c_{t'}(\delta_i)
\le-\frac14\log\det(G+\epsilon I)
                              +\frac d2\log(1+\epsilon).   (4)
```

For the factor involving `2-lambda/t`, use Jensen on the child values
and `delta_i/t'<=(bar(lambda)/t+epsilon)/(1+epsilon)` on the parents.
Their log ratio is at most `d log(1+epsilon)`. The remaining determinant
and temperature ratio give the other `d log(1+epsilon)` before division
by four. This proves (4) uniformly over all precision choices and
elimination orderings.

The finite-source Gaussian-mixture dual makes (4) a source inequality:

```math
B_M T_t(\nu)\le T_t(\nu)+\epsilon t m_2(\nu)+D_\epsilon(M),
\qquad
D_\epsilon(M)=-\frac1{4d}\log\det(G+\epsilon I)
                                      +\frac12\log(1+\epsilon). (5)
```

Here is the dual step, to specify the claim completely. Finite-source
kernel-mixture maximization has an optimal kernel vector `k_i>0`.
Its supporting hyperplane supplies `a_i=1/k_i` with
`sup_(y,lambda) sum_i nu_i a_i kernel_(t')(x_i;y,lambda)<=1` and
`-sum_i nu_i log a_i=T_(t')(nu)`. Multiply these dual weights by
`exp(-epsilon t x_i^2)`, factor each child Gaussian kernel with
`P_epsilon`, and integrate parent coordinates in the triangular order.
Each integral is at most one; (4) pays the uniform prefactor. The entropy
variational inequality and the average-marginal constraint give
`B_M T_t<=T_(t')(nu)+epsilon t m_2(nu)+D_epsilon`.
Finally `T_(t')<=T_t`, proving (5). Compactness of the finite kernel-vector
body and its positive feasible point justify the supporting hyperplane;
no minimax exchange involving the original cap is used.

Unlike a raw determinant, (5) is insensitive to `o(d)` very small
eigenvalues at any fixed `epsilon`. It therefore removes the twin-pivot
pathology at the ANALYTIC source-envelope level when the remaining
regularized bulk defect is small.

## 4. The residual for the ACTUAL cap is exactly a Gram quadratic problem

The tilted function in (1) has an additional output-energy term. Define

```math
\Gamma_M(\nu)=\frac1d\sup_\pi
                  \mathbb E X^T(G-I)X,
```

with the same average absolute-marginal constraint. Combining this exact
energy identity with (5) yields only

```math
B_M(t m_2+T_t)(\nu)
\le t m_2(\nu)+T_t(\nu)
       +t\Gamma_M(\nu)+\epsilon t m_2(\nu)+D_\epsilon(M).   (6)
```

For the symmetric Boolean source this residual is precisely

```math
\Gamma_M(\nu_1)=\frac1d\max_{x\in\{\pm1\}^d}x^T(G-I)x.
```

For the ternary source, put
`F_M(r)=max_{z in {0,+1,-1}^d, ||z||_0=r}||Mz||^2`.
Then it is the one-moment concave-hull problem

```math
\Gamma_M(\nu_p)
=\frac{\operatorname{cav}(F_M)(pd)}{dp}-1.                  (7)
```

Indeed global sign averaging makes every chosen configuration symmetric,
and the only remaining average-marginal constraint is that its number
of nonzero coordinates average to `pd`. Their common amplitude is
`1/sqrt(p)`. This proves (7), including nonintegral `pd`.

Independent inputs show `Gamma_M>=0`. Also `tr G=d`, so determinant
concavity gives `D_epsilon(M)>0`. Thus (6) cannot improve the existing
universal scalar-`T` bound; it can only add Gram costs. The planted-clique
construction shows that this Gram maximum can even diverge on some
asymptotically minimizing seeds, despite a vanishing normalized cap
perturbation.

The exact criterion (1) retains more seed information than (6), since
its output source laws and entropy are optimized jointly. But (3) is
already an unavoidable floor for that exact criterion. No additional
state or regularized determinant alone changes either fact.
