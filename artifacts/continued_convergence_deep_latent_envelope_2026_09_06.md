# A universal latent-variable lower envelope for deep recursive pairing

Date: 2026-09-06. This concerns the exact Bellman operator from
`continued_convergence_recursive_orbit_bound_2026_09_06.md`, not an upper
bound on the original minimax problem. It is a potential falsifier of the
finite recursive-weave certificate and of oversimplified deep-limit claims.

Write `f_r(nu)=B^r Phi_t(nu)` and

`g_t(v)=Phi_G(tv)=-tv(1-rho)+(1/4)log(1-rho^2)`,
`2tv=rho/(1-rho^2)`.

The convention is `g_t(0)=0`. Here `Phi_G` is the exact centered Gaussian
self-transport potential, equivalently the folded half-normal potential.

## 1. Gaussian lower bound at every finite depth

For a finite symmetric law `nu` of variance `v`, choosing iid pairs at
every level gives the normalized convolution law of `2^s` iid copies.
The central limit theorem gives convergence to `N(0,v)` in `W_2`.
The entropic self-transport cost is continuous in `W_2`: pushing a coupling
through a transport kernel contracts relative entropy, and its quadratic
cost changes by at most a constant times
`W_2(mu,nu)(sqrt(m_2(mu))+sqrt(m_2(nu)))`, with the corresponding reverse
comparison. Thus the selected iid terminal potential tends to `g_t(v)`.

Since `B Phi_t<=Phi_t` and `B` is monotone, the iterates decrease.
Appending arbitrarily many iid levels therefore proves

`f_r(nu) >= g_t(v)` for every finite `r`.

The earlier full-condensation policy also gives `f_r(nu)>=-H(nu)`.
Consequently `f_infinity(nu)=lim_r f_r(nu)` exists and is finite.

## 2. Conditional-copy construction and telescoping information cost

Let `X` have finite symmetric law `nu`, and let `L` be any finite latent
variable jointly distributed with `X`. We may assume reflection equivariance
without changing `I(X;L)` or `E Var(X|L)`: independently reflect `(X,L)` by
a fair sign and include that sign in the latent label. Symmetry of `X`
makes the additional sign independent of the reflected source.

Set `V=E Var(X|L)`. At level `d`, define `X_d` by conditionally summing
`2^d` iid copies of `X` given the **same** latent label and dividing by
`2^(d/2)`. Thus

`E Var(X_d|L)=V`,
`E[X_d|L]=2^(d/2) E[X|L]`.

In the one active plus branch, pair conditionally independent copies
`A_d,B_d` of `X_d` given `L`. The plus output has exactly the law of
`X_(d+1)`. The minus output
`Z_d=(A_d-B_d)/sqrt(2)` has unconditional mean zero and variance exactly
`V`, even when the different latent classes have different conditional
variances. Send this minus branch through any remaining recursion, using
the universal Gaussian lower bound from Section 1.

Let `J_d=I(X_d;L)` and `c_d=I(A_d;B_d)`. Conditional independence gives

`I(A_d,B_d;L)=2J_d-c_d`.

Data processing through the plus output yields

`c_d <= 2J_d-J_(d+1)`.

The actual cost at level `d` is `2^(-d-1)c_d`. Hence the total cost along
the active branch telescopes:

`sum_(d=0)^(r-1) 2^(-d-1)c_d <= J_0-2^(-r)J_r`.

This argument does not replace the conditional variance by a common
classwise variance. Only its unconditional average enters each minus
branch's Gaussian lower bound, and that average is exactly preserved.

The active terminal entropy is `O(r)`. Indeed, if the source alphabet has
`k` points, a sum of `2^r` source values has at most `(2^r+1)^k` possible
count vectors, independently of any additive relations among those points.
Thus `2^(-r) Phi_t(law X_r)>=-2^(-r-1)H(X_r)->0`.

Combining the branch contributions and taking `r` to infinity proves

`f_infinity(nu) >= g_t(V)-I(X;L)`.                            (1)

Since finite iterates lie above their limit, (1) is also a lower bound for
every finite-depth Bellman value. Optimizing over latent variables gives
the rate-distortion lower envelope

`f_infinity(nu) >= sup_L { g_t(E Var(X|L))-I(X;L) }`.          (2)

Equation (2) is an original conditional-copy derivation; no imported
rate-distortion theorem is needed to define or justify its right side.

## 3. Why Gaussian spreading plus total-entropy condensation is incomplete

The tempting candidate

`max{g_t(Var X),-H(X)}`

does not describe the deep limit for general finite symmetric input laws.
Take `X=S+epsilon T`, with independent uniform signs `S,T`,
`epsilon=1/100`, and `t=4`. Then `H(X)=log4`. Taking `L=S` in (1) gives

`f_infinity(law X) >= g_4(1/10000)-log2`.

This is about `-0.693547`, whereas the Gaussian candidate is about
`-0.777695`; the `-log4` alternative is smaller still. The coarse sign can
be concentrated while its small residual noise is spread over the discarded
branches. This is an actual lower-envelope counterexample, not merely a
failure to prove a proposed supersolution.

This four-atom example is not the original ternary input `nu_p`. Its role
is to expose the extra obligation in any universal deep-limit theorem:
coarse latent structure cannot be replaced by total input entropy alone.

The counterexample is already visible at a finite depth. With `n=2^r`,
the active law is `sqrt(n)S+(epsilon/sqrt(n))sum_(i=1)^n T_i`. Its entropy
is at most `log2+log(n+1)`, while every active pairing has information cost
exactly `log2` because `S` remains recoverable from the sign. Consequently

`f_r(law X) >= (1-2^(-r))[g_4(1/10000)-log2]`
`             -2^(-r-1)[log2+log(2^r+1)]`.

At `r=4` this is above `-0.760399099`, exceeding the Gaussian candidate
by more than `0.0172954654`. At infinite depth, the gap is greater than
`0.0841475430`. These signs are independently reproducible using the
60-digit interval checker

```
.venv/bin/python computations/continued_convergence_deep_latent_envelope_2026_09_06.py --output computations/continued_convergence_deep_latent_envelope_2026_09_06.json
```

No four-atom numerical optimizer is used in either asserted gap.

## 4. Explicit ternary latent channels

For the original input `nu_p`, the exact latent sign channel assigns the
sign at a nonzero source and a fair latent sign at zero. It has

`I(X;L)=p log2`, `E Var(X|L)=1-p`.

Thus (2) recovers the aligned-column residual lower bound, now directly
inside the recursive Bellman framework. A binary symmetric degradation
with correlation `c` gives

`I(X;L)=p[log2-h((1-c)/2)]`, `V=1-p c^2`.

More general symmetric three-label channels can be parameterized by
`w,c in [0,1]`, `b in [0,1/2]`. At zero the probabilities of labels
`(+,-,0)` are `(b,b,1-2b)`. At the positive source they are
`(w(1+c)/2,w(1-c)/2,1-w)`, and at the negative source the first two swap.
Put `q=[pw+2(1-p)b]/2`. Then

`V=1-p w^2 c^2/(2q)`,

`I=H(q,q,1-2q)-(1-p)H(b,b,1-2b)`
`  -p H(w(1+c)/2,w(1-c)/2,1-w)`.

These are concrete admissible lower tests. Optimizing them is not a proof
that all latent channels have been covered, and no negative optimizer
value is an upper certificate for the original Bellman problem.

## 5. The envelope is a subsolution; each fixed-temperature branch is a supersolution

Denote the right side of (2) by `E_t(nu)`. It satisfies the exact inequality

`B E_t >= E_t`.                                               (3)

To prove this, take any latent variable `L` for the input and the same
conditionally independent pair as in Section 2. In the plus child use `L`
again; in the minus child use the unconditioned Gaussian lower bound.
Both relevant residual variances are `V`. The resulting Bellman value is
at least

`g_t(V)-[I(U;L)+I(A;B)]/2 >= g_t(V)-I(X;L)`.

Taking the supremum over `L` proves (3). Thus this envelope is not an upper
bound merely because it is a natural limit candidate. A strict example
`B E_t(nu)>E_t(nu)` would refute the proposed identification with the deep
Bellman value while leaving (2) intact. The audit agent independently
reconstructed (3).

There is also an exact opposite-direction statement at a fixed auxiliary
temperature. Define the quadratic rate-distortion Lagrangian directly by

`J_lambda(nu)=inf_L [I(X;L)+lambda E Var(X|L)]`.

The Gaussian potential has the supporting-line representation

`g_t(v)=sup_(0<lambda<=t) [c_t(lambda)-lambda v]`,
`c_t(lambda)=(1/4)log[lambda(2t-lambda)/t^2]`.

Therefore

`E_t(nu)=sup_(0<lambda<=t) [c_t(lambda)-J_lambda(nu)]`.         (4)

For **each fixed** `lambda`, the branch
`e_lambda(nu)=c_t(lambda)-J_lambda(nu)` is a Bellman supersolution:

`B e_lambda <= e_lambda`.                                    (5)

Indeed, for a joint source `(A,B)` the Lagrangian satisfies

`J_lambda(A)+J_lambda(B)-I(A;B)`
` <= J_lambda(A,B) <= J_lambda(A)+J_lambda(B)`.

The lower bound follows from
`I(A,B;L)>=I(A;L)+I(B;L)-I(A;B)` and additivity of conditional squared
error. For the upper bound, use independent conditional reproduction
channels for the two marginals; their joint information is no greater
than the sum of the two marginal informations, and conditioning on both
labels can only improve squared prediction error. The vector Lagrangian
is invariant under orthogonal changes of coordinates. Hence

`J_lambda(U)+J_lambda(V)+I(A;B) >= 2J_lambda(nu)`,

which is (5) after multiplying by `-1/2` and adding the constant.
Here one first uses the always-safe simultaneous-reversal and input-swap
averages, so both signed input marginals are exactly `nu`; the raw
average-absolute-marginal constraint alone would not justify `2J_lambda(nu)`.

Equations (3)--(5) isolate a precise issue: the children in `B E_t` may
choose **different** auxiliary temperatures. A supremum of supersolutions
need not be a supersolution. Even proving `B E_t<=E_t` would not alone
prove that the decreasing terminal iterates converge to it; terminal-gap
decay or an appropriate fixed-point uniqueness theorem would still be needed.

## 6. Bounded general-coupling tests beyond ternary sources

`computations/continued_convergence_latent_supersolution_test_2026_09_06.py`
tests (only as a necessary condition) whether `B E_t<=E_t`. Its child
values are explicit feasible latent lower bounds. Its source envelope
has a separate upper approximation, as follows.

For a finite source `x_i` with probabilities `p_i`, and a finite grid of
reproduction points `y_j`, write

`J_lambda,grid = min_q -sum_i p_i log sum_j q_j exp[-lambda(x_i-y_j)^2]`.

Any positive numbers `a_i` give a finite-grid dual lower bound

`sum_i p_i log a_i - log max_j sum_i p_i a_i exp[-lambda(x_i-y_j)^2]`.

If the reproduction grid has mesh `h`, an optimal continuous reproduction
can be randomly rounded to its two neighboring grid points, with conditional
mean unchanged. This contracts mutual information and increases mean-square
error by at most `h^2/4`. Therefore

`J_lambda,grid-lambda h^2/4 <= J_lambda <= J_lambda,grid`.

The reproduction range may be restricted to the convex hull of the source,
by projection or by replacing a reproduction with its posterior mean.
Since `J_lambda` and `c_t(lambda)` are increasing in `lambda`, a source
envelope upper bound over `[a,b]` is
`c_t(b)-J_a`; intervals beginning at zero use `c_t(b)`.
This gives a controlled one-dimensional envelope search. Current reported
numbers use floating arithmetic, not interval-certified global inequalities.

Three bounded tests used four-atom symmetric sources: magnitudes
`(0.99,1.01)` with equal weights at `t=4`, magnitudes `(0.1,2)` with
weights `(0.9,0.1)` at `t=4`, and magnitudes `(0.4,1.2)` with equal
weights at `t=8`. No positive source-upper/child-lower gap was found.
Explicit boundary policies are included, since a smooth optimizer alone
missed the obvious conditional-copy policy in the close-cluster example.
These non-falsifications do not establish (5) for the supremum envelope.
