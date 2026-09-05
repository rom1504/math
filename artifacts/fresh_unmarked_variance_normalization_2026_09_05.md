# A dimension-free mean-standard-deviation inequality for odd transports

Date: 2026-09-05. The trace inequality was independently checked by the
root and algebra tracks immediately after derivation. It is distinct from
the matrix-fidelity expression, whose inequality has the opposite direction.

## 1. Exact Schur-power trace inequality

Let `Q` be any real positive semidefinite correlation matrix of order `n`.
For every odd integer `r=2k+1 >= 1`,

`tr[Q sqrt(Q^(circ r))]`

`>= tr[(Q^(circ(k+1)))^(3/2)] >= n`.                 (1)

For `k=0`, the assertion is the scalar convexity bound on `tr(Q^(3/2))`.
For `k>=1`, define the unital completely positive Schur map

`Psi(X)=Q^(circ k) circ X`,

and put `P=Q^(circ(k+1))`. Then `Psi(Q)=P` and
`Psi(P)=Q^(circ(2k+1))`. The map is self-adjoint for the trace inner
product. Operator concavity of the square root gives

`sqrt(Psi(P)) >= Psi(sqrt(P))`.

Multiplying by the positive semidefinite `Q` and taking trace yields

`tr[Q sqrt(Psi(P))] >= tr[Psi(Q) sqrt(P)] = tr(P^(3/2))`.

Finally `tr P=n`; convexity of `x^(3/2)` on the nonnegative eigenvalues
proves the final bound in (1). In the cubic case, the important map is
`Psi=Schur(Q)`, with intermediate matrix `P=Q^(circ2)`.

The Jensen step can also be checked without importing a general Jensen
theorem. The block matrix `[[P,sqrt(P)],[sqrt(P),I]]` is positive
semidefinite. Applying the completely positive map entrywise to its
blocks and taking the Schur complement gives
`Psi(P) >= Psi(sqrt(P))^2`. Monotonicity of the positive square root then
gives the displayed inequality. Square-root monotonicity follows directly
from its resolvent integral and order reversal under positive inversion.

## 2. Odd mixtures

Let `w_r >= 0`, `sum_r w_r=1`, with odd positive integer indices, and set

`R=sum_(r odd) w_r Q^(circ r)`.

Concavity of the square root and (1) imply

`tr[Q sqrt(R)] >= sum_r w_r tr[Q sqrt(Q^(circ r))] >= n`.       (2)

Finite mixtures suffice for polynomial responses. Infinite mixtures follow
by finite approximation, since every Schur power has operator norm at
most `||Q||op` and the omitted weight tends to zero. The square root is
continuous on positive semidefinite matrices, including singular ones.

In particular, (2) applies to the covariance kernel of every normalized
odd Gaussian response with Hermite coefficients `h_r`, using `w_r=h_r^2`.

## 3. The mean local standard deviation

Let `B` be any real symmetric square root with `B^2=Q`; since `Q_ii=1`,
each column `b_i=B e_i` has Euclidean norm one. Define

`v_i=(B R B)_ii=||sqrt(R) b_i||_2^2`.

Columnwise Cauchy--Schwarz gives

`sqrt(v_i) >= b_i^T sqrt(R) b_i`.

Sum over `i`, then apply (2), to obtain the exact inequality

`boxed: n^-1 sum_i sqrt(v_i) >= 1`.                         (3)

This does not assert `v_i>=1` pointwise. That stronger statement is false,
including for full symmetric sign seeds that give scalable tensor families.
No flatness, coherence, or asymptotic assumption is needed in (1)--(3).

## 4. Bounded deterministic variance normalization

Fix `0<eta<1` and put

`d_i=1/sqrt(max(v_i,eta))`, `Ztilde_i=d_i Z_i`.

Then `0<d_i<=eta^-1/2`; if the old unmarked field has local Gaussian
variance `v_i`, the normalized field has variance at most one. Moreover,

`n^-1 sum_i d_i v_i >= n^-1 sum_i sqrt(v_i)-sqrt(eta)`

`                         >= 1-sqrt(eta)`.                    (4)

The scalar loss is zero when `v_i>=eta` and at most `sqrt(eta)` otherwise.
Thus small or zero variances require no division by zero and no exceptional
root deletion.

The weighted projection proof tolerates any uniformly bounded deterministic
row multipliers `d_i`: insert their diagonal matrix next to the local
weight. Every operator or Frobenius estimate acquires at most the fixed
factor `max|d_i|`; input derivatives never hit the deterministic multiplier.
The surviving coefficient sum becomes

`sum_(i,k) d_i B_ik (BR)_ik = sum_i d_i (BRB)_ii = sum_i d_i v_i`.

Accordingly, for a fixed cubic response and fixed local old-field functions
`M,F`, the variance-normalized version of the weighted theorem is

`n^-1 E [M(X) circ Ztilde]^T B F(X)`

`= EM(X) E[F(X)h3(G)] * n^-1 sum_i d_i v_i + o(1)`.            (5)

Here the full weighted theorem is invoked separately; the trace inequality
does not itself prove an energy identity. For a nonnegative `M` and a
positive Hermite coefficient, (4) makes its main coefficient bounded below
independently of the operator bound `L`.

## 5. Quantifier significance and remaining loss

For each fixed `L`, the normalized field has limiting local Gaussian tails
bounded by those of a standard Gaussian, uniformly in the root. Thus
cutoff radii no longer need to grow proportionally to `L`. A crude
operator-norm energy error still contains a factor `L`; a Gaussian tail
bound can absorb it with a radius of order `sqrt(log(L/t))` when the
local modification has Gaussian probability of order `t`.

This is substantially better than a radius proportional to `L`, but the
localized gain, its smoothing loss, and the spectral-deletion loss must
still be compared explicitly in the separate gain theorem. No universal
numerical improvement is asserted solely from (3)--(5).
