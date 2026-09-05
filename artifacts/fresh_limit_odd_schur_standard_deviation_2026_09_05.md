# A lower bound for the mean standard deviation of odd Schur transports

Date: 2026-09-05. The cubic case was found by the literature agent and
independently checked by the algebra agent and director. The following
generalization includes every odd Schur power and their convex mixtures.

## Theorem

Let `Q` be any real positive semidefinite correlation matrix of order `n`.
For nonnegative weights `w_r` indexed by positive odd integers with
`sum_r w_r=1`, set

`R=sum_{r odd} w_r Q^{circ r}`.

Then

`Tr[Q sqrt(R)] >= n`.                                           (1)

If `B` is any real symmetric square root of `Q` and
`v_i=(B R B)_{ii}`, then

`sum_i sqrt(v_i) >= n`.                                         (2)

Neither statement assumes sign-flat entries, bounded operator norm, or
any coherence condition.

## Proof

First let `r=2k+1` be a fixed positive odd integer. Schur multiplication
by `Q^{circ k}` is a unital completely positive, trace-self-adjoint map
`Psi`; for `k=0` use the all-ones Schur multiplier. Set
`P=Q^{circ(k+1)}`. The Schur product theorem shows that `P` is itself a
correlation matrix. The important identities are

`Psi(Q)=P`, `Psi(P)=Q^{circ(2k+1)}`.

Operator concavity of the square root gives

`sqrt(Psi(P)) >= Psi(sqrt(P))`.

One can verify this instance directly: the positive block matrix
`[[P,sqrt(P)],[sqrt(P),I]]` stays positive under the blockwise completely
positive map. Its Schur complement gives
`Psi(P)>=Psi(sqrt(P))^2`; operator monotonicity of square root and
positivity of `Psi(sqrt(P))` give the displayed inequality.

Take trace against the positive matrix `Q` and use trace self-adjointness:

`Tr Q sqrt(Q^{circ r})
 >= Tr Q Psi(sqrt(P))
 = Tr Psi(Q) sqrt(P)
 = Tr P^{3/2}
 >= n`.

The final inequality is scalar Jensen applied to the nonnegative
eigenvalues of `P`, whose sum is `n`. Concavity of square root and trace
against `Q` now give (1) for finite convex mixtures. Countable mixtures
follow by finite-dimensional continuity, for example by replacing the
omitted total weight with that weight times `I` and taking the limit.

For (2), let `b_i` be column `i` of `B`. Since `B^2=Q` and `Q_ii=1`,
`||b_i||_2=1`. Cauchy--Schwarz and positivity therefore give

`sqrt(v_i)=||sqrt(R)b_i||_2 >= b_i^T sqrt(R)b_i`.

Summing yields `sum sqrt(v_i)>=Tr[B sqrt(R)B]=Tr[Q sqrt(R)]>=n`.

## Variance-normalized gain and cutoff parameters

For any fixed `eta>0`, let

`d_i=1/sqrt(max(v_i,eta))`.

Then `d_i^2 v_i<=1`, `d_i<=eta^{-1/2}`, and

`n^{-1}sum_i d_i v_i >= 1-sqrt(eta)/4`.                          (3)

Indeed the difference between `sqrt(v)` and `v/sqrt(max(v,eta))` is zero
when `v>=eta`, and at most `sqrt(eta)/4` when `0<=v<=eta`; maximize
`x-x^2/sqrt(eta)` at `x=sqrt(eta)/2`. Combining with (2) proves (3).
The weaker floor `1-sqrt(eta)` is therefore also valid.

For the unmarked odd response `h(G)` with Gaussian norm one, its formal
Gaussian covariance kernel is exactly this mixture with weights equal
to its squared odd Hermite coefficients. Thus the theorem supplies a
deterministic lower bound for the average standard deviation whenever
the corresponding local Gaussian transport limit has been established.
It does not itself establish that probabilistic limit or the weighted
Boolean energy identity; those are separate theorems.
