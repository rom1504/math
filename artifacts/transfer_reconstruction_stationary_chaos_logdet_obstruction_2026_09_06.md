# Stationary centered chaos and two raw log-determinant lower-bound obstructions

Date: 2026-09-06. The stationary identity is exact. The proposed raw Gaussian
log-determinant lower bound is false under the condition
`beta C max_i v_i<1`, including at nonzero globally optimal product means
with `beta C>1`. Neither counterexample below settles an additional,
genuinely small-variance hypothesis. A variance-constrained Gaussian
replacement is recorded at the end as an unproved target, not a theorem.

## 1. Exact cancellation at an interior stationary product law

Let `A` be symmetric and hollow, `b=beta/sqrt N`, and

```math
P(m)=\sum_i h\!\left(\frac{1-m_i}{2}\right)+\frac b2m^TAm.
```

Suppose `m in (-1,1)^N` is stationary, so
`a_i=atanh(m_i)=b(Am)_i`. Under independent spins of means `m_i`, put
`xi_i=X_i-m_i`. Their product mass function is
`mu_m(X)=exp(sum_i a_i X_i)/prod_i 2 cosh(a_i)`. Expanding the quadratic
and cancelling the linear term gives exactly

```math
\log Z_A(b)-P(m)
=\log\mathbb E_{\mu_m}
   \exp\left\{\frac\beta{2\sqrt N}\xi^TA\xi\right\}.       \tag{1}
```

The constant cancellation uses
`H(mu_m)=sum log(2cosh a_i)-sum a_i m_i`.
Global optimality is unnecessary for (1), though both examples below use
globally optimal product means.

Write `D=diag(v_i)`, `v_i=1-m_i^2`, and
`J=beta D^{1/2} A D^{1/2}/sqrt N`. Replacing the centered spins by centered
Gaussians with these variances gives

```math
G_A(m)=-\frac12\log\det(I-J),                             \tag{2}
```

when `||J||op<1`. Existence of this Gaussian integral does **not** imply
`log Z_A-P(m)>=G_A(m)-o(N)`.

## 2. Explicit nonzero optimal means with `beta C>1`

For odd positive `r`, let

```math
N=4^r,\quad H=(J_4-2I_4)^{\otimes r},\quad A=H+I_N.
```

Then `A` is a hollow sign matrix, `H^2=NI`, `H1=sqrt N 1`, and `H`
has diagonal `-1`. In particular
`||A||op=sqrt N+1`, with the positive endpoint attained at `1`.
Set `beta=101/100`. The globally optimal product means are attained at
the constant vectors `m=+q_N 1` and `m=-q_N 1`, where

```math
\frac{\operatorname{atanh}q_N}{q_N}
=\beta(1+N^{-1/2}),\qquad 0<q_N<1.                       \tag{3}
```

To see global optimality, the function

```math
\varphi(s)=h\!\left(\frac{1-\sqrt s}{2}\right)
=\log2-\sum_{j\ge1}\frac{s^j}{2j(2j-1)}
```

is strictly concave for `0<s<1`. For any vector of means, Jensen and the
spectral upper bound give

```math
P(m)\le N\sup_{0\le s\le1}
\left\{\varphi(s)+\frac\beta2(1+N^{-1/2})s\right\}.
```

Constant vectors attain equality. Since the derivative at `s=0` is
positive and tends to negative infinity at `s=1`, the unique maximizing
squared mean is the solution in (3).

Let `q_N->q` and `v_N=1-q_N^2->v=1-q^2`. The expansion

```math
\frac{\operatorname{atanh}q}{q}
=\sum_{j\ge0}\frac{(q^2)^j}{2j+1}
```

gives the rational bounds

```math
\frac1{40}<q^2<\frac3{100},\qquad
\frac{97}{100}<v<\frac{39}{40}.                          \tag{4}
```

Indeed, at `q^2=1/40` the tail beyond the constant term is at most
`(1/40)/(3(1-1/40))=1/117<1/100`; at `q^2=3/100`, its first two
terms already exceed `1/100`. Also `q_N>q`, hence `v_N<v`.

Take `C=101/100`; the operator bound holds once `N>=10000`. Nevertheless

```math
\beta C>1,\qquad
\beta C v_N<\frac{10201}{10000}\frac{39}{40}
=\frac{397839}{400000}<1.                               \tag{5}
```

Thus the stated weighted-norm condition holds with fixed strict slack.

The proportions of the two eigenvalues of `H/sqrt N` tend to one half,
so the raw Gaussian quantity satisfies

```math
\lim_{N\to\infty}\frac{G_A(q_N1)}N
=-\frac14\log(1-\beta^2v^2)
> -\frac14\log\left(\frac{4018791}{100000000}\right)
>\frac45.                                              \tag{6}
```

The Gaussian integral is well-defined: if `a_N=beta(1+N^{-1/2})`,
(3) and the integral representation of `atanh(q_N)/q_N` give
`a_N v_N<1` strictly. This also bounds the other spectral endpoint.

On the other hand, the spectral cap gives
`log Z_A(b)<=N log2+beta(1+N^{-1/2})N/2`, while `P(q_N1)>=N log2`.
Consequently

```math
\limsup\frac{\log Z_A(b)-P(q_N1)}N
\le\frac\beta2=\frac{101}{200}.                         \tag{7}
```

Equations (6)--(7) violate the raw Gaussian lower bound by more than
`N/4` for all sufficiently large members of this explicit sequence.
This is a hollow-sign, bounded-norm, nonzero-optimal-product counterexample,
not merely a mismatch of low-order formal cumulants.

## 3. Failure even in strict ordinary high temperature at `m=0`

There is also a simple existence counterexample with
`beta=1/100`, `C=10`, hence `beta C<1/2`.

Sample the off-diagonal signs of a hollow symmetric `A` independently.
For every fixed spin,

```math
\mathbb E_A e^{\pm(\beta/\sqrt N)H_A(x)}
=\cosh(\beta/\sqrt N)^{N(N-1)/2}.
```

Therefore Markov's inequality applied to the sum of the two partitions
shows that, with probability tending to one, simultaneously

```math
\log Z_A(\pm\beta/\sqrt N)-N\log2
\le \frac{\beta^2(N-1)}4+\varepsilon N+\log2,            \tag{8}
```

for any fixed `epsilon>0`. Only `log cosh z<=z^2/2` is used.

The operator norm is at most `10 sqrt N` with probability tending to one.
Here is an elementary bound: for every unit vector `u`, the independent
sum `u^TAu=2sum_{i<j}A_ij u_i u_j` has squared coefficients summing to
at most two, so its moment generating function gives
`P(|u^TAu|>=z)<=2exp(-z^2/4)`. A `1/4` net of the unit sphere has at most
`9^N` points, and the symmetric quadratic-form norm is at most twice its
maximum on that net. The probability of norm above `10 sqrt N` is thus
at most `2exp[(log9-25/4)N]`.

Select a realization satisfying both properties. At either objective sign,
the globally optimal product mean is zero, because

```math
P(m)\le N\log2-\frac{1-\beta C}{2}\|m\|^2.
```

Put `B=A/sqrt N`. Averaging the two Gaussian quantities removes odd
traces, and gives

```math
\frac{G_A(0)+G_{-A}(0)}2
=-\frac14\log\det(I-\beta^2B^2)
\ge\frac{\beta^2}4\operatorname{tr}B^2
   +\frac{\beta^4}8\operatorname{tr}B^4
\ge\frac{\beta^2(N-1)}4
   +\frac{\beta^4(N-1)^2}{8N}.                          \tag{9}
```

The last line uses `tr B^2=N-1` and Cauchy--Schwarz for its eigenvalues.
Choose the objective sign with the larger Gaussian quantity, and take
`epsilon=beta^4/32` in (8). Equations (8)--(9) show a discrepancy at least
`beta^4 N/16` for all sufficiently large orders. Thus even strict ordinary
high temperature cannot justify the raw Gaussian lower comparison.

## 4. A variance-constrained Gaussian target, still unproved for spins

The raw Gaussian integral allows every coordinate variance to change.
A different functional fixes those variances and subtracts the Gaussian
entropy cost of creating correlations. For symmetric hollow `J`, define

```math
\mathcal S(J)=\sup_{R\succ0,\ \operatorname{diag}R=1}
\left\{\frac12\operatorname{tr}(JR)+\frac12\log\det R\right\}.
                                                                  \tag{10}
```

Its elementary Gaussian covariance dual is

```math
\mathcal S(J)=\inf_{\Lambda\text{ diagonal},\ \Lambda-J\succ0}
\frac12\{\operatorname{tr}\Lambda-N-\log\det(\Lambda-J)\}.
                                                                  \tag{11}
```

This follows from maximizing over positive definite `R` after introducing
diagonal Lagrange multipliers; at the optimum
`R=(Lambda-J)^{-1}` and its diagonal equals one. The covariance feasible
set is bounded, has the interior point `I`, and its logarithm tends to
negative infinity at the boundary, so the finite optimum and duality cause
no compactification issue.

For a flat symmetric spectrum with endpoints `+-a` and asymptotically
constant resolvent diagonal, (11) instead gives the scalar expression

```math
\frac{\mathcal S}N=\frac{\lambda-1}{2}-\frac14\log\lambda,
\qquad \lambda=\frac{1+\sqrt{1+4a^2}}2.
```

Its small-`a` expansion is `a^2/4-a^4/8+O(a^6)`, unlike the raw
Gaussian expression `a^2/4+a^4/8+O(a^6)`. This removes the specific
radial/backtracking discrepancy exposed above.

It is **not proved here** that the spin chaos in (1) dominates
`S(beta D^{1/2} A D^{1/2}/sqrt N)-o(N)`, even under genuinely small
variances. Passing this spectral sanity check is not a discrete-spin
variational construction. The existing fixed-marginal clipped-bipartite
law gives a rigorous first correlation credit at current large `beta`
and tiny variance, but not the resummation (10).

Reproduction: `computations/transfer_reconstruction_stationary_logdet_exact_checks_2026_09_06.py`
checks the rational root brackets and the strict determinant discrepancy,
and enumerates a finite stationary example to check the exact identity (1).
