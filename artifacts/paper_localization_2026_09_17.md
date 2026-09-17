# Gaussian localization: reconstruction and exact-sign transfer audit

2026-09-17. Main paper: El Alaoui--Montanari,
[*An Information-Theoretic View of Stochastic Localization*, v2](https://arxiv.org/abs/2109.00709).
Original theorem: Eldan,
[*Taming correlations through entropy-efficient measure decompositions with
applications to mean-field approximation*](https://arxiv.org/abs/1811.11530).

Status: the primary proof and the new obstruction below have been reconstructed.
The root agent independently derived the central entropy and greedy inequalities.
This is not a convergence/nonconvergence proof and does not improve the original
numerical interval. The initial mechanism was frozen before reading the earlier
BH synthesis, in `tmp/paper_portfolio_2026_09_17/localization/initial_mechanism_frozen.md`.

## 1. Exact primary theorem, with safe hypotheses

The finite-support case suffices for both spin and coefficient-sign uses. Let X
have law mu on a finite subset of R^d, mean m_0, and covariance Sigma. Fix a
positive definite precision matrix L, take tau uniform on [1,2], independent of
X and a standard Gaussian Z, and observe

```math
Y=\sqrt\tau X+L^{-1/2}Z,\qquad \Theta=(\tau,Y).
```

Write mu_theta for the posterior law and C_theta for its covariance. Theorem 1
of the main paper, with its noise matrix Q=L^{-1}, gives

```math
\mu=\mathbb E\mu_\Theta,\quad
\mathbb E C_\Theta\preceq L^{-1},\quad
I(X;\Theta)\leq\tfrac12\log\det(I+2L^{1/2}\Sigma L^{1/2}),\quad
\mathbb E[C_\Theta L C_\Theta]\preceq\Sigma.                 \tag{1}
```

The support statement is exact: every posterior is supported on the original
finite support. It does NOT say that independent samples from its coordinate
marginals remain in that support. Nor does small covariance imply proximity to
a product law in total variation: the theorem deliberately controls a weaker
quantity. The information cost is mutual information, not the differential
entropy of the continuous observation or the logarithm of a finite component
count.

The positive-definite formulation avoids inverse ambiguities. Singular
precisions are handled by observing the corresponding image of X, or by a
limit. In contrast, singular noise can reveal a projection exactly and must
not be treated by an unqualified ordinary inverse.

### 1.1 Decisive proof, reconstructed

Use the observation process Z_t=t L^{1/2}X+B_t. At fixed t its terminal value
is a sufficient statistic for the whole path. Bayes' formula gives

```math
\mu_t(x)=\frac{\mu(x)\exp(\langle L^{1/2}Z_t,x\rangle
                  -t x^T Lx/2)}{\sum_u\mu(u)\exp(\langle L^{1/2}Z_t,u\rangle
                  -t u^TLu/2)}.
```

Let m_t be its mean and let
W_t=Z_t-integral_0^t L^{1/2}m_s ds. This is a Brownian motion in the observation
filtration: its conditional increments have zero drift and quadratic variation
t I. Direct Itô differentiation of the finite Bayes sum gives

```math
d\mu_t(x)=\mu_t(x)(x-m_t)^T L^{1/2}dW_t,
\qquad dm_t=C_t L^{1/2}dW_t.                              \tag{2}
```

All integrands are bounded on a finite time interval, so these are genuine
martingales. Since the posterior second moment is a martingale, Itô isometry
applied to m_t m_t^T yields the crucial matrix identity

```math
\frac d{dt}\mathbb E C_t=-\mathbb E[C_t L C_t].             \tag{3}
```

Integrating over [1,2] proves the last assertion in (1), because
E C_1-E C_2 <= Sigma in Loewner order. The covariance assertion follows even
more directly: the Bayes estimator is at least as good, in every positive
semidefinite quadratic loss, as Y/sqrt(t), whose error covariance is L^{-1}/t.
Thus E C_t <= L^{-1}/t for every t>0; averaging gives the slightly sharper
factor log(2) on [1,2].

Finally,
I(X;Y|tau=t)=h(Y|tau=t)-h(Y|X,tau=t). The second term is Gaussian entropy, and
the first is at most the Gaussian entropy with covariance t Sigma+L^{-1}.
This gives the log determinant in (1), since t<=2. Equivalently, differentiating
posterior entropy using (2) gives the finite-support I-MMSE identity

```math
\frac d{dt}I(X;Z_{[0,t]})=\tfrac12\mathbb E\operatorname{Tr}(L C_t). \tag{4}
```

This reconstructs the decisive argument of both papers. The main paper's
Gaussian-channel proof obtains the improved entropy constant directly; Eldan's
original proof derives (3)--(4), bounds one-dimensional covariance quadratic
forms by the reciprocal solution of f'=-f^2, and integrates over [0,2].

Audit cautions: the main paper's remark about eliminating tau by a single
matrix-valued Markov inequality is not used here; such an inference requires
care beyond the scalar Markov inequality. In the original paper's mean-field
proof, the displayed invocation `L=|J|^(1/2)` does not match the precision
convention of its Theorem 2; the subsequent displayed inequalities require
L=|J|. This is an easily repaired parameter label, not a failure of the theorem.

## 2. Exact spin transfer and the normalization obstruction

For a hollow symmetric full sign matrix A, put H_A(x)=x^T A x/2 and
Q(A)=max_{x in {+/-1}^n}|H_A(x)|. For ANY spin prior and ANY observation,

```math
\mathbb E H_A(X)=\mathbb E H_A(m_\Theta)
                    +\tfrac12\mathbb E\operatorname{Tr}(A C_\Theta). \tag{5}
```

Independently rounding the coordinates of m_theta into signs preserves
H_A(m_theta) in conditional expectation because A has zero diagonal. It does
not preserve the original spin support or its near-extremality. The localization
bound (1), by weighted Hilbert--Schmidt Cauchy--Schwarz, implies

```math
\mathbb E|\operatorname{Tr}(A C_\Theta)|
\leq\sqrt{\operatorname{Tr}(A L^{-1}A)\operatorname{Tr}\Sigma}. \tag{6}
```

Indeed tr(AC)=<L^{-1/2}A,L^{1/2}C>_HS, and
E||L^{1/2}C||_HS^2=E tr(CLC)<=tr Sigma. For L=rI and tr Sigma<=n, the energy
loss certificate is n sqrt(n-1)/(2 sqrt(r)). The generic information certificate
is n log(1+2r)/2, by concavity and tr Sigma<=n. The two certificates do not give
o(n) information together with o(n^{3/2}) energy error. This observation alone
would only be a limitation of the estimates, not an impossibility theorem.

The gradient-complexity issue is exact, not conjectural. For
f_beta(x)=beta H_A(x)/sqrt(n), beta>=0, and g standard Gaussian,

```math
\mathbb E\sup_{x\in\{-1,1\}^n}\langle g,\nabla f_\beta(x)\rangle
=\frac\beta{\sqrt n}\mathbb E\|Ag\|_1
=\beta\sqrt{\frac2\pi}\,n\sqrt{1-\frac1n}.                \tag{7}
```

Each row has Euclidean norm sqrt(n-1). Thus the Gaussian width is Theta(n) at
fixed beta for EVERY full signing. A low-gradient-complexity result cannot be
applied at this normalization without a genuinely new restricted-gradient
hypothesis. Replacing sqrt(n) by n makes the width small but changes the energy
scale.

Eldan's mean-field theorem controls a product-measure free-energy deficit by
3 log det(I+Sigma |J|). With J=beta A/(2sqrt(n)), the natural generic bound is
order n at fixed beta, not o(n). Sending beta to infinity and allowing
deterministic product measures simply recovers the original maximum; it does
not by itself supply a composition or convergence theorem.

## 3. New exact-minimizer obstruction to low-information edge rounding

This section is a deduction for this project, not a theorem attributed to the
papers and not a claim of external novelty.

### 3.1 The prior and all allowed observations

Fix ANY hollow full sign matrix A_0 of order n. Take independent uniform signs
xi_1,...,xi_n and define a random coefficient signing

```math
X_{ij}=(A_0)_{ij}\xi_i\xi_j,\qquad i<j.                     \tag{8}
```

Every point of this prior has cap Q(A_0), by vertex switching. Its edge means
are zero, edge covariance is the identity of dimension binom(n,2), and entropy
is (n-1)log(2). Let Y be ANY randomized observation of X, with information cost
I=I(X;Y), measured in nats. This includes arbitrary adaptive/anisotropic Gaussian
observations, exact revealed edges, and any mixture of such operations.

Conditional on Y=y, set b_ij=E[X_ij|y]. Independently draw coefficient signs
R_ij with these means. R is still a hollow FULL sign matrix, but its edges no
longer satisfy the correlations of the original switching orbit.

### 3.2 Entropy controls the squared posterior biases

Lift the posterior back to xi in {+/-1}^n, assigning equal probability to the
two global-flip preimages of each coefficient signing. This is exactly the
posterior of the original uniform xi, since Y depends only on X. Every marginal
xi_i remains uniform. Write

```math
D_y=D(\mathcal L(\xi|y)\|U_n)=n\log2-H(\xi|y),\qquad
c_{ij}=\mathbb E[\xi_i\xi_j|y].
```

For each fixed i, conditional entropy subadditivity gives

```math
H(\xi|y)\leq H(\xi_i|y)+\sum_{j\ne i}H(\xi_j|\xi_i,y)
 =n\log2-\sum_{j\ne i}I_y(\xi_i;\xi_j).
```

The joint law of two unbiased signs with correlation c has mutual information
[(1+c)log(1+c)+(1-c)log(1-c)]/2 >= c^2/2. Summing over i and dividing by two
therefore gives the POINTWISE inequality

```math
\sum_{i<j}b_{ij}^{\,2}=\sum_{i<j}c_{ij}^{\,2}\leq nD_y.
\quad\text{Consequently}\quad
\mathbb E\sum_{i<j}b_{ij}^{\,2}\leq nI.                   \tag{9}
```

Here E D_y=I(xi;Y)=I(X;Y), since X is a deterministic function of xi and
Y is conditionally independent of xi given X. No product-posterior premise is
used, and no particular choice of A_0 matters.

### 3.3 A uniform greedy lower bound for biased independent edges

Let b be an arbitrary deterministic edge-bias array. Set s_1=1 and, in order,
choose s_i as the sign of sum_{j<i} R_ij s_j (break a zero arbitrarily). Then

```math
H_R(s)=\sum_{i=2}^n\left|\sum_{j<i}R_{ij}s_j\right|.        \tag{10}
```

Given all previously used edges, the current row edges are independent, and
its centered sum has variance v_i=(i-1)-sum_{j<i}b_ij^2. Its conditional mean
may be arbitrary. Uniformly in that mean and all biases,

```math
\mathbb E\left|\sum_{j<i}R_{ij}s_j\right|
\geq\kappa\sqrt{v_i}-2(i-1)^{1/3},\qquad
\kappa=\sqrt{2/\pi}.                                    \tag{11}
```

Self-contained proof: smooth absolute value by phi_a(z)=sqrt(z^2+a^2).
Its third derivative has absolute value <=a^{-2}. Replace each centered edge
by an independent Gaussian with matching variance. Taylor expansion to order
two and cancellation of mean/variance terms bound the total replacement error
by (2+sqrt(8/pi))k/(6a^2), where k=i-1: centered signs have third absolute moment
at most 2 times their variance, and centered Gaussians have third absolute
moment sqrt(8/pi) times variance^{3/2}. Also
|z|<=phi_a(z)<=|z|+a. A shifted Gaussian has expected absolute value minimized
at zero, so it contributes at least kappa sqrt(v_i). Choose a=k^{1/3} and note
1+(2+sqrt(8/pi))/6<2. This proves (11), including highly biased or zero-variance
rows, with no asymptotic normality assumption.

Write D(b)=sum_{i<j}b_ij^2. Since
sqrt(k)-sqrt(k-d)<=sqrt(d), summation and Cauchy--Schwarz give

```math
\mathbb E_R Q(R)
\geq\frac{2\kappa}{3}(n-1)^{3/2}
 -\kappa\sqrt{(n-1)D(b)}-\frac32 n^{4/3}.                 \tag{12}
```

The integral comparisons used here are
sum_{k=1}^{n-1}sqrt(k)>=(2/3)(n-1)^{3/2} and
sum_{k=1}^{n-1}k^{1/3}<=(3/4)(n^{4/3}-1).

Combining (9), (12), and Jensen gives the central transfer obstruction:

```math
\boxed{\quad
\mathbb E Q(R)\geq\frac{2\kappa}{3}(n-1)^{3/2}
   -\kappa\sqrt{n(n-1)I}-\frac32 n^{4/3}.
\quad}                                                   \tag{13}
```

In particular, I=o(n) implies

```math
\liminf\frac{\mathbb E Q(R)}{n^{3/2}}
\geq\gamma:=\frac23\sqrt{\frac2\pi}=0.5319230405\ldots.   \tag{14}
```

### 3.4 High-probability statement, not just an average

If I=o(n), choose eta_n down to zero with I/(n eta_n)->0, for example
eta_n=sqrt(I/n)+n^{-1/10}. Markov's inequality gives D_y<=n eta_n with probability
1-o(1). For every such y, (9) and (12) put the conditional expected cap above
[gamma-kappa sqrt(eta_n)-o(1)]n^{3/2}. Conditional on y, changing one rounded
edge changes Q by at most two. The bounded-differences lower-tail bound is

```math
\mathbb P\{Q(R)<\mathbb E[Q(R)|y]-t\mid y\}
\leq\exp\{-t^2/[2\tbinom n2]\}.
```

Taking t=n^{5/4} proves Q(R)>=[gamma-o(1)]n^{3/2} with probability 1-o(1) under
the JOINT observation-and-rounding law.

### 3.5 Application to actual exact minimizers

At every order n=2^k, k>=1, let S_n be the symmetric Sylvester Hadamard matrix.
It has S_n^2=nI and trace zero. Its hollowing A=S_n-diag(S_n) is a full signing,
and for Boolean x,

```math
H_A(x)=\tfrac12x^TS_nx-\tfrac12\operatorname{Tr}S_n,
\qquad Q(A)\leq\tfrac12n^{3/2}.                          \tag{15}
```

Therefore an EXACT minimizer A_0 at any such order satisfies
Q(A_0)=M_n<=n^{3/2}/2. Take the prior (8) over its switching orbit. Every original
sample and every correlated posterior sample is still an exact minimizer.
But every o(n)-information channel followed by independent marginal edge
rounding has, with probability tending to one,

```math
Q(R)-M_n\geq(\gamma-\tfrac12-o(1))n^{3/2}
 =(0.0319230405\ldots-o(1))n^{3/2}.                       \tag{16}
```

This uses no reported decimal frontier or unverified minimizer structure. The
same failure already holds for the explicit Hadamard switching orbit if one
prefers a constructive prior. It defeats a universal entropy-efficient
independent-product rounding theorem even on priors supported entirely on
actual exact minimizers.

Scope: (16) does NOT rule out retaining a correlated posterior sample,
information cost of order n, a theorem selecting an exceptionally rare
posterior component, or a construction changing the prior/support. It is not
an obstruction to the original convergence problem itself.

## 4. Optional sharper information--bias tradeoff

The bias loss in (12) can be sharpened without improving the Gaussian
replacement error. Since
sqrt(k)-sqrt(k-d)<=d/sqrt(k), maximize sum_k d_k/sqrt(k) under 0<=d_k<=k and
sum d_k=D. A fractional-knapsack exchange fills the smallest k first. If
D=m(m+1)/2+r with 0<=r<m+1, the result is

```math
\sum_k\bigl(\sqrt k-\sqrt{k-d_k}\bigr)
\leq\sum_{k=1}^m\sqrt k+\frac r{\sqrt{m+1}}
\leq\frac23(2D)^{3/4}+O(D^{1/4}+1).                     \tag{17}
```

The last bound follows by integral comparison; its constant is absolute.
Since D<=nD_y and E D_y=I, concavity yields

```math
\frac{\mathbb E Q(R)}{n^{3/2}}
\geq\gamma\left[1-\left(\frac{2I}{n}\right)^{3/4}\right]-o(1), \tag{18}
```

uniformly for I<=n log 2. Thus achieving expected normalized cap at most
1/2+o(1) on the minimizer-orbit prior requires

```math
\liminf\frac In\geq
\frac12\left(1-\frac{1}{2\gamma}\right)^{4/3}>0.           \tag{19}
```

The simpler bound (13) is the principal audited statement; (17)--(19) are an
optional strengthening with the same scope.

## 5. Positive extension: covariance-adapted scalar water filling

The following finite-support theorem is a proved deduction, not a claim of
external novelty. It uses the linear-MMSE comparison behind the primary paper,
but targets one signed quadratic expectation rather than every covariance
direction.

Let mu be ANY spin law, with mean m and covariance Sigma, and set
B=Sigma^{1/2} A Sigma^{1/2}. List its positive eigenvalues as lambda_1,...,lambda_r.
For every theta>0 there is a Gaussian observation Y of the original spins,
whose posteriors have the exact original support, such that

```math
I(X;Y)\leq\frac12\sum_{j=1}^r\log_+(\lambda_j/\theta),
\qquad
\mathbb E H_A(X)-\mathbb E H_A(m_Y)
\leq\frac12\sum_{j=1}^r\min(\lambda_j,\theta).             \tag{20}
```

Here log_+(u)=max(log u,0). Independent rounding of the SPINS with means m_Y
produces a sign law with expected Hamiltonian E H_A(m_Y). Thus (20) is an
energy-preserving productization theorem whenever its positive-covariance
spectral profile is favorable. It is not edge-coefficient rounding and does
not contradict Section 3.

Proof. Work on the range of Sigma and whiten U=Sigma^{dagger/2}(X-m), so that
Cov(U)=I there. Let K share eigenvectors with B, giving a positive eigenvector
of B eigenvalue k_j=(lambda_j/theta-1)_+, and zero on all other directions.
Observe Y=K^{1/2}U+Z with independent standard Gaussian noise. The Gaussian
entropy upper bound gives I(U;Y)<=log det(I+K)/2. The Bayes estimator is at
least as good in every positive semidefinite quadratic loss as the best LINEAR
estimator. The latter has error covariance (I+K)^{-1}, whence

```math
\mathbb E\operatorname{Cov}(U|Y)\preceq(I+K)^{-1}.
```

Since the posterior covariance is positive semidefinite,
tr(B C)<=tr(B_+ C), even when B is indefinite. Applying the matrix inequality
to B_+ gives the risk bound in (20), because
lambda_j/(1+k_j)=min(lambda_j,theta). Formula (5) completes the proof.
Pseudoinverses cause no loss: X-m lies in the range of Sigma almost surely.

For a prescribed risk budget d, the k_j above also minimize the Gaussian
log-determinant certificate sum log(1+k_j)/2 subject to
sum lambda_j/(1+k_j)<=d. This follows either by a Lagrange multiplier or by
the scalar water-filling inequalities, with inactive coordinates k_j=0.
This optimality is only for these spectral Gaussian certificates, not for all
possible observations of a discrete prior.

For example, a sufficient condition for o(n)-information and o(n^{3/2}) scalar
loss is a sequence theta_n for which

```math
\sum_j\min(\lambda_j,\theta_n)=o(n^{3/2}),\qquad
\sum_j\log_+(\lambda_j/\theta_n)=o(n).                    \tag{21}
```

The signs and normalization are important: it is the positive spectrum of
Sigma^{1/2} A Sigma^{1/2}, not automatically the spectrum of A or |A|. For the
negative extremum, apply the theorem to -A. For a prior supported on
H_A>=M-delta n^{3/2}, the output expected energy is at least
M-delta n^{3/2}-the risk in (20).

No appropriate canonical near-extremal prior satisfying (21) is presently
constructed for exact minimizers. Choosing a prior concentrated on an already
known maximizing spin trivially satisfies the condition and does not solve the
original optimization or convergence problem. For a genuinely flat positive
covariance spectrum, (21) fails at the desired scales. Thus (20) isolates a
testable extra hypothesis; it does not silently assume one.

## 6. Adaptive extension and surviving open target

The natural extension is to retain the posterior correlations and control a
specific scalar objective, rather than replace the posterior by a product.
For an adaptive observation dZ_t=L_t^{1/2}X dt+dB_t, with predictable bounded
positive-semidefinite L_t, the same finite-support filtering calculation gives

```math
dm_t=C_tL_t^{1/2}dW_t,\qquad
\frac d{dt}\mathbb E C_t=-\mathbb E[C_tL_tC_t],\qquad
\frac d{dt}I(X;Z_{[0,t]})=\tfrac12\mathbb E\operatorname{Tr}(L_tC_t).
```

For a fixed spin Hamiltonian the rate of change of its covariance contribution
is -(1/2)E tr(A C_t L_t C_t). A rank-one choice v v^T has instantaneous
energy-improvement/information ratio

```math
\frac{v^TC_t A C_t v}{v^TC_t v},
```

the Rayleigh quotient of C_t^{1/2} A C_t^{1/2}. This supplies an exact control
problem, but no uniform favorable rate is presently proved. In particular,
it does not bypass (16) for edge-marginal independent rounding.

For a spin prior supported on a single maximizing level H_A=M, the posterior
still has that support, and its independent-spin rounding deficit is exactly
M-H_A(m_theta)=(1/2)tr(A C_theta)>=0. There is no signed cancellation in this
specific case. At zero deficit, the posterior product law must itself be
supported on maximizers. Because A has every off-diagonal entry nonzero, a
product law supported on maximizers has at most ONE nondeterministic coordinate:
two free coordinates would give four assignments on which the nonzero bilinear
coefficient cannot be constant. Consequently exact, zero-loss productization
of a ground-state prior requires I>=H(mu)-log2. This exact statement has no
useful o(n^{3/2}) stability estimate yet.

The remaining candidate is therefore an exact correlated coupling that keeps
the sign support and extremal scalar information, or an optimizer-specific
near-extremal spin compression theorem. Neither has been proved here. A generic
small-covariance or low-information-to-independent-rounding bridge is now
falsified quantitatively on the actual minimizer class.

## 7. Further transfer audits

### 7.1 Universal random-model lower envelope for independent rounding

A second, independently derived transfer bound strengthens the qualitative
content of (13). Let R_0 have independent unbiased coefficient signs and R_b
have independent signs of means b. For every beta>0, with D(b)=sum b_e^2,

```math
\mathbb E Q(R_b)\geq\mathbb E Q(R_0)
 -\frac{n^{3/2}\log2}{\beta}
 -\frac{\beta D(b)}{2\sqrt n}-\beta^2 n.                  \tag{22}
```

Here is a direct proof independent of the greedy algorithm. Use a soft maximum
over the 2^n indexed witnesses sigma chi_x, where x_1=1 and sigma=+/-1:
F(z)=sqrt(n) log sum exp(beta <z,sigma chi_x>/sqrt(n))/beta. Then
Q<=F<=Q+n^{3/2}log2/beta. Each second coordinate derivative is at most
beta/sqrt(n); each third derivative is at most beta^2/n because its remaining
factor is the third central moment of a +/-1 variable, with magnitude below
one. Replace the centered biased coordinates by matching Gaussians and later
replace full-variance Gaussians by unbiased signs. The sum of the two Lindeberg
errors is at most beta^2 n, by the same third-moment calculation as in (11).
Restoring the variance lost in edge e costs at most beta b_e^2/(2sqrt(n)), by
Gaussian interpolation. Finally, deterministic mean b cannot decrease the
expected even convex soft maximum under centered Gaussian noise: convexity at
the midpoint and Gaussian sign symmetry prove this. These facts give (22).

Averaging (22) over a channel and using (9) shows that I=o(n) forces

```math
\mathbb E Q(R)\geq\mathbb E Q(R_0)-o(n^{3/2}).             \tag{23}
```

For example take beta=min((n/I)^{1/2},n^{1/8}), with the first term omitted
when I=0. Thus the full random-sign model, not just the elementary constant
gamma, is an asymptotic lower envelope. No Parisi formula or numerical random
model constant is needed. The root agent independently obtained a stronger
uniform version with loss sqrt(2n log2 D(b))+O(n^{4/3}) by coupling missing
Gaussian variance and applying a union bound; its proof and the extension to
correlated/mean-retuned output laws belong to the root's separate artifact.

### 7.2 The rare-component loophole is real

An existential conclusion about one posterior component is insufficient:
reveal the entire finite input with probability epsilon and reveal nothing
otherwise. The information cost is epsilon H(X), yet the rare revealed
components are point masses and independently round with zero loss. Therefore
the obstruction in Section 3 is deliberately an expectation/high-probability
statement; it must not be promoted to a claim about every component. A positive
decomposition-to-optimization theorem needs a quantitative selection-weight or
componentwise-cost condition.

### 7.3 Fixed-law Bernoulli geometry still forgets the signing on an orbit

After freezing the Gaussian reconstruction, the definitions of Liu--Zadik,
[*A Bayesian Proof of the Bernoulli Theorem*](https://arxiv.org/abs/2608.11031),
were checked for a possible combination. Let mu_A be the uniform switching
orbit law in (8). Multiplication by A is a coordinatewise sign isometry sending
mu_A to the uniform cut-sign law. Hence the Gaussian/Cauchy channel geometry,
capped-quadratic rate--distortion functionals, and prescribed-law Gaussian and
Bernoulli maximum-coupling functionals are IDENTICAL for all full signings A.

In fact the Bernoulli fixed-law functional of mu_A equals
E_epsilon max_x H_epsilon(x). Choose a maximizer under independent edge noise,
randomize ties equivariantly under vertex switching, and use that this action
is transitive on the orbit: the maximizing index has exactly the uniform orbit
law. This attains the trivial upper bound by the expected supremum. On the
uniform union of the positive and negative switching orbits the analogous
value is E Q(epsilon). These observations concern the paper's definitions and
elementary symmetry, not a reconstruction of its full proof, which belongs to
the Bernoulli track.

A useful combined theorem must therefore retain the orbit's position relative
to the FIXED all-ones direction, or another scalar extremal constraint. Isometry-
invariant channel geometry alone cannot distinguish the minimizing signing
from an arbitrary full signing.

## 8. Scalar noise-response audit

For independent signs epsilon_e of common mean t, write
F_A(t)=E Q(A times epsilon), with entrywise multiplication. This is the Boolean
noise semigroup T_t applied to the coefficient-cube function Q. It retains the
fixed scalar cap that coordinatewise-isometry-invariant widths forget.

### 8.1 Independent exact checks through order seven

The coefficient-cube Walsh support vanishes exactly outside Eulerian edge sets,
by vertex-switching invariance, and outside even edge counts, by Q(-A)=Q(A).
Thus the first possible nonconstant degree is four, where only four-cycles
occur. The same coefficient is assigned to every four-cycle by permutation
invariance. Exact full-edge enumeration independently confirms the root's
gauge-reduced census through n=7. All exact-minimizer noise polynomials through
n=7 are nonincreasing on [0,1]; this was checked by rational real-root isolation
of their derivatives, not merely by a numerical grid. For n=4,5,6 they are

```math
F_4(t)=\frac92-\frac12t^4,
\quad F_5(t)=\frac{51}8-\frac{15}8t^4-\frac58t^6+\frac18t^{10},
```

```math
F_6(t)=\frac{585}{64}-\frac{195}{64}t^4-\frac{35}{16}t^6
       +\frac{15}{64}t^8+\frac{15}{16}t^{10}-\frac5{64}t^{12}.
```

At n=7 there are exactly three distinct minimizer noise polynomials. Their
four-cycle sums are -15,-15,-7; all three decrease. Their complete rational
coefficients are in `noise_polynomials_n7.json`. Mixed signs in higher-degree
coefficients already rule out an argument requiring every nonconstant Fourier
level to be nonpositive. These are finite diagnostics only.

### 8.2 A complete small graph counterexample to generic monotonicity

Fixed-minimizer monotonicity does not follow from convexity, a quadratic spin
cap, vertex-switching symmetry, or the absence of Fourier degrees below four.
Take a theta graph consisting of three internally disjoint paths of lengths
2,2,4 between the same two endpoints. It has seven vertices and eight edges.
The three simple cycles have lengths 4,6,6. Give one edge of the long path sign
-1 and all other edges sign +1. The four-cycle is positive and both six-cycles
are negative.

The graph's cut code has codimension two. Every nonzero cycle syndrome can be
corrected by changing one edge (each nonzero parity-check column occurs), so
every unbalanced signing has cap 8-2=6 and every balanced signing has cap 8.
In particular the chosen signing is an EXACT minimizer. The graph is bipartite,
so global reversal is a vertex switching and causes no separate issue.

The indicator of a balanced signing is one quarter of one plus the products
around the three cycles. Independent edge noise therefore gives EXACTLY

```math
F_A(t)=\frac{13}{2}+\frac12t^4-t^6,
\qquad F_A'(t)=2t^3(1-3t^2).                             \tag{24}
```

It increases up to t=1/sqrt(3), where its value is 176/27>13/2, then decreases
to the minimum cap six. This is a sparse-graph counterexample, NOT a full
complete-graph counterexample. The two other nonzero syndromes have decreasing
curve 13/2-t^4/2, so it also does not defeat existence of some favorably selected
minimizer branch. Completeness/density or another original-class hypothesis is
essential for any stronger monotonicity claim.

The first random-code counterexample, its refinement to dual distance four,
and this minimized length-eight certificate are all retained in the temporary
working directory. The graph proof above does not depend on the search.

### 8.3 A positive monotonic envelope and its precise limit obligation

For EVERY real function f on a finite cube, its optimized noise envelope
g(t)=min_a T_t f(a) is nonincreasing for t in [0,1]. Indeed, for 0<=s<=t,
T_s=T_{s/t}T_t and a Markov average is at least the global minimum of its
argument. This does not imply monotonicity at one fixed minimizer.

For the cap, let g_n(t)=min_A F_A(t), so g_n(1)=M_n. An elementary Bernstein
union bound applied to the centered edge noise gives uniformly in A

```math
F_A(t)\leq t Q(A)+C\sqrt{1-t^2}\,n^{3/2}+Cn.
```

To check the scale, the centered residual has independent coordinates, absolute
bound two, and variance 1-t^2. Every signed spin witness therefore has variance
binom(n,2)(1-t^2), and there are at most 2^n witnesses. The maximum has expected
value at most C[sqrt(n binom(n,2)(1-t^2))+n]. Choosing an exact minimizer gives

```math
0\leq g_n(t)-M_n\leq C\sqrt{1-t}\,n^{3/2}+Cn.             \tag{25}
```

Consequently convergence of g_n(t_j)/n^{3/2} for any fixed sequence t_j up to one
would imply convergence of M_n/n^{3/2}, by sending j to infinity last. Uniform
convergence in t is unnecessary. This is a quantitatively correct regularized
proxy, but its cross-order convergence is NOT proved; without that input it is
not a strict reduction of the original difficulty.

## 9. Provenance and reproducibility

Primary PDF and TeX sources are retained under
`tmp/paper_portfolio_2026_09_17/localization/`. The initial frozen mapping,
failed extraction attempt, and all later exploratory work are preserved there.
The earlier artifact consulted after freezing was
`artifacts/bh_2026_09_16_final_synthesis.md`; no result from its numerical frontier
is used in the new obstruction. No ledger, STEERING, or ACTIVE_STATE edits and
no commits were made by this agent.

The independent numerical checks are implemented in
`computations/paper_localization_2026_09_17_checks.py` and retained in
`tmp/paper_portfolio_2026_09_17/localization/checks.json`, seed 20260917. They
include 900 random posterior entropy checks, exact small switching/cap checks,
finite biased-row enumerations, and a 20,000-sample anisotropic Gaussian
posterior diagnostic. The latter is numerical only; the samplewise quadratic
identity had error below 1.7e-15. The sharper information lower constant in
(19) is 0.011748360567219888 nats per spin. No simulation is a proof input.

## 10. Later joint extensions from the same primary mechanism

The original reconstruction and scalar-noise work above remain frozen.
Subsequent independently audited results are maintained separately:

- [Bounded-cap critical-block embedding](paper_localization_bounded_cap_embedding_2026_09_17.md):
  matching block covariance and a uniform subGaussian bound do not force
  o(N^(3/2)) cap universality at linear block size, even for actual full
  sign parents with expected cap O(N^(3/2)). The family is not minimizing.
- [Balanced mode construction](paper_localization_balanced_modes_2026_09_17.md):
  exact full-sign bridges with controlled scalar responses on an intrinsic
  low-effective-mode code, plus an all-order near-level-cover extension.
- [Symmetric frame universality](paper_symmetric_frame_universality_2026_09_17.md):
  common-Gibbs fourth-moment comparison, strengthened by explicit querywise
  truncation to an unconditional O(W4^(1/4)log(K)^(3/4)) all-offset bound.
- [Anisotropic flat-mode entropy](paper_localization_anisotropic_entropy_2026_09_17.md):
  a paid coarse mode-energy revelation and Gaussian noise with constant
  physical diagonal give the entropy envelope Psi(theta)~theta log(1/theta).
- [General sign-frame information](paper_director_flat_diagonal_information_2026_09_17.md):
  the director's extension to any isotropic sign-mode law, combined with
  the support-free deterministic realization in Bernoulli Section 22.
- [Critical support-free realization](paper_localization_critical_realization_2026_09_17.md):
  nuclear-norm empirical control simultaneously for every response matrix,
  sharp aggregate covariance, and an anchored center-cover extension for
  arbitrary isotropic sign laws with k=O(sqrt(n)).
- [All-law Hadamard center obstruction](paper_localization_all_law_hadamard_obstruction_2026_09_17.md):
  an explicit quadratic absolute-overlap dual forces Hamming radius at
  least 1/72 from every vanishing-response center code under every
  isotropic physical sign law, on the actual near-half Hadamard family.
- [External-field regularization](paper_localization_external_field_regularization_2026_09_17.md):
  Gaussian near-level width control, its all-pattern exact-sign transfer,
  and I--MMSE entropy conversion; independently audited ordinary random-star
  resampling gives same-order microscopic regularity at o(n^(3/2)) cap cost.
- [All-energy Gaussian-edge stability](paper_localization_gaussian_edge_stability_2026_09_17.md):
  actual-sign preparation yields one deterministic subexponential witness
  for every ground state in a continuum of independent Gaussian edge
  perturbations. Deficit shells exclude all outside states, and the same
  witness approximates the full low-temperature partition function.
  The perturbed neighbors are weighted matrices, not new full signings.
- [Actual-sign noise stability](paper_localization_actual_sign_noise_stability_2026_09_17.md):
  cloned preparation plus all-energy positive-part replacement and
  Bernoulli martingale concentration give a deterministic subexponential
  witness for random full-sign edge-flip neighbors, with exponential
  escape probability. The independently audited
  [mixed-tail improvement](paper_director_mixed_tail_sign_stability_2026_09_17.md)
  widens the shrinking-window regime by paying Gaussian and jump metrics
  separately. Neither theorem is a uniform claim over all flip patterns.
- [Gaussian-sign scalar response](paper_localization_gaussian_sign_response_2026_09_17.md):
  exact isotropic two-polarity covariance mixtures retain a nonlinear
  absolute-response discount, but their sharp near-half floor
  0.750855124 exceeds the needed 3/4 slope. The uniform singular-endpoint
  scalar comparison and the full all-state parent certificate are paid.
- [Nonlocal sign response](paper_localization_nonlocal_sign_response_2026_09_17.md):
  explicit corrected Hadamard eigenlaws realize near-unit radial
  covariance, escape the 2/pi Gaussian-angle loss, and have full
  eta-nearcode response at most 1/sqrt(pi)+4sqrt(eta)+O(n^(-1/8)).
  Matching columns have uniform all-query exponential moments; scheduled
  identity repair restores aggregate isotropy at O(n^(5/4)) full-bridge
  cost. A finite parent certificate still retains every entropy/escape
  cost, so this is not an unconditional extension recurrence.
- [Matching increments and full-parent escape](paper_localization_matching_escape_2026_09_17.md):
  exact centered-absolute mixed tails yield an all-old-word chained
  certificate for the nonlocal bridge. An explicit rare-jump pair of
  actual eigenwords rules out a dimension-free subGaussian Hamming
  increment bound. No balanced-increment assumption is inserted.
- [Gaussian--Poisson matching comparison](paper_director_matching_gaussian_poisson_2026_09_17.md):
  the director's light-coefficient/Palm construction, with this track's
  weighted Stein sensitivity refinement, gives Wasserstein error
  O_D(p^(-1/4)) to a same-variance Gaussian plus compound-Poisson law.
  Its absolute response is below the Gaussian value even when jumps
  persist; no Gaussian CLT is claimed.
- [Adaptive actual-energy response](paper_director_adaptive_energy_response_2026_09_17.md):
  the director's high-covariance projector / low-covariance signed-energy
  truncation dichotomy and minimax give ONE exactly isotropic, uniformly
  subGaussian physical sign law with a fixed absolute-response discount
  on every |H_A|>=c n^(3/2) word, for any actual full signing with
  ||A||op<=L sqrt(n). This track independently reconstructed the full
  proof and verified that the archived scalar Gaussian-sign comparison
  applies without flat off-diagonals. The spectral hypothesis remains
  substantive; the explicit discount is not a convergence slope.
- [Operator-free low-cap response](paper_director_low_cap_uniform_response_2026_09_17.md):
  the thin signed-energy alternative, Grothendieck coordinate localization,
  and the inherited complementary half-range budget remove that spectral
  hypothesis for sufficiently low-cap actual signings. This track's
  [explicit quantitative audit](paper_localization_low_cap_response_quantitative_2026_09_17.md)
  gives a sixth-power margin-to-discount formula. One isotropic uniformly
  subGaussian physical law protects the entire macroscopic energy code.
- [Its full-parent deployment cost](paper_localization_low_cap_parent_2026_09_17.md):
  a global affine response envelope and centered absolute-increment
  chaining give actual-sign iid-parent certificates with every old word,
  new spin, and child cap paid. The current Gaussian-pair class has a
  sharp response floor0.786393873897, above the required parent slope;
  the cloned preparation also retains a non-negligible extension cost.
- [Exact quartic physical laws](paper_localization_quartic_isotropic_laws_2026_09_17.md):
  pure degree-four densities are exactly isotropic. A fixed-rank feature
  law gives a fixed scalar-response discount at information cost below
  log3, while a true marked-energy quartic changes energy variance but
  not leading responses. An isotropic query dual obstructs the entire
  convex mixture class of these quartic feature laws.
- [Smooth finite energy witnesses](paper_localization_smooth_energy_witnesses_2026_09_17.md):
  any bounded-Lipschitz likelihood of a fixed finite list of normalized
  cap-bounded actual quadratic energies leaves all leading Boolean
  responses at kappa. Rare, nonsmooth, and growing-complexity likelihoods
  are explicitly outside the proved scope.
- [Growing-rank cold feature tilts](paper_localization_growing_rank_cold_tilt_2026_09_17.md):
  a direct Fourier proof, strengthened by a positive convex extension,
  gives operator covariance error O(r/n), Frobenius error O(r^(3/2)/n),
  all-query response error O((r^2/n)^(1/4)), and Gaussian variance-change
  information cost with o(1) error for rank r=o(sqrt(n)). A separately
  audited external-field argument yields actual cold proxy1+o(1).
- [Sharp block-code information price](paper_director_block_code_response_information_2026_09_17.md):
  the director's entropy lower bound and this track's exact shared-phase
  construction (AppendixA) show cost r log(1/epsilon)+O(r) for all
  2^r block-constant queries, with exact isotropy, full physical support,
  and proxy O(epsilon^(-2)). The finite upper law has NO rank restriction
  when block size b>=4096 epsilon^(-2).
- [Sharp block-code subGaussian tradeoff](paper_localization_block_code_subgaussian_tradeoff_2026_09_17.md):
  exact isotropy and fixed proxy K force leading large-rank response at
  least kappa/sqrt(K). An actual common-phase slice law attains this
  bound for K>=2, and the director's physical Gibbs construction extends
  sharpness to K>1. The
  [joint information/tail frontier](paper_director_block_response_information_tail_frontier_2026_09_17.md)
  has leading per-feature price (1-1/K)log(1/delta) for approaching this
  floor within delta. This track independently audited the full theorem.

These are proved finite tools, unconditional structural results, and
conditional extension certificates. The low-cap theorem now DOES give
a fixed uniform response discount on macroscopic high-energy codes of
actual nearoptimal signings. None establishes a response below the needed
parent slope, a favorable full-parent recurrence, a new cap constant,
or convergence of M_n/n^(3/2).
Moreover the finite whole-nearlevel entropy-slope shortcut is now ruled
out universally by a fixed-radius flip count. The center-cover versions,
which retain and pay correlated deviations, are not subject to that
particular obstruction. The entropy theorems remain valid as finite tools.
