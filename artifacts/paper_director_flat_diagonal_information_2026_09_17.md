# Small absolute response forces low entropy: a sign-frame information theorem

2026-09-17. **PROVED and independently reconstructed by both the discrepancy
and localization researchers.** Derived jointly from the Gaussian-observation paper mechanism,
balanced sign-mode construction, and a weighted matrix Cauchy--Schwarz
inequality. No external novelty claim is made. Classical ingredients and
the optimizer-facing gap are identified explicitly below.

**Applicability correction,18:35UTC:** Sections1--5 remain valid. The
whole-nearlevel profile in Section6 is IMPOSSIBLE, not merely unproved:
every isotropic sign law has maximal response at least
sqrt(2eta/(pi c)) on the full eta-nearlevel code of a sequence with
normalized cap tending to c>0. Hence mu(eta)^2=o(eta), required by
that proposed application, cannot hold. The conditional implication is
retained to preserve the failed derivation; do not count it as an
original-problem reduction. Low-response CENTER codes and correlated
neighborhood control remain separate. See the response-rigidity proof.

## 1. Finite theorem, without an atom-count or orthogonality hypothesis

Let nu be ANY probability law on {+-1}^k satisfying

```math
\mathbb E_\nu hh^T=I_k.
```

A zero mean is not required. View a Boolean word of length n=kp as a
k by p sign matrix X, and define

```math
\Phi_\nu(X)=\frac{\mathbb E_\nu\|h^TX\|_2}{\sqrt n},\qquad
\mathcal C_\mu=\{X\in\{\pm1\}^{k\times p}:\Phi_\nu(X)\le\mu\}.
```

For every zeta,epsilon,delta>0, with d=k(k+1)/2, the following bound holds:

```math
\boxed{\begin{aligned}
\log|\mathcal C_\mu|\le{}&d\log(1+2/\zeta)\\
&+\frac n{2\delta}\left[
 \mu^2+\epsilon\mu+\frac{\zeta(\mu+\epsilon)}{\epsilon}\right]
 +n h\left(\overline\Phi(1/\sqrt\delta)\right).
\end{aligned}}                                                   \tag{1}
```

Here overline-Phi is the standard normal upper tail, not the response
functional, and h is binary entropy in nats. Empty codes are harmless.
The theorem is uniform over the number, location, and probabilities of
nu's atoms; it neither stores all responses nor assumes a product prior.

## 2. Flat-diagonal noise from weighted frame Cauchy--Schwarz

For a positive semidefinite k by k matrix S of trace one put

```math
a(h)=\sqrt{h^TSh},\quad \phi=\mathbb E_\nu a(h),\quad
w(h)=\frac{a(h)+\epsilon}{\phi+\epsilon},\quad
K=\mathbb E_\nu[w(h)hh^T].
```

Because every h_i^2=1 and E w=1,

```math
K_{ii}=1,\qquad K\succeq\frac\epsilon{\phi+\epsilon}I.
```

The block matrix

```math
\begin{pmatrix}
K&I\\ I&\mathbb E_\nu[hh^T/w(h)]
\end{pmatrix}
=\mathbb E_\nu
\binom{\sqrt{w(h)}h}{h/\sqrt{w(h)}}
\binom{\sqrt{w(h)}h}{h/\sqrt{w(h)}}^T
```

is positive semidefinite. Its Schur complement proves
K^(-1)<=E[hh^T/w(h)]. Consequently

```math
\operatorname{tr}(SK^{-1})
\le(\phi+\epsilon)\mathbb E\frac{a(h)^2}{a(h)+\epsilon}
\le\phi(\phi+\epsilon).                                  \tag{2}
```

This is the decisive combination: sign-flat coordinate magnitudes give
EXACTLY constant noise diagonals, while isotropy gives the identity
off-diagonal block and hence the inverse-energy bound. Neither scalar
variance matching nor covariance trace alone gives both properties.

## 3. Proof of the finite entropy theorem

For each X in the code set S_X=XX^T/n. It is positive semidefinite,
has trace one, and phi(S_X)=Phi_nu(X)<=mu. Choose a maximal
zeta-separated collection of S_X in trace norm, with representatives
FROM THE CODE. The unit trace-norm ball in the d-dimensional symmetric
matrix space has a packing bound (1+2/zeta)^d, by disjoint translated
balls and volume scaling. Assign each codeword to a representative
S_L within trace distance zeta. No dimension-dependent operator-net
conversion is used.

Take X uniformly from the code and condition on its label L. Build
K_L from S_L as in Section 2. Let G, conditionally independently of X,
be Gaussian with covariance delta K_L tensor I_p, and observe Y=X+G.
Because ||K_L^(-1)||op<=(mu+epsilon)/epsilon and the conditional
mean of S_X is still within trace distance zeta of S_L,

```math
\mathbb E[\operatorname{tr}(S_XK_L^{-1})\mid L]
\le\mu^2+\epsilon\mu+\zeta(\mu+\epsilon)/\epsilon.          \tag{3}
```

Gaussian-channel mutual information is bounded by the mean KL divergence
to the zero-mean noise law. This yields

```math
I(X;Y\mid L)\le\frac n{2\delta}
\left[\mu^2+\epsilon\mu+\zeta(\mu+\epsilon)/\epsilon\right].
```

This direct KL proof avoids any centering or invertibility ambiguity;
epsilon>0 makes the noise positive definite.

Every noise coordinate has variance EXACTLY delta, regardless of L.
Since X_i is a sign, the explicit decoder sign(Y_i) has error probability
overline-Phi(1/sqrt(delta)), even though the noise coordinates are
correlated. This decoder is not asserted to be MAP. Binary conditional
entropy and coordinatewise subadditivity give

```math
H(X\mid Y,L)\le n h(\overline\Phi(1/\sqrt\delta)).
```

Finally L is a deterministic function of X, so
H(X)=H(L)+I(X;Y|L)+H(X|Y,L). The packing bound and the two preceding
estimates prove (1).

## 4. Asymptotic entropy envelope and sharp declared order of limits

Suppose k=k_n, k_n^2=o(n), with arbitrary isotropic nu_n. For each
FIXED mu>0, first let n tend to infinity in (1), then zeta down to zero,
then epsilon down to zero. For every delta>0 this proves

```math
\limsup_n\frac1n\log|\mathcal C_\mu|
\le \mathcal E(\mu^2):=
\inf_{\delta>0}\left\{
\frac{\mu^2}{2\delta}+h(\overline\Phi(1/\sqrt\delta))\right\}.
                                                               \tag{4}
```

In particular

```math
\mathcal E(\theta)\le(1+o(1))\theta\log(1/\theta)
\quad(\theta\downarrow0).                                \tag{5}
```

Choose 1/(2delta)=log(1/theta). The first term is theta log(1/theta).
The elementary Gaussian tail bound gives tail probability
O(theta/sqrt(log(1/theta))); its binary entropy is
O(theta sqrt(log(1/theta)))=o(theta log(1/theta)). This suffices for
(5); no equality or optimal lower asymptotic is needed here.

If mu_n->0, (4) for every fixed mu eventually containing the code
implies log|C_(mu_n)|=o(n). A uniform finite rate is already provided
by (1); one must not silently discard the label entropy when k^2 is
comparable to n. The earlier Walsh-affine counterexample at p=1 is
outside the vanishing-label-cost regime.

For a prescribed Hadamard mode frame, its energy histogram replaces
the general matrix label and only k=o(n) is needed. That stronger
special-case range is proved separately by the localization researcher.

## 5. Support-free realization and exact-sign bridge consequence

Every declared sequence of isotropic laws nu_n with k_n^2=o(n) admits
deterministic mode labels at q=floor(epsilon n) columns such that

```math
\sum_h q_hhh^T\preceq(q+o(n))I,
```

and the count discrepancy contributes o(n^(3/2)) uniformly to mean
absolute responses. The support-free proof in
[the Bernoulli track, Section 22](paper_bernoulli_2026_09_17.md)
samples q labels and then fixes one list. The empirical covariance has
mean squared Frobenius error (k^2-k)/q. Symmetrization, sign-to-Gaussian
comparison and vector Gaussian comparison give a uniform expected
response error at most 2 sqrt(pi/2) k/sqrt(q), over ALL positive
semidefinite S of trace at most one. Markov bounds yield one list with
covariance error at most 4k/sqrt(q) and response error at most
8 sqrt(pi/2) k/sqrt(q), simultaneously. Independent auditors checked
the full proof. No bounded support or rational atom-weight hypothesis
is needed.

Sampling columns h tensor g, with independent
fair p-signs g, gives a full-sign bridge. Its exact maximization over
new spins has bounded-difference proxy (q+o(n))n and mean at most
q sqrt(n) Phi_nu(X)+o(n^(3/2)). For arbitrary n, use n0=k floor(n/k)
array coordinates and ell=n-n0<k leftovers, filled independently by
fair signs. Their mean cost is at most q sqrt(ell)=o(n^(3/2)); their
entropy cost is ell log2=o(n). Thus this is an all-order realization.

For actual old sign matrices A_n, assume their FULL absolute near-level
codes satisfy Phi_(nu_n)(X)<=mu(eta)+o_n(1), with mu(eta)->0, and

```math
K=\limsup_{\eta\downarrow0}
\frac{\mu(\eta)^2\log(1/\mu(\eta)^2)}{\eta}<\infty.
```

The direct entropy-profile augmentation theorem then applies with
every slope tau>K/2. The required near-level entropy is now PROVED
from the response hypothesis by (4)--(5), not separately assumed.
For a liminf-realizing minimizing sequence this implies K>=3c_*.
All hierarchy levels and error tolerances are fixed before n tends
to infinity, as in the earlier direct criterion.

This removes BOTH a separate entropy hypothesis and a finite-support
realization hypothesis from the conditional construction. It does
not prove favorable response control for actual minimizing children.
The sampled list has at most q modes; its existence is proved uniformly
over the original law, not assumed from a numerical quadrature.

## 6. Classical inputs and what the combination adds

The matrix Schur complement, trace-norm packing, Gaussian-channel KL
bound, and binary entropy bound are classical. The reusable statement
here is the sign-frame response-to-entropy implication with its exact
finite complexity cost and its actual-sign augmentation consequence.
It combines the localization and rounding tracks without replacing
the original word distribution by a product law. It applies to any
declared Boolean code of block responses, not only energy maximizers.
An external novelty search is required before claiming priority.
