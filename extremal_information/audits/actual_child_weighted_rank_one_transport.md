# Weighted rank-one transport exposes the posterior-retuning charge

Status: **rigorous finite theorem and exact scope no-go**.  The support-
cardinality step in Theorem 37.52 admits a weighted version for the actual
child prior.  The correct replacement for `log |supp mu|` is not the prior
Shannon entropy: it is an interpolation between prior Renyi entropy and the
prior surprisal measured under the negative-path-retuned posterior.

Consequently, a genuine effective-entropy improvement requires control of
that posterior retuning.  This is precisely the child-factor resource
isolated in Theorem 37.55, not a consequence of positive-temperature full
support or scalar child entropy.  The theorem below neither proves that the
cardinality bound is asymptotically sharp for optimizing children nor
produces a Level-6 recurrence.

## 1. Exact setup

Let `mu` be any strictly positive law on a finite set
`mathcal Q subseteq {+-1}^d`.  In the application, `mu` is the exact
orientation-conditioned actual-child prior on

```math
 Q_{ij}=SX_iY_j,
 \qquad d=mn,
 \qquad |\mathcal Q|=2^{m+n-1}.                   \tag{WT.1}
```

Let `U` be the fair law on `B in {+-1}^d`, and define the binary-channel
output density and its ordinary forward posterior by

```math
 \begin{aligned}
 p(B)&={\mathbb E_\mu e^{t\langle B,Q\rangle}
              \over(\cosh t)^d},\\
 \mu_B(Q)&={\mu(Q)e^{t\langle B,Q\rangle}
              \over\mathbb E_\mu e^{t\langle B,Q\rangle}},\\
 m(B)&=\mathbb E_{\mu_B}Q.
 \end{aligned}                                     \tag{WT.2}
```

For an arbitrary bridge law `q`, put

```math
 \bar\mu_q=\mathbb E_{B\sim q}\mu_B,
 \qquad
 \mathcal C_q
 =\mathbb E_{\bar\mu_q}[-\log\mu(Q)].             \tag{WT.3}
```

The result therefore applies to every actual negative disorder law
`q=q_s proportional p^sU`, `-lambda<=s<=0`.

## 2. Weighted-max transport

For `0<=alpha<=1`, let

```math
 H_\alpha(\mu)
 ={1\over1-\alpha}\log\sum_Q\mu(Q)^\alpha          \tag{WT.4}
```

with the usual continuous interpretations at `alpha=0,1`.

**Theorem WT.1 (weighted rank-one posterior transport).**  For every
`0<=alpha<=1` and every `q` absolutely continuous with respect to `U`,

```math
 \boxed{
 \mathbb E_q\langle B,m(B)\rangle
 \le
 \sqrt{2d\left\{
 D(q\Vert U)+(1-\alpha)H_\alpha(\mu)
 +\alpha\mathcal C_q\right\}}.}                  \tag{WT.5}
```

At the two endpoints this reads

```math
 \begin{aligned}
 \alpha=0:&\quad
 E_q\langle B,m(B)\rangle
 \le\sqrt{2d\{D(q\Vert U)+\log|\mathcal Q|\}},\\
 \alpha=1:&\quad
 E_q\langle B,m(B)\rangle
 \le\sqrt{2d\{D(q\Vert U)+\mathcal C_q\}}.
 \end{aligned}                                    \tag{WT.6}
```

Thus the first line is exactly the support transport used in Theorem 37.52.

*Proof.*  Fix `theta>0` and define the weighted maximum

```math
 X_{\alpha,\theta}(B)
 =\max_{Q\in\mathcal Q}\left\{
 \langle B,Q\rangle+{\alpha\over\theta}\log\mu(Q)
 \right\}.                                        \tag{WT.7}
```

Exponentiating the maximum and then summing gives

```math
 \begin{aligned}
 \mathbb E_Ue^{\theta X_{\alpha,\theta}}
 &\le\sum_Q\mu(Q)^\alpha
       \mathbb E_Ue^{\theta\langle B,Q\rangle}\\
 &=(\cosh\theta)^d\sum_Q\mu(Q)^\alpha
 \le e^{d\theta^2/2}\sum_Q\mu(Q)^\alpha.        \tag{WT.8}
 \end{aligned}
```

For each `Q`, (WT.7) also gives

```math
 \langle B,Q\rangle
 \le X_{\alpha,\theta}(B)
    +{\alpha\over\theta}[-\log\mu(Q)].            \tag{WT.9}
```

Average first under `mu_B`, then under `q`.  Donsker--Varadhan transport
and (WT.8) yield

```math
 \begin{aligned}
 E_q\langle B,m(B)\rangle
 &\le {1\over\theta}\left\{
 D(q\Vert U)+\log\sum_Q\mu(Q)^\alpha
 +\alpha\mathcal C_q\right\}+{d\theta\over2}.    \tag{WT.10}
 \end{aligned}
```

For `0<=alpha<=1` the expression in braces is nonnegative.  Optimize in
`theta` and use
`log sum_Q mu(Q)^alpha=(1-alpha)H_alpha(mu)` to obtain (WT.5). `square`

## 3. The missing term is exactly posterior retuning

The cross-entropy in (WT.3) has the exact decomposition

```math
 \boxed{
 \mathcal C_q
 =H(\bar\mu_q)+D(\bar\mu_q\Vert\mu).}              \tag{WT.11}
```

For any `0<gamma<1`, entropy transport on the latent alphabet gives

```math
 \boxed{
 \mathcal C_q
 \le H_{1-\gamma}(\mu)
    +{1\over\gamma}D(\bar\mu_q\Vert\mu).}         \tag{WT.12}
```

Indeed, apply the entropy variational inequality to
`f(Q)=-log mu(Q)`:

```math
 \begin{aligned}
 \gamma E_{\bar\mu_q}f
 &\le D(\bar\mu_q\Vert\mu)
   +\log E_\mu e^{\gamma f}\\
 &=D(\bar\mu_q\Vert\mu)
   +\log\sum_Q\mu(Q)^{1-\gamma},                 \tag{WT.13}
 \end{aligned}
```

and divide by `gamma`.

There is also an exact, but generally nonclosing, bridge-level bound.  Let
`Pi` be the forward output law with density `p` from (WT.2).  The posterior
kernel `B mapsto mu_B` sends `Pi` to `mu` by Bayes calibration and sends
`q` to `bar mu_q`.  Data processing therefore gives

```math
 \boxed{
 D(\bar\mu_q\Vert\mu)\le D(q\Vert\Pi).}           \tag{WT.14}
```

Combining (WT.5) and (WT.12) gives the valid prior-Renyi form

```math
 \boxed{
 \begin{aligned}
 E_q\langle B,m(B)\rangle
 \le\sqrt{2d\biggl\{
 &D(q\Vert U)+(1-\alpha)H_\alpha(\mu)\\
 &+\alpha H_{1-\gamma}(\mu)
 +{\alpha\over\gamma}D(\bar\mu_q\Vert\mu)
 \biggr\}}.
 \end{aligned}}                                    \tag{WT.15}
```

Replacing the last divergence by (WT.14) is legal but does not make the
negative path simpler: `q_s` and `Pi=q_1` lie on opposite sides of the
exponential family, and no sublinear bound on `D(q_s||Pi)` is known.

For the exact sector-factorized child law, Theorem 37.55 further splits
`D(bar mu_q||mu)` into sector retuning, the two child-factor KL retunings,
and induced cross-child dependence.  Therefore the additional term in
(WT.15) is exactly the missing directional resource, not a new scalar
consequence of child entropy.

## 4. Refined overlap floor

Use the notation of Theorem 37.52:

```math
 S_s=E_{q_s}\sum_e r_e^2,
 \qquad
 A_\rho={1\over1+\rho},
 \qquad
 C_{\rho,\delta}=A_\rho+{\delta\over1-\rho^2}.
```

The pointwise posterior-energy argument there proves

```math
 {S_s\over d}\ge {1\over C_{\rho,\delta}}
 \left\{A_\rho-{E_{q_s}\langle B,m(B)\rangle\over\rho d}\right\}
 \qquad(-\delta\le s\le0).                         \tag{WT.16}
```

Theorem WT.1 therefore gives, for every `0<=alpha<=1`,

```math
 \boxed{
 {S_s\over d}\ge {1\over C_{\rho,\delta}}
 \left[
 A_\rho-\sqrt{{2\mathcal E_{\alpha,s}\over\rho^2d}}
 \right]_+,}                                      \tag{WT.17}
```

where

```math
 \mathcal E_{\alpha,s}
 =D(q_s\Vert U)+(1-\alpha)H_\alpha(\mu)
  +\alpha\mathcal C_{q_s}.                         \tag{WT.18}
```

At `alpha=0`, (WT.17) is exactly the finite support-cardinality theorem.
An effective-entropy improvement is real only if (WT.18) is smaller; prior
entropy alone does not establish this because `mathcal C_(q_s)` is measured
under the retuned posterior.

## 5. Tiny weights and fair-base transitivity

Two exact observations delimit possible improvements.

First, if a latent set `R` satisfies

```math
 \mu(R)\le e^{-cN},
 \qquad
 \bar\mu_q(R)\ge\eta>0,                            \tag{WT.19}
```

then every `Q in R` has surprisal at least `cN`, and hence

```math
 \boxed{
 \mathcal C_q\ge\eta cN,
 \qquad
 D(\bar\mu_q\Vert\mu)\ge\eta cN-h(\eta).}       \tag{WT.20}
```

The second inequality is binary KL data processing.  Thus fixed posterior
mass on exponentially rare child states either restores a linear weighted-
max charge or is itself a linear posterior-retuning certificate.  Positive
finite temperature gives full support but does not prohibit exponentially
unequal Gibbs weights.

Second, the unweighted maximum used in Theorem 37.52 has an exact fair-base
symmetry.  Let

```math
 Q^*(B)\sim\operatorname {Unif}
 \arg\max_{Q\in\mathcal Q}\langle B,Q\rangle.     \tag{WT.21}
```

At positive temperature `mathcal Q` is the complete rank-one sign orbit.
Simultaneous row and column sign switches act transitively on `mathcal Q`,
preserve `U`, and carry the argmax set equivariantly to the switched argmax
set.  The law of `Q^*(B)` is therefore invariant under a transitive group
and hence uniform on `mathcal Q`.  Consequently,

```math
 \boxed{
 E_U[-\log\mu(Q^*)]
 ={1\over|\mathcal Q|}\sum_Q[-\log\mu(Q)]
 =\log|\mathcal Q|+D(U_\mathcal Q\Vert\mu)
 \ge\log|\mathcal Q|.}                            \tag{WT.22}
```

Thus simply weighting the fair-base maximizing word by its actual Gibbs
prior cannot replace the cardinality charge by Shannon entropy; it pays at
least cardinality on average, and more when the prior is nonuniform.

Equation (WT.22) concerns the unweighted maximum, not the softer posterior
mean.  It does not prove that the original support lower bound is
asymptotically sharp for every optimizing-child sequence.  The exact route
from the maximum to the posterior mean is (WT.5), where the unresolved
quantity is precisely `D(bar mu_(q_s)||mu)`.

## 6. Scope verdict

The weighted theorem gives a genuine finite refinement of the support
argument, but its implication is negative for an entropy-only route:

1. `alpha=0` recovers support cardinality with no child-weight information;
2. every `alpha>0` introduces posterior surprisal;
3. controlling posterior surprisal by prior Renyi entropy requires the
   posterior-retuning divergence in (WT.12);
4. linear mass on rare states makes that divergence linear by (WT.20);
5. the fair unweighted extremizer sees uniform support rather than the Gibbs
   entropy profile by (WT.22).

Therefore exact sector factorization and positive-temperature full support
do not by themselves improve `log |supp mu|` to child Shannon or Renyi
entropy in the negative-path overlap theorem.  The necessary missing
observable is posterior retuning.  Showing it is sublinear would permit an
effective-entropy refinement; showing it is linear would select a coherent
retuning/dependence alternative.  Neither conclusion follows here.
