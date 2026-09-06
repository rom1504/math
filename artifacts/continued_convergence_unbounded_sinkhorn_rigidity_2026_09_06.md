# Strict Gaussian self-transport concavity and exact zero-drift rigidity

Date: 2026-09-06. The extension and equality argument below were developed
jointly with the director. They prove the zero-drift hypothesis (B) in
`continued_convergence_terminal_gap_reduction_2026_09_06.md`. The separate
temperature-alignment hypothesis `B E<=E` remains open. No original minimax
convergence claim follows from this artifact alone.

## 1. The compact theorem and the unbounded issue

For a probability law `alpha` on `R^d` with finite second moment, define

`F_t(alpha)=inf_(gamma in Pi(alpha,alpha))`
` { D(gamma || alpha tensor alpha)+t E_gamma |X-Y|^2 }`,
`Phi_t(alpha)=-F_t(alpha)/2`, `k(x,y)=exp(-t|x-y|^2)`, `t>0`.

Feydy et al., *Interpolating between Optimal Transport and MMD using
Sinkhorn Divergences* (AISTATS 2019), Proposition 3 and Proposition 4,
give a variational representation and strict convexity of `Phi_t` on
compact spaces with a positive universal kernel. Their supplement,
Sections B.3--B.4, contains the measure-variable proof. We use their
compact result with Gaussian kernel and extend it below; the published
compact statement is not being asserted for unbounded support without
proof. [Primary paper](https://proceedings.mlr.press/v89/feydy19a/feydy19a.pdf),
[primary supplement](https://proceedings.mlr.press/v89/feydy19a/feydy19a-supp.pdf).

The auxiliary measure can have infinite total mass on an unbounded space.
The correct domain is positive Radon measures of finite Gaussian energy,
not merely finite measures. This matters already for some heavy-tailed
two-dimensional sources with finite second moment.

## 2. Constructing the scaling measure by compact truncation

Let `alpha_R` be `alpha` conditioned on a ball of radius `R`, along radii
with positive conditioning mass tending to one. The compact theorem
provides a positive measure `mu_R`, mutually absolutely continuous with
`alpha_R`, satisfying

`alpha_R=(k mu_R) mu_R`, `||mu_R||_k^2=1`.                   (1)

Here `(k mu)(x)=integral k(x,y)dmu(y)` and the squared energy is the
double kernel integral. For every unit cube `Q`, positivity of the kernel
and its lower bound on `Q tensor Q` give

`mu_R(Q)<=exp(td/2)`.                                      (2)

Thus a subsequence converges vaguely to a positive Radon measure `mu`,
with the same uniform local-mass bound. The Gaussian tails together with
(2) imply that `k mu_R` converges to `k mu` locally uniformly: partition
space into unit cubes, discard the uniformly summable Gaussian tails,
and use vague convergence on the remaining compact region. The same
argument applies to spatial derivatives if desired, although derivatives
are not needed here.

Passing (1) against compactly supported continuous test functions yields

`alpha=(k mu) mu`.                                         (3)

In particular `mu` is nonzero. Integrating (3) and using Tonelli gives
`||mu||_k^2=1`. Gaussian-kernel Cauchy--Schwarz, first for compact
restrictions of `mu` and then by monotone convergence, gives
`0<(k mu)(x)<=1` for every `x`. Positivity follows also by retaining any
compact ball with positive `mu` mass. Consequently

`f(x)=-log(k mu)(x)>=0`,
`f(x)<=t(|x|+R_0)^2-log mu(B_(R_0))`                       (4)

for a fixed such ball. Thus `f` is integrable against every source
`alpha` under consideration. Equation (3) says `mu=exp(f) alpha`.

Define the probability coupling

`gamma(dx,dy)=k(x,y) mu(dx) mu(dy)`.

Its marginals are `alpha`, and

`log[dgamma/d(alpha tensor alpha)] = f(x)+f(y)-t|x-y|^2`.

All displayed terms are integrable: (4) controls the potentials and the
finite second moments control the quadratic cost. For any competing
self-coupling `theta` with finite relative entropy, the chain rule gives

`D(theta||alpha tensor alpha)+t E_theta|X-Y|^2`
` =2 integral f dalpha + D(theta||gamma)`.

Couplings of infinite relative entropy cannot improve the infimum.
Therefore

`F_t(alpha)=2 integral f dalpha`,                           (5)

and `gamma` is the unique optimal self-coupling. This establishes the
required unbounded scaling representation directly, without assuming
existence or integrability of unbounded Schrödinger potentials.

## 3. Strict convexity on all finite-second-moment laws

The Gaussian energy is strictly positive on nonzero differences of
positive finite-energy Radon measures satisfying (2). To see this, its
quadratic form is a positive constant times

`||exp(-2t|.|^2) * (mu_0-mu_1)||_(L^2)^2`.

The local-mass bound makes the difference a tempered distribution.
Convolution by a Gaussian is injective on such distributions: after
Fourier transformation its multiplier never vanishes; testing against
compactly supported smooth functions permits division by that multiplier.
The convolution-square identity follows first by Tonelli for positive
measures and then by polarization; finite energies control the cross term.

Let `alpha_0 != alpha_1`, let their scaling measures from Section 2 be
`mu_0,mu_1`, and take `0<q<1`. Put
`alpha_q=(1-q)alpha_0+q alpha_1` and `mu_q=(1-q)mu_0+q mu_1`.
The measures are mutually absolutely continuous in corresponding pairs.
Write `I(alpha||mu)=integral log(dalpha/dmu)dalpha`, without the mass
correction used for generalized KL divergence. This expression is finite
for these pairs and their mixtures. Indeed `mu_i>=alpha_i`, while (4)
bounds `log(dmu_i/dalpha_i)` above by a quadratic; the mixture density
inherits the maximum of the two quadratic bounds.

Joint convexity of this relative-entropy expression and the energy norm
give

`I(alpha_q||mu_q)+(1/2)||mu_q||_k^2`
` <= (1-q)[I(alpha_0||mu_0)+1/2]`
`    +q[I(alpha_1||mu_1)+1/2]`,                             (6)

with strict inequality if `mu_0!=mu_1`. If instead the measures `mu_i`
coincide, strict convexity in the first argument gives strict inequality
because the source laws differ. Joint convexity is valid also for these
possibly infinite measures: it follows pointwise from the log-sum
inequality. Alternatively replace each auxiliary measure `mu` by the
finite measure `exp(-|x|^2)mu`, retaining the source `alpha`; the change
in `I(alpha||mu)` is the linear source term `-integral |x|^2 dalpha`.

For completeness, the needed variational upper bound for `Phi` requires
no unbounded dual-attainment theorem. Put `h=log(dmu_q/dalpha_q)`. For any
self-coupling of `alpha_q`, compare its entropy to the finite measure
`k mu_q tensor mu_q`. The usual entropy variational inequality gives

`Phi_t(alpha_q) <= I(alpha_q||mu_q)+(1/2)||mu_q||_k^2-1/2`.

One can use the stronger term `(1/2)log ||mu_q||_k^2` instead; the
displayed form follows from `log z<=z-1`. All potential integrals are
finite by the preceding quadratic bounds. For `i=0,1`, (5) makes this
bound an equality with `||mu_i||_k^2=1`. Combining it with the strict
version of (6) proves:

**Strict-concavity theorem.** For every `t>0`, `F_t` is strictly concave
on the convex set of probability laws on `R^d` with finite second
moment. Equivalently `Phi_t` is strictly convex there.

## 4. Strict tensorization equality

The already established tensorization inequality extends by `W_2`
continuity to these laws:

`F_t(law(X,Y)) <= F_t(law X)+F_t(law Y)`.                  (7)

Equality holds for the product law. This also follows directly by taking
a product of optimal self-couplings and applying the information chain
rule in the reverse direction. Thus the product law is a global maximum
of `F_t` on the convex set of joint laws with the specified marginals.

Strict concavity now gives the exact equality case in (7): if a distinct
joint law also attained this maximum, their midpoint would strictly
exceed it while retaining the same marginals. Hence

`F_t(law(X,Y))=F_t(law X)+F_t(law Y)`
` if and only if X and Y are independent`.                (8)

No extension of the differentiability theorem or product-potential
gradient is required for this argument.

## 5. Exactly Gaussian laws have zero one-step Bellman drift

Let `nu` be symmetric on `R`, with finite second moment. The safe global
reversal and input-swap averaging permits a Bellman pair law `pi` whose
two signed marginals are both `nu`. Write `I_pi=I(A;B)`. For any
self-coupling `((A,B),(A',B'))` of `pi`, the information identity is

`I(AB;A'B')-I(A;A')-I(B;B')`
` = I(AA';BB')-2 I_pi >= -I_pi`.

It follows after adding the separated quadratic costs and taking the
infimum that

`F_t(pi)>=2F_t(nu)-I_pi`.                                  (9)

The identity uses finite `I_pi`; policies of infinite `I_pi` have value
minus infinity and need not be considered. Combining (9) with (7) on
the orthogonal outputs `U=(A+B)/sqrt(2)`, `V=(A-B)/sqrt(2)` yields the
quantitative policy bound

`Phi_t(nu)-{[Phi_t(nu_+)+Phi_t(nu_-)]/2-I_pi/2}`
` >= I_pi/4`.                                             (10)

The compactness proof in the terminal-gap artifact shows that `B Phi`
attains its maximum. Therefore `Gamma(nu)=Phi_t(nu)-B Phi_t(nu)=0`
forces the maximizing pair to be iid. Its children have the same law
`nu'=law((A+B)/sqrt(2))`. Zero drift then says `F_t(nu')=F_t(nu)`.
Orthogonal invariance and product additivity give

`F_t(law(U,V))=2F_t(nu)=2F_t(nu')`.

By (8), `U` and `V` are independent. An elementary characteristic-function
argument identifies the source. If `var(nu)=v` and its real characteristic
function is `phi`, independence gives

`phi(a+b)phi(a-b)=phi(a)^2 phi(b)^2`.

Taking `a=b` gives `phi(2a)=phi(a)^4`. Hence for every real `z`,

`phi(z)=phi(z/2^j)^(4^j) -> exp(-v z^2/2)`

by the finite-variance expansion at zero. Thus `nu=N(0,v)`, including
`delta_0` when `v=0`. Conversely a centered Gaussian is unchanged by iid
pairing and has zero drift. We have proved the exact classification

`Gamma(nu)=0 if and only if nu is a centered Gaussian`.    (11)

On every such law the latent envelope equals the terminal potential:
the trivial label gives the Gaussian value, and the general lower-envelope
inequality gives the reverse comparison. In particular

`Gamma(nu)=0 => Phi_t(nu)=E_t(nu)`.

This discharges hypothesis (B) of the conditional terminal-gap theorem.
The remaining hypothesis is solely the temperature-alignment inequality
`B E_t<=E_t`; it is not proved by the strict-concavity argument above.
