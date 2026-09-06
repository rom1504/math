# Independent audit: unbounded Gaussian self-transport and terminal drift

Date: 2026-09-06. Full independent reconstruction of
`continued_convergence_unbounded_sinkhorn_rigidity_2026_09_06.md`
and `continued_convergence_terminal_gap_reduction_2026_09_06.md`.
Verdict: the strict-concavity and zero-drift theorem passes. Identification
of the deep value with E still requires `B E<=E`; the further Gaussian-
boundary replacement in Section 6 passes UNCONDITIONALLY.
No original minimax convergence conclusion is established.

## 1. Primary compact theorem and its exact scope

The auditor retrieved the primary [Feydy et al. paper, Propositions
3--4](https://proceedings.mlr.press/v89/feydy19a/feydy19a.pdf) and read
the complete [supplement Sections B.3--B.4](https://proceedings.mlr.press/v89/feydy19a/feydy19a-supp.pdf).
Their compact-space result supplies a symmetric scaling measure and a
strictly convex negative half self-cost. The present source correctly
does not cite that theorem as an unbounded-support result. Its Radon
extension and equality analysis are separately reconstructed below.

## 2. No loss at infinity in the scaling construction

Write `k(x,y)=exp(-t|x-y|^2)`. For the compact truncations alpha_R,
the symmetric scaling measure satisfies
`alpha_R=(k mu_R)mu_R` and Gaussian energy one. The kernel is at least
exp(-td) on a unit cube squared, so the cube mass is at most exp(td/2),
uniformly in its location and in R. This gives vague sequential
compactness on all bounded sets. The same cube bound makes the Gaussian
tails uniformly summable, proving local-uniform convergence of k mu_R
along a vaguely converging subsequence. On the compact remainder,
equicontinuity of the kernel family upgrades the pointwise vague limit
to the required uniform convergence.

Testing with compactly supported continuous functions therefore gives
`alpha=(k mu)mu`. This is the crucial step: integration of this exact
identity forces the LIMITING energy to be one, so a lower-semicontinuity
argument has not silently lost energy at infinity. In particular mu is
nonzero. Gaussian-kernel Cauchy--Schwarz on compact restrictions gives
`0<k mu<=1`. A fixed ball with positive mu mass gives the lower Gaussian
bound on k mu and hence a quadratic upper bound on `f=-log(k mu)`.
Thus f is nonnegative and integrable against every source in the theorem.

The coupling `gamma=k mu tensor mu` is a probability with the correct
marginals. Its logarithmic density against alpha tensor alpha is
`f(x)+f(y)-t|x-y|^2`, with all terms integrable. For any competing
coupling of finite entropy the chain rule yields its objective equal
to `2 integral f dalpha + D(theta||gamma)`. Infinite-entropy couplings
cannot improve it. This proves the claimed scaling representation and
unique optimum without an imported unbounded dual-attainment assertion.

## 3. Strict Gaussian energy and the variational inequality

The cube bound also makes each scaling measure translation bounded and
tempered. Gaussian energy is a positive constant times the L2 square of
convolution by exp(-2t|x|^2). Tonelli proves this first for positive
measures; finite energy and polarization control the cross term and the
difference. A zero convolution implies a zero tempered distribution:
the Fourier Gaussian multiplier never vanishes, and division is made
only on the compact support of a smooth test function. No exponentially
growing multiplier is asserted to act globally on tempered distributions.
Thus distinct such Radon measures have strictly positive difference energy.

For mixtures of the two scaling pairs, the auxiliary measure remains
mutually absolutely continuous with its source. Its log density lies
between zero and the maximum of two quadratic bounds, so the uncorrected
relative-entropy integral is finite even if the auxiliary measure has
infinite mass. Joint convexity follows from the pointwise log-sum
inequality. Its combination with strictly convex Gaussian energy gives
the needed strict inequality for the mixture.

The variational direction is correct. For a candidate auxiliary measure
mu of energy z, entropy relative to the finite measure k mu tensor mu
gives

```math
F_t(\alpha)\ge
2\int\log(d\mu/d\alpha)\,d\alpha-\log z.
```

Therefore Phi=-F/2 is bounded ABOVE by
`I(alpha||mu)+(log z)/2`, and then by
`I(alpha||mu)+(z-1)/2`. Equality holds at each constructed scaling
measure, where z=1. This proves strict convexity of Phi, equivalently
strict concavity of F, for ALL finite-second-moment laws in R^d.

## 4. Product equality and the exact zero-drift classification

The already proved self-cost tensorization says F(joint)<=sum of its
marginal self-costs, with equality for the product. The product is thus
a global maximum over the convex set of laws with those marginals.
Strict concavity makes it the UNIQUE maximizer: a second maximizer would
have a midpoint of larger value. This avoids requiring differentiability
or a product-potential gradient on an unbounded space.

For a Bellman policy first take the always-safe reversal and swap
averages, so both signed input marginals equal nu. With finite input
mutual information I_pi, every finite-cost self-coupling of the pair
satisfies

```math
I(AB;A'B')-I(A;A')-I(B;B')
 =I(AA';BB')-2I_\pi\ge-I_\pi.
```

All informations are well-defined here; data processing gives the last
inequality, and finite full information together with finite I_pi
controls the chain-rule components. Infinite-cost couplings and
infinite-information policies are irrelevant to the optimum. Adding
quadratic costs proves `F(pi)>=2F(nu)-I_pi`. Tensorization on the rotated
outputs then gives policy drift at least I_pi/4, with the factor one
quarter correct.

The attained Bellman maximum therefore has zero drift only if its input
pair is iid. Symmetry makes the plus and minus output laws identical.
Zero drift and strict product equality force those outputs to be
independent. Their characteristic functions obey
`phi(a+b)phi(a-b)=phi(a)^2phi(b)^2`, hence `phi(2a)=phi(a)^4`.
Iteration toward zero and the finite-variance expansion identify exactly
`phi(z)=exp(-v z^2/2)`. This includes the zero point mass. Conversely
Gaussian iid pairing has zero drift. Thus zero drift is exactly Gaussian,
and on this set Phi equals the latent envelope E.

## 5. Terminal-gap reduction with its remaining hypothesis

The terminal-gap source was also read completely. Its source-mixture
inequalities pass: concavity supplies the lower bound, and revealing a
shared binary mixture label costs at most its entropy for the upper
bound. Both F and J_lambda are bounded by the logarithmically growing
Gaussian self-cost at the source second moment. Removing a tail of mass
q therefore costs O(q log(1/q)), uniformly on bounded-moment sets and
uniformly for lambda<=t. On a fixed compact interval, W2 source transport
controls the costs, including J uniformly in lambda. Hence Phi and E
are genuinely weakly continuous on bounded-second-moment compact sets,
even when second-moment mass itself can escape to infinity.

Tight policy marginals, continuous child potentials, and lower
semicontinuity of relative entropy make B Phi upper semicontinuous.
The maximum is attained and Gamma=Phi-BPhi is lower semicontinuous.
The zero-drift theorem above now proves a strictly positive drift floor
on each compact set where the terminal gap D=Phi-E is bounded away
from zero.

The stopped-tree proof has the correct budgets. A near-optimal policy's
expected accumulated Gamma is at most the original D plus its
near-optimality error. A surviving branch of length r inside the compact
bad-gap set pays r times the drift floor. Along a uniformly chosen
branch the second moment is a nonnegative martingale, and bounded-time
optional stopping controls the contribution of a high-energy exit by
the original moment times `sup_(s>C) K(s)/s`. Small-gap exits pay at
most epsilon. Under the REMAINING assumption `B E<=E`, the stopped
envelope minus prefix information costs is a supermartingale. These
facts give exactly the source's quantitative terminal-gap bound.

Accordingly the previous zero-drift hypothesis is now discharged, but
temperature alignment is not. Fixed-temperature supersolutions do not
establish it, since different children can select different temperatures.
Without that inequality the conditional Bellman convergence and the
negative restricted-weave conclusion must not be asserted.

## 6. Further unconditional Gaussian-boundary replacement

The newly appended Section 6 of the terminal-gap source was then read
in full. It proves more without assuming temperature alignment, and
this further proof also passes independently.

Set `G(nu)=g_t(m2(nu))`, `h_r=B^rG`, and `H=sup_r h_r`. Iid pairing
shows BG>=G, so h_r increases. For EVERY r,k, h_r<=f_k: when r>=k
apply B^k to h_(r-k)<=Phi; otherwise use h_r<=h_k<=f_k. Thus H is
below the decreasing upper limit. The two-child sum commutes with an
increasing limit, and the policy supremum commutes with the supremum
over r, proving BH=H exactly.

Conditional copying gives H>=E without a finite-alphabet restriction.
The active plus state's second moment is at most 2^r times the original
moment; its weighted Gaussian terminal penalty is therefore O(r2^-r).
All departing minus branches have conditional-copy residual variance V.
The same information telescope then gives H>=g_t(V)-I(X;L), and taking
the latent supremum proves the comparison.

Now use H as the supersolution in the stopped tree, but KEEP the
stopping set defined by D_E=Phi-E. Its compact drift constant is already
proved positive by the Gaussian zero-drift classification. Since H>=E,
the small-gap, high-energy and surviving-terminal excess bounds are
unchanged. A near-optimal upper tree has total drift at most Phi-H plus
its approximation error, bounded by D_E. Therefore the same quantitative
estimate proves

```math
\lim_{r\to\infty}B^r\Phi_t(\nu)
 =\lim_{r\to\infty}B^rG(\nu)
```

for every symmetric finite-second-moment source, UNCONDITIONALLY.
No continuity of H is required: each finite policy tree has finitely
many deterministic source states, and bounded stopping uses only finite
sums. This identifies the decreasing and increasing terminal limits;
it does not make the compact drift constant effective, prove H=E, or
turn the negative E certificate into a negative H upper bound. Those
distinctions remain essential for the original signing problem.
