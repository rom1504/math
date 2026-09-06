# Independent primary-source and unbounded-rigidity audit

Date: 2026-09-06. Full fresh read of
`continued_convergence_unbounded_sinkhorn_rigidity_2026_09_06.md` and
`continued_convergence_terminal_gap_reduction_2026_09_06.md`.
Verdict: the unbounded strict-concavity theorem, its tensorization equality
case, and Gaussian zero-drift classification pass independent reconstruction.

## 1. Published input and its exact scope

I checked Proposition 3 and Proposition 4 in the primary paper and read
the complete proofs in supplement Sections B.3--B.4. The published
measure variational formula and strict convexity use a compact space
and a positive universal kernel. They do not themselves license an
unbounded-support statement. The Gaussian kernel restricted to each
compact ball meets their hypotheses. The source artifact uses only
that compact input before supplying its own unbounded extension.
[Primary paper](https://proceedings.mlr.press/v89/feydy19a/feydy19a.pdf),
[primary supplement](https://proceedings.mlr.press/v89/feydy19a/feydy19a-supp.pdf).

The paper's negentropy is the source artifact's Phi, not its positive
self-transport cost F. Thus strict convexity of the former corresponds
to strict concavity of the latter; no sign reversal is missing.

## 2. Unbounded scaling construction

The compact scaling measures have Gaussian energy one. Since a unit
cube has squared diameter d, each cube mass is at most exp(td/2).
This is a uniform LOCAL bound, not an unjustified total-mass bound.
It gives vague compactness, and summable Gaussian tails over cubes
turn vague convergence into locally uniform convergence of k*mu_R.
Passing the scaling equation against compact test functions therefore
gives `alpha=(k*mu)mu` for the limiting positive Radon measure.

Integrating this identity with Tonelli forces the limit's energy to
equal one, so no energy-normalization loss remains. Kernel
Cauchy--Schwarz, first for finite restrictions, gives `k*mu<=1`.
A fixed ball of positive mu mass gives a Gaussian lower bound for
k*mu and hence a quadratic upper bound for `f=-log(k*mu)`.
Consequently f is integrable against alpha under the stated second-
moment assumption. The candidate self-coupling has the exact marginal
alpha and finite entropy and quadratic cost. The relative-entropy
chain rule identifies its cost as `2 integral f d alpha` and proves
optimality without an unbounded dual-attainment assumption.

## 3. Strictness and the variational direction

For two scaling measures, finite Gaussian energy makes their Gaussian
convolutions square integrable. The difference is tempered by the
uniform cube-mass bound. Gaussian convolution is injective: its Fourier
multiplier is nowhere zero, and division on compactly supported smooth
test functions suffices, without treating an exponentially growing
inverse multiplier as a global Schwartz multiplier.

Thus the energy quadratic form is strictly convex on distinct scaling
measures. The relative-entropy term is jointly convex pointwise by the
log-sum inequality even if these positive Radon measures have infinite
total mass. The mixture log-density retains a quadratic integrable bound.
Finally the entropy variational comparison gives

`Phi(alpha_q)<=I(alpha_q||mu_q)+(1/2)||mu_q||_k^2-1/2`.

The direction and the final constant are correct; the stronger logarithmic
energy bound implies it by log z<=z-1. At each actual scaling minimizer
the inequality is equality. Strictness for the mixture proves the claimed
strict concavity of F on all finite-second-moment probability laws.

## 4. Tensorization and zero drift

The product law attains the sum of the marginal self-transport costs.
The already established tensorization inequality says it is a global
maximum over the convex fixed-marginal set. Strict concavity makes it
the UNIQUE such maximizer. This obtains the equality case without
requiring differentiability of unbounded potentials.

For a Bellman input pair pi with finite information, the exact chain
identity and data processing give

`I(AB;A'B')-I(A;A')-I(B;B')`
`=I(AA';BB')-2I_pi >= -I_pi`.

It follows that every policy's terminal drift is at least I_pi/4.
Attainment of the Bellman maximum was proved using the bounded-moment
compactness and continuity in the companion terminal-gap note. Thus
zero drift forces the maximizing pair to be iid. Equality of its input
and output costs then invokes strict tensorization to make sum and
difference independent.

For a symmetric finite-variance source, the resulting characteristic-
function identity implies `phi(2a)=phi(a)^4`. Iterating toward zero and
using the exact quadratic expansion there gives
`phi(z)=exp(-v z^2/2)`. This includes the degenerate variance-zero law.
Conversely a centered Gaussian is unchanged under iid pairing and has
zero drift. The classification is therefore exact.

## 5. Scope after the audit

This discharges the zero-drift rigidity hypothesis in the conditional
terminal-gap theorem. Its stopped-tree proof also reconstructs: the
near-optimal policy pays at most the initial gap in summed drift, the
branch second moment is a nonnegative martingale, and the three stopped
classes have the displayed small-gap, high-energy, and surviving-branch
costs. The separate temperature-alignment inequality `B E<=E` is still
an explicit open hypothesis. No unconditional Bellman upper certificate
or original minimax convergence theorem is claimed by this audit.
