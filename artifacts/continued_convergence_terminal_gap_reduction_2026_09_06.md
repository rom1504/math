# A conditional terminal-gap theorem for recursive pairing

Date: 2026-09-06. Sections 1--5 prove an implication with two explicit
hypotheses; the companion rigidity proof subsequently discharged (B).
Section 6 proves an unconditional two-sided Gaussian-boundary Bellman
limit. Temperature alignment (A) remains open. This is not a negative
Bellman certificate, a strict minimax upper bound, or a convergence theorem
for the original sequence `M_n/n^(3/2)`.
The purpose is to replace an unspecified uniform deep-CLT obligation by a
local equality-case question plus the already isolated temperature-alignment
question.

Use the Bellman operator and latent envelope from
`continued_convergence_deep_latent_envelope_2026_09_06.md`. Fix `t>0`, and
let the state space be all symmetric probability laws on the real line with
finite second moment. Write

`Phi(nu)=-F_t(nu)/2`, `E(nu)=sup_L[g_t(E Var(X|L))-I(X;L)]`,
`f_r=B^r Phi`, `D=Phi-E`, `Gamma=Phi-B Phi`, `K(s)=-g_t(s)`.

Here `g_t(s)` is the centered Gaussian terminal potential. It satisfies
`K(0)=0`, `K(s)>=0`, and `K(s)=O_t(log(1+s))` at infinity. Known inequalities
are

`E<=f_r<=Phi<=0`, `B Phi<=Phi`, `B E>=E`,
`0<=D(nu)<=K(m_2(nu))`, `Gamma>=0`.

The inequalities initially proved for finite laws extend to finite-second-
moment laws by the continuity argument below. For the lower inequality
`E<=f_r`, it suffices to extend `E<=Phi` and use `B E>=E` and monotonicity.

## 1. Uniform weak continuity despite energy escaping to infinity

For `Q=F_t`, or for any of the source rate-distortion Lagrangians
`Q=J_lambda` with `0<lambda<=t`, a binary source mixture satisfies

`(1-q)Q(mu)+q Q(eta) <= Q((1-q)mu+q eta)`
` <= h(q)+(1-q)Q(mu)+q Q(eta)`.                         (1)

The lower inequality is source concavity. For `F_t`, this is the raw
entropic-self-transport concavity used in the recursive-orbit proof. For
`J_lambda`, it follows also from its representation as an infimum of
linear source expectations over reproduction distributions. The upper
inequality follows by giving the mixture label to both coordinates of
a self-coupling, or by including it in the rate-distortion label. This
charges at most the binary entropy `h(q)`.

Both costs are nonnegative and

`J_lambda(nu)<=F_lambda(nu)<=F_t(nu)<=2K(m_2(nu))`.       (2)

The first comparison uses the second member of a self-coupling as a
reproduction label. The last comparison is the Gaussian upper bound on
the self-transport cost, obtained by appending iid pairing levels and
passing to the ordinary `W_2` central limit theorem.

If `m_2(nu)<=C`, split off the source tail `|X|>R`, of mass `q<=C/R^2`.
Its conditional second moment is at most `C/q`; the complementary
conditional moment is at most `C/(1-q)`. Equations (1)--(2) imply

`|Q(nu)-Q(nu conditioned on |X|<=R)|`
` <= h(q)+2q[K(C/q)+K(C/(1-q))]`.                       (3)

The right side tends to zero uniformly over these laws as `R` tends to
infinity, uniformly also over `0<lambda<=t`. Thus escaping energy is
harmless here even though the second moment itself is not weakly
continuous: the small tail pays only `O(q log(1/q))` information cost.

On a fixed compact source interval, `F_t` is `W_2` continuous. The
Lagrangians `J_lambda` are uniformly `W_2` continuous for `lambda<=t`:
transport the source through a `W_2` coupling and retain the old
reproduction label. Data processing decreases information, while the
change in squared prediction error is at most

`W_2(mu,nu)[2 sqrt(m_2(mu))+W_2(mu,nu)]`,

and use the reverse comparison as well. Consequently (3), with truncation
radii at continuity points of the limiting law, proves:

**Continuity lemma.** On `{nu symmetric: m_2(nu)<=C}`, `Phi` and `E` are
weakly continuous. The latter uses
`E=sup_(0<lambda<=t)[c_t(lambda)-J_lambda]` and the uniformity in `lambda`.

These bounded-moment sets are weakly compact. Indeed they are tight and
are weakly closed by lower semicontinuity of the second moment.

## 2. The one-step drift is lower semicontinuous

Before evaluating a Bellman policy, average it under simultaneous input
reversal and under input swap. This preserves the two child absolute laws
and cannot increase the relative-entropy cost. Both signed input
marginals can therefore be taken to equal the source law `nu`.

For sources `nu_j` with `m_2<=C` converging weakly, almost-maximizing pair
couplings have a weakly convergent subsequence. Its marginals are the
limiting source law. Child second moments are at most `2C`, so their
terminal potentials converge by Section 1. Relative entropy of the pair
with respect to the product of its marginals is jointly lower
semicontinuous. Hence `B Phi` is weakly upper semicontinuous on a
bounded-moment set, and

`Gamma=Phi-B Phi` is weakly lower semicontinuous.             (4)

This also supplies an attained maximizing pair law for `B Phi` at any
fixed source. No finite-support restriction is used in this compactness
step.

## 3. Two precise hypotheses; rigidity subsequently proved

The terminal-gap theorem below assumes:

**(A) Temperature alignment:** `B E<=E` on the whole state space.

**(B) Zero-drift rigidity:** whenever `Gamma(nu)=0`, one has `D(nu)=0`.

Hypothesis (A) is not implied by the fixed-auxiliary-temperature
supersolutions: different children may choose different temperatures.
Hypothesis (B) is a one-step equality-case assertion, not a reformulation
of convergence of the Bellman iterates. The stronger assertion that only
centered Gaussian laws and the zero point mass have `B Phi=Phi` suffices,
since `Phi=E` on those laws. That stronger assertion is now proved in
`continued_convergence_unbounded_sinkhorn_rigidity_2026_09_06.md` and
independently audited; hypothesis (B) is therefore no longer open.

By the continuity and compactness above, (B) implies that for every
`epsilon>0` and finite `C`,

`kappa(epsilon,C)=inf{Gamma(nu): m_2(nu)<=C, D(nu)>=epsilon}>0` (5)

whenever the set is nonempty. An empty set causes no difficulty and may
be assigned `kappa=+infinity`.

## 4. A quantitative stopped-tree bound

**Conditional theorem.** Assume (A) and (B). For a symmetric input `nu_0`
of second moment `E_0`, put `D_0=D(nu_0)`. For every `r>=1`, `epsilon>0`
and `C>E_0`,

`f_r(nu_0)-E(nu_0)`
` <= epsilon + E_0 sup_(s>C) K(s)/s`
`    + K(C) D_0/[r kappa(epsilon,C)]`.                       (6)

In particular `f_r(nu_0)` decreases to `E(nu_0)`.

Proof. Take a depth-`r` policy whose value is at least
`f_r(nu_0)-zeta>=E(nu_0)-zeta`. Follow a uniformly chosen branch of its
binary tree, writing `nu_d` for its random state. Its second moment is a
nonnegative martingale: the average of the two child second moments is
the parent's second moment, by the parallelogram identity. Pair
correlations and energy condensation do not change this fact.

At each node the actual terminal-potential drift is

`Phi(nu)-[Phi(nu_+)+Phi(nu_-)]/2 + I(A;B)/2 >= Gamma(nu)`.

Summing with tree weights telescopes. Consequently the expected sum of
`Gamma` over the full branch is at most `D_0+zeta`.

Stop the branch upon first reaching `D(nu_d)<=epsilon` or `m_2(nu_d)>C`,
and otherwise at depth `r`. Before stopping, every visited node has
`Gamma>=kappa(epsilon,C)`. The probability of surviving to depth `r`
without either exit is therefore at most

`(D_0+zeta)/[r kappa(epsilon,C)]`.                            (7)

At a small-gap stop, any remaining continuation is bounded above by
`Phi<=E+epsilon`. At a high-energy stop with second moment `s`, its
continuation is bounded above by `Phi<=E+K(s)`. Optional stopping for the
nonnegative moment martingale gives expected moment at the stopping
time exactly `E_0`, so the expected high-energy excess is at most

`E_0 sup_(s>C) K(s)/s`.

On surviving terminal branches the excess is at most `K(C)`, and (7)
bounds their probability. Finally (A), iterated only up to this bounded
stopping time, bounds the expected stopped envelope minus the accumulated
prefix information costs by `E(nu_0)`. Combining the three kinds of
stopped branches and sending `zeta` to zero proves (6).

For convergence, first choose `C` large enough that the second term is
small, since `K(s)/s` tends to zero; then choose `epsilon` small and let
`r` tend to infinity. This proves the stated conditional limit without
any uniform CLT assertion for controlled policies.

## 5. Exact remaining obligations and scope

This result provides a state-dependent terminal-majorant correction. It
does not seek a globally small additive correction to any one fixed
temperature branch, so fixed-temperature additive-defect counterexamples
do not apply to it.

The two questions isolated by the initial reduction were separate and local:

1. Can heterogeneous child temperatures violate `B E<=E`?
2. Can a non-Gaussian source have exactly zero one-step terminal drift
   while retaining `Phi>E`? The companion rigidity proof answers no.

A counterexample to (1) remains a legitimate falsifier of the proposed
identification with `E`. A proof of (1), together with the proved rigidity,
would first give
a recursive Bellman limit. Obtaining a negative restricted-weave
certificate would still require an explicit negative value of the
resulting envelope at an allowed `(p,t)`, and an original-problem
convergence theorem would remain a further, different obligation.

## 6. Unconditional Gaussian-boundary replacement after zero-drift rigidity

The companion proof
`continued_convergence_unbounded_sinkhorn_rigidity_2026_09_06.md` proves
that `Gamma=0` holds exactly at centered Gaussian laws. This discharges
hypothesis (B). More strongly, the stopping method gives a deep-limit
identification that does not assume temperature alignment at all.

Define `G(nu)=g_t(m_2(nu))`, and

`h_r=B^r G`, `H(nu)=sup_r h_r(nu)`.

Choosing an iid input pair preserves the variance in both children and
costs zero information, so `B G>=G`. Thus `h_r` increases. Since
`G<=Phi` and `B Phi<=Phi`, it stays below every sufficiently earlier
upper iterate, and in particular `H<=f_infinity`. The finite two-child
sum commutes with an increasing limit, as does the supremum over pair
policies. Consequently

`B H=H`.                                                   (12)

Also `H>=E`. Repeat the conditional-copy construction in the latent
envelope artifact, but terminate discarded branches at `G` instead of
`Phi`. Their variance is exactly `V=E Var(X|L)`. The active branch at
depth `r` has second moment at most `2^r m_2(nu)`, so its weighted
Gaussian value tends to zero because `2^(-r)K(2^r m_2(nu))->0`.
The same telescoping information inequality therefore gives
`H>=g_t(V)-I(X;L)` for every latent channel and hence `H>=E`.
This version does not require a finite source alphabet or an entropy
bound on its growing sums.

Apply the stopped-tree proof with `H` as the supersolution, using (12),
but retain the stopping criterion `Phi-E<=epsilon`. Its compactness
constant `kappa(epsilon,C)` is available because zero drift implies
`Phi=E`. At every stop, `Phi-H<=Phi-E`, so all excess estimates remain
valid. A near-optimal policy has total drift budget at most
`Phi(nu_0)-H(nu_0)<=D_0`. Thus, with no hypothesis (A),

`0<=f_r(nu_0)-H(nu_0)`
` <=epsilon+E_0 sup_(s>C)K(s)/s`
`   +K(C)D_0/[r kappa(epsilon,C)]`.                         (13)

We conclude the unconditional equality

`lim_r B^r Phi_t(nu)=lim_r B^r G(nu)`                       (14)

for every symmetric finite-second-moment source. The left side decreases
and the right side increases, so these are genuine two-sided Bellman
approximations to the same value. There is no asserted effective lower
bound on `kappa`, and no uniform growing-depth Hadamard type theorem is
being imported into this analytic statement.

Equation (14) does not say that the common value equals `E`. Repeated
pairings applied to Gaussian boundary data could still exceed the
single-latent-channel envelope. Proving `B E<=E` would identify the
common value with `E` using Section 4; a strict opposite example would
disprove that simpler formula while leaving (14) intact. No negative
restricted-weave upper certificate is asserted here.

There is nevertheless a genuine fixed-point uniqueness conclusion:
`H` is the unique Bellman fixed point in the order interval `G<=u<=Phi`.
Indeed iteration gives `B^r G<=u<=B^r Phi`, and (14) squeezes both sides.
More generally every terminal potential `psi` satisfying `G<=psi<=Phi`
has `B^r psi->H`, whether or not its own iterates are monotone. These
statements concern only the analytic Bellman problem, not the original
minimax sequence.
