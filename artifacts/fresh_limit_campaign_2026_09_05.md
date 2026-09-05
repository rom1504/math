# Fresh limit campaign, 2026-09-05

Started at 18:29 UTC from commit 3d05fddeda74e4987608e57257ce564f0a331361.
The user authorized six hours of substantive research on convergence or
nonconvergence of M_n/n^(3/2), with freedom to replace prior architecture.
No active-time claim is inferred solely from elapsed wall time.

## Initial operational check

Local write access verified by creation of this research record. The local
Python environment is Python 3.9.2; git and rg are available. There are eight
CPU cores and approximately 66 GiB free disk. Lean and Lake are not on PATH.
Three independent researchers are assigned algebraic/combinatorial,
variational/structural, and primary-literature attacks. Their initial prompts
contain the original problem and numerical frontier, not the archived route
judgments. Existing untracked files from earlier sessions are preserved.

## First-principles notebook

Only the problem, rigorous interval, small exact values, and exact augmented
cut-code identity were read from ACTIVE_STATE before this analysis. The
original objective remains unsolved. Concrete derivations and their later
archive comparisons will be recorded below.

### Frozen director candidates before historical assessment

1. **A direct universal endpoint inequality.** Seek a universal bound
   Q(A) >= (1/2)n^(3/2)-o(n^(3/2)), possibly via the cut-Laplacian at a
   maximizing spin and its negative-energy counterpart. This is stronger
   than limit existence, but might admit a shorter proof. At a positive
   maximizer, switching makes every cut weight nonnegative; bounding its
   maximum by Q imposes a simultaneous upper cut constraint. The concrete
   first question is whether those two cut constraints and |a_ij|=1 force
   sum_i h_i >= n^(3/2)-o(n^(3/2)). No such implication is asserted.
2. **Scale transfer by nonlocal orthogonal lifting.** For an arbitrary good
   child A, seek sign realizations at larger orders whose Boolean cap is
   controlled by Q(A), rather than its operator norm. A naive A tensor H
   replaces scalar spins by coupled positive/negative eigenspace vectors.
   The exact test is whether this vector relaxation has a vanishing, rather
   than fixed, gap when its inner Hadamard dimension grows. A fixed gap kills
   this version. This is a candidate, not an imported construction theorem.
3. **Direct covering threshold.** Express the optimum as the deficit of
   the augmented cut code, and seek an all-order variational characterization
   of its worst-coset radius. Unlike ordinary radial density, any candidate
   must distinguish an empty covering hole from exponentially small positive
   mass. The concrete question is whether cut-code algebra supplies a
   nonperturbative inequality for the smallest covering multiplicity; an
   entropy limit alone does not answer it.

These are initial hypotheses, not progress claims. The team is separately
developing additional candidates. Historical comparison may now be used to
test the exact mechanisms above.

**Subsequent correction:** the last sentence of candidate 3 was too strong.
An actual ordinary entropy limit *does* imply convergence, using uniform
positive-entropy thickening. The support audit below corrects this initial
assessment as well as the corresponding archived interpretation.

## First substantive findings (campaign continuing)

### Proved: polar-Gram tradeoff

For every hollow signing of order n>=2, with undoubled Boolean cap Q,

```math
Q(A)\ge\frac{n(n-1)}{\pi}
\arcsin\!\left(\frac{n}{\|A\|_*}\right).
```

The associated absolute vector relaxation satisfies
V(A)||A||_* >= n^2(n-1)/2. Its exact minimum n sqrt(n-1)/2 is attained
precisely by conference matrices; a relative epsilon excess forces mean
squared spectral distance from the two-point spectrum at most
2 epsilon/(1+epsilon). The polar-Gram proof was derived by the director
and independently checked by the variational and literature researchers.
It is not an assertion that Boolean minimizers minimize the relaxation.

For Q(A_n)/n^(3/2)<=c+o(1), the Boolean theorem implies
||A_n||_*/(n sqrt(n-1)) >= 1/(pi c)-o(1). This supplies a quantitative
constraint, not a better universal lower constant than the existing bound.
Full proof: [variational artifact](fresh_limit_variational_2026_09_05.md).

### Proved: a tensor-stable minimum has a limit

The algebra researcher constructed explicit symmetric regular Hadamards of
orders 4 and 144. For their tensor semigroup, the same-spin regularization

```math
\mathcal R(B)=\sup_{s=4^a144^b}
\frac{Q(H_s\otimes B)}{s^{3/2}}
```

is genuinely tensor stable. The director's full-sign-seed observation
removes the leading completion payment: allow diagonal entries +/-1 in
the finite seed and erase only the final parent's diagonal. Thus

```math
\limsup_n\frac{M_n}{n^{3/2}}
\le c_{\mathcal R}
=\inf_{B\text{ full symmetric sign}}
\frac{\mathcal R(B)}{(\operatorname{ord}B)^{3/2}},
\qquad
\frac1{\sqrt{2\pi}}\le c_{\mathcal R}\le\frac12.
```

The minimum regularized ratio itself converges to c_R. This is **not** a
lower bound on the original liminf. A finite full seed with a rigorous
regularized ratio below 1/2 would improve the original asymptotic upper
bound; none has been certified. The natural PSD-majorant upper certificate
has an exact half-scale floor for full sign seeds, so it cannot certify
such an improvement. Proof and exact tests are in the
[algebra artifact](fresh_limit_algebra_2026_09_05.md).

One quantifier correction emerged in the director audit: for a *fixed*
seed order, the finite minimum over seeds does commute with the increasing
directed supremum over outer factors. The unproved interchange concerns
unbounded seed orders, or removing the tensor-family restriction.

### Exact/numerical checks and research decision

The independent `computations/audit_fresh_nuclear_2026_09_05.py` audit
enumerates all 2,131,019 root-gauged signings through order eight.
Energy and bilinear-cap calculations are exact integers; eigensolver and
arcsine comparisons are floating-point diagnostics of the analytic proof.
It also verifies the archived order-five cut-Laplacian obstruction using
the exact negative principal determinant -1, without trusting eigenvalues.

The checked coding, Hilbertian, and synchronous-XOR sources do not give an
imported convergence theorem. In particular, independent row/column
Hadamard excess is not same-switch quadratic excess. See the
[source and adversarial audit](fresh_limit_literature_2026_09_05.md).

The original interval is unchanged. The team continues with the actual
regularized Boolean norm and its realizability constraints, while checking
whether the spectral tradeoff supports a quantitatively useful Boolean
argument. Neither auxiliary result is counted as solving convergence.

## Entropy support audit and reopened obligation (19:35 UTC)

The director reconstructed, and the source researcher independently checked,
the following implication. If

```math
S_n(c)=n^{-2}\log\bigl(1+\#\{A:Q(A)\le c n^{3/2}\}\bigr)
```

converges for a dense set of fixed thresholds c, then M_n/n^(3/2)
converges. No formula for S, positivity at its endpoint, or refined entropy
scale is needed. Uniform independent edge noise turns any one signing of
cap c into exp(eta(epsilon)n^2) signings of cap c+epsilon. Consequently an
oscillating support minimum forces an oscillating ordinary entropy.

Also, an empty lower-tail event has exact rate +infinity, whereas a single
signing has rate at most (log 2)/2 at speed n^2. A genuine lower-tail LDP
therefore distinguishes them. Two archived entropy assessments had asserted
otherwise; they are corrected with links to the full proof. This audit
removes a false restriction, not the mathematical task of proving a limit.

The primary-source check found no applicable Bernoulli Boolean lower-tail
LDP. Gaussian upper-tail LDPs have the wrong event and speed. Gaussian
small-operator-norm lower tails cannot be transferred at speed n^2: a
normalized hollow sign matrix has norm at least sqrt(1-1/n), while Gaussian
small-norm events below one still have finite n^2 rate. These are precise
scope distinctions, not an impossibility claim about correlated counting.

Proof: [entropy support audit](fresh_entropy_support_audit_2026_09_05.md).
Independent audit and primary sources:
[entropy/LDP literature](fresh_limit_entropy_ldp_literature_2026_09_05.md).

The campaign continues with actual counting/entropy mechanisms in parallel
with nonquadratic tensor certificates. The original interval is unchanged.

## Universal lower-bound improvement (20:25 UTC)

An analytic theorem, independently checked by the director and two agents,
now proves

```math
\liminf_n M_n/n^{3/2}>0.3396496212118657.
```

For B=A/sqrt(n-1), a bounded-normalized-cap sequence, uniform independent
signs S, and a fixed smooth even h, the rooted field B[S h(BS)] has a
coordinatewise Gaussian limit of variance E h(Z)^2, uniformly in the
coordinate. The low-cap spectral bootstrap, exact endpoint-spin transport,
finite-degree invariance, and four explicit Gaussian contraction cases
prove this without a generic AMP assumption. Paired biased rounding and
a partial best-response update turn it into the improved lower bound.

At the rational parameters t=7/8 and p=8/125 an exact outward-rounded
Fraction calculation bounds the resulting constant between
0.339649621211865756397462117164309583060688819678602946169260 and
0.339649621211865756397462117164309583060688819678602946169325.
The decimal optimization is not needed for the theorem.

Proof: [rooted Gaussian theorem](fresh_limit_rooted_gaussian_2026_09_05.md),
with [transport details](fresh_limit_rooted_response_2026_09_05.md).
Audits: [variational reviewer](fresh_limit_rooted_independent_audit_2026_09_05.md)
and [literature/adversarial reviewer](fresh_rooted_literature_adversarial_audit_2026_09_05.md).
The finite stress test includes random, planted, gauged, regular-Hadamard,
and nonconference tensor families; these are diagnostics, not asymptotic
proof evidence.

The upper bound remains 1/2 and convergence remains open. The next live
theorem concerns the positive oriented energy of a custom second response.
That claim is being tested independently; this checkpoint does not stop the
campaign and does not import its unverified numerical consequences.

## Positive second-response energy and optimized mask (20:30 UTC)

The next theorem has now passed independent audits by the director and two
agents. For G=BS and Y=B[S h(G)], with h fixed smooth even, the exact
asymptotic identity is

```math
\frac1n\mathbb E F(Y)^\mathsf T B[S H(Y)]
\longrightarrow \mathbb E F'(Y_*)\,
\mathbb E[h(G_*)H(Y_*)],
```

where (G*,Y*) is Gaussian with variances 1,E h² and covariance E h.
All indirect paths are summed before applying the actual bilinear cap.
Their dependence on any one seed spin vanishes uniformly by an explicit
cube-gradient/high-moment bound. The direct term is factored using the
proved two-coordinate rooted Gaussian law on all but a vanishing fraction
of pairs. This is not an assertion of generic AMP universality.

At fixed E h=rho, E h²=1, the exact best smooth mask for an outer threshold
alpha is the normalized projection of
k(g)=Pr(|rho g+sqrt(1-rho²)Z|<=alpha) onto constants and its centered part.
It is an actual bounded smooth even function. Cauchy--Schwarz proves this
optimization, and a positive Hermite series with an explicit geometric tail
gives an exact arithmetic certificate. At rho=47/50, alpha=81/100,

```math
\liminf_n M_n/n^{3/2}
\ge 0.385785876908778466066127787300791530025224664775532506435400.
```

The scalar certificate's upper endpoint is 0.385785876909691211.
Proofs: [energy identity](fresh_limit_custom_response_2026_09_05.md),
[smooth-mask optimization](fresh_limit_response_variational_2026_09_05.md).
Audits: [independent counterexample reviewer](fresh_custom_response_independent_audit_2026_09_05.md)
and [independent analytic reviewer](fresh_custom_response_literature_audit_2026_09_05.md).

An exact six-vertex negative-response example falsifies a finite unsmoothed
positivity claim. It does not falsify the theorem, whose dimension limit is
taken with smoothing fixed. Saved tests of the explicit optimized mask on
random, Gram-sign, planted, Hadamard, and tensor families are diagnostics.
Convergence remains open. The next discriminating task is joint two-field
rounding and a proved description of what this rounding family cannot do.

## Hierarchical tree fixed point: a bound above 0.426 (21:15 UTC)

The new analytic result, independently reconstructed by three reviewers, is

```math
\liminf_n M_n/n^{3/2}\ge0.4260903547524242.
```

The director's elementary tree-moment argument supplies the key extension.
For an injectively labeled rooted tree with external-root degree one and
all other degrees odd, every leading moment quotient either contains a
nonempty parity graph, killed by the actual bilinear cut bound, or pairs
whole tree copies. The latter are exactly Gaussian Wick moments. A parallel
energy expansion pairs the physical bridge with one top tree edge; its
branches give normalized Hermite polynomials. Own-spin-free Sobolev
transport extends this identity to bounded smooth Boolean responses.

These child Hermite polynomials form an orthonormal basis of the even
countable Gaussian space. The resulting isometry U maps this space to its
first Gaussian chaos. For an even unit-norm response g with E g'^2<1,
V_(t+1)=U g(V_t) converges in the same Gaussian space to a standard normal
V satisfying V=U g(V). For H=1{|V|<=alpha}, W=U H is jointly Gaussian with
V, with variance p=2Phi(alpha)-1 and covariance w=E g(Z)H(Z). The exact
Boolean certificate is E|W|1{|V|>alpha}. The passage back to sign matrices
uses fixed finite tree/response approximations before the dimension limit.

A completely finite choice of g has Hermite coefficients proportional to
c_r/(97/10+r), for even r=0,...,200, where c_r are the coefficients of
1{|Z|<=37/50}. Exact arithmetic gives derivative energy below
0.997762540769428 and a lower-bound interval

```math
[0.426090354752424224090748824426603635740460994016689775818265,
 0.426090354752424224090748824426603635740460994016689775819420].
```

Proofs: [tree moments](fresh_limit_second_rooted_tree_2026_09_05.md),
[hierarchical energy and fixed point](fresh_limit_hierarchical_tree_energy_2026_09_05.md),
[own-free transport](fresh_ownfree_sobolev_transport_audit_2026_09_05.md).
Independent audits: [tree moments](fresh_tree_chaos_literature_audit_2026_09_05.md),
[energy combinatorics](fresh_tree_energy_independent_audit_2026_09_05.md),
[complete hierarchy](fresh_hierarchical_energy_literature_audit_2026_09_05.md).
The certificate and saved exact output are under `computations/`.

There is also a proved limitation: the paired Gaussian certificate alone
is bounded by the scalar envelope
2 sqrt(p) phi(Phi^{-1}((1+p)/2)), whose numerical supremum is about 0.44955.
An elementary analytic bound already puts this below 1/2. This does not
bound the actual unoriented output energy, nor the original optimum.
The campaign continues with fixed-point improvements and that discarded
energy; no convergence or nonconvergence theorem is claimed.

## Anchored fixed point and midpoint reassessment (21:28 UTC)

The independently audited anchored extension gives the stronger bound

```math
\liminf_n M_n/n^{3/2}\ge0.4283764176561877.
```

The constant Gaussian coordinate U1 is retained as an anchor. Only the
innovation degree is charged in the fixed-point contraction; this is a
strictly larger admissible response class than the previous scalar
construction. An explicit centered bivariate polynomial of total degree
200 has conditional derivative energy below 0.992418. Exact rational
interval arithmetic yields lower endpoint
0.428376417656187756438676313934781469275119835457995857473025.
The director and algebra reviewer independently checked the proof and
certificate, and the director reran it with identical output.

[Complete theorem and certificate](fresh_anchored_fixed_point_independent_2026_09_05.md).
No new matrix-limit theorem is assumed: the existing finite tree-energy
identity realizes the bounded final mask after Gaussian-space approximation.

At the campaign midpoint, this is substantial movement of the original
rigorous interval, but not of an all-order convergence recurrence. The
paired-certificate ceiling remains below 1/2. Work continues in parallel on
a rigorously controlled mask improvement and on new unmarked internal-tree
fields that can see the quadratic energy discarded by pairing. These new
fields are not yet a universal theorem. The original question is open.
