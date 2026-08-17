# Independent audit: rare-event spectral-spike state

**Verdict: PASS for Theorems 1--2 after scope repairs; qualify the state and
composition claims.**

This audit freezes `rare_event_spectral_spike_state.md` at SHA-256
`71935285ab5752431f83c15c47f26f8e3846d0264ec76128f8ff51874cb5ac9d`.
The rank-one limit, its phase-transition constant, and the fixed-depth moment
falsifier are correct.  The finite-rank multiset law is also correct for a
jointly generic, fixed-rank family.  It is not a deterministic congruence for
arbitrary spike continuations, and `(mu,Theta)` is a **presented**
bulk--perturbation state rather than an intrinsic statistic of the total
matrix.  Those distinctions are essential to the benchmark's interpretation.

## 1. Normalization and exact rank-one carrier: PASS

For `x in S^(N-1)`, the Rayleigh--Ritz normalization is exact:

```math
\max_x x^TCx=\lambda_{\max}(C).
```

The rank-one determinant identity has the correct sign,

```math
\det(zI-B-\theta uu^T)
=\det(zI-B)(1-\theta u^T(zI-B)^{-1}u).
```

For a positive spike and Haar `u`, the weighted measure `nu_(N,u)` almost
surely charges every eigenspace.  Its support and weights determine the
unique secular root above `lambda_max(B_N)`, so it is an exact carrier for
the perturbed **top** eigenvalue.  For a general non-Haar vector orthogonal to
an eigenspace, one should also retain uncharged eigenvalues; the draft's Haar
hypothesis avoids that null event.

The concentration claim is correct.  For every bounded continuous `f`,

```math
u_N^Tf(B_N)u_N-N^{-1}\operatorname{tr}f(B_N)\to0
```

almost surely by spherical concentration and Borel--Cantelli.  Uniform
operator bounds and a countable dense test family give
`nu_(N,u_N) => mu`.  Applying the same argument to resolvents and a compact
exhaustion gives local uniform convergence of `G_(N,u_N)` on `(b,infinity)`.

## 2. Theorem 1, rank one: PASS

For `z>lambda_max(B_N)`, `G_(N,u_N)(z)` is positive and strictly decreasing.
Almost surely it diverges at the finite top eigenvalue, so the positive
rank-one perturbation has one secular root above that eigenvalue.  If

```math
\theta G_\mu(b+)>1,
```

strict monotonicity gives a unique `rho>b` with
`G_mu(rho)=1/theta`; local uniform convergence sends the finite root to
`rho`.  If `theta G_mu(b+)<=1`, then for every `epsilon>0`,
`theta G_mu(b+epsilon)<1`, which puts the root below `b+epsilon` eventually;
positive interlacing and `lambda_max(B_N)->b` give the lower bound.  Thus
(7)--(8), including the threshold equality case, are correct.

A useful normalization check is the semicircle law: `b=2`, the threshold is
`theta=1`, and the supercritical root is `theta+theta^(-1)`, as expected.

The cited primary source is appropriate.  The additive finite-rank phase
transition is Theorem 2.1, and convergence of the matrix-valued weighted
spectral measures used in its proof is Proposition 9.3, in
[Benaych-Georges and Nadakuditi, *The eigenvalues and eigenvectors of finite,
low rank perturbations of large random matrices*](https://arxiv.org/abs/0910.2120).
The draft's rank-one proof is sufficiently self-contained and does not rely
on a stronger random-base theorem.

## 3. Required repair 1: restrict or extend the spike signs

The model introduces an unrestricted multiset `Theta`, while Theorem 1 and
the displayed response law treat only `theta>0` and the upper edge.  For the
claims as written, require

```math
Theta subset (0,infinity).
```

If negative spikes are intended, the state must also retain the lower bulk
edge and use the lower-edge Cauchy-transform branch.  A negative finite-rank
deformation can produce a lower outlier but does not obey the upper-edge rule
(8).  “All outlying extremal responses” should therefore be replaced by
“all upper outliers from fixed positive spikes,” unless the two-sided formula
is supplied.

## 4. Finite-rank and repeated-spike composition: PASS only under joint
genericity

For a fixed Haar orthonormal frame `U_N`, the matrix secular equation is

```math
\det(I-\Theta U_N^T(zI-B_N)^{-1}U_N)=0.
```

Joint isotropic concentration gives

```math
U_N^T(zI-B_N)^{-1}U_N\to G_\mu(z)I_r,
```

including vanishing off-diagonal entries.  Hence the limiting secular
determinant factors as `prod_j(1-theta_j G_mu(z))`, with the asserted
multiplicities.

If a second independent Haar frame is added, its cross Gram and cross
resolvent blocks vanish, while each diagonal block converges to the same
scalar resolvent.  Orthogonalizing the combined fixed-size frame perturbs its
finite spike eigenvalues by `o(1)`.  Thus multiset union is an
**asymptotically exact** update for jointly independent generic frames:

```math
(\mu,\Theta)\star\Theta'
=(\mu,\Theta\uplus\Theta')+o(1)
```

at the response level.  It is not an exact finite-`N` identity.

The generic-direction hypothesis cannot be weakened silently.  Take
`B_N=0` and two unit spikes of strength one.  Two asymptotically orthogonal
directions give limiting top eigenvalue `1`, whereas two copies of the same
direction give `2`.  Both would be labelled by the formal multiset
`{1,1}` if relative orientation were discarded.  Thus the union law is
false for correlated/adversarial continuations; the missing datum is their
relative Gram geometry.  This is a decisive falsifier outside the stated
generic regime.

## 5. Required repair 2: this is a presented state, not an intrinsic
renormalized spectrum

The notation `S(B_N,P_N)` correctly exposes a decomposition.  The prose
calling `(mu,Theta)` “equivalently” a renormalization of finite-rank spectral
mass is too strong:

* `theta_j` is an eigenvalue of the **perturbation**, not generally an
  eigenvalue of the total matrix `B_N+P_N`;
* a supercritical total-matrix outlier is at `rho_j`, related by
  `G_mu(rho_j)=1/theta_j`;
* a subcritical spike leaves no separated atom in the total empirical
  spectrum from which `theta_j` could be recovered;
* one total matrix may admit more than one bulk-plus-spike presentation.

Therefore the state is a strict quotient of a **presented composition
history**, not automatically a statistic computable from the landscape
matrix alone.  This is legitimate in the repository's presented-carrier
language, but it must be declared.

Nor is the state finite-dimensional over an arbitrary bulk class: `mu` may
be any compactly supported probability measure and its full Cauchy transform
is needed.  The rare mark `Theta` has fixed dimension, and the whole state is
size-independent; call it “finite-dimensional” only when `mu` is fixed or
belongs to a finite-parameter family.

For the single top-eigenvalue query under further independent positive
spikes, `(mu,Theta)` is also not minimal.  Existing upper outliers merely
persist, so `(mu,current top)` updates by

```math
L\mapsto\max\{L,\mathcal R(\mu,\theta')\}.
```

The full multiset is justified only if the declared outputs include all
outlier locations/multiplicities, or if a later bulk-changing composition
makes latent strengths relevant.  No minimality theorem is proved here.

## 6. Deterministic-versus-random quantifiers: correct after a sharper
declaration

The almost-sure theorem implies existence of deterministic good sequences:
one may intersect the probability-one isotropic-resolvent events for a
countable collection of deterministic bulk sequences and a countable family
of predeclared generic direction frames, then freeze a realization.  The
same rooted-resolvent event actually handles all positive strengths, so the
rational-strength restriction is harmless but unnecessary.

Freezing does **not** convert the result into a uniform theorem over future
directions.  An adversary may choose a direction after seeing the frozen
matrices and align it with a bulk eigenvector or a previous spike.  For a
concrete example, let half the eigenvalues of `B_N` be `0` and half `-1` and
take `theta=1`.  Then `mu=(delta_0+delta_{-1})/2`, and the generic root is
`rho=1/sqrt2`; a spike inside the top eigenspace gives top eigenvalue `1`,
while a spike inside the bottom eigenspace leaves the top at `0`.  The same
`(mu,theta)` therefore has three different behaviors according to rooted
orientation.

The correct deterministic formulation is consequently:

> for each countable, predeclared family of generic direction sequences,
> there exists a deterministic frozen realization on which the state law
> holds simultaneously.

It is not a contextual congruence for all deterministic spike queries.

Also replace “conjugating by **any** orthogonal matrix makes both matrices
dense” by “a generic common orthogonal conjugation makes them dense almost
surely.”  The identity matrix is an immediate counterexample to “any.”

## 7. Theorem 2 and fixed-depth moments: PASS

A rank-one Hermitian perturbation changes the empirical distribution
function by at most `1/N`, and both matrix sequences have uniformly bounded
operator norm.  Equivalently, telescoping powers gives the quantitative
bound

```math
\left|N^{-1}\operatorname{tr}[(B_N+P_N)^k-B_N^k]\right|
\le {C_{k,\theta}\over N}
```

for every fixed `k`.  Thus (11) holds simultaneously for every fixed finite
depth.  Theorem 1 and the assumed edge convergence give (12).  All
normalizations are correct.

The conclusion should be phrased as follows: no summary whose limiting value
is a continuous function **only of the weak empirical spectral law of the
already composed matrix** can recover its top response uniformly on this
class.  By contrast, `mu` is sufficient for a known generic spike supplied
as an external query via (8).  The mark is needed to make weak-bulk
compression reusable after composition.  Declaring that input/query split
removes an apparent contradiction between Theorems 1 and 2.

The adjective “sharp” is not established in an information-theoretic sense.
The theorem gives a decisive fixed-depth falsifier but no matching transition
for growing moment depth (which can begin to see an outlier around logarithmic
depth in simple separated models) and no bit lower bound.  “Fixed-depth
insufficiency theorem” is accurate.

## 8. Collins--Male citation: correct theorem, missing hypotheses

[Collins and Male, *The strong asymptotic freeness of Haar and deterministic
matrices*](https://arxiv.org/abs/1105.4345) proves norm as well as trace
convergence when the deterministic tuple already has a **strong limiting
distribution**.  Their Corollary 2.2 then gives Hausdorff convergence of the
spectrum of independently unitarily invariant sums to the support of free
additive convolution.

Weak empirical convergence and a bound on `||B_N||` alone are not the cited
hypotheses.  To claim edge control for composed bulks, require strong
convergence of each bulk sequence (for a single self-adjoint sequence,
equivalently no spectral outliers and convergence of the relevant extremal
support data).  The draft correctly says this extension is unused, but its
one sentence should state this condition explicitly.  The paper's strongest
deterministic-matrix theorem is formulated over complex/unitary Haar
matrices; its orthogonal result in that paper covers independent Haar
orthogonal matrices, so a real orthogonal-conjugation application should be
cited/formulated with care rather than inferred without hypotheses.

## 9. Archive novelty and theory judgment

The repository already contains stronger Boolean-specific warnings that
complete spectra, fixed trace polynomials, and bounded-depth observables can
miss extremal responses (`obstruction_atlas_report.md`, especially rows S
and T).  Hence Theorem 2 is not a new general obstruction principle.  It is a
particularly clean zero-density spectral instance.

What is new to the archive is the **positive** benchmark:

```text
presented zero-density spike mark
 + Haar synchronization of the rooted resolvent
 + scalar secular equation
 -> reusable asymptotic response law under generic finite-rank additions.
```

This is not the Gaussian tangent-state theorem and is not merely the earlier
spectral anti-pin calculation.  Its mathematics is classical random-matrix
theory; the new contribution is its placement as a strict presented-carrier
example and its explicit synchronization boundary.

Under the near-minimizer campaign's existing taxonomy, this is best called a
**Level-3 positive benchmark** or a scoped rare-event portfolio example, not
a Level-4 exact-sign result: the entries are real, the configuration space is
the sphere, and the future directions are generic/predeclared rather than
adversarial.  It makes no frontier change for dense Boolean signings.

## 10. Final classification

| Claim | Verdict |
|---|---|
| determinant/secular normalization | PASS |
| Haar rooted-measure synchronization | PASS |
| rank-one threshold and response law | PASS |
| positive fixed-rank extension | PASS under joint Haar genericity |
| multiset-union composition | asymptotically correct, false for correlated directions |
| deterministic freezing | existence for countable predeclared generic contexts only |
| fixed-depth normalized-moment insufficiency | PASS |
| intrinsic/finite-dimensional state claim | REPAIR: presented state; arbitrary `mu` is infinite-dimensional |
| Collins--Male extension | REPAIR: add strong-convergence/no-outlier hypotheses |
| archive novelty | new positive benchmark, classical imported theorem, no Boolean frontier advance |

No counterexample was found within the repaired positive-spike, fixed-rank,
jointly Haar-generic hypotheses.  The explicit aligned-spike examples above
show that each genericity qualifier is mathematically necessary.
