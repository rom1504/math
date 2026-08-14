# Strategic steering

Evidence cutoff: distributed-irregularity/Frobenius-stability checkpoint (2026-08-14).
Status: **finite-temperature averaged-pressure campaign active; coordinate puncturing closed**.
Next mandatory blank-slate refresh: Wave 61 if ordinary waves resume, or after
five substantial campaign checkpoints from this one, whichever comes first.

## User-stated objective and workflow directives

The research objective is to determine whether `M_n/n^(3/2)` converges. The conjectural value `1/2` is not an additional user objective.

Reproducible artifacts, independent verification, regular Git checkpoints, and the README consolidation/stopping rule remain workflow directives. Mathematical conjectures and route rankings below are agent-authored.

The latest user input supplied the fixed-temperature pressure reduction as a strategic lead, recommended a quadratic-Littlewood literature audit, and requested that the Walsh/bent basin line end unless its final bounded experiment found a persistent basin. It did not direct the mathematical conclusion or make any route a user objective.

## Agent-authored assessment

The rigorous ground-state frontier is unchanged:

```math
0.336493364431\ldots
\le \liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
\le \limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
\le \frac12.                                      \tag{S1}
```

For

```math
P_n(\beta)=\min_A\log\left(2^{-n}\sum_x
\cosh\!\left({\beta H_A(x)\over\sqrt n}\right)\right),
```

convergence of `P_n(beta)/n` for every fixed `beta>0` would imply convergence in (S1) through the exact soft-max squeeze. Scalar bridge subadditivity still contracts `beta`, and the analytic countermodel plus the exact order-four falsifier show that centering, convexity, restriction monotonicity, and scalar interpolation do not repair this.

## New fixed-temperature theorems

The July 2026 theorem of Fan--Misiakiewicz--Wang--Wen applies exactly to normalized symmetric conference matrices. Define

```math
\psi(c)=\frac14\left[
\sqrt{1+4c^2}-1-
\log\frac{1+\sqrt{1+4c^2}}2\right].                \tag{S2}
```

Along Paley conference orders,

```math
\frac1r\log\overline Z_r(A_r,\beta/\sqrt r)
\longrightarrow\psi(\beta),
\qquad 0<\beta<\frac12.                            \tag{S3}
```

Nearby Paley orders plus restriction monotonicity give the unconditional bound

```math
\limsup_n\frac{P_n(\beta)}n\le\psi(\beta),
\qquad 0<\beta<\frac12.                            \tag{S4}
```

There is also a new exact spectral theorem. If `mu` is any compact mean-zero law of variance `c^2` whose real R-transform branches reach one, and

```math
J(\mu)=\frac12\int_0^1R_\mu(u)\,du,
```

then

```math
J(\mu)+J(-\mu)\ge2\psi(c),
\qquad
\max\{J(\mu),J(-\mu)\}\ge\psi(c).                 \tag{S5}
```

Equality in the maximum is unique to the symmetric Bernoulli law. The proof is a variational Cauchy-transform identity followed by product Jensen and AM--GM. The `cosh` partition selects the maximum of the two orientations, so no spectral-symmetry assumption is hidden.

Consequently, every signing sequence with scaled operator norm uniformly below `1/2` and the fixed-power delocalization of Assumption 2.9 has limiting pressure at least `psi(beta)`. A sequence below `psi(beta)-delta` must approach the norm threshold or have a persistent fixed-power entrywise/traffic irregularity.

The available signed-permutation universality theorem weakens uniform diagonal control to convergence of the averaged diagonal-monomial distribution, but retains uniform `n^(-1/2+o(1))` bounds for every off-diagonal diagonal monomial. Frobenius-near conference structure supplies the averaged part. An explicit `n^(3/4)`-edge conference star perturbation remains Frobenius-near and norm-near while creating an `n^(-1/4)` entry of the squared scaled matrix, so Frobenius closeness alone cannot supply the missing uniform condition.

That obstruction is now scalable across coordinates: a distributed conference-star construction has conference pressure, the full averaged diagonal-monomial law, and `n^(-1/2+eta)` entries of `X^2` on a linear matching. Every endpoint cover is linear. Thus `o(n)`-coordinate puncturing is false even for correct-pressure sequences, and maximum-entry delocalization is not necessary for the thermodynamic answer.

The positive replacement is a dimension-free high-temperature interpolation inequality. Under a common operator margin `kappa<1/2`, `|log Zbar(X)-log Zbar(Y)| <= (K_kappa/2)||X-Y||_*`; hence normalized-Frobenius `o(1)` distance costs `o(n)` pressure. If `X^2` is normalized-Frobenius close to `beta^2 I`, pressure transfers to the polar involution `beta sign(X)`. This is a verified reduction unavailable in the previous assessment.

## Uniform joint reverse-KL route

For two conference children and a uniform sign bridge, the random parent law is

```math
\tfrac12(\delta_{-1/\sqrt2}+\delta_{1/\sqrt2})
\boxplus {\rm SC}(1/2).
```

The parent also satisfies the deterministic theorem pathwise. For every fixed `0<beta<sqrt(2)/6`, the exact full joint output obeys

```math
\frac{D_{\rm KL}(U\Vert\Pi)}r\to0,
\qquad
\frac{G_{r,r}(\beta)}r\to
\gamma(\beta)=\frac{\beta^2}{4}-2\psi(\beta)
+2\psi(\beta/\sqrt2)>0.                            \tag{S6}
```

Moreover, a uniform bridge reaches the same-temperature child target with probability at most `exp(-c_beta r)`. This converts the earlier nonuniform fourth-order tangent into a scalable fixed-temperature obstruction. It closes uniform averaging, fixed quantiles, polynomial sampling, and every compressed state that merely lower-bounds the same reverse KL on conference children.

This does not falsify the minimizer-optimized criterion: conference signings are not known to be exact pressure minimizers, and an exponentially rare algebraic bridge may exist. Do not spend further effort lower-bounding this KL on conference children; even the full divergence is sublinear.

## Leading target and implication chain

The strongest current target is the **polar-involution pressure theorem**. For fixed `beta<1/2`, prove that every symmetric involution `U_n^2=beta^2 I` arising at normalized-Frobenius distance `o(1)` from a flat hollow signing interaction satisfies

```math
\log\overline Z(U_n)\ge n\psi(\beta)-o_\beta(n).       \tag{S7}
```

If exact pressure minimizers can additionally be moved at `o(n)` pressure cost into the strict-norm Frobenius-near-conference class, Frobenius stability and (S7) supply the missing liminf, while (S4) supplies the limsup. This proves `P_n(beta)/n->psi(beta)` in that high-temperature range.

This target is narrower than entrywise universality: the proved Frobenius inequality removes every maximum-entry repair obligation and permits distributed edits. It is not an eigenvalue-only statement. A padded active edge makes the exact rms, spherical, and determinant lower bounds false, although only by `O(1)`, while the best isotropic Gaussian trial for an involution misses `psi(beta)` by `beta^4/12+O(beta^6)`. Any proof of (S7) must use frame-dependent dense flat-hollow information or allow an `o(n)` localization correction. A separate minimizer-to-Frobenius bridge remains open, so the implication chain has two genuine obligations.

Even complete success for `beta<1/2` would not prove (S1); the zero-temperature squeeze needs fixed `beta` arbitrarily large. Extension beyond the replica-symmetric high-temperature theorem is a separate bridge.

## Ranked alternatives and falsification criteria

1. **Averaged/polar pressure:** prove (S7), or construct a dense flat-derived polar family with pressure below `psi(beta)` by a linear amount. An `O(1)` padded-core deficit does not falsify it.
2. **Minimizer landing:** put exact pressure minimizers into the strict-norm Frobenius-near-conference class at `o(n)` pressure cost. Coordinate puncturing and maximum-entry repair are closed; a global Frobenius edit is allowed.
3. **Exponentially rare bridge:** give a compressed algebraic selector and a proved pressure/composition bound. Uniform or polynomial bridge search is closed by (S6).
4. **Large-temperature/direct composition:** find a theorem beyond operator norm `1/2`, or an explicit constructor with summable normalized defect.
5. **Genuine nonconvergence:** produce fixed `epsilon>0` and two infinite subsequences separated by `epsilon`, or prove `liminf<limsup`. Route-specific failures do not count.

Walsh basin, scalar finite-temperature recurrences, conference KL lower bounds, selected-prior/common-active-face, scalar atoms, separately paid channels, canonical same-map Krivine rounding, and fixed-level SOS remain inactive.

## Checkpoint decision

The Frobenius pressure inequality is primary theorem-level progress and the distributed matching is a scalable falsifier of the coordinate architecture. They do not improve a ground-state bound, complete a recurrence, or prove convergence/nonconvergence. Continue only the bounded polar-pressure test and the minimizer-to-Frobenius bridge. If neither yields a theorem or a linear dense-family counterexample at the next substantive checkpoint, consolidate the finite-temperature campaign rather than returning to maximum-entry traffic repairs.
