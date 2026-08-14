# Strategic steering

Evidence cutoff: fixed-temperature spectral-extremality checkpoint (2026-08-14).
Status: **finite-temperature regularization campaign active; uniform conference reverse-KL closed**.
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

The strongest current target is **pressure-preserving regularization**. For each fixed sufficiently small `beta`, transform an exact pressure minimizer `A_n^*` into a signing `A_tilde_n` such that

```math
\log\overline Z(\beta\widetilde A_n/\sqrt n)
\le P_n(\beta)+o_\beta(n),                            \tag{S7}
```

while `||beta A_tilde_n/sqrt(n)||_op<=1/2-eta_beta`, the averaged diagonal-monomial law is conference-like, and every off-diagonal diagonal monomial is `n^(-1/2+o(1))` uniformly. The stronger original fixed-power Assumption 2.9 also suffices. Then (S5) supplies the missing liminf, (S4) supplies the limsup, and `P_n(beta)/n->psi(beta)` in that high-temperature range.

This is demonstrably narrower than optimizing over all spectral laws: (S5) has removed that obligation, and the literature audit removes uniform diagonal matching under Frobenius-near structure. It is not yet known to be substantially easier than the full fixed-temperature minimization. Norm-only and Frobenius-only penalties are both falsified by localized flat-sign perturbations. A successful theorem must puncture/edit exceptional coordinates at `o(n)` pressure cost, or directly control off-diagonal open-cactus corrections.

Even complete success for `beta<1/2` would not prove (S1); the zero-temperature squeeze needs fixed `beta` arbitrarily large. Extension beyond the replica-symmetric high-temperature theorem is a separate bridge.

## Ranked alternatives and falsification criteria

1. **Regularization/traffic:** prove (S7) with uniform off-diagonal diagonal-monomial control, or construct a scalable signing family with pressure below `psi(beta)` and a quantified persistent irregularity. Stop on norm-only/Frobenius-only estimates or a merely spectral reformulation.
2. **Exponentially rare bridge:** give a compressed algebraic selector and a proved pressure/composition bound. Uniform or polynomial bridge search is closed by (S6).
3. **Large-temperature thermodynamic limit:** find an interpolation or variational theorem valid beyond operator norm `1/2`; a small-`beta` analytic continuation claim is insufficient.
4. **Direct composition:** retain an explicit complete-support constructor with geometrically summable normalized defect. No such landing theorem is known.
5. **Genuine nonconvergence:** produce fixed `epsilon>0` and two infinite subsequences separated by `epsilon`, or prove `liminf<limsup`. Route-specific failures do not count.

Walsh basin, scalar finite-temperature recurrences, conference KL lower bounds, selected-prior/common-active-face, scalar atoms, separately paid channels, canonical same-map Krivine rounding, and fixed-level SOS remain inactive.

## Checkpoint decision

The imported theorem and (S5) are material theorem-level progress: they remove the conference uniform-remainder question and the high-temperature spectral optimization. They do not improve a ground-state bound, complete a recurrence, or prove convergence/nonconvergence. Continue only the bounded regularization/traffic test and the exact literature boundary. If neither yields a pressure-preserving theorem nor a scalable irregular counterexample, consolidate the finite-temperature campaign rather than returning to rapid scalar variants.
