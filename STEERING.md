# Strategic steering

Evidence cutoff: finite-temperature minimax checkpoint (2026-08-14).
Status: **finite-temperature joint-cancellation campaign active; Walsh basin archived**.

## User-stated objective and workflow directives

The research objective is to determine whether `M_n/n^(3/2)` converges. The conjectural value `1/2` is not an additional user objective.

Reproducible artifacts, independent verification, regular Git checkpoints, and the README consolidation/stopping rule remain workflow directives. Mathematical conjectures and route rankings below are agent-authored.

The latest user input supplied three strategic tests: evaluate the finite-temperature pressure as a possible convergence mechanism, audit the quadratic-Littlewood literature, and perform one final bounded full-basin experiment at `k=4,5,6`. It recommended ending the Walsh/bent/local-search line unless that experiment found a persistent basin. The route judgment below is the agent's assessment of those tests.

## Agent-authored assessment

The rigorous frontier is unchanged:

```math
0.336493364431\ldots
\le \liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
\le \limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
\le \frac12.                                      \tag{S1}
```

## Finite-temperature reduction and scalar boundary

With `H_A(x)=sum_(i<j) a_(ij)x_ix_j`, define

```math
\Phi_n(\beta)=\frac1n\min_A\log\sum_x
2\cosh\!\left(\frac{\beta H_A(x)}{\sqrt n}\right).
```

The exact soft-max squeeze implies that convergence of `Phi_n(beta)` for every fixed `beta>0` would force convergence in (S1), by taking `beta -> infinity` after the `n` limit.

This pressure is the repository's earlier soft-cap pressure plus its entropy normalization. For

```math
F_n(t)=\min_A\log\left(2^{-n}\sum_x\cosh(tH_A(x))\right),
\qquad R_n(t)=F_n(t)-{n\choose2}\log\cosh t,
```

bridge signs plus one child orientation give

```math
R_{m+n}(t)\le R_m(t)+R_n(t).
```

On the required diagonal this becomes

```math
r_{m+n}(\beta)
\le\theta r_m(\beta\sqrt\theta)
 +(1-\theta)r_n(\beta\sqrt{1-\theta}),           \tag{S2}
```

so balanced composition contracts `beta` to `beta/sqrt(2)`. Exact centering does not fix that. A proved analytic countermodel satisfies the entropy squeeze, convexity, correct variance and zero-temperature slope, restriction monotonicity, uniform diagonal Lipschitz control, and centered subadditivity, while its diagonal pressure oscillates. The natural quadratic scale-transport inequality is also false for an actual order-four signing. Scalar Fekete/interpolation variants are therefore inactive.

## Leading joint target

Put `P_n(beta)=F_n(beta/sqrt(n))`. For contracted-temperature child minimizers, let `T_(m,n)` be the two child temperature increments and let `D_(m,n)=D_KL(U||Pi)` be the reverse KL divergence of the joint orientation-and-noisy-rank-one bridge output from uniform. The exact quenched identity gives

```math
P_{m+n}(\beta)\le P_m(\beta)+P_n(\beta)+G_{m,n}(\beta),
```

```math
G_{m,n}(\beta)=mn\log\cosh\!\left(\frac\beta{\sqrt{m+n}}\right)
-T_{m,n}(\beta)-D_{m,n}(\beta).                  \tag{S3}
```

The exact sufficient lemma is: for every fixed `beta>0`, uniformly on `N/4 <= m,n <= 3N/4`,

```math
G_{m,n}(\beta)\le C_\beta N^{1-\delta_\beta}.    \tag{S4}
```

A balanced-tree Hammersley argument makes (S4) imply convergence of `P_n(beta)/n`; the soft-max squeeze then proves convergence in (S1). This state is genuinely joint and uses an averaged logarithm rather than a parent extremal spin. It is nevertheless stronger than existence of one good bridge and still has `2^(mn+1)` output sectors, so (S4) is a precise sufficient reformulation, not yet a verified reduction.

Exact enumeration through total order seven finds positive `G` in every tested split. For `beta<=1`, the balanced values are nearly order-independent and remain compatible with (S4); at larger `beta`, exceptional parents already beat the uniform average, showing that (S4) is strictly stronger than scalar same-temperature subadditivity.

Conference children give the main obstruction. Their average-log compensation margin has Taylor expansion

```math
-\frac{\beta^2}{4}
-\frac{9r^2-25r+15}{48r}\beta^4+O_r(\beta^6),    \tag{S5}
```

while reverse KL first appears at order `beta^8` with bounded coefficient. Thus the adverse fourth-order coefficient is linear in the parent order. The remainder is not uniform in `r`, and conference children are not proved pressure minimizers along a fixed-`beta` infinite family, so (S5) falsifies only a conference/perturbative implementation, not (S4).

## Closed and alternative tracks

- The final uniform full two-block experiment found zero positive defects and zero nonzero unmatched cores in `10,000`, `2,000`, and `250` samples at `N=512,2048,8192`. This favors basin decay but is finite evidence. The Walsh/bent/greedy-basin line is now archived.
- The exact disorder-counting product theorem constructs exponentially many parents, but retains the same contracted child parameters and misses an order-`N` endpoint term. A bare speed-`n^2` disorder LDP does not locate the needed support edge.
- Quadratic Littlewood, Boolean BH, Sidon, KSZ, and current SK theorems supply no `1+o(1)` complete-support cross-order law. Guerra--Toninelli averages Gaussian disorder; it does not commute with the adversarial signing minimum.
- Selected-prior/common-active-face, scalar atoms, separately paid channels, canonical same-map Krivine rounding, fixed-level SOS, and full bridge optimization remain inactive.

## Targets and falsification criteria

1. **Leading, falsification first:** determine the scale of `G_(m,n)(beta)` for exact low-temperature minimizers. Prove (S4), or prove for some fixed `beta,c_beta>0` and balanced infinite sequence that `G>=c_beta N`. Conference calculations count only after a uniform remainder theorem and a verified minimizer link.
2. **If the uniform-output state fails:** seek a bounded-complexity overlap or disorder-support state that controls a lower quantile of bridge log pressure, not its full minimum and not merely its uniform mean.
3. **Direct composition alternative:** retain an explicit complete-support constructor with a geometrically summable normalized defect. No such constructor is known.
4. **Genuine nonconvergence:** produce fixed `epsilon>0` and two infinite subsequences separated by `epsilon`, or prove `liminf<limsup`. Route-specific failures do not count.

## Checkpoint decision

No rigorous bound, recurrence step, or convergence/nonconvergence mechanism improved, so this is not primary progress under the README definition. Continuing is justified only as a bounded test of (S4), because it is the one surviving joint, same-temperature statement with a correct leading scale and an immediate scalable falsifier. Do not resume Walsh-basin work or scalar finite-temperature variants. If (S4) acquires a linear fixed-`beta` obstruction without a compressed replacement state, consolidate rather than generate another rapid wave.
