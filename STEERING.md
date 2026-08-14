# Strategic steering

Evidence cutoff: regular-Hadamard Walsh/basin checkpoint (2026-08-14).
Status: **trajectory/basin campaign active; bent-only inverse inactive**.

## User-stated objective and workflow directives

The research objective is to determine whether `M_n/n^(3/2)` converges. The conjectural value `1/2` is not an additional user objective.

The user authorized sustained computational--composition research, then focused the active campaign on a joint same-switch certificate. Reproducible artifacts, independent verification, regular Git checkpoints, and the README consolidation/stopping rule remain workflow directives; they do not make mathematical conjectures or rankings user directives.

The latest directive requests a Walsh/bent-coordinate audit of the regular-Hadamard obstruction, a search for applicable bent/plateaued, automorphism, and Reed--Muller results, and symmetry-aware measurement of the explicit bad core's basin. The requested decision rule is to stop a Walsh-structure route if large cores are generic and continue it if they correlate strongly with bent/near-bent structure. The interpretation below is agent-authored.

## Agent-authored assessment

The rigorous frontier is unchanged:

```math
0.336493364431\ldots
\le \liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
\le \limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
\le \frac12.                                      \tag{S1}
```

The preceding joint same-switch certificate remains valid:

```math
Q(A)\ge \mathbb E\|AX\|_1-\mathbb E\Delta_{\rm gr}(A,X),\qquad
\frac{\mathbb E\Delta_{\rm gr}}{n^{3/2}}
\le C\alpha\sqrt{\log(2e/\alpha)}
+K\Pr\{\kappa_*>\alpha n\}.                       \tag{S2}
```

Thus `\kappa_*/n\to0` in probability would improve the lower constant to `1/\sqrt{2\pi}=0.398942\ldots`. It would improve the interval, not prove convergence.

## Walsh audit

For the quadratic Walsh matrix `H`, put `g=fX`. The exact row/energy dictionary is

```math
D_XCD_X=D_gHD_g-I,\qquad
X_u(CX)_u=g_u(Hg)_u-1,\qquad X^TCX=g^THg-m.       \tag{S3}
```

Masked forms of (S3) give both shores, the augmented field, initialization, every flip gain, terminal margins, unmatched mass, and `\kappa`; no parent Boolean maximization is hidden.

For a regularizing base `g` and anti-regularizing terminal `z`,

```math
\Delta=[E(g)+E(z)-2m]_+=[A(z)-B(g)-2m]_+,         \tag{S4}
```

where `A,B` are squared deviations of signed Walsh amplitudes from `\sqrt m`. A leading defect forces increased anti-flatness and fourth moment, not near-bentness.

The original input is self-dual bent and its displayed endpoint is anti-self-dual 2-plateaued, but that shape is not necessary:

- verified static project-scale families have simultaneous linear cores with full-support nonflat spectra;
- the prescribed least-index trajectory from a three-level non-bent tensor seed has certified `\kappa/m=0.0930862` at `k=9`, ending with Walsh support `262126/262144` and 440 absolute magnitudes;
- therefore “large core implies bent/plateaued” is false as a finite uniform mechanism. An asymptotic inverse would require a genuinely dynamical hypothesis.

The trajectory shadows an explicit tensor terminal. If all deviations lie in three medium fibres with density `\eta<1/8`, then

```math
\frac{\kappa}{m}\ge
\frac{1/8-\eta-2/\sqrt m}{3-2/\sqrt m}.           \tag{S5}
```

Its endpoint condition reduces to one self-consistent signed-Walsh threshold set on `m/16` tail coordinates. Proving fibre confinement with `\limsup\eta<1/8` is the strongest local target. It gives a scalable pointwise core, but one seed does not falsify the probability statement in (S2).

## Basin and literature boundary

At `k=2`, exact projective enumeration gives basin `70/32768` for the explicit endpoint and `571/32768` for its eight-point affine-orthogonal orbit. Least-index tie breaking is not equivariant; endpoint basin sizes range from 67 to 74. A symmetry-invariant random-tie control gives orbit probability `0.0140747`.

There were no explicit-orbit or anti-self-dual 2-plateaued hits in 50,000 local samples at `k=3` and no 2-plateaued hit in 10,000 at `k=4`. Uniform full two-block samples had two diffuse events in 20,000 at `k=2` and none in 3,000 at `k=3`. This favors basin decay in the tested family but is not an asymptotic theorem.

Primary literature classifies and counts exact self-dual and plateaued functions; affine-support 2-plateaued vectors and fixed small Hamming neighborhoods are exponentially sparse. No checked theorem converts a Walsh `L^1` or fourth-moment defect into proximity to that class, and endpoint sparsity does not bound deterministic preimages. Reed--Muller covering-radius results control the largest Walsh coefficient, not (S4) or greedy basins.

## Exact targets and falsification criteria

1. **Leading: trajectory/basin entropy.** Prove `\Pr\{\kappa_*>\alpha n\}\to0` uniformly for every fixed `\alpha>0` using a state retaining tensor-fibre multiplicities and lexicographic interleaving, or construct a project-scale family with nonvanishing bad probability. A single seed, endpoint count, or moment bound does not decide this.
2. **Local structured falsifier.** Prove or falsify fibre confinement and `\limsup\eta_k<1/8` for the non-bent tensor seed. This tests the strongest observed obstruction but cannot alone settle target 1.
3. **Convergence composition.** The clean interface remains a constructor `T_k`, defined without `M_n`, satisfying

   ```math
   \operatorname{cap}(T_k(A))
   \le(1+\eta_k)k^{3/2}\operatorname{cap}(A)+Cnk^{3/2},
   \qquad \eta_k\to0.                              \tag{S6}
   ```

   Applying (S6) to a near-liminf seed and filling a fixed remainder would prove `\limsup\le\liminf`. Scalar atoms, separately paid channels, canonical same-map Krivine rounding, fixed-level SOS, and exact full bridge optimization remain inactive.
4. **Genuine nonconvergence.** This requires fixed `\epsilon>0` and two asymptotically separated subsequences, or strict `\liminf<\limsup`. No candidate subsequences are known; route-specific cores do not count.

## Checkpoint decision

Stop the bent/plateaued stability route: exact spectral shape is neither necessary for static cores nor for the tested prescribed trajectory. Do not reopen moment-only, exact-plateaued, or Reed--Muller-deep-hole variants.

Continue only the dynamics-specific basin question, with the reduced tensor threshold law as a controlled adversarial test. No primary progress under the README definition occurred: (S1), the recurrence, and convergence status are unchanged. Wave 61 remains the next scheduled blank-slate boundary if ordinary waves resume.
