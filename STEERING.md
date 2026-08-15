# Strategic steering

Evidence cutoff: Section 10.138 (2026-08-15). Status: **bounded high-rank
moving-representation attempt active; scalar fibers and small tilts closed**.

Next mandatory blank-slate refresh: Wave 61 if ordinary waves resume, or five
substantial campaign checkpoints after the last full refresh; earlier after
a decisive theorem, counterexample, or route change.

## User-stated objective and workflow directives

The objective is to determine whether `M_n/n^(3/2)` converges; the
conjectural value `1/2` is not another user objective. Reproducibility,
verification, Git checkpoints, and README stopping rules remain directives.
The latest input supplied external research feedback, not mathematical
directives. Route rankings below are agent-authored; `README.md` is unchanged.

## Agent-authored assessment

The rigorous frontier is unchanged:

```math
0.336493364431\ldots
\le \liminf_{n\to\infty}{M_n\over n^{3/2}}
\le \limsup_{n\to\infty}{M_n\over n^{3/2}}
\le {1\over2}.                                      \tag{S1}
```

## New moving-representation theorem

For `E=binom(n,2)`, augmented cut code `C_n^+`, and `mu(a)=Q(a)/E`, the
imported moving projections have overlap `kappa>=0` and unsplit remainder

```math
F=(\tau-\lambda)\kappa\succeq_{\rm PD}0.
```

```math
T_a=\sum_{c\in C_n^+}\kappa(ac),
\qquad J_n=\sum_{c\in C_n^+}F(c),
```

subgroup Fourier positivity proves

```math
\boxed{(\lambda-\mu(a))_+T_a\le J_n.}               \tag{S2}
```

The published packing theorem alone is coset-blind: `C_n^+` has minimum
distance `n-1`, and every translate has identical internal distances.  The
new information is the rooted cross term in (S2), not the packing bound.

The exact convergence lemma is a family defined without `M_n` satisfying

```math
\lambda_n={1-o(1)\over\sqrt n},
\qquad
\sup_a{J_n\over T_a}=o(n^{-1/2}).                    \tag{S3}
```

(S2) then gives `Q(a)>=(1/2-o(1))n^(3/2)` uniformly, and conference upper
bounds prove convergence. The remaining root-mass bound is not demonstrably
simpler than the bare tail. Canonical rank-one kernels give only RMS through
tested order 18; their higher-level root mass collapses.

Any rank-one partial dual-transversal has support at most `2^n`, so the
Bollobás--Lee--Letzter cube spectral theorem forces

```math
\lambda=O((n\log n)^{-1/2})=o(n^{-1/2}).             \tag{S4}
```

Thus scalar transversal fibers cannot reach (S3), even with perfect
quotient alignment.  Correct scale needs `exp(Omega(n log n))` hidden
support or a genuinely high-rank/nonabelian representation.

## Leading route and falsification criterion

The leading target is a **cut-specific high-rank moving fiber**: an
`S_n`-equivariant family with possibly `exp(Omega(n log n))` hidden dimension,
but polynomially closed `J_n,T_a` satisfying (S3).

The state is genuinely weaker than parent maximization only if its
definition is independent of `a` and its root-mass estimate follows from an
algebraic identity or uniform annular bound, not from enumerating the signed
Eulerian/coset histogram.  Reject a candidate if:

- it is radial rank one or an `exp(O(n))` scalar transversal;
- its `T_a` collapses on a scalable signing;
- `J_n/T_a` remains above `lambda_n` by fixed `n^(-1/2)` scale; or
- its purported compressed state invertibly reconstructs the full coset
  energy distribution.

The focused literature check found no Delsarte, Terwilliger, cocycle-code, or
signed-graph theorem supplying this root-mass bound.

## Entropy-tilted bridge alternative

For bridge pressure `L` and uniform output law `U`, the exact soft minimum is

```math
\mathcal R_\lambda
=-{1\over\lambda}\log\mathbb E_Ue^{-\lambda L}
=\log\mathbb E_Ue^L-D_{1+\lambda}(U\Vert\Pi).        \tag{S5}
```

Its reveal recursion is an entropy-weighted martingale. An `exp(-O(N))`
target basin with growing tilt gives a summable composition defect.

On conference children, every small fixed tilt retains the full linear
defect, while the known bridge has `Theta(r^2)` entropy. Keep only a diffuse
law with

```math
D(q_r\Vert U)=O(r)
```

and a linear pressure gain.  Uniform bridges, polynomial sampling, small
fixed tilts, and isolated algebraic orbits are closed.

## Ranked alternatives

1. **High-rank rooted moving representation:** prove (S3), or establish a
   scalable obstruction for every algebraically closed `S_n` family.
2. **Diffuse entropy-tilted bridge law:** prove the `O(r)`-entropy pressure
   gain, or a matching entropy-rate no-go.
3. **Polar pressure:** return only with a frame-dependent involution theorem
   or dense linear counterexample; landing and large `beta` remain open.
4. **Direct composition/large temperature:** find a summable cross-order
   defect without scalar temperature contraction.
5. **Genuine nonconvergence:** produce fixed `epsilon>0` and two infinite
   subsequences separated by `epsilon`, or prove `liminf<limsup`.

Walsh basins, selected-prior/common-active-face, scalar atoms, separately
paid channels, canonical same-map Krivine rounding, fixed-level SOS, radial
rank-one moving kernels, scalar transversals, and small fixed disorder tilts
remain inactive.

## Checkpoint decision

No ground-state bound, recurrence, convergence theorem, or nonconvergence
mechanism improved. The root-mass lemma may be equivalent in difficulty, so
this is the first no-primary-progress checkpoint after the Frobenius reset.

Continue one bounded high-rank construction attempt.  If the next
substantive checkpoint also yields no primary progress, apply the README
consolidation rule and a blank-slate diagnostic before opening further
variants.
