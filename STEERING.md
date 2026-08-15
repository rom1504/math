# Strategic steering

Evidence cutoff: Section 10.140 (2026-08-15). Status: **autonomous route
generation paused after consolidation diagnostic; external review advised**.

If ordinary waves resume, keep Wave 61 as the next blank-slate boundary;
refresh earlier after a decisive result or architecture change.

## User-stated objective and workflow directives

The objective is to determine whether `M_n/n^(3/2)` converges; `1/2` is a
conjectural assessment, not another user objective. Reproducibility,
verification, Git checkpoints, both convergence and nonconvergence, and the
README stopping rule remain user directives.

The user's latest input supplied external research suggestions, not
mathematical directives. They were independently evaluated; route selection
and conjectures remain agent-authored. `README.md` is unchanged.

## Agent-authored frontier

```math
0.336493364431\ldots
\le \liminf_{n\to\infty}{M_n\over n^{3/2}}
\le \limsup_{n\to\infty}{M_n\over n^{3/2}}
\le {1\over2}.                                      \tag{S1}
```

No bound, recurrence step, convergence theorem, or genuine nonconvergence
mechanism improved at the last two substantive checkpoints.

## Strongest new theorem and its limit

For a pointwise-positive operator kernel on the edge cube whose Fourier
coefficients and complete remainder coefficients are Loewner-positive, the
full joint same-switch theorem is

```math
\boxed{(\lambda-\mu(a))_+T_a\preceq J.}             \tag{S2}
```

It preserves matrix cancellation. If the root is constant, however, a least
generalized eigenvector scalarizes the certificate with exactly the same
bound. Correct scale requires Fourier support `exp(Omega(n log n))`; rank
alone supplies no strength.

Partial matchings have that hidden size, but their direct uses are closed:

```math
\lambda-{J\over T_a}=O(n^{-1})
```

for radial matching squares, and

```math
\lambda\le {1\over n-1}
```

for arbitrary-rank kernels supported directly on matchings. Every pure
matching add/delete representation also has coefficient at most

```math
{4\over3\sqrt3}=0.769800\ldots<1.                   \tag{S3}
```

Broader difference-of-matchings blocks remain formally open, but no uniform
root theorem emerged and pure matching transitions cannot reach the
conference coefficient.

## Exact sufficient convergence statement

A moving-kernel family defined independently of `M_n` would prove convergence
to `1/2` if

```math
\lambda_n={1-o(1)\over\sqrt n},
\qquad
\sup_a\lambda_{\min}(T_a^{-1/2}J_nT_a^{-1/2})
=o(n^{-1/2}).                                      \tag{S4}
```

Then (S2) gives `Q(a)>=(1/2-o(1))n^(3/2)` uniformly and conference orders
give the upper subsequence. The root estimate can encode the full signed
coset histogram, so (S4) is not presently a strict reduction.

## Blank-slate abstraction audit

1. **Uniform block completion.** The proposed inequality was

   ```math
   D(A*_RB)\le\sqrt{N/p}D(A)+\sqrt{N/q}D(B)+O(N^{3/2-\delta}).
   ```

   It would make `M_n/sqrt(n)` almost subadditive. Comparison showed that it
   is exactly the existing state-dependent bridge objective, with a stronger
   universal child quantifier and the same rectangular norm obstruction.

2. **Bounded dual-cycle pressure.** Fixed-degree Eulerian cumulants were
   proposed to determine the minimized finite-temperature pressure. The exact
   pressure-to-ground-state mapping is valid. Fixed degree, however, misses
   planted zero-entropy resonances; degree `Theta(n)` restores the complete
   signed coset histogram. This is not a strict compression.

3. **Symplectic/Witt nonconvergence.** Orders `4^k` and `2*4^k` were proposed.
   The symplectic family has exact cap `(n/2)(sqrt(n)+1)`, so separation would
   require a bad subsequence above `1/2`. Dense Paley conference restrictions
   give `M_n<=(1/2+o(1))n^(3/2)` on every order and rigorously falsify it.

The candidates were respectively an equivalent obligation, an obstructed
bounded-state mechanism, and a falsified construction. None justifies a new
campaign. Details are in Section 10.140 and the diagnostic artifact.

## Inactive mechanisms

Selected-prior/common-active-face, scalar atoms, separately paid channels,
ordinary polarization, same-map Gaussian response, fixed-level SOS and
replicas, Walsh basin classification, radial moving kernels, scalar
transversals, direct matching support, pure matching add/delete modules,
uniform bridge sampling, fixed small disorder tilts, isolated algebraic
bridges, and eigenvalue-only polar variants remain inactive.

## Conditions for a justified restart

External review or a new architecture should supply at least one of:

1. a coefficient-one graph-orbit moving representation with an algebraic
   rooted mass theorem satisfying (S4);
2. a diffuse bridge law with `O(n)` entropy and a proved linear pressure gain;
3. a strict reduction of the existential bridge state to demonstrably less
   information than full parent optimization; or
4. a genuine nonconvergence mechanism producing two constants below `1/2` on
   multiplicatively separated epochs.

Suggestions must be checked against the ledger before becoming targets.
Vocabulary-only reformulations and additional finite caps do not qualify.

## Research-director decision

This is the second consecutive substantive checkpoint without primary
progress, followed by the required bounded diagnostic. The recent external
moving-representation idea was tested with executor, adversarial-verifier,
and director passes and yielded rigorous but route-limiting theorems.

Pause autonomous route generation and seek external mathematical review. This
does not claim the problem is solved or blocked. Resume only after explicit
authorization tied to a genuinely new input; do not open a nearby variant.
