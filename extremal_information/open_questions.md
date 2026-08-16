# Minimal open questions

These are selected theorem targets, not a catalogue of nearby ideas.  The
first two syndrome checkpoints produced Theorem 8.3; the remaining question
explicitly asks for complexity beyond that block construction.

## 1. Macroscopic syndrome-response rate--distortion

Fix `G=F_2^w`.  A full-rank parity-check fragment is operationally represented
by its nonzero column-type support `S subset G\{0}`.  Its response to an
appended fragment `E` is

```math
\mathcal R_S(E)=\rho(\ker[H_S\ E]).
```

Theorem 8.3 answers the first nontrivial normalized case.  On a direct-sum
block family, subset-selecting future fragments turn `Theta(w)` latent bits
into coherent radius differences.  For every fixed `epsilon<1/8`, error
`epsilon*w` has complexity `Theta_epsilon(w)` on that family, using only
linear-length fragments.

Determine the metric entropy for **arbitrary** spanning supports.

**Positive target:** a composition-stable quotient with
`exp(o(2^w))` possible states and uniform response error `o(w)` under repeated
support union/min-plus convolution.

**Negative target:** a realizable family with superlinear response-information
rate—and ultimately, if true, `exp(Omega(2^w))` supports—whose complete
future-response maps are pairwise `Omega(w)` apart without a supplied
direct-sum decomposition.

**Stop condition:** do not repeat the block subset-count construction or the
old unit-gap support probes.  A further lower bound must increase the source
complexity or expose a new interaction among blocks; an upper bound must
propagate its error under union.

## 2. Derive, rather than assume, deterministic linkage

Theorem 9.2 proves uniform zero-temperature synchronization from mixture
ultrametricity and `(D,tau)` monotone linkage of total-overlap pair labels.
The rare matching example proves that PSD, exchangeability, ultrametricity of
all nonnegative species mixtures, and vanishing conditional variance are not
enough.

Find a natural finite hypothesis, stated without inspecting the hidden species
values, that forces `tau+3D*eta=o(1)` on every pair label exposed by a declared
zero-temperature query.

**Success:** a model class in which the synchronized scalar profile is a
strict quotient and the theorem controls every declared maximum.

**Falsifier:** a family satisfying the proposed visible hypotheses with a
fixed response gap on an exposed fibre.  Average-only hypotheses are already
closed by Example 7.

## 3. Effective-rank composition beyond fixed dimension

Shapley--Folkman bounds the loss from replacing a nonconvex reachable sum by
its convex roof by the `r` largest component diameters, where `r` is the
effective affine-difference rank.  Determine whether **exposed** or
query-dependent dimension can replace `r`.

**Success:** a bound controlled by the dimension or Gaussian width of faces
actually reached by a declared Lipschitz query, with subextensive error in a
model where ambient rank grows.

**Falsifier:** a same-roof pair with low exposed complexity for every local
component but a leading composed response gap.  A representation of an
arbitrary convex body is not automatically succinct.

## 4. A normalized massive-witness theorem, or its impossibility

The query-weighted four-cell theorem is rigorous, but every matching-based
certificate on the canonical code transversal is exponentially small under
diffuse state-pair sampling.  The graph-code distance table even has
exponential exact tropical rank and a rank-one normalized approximation with
mean-square error `1/(16t)`.

Ask whether another natural conditional-response family admits a jointly
charged witness involving a positive fraction of its table, with no separate
payment of scalar cells or channels.

**Success:** a polynomially described witness yielding a nonvanishing
normalized average-error lower bound.

**Stop:** do not try more pairwise graphs on the same transversal.  If every
joint witness is equivalent to the full table, record an information lower
bound instead.

## 5. Constrained compactness and realization

Fixed-interface response bodies have unrestricted finite recovery sequences.
Characterize which limiting response states are realized at all large sizes
inside a constrained family such as linear codes, dense CSPs, or bounded-width
factor graphs.

**Success:** a Gamma-limsup/recovery theorem preserving the declared response
with vanishing normalized loss and without storing a target optimizer.

**Stop:** finite-state approximation outside the constrained model is already
known and does not answer this question.

## Reconnection rule

The program now has strict composable quotients in nontrivial restricted
models, but none closes the dense-sign bridge interface.  Do not return to the
motivating signing problem unless one of Questions 1--3 produces a state whose
complexity and error remain controlled when its interface rank grows.
