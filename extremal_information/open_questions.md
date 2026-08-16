# Minimal open questions

These questions are ordered by their ability to distinguish a real theory
from a renamed obstruction.  Each includes a stopping or falsification test.

## 1. Minimal rooted repair of pair entropy

The code pair `C^r,D^r` has identical unrooted energy--energy--distance data
and linearly separated worst-root response.  Determine the least rooted
augmentation that fixes every query

```math
\max_x\{H(x)-\lambda d(x,z)\}.
```

Concrete candidates are the multiset of rooted distance profiles, the upper
roof of `(H(x),d(x,z))`, or a finite association-scheme module.

**Success:** an exact sufficiency theorem with a state smaller than the full
root-by-state distance table, plus a tensor composition law.

**Stop:** an injectivity or packing theorem showing that uniform rooted
queries recover the complete landscape/code at the requested accuracy.

## 2. Can structure rescue a finite global replica hierarchy?

The universal version is false: for every fixed `k`, Theorem 3.3 constructs
code landscapes with identical complete unrooted data through `k` points and
a positive normalized covering-radius gap after tensorization.  The remaining
question is whether a meaningful structured class—dense quadratic signings,
models obeying synchronization identities, or a closed association-scheme
family—can be controlled by finite `k`.

For dense quadratic signs, retain the exact energy labels and all entries of
the global overlap matrix of `k` replicas.  Can two sequences have identical
such data but a separated intrinsic `(k+1)`-replica response, or does a
verifiable structural hypothesis force closure?

The response must be permutation-invariant; otherwise the fixed-half example
already answers the easier rooted question.

**Success:** a scalable dense-quadratic construction, or a theorem that a
meaningful nontrivial class synchronizes at finite `k`.

**Stop:** do not infer sufficiency from the absence of collisions at small
orders.  The exact dense-sign census currently reaches only order eight.

## 3. Query-generated feature-algebra growth

Starting from a feature family `F_0`, let `F_t` contain all observables created
by `t` allowed compositions and contractions.  Bound the response metric
entropy of its exposed roofs.

**Success:** for a nontrivial model, prove

```math
\log\operatorname{Cov}(\mathcal R_{n,t},\epsilon_n)
=o(|\Omega_n|)
```

at an error that preserves the normalized extreme, together with a closed
composition law.

**Stop:** if the algebra separates every Boolean state, prove this explicitly
and classify the resulting full-landscape rate rather than carrying the
hierarchy further.

## 4. Bounded-strength extremal rate--distortion

The quadratic pinned-query lower bound uses external fields of magnitude
`Theta(n)`, or bounded per-edge rank-one coupling interventions.  Determine
the rate for bounded coordinate fields, sparse interventions, and low-rank
bridges.

**Success:** a nonzero asymptotic rate function with matching upper and lower
bounds on a natural quadratic ensemble.

**Falsifier for a quadratic lower rate:** construct sketches of `o(n^2)` bits
that answer the declared bounded query family to the target normalized error.

## 5. A unified extremal information inequality

Current lower bounds use either response packing or posterior sign
polarization.  Seek an inequality relating query distortion, posterior
feature variance, and mutual information, for example

```math
I(H;Z)\ge \Psi
\bigl(\mathbb E\,d_\Theta(H,\widehat H_Z),
      \text{exposed-face geometry}\bigr).
```

It should specialize both to binary Hamming rate--distortion and to the
sign-near barycenter bound.

**Success:** a theorem with a nontrivial application to at least two of
quadratic forms, code covering radius, Max-Cut, or random CSPs.

**Stop:** a restatement of the definition of mutual information or a packing
bound with no model-specific geometry is not the requested unification.

## 6. Deterministic synchronization

The Curie block pair fails because total overlap does not determine species
overlaps.  Random multi-species spin glasses can obtain synchronization from
strong distributional identities.  Identify deterministic hypotheses under
which omitted block features are functions, approximately and uniformly, of
a retained global order parameter.

**Success:** a theorem with checkable finite hypotheses and an error bound
stable at zero temperature.

**Falsifier:** a deterministic family satisfying the proposed hypotheses but
retaining an order-one block-response gap.

## 7. Constrained compactness and recovery

Fixed-interface response bodies have unrestricted finite recovery sequences.
Which closed subsets arise from a constrained family—complete sign matrices,
linear codes, dense CSPs—at every sufficiently large size?

**Success:** a Gamma-limsup/recovery theorem that preserves the relevant
response roof with vanishing normalized loss and does not assume the target
finite optimum.

**Stop:** arbitrary finite-state approximation is already proved and does not
address the constraint.  A construction that stores a target optimizer is
circular for extremal compression.

## 8. Dense-sign pair-overlap collision

Find two nonisomorphic hollow `+/-1` quadratic sequences with identical
energy--energy--global-overlap data at a declared asymptotic resolution but a
separated anchored response.

**Success:** a scalable construction with an exact normalization and response
gap.

**Finite guide:** the exhaustive rooted-gauge census found no collision with
different one-vertex response multisets through order eight.  This is only a
search boundary, not evidence of impossibility.

## 9. Known-model validation

For each of the following, identify the query interface and determine whether
the standard order parameter is exactly an upper response roof, a quotient of
one, or something genuinely different:

- REM/GREM and the Parisi zero-temperature limit;
- random Max-Cut and dense CSP interpolation;
- code/coset covering radius;
- discrepancy and vector balancing; and
- Littlewood/unimodular polynomial sup norms.

**Success:** a known theorem becomes a short corollary of a common response or
rate--distortion statement.

**Stop:** shared terminology without a checked implication is not validation.

## 10. Conditions for returning to the motivating problem

Only after Questions 3, 6, or 7 produces a closed state with demonstrably
sub-landscape information should the theory be tested on the original bridge
composition.  The first target would then be a theorem saying that the child
states determine the parent cap with `o(n^(3/2))` error.

If the required interface separates every spin state or has quadratic edge
information rate, record that as an explanatory lower bound and do not market
it as a convergence reduction.
