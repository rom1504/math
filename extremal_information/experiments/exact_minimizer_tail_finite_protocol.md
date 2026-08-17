# Frozen protocol: finite audit of `L_tail`

**Status:** preregistered finite protocol; no asymptotic inference.

This protocol was fixed before running the task-local tail script.  The only
matrix corpus used for low-cap objects is

```text
extremal_information/experiments/nearmin_blind_structural_results.json
sha256 2c086cf7523ead804942948e800c6231eac33d954e049b5aa113c9fb0cca47a5
```

The target statement is Section 5 of
`drafts/exact_minimizer_switching_broadcast.md`, frozen for this audit at
SHA-256
`87f517bdf945f71e2c45a82bcdb75bee4c55d5498b067a86845339b0b7a5c5ea`.
No new minimizer or near-minimizer search is permitted.

**Post-freeze traceability note.**  During the run, that source received a
scope-only repair and now has SHA-256
`37cc9807b847ab1b2935dd5f5d1f9c2dfd9a0630e536c37d761154e49c3bcd18`.
The formulas SB.8--SB.9 audited here did not change.  This note records the
new source identity; it changes none of the predeclared thresholds, strata,
orientation rules, or output metrics below.

After its independent proof audit, the source received only a status-line
update and has final SHA-256
`5278e6cb96a3a554141fe52cfd31dbb1ca38cf7b1260a33a554c116bf6074e8f`.
Again, SB.8--SB.9 and every input to this finite protocol are unchanged.

## 1. Predeclared thresholds

For a hollow order-`n` signing, put

```math
H_A(x)=\sum_{i<j}A_{ij}x_ix_j,
\qquad Q(A)=\max_x|H_A(x)|.
```

The complete, fixed threshold grid is

```text
d0 in {1/64, 1/32, 1/16, 1/8}.
```

No post-hoc threshold will be added to the principal tables.  The largest
threshold `1/8` is the one used in the archived Walsh application of
Theorem 21.8; the smaller values diagnose the severe lattice effect at
orders at most fourteen.

## 2. Orientation rule

Let `P_+=max_x H_A(x)` and `P_-=-min_x H_A(x)`, so
`Q(A)=max(P_+,P_-)`.  The matrix is globally negated when only `P_-` equals
`Q(A)`.  When `P_+=P_-=Q(A)`, both global orientations are admissible.  To
search for counterexamples rather than flatter `L_tail`, the principal
record uses the orientation with the **larger** upper-tail density at each
`d0` (equivalently the smaller entropy-deficit exponent).  Both directional
counts are retained.

The full-cube and projective-cube densities coincide because every quadratic
energy is invariant under `x -> -x`.  Computation may therefore enumerate
the `2^(n-1)` projective spins with the first coordinate fixed.

For the selected orientation define

```math
p_(A,d0)=2^(-(n-1))
 #\{x\text{ projective}:Q(A)-H_(+-A)(x)<d0 n^(3/2)\},
```

and report

```math
I_bits(A,d0)=-{1\over n}\log_2p_(A,d0),
\qquad
I_nats(A,d0)=I_bits(A,d0)\log2.
```

The finite analogue of the `kappa` in `L_tail` is `I_nats`: indeed the
full-cube count is exactly `exp((log 2-I_nats)n)`.  The strict inequality in
the shell definition is evaluated exactly for rational `d0=p/q` through

```math
(q\,[Q-H])^2<p^2n^3,
```

so no floating threshold decision is allowed.

The spectral diagnostic is

```math
S(A)={\|A\|_(2 to2)\over\sqrt n},
```

where the norm is the maximum absolute eigenvalue of the real symmetric
matrix.

## 3. Frozen strata

The principal nonexhaustive strata are:

1. `repository_exact`: every byte-distinct matrix in the frozen
   `repository_exact_representatives` list;
2. `repository_one_step_near`: every byte-distinct matrix in the frozen
   `repository_one_step_near_representatives` list;
3. `adversarial_low_cap`: every byte-distinct matrix in the frozen
   `cap_constrained_adversarial_samples` list;
4. `greedy_low_cap`: the separately labelled frozen
   `independently_generated_greedy_low_cap` list;
5. `uniform_random`: the original 48 unconditioned controls at every order
   `3,...,14`, reconstructed from the frozen seed `20260817` and sample
   count in the JSON protocol; their cap histograms must exactly reproduce
   the frozen `population_summaries.uniform_random` histograms;
6. `cyclic_structured`: all cyclic-distance matrices from the original
   deterministic control family, reconstructed at every order and checked
   against the frozen cyclic cap histograms.

Reconstructing the two control families is permitted because it performs no
optimization and uses their already frozen generator, seed, and size.  A
control cross-check failure aborts the run.  Uniform controls are not
conditioned on cap.  The frozen low-cap random subset is also reported as
`uniform_random_low_cap` to separate conditioning from the true control.

Records are deduplicated by the SHA-256 hash of the integer matrix.  This
removes byte-identical repetitions but is not claimed to quotient switching
or permutation symmetry.  Per-order tables are primary; pooled summaries
are descriptive only.

## 4. Exhaustive small orders

For every `3<=n<=7`, enumerate the complete root-gauged population obtained
by fixing the first row positive and varying the remaining
`binom(n-1,2)` edges.  Partition it by exact cap.  Report at least the exact
stratum `Q=M_n` and one-step stratum `Q=M_n+2`, including counts, tail-rate
ranges for every predeclared `d0`, and spectral ranges.  The exhaustive
minimizer counts must reproduce the frozen counts

```text
2, 6, 12, 12, 3240.
```

This enumeration is classification, not a new search for larger-order
minimizers.

## 5. Counterexample-first reporting

For every stratum/order/threshold, report the smallest `I_nats` (largest
tail) and its matrix hash.  For the spectral condition report the largest
`S(A)` and its hash.  Preserve an explicit matrix for each global worst
exact and one-step witness.

The audit may say only:

- whether either candidate is already falsified at a finite order;
- how broad the observed finite constants are;
- whether exact representatives differ descriptively from controls; and
- what larger-order certificate would discriminate the hypotheses.

It may **not** infer a uniform positive asymptotic `kappa`, a uniform
operator constant, physical response packing, or convergence of any
sequence from orders at most fourteen.
