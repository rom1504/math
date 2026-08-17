# Audit: Boolean-port convolution and reuse

**Object audited.** `boolean_port_convolution_reuse.md`, Theorems CR.1--CR.3.

**Verdict.** PASS.  The positive sampling theorem is correctly limited to a
declared occurrence tree; the counterexample kills diagonal reuse of one
sample bank, not every possible self-convolution algorithm.

## Semantic algebra

Expanding definitions gives

```math
R_{\mu*\lambda}(\epsilon)
=\sum_{s,t}\mu(s)\lambda(t)K(st\epsilon)
=\sum_t\lambda(t)R_\mu(t\epsilon).
```

This proves CR.2 with no normalization factor.  A probability average is an
`ell_infty` contraction.  Replacing product factors one at a time therefore
gives CR.4.

On the two-point subgroup `{e,a}`, the convolution power is precisely the
parity law.  Since the response difference of a two-point mixture from
`delta_e` is a positive scalar multiple of one fixed response vector, its
supremum norm scales by exactly that scalar.  Thus `q_L(t)/t -> L`; the
generic telescoping coefficient is genuinely locally sharp.

## Sampling and reuse distinction

In a declared occurrence tree, for fixed replica `ell` all leaf-occurrence
samples are independent.  Their product has the required convolution law.
Different replicas use disjoint randomness, so each node marginal consists
of `k` iid samples.  RC.1 and its bounded-difference tail apply to that
marginal, and a union bound does not require different nodes to be
independent.  CR.11 is therefore valid.

The state represents a normalized histogram.  For actual tensor products,
storing the row count separately is sufficient: counts multiply exactly and
the unnormalized response is `pN R`.  No row-mass normalization is silently
lost.

If a semantic leaf is used twice, the theorem explicitly requires two leaf
**occurrences** with independent banks.  Feeding one bank into both inputs
instead produces `S_ell^2=e`.  Uniform measure on `{e,a}` is convolution
idempotent but lies at response distance
`d(delta_e,delta_a)/2` from `delta_e`; this is a fixed error whenever `a` has
linear projective Hamming weight.  The full-uniform calculation is also
correct: its response is constant
`c_p=E|sum_iX_i|/p`, and for `p>=2`, `c_p<=1/2`, so the maximum deviation from
the point response is `1-c_p`.

## Doeblin contraction

Uniform averaging of every translated kernel is the same, hence a difference
of two response roofs has uniform mean zero.  Decomposing

```math
\lambda=\alpha u_p+(1-\alpha)\lambda'
```

annihilates the first component and leaves a probability averaging operator
of norm at most one.  This proves the exact factor `1-alpha`.  Triangle
inequality then gives recurrence

```math
e_j\le\eta_j+(1-\alpha_j)e_{j-1},
```

whose iteration is CR.19.  The uniform bound `eta/alpha` in CR.20 has the
correct indexing and geometric-series constant.

## Scope

The suite proves three different claims and does not conflate them:

1. exact semantic convolution is nonexpansive;
2. independent occurrence sampling is reusable along a fixed tree;
3. semantic/DAG self-reuse needs fresh independence, finer information, or a
   contraction mechanism.

Neither the declared-tree compiler nor Doeblin forgetting applies to an
arbitrary adversarial dense interaction, and no such claim is made.
