# Independent-style self-audit: dimension-free Boolean-port coresets

**Object audited.** `boolean_port_dimension_free_coreset.md`, Theorem RC.1.

**Verdict.** PASS, with the theorem scoped to one-shot response approximation
and the fixed-error entropy conclusion scoped to `0<eta<1/2`.

## Constant audit

For a fixed sample, symmetrization contributes the first factor `2`:

```math
\mathbb E\sup_f|(P_k-P)f|
\le2\mathbb E\sup_f|k^{-1}\sum_j\sigma_jf(S_j)|.
```

The safe absolute-value form of scalar Rademacher contraction contributes a
second factor `2`.  For the underlying linear class, exact optimization over
the Boolean endpoint gives

```math
\mathbb E_\sigma\sup_\epsilon
\left|{1\over kp}\sum_j\sigma_jS_j\cdot\epsilon\right|
={1\over kp}\sum_{i=1}^p
 \mathbb E_\sigma\left|\sum_j\sigma_j(S_j)_i\right|
\le {1\over\sqrt k}.
```

Therefore the displayed `4/sqrt(k)` bound is valid.  Taking
`k=ceil(16/eta^2)` indeed makes this at most `eta`; no hidden dependence on
the `2^(p-1)` queries remains.

The proof is unchanged for generators in `[-1,1]^p`, because
`sum_j(S_j)_i^2<=k` is all that the last estimate uses.

## Counting and metric audit

There are `|G_p|=2^(p-1)` projective row types.  Ordered empirical samples
therefore give at most `2^{(p-1)k}` response centres.  With **natural**
logarithm this is exactly

```math
\log \operatorname {Cov}_p(\eta)
\le k(p-1)\log2.
```

Repeated sample values and response collisions only reduce this count.

For the lower bound one must choose `eta<theta<1/2`, not `theta=eta`.
The projective code points have response distance at least `2theta>2eta`;
hence an arbitrary radius-`eta` ball contains at most one code point.  The
greedy code exponent `1-H_2(theta)` is positive, yielding the claimed
`Omega_eta(p)` natural-log lower bound.

## Scope audit

RC.1 proves existence of a small empirical representative of each static
response roof.  It does **not** make independently chosen representatives a
congruence under repeated merges.  A fresh error can enter at every
resparsification.  Consequently the result closes the fixed-scale metric
entropy gap but does not by itself solve the dynamic-memory problem.

The finite verifier checks the two conditional identities for all sample
matrices in its stated range.  It is a normalization diagnostic, not a
replacement for the analytic contraction proof.
