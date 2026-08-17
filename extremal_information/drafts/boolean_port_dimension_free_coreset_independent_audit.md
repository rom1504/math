# Independent audit: dimension-free Boolean-port response coresets

**Verdict: PASS.**  The symmetrization and contraction constants in RC.1
are safe, the projective conventions are consistent, the cover count is
correct with natural logarithms, and the packing argument proves the stated
fixed-`eta` linear law.

## 1. Empirical-process constants

For

```math
f_\epsilon(s)=|s\cdot\epsilon|/p,
```

standard ghost-sample symmetrization gives

```math
\mathbb E\sup_\epsilon|(\mu-\mu_k)f_\epsilon|
\le {2\over k}\mathbb E\sup_\epsilon
 \left|\sum_{j=1}^k\sigma_jf_\epsilon(S_j)\right|.
```

Apply the absolute-supremum contraction lemma to
`phi(t)=|t|` and `t=(s dot epsilon)/p`.  Its safe form has factor two, so the
coefficient becomes `4/k`, exactly as in RC.8.  Conditional on the sample,

```math
\sup_\epsilon\left|\sum_j\sigma_j{S_j\cdot\epsilon\over p}\right|
={1\over p}\sum_{i=1}^p\left|\sum_j\sigma_j(S_j)_i\right|.
```

Each inner Rademacher sum has second moment `k`; hence its expected absolute
value is at most `sqrt(k)`.  Therefore

```math
\mathbb E d_p(\mu,\mu_k)\le4/\sqrt{k}.
```

Choosing `k=ceil(16/eta^2)` really does make this at most `eta`.  No endpoint
union bound or hidden dependence on `p` remains.

## 2. Quotient and cover conventions

Both `s -> -s` and `epsilon -> -epsilon` leave the response unchanged, so
sampling canonical representatives of the projective cube loses no query or
measure information.  There are `2^(p-1)` possible sampled atoms.  Enumerating
all ordered empirical samples gives at most

```math
2^{(p-1)k}
```

centres, and hence, when `log` is the natural logarithm,

```math
\log\operatorname {Cov}_p(\eta)
\le(p-1)k\log2.
```

This is RC.5.  Ordered samples overcount empirical measures, which is harmless
for an upper bound.

## 3. Packing-radius check

For point masses, if `h(s,t)` is projective Hamming distance, the exact metric
identity is

```math
d_p(\delta_s,\delta_t)=2h(s,t)/p.
```

The reverse triangle inequality gives the upper bound, and querying at
`epsilon=s` gives equality.  For fixed `eta<theta<1/2`, a greedy code in the
projective cube with distance at least `theta p` has

```math
2^{(1-H_2(\theta)-o(1))p}
```

words: a projective ball below the antipodal threshold has the usual binary
entropy exponent.  Its response points are separated by at least `2theta`,
which is strictly greater than `2eta`; consequently radius-`eta` cover balls
contain at most one codeword.  This supplies the positive linear lower
exponent required for RC.6.

## 4. Scope

RC.1 is a static existence theorem.  The empirical coreset is allowed to
depend on the whole input measure, and repeated resampling can accumulate
error.  Thus the draft correctly refrains from claiming an all-depth merge
scheme or a dynamic congruence.

The finite verifier also passes:

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_boolean_port_dimension_free_coreset.py
```

It checks the conditional contraction and `ell_1` identities exactly on
5,054 small sample matrices.  Those checks are diagnostic; the argument
above proves the general statement.
