# Dimension-free sampling for Boolean-port response roofs

**Status.** Rigorous theorem.  This closes the fixed-distortion gap in
Theorem PF.2 of `boolean_port_fourier_feature_algebra.md`: the logarithmic
covering number is linear, not merely quadratic, in the number of ports.
It is a static response theorem and makes no all-depth mergeability claim.

## 1. Setting

Let

```math
G_p=\{+-1\}^p/\{s\sim-s\}
```

and, for a probability measure `mu` on `G_p`, put

```math
R_\mu(\epsilon)
=\mathbb E_{s\sim\mu}{|s\cdot\epsilon|\over p},
\qquad \epsilon\in G_p.                            \tag{RC.1}
```

The normalized contextual metric is

```math
d_p(\mu,\nu)=\|R_\mu-R_\nu\|_\infty.              \tag{RC.2}
```

Write `Cov_p(eta)` for the smallest number of response functions needed to
cover `{R_mu:mu in P(G_p)}` in this metric at radius `eta`.  Every unlabelled
`log` below is the natural logarithm.

## 2. Dimension-free empirical coreset

### Theorem RC.1

For every `p>=1`, every probability measure `mu` on `G_p`, and every
`0<eta<1`, there is an equally weighted empirical measure `nu` supported on
at most

```math
k=\left\lceil {16\over\eta^2}\right\rceil          \tag{RC.3}
```

projective sign vectors such that

```math
d_p(\mu,\nu)\le\eta.                               \tag{RC.4}
```

Consequently

```math
\boxed{
\log \operatorname {Cov}_p(\eta)
\le (p-1)\log2\left\lceil {16\over\eta^2}\right\rceil .}
                                                               \tag{RC.5}
```

In particular, for every fixed `0<eta<1/2`, the point-mass packing from
PF.2 and (RC.5) give

```math
\boxed{\log\operatorname {Cov}_p(\eta)=\Theta_\eta(p).}       \tag{RC.6}
```

Thus fixed-error Boolean-port response compression costs a linear number of
bits in the port count.  The earlier `O(p^2/eta^2)` estimate was an artifact
of applying a union bound separately to `2^(p-1)` endpoint queries.

### Proof

Choose representatives in `{+-1}^p`.  For `epsilon in {+-1}^p`, define

```math
f_\epsilon(s)={|s\cdot\epsilon|\over p},
\qquad
g_\epsilon(s)={s\cdot\epsilon\over p}.             \tag{RC.7}
```

Let `S_1,...,S_k` be independent samples from `mu` and let `mu_k` be their
empirical measure.  Symmetrization followed by the scalar Rademacher
contraction principle for the map `t mapsto |t|` gives

```math
\begin{aligned}
\mathbb E_S d_p(\mu,\mu_k)
&\le {2\over k}\mathbb E_{S,\sigma}
   \sup_\epsilon\left|\sum_{j=1}^k
      \sigma_j f_\epsilon(S_j)\right|\\
&\le {4\over k}\mathbb E_{S,\sigma}
   \sup_\epsilon\left|\sum_{j=1}^k
      \sigma_j g_\epsilon(S_j)\right|.             \tag{RC.8}
\end{aligned}
```

We use the safe factor-two absolute-value form of contraction: for any
`T subset R^k` and coordinatewise contractions `phi_j` with `phi_j(0)=0`,

```math
\mathbb E_\sigma\sup_{t\in T}
 \left|\sum_j\sigma_j\phi_j(t_j)\right|
\le2\mathbb E_\sigma\sup_{t\in T}
 \left|\sum_j\sigma_jt_j\right|.                  \tag{RC.9}
```

This is the elementary Rademacher contraction lemma (apply the usual
one-coordinate contraction induction to `T union (-T)`).  No cardinality
of the query class enters it.

For fixed sampled rows, set

```math
V_i=\sum_{j=1}^k\sigma_j(S_j)_i.
```

The remaining supremum is explicit:

```math
\sup_{\epsilon\in\{+-1\}^p}
 \left|\sum_j\sigma_jg_\epsilon(S_j)\right|
={1\over p}\sup_\epsilon\left|\sum_i\epsilon_iV_i\right|
={1\over p}\sum_i|V_i|.                           \tag{RC.10}
```

By Cauchy--Schwarz in the Rademacher variables,

```math
\mathbb E_\sigma|V_i|
\le\sqrt{\mathbb E_\sigma V_i^2}=\sqrt{k}.         \tag{RC.11}
```

Equations (RC.8)--(RC.11) therefore yield the dimension-free estimate

```math
\mathbb E_S d_p(\mu,\mu_k)\le {4\over\sqrt{k}}.    \tag{RC.12}
```

For `k` in (RC.3), some sample has discrepancy at most `eta`.  There are at
most `|G_p|^k=2^{(p-1)k}` ordered empirical samples, so their response
functions form the cover (RC.5).

For the lower bound, fix **strictly** `eta<theta<1/2`.  A projective binary code of
relative Hamming distance at least `theta` has
`2^{(1-H_2(theta)-o(1))p}` words, while PF.2 gives exact pairwise response
distance at least `2theta`.  Since no radius-`eta` ball can contain two
points at mutual distance greater than `2eta`, every radius-`eta` cover has
at least this many members, proving the lower half of (RC.6). `square`

## 3. General zonoid form

The proof did not use that sampled generators are exactly Boolean, only
that every coordinate lies in `[-1,1]`.  If `mu` is a probability measure
on `[-1,1]^p` and

```math
R_\mu(\epsilon)=\mathbb E_\mu|s\cdot\epsilon|/p,
```

then the same empirical support bound `ceil(16/eta^2)` holds for the
restriction of its zonoid support function to Boolean directions.  Thus
RC.1 is a uniform coreset theorem for this restricted support-function
query class, rather than a peculiarity of the finite projective group.

## 4. What this does and does not compress

The exact Fourier state still has `2^(p-1)-1` independent coordinates.
RC.1 says that its **fixed-scale response image** has only `Theta_eta(p)`
bits of metric entropy.  It does not give a reusable bounded-error merge
scheme: replacing a measure by a fresh empirical coreset after every
composition can still add error.  Dynamic congruence remains a separate
resource.

## 5. Finite diagnostic

The companion verifier exhausts all projective sample matrices for
`p<=4`, `k<=4`.  In exact integer arithmetic it checks the conditional
contraction inequality used in (RC.8) and the explicit `ell_1` identity
(RC.10):

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_boolean_port_dimension_free_coreset.py
```

The finite check is a normalization diagnostic; the proof above is the
general certificate.
