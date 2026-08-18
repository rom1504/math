# Independent external-theorem audit: constant-density projected sharp edge

**Frozen source:**
`extremal_information/drafts/constant_density_row_sharp_edge.md`

**SHA-256:**
`1d12a93158191635139175cd9818ccddeaaa93745dbc78290c4d908a0e758690`

**Verdict:** **PASS.**  The Chafaï--Tikhomirov STP hypothesis and expectation
edge theorem, the Yaskov quadratic-form hypothesis and Marchenko--Pastur
theorem, all dimension/aspect normalizations, and the final
expectation-plus-MP implication are used correctly.  No mathematical repair
is required.

The audit reads the actual primary sources:

- D. Chafaï and K. Tikhomirov, [*On the convergence of the extremal
  eigenvalues of empirical covariance matrices with
  dependence*](https://arxiv.org/abs/1509.02231), Definition 1.2 and
  Theorem 1.7.
- P. Yaskov, [*The necessary and sufficient conditions in the
  Marchenko--Pastur theorem*](https://arxiv.org/abs/1511.02711), assumption
  (A1), the definition of (MP), and Theorem 3.3.

For reproducibility, the arXiv v1 PDFs inspected locally had SHA-256 hashes

```text
1509.02231  67f851a4299465fc4d44eba6cca610f7b35f856dae06d2235f4bb590ed6a5dfc
1511.02711  df932f14cc0c49d532e96e037f64ed7b0ef48d8392fffc54cd496cf771c9bb2e
```

## 1. Chafaï--Tikhomirov: exact STP quantifiers

Definition 1.2 concerns a sequence of centered isotropic vectors, one in
each ambient dimension.  It requires two functions `f:N->[0,1]` and
`g:N->R_+` which both tend to zero with the *projection rank*.  For every
ambient dimension, every nonzero orthogonal projection `Q`, and every

```math
t\ge f(k)k,\qquad k=\operatorname{rank}Q,
```

the required one-sided estimate is

```math
P\{\|QX\|_2^2-k\ge t\}\le {g(k)k\over t^2}.
\tag{EA.1}
```

Thus the definition requires uniformity in ambient dimension and in the
projection, but only an upper tail.  CE.13--CE.16 establish precisely this
property.

After whitening, `X_r` is centered because `E_r=-E_r`, and it is exactly
isotropic on the `d_r`-dimensional regular subspace.  For a rank-`k`
projection, Hanson--Wright after conditioning gives, uniformly in `r,Q`,

```math
P\{\|QX_r\|_2^2-k\ge t\}
\le {C\over p_0}\exp\{-c\min(t^2/k,t)\}
\qquad(t\ge4k^{3/4}).
\tag{EA.2}
```

The trace recentering used here is correct.  Since
`delta_r=r^(-1/4)`, `k<=d_r<=r`, and the regular covariance eigenvalues lie
in `[1-delta_r,1+delta_r]`,

```math
|\operatorname{tr}M-k|
\le {\delta_r\over1-\delta_r}k
\le2k^{3/4}.
```

Thus `t>=4k^(3/4)` leaves at least `t/2` for Hanson--Wright.

Taking

```math
f(k)=\min(1,4k^{-1/4})
```

is valid.  For all sufficiently large `k`, multiplying (EA.2) by
`t^2/k` and taking the supremum over `t>=4k^(3/4)` gives a function tending
to zero.  The finitely many smaller ranks can be assigned finite larger
values of `g(k)`: on the bounded interval before (EA.2) applies use the
trivial probability bound, and afterward use (EA.2).  These values do not
affect `g(k)->0`.  If the monotonicity used without loss of generality in
the paper's proof is desired, replace `g(k)` by its decreasing tail
envelope; Definition 1.2 itself does not require monotonicity.

This verifies the full STP property, not merely a projection estimate for
ranks growing proportionally to the ambient dimension.

## 2. Chafaï--Tikhomirov: dimensions, aspect, and conclusion

In the primary paper, `X_n in R^n`, `m_n` independent samples, and

```math
A_n=\sum_{j=1}^{m_n}X_n^{(j)}(X_n^{(j)})^T.
```

Theorem 1.7 assumes

```math
0<\liminf n/m_n\le\limsup n/m_n<\infty
```

and STP, and concludes

```math
\limsup {E\lambda_{max}(A_n)\over(\sqrt{m_n}+\sqrt n)^2}\le1.
\tag{EA.3}
```

Here the ambient dimension is `d_r`, the number of independent samples is
`r`, and the sample sum is exactly
`mathbb X_r^T mathbb X_r`.  Since `d_r/r->1`, both aspect hypotheses hold,
and (EA.3) is exactly CE.18.  There is no missing factor `1/r`: the theorem
is stated for the unnormalized sample sum, while Yaskov below uses the
normalized covariance.

The paper writes the ambient dimension itself as the sequence index.  The
present triangular notation with dimension `d_r` is harmless.  Formally,
if a claimed conclusion failed, choose a violating subsequence and then a
further subsequence on which `d_r` is strictly increasing (possible because
`d_r->infinity`).  Index Chafaï--Tikhomirov by those ambient dimensions and
fill missing dimensions with standard isotropic Rademacher vectors, using
`m_n=n` there.  A common STP envelope works for both families and the aspect
ratio remains one.  The theorem then applies to the violating subsequence.
This also proves the claimed uniformity over all admissible row fibres by
the same worst-case subsequence argument.

The conclusion imported from Theorem 1.7 is only an expectation upper
edge.  The source correctly does not claim that it alone gives convergence
in probability.

## 3. Yaskov: exact (A1) and aspect mapping

Yaskov's assumption (A1) says that for every sequence of real symmetric
positive-semidefinite matrices `A_p` with uniformly bounded operator norm,

```math
{x_p^TA_px_p-\operatorname{tr}A_p\over p}
\longrightarrow0
\quad\text{in probability}.
\tag{EA.4}
```

This is exactly the condition checked in CE.17, with `p=d_r`.  For the
pulled-back Rademacher quadratic form, the operator norm is uniformly
bounded and the Frobenius norm squared is `O(d_r)`.  Hanson--Wright and
conditioning by an event of probability at least `p_0` therefore make a
fixed `epsilon d_r` deviation `exp(-Omega_epsilon(d_r))`.  The difference
between the unconditioned trace center and `tr A_r` is
`O(delta_r d_r)=o(d_r)`.  This verifies (A1) for every allowed matrix
sequence, not only for projections.

Yaskov uses a `p by n` data matrix whose columns are independent copies of
the isotropic `p`-vector and assumes `p/n->rho>0`.  Transposing the source's
`r by d_r` row matrix gives

```math
p=d_r,\qquad n=r,\qquad p/n\to1.
```

Thus its normalized covariance is exactly

```math
{1\over r}\mathbb X_r^T\mathbb X_r.
```

Theorem 3.3 says that (A1) implies the paper's property (MP), which is
almost-sure weak convergence for every allowed orthonormal row compression.
Taking the identity compression gives almost-sure convergence of this
empirical spectral law to the parameter-one MP law.  The source only uses
the weaker in-probability consequence.  No almost-sure edge convergence is
being imported from Yaskov.

As for Chafaï--Tikhomirov, a violating-subsequence/fill-in argument removes
the notational difference between ambient dimensions `d_r` and a sequence
indexed by every integer.  No cross-order independence is needed for the
in-probability conclusion used here.

## 4. MP bulk plus expectation edge really gives operator convergence

The parameter-one MP law has positive mass in every interval immediately
below its upper support edge `4`.  Weak empirical-law convergence therefore
implies, for every `delta>0`,

```math
P\left\{{\lambda_{max}(\mathbb X_r^T\mathbb X_r)\over r}
        \ge4-\delta\right\}\longrightarrow1.
\tag{EA.5}
```

This is only a lower edge bound for the largest eigenvalue; MP convergence
does not exclude outliers above four.

On the other hand, CE.18 and `d_r/r->1` give

```math
\limsup E{\lambda_{max}(\mathbb X_r^T\mathbb X_r)\over r}\le4.
\tag{EA.6}
```

Suppose an upper deviation `4+epsilon` retained probability at least
`eta>0` along a subsequence.  Combining it with (EA.5) gives

```math
\liminf E{\lambda_{max}\over r}
\ge(4+\epsilon)\eta+(4-\delta)(1-\eta)
>4
```

after taking, for example, `delta<eta epsilon/2`.  This contradicts
(EA.6).  Hence the normalized largest eigenvalue converges to four in
probability, and the largest singular value converges to two.

The phrase “probability-one lower bound” in the explanatory sentence after
CE.19 should be read as “a lower bound holding with probability tending to
one”; CE.19 itself states the correct mode.  This is a wording clarification,
not a proof defect.

Finally, because every eigenvalue of `S_r` lies in
`[1-delta_r,1+delta_r]`,

```math
\sqrt{1-\delta_r}\,\|\mathbb X_r\|_{op}
\le\|\mathbb X_rS_r^{1/2}\|_{op}
\le\sqrt{1+\delta_r}\,\|\mathbb X_r\|_{op}.
```

This transfers the limit two to `B_rV_r` and verifies both halves of CE.4.

## 5. Precise corrections

None required.  Two optional clarifications would make the source more
self-contained:

1. identify Yaskov's “main MP theorem” explicitly as Theorem 3.3 and state
   that its empirical-law conclusion is almost sure, although only
   convergence in probability is used;
2. replace “probability-one lower bound” in the expectation argument by
   “lower bound with probability tending to one.”

Neither changes the theorem or any quantitative conclusion.
