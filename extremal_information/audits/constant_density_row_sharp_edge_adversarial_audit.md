# Adversarial audit: sharp edge for constant-density row fibres

**Frozen source:**
`extremal_information/drafts/constant_density_row_sharp_edge.md`

**SHA-256:**
`1d12a93158191635139175cd9818ccddeaaa93745dbc78290c4d908a0e758690`

**Verdict:** **PASS.**  The whitening, uniform Strong Tail Projection
verification, weak quadratic-form concentration, both imported random-matrix
applications, edge matching, matrix orientation, and pressure implication
are correct.  I found no counterexample and no missing substantive
hypothesis.

Two optional wording repairs would make the imported-theorem application
fully explicit:

1. In the expectation-to-probability argument, “probability-one lower bound”
   should read “lower bound holding with probability tending to one.”  The
   displayed argument itself is valid.
2. The cited theorems index the row law by its ambient dimension, whereas the
   source uses the pair `(d_r,r)` and `d_r` need not be strictly increasing.
   Add the standard sentence that one applies each theorem after arbitrary
   subsequence extraction and reindexing/filling unused dimensions.  Since
   `d_r -> infinity`, `d_r/r -> 1`, and all tail constants are uniform, this
   is routine and changes no conclusion.

I read the primary sources rather than relying on their abstracts:

- D. Chafaï and K. Tikhomirov,
  [*On the convergence of the extremal eigenvalues of empirical covariance
  matrices with dependence*](https://arxiv.org/abs/1509.02231), Definition
  1.2 and Theorem 1.7.
- P. Yaskov,
  [*The necessary and sufficient conditions in the Marchenko--Pastur
  theorem*](https://arxiv.org/abs/1511.02711), Assumption (A1), Definition
  (MP), and Theorem 3.3.

## 1. Covariance peel and whitening

Central symmetry of `E_r` gives `E R=0`.  The companion Parseval estimate

```math
\|\Sigma_r-I\|_F^2\le2(p_0^{-1}-1)
```

is correct: the ordered off-diagonal covariance sum is twice the sum of the
squared level-two Fourier coefficients and is bounded by twice the complete
nonconstant Parseval mass.

Every eigenvalue removed by `P_r` differs from one by more than
`delta_r=r^(-1/4)`, so

```math
k_r\delta_r^2\le\|\Sigma_r-I\|_F^2,
\qquad k_r=O_{p_0}(r^{1/2}),
\qquad d_r=r-k_r=r-o(r).
```

On `ran(V_r)`, `S_r` has spectrum in
`[1-delta_r,1+delta_r]` and is invertible for all sufficiently large `r`.
Because `V_r` is a spectral projection of `Sigma_r`, it commutes with that
matrix.  Consequently

```math
X_r=S_r^{-1/2}V_rR
```

is centered and exactly isotropic:

```math
\mathbb E X_rX_r^T
=S_r^{-1/2}(V_r\Sigma_rV_r)S_r^{-1/2}=I_{d_r}.
```

In an orthonormal coordinate system on `ran(V_r)`, a row of `B_rV_r` is
`X_r^TS_r^(1/2)`, hence

```math
B_rV_r=\mathbb X_rS_r^{1/2}.
```

This confirms the orientation in CE.8.

## 2. Uniform projection tails

Let `Q` be any rank-`k` orthogonal projection on the whitened space and

```math
M=V_rS_r^{-1/2}QS_r^{-1/2}V_r.
```

Then `M` is positive semidefinite, has rank at most `k`, and

```math
\|M\|_{op}\le(1-\delta_r)^{-1},
\qquad
\|M\|_F^2\le k(1-\delta_r)^{-2}.
```

Under the conditioned law,

```math
\mathbb E_{\mu_r}R^TMR
=\operatorname {tr}(\Sigma_rM)=\operatorname {tr}Q=k,
```

whereas the unconditioned Rademacher center is `tr M`.  Functional calculus
on `S_r` gives

```math
|\operatorname {tr}M-k|
=|\operatorname {tr}(Q(S_r^{-1}-I))|
\le {\delta_r\over1-\delta_r}k.
```

Since `k<=d_r<=r`, for all large `r`

```math
{\delta_r\over1-\delta_r}k\le2k^{3/4}.
```

This is the key uniformity calculation: it works simultaneously for rank
one, intermediate ranks, and ranks comparable to the ambient dimension.

Rademacher Hanson--Wright yields

```math
\Pr_U\{|W^TMW-\operatorname {tr}M|\ge u\}
\le C\exp[-c\min(u^2/k,u)],
```

with constants independent of `r`, `E_r`, and `Q`.  Conditioning costs only
`p_r^{-1}<=p_0^{-1}`.  For `t>=4k^(3/4)`, the centering discrepancy is at
most `t/2`, proving CE.13 after changing absolute constants.

This really is Chafaï--Tikhomirov STP, not merely a fixed-rank tail.  For
large `k`, take

```math
f(k)=4k^{-1/4},
\quad
g(k)={C\over p_0}
\sup_{t\ge4k^{3/4}}{t^2\over k}
e^{-c\min(t^2/k,t)}.
```

At the lower endpoint the exponent is `Omega(sqrt(k))`; beyond it the
displayed expression decreases up to harmless changes of constants.
Therefore `g(k)->0`.  For the finitely many small `k`, set `f(k)=1` and
enlarge `g(k)`: on `k<=t<4k^(3/4)` the trivial probability bound is
dominated by choosing `g(k)>=16sqrt(k)`, and CE.13 handles larger `t`.
Finitely many small ambient orders can be absorbed in the same values.
Thus the functions are uniform over every ambient dimension and every
projection, exactly as Definition 1.2 requires.

## 3. Chafaï--Tikhomirov Theorem 1.7

The primary theorem assumes centered isotropic vectors, iid sample rows,
STP uniform over dimensions, and sample-to-dimension ratios bounded away
from zero and infinity.  Here:

- centering and isotropy were checked above;
- the `r` rows of `mathbb X_r` are iid;
- Section 2 verifies STP uniformly;
- ambient dimension is `d_r` and sample size is `r`, with `d_r/r->1`.

The theorem therefore gives

```math
\limsup_r
{\mathbb E\lambda_{max}(\mathbb X_r^T\mathbb X_r)
 \over(\sqrt r+\sqrt{d_r})^2}\le1,
```

exactly CE.18.  The cited result is an expectation upper bound, so the source
correctly does not infer convergence from it alone.

The papers write a sequence `X_n in R^n`; our ambient dimensions `d_r` may
skip or repeat integers.  This is harmless but worth stating.  Given any
subsequence in `r`, extract a further subsequence on which `d_r` is strictly
increasing, place these laws at the corresponding ambient dimensions, and
fill unused dimensions with standard Rademacher laws and sample sizes of
ratio tending to one.  Uniform STP persists.  Theorem 1.7 applies to the
filled sequence.  Since every original subsequence has such a further
subsequence, CE.18 holds along the full sequence.

## 4. Yaskov (A1) and the MP law

For a deterministic positive semidefinite `d_r`-by-`d_r` matrix `A_r` with
`||A_r||op<=L`, use

```math
M=V_rS_r^{-1/2}A_rS_r^{-1/2}V_r.
```

Then

```math
\|M\|_{op}\le {L\over1-\delta_r},
\qquad
\|M\|_F^2\le {L^2d_r\over(1-\delta_r)^2},
```

and

```math
|\operatorname {tr}M-\operatorname {tr}A_r|
\le {L\delta_r\over1-\delta_r}d_r=o(d_r).
```

Conditioned Hanson--Wright at deviation `epsilon d_r`, after subtracting
this deterministic centering error, has probability at most
`C p_0^{-1}e^{-c_{epsilon,L}d_r}`.  Thus

```math
{X_r^TA_rX_r-\operatorname {tr}A_r\over d_r}
\longrightarrow0
```

in probability for every such matrix sequence.  This is precisely Yaskov's
Assumption (A1); positivity and the normalization by ambient dimension both
match the primary source.

Yaskov Theorem 3.3 says that (A1), under `p/n->rho`, is sufficient for the
MP property of the sample covariance.  Apply it to the transpose of
`mathbb X_r`: in Yaskov's notation `p=d_r`, `n=r`, and the iid vectors are
the columns of the `d_r`-by-`r` matrix `mathbb X_r^T`.  Therefore

```math
\mu_{r^{-1}\mathbb X_r^T\mathbb X_r}
\Longrightarrow \operatorname {MP}_1.
```

The same subsequence/reindexing argument as above resolves the cosmetic
`d_r` indexing issue.  No coordinate independence is assumed by Yaskov's
sufficient direction.

## 5. Matching the upper edge

The parameter-one MP law has support `[0,4]` and positive mass in every
interval `(4-eta,4)`.  Weak empirical-law convergence therefore implies

```math
\Pr\{\lambda_{max}(\mathbb X_r^T\mathbb X_r)/r\ge4-\eta\}
\longrightarrow1
```

for every fixed `eta>0`.  This is CE.19.  It is a lower bound on the largest
eigenvalue obtained from mass near the MP *upper* support edge; it is not a
claim that MP convergence itself excludes outliers.

Let `Z_r=lambda_max/r`.  CE.18 and `d_r/r->1` give
`limsup E Z_r<=4`.  If along a subsequence
`Pr(Z_r>=4+epsilon)>=eta`, then for any `delta>0`,

```math
\mathbb E Z_r
\ge(4-\delta)\Pr(Z_r\ge4-\delta)
 +(\epsilon+\delta)\Pr(Z_r\ge4+\epsilon).
```

Taking the subsequential limit and then choosing
`delta<eta epsilon/2` forces `liminf E Z_r>4`, a contradiction.  Hence the
upper tail also vanishes and `Z_r->4` in probability.  This verifies the
nontrivial expectation-plus-MP inference in the source.  Its phrase
“probability-one lower bound” should only be read as “with probability
tending to one.”

Taking square roots yields `||mathbb X_r||op/sqrt(r)->2`.

## 6. Reverse inequality and matrix orientation

The upper half of CE.4 follows immediately from

```math
\|\mathbb X_rS_r^{1/2}\|_{op}
\le\sqrt{1+\delta_r}\|\mathbb X_r\|_{op}.
```

The reverse inequality in the source is also correct, although operator
norm is not generally monotone under arbitrary right multiplication.  Since
`S_r^(1/2)` is invertible,

```math
\|\mathbb X_r\|_{op}
=\|(\mathbb X_rS_r^{1/2})S_r^{-1/2}\|_{op}
\le {1\over\sqrt{1-\delta_r}}
 \|\mathbb X_rS_r^{1/2}\|_{op}.
```

Therefore

```math
\|B_rV_r\|_{op}
=\|\mathbb X_rS_r^{1/2}\|_{op}
\ge\sqrt{1-\delta_r}\|\mathbb X_r\|_{op},
```

exactly as claimed.

## 7. Pressure consequence and counterexample search

The companion reduction supplies a deterministic projection of rank
`O(sqrt r)`, expected removed nuclear cost `O(r)`, and a coupled projected
Frobenius cost `O(r^(3/4))`.  CE.5 supplies the remaining regular-bulk edge.
Intersect its probability-one-asymptotic event with the standard iid edge
event to obtain the `G_r` required by the repaired projected-coupling
criterion.  That criterion only requires `Pr(G_r^c)=o(1)`, so no unproved
tail rate is being used.  CE.21 follows.

I tried the natural hostile constant-density events: a parity fibre, a union
of distant Hamming layers, a tail band producing a large covariance spike,
and an event aligned with a fixed quadratic form.  Parity and high-order
dependence are already covered by conditioning Hanson--Wright; covariance
spikes are removed by `P_r`; and a quadratic alignment still costs only the
fixed likelihood factor `p_0^{-1}` in every quadratic tail.  None escapes
STP after whitening.

## 8. Required repairs

No mathematical repair is required.  The two optional precision edits are:

- replace “probability-one lower bound” by “lower bound with probability
  tending to one” in Section 4;
- add one sentence explaining subsequence/reindexing from `(d_r,r)` to the
  ambient-dimension indexing used in the cited theorems.

The frozen theorem remains valid without these editorial additions.

## 9. Repair and Corollary CE.2 verification

The strengthened source with SHA-256
`e12b1550da75ac2e9e9cc54c19c3f9b0b421ad3298583f53d9264a7b1cafa3dc`
implements both optional precision edits above and adds Corollary CE.2.
I independently checked the extension.  **Verdict for the strengthened
source: PASS.**

Let `g_r=dmu_r/dU_r` satisfy `E_U g_r=1` and `E_U g_r^2<=K`.  Parseval gives

```math
\sum_{S\ne\varnothing}\langle g_r,\chi_S\rangle_U^2
=\mathbb E_Ug_r^2-1\le K-1.
```

The ordered off-diagonal entries of `Sigma_r-I` use each level-two
coefficient twice, proving CE.22 with the factor `2(K-1)`.  Thus the same
`delta_r=r^(-1/4)` spectral peel has `k_r=O_K(sqrt r)`.

The KL step is exact:

```math
D(\mu_r\|U_r)=\mathbb E_{\mu_r}\log g_r
\le\log\mathbb E_{\mu_r}g_r
=\log\mathbb E_Ug_r^2\le\log K,
```

where Jensen is under `mu_r`.  Hence the companion entropy-transport
argument still couples one row to a uniform Rademacher row with expected
Hamming cost `O_K(sqrt r)`; independent rowwise copies give projected
Frobenius cost `O_K(r^(3/4))`.

For every event `A`, Cauchy--Schwarz gives

```math
\mu_r(A)=\mathbb E_U[g_r1_A]
\le\sqrt K\,U_r(A)^{1/2}.
```

Applied after the same whitening and centering calculation, this changes a
Rademacher Hanson--Wright tail `Ce^{-ca}` into at most
`sqrt(KC)e^{-ca/2}`.  The exponent remains uniform over all projection
ranks and dimensions.  Therefore the STP functions in CE.14--CE.16 still
tend to zero, and a fixed `epsilon d_r` quadratic-form deviation still has
probability `e^{-Omega_{K,epsilon}(d_r)}`.  Both the
Chafaï--Tikhomirov and Yaskov hypotheses remain valid.

Finally,

```math
\begin{aligned}
\mathbb E\|B_rP_r\|_*
&\le\sqrt{k_r}\,(\mathbb E\|B_rP_r\|_F^2)^{1/2}\\
&=\sqrt{k_r}\,(r\operatorname {tr}(P_r\Sigma_r))^{1/2}\\
&\le r\sqrt{k_r}=O_K(r^{5/4})=o(r^{3/2}),
\end{aligned}
```

because `Sigma_r` is positive semidefinite with
`tr Sigma_r=E||R||_2^2=r`.  The iid removed component remains smaller,
`O_K(r)`.  Thus every hypothesis of the repaired projected-coupling theorem
is supplied.

The pre-existing CE.1 proof remains unchanged mathematically.  Its added
subsequence/reindexing sentence and corrected probability wording are
accurate.  No further repair is required.
