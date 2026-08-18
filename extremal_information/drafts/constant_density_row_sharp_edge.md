# Sharp projected edge for every constant-density row fibre

**Status.** Task-local theorem report.  The only imported ingredients are
Theorem 1.7 of Chafaï--Tikhomirov and the Marchenko--Pastur sufficiency
theorem of Yaskov; their hypotheses are verified below.  Together with the
separate projected-coupling theorem, this closes the constant-density,
centrally symmetric row-product class throughout the strict conference
interval.

## 1. Setup and statement

Let `U_r` be uniform on `{+-1}^r`.  Let `E_r=-E_r` and

```math
p_r=U_r(E_r)\ge p_0>0,
\tag{CE.1}
```

and put `mu_r=U_r(.|E_r)`.  For `R~mu_r`, let

```math
\Sigma_r=\mathbb E RR^T.
\tag{CE.2}
```

Fix `delta_r=r^(-1/4)`, and let `P_r` be the spectral projection of
`Sigma_r` onto the eigenvalues outside
`[1-delta_r,1+delta_r]`.  Write `V_r=I-P_r`, and let `B_r` be the square
matrix whose `r` rows are independent copies of `R`.

The Fourier--Parseval argument in the companion reduction gives

```math
\|\Sigma_r-I\|_F^2\le 2(p_0^{-1}-1),
\qquad
k_r:=\operatorname{rank}P_r=O_{p_0}(r^{1/2}).
\tag{CE.3}
```

### Theorem CE.1 (sharp regular-bulk edge)

Under (CE.1),

```math
\boxed{
{\|B_rV_r\|_{op}\over\sqrt r}\longrightarrow 2
\quad\hbox{in probability}.}
\tag{CE.4}
```

In particular, the upper bound required by the projected-coupling theorem
holds uniformly over this class:

```math
\|B_rV_r\|_{op}\le(2+o_{\Pr}(1))\sqrt r.
\tag{CE.5}
```

## 2. Whitening the regular subspace

Identify the range of `V_r` with `R^{d_r}`, where
`d_r=r-k_r=r-o(r)`, and let

```math
S_r=V_r\Sigma_rV_r|_{\operatorname{ran}V_r}.
\tag{CE.6}
```

All eigenvalues of `S_r` lie in
`[1-delta_r,1+delta_r]`.  The random vector

```math
X_r=S_r^{-1/2}V_rR\in\mathbb R^{d_r}
\tag{CE.7}
```

is centered and isotropic.  Let `mathbb X_r` be the `r by d_r` matrix with
independent rows distributed as `X_r`.  In compatible coordinates,

```math
B_rV_r=\mathbb X_rS_r^{1/2},
\qquad
\|B_rV_r\|_{op}\le\sqrt{1+\delta_r}\,\|\mathbb X_r\|_{op}.
\tag{CE.8}
```

It remains to prove the sharp square-edge law for `mathbb X_r`.

## 3. Uniform quadratic-form concentration

Let `Q` be any orthogonal projection on `R^{d_r}` of rank `k>=1`.  In the
original coordinates,

```math
\|QX_r\|_2^2=R^TMR,
\qquad
M=V_rS_r^{-1/2}QS_r^{-1/2}V_r.
\tag{CE.9}
```

The matrix `M` is positive semidefinite and

```math
\|M\|_{op}\le(1-\delta_r)^{-1},
\quad
\|M\|_F^2\le k(1-\delta_r)^{-2},
\quad
|\operatorname{tr}M-k|
\le {\delta_r\over1-\delta_r}k.
\tag{CE.10}
```

The last identity uses `E_mu R^TMR=tr Q=k`, whereas `tr M` is the
expectation under the unconditioned Rademacher law.

For an independent Rademacher vector `W`, Hanson--Wright gives absolute
constants `c,C>0` such that, for all `u>=0`,

```math
U_r\{|W^TMW-\operatorname{tr}M|\ge u\}
\le C\exp\{-c\min(u^2/k,u)\},
\tag{CE.11}
```

after harmless adjustment of the constants for the finitely many small
orders.  Conditioning costs at most `p_0^{-1}`.

Because `k<=d_r<=r`,

```math
{\delta_r\over1-\delta_r}k\le2k^{3/4}
\tag{CE.12}
```

for all large `r`.  Hence, whenever `t>=4k^{3/4}`, (CE.11) yields

```math
\Pr\{\|QX_r\|_2^2-k\ge t\}
\le {C\over p_0}\exp\{-c\min(t^2/k,t)\}.
\tag{CE.13}
```

This verifies the **Strong Tail Projection** property of
Chafaï--Tikhomirov.  Indeed, take

```math
f(k)=\min(1,4k^{-1/4})
\tag{CE.14}
```

and enlarge a finite number of values of `g(k)` so that, for all dimensions,
all rank-`k` projections, and all `t>=f(k)k`,

```math
\Pr\{\|QX_r\|_2^2-k\ge t\}
\le {g(k)k\over t^2}.
\tag{CE.15}
```

For large `k`, (CE.13) permits

```math
g(k)={C\over p_0}\sup_{t\ge4k^{3/4}}
{t^2\over k}\exp\{-c\min(t^2/k,t)\}=o(1).
\tag{CE.16}
```

Thus both defining functions tend to zero as required.  Notice that this
checks the projection-tail hypothesis itself; coordinate independence of
`X_r` is neither asserted nor used.

The same computation verifies Yaskov's weak quadratic-form condition.  If
`A_r` is any positive semidefinite `d_r by d_r` matrix with uniformly
bounded operator norm, apply (CE.11) to
`M=V_rS_r^{-1/2}A_rS_r^{-1/2}V_r`.  Then

```math
{X_r^TA_rX_r-\operatorname{tr}A_r\over d_r}
\longrightarrow0
\quad\hbox{in probability},
\tag{CE.17}
```

because the uniform-law centering error is
`O(delta_r d_r)=o(d_r)` and a fixed `epsilon d_r` deviation has probability
`O_{p_0}(exp(-c_epsilon d_r))`.

## 4. Matching upper and lower edges

Theorem 1.7 of Chafaï--Tikhomirov applies to the independent isotropic rows
`X_r`, since (CE.15) is their STP property and `d_r/r -> 1`.  Formally, if
an asserted limit failed, pass to a subsequence on which the integer
dimensions `d_r` are strictly increasing and fill any skipped dimensions by
standard Rademacher vectors; the uniform STP functions in (CE.14)--(CE.16)
are unchanged.  It gives

```math
\limsup_{r\to\infty}
{\mathbb E\lambda_{\max}(\mathbb X_r^T\mathbb X_r)
 \over(\sqrt r+\sqrt{d_r})^2}\le1.
\tag{CE.18}
```

Yaskov's Theorem 3.3 (the Marchenko--Pastur theorem) applies by (CE.17), again with aspect
ratio `d_r/r -> 1`.  Therefore the empirical spectral law of
`r^{-1}mathbb X_r^Tmathbb X_r` converges to the parameter-one
Marchenko--Pastur law.  Since that law has positive mass in every interval
immediately below its upper edge four,

```math
\liminf_{r\to\infty}
{\lambda_{\max}(\mathbb X_r^T\mathbb X_r)\over r}\ge4
\quad\hbox{in probability}.
\tag{CE.19}
```

Equations (CE.18)--(CE.19) imply convergence in probability to four.  For
completeness, if an upper deviation by `epsilon` retained probability
`eta>0`, combine it with the lower bound `4-delta` holding with probability
tending to one, with
`delta<eta epsilon/2`, to force the expectation above `4+o(1)`, contrary
to (CE.18).  Thus

```math
{\|\mathbb X_r\|_{op}\over\sqrt r}\longrightarrow2
\quad\hbox{in probability}.
\tag{CE.20}
```

Equation (CE.8) proves the upper half of (CE.4); its lower half follows
from (CE.20) and the reverse inequality
`||B_rV_r||_op >= sqrt(1-delta_r)||mathbb X_r||_op`.  This proves Theorem
CE.1.  `square`

### Corollary CE.2 (bounded Renyi-two row laws)

The event-conditioning hypothesis can be replaced by the following strictly
more general assumption.  Let `mu_r` be any centrally symmetric law on the
row cube, write `g_r=dmu_r/dU_r`, and suppose

```math
\boxed{\mathbb E_{U_r}g_r^2\le K<\infty}
\tag{CE.21}
```

uniformly in `r`.  Then Theorem CE.1 and its pressure consequence remain
valid, with constants depending only on `K`.

**Proof.**  Parseval now gives

```math
\|\Sigma_r-I\|_F^2\le2(K-1),
\tag{CE.22}
```

so the same covariance peel has rank `O_K(sqrt r)`.  Renyi monotonicity (or
Jensen under `mu_r`) gives

```math
D(\mu_r\|U_r)\le\log\mathbb E_Ug_r^2\le\log K,
\tag{CE.23}
```

which is exactly what the Hamming transport proof uses.  Finally,
Cauchy--Schwarz replaces the bounded-likelihood-ratio step:

```math
\mu_r(A)=\mathbb E_U[g_r1_A]
\le\sqrt K\,U_r(A)^{1/2}.
\tag{CE.24}
```

Applying (CE.24) to the Hanson--Wright event only halves its exponential
constant.  Equations (CE.13)--(CE.17), including STP and (A1), therefore
continue to hold.  For the removed component one no longer needs the sharper
constant-density estimate: since `tr(P_rSigma_r)<=r`,

```math
\mathbb E\|B_rP_r\|_*
\le\sqrt{k_r}\,(r\operatorname{tr}(P_r\Sigma_r))^{1/2}
\le r\sqrt{k_r}=O_K(r^{5/4})=o(r^{3/2}).
```

Thus the projected-coupling hypotheses still all hold.  `square`

## 5. Pressure consequence

The companion reduction already proves, uniformly over (CE.1):

1. `rank(P_r)=O_{p_0}(sqrt r)`;
2. the removed nuclear costs are `O_{p_0}(r)=o(r^(3/2))`;
3. there is a rowwise coupling to an iid Rademacher bridge whose projected
   Frobenius cost is `O_{p_0}(r^(3/4))=o(r)`.

Theorem CE.1 supplies the last hypothesis of the projected-coupling
criterion.  Consequently, for every fixed `0<beta<sqrt(2)/6`, every
centrally symmetric row-product bridge law satisfying either (CE.1) or the
bounded-Renyi-two condition (CE.21) satisfies

```math
\boxed{
\mathbb E[(h_\beta-f(B_r)/r)_+]\longrightarrow0.}
\tag{CE.25}
```

Thus no such law can lower conference pressure by a positive amount per
spin.  A favorable speed-`r` law, if one exists, must use non-product
dependence between rows, row events whose mass vanishes with `r`, or a
mechanism outside this bridge-law model.

## 6. Imported results and hypothesis map

- D. Chafaï and K. Tikhomirov, *On the convergence of the extremal
  eigenvalues of empirical covariance matrices with dependence*,
  arXiv:1509.02231, Definition 1.2 and Theorem 1.7.  Definition 1.2 is
  exactly (CE.15); Theorem 1.7 is exactly the expectation bound (CE.18).
- P. Yaskov, *The necessary and sufficient conditions in the
  Marchenko--Pastur theorem*, arXiv:1511.02711, Assumption (A1) and the main
  MP theorem.  Assumption (A1) is exactly (CE.17); its sufficient direction
  gives the empirical-law statement used in (CE.19).

The first paper controls the upper edge but does not by itself supply the
matching lower edge.  The second controls the empirical law but not the
upper edge.  Both are needed, and neither theorem is being applied merely
from a bounded-subgaussian slogan.

## 7. Frontier movement

The sharp-edge lemma left open in the constant-density reduction is now a
theorem.  This closes an entire class rather than one named fibre.  It also
identifies the next irreducible question: whether `Theta(r)` total bridge
information can create a pressure gain through **cross-row dependence**
while every individual row remains spectrally regular after its low-rank
covariance peel.
