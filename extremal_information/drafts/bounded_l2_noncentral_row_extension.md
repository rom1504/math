# Bounded-L2 row laws: noncentral sharp edge and pressure no-gain

**Status.** Task-local theorem draft.  This removes central symmetry and the
hard-event form from the constant-density row-product theorem.  The row law
may be any exact-sign law whose density with respect to the uniform cube has
uniformly bounded `L2` norm.  The mean is added to the exceptional response
space; the resulting joined projection need not commute with the second
moment, and that point is handled explicitly below.

The random-matrix inputs are the same Strong Tail Projection upper-edge
theorem of Chafaï--Tikhomirov and Marchenko--Pastur sufficiency theorem of
Yaskov used in the companion sharp-edge report.

## 1. Setup and statements

Let `U_r` be uniform on `{+-1}^r`, and let `mu_r` be an arbitrary probability
law on that cube with density

```math
g_r={d\mu_r\over dU_r},
\qquad
g_r\ge0,
\qquad
\mathbb E_{U_r}g_r=1,
\qquad
\mathbb E_{U_r}g_r^2\le K
\tag{BL.1}
```

for a fixed finite `K`.  If `R_r~mu_r`, write

```math
m_r=\mathbb E R_r,
\qquad
\Sigma_r=\mathbb E R_rR_r^T.
\tag{BL.2}
```

Thus `Sigma_r` is the uncentered second moment.  Let
`delta_r=r^(-1/4)`.  Let `P_r^0` be the spectral projection of `Sigma_r`
onto eigenvalues outside `[1-delta_r,1+delta_r]`, and let `P_r` be the
orthogonal projection onto

```math
\operatorname{ran}P_r^0+\operatorname{span}\{m_r\}.
\tag{BL.3}
```

Put `V_r=I-P_r`.  If `m_r=0`, the second summand is omitted.

Let `B_r` be the square matrix with independent rows distributed as `R_r`.

### Theorem BL.1 (noncentral bounded-L2 sharp projected edge)

Under (BL.1),

```math
\boxed{
{\|B_rV_r\|_{op}\over\sqrt r}\longrightarrow2
\quad\hbox{in probability}.}
\tag{BL.4}
```

The assertion is uniform over all triangular sequences satisfying (BL.1)
with the same `K`.

### Theorem BL.2 (bounded-L2 product rows cannot lower pressure)

For every fixed

```math
0<\beta<{\sqrt2\over6},
\tag{BL.5}
```

the corresponding conference pressure satisfies

```math
\boxed{
\mathbb E\left[
\left(h_\beta-{f_r(B_r)\over r}\right)_+
\right]\longrightarrow0.}
\tag{BL.6}
```

Here `f_r` and `h_beta` have the normalization of the companion projected-
coupling theorem.  No symmetry of `mu_r` is assumed.

## 2. Fourier control of the mean and second moment

For nonempty `S subseteq [r]`, let

```math
a_S=\langle g_r,\chi_S\rangle_{U_r}.
\tag{BL.7}
```

Normalized Parseval gives

```math
\sum_{S\ne\varnothing}a_S^2
=\mathbb E_Ug_r^2-1
\le K-1.
\tag{BL.8}
```

The singleton coefficients are the coordinates of `m_r`.  Since the rows
are exact signs, `Sigma_r` has diagonal one, and for `i != j` its `(i,j)`
entry is `a_{\{i,j\}}`.  Hence

```math
\boxed{
\|m_r\|_2^2+{1\over2}\|\Sigma_r-I\|_F^2
\le K-1.}
\tag{BL.9}
```

In particular,

```math
\operatorname{rank}P_r^0
\le {2(K-1)\over\delta_r^2},
\qquad
k_r:=\operatorname{rank}P_r
\le {2(K-1)\over\delta_r^2}+1
=O_K(\sqrt r)=o(r).
\tag{BL.10}
```

Also

```math
\|\Sigma_r\|_{op}
\le1+\|\Sigma_r-I\|_F
\le 1+\sqrt{2(K-1)}=:L_K.
\tag{BL.11}
```

## 3. Why the noncommuting joined projection is harmless

In general `P_r` is not a spectral projection of `Sigma_r`, because the
mean need not be invariant under `Sigma_r`.  It would therefore be invalid
to claim that `P_r` commutes with `Sigma_r`.

What is needed is only a compression bound.  Since

```math
\operatorname{ran}V_r
\subseteq\operatorname{ran}(I-P_r^0),
\tag{BL.12}
```

every `x in ran(V_r)` lies in the regular spectral subspace of `Sigma_r`.
Thus the spectral theorem gives the Rayleigh-quotient inequalities

```math
(1-\delta_r)\|x\|_2^2
\le x^T\Sigma_rx
\le(1+\delta_r)\|x\|_2^2.
\tag{BL.13}
```

Moreover, by construction,

```math
V_rm_r=0.
\tag{BL.14}
```

Let `d_r=rank(V_r)=r-k_r` and, on `ran(V_r)`, define

```math
S_r=V_r\Sigma_rV_r|_{\operatorname{ran}V_r}.
\tag{BL.15}
```

Equation (BL.13), not commutation, implies

```math
(1-\delta_r)I\preceq S_r\preceq(1+\delta_r)I.
\tag{BL.16}
```

Therefore

```math
X_r=S_r^{-1/2}V_rR_r\in\mathbb R^{d_r}
\tag{BL.17}
```

is centered by (BL.14) and isotropic because

```math
\mathbb E X_rX_r^T
=S_r^{-1/2}(V_r\Sigma_rV_r)S_r^{-1/2}=I_{d_r}.
\tag{BL.18}
```

This is the only place the mean direction is required.

## 4. L2 change of measure preserves projection tails

For every event `A` in the cube, Cauchy--Schwarz and (BL.1) give

```math
\boxed{
\mu_r(A)\le\sqrt K\,U_r(A)^{1/2}.}
\tag{BL.19}
```

This replaces the bounded-likelihood-ratio step in the constant-density
proof.  It halves exponential tail rates but preserves their scale.

Let `Q` be an orthogonal projection on `R^{d_r}` of rank `ell>=1`.  In the
original coordinates,

```math
\|QX_r\|_2^2=R_r^TMR_r,
\quad
M=V_rS_r^{-1/2}QS_r^{-1/2}V_r.
\tag{BL.20}
```

By (BL.16),

```math
\|M\|_{op}\le(1-\delta_r)^{-1},
\qquad
\|M\|_F^2\le\ell(1-\delta_r)^{-2}.
\tag{BL.21}
```

Under `mu_r`, isotropy gives

```math
\mathbb E_{\mu_r}R_r^TMR_r=\ell.
\tag{BL.22}
```

Under the uniform law the expectation is `tr M`.  In coordinates on
`ran(V_r)`,

```math
\operatorname{tr}M=\operatorname{tr}(S_r^{-1}Q),
\tag{BL.23}
```

so (BL.16) implies

```math
|\operatorname{tr}M-\ell|
\le {\delta_r\over1-\delta_r}\ell.
\tag{BL.24}
```

For a uniform Rademacher vector `W_r`, Hanson--Wright and (BL.21) give

```math
U_r\{|W_r^TMW_r-\operatorname{tr}M|\ge u\}
\le C\exp\{-c\min(u^2/\ell,u)\}.
\tag{BL.25}
```

Apply (BL.19), absorb the square root by changing absolute constants, and
use

```math
{\delta_r\over1-\delta_r}\ell
\le2\ell^{3/4}
\qquad(1\le\ell\le d_r\le r)
\tag{BL.26}
```

for all large `r`.  Whenever `t>=4 ell^(3/4)`, this yields

```math
\Pr\{\|QX_r\|_2^2-\ell\ge t\}
\le C_K\exp\{-c_K\min(t^2/\ell,t)\}.
\tag{BL.27}
```

Consequently the array `X_r` satisfies the Strong Tail Projection property
with, for example,

```math
f(\ell)=\min(1,4\ell^{-1/4})
\tag{BL.28}
```

and a function `g(ell)->0` obtained by dominating the exponential right
side of (BL.27) by `g(ell) ell/t^2`.  The functions depend only on `K`, not
on the particular density.  Coordinate independence of `X_r` is neither
asserted nor needed.

## 5. Marchenko--Pastur condition

Let `A_r` be any positive semidefinite `d_r by d_r` matrix with uniformly
bounded operator norm, and set

```math
M_r=V_rS_r^{-1/2}A_rS_r^{-1/2}V_r.
\tag{BL.29}
```

Then `||M_r||op=O(1)`, `||M_r||F^2=O(d_r)`, and

```math
|\operatorname{tr}M_r-\operatorname{tr}A_r|
\le {\delta_r\over1-\delta_r}\operatorname{tr}A_r
=o(d_r).
\tag{BL.30}
```

Hanson--Wright followed by (BL.19) shows, for every fixed `epsilon>0`,

```math
\Pr\{|R_r^TM_rR_r-\operatorname{tr}M_r|>\epsilon d_r\}
\le C_{K,\epsilon}e^{-c_{K,\epsilon}d_r}.
\tag{BL.31}
```

Combining (BL.30)--(BL.31),

```math
{X_r^TA_rX_r-\operatorname{tr}A_r\over d_r}
\longrightarrow0
\quad\hbox{in probability}.
\tag{BL.32}
```

This is Yaskov's weak quadratic-form condition.

## 6. Proof of the sharp edge

Let `mathbb X_r` be the `r by d_r` matrix with independent rows distributed
as `X_r`.  Since `d_r/r->1`, the Strong Tail Projection theorem gives

```math
\limsup_{r\to\infty}
{\mathbb E\lambda_{\max}(\mathbb X_r^T\mathbb X_r)
 \over(\sqrt r+\sqrt{d_r})^2}\le1.
\tag{BL.33}
```

The Marchenko--Pastur theorem applied through (BL.32) gives empirical
spectral convergence at aspect ratio one.  Its positive mass immediately
below the upper edge four implies

```math
\liminf_{r\to\infty}
{\lambda_{\max}(\mathbb X_r^T\mathbb X_r)\over r}
\ge4
\quad\hbox{in probability}.
\tag{BL.34}
```

The expectation upper bound and probability lower bound together force

```math
{\|\mathbb X_r\|_{op}\over\sqrt r}\longrightarrow2
\quad\hbox{in probability}.
\tag{BL.35}
```

In compatible coordinates,

```math
B_rV_r=\mathbb X_rS_r^{1/2}.
\tag{BL.36}
```

The two singular-value inequalities from (BL.16) sandwich its norm between
`sqrt(1-delta_r)||mathbb X_r||op` and
`sqrt(1+delta_r)||mathbb X_r||op`.  Equations (BL.35)--(BL.36) prove BL.4.
`square`

Every tail and centering estimate above is uniform in the density subject
to the same `K`.  Equivalently, if uniform convergence failed, one could
choose a violating density along a subsequence and fill the remaining
orders arbitrarily, contradicting BL.1 for that triangular sequence.

## 7. Nuclear cost and Hamming transport

The joined projection has the same subcritical cost as the earlier spectral
peel.  By (BL.11),

```math
\mathbb E\|P_rR_r\|_2^2
=\operatorname{tr}(P_r\Sigma_r)
\le L_Kk_r.
\tag{BL.37}
```

Therefore, using `rank(B_rP_r)<=k_r`,

```math
\mathbb E\|B_rP_r\|_*
\le\sqrt{k_r}(rL_Kk_r)^{1/2}
=\sqrt{L_K}\,k_r\sqrt r
=O_K(r)=o(r^{3/2}).
\tag{BL.38}
```

For an iid Rademacher bridge `W_r`, the corresponding expectation is at
most `k_r sqrt(r)=O_K(r)`.

The `L2` hypothesis also controls entropy.  Jensen's inequality under
`mu_r` gives

```math
D(\mu_r\|U_r)
=\mathbb E_{\mu_r}\log g_r
\le\log\mathbb E_{\mu_r}g_r
=\log\mathbb E_{U_r}g_r^2
\le\log K.
\tag{BL.39}
```

The cube `T1` inequality thus supplies a coupling `(R_r,W_r)` with

```math
\mathbb E d_H(R_r,W_r)
\le\sqrt{{r\log K\over2}}.
\tag{BL.40}
```

Couple the rows independently.  Exact signs then give

```math
\mathbb E\|(B_r-W_r)V_r\|_F
\le\left(4r\,\mathbb E d_H(R_r,W_r)\right)^{1/2}
=O_K(r^{3/4})=o(r).
\tag{BL.41}
```

Equations (BL.4), (BL.38), and (BL.41) verify every hypothesis of the
audited projected-coupling theorem.  This proves BL.6.  Notice that the
possibly macroscopic mean response lies in `P_r` and is restored by the
one-sided convex supporting inequality; the full endpoint need not be
operator-regular.  `square`

## 8. Latent-mixture corollary

The uniformity gives a useful extension beyond a single product law.

### Corollary BL.3 (uniformly L2-bounded latent products have no gain)

Let

```math
q_r=\int \mu_{r,z}^{\otimes r}\,\pi_r(dz)
\tag{BL.42}
```

be any mixture of row-product exact-sign bridge laws.  Suppose every
component has a density `g_{r,z}` relative to `U_r` satisfying

```math
\mathbb E_U g_{r,z}^2\le K
\tag{BL.43}
```

with the same fixed `K`, for `pi_r`-almost every `z`.  Then, for every fixed
`beta<sqrt(2)/6`,

```math
\boxed{
\mathbb E_{C_r\sim q_r}\left[
\left(h_\beta-{f_r(C_r)\over r}\right)_+
\right]\longrightarrow0.}
\tag{BL.44}
```

**Proof.**  Apply the uniform form of BL.2 conditionally on `z` and
integrate.  The exceptional projection may depend on `z`; only the uniform
error bound is integrated, so no common projection is required.  `square`

This corollary can cover bridge laws with linear row total correlation and
arbitrarily large latent support.  It does not contradict the separate
total-correlation lower bound: a reusable latent variable is harmless here
because every conditional product component itself has bounded `L2`
density.  A favorable latent mixture must contain components whose row
likelihood becomes increasingly concentrated, or must retain dependence
even after conditioning on the proposed latent state.

## 9. Limitations and stress tests

1. **No central symmetry is used.**  The mean may be nonzero, but (BL.9)
   bounds its squared norm by `K-1`, and one additional response direction
   removes it exactly.

2. **The joined projection is not assumed invariant.**  The proof uses only
   the subspace inclusion (BL.12) and compression inequalities (BL.13).
   Replacing those by a false commutation assertion would create a gap.

3. **A biased half-cube is covered.**  For example, conditioning on a
   majority halfspace gives a nonzero diffuse mean.  Its span is peeled;
   the remaining centered rows have sharp edge two.  Thus a mean spike is
   not a counterexample.

4. **Arbitrary weighted laws are covered.**  The proof uses neither a hard
   event nor a pointwise density bound.  The only change-of-measure step is
   Cauchy--Schwarz in (BL.19).

5. **Uniform `L2` control is essential to this argument.**  If
   `E_U g_r^2` diverges, the Fourier rank, change-of-measure tails, entropy
   transport, and mixture uniformity all change scale.  The theorem makes
   no claim in that regime.

6. **The theorem remains one-sided and high-temperature.**  It excludes a
   lower pressure phase only for fixed `beta<sqrt(2)/6`; it does not identify
   the full pressure law, cover the endpoint, or prohibit pressure increase.

## 10. Frontier movement

Central symmetry and literal constant-density fibres are not structural
requirements of the no-gain mechanism.  The real row-level requirements are
bounded quadratic density information, an `o(r)`-rank mean/covariance peel,
and projection-tail concentration.  Uniform `L2` density supplies all three.

Moreover, cross-row total correlation alone is not the final obstruction:
even a highly correlated latent mixture is harmless when it decomposes into
uniformly bounded-L2 product components.  A surviving favorable law must
therefore have either increasingly singular conditional row components or
irreducible dependence not captured by such a latent product decomposition.
