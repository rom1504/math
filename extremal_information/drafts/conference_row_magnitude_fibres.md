# Constant-probability row-magnitude fibres retain conference pressure at sufficiently high temperature

**Status.** Task-local theorem report; no canonical files are edited.  This
tests the sign-invariant rowwise magnitude candidate left open by the
halfcube/gauge audit.  The result is a genuine quenched pressure theorem in
a nonempty strict-high-temperature interval.  It deliberately distinguishes
the theorem from the weaker observation that the conditioned row covariance
has only a finite-rank defect.

## 1. Magnitude fibres and exact entropy

Fix sign vectors `v_r in {+-1}^r`.  Let `I_r` be an arbitrary subset of the
possible absolute row sums and define

```math
E_r=\{R\in\{+-1\}^r:|\langle R,v_r\rangle|\in I_r\}.
\tag{RM.1}
```

Assume that its uniform-row probability satisfies

```math
p_r=2^{-r}|E_r|\ge p_0>0.
\tag{RM.2}
```

The row-product fibre is

```math
\mathcal M_r(E_r)=\{B:B_{i,*}\in E_r\text{ for all }i\}.
\tag{RM.3}
```

It has the exact cardinality

```math
\boxed{
|\mathcal M_r(E_r)|=|E_r|^r
=2^{r^2}p_r^r.}
\tag{RM.4}
```

Thus every fixed nontrivial row-magnitude constraint has precisely speed
`r` entropy loss.  It is invariant under every individual row sign and is
not a switching cross-section.

Let `A_r` be a symmetric conference signing, fix an orientation `epsilon`,
and retain the pressure notation

```math
f_{\epsilon,r}(B)=
\log\left[2^{-2r}\sum_{x,y}
\cosh\left\{{\beta\over\sqrt{2r}}
\big(H_A(x)+\epsilon H_A(y)+x^TBy\big)\right\}\right].
\tag{RM.5}
```

## 2. What covariance says, and what it does not say

After the column gauge `v_r -> 1`, a uniform row conditioned on `E_r` is
exchangeable and centrally symmetric.  Put

```math
\alpha_r={1\over r}
 \mathbb E\left[\left(\sum_jR_j\right)^2\middle|E_r\right].
\tag{RM.6}
```

Its covariance has eigenvalue `alpha_r` in the all-ones direction and the
repeated eigenvalue

```math
{r-\alpha_r\over r-1}
\tag{RM.7}
```

on its orthogonal complement.  Equation (RM.2) and the unconditioned second
moment give `alpha_r<=1/p_0`.  Hence the covariance empirical law converges
to `delta_1`; the conditioning changes at most one population direction.

This is only bulk-spectral evidence.  A finite-rank direction can change an
Ising free energy by order `r` after a mean-field instability, and empirical
spectral convergence by itself does not verify Boolean pressure.  The proof
below instead couples the actual landscapes and invokes a pressure stability
theorem under a verified strict operator bound.

## 3. A row-layer coupling

### Lemma RM.1 (uniform and conditioned rows couple within `O(sqrt r)` edits)

There is a coupling `(W,R)` in which `W` is a uniform sign row, `R` is
uniform on `E_r`, and

```math
\mathbb E d_H(W,R)\le C_{p_0}\sqrt r.
\tag{RM.8}
```

The coupling can be chosen independently from row to row.

**Proof.**  Apply the column gauge and let `K,K'` be the numbers of plus
signs in `W,R`.  Give `K` its binomial law and `K'` its law conditioned on
`|2K'-r| in I_r`.  They may initially be coupled independently.  Since

```math
\mathbb E|2K-r|\le\sqrt r,
\qquad
\mathbb E[|2K'-r|]
\le {\mathbb E|2K-r|\over p_r}
\le {\sqrt r\over p_0},
\tag{RM.9}
```

one has `E|K-K'|<=C_(p_0)sqrt(r)`.

Conditional on `K=k,K'=k'`, couple the two uniform plus sets by nesting the
smaller uniformly inside the larger.  Symmetry preserves both uniform
layer marginals and makes their Hamming distance exactly `|k-k'|`.  Undo
the column gauge. `square`

For `r` independent copies, write the coupled matrices as `W_r,R_r`.  Then

```math
\mathbb E\|R_r-W_r\|_F^2
=4\mathbb E d_H(R_r,W_r)
\le C_{p_0}r^{3/2},
\tag{RM.10}
```

and consequently

```math
\mathbb E\|R_r-W_r\|_F=O_{p_0}(r^{3/4})=o(r).
\tag{RM.11}
```

## 4. A uniform operator bound for the conditioned rows

### Lemma RM.2 (crude but uniform dependent-row norm bound)

There are constants `L_(p_0),c_(p_0)>0` such that

```math
\Pr\{\|R_r\|_{op}>L_{p_0}\sqrt r\}
\le2e^{-c_{p_0}r}.
\tag{RM.12}
```

**Proof.**  Central symmetry of `E_r` gives `E R=0`.  For every deterministic
unit vector `z`, Hoeffding's inequality before conditioning gives

```math
\Pr\{|\langle R,z\rangle|>u\mid E_r\}
\le {2\over p_0}e^{-u^2/2}.
\tag{RM.13}
```

Thus all row linear forms have a subgaussian norm bounded only in terms of
`p_0`.  For fixed `z`, Bernstein's inequality for the independent
subexponential variables `\langle R_i,z\rangle^2` bounds
`\sum_i\langle R_i,z\rangle^2` by `C_(p_0)r` outside an
`e^{-c_(p_0)r}` event.  A fixed `1/4`-net of the unit sphere has at most
`9^r` points.  Increasing `C_(p_0)` so the Bernstein exponent dominates
`r log 9`, and using the standard net comparison for operator norm, proves
(RM.12). `square`

The constant is intentionally not claimed sharp.  A Bai--Yin edge `2` for
these dependent-coordinate rows would enlarge the temperature interval,
but is not needed for the theorem below.

## 5. Quenched pressure theorem

### Theorem RM.3 (row-magnitude fibres are not favorable at sufficiently high temperature)

For every `p_0>0` there is `beta_0(p_0)>0` such that, whenever

```math
0<\beta<\beta_0(p_0),
\tag{RM.14}
```

the following holds uniformly over all choices (RM.1)--(RM.2).  If `B_r`
is uniform on `M_r(E_r)`, then, for both orientations,

```math
\boxed{
{f_{\epsilon,r}(B_r)\over r}\longrightarrow
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4}}
\tag{RM.15}
```

in probability and in `L^1`.  In particular, these exact speed-`r` fibres
do not reach the smaller same-temperature target
`tau_beta=h_beta-gamma(beta)`.

One may take any `beta_0` small enough that, for some `kappa<1/2`,

```math
{\beta_0\over\sqrt2}
\left(1+\max\{3,L_{p_0}\}\right)<\kappa.
\tag{RM.16}
```

**Proof.**  Use the coupling in RM.1.  The uniform bridge satisfies
`||W_r||_op<3sqrt(r)` outside an `e^{-cr}` event.  Lemma RM.2 supplies the
analogous bound for `R_r`.  By (RM.16), on their intersection the two
scaled parent interactions

```math
X_B={\beta\over\sqrt{2r}}
\begin{pmatrix}A_r&B\\B^T&\epsilon A_r\end{pmatrix}
\tag{RM.17}
```

have operator norm at most `kappa`.

The audited strict-high-temperature pressure stability theorem gives

```math
\begin{aligned}
|f(R_r)-f(W_r)|
&\le {K_\kappa\over2}\|X_{R_r}-X_{W_r}\|_*\\
&\le {K_\kappa\beta\over\sqrt2}\|R_r-W_r\|_F.
\end{aligned}
\tag{RM.18}
```

Indeed the symmetric off-diagonal block difference has nuclear norm twice
that of the bridge difference, and
`||R-W||_*<=sqrt(r)||R-W||_F`.  Equations (RM.11) and (RM.18) make the
expected good-event pressure difference `O(r^(3/4))=o(r)`.

On the complement, both sign-parent pressures lie between `0` and
`C_beta r^(3/2)`, while the bad operator event has probability `e^{-c r}`.
Thus

```math
\mathbb E|f(R_r)-f(W_r)|=o(r).
\tag{RM.19}
```

The uniform conference theorem gives `f(W_r)/r -> h_beta` in probability
and mean.  Equation (RM.19) transfers both conclusions to `R_r`, proving
(RM.15). `square`

## 6. What this resolves and what remains

The theorem covers every rowwise constraint determined by
`|<R,v_r>|` whose one-row probability stays bounded below: central bands,
tail bands with a fixed normalized threshold, finite unions of magnitude
intervals, and exact parity-compatible variants having constant total row
mass.  It is a pressure theorem, not an ESD conjecture.

There are two precise limitations.

1. The interval `beta<beta_0(p_0)` comes from the crude net norm constant.
   Extending all the way to the campaign range `beta<sqrt(2)/6` requires a
   sharp operator-edge theorem (or another pressure argument) for the
   independent rows with dependent coordinates.  Bulk covariance and an
   anticipated Marchenko--Pastur law do not supply this by themselves.
2. The result assumes one-row mass bounded below.  If `p_r` decays with
   `r`, the fibre entropy can exceed speed `r`, the subgaussian constant
   deteriorates, and the theorem makes no assertion.

The candidate is therefore closed rigorously in a nonempty uniform
high-temperature regime and reduced, in the larger regime, to a concrete
random-matrix/operator-margin question.  It does not exhibit a speed-`r`
favorable basin.  The next genuinely different candidate would need
non-product correlations between rows or a row constraint not determined
by a single magnitude statistic.
