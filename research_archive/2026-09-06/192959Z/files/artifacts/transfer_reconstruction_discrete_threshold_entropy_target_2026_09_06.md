# Discrete threshold entropy: exact curvature and an essential flatness obstruction

Date: 2026-09-06. The requested **flat-sign, bounded-operator-norm** entropy
upper bound remains unproved and is not falsified here. This note proves its
exact infinitesimal coefficient and a scalable obstruction to replacing its
flatness assumption by operator/Frobenius bounds alone. The obstruction
isolates a genuinely non-Gaussian rare-threshold effect.

## 1. The exact target and what it would give

Let `A_N` be a symmetric hollow sign matrix with
`||A_N||op<=C sqrt N`. Fix `delta in (0,1/2)`, `s`, and `C` before
taking the order limit, with `|s|C<1`. Put

```math
R_N(s)=I+sA_N/\sqrt N,\quad
G\sim N(0,R_N(s)),\quad X_i=\operatorname{sign}(G_i-t),
\quad t=\Phi^{-1}(\delta),
```

and `m=1-2delta`, `v=1-m^2`, `a=2phi(t)`, `eta=a^2/v`.
The precise open bound is of the form

```math
D(\mathcal L(X)\Vert\mu_m)
\le F(\delta,s,C)\frac{\eta^2s^2}{4}N+o(N),             \tag{1}
```

uniformly over the actual sign matrices, with an explicit factor small
enough to improve the existing clipped-law credit. An unspecified factor
depending arbitrarily on the rare probability would not meet this target.

For these explicit matrices `R_N(s)`, the Hermite remainder from the
threshold-resummation proof contributes only `O(s^2 sqrt N)` to energy.
Thus (1) would give the correlation credit per spin

```math
\frac{\beta a^2s}{2}-F\frac{a^4s^2}{4v^2}.              \tag{2}
```

This identifies the quantitative burden. At the current flip density,
`beta v^2/a^2` is of order one, so a fixed spectral-margin choice of `s`
could be useful. The exact entropy, rather than the Gaussian data-processing
upper bound on it, is required.

## 2. Exact infinitesimal coefficient

For every fixed finite matrix, let `p_s` be the thresholded law and
`p_0=mu_m`. The Gaussian density score at `s=0` is

```math
\left.\partial_s\log\frac{dN(0,R_N(s))}{dN(0,I)}\right|_{s=0}
=\sum_{i<j}\frac{A_{ij}}{\sqrt N}G_iG_j.
```

Conditional expectation on the threshold cells gives the discrete score.
Indeed

```math
\mathbb E[G_i\mid X_i]=\frac a v(X_i-m).
```

Writing `U_i=(X_i-m)/sqrt v`, the exact threshold score is therefore

```math
L_1(X)=\left.\partial_s\frac{p_s(X)}{p_0(X)}\right|_{s=0}
=\eta\sum_{i<j}\frac{A_{ij}}{\sqrt N}U_iU_j.             \tag{3}
```

The products `U_iU_j` are orthonormal under the independent product law.
Consequently

```math
D(p_0\Vert p_0)=D'(0)=0,\qquad
D''(0)=\mathbb E L_1^2=\eta^2\frac{N-1}{2}.              \tag{4}
```

All differentiations are justified at fixed `N` by Gaussian domination in
a sufficiently small covariance neighborhood and by the finite positive
threshold-cell probabilities. In particular

```math
D(p_s\Vert p_0)=\frac{\eta^2(N-1)}4s^2+O_{N,\delta,A}(s^3).
```

The dependence of the remainder on `N` is deliberately displayed. This
calculation does not prove (1) at any fixed nonzero `s` as `N` grows.

## 3. Scalable failure of an operator-only entropy contraction

Suppose one drops entrywise flatness and only assumes a bounded covariance
perturbation norm and extensive squared Frobenius norm. For even `N`, take
`R_N` to be `N/2` independent copies of

```math
\begin{pmatrix}1&\rho\\\rho&1\end{pmatrix},\qquad 0<\rho<1.
```

Then `||R_N-I||op=rho` and `tr(R_N-I)^2=rho^2 N`. Set
`X_i=1-2*1{G_i>=z}` at `z>0`, so that its mean remains `m=1-2delta`,
and write `delta=Phi(-z)` for the rare probability. If `I_delta(rho)` is
the mutual information of one thresholded pair, exactly

```math
D(\mathcal L(X)\Vert\mu_m)=\frac N2 I_\delta(\rho).      \tag{5}
```

For every fixed `rho>0`,

```math
\frac{I_\delta(\rho)}{\eta^2\rho^2}\longrightarrow\infty
\quad\text{as }\delta\downarrow0.                       \tag{6}
```

Here `eta=a^2/v=phi(z)^2/[delta(1-delta)]`, exactly as above.
Thus no factor depending only on a fixed spectral margin can supply the
proposed rare-threshold contraction for this larger matrix class.

An explicit elementary proof is enough. Let
`r_z=P(G_1>=z,G_2>=z)`. Integrating the bivariate Gaussian density on
`[z,z+1/z]^2` gives, for all sufficiently large `z`,

```math
r_z\ge c_\rho z^{-2}\exp\{-z^2/(1+\rho)\},             \tag{7}
```

where `c_rho>0` is fixed. The Gaussian tail upper bound is
`delta<=phi(z)/z`. Therefore `log(r_z/delta^2)` grows at least as a
positive constant times `z^2`. Applying data processing to the event
that both bits are rare gives

```math
I_\delta(\rho)\ge
D(\operatorname{Ber}(r_z)\Vert\operatorname{Ber}(\delta^2))
\ge r_z\log(r_z/\delta^2)-r_z.
```

The lower Mills bound `delta>=phi(z)z/(1+z^2)` implies
`eta^2=O(z^2 exp(-z^2))`. Combining these estimates proves (6).

For a fully explicit choice `rho=1/2` and `z>=4`, the rectangle calculation
can use `c_rho=e^{-2}/(pi sqrt3)`. It yields

```math
I_\delta(1/2)\ge\frac{e^{-2}}{12\pi\sqrt3}e^{-2z^2/3},
\qquad
\frac{I_\delta(1/2)}{\eta^2(1/2)^2}
\ge\frac{e^{-2}}{24\sqrt3}\frac{e^{z^2/3}}{z^2}.         \tag{8}
```

The perturbation in this example has fixed nonzero entries, not entries
of size `s/sqrt N`. Expressing it as `s A_N/sqrt N` would require entries
of `A_N` of order `sqrt N`, not a signing. Thus (5)--(8) are **not** a
counterexample to the actual target (1). They rule out a factor uniform
as `delta` tends to zero using only spectral stability and Frobenius norm.
They do not rule out an explicit `delta`-dependent operator-only bound,
nor prove that every successful proof of (1) must use flatness.

## 4. Primary literature checked and what it does not settle

Macke, Opper, and Bethge analyze higher-order correlations and entropy in
a thresholded common-Gaussian-input model. Their homogeneous model keeps
the off-diagonal Gaussian correlation fixed as the population grows; its
largest covariance eigenvalue consequently grows linearly. The explicit
derivations are useful evidence of rare-threshold nonlinear effects, but
do not meet the bounded-norm, vanishing-entry hypotheses in (1).
[Original paper](https://arxiv.org/abs/1009.2855),
[authors' detailed supplement, Section I.B](https://www.mackelab.org/pubs/Macke_Opper_2011_supplement.pdf).

No primary theorem found in this bounded search supplies the uniform
discrete-entropy remainder required by (1). The previously banked
Gaussian-threshold lower bound and clipped correlated-law lower bound
remain valid; neither is weakened by the present obstruction.

## 5. Reproduction and current stopping point

`computations/transfer_reconstruction_threshold_entropy_exact_checks_2026_09_06.py`
replays (3)--(4) exactly for all order-four signings and separately
computes labeled numerical pair-resonance examples. The fixed-flat-sign
entropy upper bound, or a counterexample satisfying all of its hypotheses,
is still the missing result. No claim of a dimension-uniform Taylor theorem
or a new cap coefficient follows from this note.
