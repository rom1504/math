# Adversarial audit: fixed negative-tilt product-value reduction

Status: **independent proof audit passed, including the cap-free sharpening;
narrow RESET, Level 5**.

The audited claims are NT.1--NT.16 in
[`actual_child_negative_tilt_product_value.md`](actual_child_negative_tilt_product_value.md).
The audit checks the variational signs, the integrated Herbst factor, the
raw/capped comparison, and the translation to alternatives (ii)--(iii).

## 1. Variational sign

For `F=L wedge CN`, put

```math
\mathcal G_F(P)=E_PF+\lambda^{-1}D(P\Vert U),
\qquad
V_C^{\rm row}=\inf_{P\ {\rm product}}\mathcal G_F(P).
```

Donsker--Varadhan gives, over all bridge laws,

```math
\inf_P\mathcal G_F(P)=-\lambda^{-1}\log E_Ue^{-\lambda F}.
```

The row products are a restriction and `U` is a row product.  Therefore

```math
-\lambda^{-1}\log E_Ue^{-\lambda F}
\le V_C^{\rm row}\le E_UF,
```

which is exactly the direction required in NT.6.  There is no illicit
exchange of a minimum and expectation.

## 2. Integrated logarithmic-Sobolev factor

Let

```math
\psi(u)=\log E_Ue^{-u(F-E_UF)},
\qquad 0\le u\le\lambda.
```

With the half-flip normalization, the fair-cube logarithmic-Sobolev and
exponential-gradient estimates give

```math
{d\over du}{\psi(u)\over u}
\le C_{\rm LS}E_{\Pi_{-u,F}}\Gamma_F.
```

Since `psi(u)/u->0` at zero, integration yields

```math
{\psi(\lambda)\over\lambda}
\le C_{\rm LS}\int_0^\lambda
 E_{\Pi_{-u,F}}\Gamma_F\,du.
```

The exact cavity estimate
`Gamma_F<=t^2 sum_a r_a^2` then gives

```math
{\psi(\lambda)\over\lambda}
\le C_{\rm LS}\lambda t^2mn
 \overline\rho_N^-(\lambda).
```

Thus NT.7a has the correct power of `lambda`; replacing the integral by a
supremum is unnecessary.

## 3. Raw/capped comparison

Because `F<=L`,

```math
V_{\lambda,F}:=-\lambda^{-1}\log E_Ue^{-\lambda F}
\le V_{\lambda,L}.
```

Hence

```math
E_UL-V_{\lambda,L}
\le(E_UL-E_UF)+(E_UF-V_{\lambda,F}).
```

The first term is `O(N^(3/2)e^(-cN))` by RT.2 with the fair product.
The same truncation theorem applies to both the raw and capped row-product
optimizers because their coordinate best responses have the common row
collision bound.  Therefore NT.8 and the raw branch statement do not assume
that the full joint escort itself has bounded collision norm.

## 4. Branch implication

The fair bridge law is an admissible competitor for the reverse product
projection, so

```math
\mathcal I_\lambda^{\leftarrow}
\le D(U\Vert q_\lambda)
=\lambda(E_UL-V_{\lambda,L}).
```

Combining Sections 2--3 proves

```math
\mathcal I_\lambda^{\leftarrow}
\le C_{\rm LS}\lambda^2t^2mn
 \overline\rho_N^-(\lambda)+O(N^{3/2}e^{-cN}).
```

At a balanced split `t^2mn=Theta(N)`.  Thus:

- `bar rho_N^-=o(1)` implies `I^leftarrow=o(N)`;
- `I^leftarrow>=eta N` implies `bar rho_N^->=c>0`;
- if also `J>=eta N`, then
  `J-I^leftarrow>=eta N-o(N)`.

This is a valid one-sided separator between alternatives (ii) and (iii).
Positive integrated overlap alone does not prove alternative (ii), and the
main note does not claim otherwise.

## 5. Cap-free sharpening on the actual negative path

The cap is not needed for the branch separator.  Put

```math
{d\widehat\Pi_s\over dU}={e^{sL}\over E_Ue^{sL}},
\qquad
\widehat\rho_N^-(\lambda)={1\over\lambda mn}
 \int_{-\lambda}^0E_{\widehat\Pi_s}\sum_a r_a^2\,ds.
```

Repeating Sections 1--2 directly with `F=L` proves

```math
\boxed{
\begin{aligned}
0\le E_UL-V_\lambda^{\rm row}
 &\le C_{\rm LS}\lambda t^2mn\widehat\rho_N^-(\lambda),\\
\mathcal I_\lambda^{\leftarrow}
 &\le C_{\rm LS}\lambda^2t^2mn\widehat\rho_N^-(\lambda).
\end{aligned}}
```

All exponential-gradient integrands are dominated on `[-lambda,0]` because
`L>=0` and `e^(sL)<=1`.  At the left endpoint,
`widehat Pi_(-lambda)=q_lambda` is exactly the actual negative-disorder law.
The raw Gaussian-smoothing computation likewise gives

```math
\widehat\rho_N^-(\lambda)
={1\over\lambda}\int_{-\lambda}^0
 {\widehat A_N(s)-1\over s-1}\,ds+O_\beta(N^{-1/2}).
```

The denominator stays away from zero.  Thus the proof uses no cap, tail
comparison, replica-one mixed derivative, or growing window.  This exact
version supersedes the capped estimate for deciding the optimized product
phase; the capped theorem remains useful only for two-sided uniform response
questions.

## 6. Information and recurrence scope

The statistic is carrier-independent and uses one fixed negative interval.
Via the raw heat identity it is one integrated smoothing secant, with no
replica-one mixed derivative or cap.  It is analytically smaller than the
carrier response table, but exact evaluation still invokes the actual-child
Gibbs law; no polynomial-time or finite-state realization is proved.

There is no Level-6 recurrence because neither
`hat rho_N^-=o(1)` nor a directional converse from positive overlap is
known.  The correct classification is a narrow RESET: the previous
bounded-row-degree cross-row optimization has been replaced, for the branch
decision, by one integrated fixed-window actual-child observable.
