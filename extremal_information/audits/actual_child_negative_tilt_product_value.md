# The fixed negative-tilt overlap controls the optimized product value

Status: **rigorous actual-child full-product value reduction**.  The common
compact-tilt response theorem can be strengthened for the quantity that is
actually optimized.  Relative entropy is already part of the row-product
objective, so Donsker--Varadhan duality cancels its entire cost at the single
physical parameter `-lambda`.  No carrier alphabet, growing tilt window, or
positive-tilt cap plateau remains.

This theorem does not prove overlap decay.  It reduces the row-product gain
to one fixed negative-tilt cavity-overlap envelope of the actual children.

## 1. Setup

Let `A,D` be actual contracted-temperature minimizing children, let `L` be
their bridge pressure, and let

```math
F=L\wedge CN                                                    \tag{NT.1}
```

be any cap satisfying the full bounded-`L^2` response theorem RT.2.  For a
row-product law `P` define

```math
\mathcal F_C(P)=E_PF+{1\over\lambda}D(P\Vert U_{mn}),
\qquad
V_{\lambda,C}^{\rm row}=\inf_{P\in\mathcal P_{\rm row}}
\mathcal F_C(P).                                  \tag{NT.2}
```

The fair law is a row product, so

```math
V_{\lambda,C}^{\rm row}\le E_UF.                 \tag{NT.3}
```

Let

```math
{d\Pi_{s,F}\over dU}={e^{sF}\over E_Ue^{sF}},
\qquad
\rho_N^-(\lambda)={1\over mn}
 \sup_{-\lambda\le s\le0}
 E_{\Pi_{s,F}}\sum_a r_a^2.                       \tag{NT.4}
```

The weaker integrated statistic is

```math
\overline\rho_N^-(\lambda)
={1\over\lambda mn}\int_{-\lambda}^0
 E_{\Pi_{s,F}}\sum_a r_a^2\,ds
\le\rho_N^-(\lambda).                            \tag{NT.4a}
```

Only negative tilts appear.  They downweight the high-pressure cap plateau,
so the fixed cap is harmless on this entire window.

All estimates below are uniform in the two actual children and in either
relative orientation.  They may therefore be applied directly in whichever
orientation is target-reaching by Theorem 37.34.

## 2. Exact entropy cancellation

### Theorem NT.1 (negative-tilt product-value theorem)

For every row product `P`, Donsker--Varadhan duality at parameter
`-lambda` gives

```math
E_PF+{1\over\lambda}D(P\Vert U)
\ge-{1\over\lambda}\log E_Ue^{-\lambda F}.       \tag{NT.5}
```

Consequently

```math
\boxed{
0\le E_UF-V_{\lambda,C}^{\rm row}
\le {1\over\lambda}\log E_U
 e^{-\lambda(F-E_UF)}.}                          \tag{NT.6}
```

The fair-cube logarithmic-Sobolev inequality, clipping contraction, and the
exact cavity gradient now imply

```math
\boxed{
0\le E_UF-V_{\lambda,C}^{\rm row}
\le C_{\rm LS}\lambda t^2mn\rho_N^-(\lambda).}   \tag{NT.7}
```

Here `C_LS` is absolute and may be taken equal to one with the present
half-flip normalization.

*Proof.*  The Gibbs variational inequality (NT.5) holds for every law, hence
also after taking the infimum over row products.  Subtract it from `E_UF` to
obtain (NT.6).  Apply the negative-sign version of CT.8 with
`s=lambda`:

```math
\log E_Ue^{-\lambda(F-E_UF)}
\le C_{\rm LS}\lambda^2
 \sup_{-\lambda\le u\le0}E_{\Pi_{u,F}}\Gamma_F.
```

Finally `Gamma_F<=t^2 sum_a r_a^2` proves (NT.7). `square`

Keeping the integral instead of replacing it by its supremum gives the
strictly weaker bound

```math
\boxed{
0\le E_UF-V_{\lambda,C}^{\rm row}
\le C_{\rm LS}\lambda t^2mn
       \overline\rho_N^-(\lambda).}              \tag{NT.7a}
```

The same argument applies when the infimum is restricted to the complete
fixed-degree square carrier, because that carrier contains the uniform row
density.  Uniform truncation RT.2 gives, after enlarging constants,

```math
\boxed{
0\le E_UL-V_\lambda^{\rm row}
\le C_{\rm LS}\lambda t^2mn\overline\rho_N^-(\lambda)
   +O_{\beta,\lambda,K,C}(N^{3/2}e^{-cN}),}      \tag{NT.8}
```

and the analogous bound for the restricted square-carrier value.

On balanced splits `t^2mn=Theta_beta(N)`.  Thus

```math
\overline\rho_N^-(\lambda)=o(1)
\quad\Longrightarrow\quad
E_UL-V_\lambda^{\rm row}=o(N),                  \tag{NT.9}
```

while `bar rho_N^-(lambda)=O(N^(-alpha))` gives the stronger value error
`O(N^(1-alpha))`.  This improves the square-root loss in the uniform
response-range theorem because the objective's entropy penalty cancels
exactly instead of being charged as an external query budget.

## 3. Converse and directed dependence

The converse is immediate from (NT.7): if

```math
E_UF-V_{\lambda,C}^{\rm row}\ge\eta N,           \tag{NT.10}
```

then

```math
\boxed{
\overline\rho_N^-(\lambda)
\ge c_{\eta,\beta,\lambda}>0.}                  \tag{NT.11}
```

Thus the fixed negative-tilt curve is a necessary certificate for a linear
optimized product gain, not merely for a declared pair of carriers.

There is a more direct branch consequence.  Since the fair bridge law is an
admissible product competitor in the reverse information projection,

```math
\mathcal I_\lambda^{\leftarrow}
\le D(U\Vert q_\lambda)
=\lambda(E_UL-V_\lambda).                        \tag{NT.11a}
```

Because `F<=L`, the raw escort gain is no larger than its capped analogue,
apart from `E_UL-E_UF`; the latter is exponentially small by RT.2.  Applying
the integrated log-Sobolev estimate therefore proves:

### Corollary NT.2 (fixed negative-overlap branch separator)

```math
\boxed{
\mathcal I_\lambda^{\leftarrow}
\le C_{\rm LS}\lambda^2t^2mn
       \overline\rho_N^-(\lambda)
 +O_{\beta,\lambda}(N^{3/2}e^{-cN}).}            \tag{NT.11b}
```

Consequently, on balanced splits,

```math
\overline\rho_N^-(\lambda)=o(1)
\Longrightarrow \mathcal I_\lambda^{\leftarrow}=o(N),        \tag{NT.11c}
```

and an extensive reverse projection forces
`overline rho_N^-(lambda)>=c>0`.  Conditional on
`mathcal J>=eta N`, (NT.11c) selects the coherent-retuning alternative

```math
\boxed{
\mathcal J-\mathcal I_\lambda^{\leftarrow}
\ge\eta N-o(N).}                                 \tag{NT.11d}
```

Thus `overline rho_N^-(lambda)` is a concrete carrier-independent analytic
separator between alternatives (ii) and (iii), conditional on the already
isolated linear-`J` phase.  It remains an expectation over the capped Gibbs
law rather than a proved finite-state evaluator.

The fair-base smoothing identity FB.14 gives an exact one-dimensional
presentation of this separator.  With `A_N,p_N,Q_N` as in FB.11--FB.14,
the whole negative interval stays away from the replica-one singularity and

```math
\boxed{
\overline\rho_N^-(\lambda)
={1\over\lambda}\int_{-\lambda}^0
 {A_N(s)-p_N(s)\over s-1}\,ds
 +O_\beta(N^{-1/2})+e^{-cN}.}                   \tag{NT.11e}
```

No mixed derivative at `s=1`, growing tilt window, or adaptive cap is needed.
Thus the SML can equivalently be stated as decay of one integrated ordinary
replicated-smoothing secant on a fixed interval.

Let

```math
V_\lambda=-{1\over\lambda}\log E_Ue^{-\lambda L},
\qquad
\mathcal I_\lambda^{\leftarrow}
=\lambda(V_\lambda^{\rm row}-V_\lambda).         \tag{NT.12}
```

The exact decomposition AC.24 is

```math
E_UL-V_\lambda
=E_UL-V_\lambda^{\rm row}
 +{1\over\lambda}\mathcal I_\lambda^{\leftarrow}.             \tag{NT.13}
```

The same integrated estimate bounds the fair escort gain itself.  Hence a
linear value in the left side of (NT.13) is incompatible with
`bar rho_N^-(lambda)=o(1)`; (NT.13) supplies no additional live branch in
that regime.  The useful implication is NT.11c--NT.11d.

There is also an exact translation of the coherent-retuning branch.  If `r`
is the canonical product of the actual escort row marginals and

```math
\mathcal G_L(r)=E_rL+{1\over\lambda}D(r\Vert U),
```

then the established identity is

```math
\mathcal J-\mathcal I^{\leftarrow}
=\lambda\{\mathcal G_L(r)-V_\lambda^{\rm row}\}.               \tag{NT.15}
```

Consequently `bar rho_N^-(lambda)=o(1)` gives

```math
\boxed{
\mathcal J-\mathcal I^{\leftarrow}
=\lambda\{\mathcal G_L(r)-E_UL\}+o(N).}          \tag{NT.16}
```

Thus, on the small-overlap branch, the global product optimization disappears
to `o(N)` and, conditional on linear `J`, alternative (iii) is forced.  If
the overlap stays positive, NT.11 only certifies the possibility of a linear
optimized product gain.  It does not by itself show that the gain occurs or
identify its branch.

## 4. The cap is unnecessary on the negative path

The preceding capped argument interfaces directly with the full-carrier
response theorem, but the optimized-value separator has a still simpler
form.  Since `L>=0`, every tilt in `[-lambda,0]` suppresses the high-pressure
tail.  Define the **raw** path

```math
{d\widehat\Pi_s\over dU}={e^{sL}\over E_Ue^{sL}},
\qquad
\widehat\rho_N^-(\lambda)
={1\over\lambda mn}\int_{-\lambda}^0
 E_{\widehat\Pi_s}\sum_a r_a^2\,ds.              \tag{NT.17}
```

### Theorem NT.3 (raw negative-disorder branch separator)

Without any truncation error,

```math
\boxed{
\begin{aligned}
0\le E_UL-V_\lambda^{\rm row}
 &\le C_{\rm LS}\lambda t^2mn
       \widehat\rho_N^-(\lambda),\\
\mathcal I_\lambda^{\leftarrow}
 &\le C_{\rm LS}\lambda^2t^2mn
       \widehat\rho_N^-(\lambda).
\end{aligned}}                                   \tag{NT.18}
```

*Proof.*  Apply Donsker--Varadhan directly to `L`, and put

```math
\psi(u)=\log E_Ue^{-u(L-E_UL)}.
```

The fair-cube logarithmic-Sobolev inequality and the exact, **unclipped**
cavity derivative give

```math
{d\over du}{\psi(u)\over u}
\le C_{\rm LS}E_{\widehat\Pi_{-u}}\Gamma_L
\le C_{\rm LS}t^2E_{\widehat\Pi_{-u}}\sum_a r_a^2.
```

Integrate from zero to `lambda`.  The first line follows as in NT.6; the
second follows from
`I^leftarrow<=D(U||q_lambda)=psi(lambda)`. `square`

At the endpoint, `widehat Pi_(-lambda)=q_lambda` is exactly the actual
negative-disorder Gibbs law.  The rare rank-one high-pressure atoms that
destroy positive all-tilt concentration are automatically downweighted and
never enter the argument.

There is also a cap-free smoothing representation.  For `s in [-lambda,0]`
define

```math
\widehat{\mathscr H}_s(u)
={1\over s}\log E_{B\sim U,G}e^{sL(B+\sqrt uG)},
\qquad
\widehat A_N(s)={2\widehat{\mathscr H}_s'(0+)
                  \over t^2mn},                 \tag{NT.19}
```

with the expectation-valued extension at zero.  The raw heat identity is

```math
\widehat A_N(s)
=1+(s-1)E_{\widehat\Pi_s}\mathcal O(B).          \tag{NT.20}
```

Consequently

```math
\boxed{
\widehat\rho_N^-(\lambda)
={1\over\lambda}\int_{-\lambda}^0
 {\widehat A_N(s)-1\over s-1}\,ds
 +O_\beta(N^{-1/2}).}                            \tag{NT.21}
```

All exponential moments and differentiated integrands are dominated on this
fixed nonpositive interval because `e^(sL)<=1`, the Gibbs first derivatives
are bounded by `t`, and their second derivatives by `t^2`.  Thus no cap,
tail correction, mixed derivative, or adaptive window appears in NT.21.

## 5. Revised smallest lemma

For the balanced product phase, the required statement is now only

> **raw fixed negative-tilt overlap classification:** prove
> `hat rho_N^-(lambda)=o(1)` (with a power rate for a summable recurrence),
> which selects coherent retuning whenever `J` is linear; or prove that
> positive `hat rho_N^-` yields an explicit favorable reverse-product
> direction, or construct an actual minimizing sequence showing why it does
> not.

This is strictly weaker than full response-range control, all-fixed-window
rigidity, a mixed derivative at replica one, or a growing-window estimate.
Existing child minimality supplies no negative-replica cross-edge comparison
at the required precision.
