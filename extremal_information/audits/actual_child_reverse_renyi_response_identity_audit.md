# Audit of the actual-child reverse-Renyi response identity

Status: **independent adversarial audit; passed after one substantive
localization correction**.

Audited source:
`drafts/actual_child_reverse_renyi_response_identity.md`.

## 1. Separate-coordinate derivation

Give the channel coordinates separate amplitudes `u_f`, write
`rho_f=tanh u_f`, and hold every `u_f`, `f ne e`, fixed.  Deleting `e`
gives

```math
p(B)=p_{-e}(B_{-e})(1+\rho_e B_e r_e(B_{-e})).
```

In particular `r_e` has no `u_e` derivative.  If `a=rho_e r_e` and
`z=arctanh a`, summing the two inverse-escort weights over `B_e` gives

```math
g_e(a)={ (1+a)^{-\lambda}+(1-a)^{-\lambda}\over2}
       =(1-a^2)^{-\lambda/2}\cosh(\lambda z).
```

Since `partial_(u_e)a=(1-rho_e^2)r_e`,

```math
{1\over\lambda}\partial_{u_e}\log g_e(a)
={(1-\rho_e^2)r_e
  [\rho_e r_e+\tanh(\lambda\operatorname{arctanh}(\rho_e r_e))]
  \over1-\rho_e^2r_e^2}.
```

The marginal of the full inverse escort on `B_{-e}` is exactly the measure
obtained after multiplying by `g_e`.  Summing these partial derivatives and
then setting all `u_e=u` proves (RR.7).  Thus no derivative of the extrinsic
response was omitted: such a derivative would arise only if one incorrectly
differentiated all channel amplitudes while treating a single-coordinate
partial derivative.

The sign also checks directly from

```math
\mathcal R_\lambda'(u)
=-\mathbb E_{q_{\lambda,u}}\partial_u\log p_u(B).
```

For the one-point prior `Q=+1`, the formula reduces per coordinate to
`tanh u+tanh(lambda u)`, which is the derivative of
`log cosh u+lambda^(-1)log cosh(lambda u)`.

## 2. Constants and Renyi orientation

For `a=rho r`, the ratio
`tanh(lambda arctanh a)/a` lies in
`[min{lambda,1},max{lambda,1}]`.  Also

```math
1-\rho^2
\le {1-\rho^2\over1-\rho^2r^2}\le1.
```

These two observations give (RR.9) with precisely
`c_lambda=1+min{lambda,1}` and
`C_lambda=1+max{lambda,1}`; integration gives (RR.10).

The divergence orientation is correct:

```math
D_{1+\lambda}(U\Vert\Pi_u)
={1\over\lambda}\log
 \mathbb E_{\Pi_u}(dU/d\Pi_u)^{1+\lambda}
={1\over\lambda}\log\mathbb E_U p_u^{-\lambda}.
```

At the physical amplitude, conditioning the joint output law (2.4) on the
chosen orientation gives `p_t=c_epsilon exp(L_epsilon)` for a constant
independent of `B`.  Hence `L_t=C_t+log p_t`,
`V_lambda=C_t-R_lambda(t)`, and
`E_U L_t-V_lambda=R_lambda(t)-D(U||Pi_t)`.  No sector-normalization factor is
missing.

## 3. Density scaling and the correction

The original draft correctly obtained a point with dense `S_lambda` from
the upper half of (RR.10), but did not locate that point away from `u=0`.
Without such localization, (RR.23) did not follow: (RR.21) loses a factor
`rho(u)^2`, and the selected amplitude could formally have been
`o(N^{-1/2})`.

The draft has been corrected.  If `R_lambda(t)>=eta N`, set

```math
\alpha=\min\{1/2,\sqrt{2\eta/(C_\lambda\beta^2)}\},
\qquad t=\beta/\sqrt N.
```

Then `mn log cosh(alpha t)<=eta N/(4C_lambda)`, so

```math
\int_{\alpha t}^t\tanh u\,S_\lambda(u)du
\ge {3\eta N\over4C_\lambda}.
```

Consequently some `u_N in [alpha t,t]` has

```math
S_\lambda(u_N)
\ge {3\eta N\over4C_\lambda\log\cosh t}.
```

For balanced splits this is a fixed positive density of the `mn`
coordinates.  Since the Bernoulli information integrand in (RR.22) is at
least `s^2/2`, (RR.20)--(RR.21) now give

```math
\sum_e I(B_e;B_{-e})
\ge {\min\{\lambda,1\}^2\over2}
     \tanh^2(\alpha t)S_\lambda(u_N)
=\Omega(N).
```

Thus the repaired (RR.23) is rigorous.  Conversely, uniform
`S_lambda(u)/(mn)=o(1)` makes (RR.10) `o(N)`, as claimed.

## 4. Escort bias and fairness

Conditioning the inverse escort on `B_{-e}` gives weights
`(1+rho b r)^{-lambda}`.  Their signed mean is exactly

```math
s_{e,u}=-\tanh(\lambda\operatorname{arctanh}(\rho r_{e,u})),
```

so (RR.20)--(RR.21) have the correct sign, factor, and orientation.

For the actual sector, global inversion of the bridge is a symmetry:
equivalently, the zero-bridge prior on `Q` is centrally symmetric (global
spin inversion sends `Q` to `-Q` without changing the child energies).
Therefore `p_u(-B)=p_u(B)`, hence also `q_{lambda,u}(-B)=q_{lambda,u}(B)`.
Each bridge bit is marginally fair and (RR.22) is valid.  This fairness is
not asserted for an arbitrary non-central prior; it is used only for the
actual child prior.

## 5. Scope verdict

The exact identity is not circular and does not differentiate an extrinsic
quantity along the wrong path.  It classifies a necessary resource: linear
reverse-Renyi compensation forces dense weak posterior response, while
uniformly vanishing response density rules it out.  Presence of dense
response at one isolated amplitude is not by itself sufficient for linear
work, so the opening language was corrected from “equivalent” to the
one-way implication actually proved.  Finally, the result is leave-one-bit,
not cross-row: it does not yet control the directed row-product excess.

## 6. Addendum: the reverse-KL endpoint

The subsequently added (RR.15a)--(RR.15b) also pass audit.  As
`lambda downarrow 0`,

```math
{1\over\lambda}\log\mathbb E_U p_u^{-\lambda}
\longrightarrow-\mathbb E_U\log p_u=D(U\Vert\Pi_u),
\qquad q_{\lambda,u}\longrightarrow U.
```

In (RR.7), `tanh(lambda arctanh(rho r))` tends to zero, leaving exactly

```math
\mathcal R_0'(u)=\sum_e\mathbb E_U
 {\rho(1-\rho^2)r_{e,u}^2\over1-\rho^2r_{e,u}^2}.
```

This agrees with a direct fair-bit average of
`-partial_u log(1+rho B_e r_e)`.  Since both `R_lambda(0)` and `R_0(0)`
vanish, the fundamental theorem of calculus gives (RR.15b) with no missing
endpoint constant.  The warning that the difference of derivatives need
not be pointwise nonnegative is necessary and correct.
