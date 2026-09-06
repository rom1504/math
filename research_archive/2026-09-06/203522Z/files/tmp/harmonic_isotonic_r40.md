# Wave 40 harmonic isotonic-defect audit

## Exact local formulas

On an `i`-edge with score gap `delta=g(y)-g(x)`, base log odds
`omega=log(nu(y)/nu(x))`, and
`p_s=logistic(omega+s delta)`,

```math
v_{i,s}=A''_{i,e}(s)=p_s(1-p_s)\delta^2,
\qquad
\partial_s v_{i,s}=A'''_{i,e}(s)
=p_s(1-p_s)(1-2p_s)\delta^3.
```

The latter is the third conditional central moment.  Because the conditional
law of `D` given `g(D)=u` is independent of interpolation time,

```math
\partial_s m_s(u)
=\mathbb E_\nu\left[\sum_i A'''_{i,D_{-i}}(s)\mid g=u\right].
```

For two independent draws `B,B'` from the time-`s` conditional bit law,

```math
v_{i,s}(e)
=\frac12\mathbb E[(g(e,B)-g(e,B'))^2].
```

At a matched restriction endpoint, the exact exclusion identity is
`delta=log(r_i(x)/r_i(y))`.  All these formulas are checked by
`tmp/harmonic_isotonic_r40_check.py`.

## Literal monotonicity is false on actual finite minimizers

The checker scans 48 Gauss--Legendre times for several `A_6,A_8,A_9`
restriction endpoints.  Literal monotonicity of `m_s` fails broadly (usually
at all 48 nodes).  The weighted isotonic defect can nevertheless be small
because the inverted score levels may have very small scalar-mixture mass.
The largest scanned ratio was on `A_9,beta=2,m=7`:

```text
K^2 = 0.0009978621442,
C_V(1) = 0.05086819103,
K^2/C_V(1) = 0.01961662335.
```

These finite decimals are numerical only.  They rule out pointwise
monotonicity as a proof target but do not disfavor a signing-specific
weighted estimate.

## Exact matched-erasure wall to endpoint or Dirichlet control of `K`

Let the state space be `{0,1}^3`, let `nu` be uniform, and, up to the
irrelevant normalizing additive constant, put

```math
g_L(x)=L\mathbf1_{\{x_1\vee x_2=1\}}.
```

The third bit is inert.  This is not merely an arbitrary exponential family:
it has an exact positive fixed-size erasure-mixture lift.  Fix
`0<epsilon<1/2` and set

```math
r_1=r_2=\epsilon e^{-g_L},
\qquad
r_3=1-2\epsilon e^{-g_L},
\qquad
b_i=e^{g_L}r_i.
```

Then `sum_i r_i=1=n-m`, `b_1=b_2=epsilon`, and
`b_3=e^{g_L}-2epsilon`.  Each `b_i` is independent of bit `i`, all are
strictly positive, and `sum_i b_i=e^{g_L}`.  Thus the `b_i` are valid
positive likelihood components indexed by the selector omitting coordinate
`i`, and their posterior omission probabilities are exactly the `r_i`.
Moreover, on every `i`-edge,

```math
g_L(y)-g_L(x)=\log\frac{r_i(x)}{r_i(y)}.
```

Hence the example satisfies component invariance, fixed-size posterior sum,
and the matched exclusion identity exactly.

Put

```math
q_s=\frac{e^{sL}}{1+e^{sL}},
\qquad
v_s=L^2q_s(1-q_s).
```

Only the two active edges out of the all-zero active state have nonzero
curvature.  Conditional on the two score levels,

```math
m_s(0)=2v_s,
\qquad
m_s(L)=\frac23v_s.
```

Thus literal monotonicity fails for every `s`.  If

```math
a_s=\lambda_s\{0\}=\frac1{1+3e^{sL}},
\qquad
w_s=\frac{a_s+a_1}{2},
```

the two-level isotonic projection pools both values and gives exactly

```math
\overline{\mathcal R}_s
=w_s(1-w_s)\left(\frac43v_s\right)^2.
```

Consequently

```math
\frac{\mathcal K_L}{\sqrt L}
\longrightarrow
\kappa
=\frac43\int_0^\infty
\frac{e^u}{(1+e^u)^2}
\sqrt{
u\,\frac1{2(1+3e^u)}
\left(1-\frac1{2(1+3e^u)}\right)
}\,du>0.
```

In particular `K_L^2=Theta(L)`.  This can also be certified without
dominated convergence by restricting the integral to `u=sL in [1,2]`.

The endpoint vertex cost, however, is

```math
C_V(1)
=2\frac{1+e^L}{1+3e^L}
D(\operatorname{Ber}(q_1)\Vert\operatorname{Ber}(1/2))
\longrightarrow\frac23\log2.
```

The parent entropy tends to `log(4/3)`.  Even the exact restoring estimate is
uniformly favorable:

```math
\frac{\operatorname{Var}_{\mu_s}(g_L)}
     {\mathbb E_{\mu_s}V_s}
=\frac{3(1+e^{sL})}{2(1+3e^{sL})}
\le\frac34.
```

Also `int_0^1 s E_{mu_s}V_s ds=O(1)` and the actual signed migration tends
to the finite negative constant

```math
-\frac43\int_0^\infty
\frac{u e^u}{(1+e^u)^2(1+3e^u)}\,du.
```

Numerically this is about `-0.05133578656`, whereas
`K_L^2/L` approaches about `0.0152`.  Thus the scalar isotonic--information
majorant loses a growing factor even though the actual adverse migration is
bounded.

This exactly falsifies any universal bound of `K^2` by endpoint vertex cost,
parent entropy, the time-weighted scalar Dirichlet integral, or those
quantities together with restoring, even after imposing the fixed-size
matched-erasure identities.  The third-moment derivative does detect the
wall, but at the cost

```math
\int_0^1 |v_s'|\,ds=v_0-v_1=\Theta(L^2),
```

so it is not controlled by the existing budgets.

The wall is **not** a quadratic-signing endpoint.  Its positive component
likelihoods are arbitrary omitted-coordinate functions, not the external
completion partition functions of restrictions of one complete signing;
the uniform base law is not the nonzero-temperature Gibbs law of an exact
quadratic-signing minimizer; and no order-minimal signing is supplied.
Therefore it does not falsify (10.1096) for actual minimizers.  It proves
that any positive theorem must use genuinely signing-specific
outside-completion/minimality structure (or replace the lossy isotonic
majorant by a sharper signed-migration argument), not only endpoint cost,
restoring, component invariance, fixed-size posterior mass, and matched
exclusion.
