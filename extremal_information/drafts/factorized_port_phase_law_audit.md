# Independent audit: factorized port phase law

**Verdict.** Pass after two notation/scope corrections incorporated into the
main draft: the majority tie rule is fixed before invoking the common
selector witness, and the factor defect in (FP.15) is consistently denoted
`delta_t`.

## 1. Support collapse

After gauging by the base pole, the row law is

```math
(1,X_1,\ldots,X_L),
```

with independent factor blocks. For fixed endpoint signs, the centered
block sum has variance at most `q_t^2`. Therefore

```math
0\le E|S|-|ES|\le\sqrt{\sum_tq_t^2}.
```

The maximum absolute mean is exactly
`1+sum_(t,j)|mu_(t,j)|`. This proves both uniform bounds (FP.10)--(FP.11),
including the normalization by `p_L`.

## 2. Relative defect and selector witness

Fixing one antipodally odd tie-broken majority selector at each arity makes
its Walsh expansion supported on odd products. Every such product belongs
to the full Cartesian pole span `U^(L)`, so the relative defect (FP.14)
bounds its quadratic loss by `e_L r_L n_L/2`. The field contribution is
exactly `m_L||z_epsilon||_1`. Thus (FP.21), followed by the support-collapse
bound, gives (FP.17)--(FP.18) with no missing factor of two.

For factor compressions, the existing Cartesian relative-synchronization
theorem gives

```math
e_L\le\sum_t\delta_t.
```

When all compressed contractions are positive semidefinite and the
represented span is the full tensor product, the smallest tensor eigenvalue
is `prod_t(1-delta_t)`, so in fact

```math
e_L=1-\prod_t(1-\delta_t).
```

The draft correctly does not infer `e_L->0` from ordinary summability of a
fixed sequence of nonzero factor defects.

## 3. Exact-sign completion

Trace zero makes hollowing invisible because

```math
x^T(H-\operatorname {diag}H)x=x^THx-\operatorname {tr}H=x^THx.
```

There are `d_L^(aux)=m_Lp_L` new vertices. Any hollow exact-sign fill on
them changes the cap by at most `binom(d_L^(aux),2)`. Under (FP.22) and the
stated squared-cost hypothesis, this is `o(r_Ln_L)`, while

```math
{r_Ln_L\over(n_L+d_L^{aux})^{3/2}}\longrightarrow\rho.
```

This verifies (FP.24).

## 4. Seed and obstruction

For the order-16 seed, the relative generator means are `1/2` and `0`.
With both generators at every factor,

```math
\theta_L={1+L/2\over2L+1}\to1/4,
```

and `m_Lp_L/r_L->1`, giving the completed-parent limit `3/4`.

With only the mean-`1/2` generator in `X` factors and the mean-zero
generator in `Y` factors,

```math
\theta_L={1+K_L/2\over L+1}.
```

Alternating blocks that dominate the entire preceding prefix make
`K_L/L` tend to one and zero on alternating endpoints. Every pole remains
an exact positive tensor pole, so `e_L=0`; the completed cap consequently
has subsequential limits `1` and `1/2`. This is a signing-sequence
counterexample, not a statement about minima.

The obstruction is distinct from the Walsh-prefix mantissa law. It occurs
at tensor endpoints through a nonconvergent empirical factor phase; the
Walsh-prefix obstruction occurs between tensor endpoints for one fixed
factor hierarchy.

The companion verifier passed all 44 finite checks during this audit.
