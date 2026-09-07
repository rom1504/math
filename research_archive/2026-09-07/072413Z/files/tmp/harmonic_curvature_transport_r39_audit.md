# Independent audit of Wave 39 harmonic curvature transport

## Verdict

**PASS.**  No sign, factor, conditioning, or bootstrap correction is needed.
The original checker was rerun under `.venv` and passes on all four displayed
finite cases.  The conclusions remain conditional: endpoint vertex cost,
restoring, orientation, and selector Hellinger are separate assumptions, and
no estimate of the curvature isotonic defect is proved.

## Identity checks

### H6

For each context, `k(t)=int_0^t s A''(s) ds`.  Therefore

```math
C_V(1)=\int_0^1sE_{\mu_1}V_sds,
\qquad
W_V=\int_0^1sE_{\mu_s}V_sds.
```

The already verified boundary identity is
`W_V=C_V(1)-int Cov_mu_t(g,L_t)dt`.  Hence

```math
\int_0^1Cov_{\mu_t}(g,L_t)dt
=\int_0^1s(E_{\mu_1}V_s-E_{\mu_s}V_s)ds,
```

with exactly the sign shown in H6.

### Conditioning on `g`

For every two interpolation times,

```math
d\mu_t/d\nu=exp(tg-psi(t)).
```

This density is constant on each exact `g`-level.  Thus the conditional law
of the full state given `g` is genuinely independent of time.  Although
`V_s` itself depends on `s`, for each fixed `s` it is one state function and
the same conditional regression `m_s(g)` is valid under both `mu_s` and
`mu_1`.  No conditioning assumption is hidden here.

### H12 and H14

Under the balanced mixture, the posterior time probability is
`eta=p_1/(p_s+p_1)`.  Direct expansion gives

```math
Cov_{(lambda_s+lambda_1)/2}(m_s,eta)
=\frac14(E_{\lambda_1}m_s-E_{\lambda_s}m_s),
```

so the factor `4` in H12 is correct.  Binary Pinsker gives
`Var(eta)<=I_s/2`.  Removing a nondecreasing isotonic fit can only hurt the
negative covariance because `eta` is nondecreasing.  Multiplication by `4`
then gives `sqrt(8 I_s Rbar_s)`, the exact factor in H14.

### H16 and H17

Equal-weight Jensen--Shannon information obeys

```math
I_s\le\frac14[D(mu_1||mu_s)+D(mu_s||mu_1)].
```

For the exponential family, the Jeffreys sum is
`(1-s)(psi'(1)-psi'(s))`, proving H16.  Since

```math
int_s^1 Var_{mu_u}(g)du<=mathscr H/s,
```

the adverse part of H6 is at most

```math
sqrt(2 mathscr H)
int_0^1sqrt{s(1-s) Rbar_s}ds.
```

The endpoint `s=0` is harmless because it is measure zero and the combined
integrand has the displayed finite form.  H17 has the correct factor and
time weight.

### H18--H20 bootstrap

Let `W_V,W_0` be weighted vertex and orientation energies.  H17 and the
endpoint assumption give

```math
W_V<=O(a_n)+sqrt(2 mathscr H) K.
```

The separately assumed restoring comparison gives
`mathscr H<=C_R(W_V+W_0)`.  With `W_0=O(a_n)`, setting
`h=sqrt(mathscr H)` yields

```math
h^2<=O(a_n)+O(K)h.
```

Thus `K^2=O(a_n)` implies `mathscr H=O(a_n)`.  This is a valid quadratic
bootstrap, not a circular use of the desired entropy estimate.  It does rely
essentially on the unproved restoring and orientation hypotheses, which the
memo labels explicitly.

## Numerical scope

The checker uses floating grouping of `g`-levels and Gauss--Legendre
quadrature.  Its residuals verify the algebra numerically but do not provide
certified signs or asymptotic evidence.  The memo labels the finite
monotonicity failures and magnitudes numerical, so no scope correction is
needed.
