# Wave 39 harmonic attack: Bregman curvature transport and scalar time information

## Status

- **Verified:** integrated signed vertex migration is exactly the transport of
  the same conditional Bregman-curvature profile from the contemporaneous law
  `mu_s` to the endpoint law `mu_1`.
- **Verified:** the interpolation parameter is conditionally independent of
  the full state given the scalar score `g`.  Migration therefore depends
  only on a one-dimensional regression of total local curvature on `g`, with
  no loss from fluctuations inside a `g`-level.
- **Verified sufficient package:** a scalar isotonic-defect bound at scale
  `n^(1/2-2c)`, together with the already separate endpoint, restoring, and
  orientation assumptions, closes a quadratic bootstrap and proves the
  parent entropy target.
- **Open:** no exact-minimizer estimate of the new scalar curvature defect is
  proved.  Endpoint vertex KL, restoring, orientation, and adjacent-selector
  Hellinger control remain separate inputs.
- **Numerical:** exact monotonicity is false as a finite shortcut.  The new
  identities and inequalities pass on `A_4,A_6,A_8,A_9`; finite magnitudes
  are reported only as mechanism diagnostics.

The checker is `tmp/harmonic_curvature_transport_r39_check.py`.

## 1. Conditional KL is integrated Bregman curvature

Retain the nonreference vertex coordinates `V_*`.  For coordinate `i` and
context `e`, let

```math
A_{i,e}(t)
=\log\mathbb E_{\nu(D_i\mid e)}e^{t g},
\qquad
v_{i,s}(e)=A''_{i,e}(s)
=\operatorname{Var}_{\mu_s(D_i\mid e)}g.
\tag{R39.H1}
```

The conditional binary KL is the Bregman divergence

```math
k_{i,e}(t)=tA'_{i,e}(t)-A_{i,e}(t),
```

so

```math
\boxed{
k_{i,e}(t)=\int_0^t s\,v_{i,s}(e)\,ds.
}
\tag{R39.H2}
```

Define the statewise total instantaneous vertex curvature

```math
V_s(d)=\sum_{i\in V_*}v_{i,s}(d_{-i}).
\tag{R39.H3}
```

At the endpoint, Fubini and (R39.H2) give

```math
\boxed{
C_V(1)=\mathbb E_{\mu_1}L_1
=\int_0^1s\,\mathbb E_{\mu_1}V_s\,ds.
}
\tag{R39.H4}
```

At the contemporaneous law, conditional variance gives

```math
\boxed{
\int_0^1s\mathcal E_{s,V}(g)\,ds
=\int_0^1s\,\mathbb E_{\mu_s}V_s\,ds.
}
\tag{R39.H5}
```

Subtracting (R39.H5) from (R39.H4), or comparing directly with (10.1082),
proves the exact curvature-transport identity

```math
\boxed{
\int_0^1\operatorname{Cov}_{\mu_t}(g,L_t)\,dt
=\int_0^1s\left{
\mathbb E_{\mu_1}V_s-\mathbb E_{\mu_s}V_s
\right}ds.
}
\tag{R39.H6}
```

This is not a new bound by itself.  It replaces the accumulated KL load by
the instantaneous conditional Fisher/Bregman curvature and moves all time
dependence into a comparison of two state laws.

## 2. The time comparison is exactly one-dimensional

Put

```math
\psi(t)=\log\mathbb E_\nu e^{tg}.
```

For `0<=s<1`,

```math
\frac{d\mu_1}{d\mu_s}(d)
=\exp\{(1-s)g(d)-\psi(1)+\psi(s)\}.
\tag{R39.H7}
```

The likelihood ratio depends only on `g`.  Consequently the conditional law
of `D` given `g(D)=u` is independent of the interpolation time.  Define the
time-independent conditional regression

```math
m_s(u)=\mathbb E_\nu[V_s(D)\mid g(D)=u].
\tag{R39.H8}
```

If `lambda_t` is the law of `g(D)` under `mu_t`, then

```math
\mathbb E_{\mu_t}V_s=\mathbb E_{\lambda_t}m_s,
\qquad t\in\{s,1\}.
\tag{R39.H9}
```

Thus within-level fluctuations of `V_s` cancel exactly; unlike (10.1085), no
approximation of the full statewise statistic is required.

There is an equivalent information statement.  Give a time label `J` equal
to `s` or `1` with probability `1/2`, and conditionally sample
`D sim mu_J`.  Then

```math
\boxed{
I(J;D\mid g(D))=0,
\qquad I(J;D)=I(J;g(D)).
}
\tag{R39.H10}
```

In other words, `g` is a sufficient statistic for distinguishing the two
interpolation times.

Equation (R39.H7) also says that `lambda_1` is an increasing likelihood-ratio
tilt of `lambda_s`, hence it stochastically dominates `lambda_s`.  Let
`Q_t(a)` be the left quantile of `lambda_t`.  The canonical quantile coupling
satisfies `Q_1(a)>=Q_s(a)` and gives the exact signed scalar formula

```math
\boxed{
\int_0^1\operatorname{Cov}_{\mu_t}(g,L_t)\,dt
=\int_0^1\int_0^1
s\{m_s(Q_1(a))-m_s(Q_s(a))\}\,da\,ds.
}
\tag{R39.H11}
```

This is an exact one-dimensional reformulation of (10.1083), preserving the
sign and all cancellation across coordinates and time.  It does not assert
that `m_s` is pointwise nondecreasing.

## 3. A scalar approximate-order inequality

Let

```math
\overline\lambda_s=\frac12(\lambda_s+\lambda_1),
```

and let

```math
\eta_s(u)=\Pr(J=1\mid g(D)=u).
```

The posterior `eta_s` is nondecreasing in `u`.  Directly from the two mixture
components,

```math
\boxed{
\mathbb E_{\mu_1}V_s-\mathbb E_{\mu_s}V_s
=4\operatorname{Cov}_{\overline\lambda_s}(m_s,\eta_s).
}
\tag{R39.H12}
```

Define the scalar isotonic defect

```math
\overline{\mathcal R}_s
=\inf_{\phi\text{ nondecreasing}}
\mathbb E_{\overline\lambda_s}[m_s-\phi]^2.
\tag{R39.H13}
```

Let

```math
\mathcal I_s=I(J;D)=I(J;g(D))
```

for the balanced time experiment.  Since `eta_s` and every admissible `phi`
are comonotone,

```math
\operatorname{Cov}(\phi,\eta_s)\ge0.
```

Furthermore binary Pinsker gives

```math
\operatorname{Var}_{\overline\lambda_s}(\eta_s)
=\mathbb E(\eta_s-1/2)^2\le\mathcal I_s/2.
```

Applying Cauchy--Schwarz only to the isotonic residual in (R39.H12) yields

```math
\boxed{
\left[
\mathbb E_{\mu_s}V_s-\mathbb E_{\mu_1}V_s
\right]_+
\le\sqrt{8\mathcal I_s\overline{\mathcal R}_s}.
}
\tag{R39.H14}
```

This is a scalar, time-aware successor to (10.1084).  It exactly discards
variation of total curvature inside a `g`-level and weights the surviving
isotonic defect by how distinguishable time `s` is from the endpoint.

For comparison, the quantile coupling also gives the direct sufficient
pressure

```math
\int_0^1\int_0^1
s[m_s(Q_s(a))-m_s(Q_1(a))]_+\,da\,ds.
\tag{R39.H15}
```

Bounding (R39.H15) at the project scale proves (10.1083), but it discards
more cancellation than the exact signed formula (R39.H11).

## 4. Bregman time information closes a restoring bootstrap

The time information in (R39.H14) has a global Bregman upper bound.  Convexity
of relative entropy in its second argument gives

```math
\mathcal I_s
\le\frac14\{D(\mu_1\Vert\mu_s)+D(\mu_s\Vert\mu_1)\}.
```

The two KL divergences are the opposite Bregman divergences of `psi`, so

```math
\boxed{
\mathcal I_s
\le\frac{1-s}{4}\{\psi'(1)-\psi'(s)\}
=\frac{1-s}{4}\int_s^1
\operatorname{Var}_{\mu_u}(g)\,du.
}
\tag{R39.H16}
```

Let the parent entropy be

```math
\mathscr H
=\operatorname{Ent}_\nu(f)
=\int_0^1u\operatorname{Var}_{\mu_u}(g)\,du.
```

For `s>0`, the final integral in (R39.H16) is at most `mathscr H/s`.
Combining (R39.H6), (R39.H14), and (R39.H16) proves

```math
\boxed{
\left[-\int_0^1
\operatorname{Cov}_{\mu_t}(g,L_t)\,dt\right]_+
\le\sqrt{2\mathscr H}\,\mathcal K,
\qquad
\mathcal K
=\int_0^1\sqrt{s(1-s)\overline{\mathcal R}_s}\,ds.
}
\tag{R39.H17}
```

This yields a concrete sufficient package.  Put
`a_n=n^(1/2-2c)` and assume uniformly over the relevant project pairs:

```math
\boxed{
C_V(1)=O(a_n),
\qquad
\int_0^1t\mathcal E_{t,\mathrm{orientation}}(g)\,dt=O(a_n),
\qquad
\operatorname{Var}_{\mu_t}(g)
\le C_R\mathcal E_t(g),
\qquad
\mathcal K^2=O(a_n).
}
\tag{R39.H18}
```

The third clause is precisely the bounded restoring comparison already
separated in (10.983), not a new unconditional fact.  To see that (R39.H18)
closes, let `W_V` and `W_0` be the weighted vertex and orientation energies.
Equations (10.1082) and (R39.H17) give

```math
W_V\le C_V(1)+\sqrt{2\mathscr H}\,\mathcal K.
```

The restoring clause gives

```math
\mathscr H\le C_R(W_V+W_0).
```

Writing `h=sqrt(mathscr H)` therefore yields a quadratic inequality

```math
h^2\le O(a_n)+O(\mathcal K)h.
```

Since `mathcal K^2=O(a_n)`, it follows that

```math
\boxed{\mathscr H=O(a_n).}
\tag{R39.H19}
```

Thus (R39.H18), together with the still separate adjacent-selector Hellinger
input, proves the parent part of the harmonic convergence package.

The new migration-side approximate-order target is

```math
\boxed{
\left(
\int_0^1\sqrt{s(1-s)\overline{\mathcal R}_s}\,ds
\right)^2
=O(n^{1/2-2c}).
}
\tag{R39.H20}
```

No estimate of (R39.H20) is currently known.

## 5. Assumptions, falsification, and finite audit

Nothing above controls the endpoint cost `C_V(1)`.  The endpoint selector-
erasure identities remain possible inputs, but endpoint erasure alone was
already shown not to control transient migration.  Likewise, (R39.H20) does
not prove restoring, orientation, or selector Hellinger control.

An unbounded actual-minimizer family with `mathcal K^2=omega(a_n)` falsifies
the specific approximate-order implementation (R39.H20), not the exact
signed target: (R39.H14) takes an adverse part at each time.  The exact scalar
identity (R39.H11) is equivalent to migration, so an integrated negative part
of order `omega(a_n)` there is an actual migration obstruction.

The checker uses 96-point Gauss--Legendre quadrature.  Direct covariance,
curvature transport, the endpoint Bregman integral, and the energy boundary
identity agree below `3e-15`.  Regression uses `g`-levels merged only across
floating symmetry roundoff.

```text
case           signed migration   quantile downward   isotonic-info bound   K^2
A4,.5,m3        0.000006998        0.000002170          0.000003624       2.96e-9
A6,.5,m3        0.000005609        0.000001685          0.000007230       7.57e-8
A8,2,m5        -0.007441065        0.013519670          0.013694017       6.70e-4
A9,2,m7        -0.001467944        0.009703035          0.005469408       9.98e-4
```

The nonzero downward pressures in the favorable `A_4/A_6` cases and the
negative `A_8/A_9` signs rule out assuming pointwise coordinate or regression
monotonicity.  On `A_9`, the exact signed magnitude is much smaller than both
one-sided sufficient bounds.  These are finite numerical warnings, not
asymptotic evidence against (R39.H20).

## Frontier

Signed harmonic migration is exactly the endpoint transport of total local
Bregman curvature.  The exponential interpolation makes `g` a sufficient
one-dimensional statistic for that transport, and the remaining structured
quantity is the isotonic defect of the scalar regression `m_s(g)`.  Time
mutual information and the restoring comparison reduce the required scale to
(R39.H20) through a closed quadratic bootstrap.

The next harmonic attack should use the matched principal-restriction form of
the exclusion probabilities to estimate this **summed curvature regression**.
It should not split external partition from exclusion shape, impose a
pointwise covariance sign, or return to coordinatewise adverse parts.
