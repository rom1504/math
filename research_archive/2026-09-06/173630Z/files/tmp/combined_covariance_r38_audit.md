# Independent audit of `combined_covariance_r38`

## Result

**PASS, with no mathematical correction required.**  The claims H2--H17
have the stated scope.  In particular, all of H2--H6 are for the
`n-1` nonreference vertex coordinates only; the global-orientation edge is
not silently included.  The memo explicitly retains endpoint cost,
restoring comparison, orientation, and selector Hellinger as separate inputs.

The original checker was rerun with `.venv/bin/python` and ended in
`PASS combined_covariance_r38_check`.  Independent quadrature at orders
`32,64,96,128` reproduced the displayed A8 and A9 signed integrals to at
least 14 decimal places.  This remains numerical evidence, not a certified
sign proof.

## Reconstruction

### H2 and H3: coordinate means

For every vertex coordinate `i`,

```math
X_i(D_{-i})=E_{\mu_t}[g\mid D_{-i}],
\qquad E_{\mu_t}X_i=E_{\mu_t}g=\psi'(t).
```

Thus, with `k_i=k_{i,D_-i}` context-measurable,

```math
\operatorname{Cov}_{M_i}(X_i,k_i)
=E_{\mu_t}[(g-\psi')k_i].
```

Summation gives `Cov_mu(g,sum_i k_i)`.  The means `E k_i` may depend on
`i`; they cause no missing between-coordinate term because every `X_i` has
the same mean.  This also proves the uniform-coordinate mixture formula H3.

### H4--H6: derivative and boundary sign

For one coordinate,

```math
\frac d{dt}E_{M_i(t)}k_i(t)
=tE_i(t)+\operatorname{Cov}_{M_i(t)}(X_i,k_i).
```

The first term uses `k'_e=tA''_e`; the second is exactly the derivative of
the moving context law.  Summing and applying H2 proves H4.  Since all
conditional binary costs vanish at `t=0`, integration proves H5 with the
displayed minus sign.  Therefore endpoint vertex cost plus
`[-integral Cov(g,L)]_+` bounds the vertex energy, proving H6.

H6 is a bookkeeping-equivalent sufficient condition once endpoint cost is
controlled, rather than independent progress on the energy.  The source
memo already states this, so there is no circular proof claim.

### Positive-part hierarchy

The inequalities in H7 have the correct direction:

```math
[-\int c]_+\le\int[-c]_+,
\qquad
[-\sum_i c_i]_+\le\sum_i[-c_i]_+,
```

followed by the already proved coordinatewise transient-flux majorant.
There is no interchange asserted as equality.  “Strict hierarchy” should be
read as a hierarchy of successively stronger hypotheses, not strict numeric
inequality; equality occurs in some displayed finite cases.

### H11 and H12: isotonic defect

For any nondecreasing `phi`, independent-copy symmetrization gives
`Cov(g,phi(g))>=0`.  Hence

```math
[-Cov(g,L)]_+
\le |Cov(g,L-\phi(g))|
\le\sqrt{Var(g)E(L-\phi(g))^2}.
```

Taking the infimum over the finite-dimensional closed isotonic cone proves
H11, and integration proves H12.  This is a valid sufficient condition; it
does not assume the desired covariance sign.

### H13--H15: matched normalization and fixed size

For an omitted vertex `i`,

```math
B_i(D_{-i})=U(D)r_i(D),
\qquad \widetilde B_i=f(D)r_i(D).
```

Thus `f=tilde B_i/r_i` on each `i`-edge and

```math
X_i=\log\widetilde B_i+E_{\mu_t(D_i\mid D_{-i})}[-\log r_i].
```

Its context average is `E_mu log f=E_mu g`, so no coordinate-dependent
centering is hidden in H14.  Pointwise posterior fixed size gives

```math
\sum_{i=1}^n\widetilde B_i
=f\sum_i r_i=(n-m)f,
```

which proves H15.  H15 includes the chart reference vertex whereas H2 uses
only the `n-1` coordinate vertices.  The memo states this distinction and
does not substitute H15 directly into `L_t`; doing so without separately
handling the reference term would be invalid.

The backward-erasure formula H17 was also checked independently at five
temperatures on A4, A6, A8, and A9.  The largest probability identity error
was `5.7e-14`, and `k_i(t)<=t q_i^leftarrow(t)` held to floating precision.

## Numerical scope

- The algebraic identities are finite-state exact identities; their checker
  evaluations are floating-point audits.
- The integrated values use 96-point Gauss--Legendre quadrature.  The tiny
  boundary residual is a consistency check, not a rigorous quadrature error
  certificate.
- The negative A8/A9 covariance and the rejection of universal exact
  monotonicity must remain labelled numerical unless an interval or exact
  sign certificate is added.  The source memo does label these values
  numerical.
- The checker deliberately omits the global-orientation edge.  This matches
  the theorem's vertex-only notation and is not an omission from a claimed
  full-coordinate calculation.
