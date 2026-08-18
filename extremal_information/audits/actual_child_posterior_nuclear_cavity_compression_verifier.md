# Independent audit of posterior nuclear cavity compression

**Disposition: PASS with one corrected bit count and a scope qualification.**
The signed posterior-replica cost in the source is now `R(N-1)`, not
`R(N-2)`.  The theorem gives querywise geometry; it does not synchronize a
common frame or convert one to product gain.

## Constant checks

Every signed rank-one word has nuclear norm `sqrt(mn)`, so convexity gives
`||M||_*<=sqrt d` exactly.  Bit insertion gives

```math
M_e={r_e+rho B_e\over1+rho B_er_e},
```

and `|M_e-r_e|<=2rho`; summation gives `4rho^2d`.  If the singular values
are decreasing, then

```math
sigma_(R+1)<=||M||_* /(R+1),
\qquad
\sum_(j>R)sigma_j^2<=d/(R+1).
```

The triangle inequality and `(a+b)^2<=2a^2+2b^2` therefore give precisely
the constants `2` and `8` in PN.11--PN.14.  At physical scale,
`t^2d=Theta(N)` and `t^2rho^2d=O(1)`.

For posterior replicas, conditional variance is exactly

```math
E[||\bar Q_R-M||_F^2\mid B]
={d-||M||_F^2\over R}.
```

Finally,

```math
||M||_F^2<=sigma_1(M)||M||_*<=sigma_1(M)sqrt d
```

verifies the overlap-to-spike conversion.  Positive expected normalized
overlap yields a spike on fixed positive path-mixture mass because
`sigma_1<=sqrt d`.

## Metric entropy

Rank truncation followed by Stiefel and singular-value nets at precision
`Theta(R^(-1/2))` gives logarithmic covering size `O(RN log(CR))` at
Frobenius scale `Theta(sqrt(d/R))`.  Conversely, equal-singular-value
constant-angle Grassmann packings give `Omega(RN)` at the same scale for
comparable sides and `R<=c min(m,n)`.  This is an ambient nuclear-ball
statement, not a claim that actual posterior means fill the packing.

## Scope of combination with the actual overlap floor

Apply Theorem 37.56 to the path-mixture law

```math
nu_N={1\over lambda}\int_{-lambda}^0q_s\,ds.
```

Then PN.4 proves a macroscopic posterior singular value under `nu_N`, and
equivalently at some tilt.  It does **not** prove the endpoint statement at
`q_lambda`, one common singular direction, or synchronization across bridge
words or tilts.  Frame response entropy is therefore useful only together
with an integrable child-generated rule and a directional conversion to
product gain or target reach.
