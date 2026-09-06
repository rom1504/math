# Independent audit of the cubic Boolean chain remainder

2026-09-06. PASS. Independently reconstructed the complete argument in
`transfer_director_boolean_cubic_chain_nuclear_remainder_2026_09_06.md`.
The conclusion is precisely the chain remainder, with the random Stein
brackets retained; it is not a Gaussian or covariance-substitution result.

## 1. Exact cancellation and constants

Conditioning off coordinate a makes both `Delta_a Z` and
`d=Delta_a V` deterministic. With `h(t)=f(V0+td)`, the conditional
chain error divided by `Delta_a Z` is

```
[h(1)-h(-1)-h'(1)-h'(-1)]/2.
```

Twice integrating by parts gives

```
integral_-1^1 (1-t^2) h'''(t) dt
 =2[h'(1)+h'(-1)-h(1)+h(-1)].
```

Thus the coefficient is exactly -1/4, and its absolute bound is
`||D3 f||op,infty ||d||^3/3`. The cubic `h(t)=t^3` gives error -2 and
attains the bound `6/3=2`. No unaveraged second-order term survives.

## 2. Elementary fourth moment and vector transport

For `P=P0+s P1`, put `a=||P0||4^2`, `b=||P1||4^2`.
The exact fourth-power expansion and Cauchy--Schwarz give

```
||P||4^4 <= a^2+6ab+b^2 <= (a+3b)^2.
```

Induction yields `||P||4^2<=sum_S 3^|S| phat(S)^2`.
For a degree-d vector query, Minkowski applied to the sum of squared
derivatives consequently gives

```
(E ||Delta_a V||^4)^(1/4)
 <= 3^((d-1)/2) sqrt(I_a(V)).
```

This uses the aggregate coordinate influence, so it remains valid for
varying vector dimension with the stated aggregate hypotheses.

## 3. Exact brackets and nuclear scale

For a pure positive-degree K Fourier chaos Z, Fourier orthogonality gives
the exact identity

```
E[Z U]=(1/K) sum_a E[Delta_a Z Delta_a U]
```

for every cube function U. Applying the conditional identity above leaves
the ACTUAL bracket `sum_a Delta_a Z Delta_a V_l` in the first-order term.
Hölder and the two fourth-moment estimates give the error entry bound

```
M/(3K) * 3^((K+3d-4)/2) * sqrt(c/n)
    * sum_a I_a(V_i)^(3/2)
 <= [M I0 sqrt(c)/(3K)] 3^((K+3d-4)/2) sqrt(kappa_n/n).
```

The exponent is `(K-1)/2+3(d-1)/2`, exactly as stated.
For the n-by-n error matrix,
`||R||_*<=sqrt(n)||R||F<=n^(3/2)max|R_ij|`, proving the asserted
normalized bound `C sqrt(kappa_n)`.

All constants and dimension factors pass. The low-influence hypothesis is
essential to obtaining o(n) nuclear scale from this estimate; diffuse Z
alone does not establish it for V. In particular an order-one diagonal
coordinate in a coherent QS query is not covered by a vanishing-kappa
claim. No bracket fluctuation, normal approximation, polynomial-density,
or high-influence closure statement was used in this audit.

## 4. Scope of the illustrative cubic query

The multilinear representation of `(n^-1/2 sum_a S_a)^3` has degree at
most three. Its a-th derivative equals
`3(n^-1/2 sum_(b!=a) S_b)^2/sqrt(n)+n^-3/2`, hence has squared L2 norm
O(1/n) and bounded total influence. Thus it is a valid low-influence
query example and converges to a cubic of a normal random variable.
The source invokes a separate banked moment counterexample for its
moment-indeterminacy comment; that comment is not a dependency of the
chain theorem proved here.
