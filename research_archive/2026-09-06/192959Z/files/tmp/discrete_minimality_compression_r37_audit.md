# Independent audit: discrete minimality compression, Wave 37

## Verdict

**PASS after two exact corrections.**

1. The incidence matrices are enumerated exactly, but the claimed fractional
   cover values `2,10,55/4` are obtained only from floating-point
   `scipy.optimize.linprog`.  The current checker does not contain rational
   primal and dual certificates.  Either label these three LP values
   **Numerical**, or add rational feasible primal/dual vectors with equal
   objectives and verify every constraint in `Fraction` arithmetic before
   calling them exact.
2. In finite diagnostic point 2, replace “the trace/Markov fact that many
   cuts have row at most the mean” by “the natural mean row cap” (or “the
   fact that at least one cut has row at most the mean”).  Mean alone gives
   at least one such cut; Markov gives at least half only at cap
   `2n(n-1)`, not at `n(n-1)`.  The checker deliberately uses the mean cap,
   so changing the prose is preferable to changing the computation.

All other signs, factors, regimes, constants, implications, and scope
statements check.

## Algebra audit

Write `E=2 sum_e s_e` and `c_S=2 sum_(e in E(S))s_e`.  Flipping the complement
of `E(S)` changes the oriented energy by

```math
E-4\sum_{e\notin E(S)}s_e
=E-2(E-c_S)=2c_S-E,
```

so (D37.1) has the correct sign and factor.  Since the perturbed complete
signing has norm at least `q_n`, some lifted state has positive oriented
energy at least `q_n`.  Substituting `E=q_n-Delta` gives

```math
c_S\ge q_n-\Delta/2,
```

exactly (D37.2).

The loss identity is

```math
h_d(S)=Q(A[S])-c_S-p_2\Delta-B_{n,m}.
```

Using `Q(A[S])<=q_n` and (D37.2) gives

```math
h_d(S)\le Q(A[S])-q_n+(1/2-p_2)\Delta-B_{n,m}.
```

Thus (D37.3) is valid precisely when `p_2>=1/2`.  Also
`B_(n,m)>0` for every proper selector because
`p_2<p^2<p^(3/2)`.  The asymptotic threshold is a ratio strictly above
`1/sqrt(2)`, with the stated finite correction.

## Fractional cover and sufficiency

(D37.4) and (D37.5) are the standard finite covering LP and its exact dual.
If `w` is feasible, averaging its row constraints over `U_m` gives

```math
1\le\sum_dw_dU_m(I_d).
```

Therefore one row-good column has uniform mass at least `1/tau_flip`.
Equation (D37.3) places all of that incidence inside `{h_d<=0}`.  With the
scales in (D37.6), (10.1023) yields (D37.8), so the claimed sufficient
package is correct.

Across a fixed normalized gap, (10.1022) bounds each row-good incidence by
`exp{-Omega(n^(3/4))}`: its three exponents are respectively at least
`n^(3/4+c),n`, and `n^(3/4)`.  Averaging any fractional cover then forces
`tau_flip>=exp{Omega(n^(3/4))}`, proving (D37.9), including the case
`tau_flip=+infinity`.

## Full edge-cube volume

For `X~Bin(N,1/2)` and
`R=N/2-q_n/4`, Hoeffding gives

```math
P(X\le R)\le
\exp\{-2(q_n/4)^2/N\}
=\exp\{-q_n^2/(8N)\}.
```

Averaging a fractional radius-`R` cover therefore proves

```math
\sum_dw_d\ge
2^N/\sum_{j\le R}\binom Nj
\ge\exp\{q_n^2/(8N)\}.
```

The constant `1/8` is correct.  Since `q_n=Omega(n^(3/2))` and
`N=Theta(n^2)`, this is `exp{Omega(n)}`.  Removing centers by a row cap can
only make a retained full-cube cover harder, so the scope statement is also
correct.

## Novelty relative to (10.733) and (10.744)

- (D37.1)--(D37.2) are a structured-selector specialization of the already
  proved full edge-cube cover (10.733), as the memo correctly acknowledges.
- The favorable sign consequence (D37.3) in the `p_2>=1/2` window and the
  row-restricted fractional LP `tau_flip` are the new packaging.
- The sphere-volume lower bound for retaining the full cube is a new direct
  consequence of (10.733), not a new minimality identity.
- (D37.9) is the row-sensitive fractional-cover analogue of the earlier
  codebook converse (10.744).  Its `n^(3/4)` exponent comes from the new
  row-sensitive tail (10.1022); it is not merely the old general
  `Omega(sqrt n)` statement restated.

