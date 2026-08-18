# Independent audit of the global latent coreset theorem

**Disposition: PASS after minor clarifications.**  Every normalization,
constant, and exponent in GC.0--GC.7 checks.  The source now records the
reverse full/deleted comparison and the simultaneous probabilistic-method
selection.

## Algebra and sampling checks

The identity

```math
(1+rho z)^2=(1+rho^2)(1+tanh(2t)z)
```

gives the powers `d` and `d-1` in GC.11a exactly.  For a deleted posterior,
the good denominator event contributes at most `16K_e/R`; Chebyshev and the
bounded response range contribute at most `16(K_e-1)/R`.  Hence the constant
`32` in GC.16 is valid, and aggregate expectation selects one global sample
without a union bound.

The full posterior differs from the deleted posterior by a density factor in
`[e^(-2t),e^(2t)]`; squaring gives the two-sided `e^(4t)` comparison.  The
half-flip derivative of the empirical likelihood is exactly
`atanh(rho r_e^(R))`, and its Lipschitz constant on `[-1,1]` is
`rho/(1-rho^2)`.  Thus the empirical channel is genuinely curl-free.

## Logarithmic lower tail and exponents

For a positive mean-one variable with second moment `K`, Chebyshev gives the
`4(K-1)/R` bad-denominator probability and Paley--Zygmund gives
`Pr{W>=1/2}>=1/(4K)`.  Together with the deterministic channel floor this
proves GC.26 and the averaged truncation in GC.29.

With `H=K N^(1/2+zeta)` and `R=16H(log N)^2`,

```math
t^2d=Theta(N),\qquad td=Theta(N^(3/2)),
```

so the cavity/gradient error is
`O(N^(1/2-zeta)/(log N)^2)` and the lower-tail pressure term is
`O(N^(1-zeta))+exp{-Omega((log N)^2)}`.  All other pressure terms are `o(N)`
when `log K=o(N)`.  A normalized-sum probabilistic argument chooses one
coreset satisfying both conclusions at a harmless factor two.

## Archive and scope

The condition `log overline K_del=o(N)` is a scalar two-temperature
posterior-collision condition and retains strictly less output information
than the complete cavity table.  It is not yet proved easier for actual
minimizers.  Existing conditional row-Renyi theorems concern the output law
in the opposite channel direction; fair-`L^2`, min-entropy, retuning, and
replica-Gram results neither imply nor falsify it.  The theorem supplies one
global integrable coreset, but not an `o(N)`-bit mergeable state, directional
target relevance, or a Level-6 recurrence.
