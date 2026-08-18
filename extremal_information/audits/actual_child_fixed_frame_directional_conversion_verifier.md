# Independent audit of integrable-frame directional conversion

**Disposition: PASS.**  The signs, constants, Walsh inversion, and physical
asymptotics in FF.1--FF.3 are correct.

Central symmetry gives `Z=log E cosh>=0`.  Direct substitution gives

```math
G=lambda A+log z,
\qquad
J=\sum_iD_i+lambda E_rZ+log z.
```

Using `log z<=0`, Jensen's `log z>=-lambda A`, and the fair-product
competitor proves FF.6 and its factor-two alternative with the displayed
directions.  The identity `G=lambda(E_UL-V_lambda)` has no missing constant.

For the roof, each `t<B,C>` is centered subgaussian with proxy
`t^2||C||_F^2`, giving `t sqrt(2d log K)`.  At physical scale this is
`O(sqrt(N log K))=o(N)` when `log K=o(N)`.

In the fixed-frame theorem, all vector-valued Walsh coefficients are
`rho^|S| E[Pi_perp Q chi_S(Q)]`.  Since `rho>0`, inversion forces every
support atom into the fixed span, hence into `{C,-C}`.  The one-row KL limit,
`O(sqrt N)` fair remainder, local-CLT `sech^a` law, and overlap integrals
`1/(1+a)` and `log(1+lambda)/lambda` are all correct.  Deleted cavities have
the same limit by `|M_e-r_e|<=2rho`.

The scope is essential: a subexponential number of real frames need not have
a low-bit coordinate description, adaptive posterior low rank is not an
integrable fair roof, and the one-row cost plus target calibration remain
separate obligations.
