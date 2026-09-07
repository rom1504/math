# Independent audit of the balanced-incidence ensemble comparison

2026-09-07. Status: PASS, using the signed Eulerian counting theorem in
`principle_director_signed_eulerian_entropy_2026_09_07.md`. This checks the
pointwise measure ratios, all-order realization, and negative-disorder
pressure normalization. It does not assert a pressure limit.

Let the connected even-degree graph have degrees `2k_v`, `e` edges, and
`D=product_v binom(2k_v,k_v)` balanced incidence arrays. The prescribed-sign
fiber count obeys

```math
\frac{D}{2^e\prod_v k_v}\le C_G(S)\le 2^e
```

on precisely the parity class `product_e S_e=(-1)^e`. The upper bound fixes
one incidence at every edge; the other is then determined. Uniform measure
on this parity class has atom mass `2^{-(e-1)}`. Consequently

```math
\frac{dP_G}{dU_G^*}(S)=\frac{2^{e-1}C_G(S)}D
\in\left[\frac1{2\prod_v k_v},\frac12\prod_v(2k_v+1)\right],
```

where the upper bound uses
`binom(2k,k)>=2^(2k)/(2k+1)`. Thus every nonnegative function, including a
single-atom indicator, has the claimed two expectation inequalities.

For odd `n`, use `K_n`; for even `n>=4`, remove a perfect matching and fill
that matching independently. The graph is connected and has positive even
degree. A fixed remaining edge exchanges the two graph-parity classes.
Changing it affects every spin energy by magnitude two, hence changes the
absolute cap by at most two. The support optimum therefore belongs to
`[M_n,M_n+2]` without a continuity assumption.

For `f(A)=exp[-lambda sqrt(n) Q(A)]`, the one-edge bijection compares the two
uniform parity-class expectations by factors between
`exp(-2 lambda sqrt(n))` and `exp(2 lambda sqrt(n))`. Averaging those two
expectations preserves these bounds. Combining this with the pointwise
density comparison gives the director's displayed estimate

```math
|\log Z_n^{\rm inc}(\lambda)-\log Z_n(\lambda)|
\le n\log(n-1)+\log2+2\lambda\sqrt n.
```

For odd `n`, the upper density logarithm is `n log n-log2`; it is indeed
covered by the displayed bound because
`n log(n/(n-1))<=2 log2` for `n>=3`. For even `n`, `2k_v+1=n-1` directly.

Finally the smallest cap term bounds the expectation from above, while one
minimizing signing, of uniform mass `2^{-binom(n,2)}`, bounds it from below:

```math
\frac{M_n}{n^{3/2}}
\le\frac{-\log Z_n(\lambda)}{\lambda n^2}
\le\frac{M_n}{n^{3/2}}+
\frac{\binom n2\log2}{\lambda n^2}.
```

Hence existence of the incidence pressure limit for every fixed positive
`lambda`, followed by `lambda -> infinity`, would imply the original cap
limit. The quantifiers and normalization are correct.

## What independent rows do and do not remove

The reference row variables are independent, but each carries order `n`
bits. The negative-disorder tilt couples them through the complete Boolean
maximum. Product-row trial laws therefore cannot be substituted for arbitrary
tilted laws without a separate argument. In particular, a product-row law
whose edge products are almost surely one prescribed signing must make each
individual incidence deterministic: independent endpoint signs with a
deterministic product are each deterministic. Such a law loses the large
conditional fiber entropy that the signed count preserves.

Removing a vertex also changes the balance condition in all the surviving
rows. An elementary coupling changing order `n` incidences only gives an
order-`n` deterministic cap error and thus an order-`n^(3/2)` logarithmic
weight error at the displayed critical tilt. This observation is a scale
warning, not a no-go theorem for a sharper cavity comparison.

No total-variation approximation, finite-state compression, pressure
subadditivity, or equality of positive-temperature product-row variational
values is inferred.
