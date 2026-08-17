# Adversarial audit: dense-bridge packing for sign quadratics

**Verdict: PROMOTE, with two scope/notation edits.**  The deterministic
pole-locking argument, the exponential random construction, and the
projective response separation are correct.  The theorem gives an intrinsic
`Omega(n)`-bit lower bound at error `Theta(n^(3/2))` for the class of all
complete sign-quadratic children.  It does **not** give such a lower bound for
bounded-cap or near-minimizing children: the planted quadratics used here have
unperturbed maximum `(n^2-n)/2`.

## 1. Deterministic checks

For

```math
H_z(x)=((x^Tz)^2-n)/2
```

and the nearer pole `sz`, put `r=d_H(x,sz)<=n/2`.  The exact loss is
`2r(n-r)>=nr`; changing `r` coordinates can improve the field by at most
`2r||h||_infinity<nr`.  Hence every non-pole loses strictly, including the
tie case `r=n/2`, and the two poles have scores
`(n^2-n)/2 +- h^Tz`.  Lemma QD.1 is therefore correct with no parity caveat.

The response is consequently

```math
(P_BH_{z_d})(y_c)=(n^2-n)/2+|z_d^TBy_c|.
```

If the diagonal term is at least `d_0 n^(3/2)` and every off-diagonal
absolute term is at most `d_1 n^(3/2)`, then for the difference of children
`c,d` the query `y_c` gives at least `(d_0-d_1)n^(3/2)` and `y_d` gives at
most its negative.  Thus the sup-norm distance is at least this amount and
the projective distance, one half of the oscillation, is also at least this
amount.  With `d_1=d_0/2`, the common choice `g=d_0/2` in QD.2 is correct.

The off-diagonal condition also proves `z_c != +-z_d`, since otherwise the
off-diagonal absolute term at `y_c` would equal its diagonal term.  Hence the
quadratics, which identify `z` and `-z`, really are pairwise distinct.

The contextual identities are exact.  Max-plus transfer is sup-norm
nonexpansive, while a sufficiently negative future outside one selected
query pins that query for both response functions.  Modulo constants,
`inf_a ||f-g-a||_infinity=osc(f-g)/2`.  A projective packing at distance
`g n^(3/2)` therefore requires distinct decoded states whenever
`2 epsilon<g`.

## 2. Probabilistic audit

Take independent Rademacher `B` and independent uniform `y_c`.

1. The standard rectangular-sign-matrix estimate gives
   `P(||B||_op>C_0 sqrt(n))<=2 exp(-c_op n)`.
2. For fixed `y`, the `n` coordinates of `By` are independent length-`n`
   Rademacher sums.  Their absolute values have mean at least `sqrt(n/2)`
   and subgaussian scale `O(sqrt(n))`; concentration of their sum at a fixed
   fraction below its mean is `exp(-Omega(n))`.  This justifies the diagonal
   estimate uniformly after a union bound over `N` queries.
3. Each coordinate of `By_c` is a length-`n` Rademacher sum, so
   `P(|(By_c)_i|>=n/2)<=2e^{-n/8}`.  The `nN` union bound gives the strict
   hypothesis needed by QD.1.
4. Conditional on `B,y_d`, `z_d` is fixed and `y_c` remains independent for
   `c!=d`.  On the operator event,
   `||B^Tz_d||_2<=||B||_op sqrt(n)<=C_0 n`.  The two-sided Rademacher tail at
   `d_1 n^(3/2)` is therefore at most
   `2 exp(-d_1^2 n/(2C_0^2))`.  Conditioning creates no dependence problem;
   summing these conditional bounds over all ordered pairs is a valid union
   bound.

With `N=floor(exp(2 gamma n))`, the four failure exponents are respectively
`c_op`, `c_diag-2gamma`, `1/8-2gamma`, and
`c_cross-4gamma`.  The stated choice

```math
gamma < (1/8) min(c_op,c_diag,1/8,c_cross)
```

makes all positive.  Thus the simultaneous event has positive probability
for all sufficiently large `n`.  No independence between the diagonal
events is being assumed.

## 3. Coefficient and normalization audit

Each planted child has all off-diagonal coefficients `z_i z_j in {+-1}` and
no diagonal coefficient, exactly the declared class.  The bridge term
`x^TBy` uses `n^2` signs, and both the diagonal/off-diagonal response signal
and its separation are on the claimed `n^(3/2)` scale.  The common child
maximum is order `n^2`, but it is exactly the same additive baseline for
every child and cancels projectively.  This is mathematically valid, while
also explaining the theorem's essential scope limitation.

Coefficient rounding changes every landscape value by at most
`p Delta/2`, hence every bridge response and every later scalar optimum by
the same amount.  QD.3 is correct.  For precision, the uniform bit bound is
best displayed as

```math
O(n^2 log(2+sqrt(n)/epsilon)),
```

rather than `O(n^2 log(n/epsilon))`; the latter is harmless for the fixed
small `epsilon` regime of the theorem but is not uniform in arbitrary
`epsilon`.

## 4. Reproducible verification

Running

```bash
python3 extremal_information/experiments/verify_quadratic_dense_bridge_response.py
```

completed successfully.  It performed 1,400 exhaustive small-order
pole-lock checks and 3,240 exact rational quantization/continuation checks.
Its seeded `n=32`, 12-query certificate has maximum field coordinate 14,
minimum directed gap 56, and minimum projective gap 68, satisfying the
strict pole-lock hypothesis and the claimed response separation.

## 5. Exact implication and non-implications

The proved statement rules out any uniform summary with `o(n)` bits that
answers all future continuations, at `epsilon n^(3/2)` accuracy for fixed
small `epsilon`, on **all** complete sign-quadratic children across the
constructed dense sign bridge.  This is stronger than an algorithm-specific
lower bound and weaker than the arbitrary-landscape doubly-exponential state
obstruction.

It does not rule out:

- an `O(n)`-bit state (the constructed pole subclass itself has only
  `n-1` independent bits);
- a subquadratic state for the full sign-quadratic class;
- compression specialized to bounded-cap or near-minimizing signings;
- restricted future/query families; or
- a compositional proof that does not uniformly preserve the full response
  function.

The main draft should state the bounded-cap limitation next to its headline
conclusion and replace the coarse logarithm in QD.3 if a uniform-in-`epsilon`
claim is intended.  Subject to those edits, the theorem is ready for the
rigorous results file.
