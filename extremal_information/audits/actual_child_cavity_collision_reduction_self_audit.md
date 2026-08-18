# Self-audit: collision--cavity reduction

**Object audited:**
[`../drafts/actual_child_cavity_collision_reduction.md`](../drafts/actual_child_cavity_collision_reduction.md).

**Verdict:** **PASS as an exact reduction, a sufficient theorem, and two
generic scalable falsifiers; no actual-minimizer closure is claimed.**  An
independent derivation checked CC.8, CC.11, CC.15, and CC.21--CC.23.

## Normalizations

The bit-channel likelihood relative to a fair row is
`exp(u<b,z>)/(cosh u)^n`, not `2^(-n)` times that expression.  Averaging over
the fair conditional sign `sX_i` changes the exponential to `cosh` and gives
exactly CR.0's row likelihood.  Hence all `(cosh u)^(-n)` factors cancel
between the numerator and denominator in (CC.8).  The sector weights already
sum to one, so no orientation factor is missing.

Substitution of `p=G prod_i z_i` gives

```math
 {dq\over dr}={G^{-\lambda}\over E_rG^{-\lambda}},
```

and direct evaluation of `D(r||q)` confirms both signs in (CC.6).

## Collision--cavity split

For fixed `(s,y)`,

```math
 k_u(b|s,x_i,y)
 ={\cosh(u<b,y>)\over(\cosh u)^n}
  [1+s x_i\tanh(u<b,y>)].
```

This verifies every factor in (CC.9)--(CC.11).  The expectation over the
single shared `(s,y)` is essential.  Replacing it by independent copies
would erase `C_D`; bounding `C_D` and the left response separately would
lose their cancellation and is intentionally avoided.

## CC.2 constant

For a function of independent coordinates whose `i`th coordinate range is
at most `c_i`, the Doob-martingale form of Hoeffding's lemma gives

```math
 log E exp(theta(f-Ef)) <= theta^2 sum_i c_i^2/8.
```

Bayes' rule proves that the row range of `log G` is exactly the projective
oscillation in (CC.14), so the factor `lambda^2/8` is correct.  The theorem
uses a supremum over all other rows; it is therefore rigorous but may be too
strong asymptotically.

For CC.3, Efron--Stein contributes
`sum_i E Var(h|B_(-i))` with no extra factor: the usual one-half cancels the
factor two from two iid resamplings of row `i`.  Popoviciu's bound then gives
`Var(h|B_(-i))<=delta_i(B_(-i))^2/4`.  This controls only the product
endpoint and is not silently extended to the dependent hybrid laws.

## Scalable example

For `Q=+-11^T`, summing the two planted words gives the two cosh formulas in
(CC.19).  Under the canonical row law, multiplication by
`exp(lambda sum_i g(V_i))` cancels `prod_i cosh(V_i)^(-lambda)` exactly.
This proves (CC.21).  Jensen gives

```math
 log E_U exp[-lambda g(sum_iV_i)]
 >=-lambda E_U g(sum_iV_i),
```

and `E_U g(sum_iV_i)<=sqrt(E_U(sum_iV_i)^2)=c sqrt(r)`, proving (CC.22).
The limiting tilted Gaussian differs from the Gaussian whenever `c,lambda`
are positive, so its KL divergence is strictly positive.  The example is a
generic central rank-one prior, not an actual optimized-child law.  For
`r=1`, `g(V)` is constant and `kappa_r=0`; the asymptotic statement and the
claim `kappa_r>0` correctly use `r>=2`.

In the block-parity example, expanding the independent block law leaves
only the constant and full `ell`-spin character, proving (CC.27).  All
moments of order below `ell` vanish.  The canonical row law is a product,
so the centered log-MGF adds exactly across blocks.  The tilted CLT and
boundedness from `|delta|<1` justify (CC.29); strictness follows because the
limiting block product is nonconstant.

## Information content

The exact function `G` still retains every even left-child correlation via
(CC.12), so calling CC.8 alone a difficulty reduction would be too strong.
The scalar `Delta_u` is genuinely less information than the full conditional
tables and CC.2 is a strict sufficient reduction.  Naively verifying its
supremum can nevertheless require exponential work, so it is not yet an
algorithmic compression.  Proving it for actual minimizers remains open.
