# Audit of the actual-child row-product global certificate

Status: **independent adversarial audit; PASS for the analytic theorems,
coarse lower certificates, and target-excess theorem after certificate
hardening**.  The imported full-product upper endpoints and complete-alphabet
spectral radii remain numerical, as explicitly classified below.

Audited source:
`drafts/actual_child_row_product_global_certificate.md` and
`experiments/actual_child_row_product_certificate.py`.

## 1. Rectangle-Hessian theorem

Put `delta_i=p_i-r_i` and `d_i=TV(p_i,r_i)`.  Along the simultaneous affine
interpolation of all marginals,

```math
g''(s)=2\sum_{i<j}
 E_{\delta_i\otimes\delta_j\otimes r_{-ij}(s)}f.
```

Writing `delta_i=d_i(alpha_i-beta_i)` turns each mixed expectation into an
average of row rectangles, so its absolute value is at most
`C_ij d_i d_j`.  Taylor's integral remainder is therefore at least
`-sum_(i<j) C_ij d_i d_j`.  The exact entropy expansion about `r_i` is

```math
D(p_i\Vert u_i)-D(r_i\Vert u_i)
=D(p_i\Vert r_i)+\langle p_i-r_i,\log(r_i/u_i)\rangle.
```

At a Gibbs fixed point its linear part cancels `g'(0)`.  Pinsker gives
`D(p_i||r_i)>=2d_i^2`, hence

```math
F_\lambda(p)-F_\lambda(r)
\ge {1\over2}d^T[(4/\lambda)I-C]d.
```

Thus the factor four in (GC.4)--(GC.5) is correct.  Since the least
eigenvalue of `(4/lambda)I-C` is
`4/lambda-lambda_max(C)`, strict positivity proves global uniqueness, not
merely local stability.  For an approximate fixed point, a zero-mass signed
measure obeys `|<delta_i,h_i>|<=d_i osc(h_i)`; completing the square gives
(GC.7) with the stated factor one half.  No simplex-boundary case is lost,
because the inequalities hold for arbitrary product laws while `r` has full
support.

## 2. Data processing and exact coarse reduction

For every row-product `p`, its deterministic rowwise image `P` is product and

```math
D(p\Vert q)\ge D(P\Vert Q).
```

Taking the infimum first over realizable images and then enlarging to all
product laws on the feature alphabet gives exactly the direction (GC.12).
It is a lower certificate for the reverse projection; it is not forward
total correlation.

For the order-eight feature map, equality of the four atom classes is an
integer equality of pressure-signature count columns.  Their four values
span exactly

```math
1,\quad y_1y_2+y_3y_4,\quad
(y_1+y_2)(y_3+y_4),\quad y_1y_2y_3y_4,
```

so (GC.20) is exact rather than a numerically truncated Walsh expansion.
Every minimizer of (GC.15) is interior.  Holding `(s_3,s_4)` fixed, subtraction
of the first two stationarity equations gives (GC.22); since
`J-Ks_3s_4>=J-K>0`, both summands have the sign of `s_1-s_2`, forcing equality.
The second twin pair is identical.  Finally `k>0` makes equal signs weakly
better for fixed absolute values, and global sign symmetry permits
`u,v in [0,1]`.  Thus the two-dimensional reduction has no omitted branch.

Feature-mask selection was post hoc but creates no statistical-validity
issue: after a deterministic map is frozen, KL data processing is a pointwise
identity for that same finite law.  No sample splitting or probabilistic
confidence statement is used.

## 3. Computer-assisted lower bound

The hardened script evaluates each entropy endpoint by `mpmath.iv` from the
endpoint's exact dyadic rational (`as_integer_ratio`), so the platform
`libm` logarithm is absent from the proof.  The remaining per-box arithmetic
has fewer than 64 binary64 basic operations and, under the asserted runtime
coefficient check and boxes in `[0,1]^2`, every intermediate has absolute
value below 16.  With IEEE-754 round-to-nearest, the conservative accumulated
absolute error bound

```math
64\cdot16\cdot2^{-52}<2.28\times10^{-13}
```

is below one quarter of the `10^{-12}` downward subtraction.  This also
covers an inward one-ulp rounding in the elementary interval products.
The final `nextafter` is downward.  The smallest pruned lower bound in the
fresh run is `1.075000000841172`, over `8.4e-10` above the exact decimal
target `1.075`; hence the fact that the Python input float is the neighboring
binary number cannot affect the claim.  All `752768` boxes are exhausted.
The feasible coarse product is enclosed in

```text
[1.075619910790213, 1.0756199107922173],
```

so (GC.24) passes independently.

## 4. Actual minimizers, order-nine uniqueness, and target excess

The fresh script first selects the child histogram and then independently
encloses every distinct histogram pressure.  In all three cases the selected
upper endpoint lies strictly below every competitor lower endpoint (the
smallest recorded contracted-temperature gap is about `.324`).  Thus the
finite applications really use optimized children; the initial
high-precision selector is not the final certificate.

For the two order-nine coarse laws, the Walsh interval bounds give rectangle
row-sum upper bounds `.570579` and `1.265695`.  Global complement symmetry
kills every singleton coefficient of `-log Q`, so the uniform product is a
coordinate Gibbs fixed point.  Since these row sums upper-bound
`lambda_max(C)` and are below four, GC.1 proves that this fixed point is the
unique global coarse minimizer.  Evaluating it gives the intervals reported
in (GC.28).

At order eight the negative-moment soft value and same-temperature target are
recomputed independently as outward intervals.  Combining the soft lower
endpoint, the exact decimal `lambda_*`, and the coarse reverse-KL lower bound
gives

```text
V_row - T >= 0.19973600675473155.
```

Therefore the positive target excess does not assume that the numerically
located threshold is an exact root.  The rounded theorem claim `>=0.19973`
passes.

## 5. Evidentiary boundary and verdict

The complete-row rectangle spectral radii use floating pressure cubes and
`numpy.linalg`; they are correctly labeled numerical falsifiers and are not
used in a positive theorem.  Their margins over the threshold are very large.
The upper endpoints for the full row-product projection are imported
feasible coordinate-Gibbs evaluations padded by `10^{-6}`.  They are
reproducible numerical feasible-law bounds, not certificates of global
optimality; none of the theorem's nonzero lower bounds or the target-excess
conclusion depends on those upper endpoints.

Subject to that explicit numerical status, the certificate passes.  In
particular, there is no post-selection fallacy, no reversed data-processing
inequality, no missing symmetry branch, and no unsupported use of the
binary64 threshold root in the rigorous lower conclusion.
