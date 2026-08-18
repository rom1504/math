# Adversarial audit: global actual-child row-product certificates

**Objects audited:**

- [`drafts/actual_child_row_product_global_certificate.md`](../drafts/actual_child_row_product_global_certificate.md)
- [`experiments/actual_child_row_product_certificate.py`](../experiments/actual_child_row_product_certificate.py)
- [`../../computations/results/actual_child_row_product_certificate.json`](../../computations/results/actual_child_row_product_certificate.json)

**Verdict:** **PASS.**  GC.1 has the correct constant four and the a
posteriori interval has the correct sign and inverse Hessian.  KL data
processing is in the required direction.  The exact N=8 twin reduction and
the N=9 uniform-global argument are valid.  The implementation now encloses
all transcendental quantities used in the certificate with `mpmath.iv`,
including binary entropy at exact dyadic subdivision endpoints, `log(2)`,
the child-minimum comparisons, negative-moment soft pressure, and the
same-temperature target.  A clean hardened rerun reproduced every reported
field except wall-clock values and certified the N=8 target excess lower
bound `0.19973600675473155`.

## 1. GC.1 constants and a posteriori bound

With `delta_i=p_i-r_i` and `d_i=TV(p_i,r_i)`, the convention

```math
\delta_i=d_i(\alpha_i-\beta_i)
```

is correct.  For the multilinear energy interpolation,

```math
g''(s)=2\sum_{i<j}E_{\delta_i\delta_jr_{-ij}(s)}f.
```

Taylor's integral contributes
`-sum_(i<j) C_ij d_i d_j`.  The exact entropy identity around a fixed point
leaves `(1/lambda)sum_i D(p_i||r_i)`, and Pinsker gives
`(2/lambda)sum_i d_i^2`.  Together these equal

```math
{1\over2}d^T[(4/\lambda)I-C]d,
```

so neither factor two is missing.  Since `C` is symmetric,
`lambda lambda_max(C)<4` is exactly positive definiteness of this Hessian
minorant; negative eigenvalues of `C` cause no problem.

For an approximate fixed point, a signed mass-zero row perturbation pairs
with the residual by at most `eta_i d_i`.  Minimizing the resulting quadratic
over all real vectors (a larger set than `d_i>=0`) gives
`-eta^T H^{-1}eta/2`.  Thus (GC.7) is a valid, possibly conservative,
a posteriori lower bound.

## 2. Data-processing direction

A rowwise map sends every original row product to a coarse row product, and

```math
D(p\Vert q)\ge D(\phi_\#p\Vert\phi_\#q).
```

The set of images is only a subset of all coarse products.  Minimizing first
over original products and then enlarging the coarse feasible set therefore
preserves the direction in (GC.12).  The theorem gives a lower certificate
for the reverse projection, not for forward total correlation.

## 3. N=8 exact reduction

The integer signature-count equalities in (GC.18) imply the zero Walsh
coefficients and the repeated `J,k` coefficients exactly, before any
transcendental evaluation.  In the four-variable product objective, the
first two stationarity equations differ by

```math
\operatorname{atanh}s_1-\operatorname{atanh}s_2
 +(J-Ks_3s_4)(s_1-s_2)=0.
```

Because `J-K>0`, both summands have the sign of `s_1-s_2`, proving
`s_1=s_2`; similarly `s_3=s_4`.  Full support excludes boundary minimizers.
After writing the twin means as `x,y`, changing one sign when `xy<0`
strictly improves the `-4kxy` term and leaves every other term fixed.
Global complement symmetry then reduces to `u,v in [0,1]`.  Equation
(GC.23) and all its signs are correct.

On a nonnegative box, `p log p+(1-p)log(1-p)` is increasing with the mean,
so evaluating it at the point closest to zero is a valid lower bound.  The
interval choices for the negative `J,k` terms and positive `K` term are also
in the correct directions.  The hardened routine reconstructs each binary64
subdivision endpoint by `as_integer_ratio` and evaluates its entropy using
outward `mpmath.iv` arithmetic.  The remaining elementary binary64
accumulation is explicitly bounded by fewer than 64 operations of magnitude
below 16, while the `1e-12` downward safety margin exceeds the resulting
worst-case roundoff bound by more than a factor four.  Thus the reported
branch-and-bound lower bound `1.075` is a valid outward certificate.

## 4. N=9 uniform-global claim

Complement symmetry eliminates all odd-degree log-density Walsh
coefficients.  Averaging the effective potential over uniform values of the
other rows therefore leaves no singleton field, so the uniform product is a
coordinate Gibbs fixed point.  The coarse rectangle matrix is nonnegative
and symmetric; its certified maximum row sum bounds its spectral radius.
The reported `.571` and `1.266` margins are far below four, so GC.1 proves
the uniform product is the unique global minimizer.  The value is then
`-a_0-4log 2`, as used in the result file.  The hardened implementation
evaluates `log 2` by outward `mpmath.iv` arithmetic, so this endpoint is also
formally enclosed.

## 5. Reproduction and scope

Running

```text
.venv/bin/python extremal_information/experiments/
actual_child_row_product_certificate.py
```

to a file under `/home/math/quadra/tmp` reproduced all probability,
coefficient, rectangle, branch-and-bound, projection, child-minimality,
soft-pressure, target, and target-excess fields exactly; only the recorded
wall times differed.  The full-alphabet rectangle spectral radii are
correctly labelled floating falsification diagnostics, not proof
certificates.  The finite lower bounds are genuine actual-child dependence
witnesses; they do not assert extensive asymptotic dependence.
