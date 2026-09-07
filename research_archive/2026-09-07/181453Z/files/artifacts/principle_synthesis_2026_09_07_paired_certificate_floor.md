# A .497 floor for the paired anti-weave Bellman certificate, not its actual cap

2026-09-07. Status: proved, with a short outward-interval endpoint check.
The actual anti-family below `.499` remains valid. This is a floor only for
the unmodified paired-row Bellman/first-moment certificate used there;
it does not constrain the root's new physical-stability-filtered operation.

## 1. Two retained source channels already force the floor

The certificate of `principle_synthesis_2026_09_07_anti_invariant_weave.md`
has the form

```
C_pair(p,t) = [t+sup_theta R_(p,t)(theta)]/(2t sqrt(p)),
R(theta) = p/2 [log2+h(theta)]
  +[E_t(mu_(p theta))+E_t(mu_(p(1-theta)))]/2,
```

where `0<p<=1`, `t>0`, and `mu_s` puts mass `1-s` at zero and mass `s/2`
at each of `+-sqrt(2/p)`. At `theta=0`, reveal the nonzero ternary child
completely. Its conditional variance becomes zero and the information cost
is `h(p)+p log2`. This admissible channel gives

```
R(0) >= -h(p)/2.
```

At `theta=1/2`, both children have second moment one. The uninformative
channel gives

```
R(1/2) >= p log2+g_t(1).
```

Consequently

```
C_pair(p,t) >=
 max{t-h(p)/2, t+p log2+g_t(1)}/(2t sqrt(p)).       (1)
```

Here the exact Gaussian boundary is

```
g_t(1)=max_(0<lambda<=t) [c_t(lambda)-lambda],
c_t(lambda)=.25 log[lambda(2t-lambda)/t^2],
lambda_t=2t/[4t+1+sqrt(16t^2+1)].
```

We prove that the right side of (1) is at least `b=497/1000` for every
`p,t`. A floating search suggests its exact infimum is about `.49747633`,
but that number is not used or asserted as a certified optimum.

## 2. Reduction to a compact one-dimensional check

Set `r=2b sqrt(p)` and `delta=1-r>0`. To have both terms in (1) below `b`
would require

```
t < T(p):=h(p)/(2delta),
F_p(t):=p log2+g_t(1)+delta t < 0.                 (2)
```

Differentiating the optimized Gaussian expression gives
`g'_t(1)=-lambda_t/t`. The function `F_p` decreases up to

```
t_*(p)=(1-delta)/[2delta(2-delta)]
```

and increases afterwards. Its unconstrained minimum is exactly

```
min_(t>0) F_p(t)=p log2+.25 log(1-4b^2 p).         (3)
```

For `0<p<=p0=93/100`, the expression on the right of (3) is nonnegative.
It is concave in `p`, is zero at `p=0`, and its value at `p0` is strictly
positive by the directed check below. This excludes (2) in that region.

For `p0<=p<1`, one has `h(p)<=h(p0)<3/10` and `r>9/10`. Therefore

```
h(p) < r/(1+r),       T(p)<t_*(p).
```

So `F_p(t)>=F_p(T(p))` throughout the time range in (2), and

```
F_p(T(p))=p log2+g_(T(p))(1)+h(p)/2.               (4)
```

The case `p=1` is immediate from the first term of (1), which is `1/2`.

## 3. The short interval verification

Partition `[.93,1]` into the 700 rational intervals
`[j/10000,(j+1)/10000]`, `j=9300,...,9999`. Since entropy is decreasing
there and `g_t` is decreasing in `t`, for an interval `[l,u]` one has

```
T(p) <= h(l)/[2(1-2b sqrt(u))] =: T_upper,
F_p(T(p)) >= l log2+g_(T_upper)(1)+h(u)/2.          (5)
```

`computations/principle_synthesis_2026_09_07_paired_certificate_floor.py`
evaluates every lower bound (5) with 50-decimal directed intervals, using
exact rational grid endpoints and the stable displayed formula for
`lambda_t`. Every lower endpoint is strictly positive. The smallest checked
lower endpoint is exactly

```
150747305007794028440063327974512410775466568163527
/11972621413014756705924586149611790497021399392059392
```

(approximately `.0125910`), on `[4881/5000,9763/10000]`.
The same run verifies the endpoint in (3) is positive and the two elementary
inequalities used for `T<t_*`. Its complete concise result is
`computations/results/principle_synthesis_2026_09_07_paired_certificate_floor.json`.

Combining the regions proves

```
inf_(0<p<=1,t>0) C_pair(p,t) >= .497.
```

This is already strictly above the original all-order upper endpoint below
`.493609`. Hence optimizing retention, temperature, recursive depth, or the
entire scalar latent channel **within this same paired Bellman certificate**
cannot perform asymptotically lossless landing of the original optimal value.

The result is not a lower bound for `min_{A,C}Q([[A,C],[-C,-A]])`, not a
selectable-child obstruction, and not a floor for stability-aware counting.
In particular the root's newer selector theorem can exclude concentrated
spectral words before counting extrema, and then impose an actual local-
stability penalty; its exponent need not remain above this Bellman envelope.
The point is the exact distinction between a successful subhalf family and
the completeness of its original certificate, not a revival of the false
Boolean half-floor.
