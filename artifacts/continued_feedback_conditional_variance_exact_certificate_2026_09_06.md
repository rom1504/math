# Exact negative conditional-variance-envelope certificate

Date: 2026-09-06. Exact integer/rational certificate for the root envelope
in `continued_director_conditional_variance_envelope_test_2026_09_06.md`.
The certificate below is unconditional as a numerical inequality. Its
signing consequence uses the separately proved/audited supersolution and
stopped-tree arguments; those are not inferred merely from this check.
The complete code and proof passed independent full read and exact replay;
the auditor obtained all three displayed rational endpoints unchanged.
The separately audited construction chain is recorded in
`continued_audit_exact_conditional_envelope_to_all_order_cap_2026_09_06.md`.

Let p=31/32 and t=4. Put

`F(z,s)=h(z)+z h((1+s)/2)+g_t((z-z^2 s^2)/p)`,

for 0<=z,s<=1, where h is binary entropy and g_t is the centered Gaussian
self-transport potential. The exact posterior reduction makes the root
offset

`-h(p)+t(1-sqrt(p))+sup_(E z=p) E F(z,s)`.                  (1)

The verified bound is

`(1) <= -91470529542342299/20460000000000000000`
`     = -0.004470700368638431... < 0`.                       (2)

## 1. Reproducible exact verifier

From the repository root run

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python computations/continued_feedback_conditional_variance_exact_certificate_2026_09_06.py --output computations/results/continued_feedback_conditional_variance_exact_certificate_2026_09_06.json
```

The run takes approximately twenty seconds in the campaign environment.
It uses rational precomputation and signed 64-bit INTEGER arrays, not
floating-point arrays or transcendental optimization. Floats in its
output JSON are displays only. The exact integer overflow conditions
are asserted before the calculation.

The imported logarithm routine encloses positive rational logarithms
by a 25-term atanh series with an explicit rational tail and outward
rounding to 10^-12. Entropy intervals are exact weighted sums of these
logarithm intervals. Square-root endpoints are obtained by integer
square root at the same scale and verified by squaring. These routines
already passed independent exact replay in the original latent-envelope
certificate.

## 2. Exact finite grid upper values

Use N=2500 and z=i/N,s=j/N. Precompute upper binary entropies at all
the z grid points and all `(1+s)/2` grid points. Also precompute upper
Gaussian potentials at variance k/1000, for 0<=k<=1033.

For positive variance v, write x=tv and

`r=(sqrt(1+16x^2)-1)/(4x)`,
`g_t(v)=-x(1-r)+(1/4)log(1-r^2)`.

Use the UPPER r endpoint in the increasing linear term and the LOWER r
endpoint in the decreasing logarithmic term. This gives a rational upper
bound for g_t(v); at zero its value is exactly zero. Round each table
entry upwards to an integer multiple of 10^-12.

At each grid pair, the variance lookup index is the exact integer

`k=floor[32000(i N^3-i^2 j^2)/(31 N^4)]`.

It is the floor of 1000 times the true variance. Since g_t decreases,
the upper potential at k/1000 is an UPPER value at the true variance.
No additional lookup error has to be appended: the optimistic lookup
is already included in every pointwise upper value.

The largest intermediate variance numerator is bounded by
`32*1000*N^4<2^63`. The entropy product integer is bounded by
`N*10^12<2^63`. The code also checks every lookup index. Thus the array
operations are exact signed 64-bit arithmetic without wraparound.
Division in the z-times-entropy product is rounded upwards.

For each z grid point retain the largest upper value over ALL 2501 s
points. Compute the exact rational least concave majorant of those
2501 values by a decreasing-slope stack. Evaluating it at p gives the
maximum of the grid upper averages over all grid laws with E z=p.
Including the outward root offset gives

`grid_offset <= -13307450618189/1056000000000000`
`             = -.012601752479345644...`.                   (3)

## 3. Passing from the grid to every posterior law

The Gaussian potential is t-Lipschitz in its variance: its derivative
is `-t(1-r)` in [-t,0]. The variance function has bounds

`|partial_z ((z-z^2 s^2)/p)|<=1/p`,
`|partial_s ((z-z^2 s^2)/p)|<=2/p`.

Round an arbitrary z RANDOMLY to its two adjacent grid points, preserving
its expectation exactly; its displacement is at most 1/N. Round s to a
nearest grid point, with displacement at most 1/(2N). The latter has no
mean constraint. These operations preserve the only posterior mixing
constraint E z=p. Binary entropy has modulus h(delta), so the pointwise
loss is at most

`omega=h(1/N)+(log2)/N+h(1/(4N))+2t/(pN)`.                  (4)

The h(1/(4N)) term is the entropy change of `(1+s)/2`, whose displacement
is at most 1/(4N). The two variance displacements each contribute at
most t/(pN). The factor z multiplying the spin entropy costs at most
(log2)/N. Thus every continuous posterior law is bounded by some mean-
preserving grid law plus omega, without assuming its support was on the
grid or that a sampled local maximum was global.

The exact outward evaluation of (4) is

`omega <= 630156538579809/77500000000000000`
`       = .008131052110707212...`.

Adding it to (3) gives exactly (2). All optimization and interpolation
errors have been paid before the strict negativity assertion is made.
The certificate covers the FULL posterior concave envelope in (1), not
only the two support points observed in a floating diagnostic.
