# Exact local-stability wind-tunnel

This is a small sampling diagnostic, **not** a certified expectation or a
new cap search. All151 per-matrix histograms are exact:144 random twists,
24 per group, plus seven known family-minimum witnesses. Random seed is
2026091904. Runtime was3.48 seconds, below the authorized five-minute ceiling.
Every input child, signed permutation, matching and parent matrix is saved in
`computations/results/twisted_chiral_stability_wind_tunnel_2026_09_19.json`.

The C++ evaluator enumerates all projective parent spins with coordinate0
fixed positive, maintaining exact integer energies and local fields by Gray
updates. A positive one-spin local maximum requires x_i(Dx)_i>0 at **every**
coordinate, including0. These full even-order sign matrices have odd fields,
so there are no ties. Six known-minimum histograms were independently replayed
with direct batched NumPy matrix products, not Gray updates, and agree exactly.

All three thresholds are strict. Integer comparisons are e²>8Q(A)²,
3e>4β(A), and e>the separately certified fixed-child family minimum.
Known witnesses are excluded from sample means. A positive energy above a
threshold exists exactly when a positive stable energy above it exists:
strictly increasing one-spin ascent terminates at a local maximum. This
equivalence was also checked for every sampled matrix and threshold.

| Child group | Q(A) | β(A) | Mean raw count above2√2Q | Mean stable count | Stable/raw |
|---|---:|---:|---:|---:|---:|
| Optimal6, class0 |5|12|60.000|7.750|0.1292|
| Optimal8, class0 |10|24|76.083|10.917|0.1435|
| Optimal8, class1 |10|24|82.750|11.042|0.1334|
| Optimal10, class0 |13|40|1338.625|35.958|0.02686|
| Optimal10, class1 |13|34|1348.375|36.208|0.02685|
| Conference10 |15|30|121.292|18.833|0.1553|

All144 sampled twists violate the2√2Q threshold. Filtering is material,
especially the roughly37-fold reduction for both optimal order10 classes,
but every corresponding sample mean remains well above1.

| Child group | Mean raw above4β/3 | Mean stable above4β/3 | Draws with violation | Mean stable above family minimum |
|---|---:|---:|---:|---:|
| Optimal6, class0 |28.167|6.875|24/24|4.708|
| Optimal8, class0 |15.375|4.833|24/24|7.500|
| Optimal8, class1 |17.042|5.458|24/24|8.042|
| Optimal10, class0 |3.583|1.542|20/24|14.625|
| Optimal10, class1 |107.500|15.208|24/24|15.208|
| Conference10 |264.708|28.458|24/24|28.458|

No requested sampled stable-count mean is below1. The easiest observed
resource threshold is4β/3 for optimal10 class0, but its sample mean1.542
does not certify a first-moment bound. The two optimal order10 classes
have different exact bilinear norms40 and34 despite equal Q13 and equal
twisted-family minimum44.

Reproduce with `.venv/bin/python computations/twisted_chiral_stability_wind_tunnel_2026_09_19.py`.
The driver auto-builds the retained C++ source if missing or stale. It is
restricted to full even-order signing matrices of order at most20, and
should not be used with the same strict-stability test on matching-free
cores, whose fields may vanish.
