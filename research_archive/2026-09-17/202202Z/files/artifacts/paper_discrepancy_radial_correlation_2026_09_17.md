# Exact radial sign-correlation laws on stored finite minimizers

2026-09-17 campaign. These are finite, full-cube linear programs with
replayed rational certificates. Global-minimizer provenance is imported
from the archived solver certificates; matrix caps and every displayed
correlation/response inequality are checked again. No asymptotic
optimizer structure is inferred.

## 1. The radial problem and exact certificate convention

For a fixed hollow full signing A of order n and direction s∈{±1},
maximize α subject to a probability law ν on physical sign vectors h
satisfying

$$
\mathbb E_\nu hh^T=I+s\alpha A.                       \tag{1}
$$

All 2^(n−1) projective sign columns are included. A projective law may
be globally symmetrized if mean zero is desired; this changes neither
its correlations nor any absolute response. Put t=α√n.

The script reconstructs rational primal probabilities and a rational
quadratic polynomial

$$
P(h)=\alpha+\sum_{i<j}b_{ij}h_ih_j\ge0
\quad\hbox{for EVERY sign }h,\qquad
-s\sum_{i<j}b_{ij}A_{ij}=1.                          \tag{2}
$$

Taking expectations in (2) under a feasible law of radial value α′
gives α−α′≥0. Matching the primal proves the exact optimum.

The finite directional energy bound is

$$
\alpha\le\frac{\max_x sH_A(x)}{\binom n2},\qquad
t\le\frac n{n-1}\frac{2\max_x sH_A(x)}{n^{3/2}}.     \tag{3}
$$

The factor n/(n−1) must not be omitted in finite comparisons.
For comparison, Gaussian sign rounding from covariance I+sρA,
ρ=−1/λ_min(sA), gives the feasible radial value
`α_G=(2/pi) arcsin(ρ)`. The linearized common operator-norm baseline
is `t_G≈(2/pi)/||A/sqrt(n)||_op`. Neither comparison characterizes
the full correlation polytope.

## 2. Exact finite results

All 30 directional optima for the 15 stored cases from n=3 through14
were certified. The even-order cases are:

| Stored case | Exact α₊=α₋ | t=α√n | Exact arcsine baseline t_G |
|---|---:|---:|---:|
| n4 | 1/3 | 0.666667 | 0.590334 |
| n6 | 1/3 | 0.816497 | 0.723009 |
| n8 class0 | 1/4 | 0.707107 | 0.506041 |
| n8 class1 | 1/3 | 0.942809 | 0.611921 |
| n10 | 1/5 | 0.632456 | 0.493183 |
| n12 | 1/4 | 0.866025 | 0.566877 |
| n14 | 3/13 | 0.863459 | 0.669429 |

Only n6 and n14 in this table saturate (3). Thus nonlocal physical
sign laws substantially exceed the Gaussian-sign baseline in these
actual finite examples, but the optimal t is not one common constant.

The replay is
[radial_correlation.py](../computations/paper_discrepancy_2026_09_17_radial_correlation.py).
The full 30 rational primal and dual witnesses are frozen in
`tmp/paper_portfolio_2026_09_17/discrepancy/radial_correlation_certificates.json`.
The `--full` option regenerates them.

## 3. Covariance cancellation does not determine absolute response

Given maximal radial laws ν₊,ν₋, mix them with weight
θ=α₋/(α₊+α₋) on ν₊. The resulting covariance is exactly I.
We separately minimized its maximum absolute response on the COMPLETE
absolute ground code, retaining both radial endpoint constraints.

| Case | Best response at fixed radial endpoints |
|---|---:|
| n4 | 4/3 |
| n6 | 5/3 |
| n8 class0 | 7/4 |
| n8 class1 | 5/3 |
| n10 | 747/335 |
| n12 | 9/4 |
| n14 | 98/39 |

Every value now has matching rational primal laws and a rational
all-sign-query lower certificate. At n12, direct rational rounding
initially failed and produced the certified interval
[9/4−33/(8·10¹¹),9/4]. Exact row reduction on its 152-variable active
dual face subsequently recovered an endpoint certificate, checked
against all columns in both radial sectors. The earlier failure is
preserved in the output status history. This final equality follows
from that exact certificate, not from the decimal LP value.

At n8 class0, the initially selected, exactly certified pair of laws
has maximum ground response 35/16. A different pair at the SAME
covariances I±A/4 has response 7/4. This is a concrete exact warning
against replacing a physical sign law's absolute response by its
covariance or by an assumed Gaussian formula.

Full witnesses are in
`tmp/paper_portfolio_2026_09_17/discrepancy/radial_correlation_mixture_certificates.json`.
Use `--even-only --optimize-mixture --full` to regenerate them.
The n12 exact endpoint recovery is separately replayed by
`computations/paper_discrepancy_2026_09_17_radial_dual_face.py`, with
its supplemental witness in
`tmp/paper_portfolio_2026_09_17/discrepancy/radial_n12_dual_face_certificate.json`.

## 4. Averaging against the radial laws' own supports can miss the obstruction

There is an explicit full-sign family showing why a universal
order-√n response lower bound cannot follow merely by averaging a
balanced radial mixture against its own sector supports.

Let H=(J₄−2I₄)^(tensor a), n=4ᵃ, d=(−1)ᵃ, and A=H−dI.
Then H²=nI, its diagonal is d, and A is a hollow full signing. The
tensor Walsh basis consists of n mutually orthogonal Boolean
eigenvectors of H. Let ν₊ and ν₋ be uniform on the respective basis
subsets with eigenvalues +√n and −√n. Their sizes are

$$
r_\pm=(n\pm d\sqrt n)/2.
$$

Writing P_±=(I±H/√n)/2 for the eigenprojections gives

$$
\mathbb E_{\nu_+}hh^T=\frac n{r_+}P_+
=I+\frac{A}{\sqrt n+d},\qquad
\mathbb E_{\nu_-}hh^T=\frac n{r_-}P_-
=I-\frac{A}{\sqrt n-d}.                             \tag{4}
$$

Both radial values are optimal: the corresponding Boolean eigenwords
attain the spectral directional cap, and (3) matches (4).
Mix the laws with weights r₊/n,r₋/n. This is uniform on the full
orthogonal Boolean basis and has covariance I. At EVERY word x in
either of those own supports,

$$
\mathbb E|h\cdot x|=1,                              \tag{5}
$$

because exactly one of the n basis vectors has overlap magnitude n
and all others have overlap zero. The entire absolute near-ground
code may contain many other words with larger response; (5) says
nothing about their maximum. It specifically defeats the proposed
own-sector averaging shortcut, without using a non-full or weighted
matrix and without pretending that covariance determines L¹ norms.
