# Quantifier audit: variance normalization and the localized gain

Date: 2026-09-05. This records a conditional implication whose analytic
inputs are listed explicitly. It is not a standalone numerical certificate.

## 1. Fixed Gaussian baseline assumptions

Suppose a fixed finite old-tree Gaussian construction supplies an even
mask `H` in `[0,1]`, a centered Gaussian preactivation `W`, and the hard
response `F_0=(1-H) sign(W)`. Assume the following fixed constants exist.

1. The old paired energy is `J_0=E[W(1-H) sign(W)]` (or the corresponding
   exact finite projected preactivation identity).
2. The cubic edge coefficient `b_3=E[F_0 h3(G)]` is nonzero.
3. For all small fixed `t>0`, a local mask `0<=M_t<=1` can be chosen with
   support in `{H=0, |W|<=a t}`, where `a<1` is fixed, and
   `c_0 t <= E M_t <= C_0 t`.
4. An odd bounded smooth sign approximation `psi` has
   `sup_(|s|<=a)|psi(s)|<=1-s_0`, with `s_0>0`, and its paired smoothing
   loss for `F_t=(1-H)psi(W/t)` is at most `C_s t^2`.

For a usual bounded transition equal to sign outside a compact interval,
assumption 4 follows from the bounded density of the nondegenerate
Gaussian `W`. Also `F_t -> F_0` in Gaussian `L2`, so after possibly
decreasing the allowed range of `t`, its cubic coefficient has the same
sign and magnitude at least `|b_3|/2`.

The support/probability hypothesis in 3 must be proved for the particular
baseline. It cannot be inferred merely from `P(H=0)>0`.

## 2. Exact paired perturbation and its feasibility

For an actual signing with fixed `||B||op<=L`, normalize the cubic transport
by `d_i=1/sqrt(max(v_i,eta))`, with one fixed `0<eta<1`, and let
`Ztilde_i=d_i sign(b_3) Z_i`. Define

`D_i=M_t(X_i) clip_R(Ztilde_i)`,

`mu_+ = F_t + S H + epsilon D`,

`mu_- = -F_t + S H + epsilon D`.

Take `0<epsilon<=s_0/R`. On the support of `D`, the old marked mean `H`
is zero and `|F_t|<=1-s_0`; off that support `D=0`. Hence both means
lie in the Boolean cube. This is a pointwise check, not an average one.

Their half energy difference is exactly

`(E mu_+^T B mu_+ - E mu_-^T B mu_-)/(4n)`

`= n^-1 E F_t^T B[S H] + epsilon n^-1 E F_t^T B D`.       (1)

Every quadratic term in `D`, and its cross term with `S H`, cancels.
No reversal of the signing and no second-derivative energy estimate are
needed. The absolute quadratic maximum bounds the left side of (1).

## 3. Uniform positive channel and the only cutoff loss

The weighted projection and mean-standard-deviation inequality give

`liminf_n n^-1 E F_t^T B[M_t circ Ztilde] >= kappa t`,      (2)

where `kappa=c_0 |b_3|(1-sqrt(eta))/2>0` is independent of `L`.
The dimension limit is taken at fixed `L,t,eta` and fixed baseline.

The local normalized field has limiting Gaussian variance at most one
and is independent of the old fields. Therefore, uniformly along roots,
the mean square cutoff error is at most

`C_0 t * E[N^2 1_(|N|>R)] + o(1)`.

Bounded-operator Cauchy--Schwarz bounds the corresponding energy error by

`L sqrt(C_0 t) * [E N^2 1_(|N|>R)]^(1/2) + o(1)`.       (3)

For `R>=1`, the standard Gaussian tail is at most a fixed constant times
`(R+1) exp(-R^2/2)`. Thus one can choose

`R <= C_1 sqrt(log(e+L/t))`

with a baseline-dependent `C_1`, large enough that (3) is at most
`kappa t/2`. The logarithmic factor absorbs `L`; the unnormalized field
would have needed a radius proportional to `L`.

Combining (1)--(3), with `epsilon=s_0/R`, gives a gain over `J_0` of at least

`c_1 t/sqrt(log(e+L/t)) - C_s t^2`,                    (4)

for a positive constant `c_1` independent of `L`.

## 4. A gain that decays only logarithmically in the operator cap

Let `ell=log(e+L)`. Choose `t=a_0/sqrt(ell)`, with a sufficiently small
fixed `a_0>0`. It is within the allowed range for every sufficiently large
`L`, and

`log(e+L/t) <= C_2 ell`.

The positive term in (4) is then at least
`c_1 a_0/(sqrt(C_2) ell)`, while its negative term is
`C_s a_0^2/ell`. Taking `a_0` small enough proves

`gain >= c_*/log(e+L)`                                (5)

for a fixed `c_*>0` and all sufficiently large fixed `L`.

This order is sufficient to dominate a polynomial spectral-deletion
loss. Precisely, if the bounded-op reduction deletes at most an `epsilon`
fraction and produces a cap `L<=C epsilon^-p`, for fixed positive `C,p`,
then the lower bound after deletion is

`(1-epsilon)^(3/2) [J_0+c_*/log(e+C epsilon^-p)]`.

As `epsilon` decreases, the gain is of order `1/log(1/epsilon)` while the
loss is `O(epsilon)`. A sufficiently small but fixed epsilon therefore
makes the displayed quantity strictly exceed `J_0`. Then fix epsilon,
its cap, the smoothing and cutoff parameters, and the finite construction
before taking the matrix dimension to infinity.

## 5. Infinite Gaussian closures

If `J_0` is an infinite-coordinate Gaussian variational closure, a gain
for each finite approximation does not alone prove a gain over `J_0`.
To use (5), the finite approximations must retain uniformly positive
constants `|b_3|,c_0,s_0` and uniformly finite `C_0,C_s`, while their paired
energies tend to `J_0`. These are distributional conditions on a finite
number of preactivations/cutoff events and may hold for an explicit
Gaussian construction, but require a separate check.

If those uniform conditions hold, first fix the small deletion fraction
and hence a positive lower bound on the gain in (5), then choose one
finite approximation whose energy loss is less than, say, one quarter
of that gain. All matrix parameters are then fixed, and the strict
improvement survives. Reversing this order without uniform constants
would be unjustified.

## 6. Final audit of the selected anchored construction

The completed proof in
`fresh_normalized_unmarked_gain_and_limit_order_2026_09_05.md`, Sections
3--7, passes this independent end-to-end audit. Its finite smooth
approximations establish the uniform constants required in Section 5:
the covariance triples converge to a nondegenerate `(V,W,G0)` triple;
the smooth masks are identically zero on a fixed positive-density
rectangle beyond their transition strip; and the cubic coefficient
converges to a nonzero value. The finite coordinate projection may be
enlarged to an ancestor-closed family before defining the projected
preactivation, without affecting the argument.

The arithmetic in its logarithmic gain is exact. With
`R=3 sqrt(log L)` and `t=a/(8 C R)`, the clipping-adjusted gain before
finite approximation is `a^2/(64 C R^2)`. Spending one half on the finite
approximation leaves `a^2/(1152 C log L)`. Each parameter is fixed before
the dimension limit. A sufficiently small fixed deletion fraction then
gives a strict unrestricted gain.

The diagonal Grothendieck regularization proof in
`fresh_range_and_spectral_regularization_2026_09_05.md`, Section 4, was
also read and checked independently. Its semidefinite dual uses Gram
vectors with unit combined norm, and the elementary Grothendieck bound
gives a dominating diagonal of trace at most `K beta(A)`. Removing the
entries above its Markov threshold gives the stated operator cap
`4 K Q(A)/(epsilon n)`. Relative to the retained order this is a fixed
cap of order `1/epsilon`, exactly what the logarithmic comparison needs.

Finally, the exact Fraction script
`computations/fresh_finite_anchor_cubic_coefficient_certificate.py` was
independently rerun. Its interval is

`[-0.063914331766535141564927573349386587491619191945279997110215,`

` -0.063914331766535141564927573349386587491619191945279997110143]`.

The closed formula was separately reconstructed from
`G0=rho V+k Z+independent residual` and the four Gaussian Hermite terms
of its conditional cubic. The boundary term, the two scalar Gaussian
integrals, and the final `2/sqrt(6)` normalization all agree with the
script. Its positive conditional variances also verify the required
nondegeneracy.

For a more quantitative slice than the proof needs, the selected
`alpha=361/500` admits

`M_t=1{3/4<=|V|<=4/5, |W|<=t/2}`, for `0<t<=1/100`.

The certified covariance lies within
`p<.531`, `w<.701`, `.039<sigma^2<.04`. The joint Gaussian density on
this rectangle is greater than `1/200`: its exponent is below `4.5`,
its prefactor is greater than `35/44`, and `exp(4.5)<100`. Its area is
`t/10`, giving `EM_t>=t/2000`; the upper bound `EM_t<=t` follows from
the Gaussian `W` density. These strict bounds persist under sufficiently
accurate finite approximation. A smooth cutoff on a slightly smaller
rectangle gives the same property with a slightly smaller fixed constant.

The resulting conclusion is a genuine strict unrestricted improvement
over the anchored Gaussian certificate. This audit does not certify a
numerical increment and does not imply convergence of the original
minimum sequence.
