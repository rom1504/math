# Independent audit of the variance-normalized localized gain

Date: 2026-09-05. This checks the proposed strict-improvement mechanism,
including its finite/infinite approximation order. It is not a new decimal
certificate and does not prove convergence of the original minimax sequence.

## 1. Weighted identity and its analytic scope

I read the complete files
`fresh_full_weighted_unmarked_projection_2026_09_05.md` and
`fresh_tree_matrix_derivative_and_weighted_third_chain_2026_09_05.md`.
The matrix-chaos import was also checked directly against Theorems 2.1
and 2.4 of the primary paper
<https://web.math.princeton.edu/~rvan/chaosconf241224.pdf>.
The former compares square-free coupled and decoupled matrix chaoses;
the latter bounds the decoupled expectation by a logarithmic factor times
the largest coefficient flattening with matrix indices on opposite sides.
The derivative module checks precisely those flattenings, and its
dimension-free fixed-moment upgrade suffices for the subsequent Holder
bounds.

The algebraic split in the weighted proof passes. Mixed derivative
allocations have one flat Hadamard factor and hence a Frobenius gain.
The second-chain remainders retain the needed powers of `Q_ij`. The
single-third-chain exception uses the exact one-level own-coordinate
exclusion before Gaussian integration by parts. In particular it does
not replace a pointwise `O(n^-1/2)` third derivative by a nonexistent
small operator norm. The all-edge surviving term has the correct trace
`sum_i v_i`, where `v_i=(B Q^{circ3} B)_{ii}`.

An arbitrary deterministic multiplier `d_i` with uniformly bounded
absolute values can be inserted at the local-weight root. It does not
receive a spin derivative, and all diagonal operator bounds gain at most
its fixed supremum. The main term becomes

`E M * E[F h3(G)] * n^{-1}sum_i d_i v_i`.                       (1)

Taking `d_i=1/sqrt(max(v_i,eta))` is legitimate at a fixed `eta>0`, even
though the multiplier depends on the signing. No random-data derivative
is hidden in this choice.

## 2. A pointwise feasible curve with exact energy cancellation

Let `H` be an even mask in `[0,1]`, let `W` be an old Gaussian
preactivation, and put

`F_t=(1-H) psi(W/t)`.

Choose a fixed smooth odd `psi`, bounded by one, agreeing with sign outside
a compact interval, and satisfying `|psi(u)|<=1-s_0` for `|u|<=a`, where
`0<a<1` and `s_0>0`. Let `0<=M_t<=1` be even, supported where
`H=0` and `|W|<=a t`. Set

`D_i=M_t(X_i) clip_R(d_i sign(b3) Z_i)`,

`mu_+=F_t+S H+epsilon D`, `mu_-=-F_t+S H+epsilon D`.

For `epsilon<=s_0/R`, both means belong to `[-1,1]^n` for every input.
Indeed on the perturbation support `H=0` and the old unmarked mean has
slack `s_0`; elsewhere the perturbation vanishes. Independent conditional
rounding is valid because the objective is hollow.

The half energy difference is exactly

`n^{-1} E F_t^T B[S H] + epsilon n^{-1} E F_t^T B D`.           (2)

There is no quadratic perturbation penalty. Both the `D^TBD` term and the
cross term with `S H` cancel. This uses two means for the same signing,
not any assumption about cancellation between blocks or child spins.

## 3. Main constants and cutoff dependence on the operator cap

Assume for the fixed target Gaussian construction that

* `b3=E[(1-H) sign(W) h3(G)]` is nonzero;
* `E M_t` lies between two positive constant multiples of `t`;
* the old soft-sign loss is at most `C_s t^2`.

These constants can be kept uniform in a small covariance neighborhood of
a target nondegenerate Gaussian tuple. Bounded density of `W` gives the
quadratic soft-sign loss. Positive conditional probability of the chosen
`H=0` region given `W=0` gives the linear local-mask mass. Continuity of
the cubic coefficient gives `|E F_t h3(G)|>=|b3|/2` for all sufficiently
small `t`. A nonzero *nonlinear* conditional edge projection does not by
itself imply `b3!=0`; the cubic coefficient is a separate condition.

The odd-Schur standard-deviation theorem gives the uniform floor

`n^{-1}sum_i d_i v_i >= 1-sqrt(eta)/4`.

Thus (1) gives an unclipped direction at least `c t`, with `c>0`
independent of the fixed operator cap `L`. In the matrix limit, the
normalized local unmarked field is Gaussian of variance at most one and
independent of the old local fields. Therefore the limiting mean-square
cutoff error is bounded by

`E M_t * tau_R`, `tau_R=E[N^2 1{|N|>R}]`.

The corresponding normalized energy error is at most

`L sqrt(E M_t * tau_R)`.

The only displayed dependence on `L` is the operator Cauchy--Schwarz
factor. Uniform integrability may have constants depending on `L` and the
finite family, but it is used only before taking the fixed-parameter
matrix limit, and introduces no surviving extra constant into this
Gaussian tail estimate.

For `R>=1`, `tau_R<=C(R+1) exp(-R^2/2)`. For all sufficiently large `L`,
take `R=3 sqrt(log L)` and `t=c_0/R`. The cutoff error is then negligible
relative to `t`. Equation (2), with `epsilon=s_0/R`, has positive gain
of order `c_0/R^2`, while the soft-sign loss is `O(c_0^2/R^2)`. Choose
`c_0>0` small once and for all. The net gain is at least

`c_*/log L`,                                                   (3)

where `c_*>0` is independent of `L`. This calculation is valid for large
fixed `L`, not for `L=L_n` inside a single matrix-limit theorem.

## 4. Finite realization of an infinite Gaussian baseline

The following concrete approximation avoids differentiating a hard mask
through a matrix limit. Suppose the target uses a first-chaos unit
Gaussian `V`, a central indicator `H_*=1{|V|<=alpha}`, and `W_*=U H_*`.

Choose finite-coordinate first-chaos `V_j -> V` in `L^2`. Use smooth even
compact-transition masks `H_j(V_j)`, valued in `[0,1]`, converging to
`H_*`, and equal to zero outside `|V_j|<=alpha+delta_j`, with
`delta_j->0`. Choose a finite coordinate projection `P_j` containing the
coordinates of `V_j` and the edge coordinate `G`, and enlarge it so that

`W_j=P_j U H_j -> W_*` in `L^2`.

The old responses and local masks are then smooth functions of finitely
many old fields. The final bounded clip is used only through its tail
error, not differentiated in the old energy theorem. The omitted
first-chaos tail `U H_j-W_j` is independent
of all response coordinates and has mean zero. Consequently its old
paired energy is exactly

`E F_{t,j} U H_j=E F_{t,j} W_j`.

Choose a fixed smooth local mask supported in
`|V_j|>=alpha+1` and `|W_j|<=a t`, so `H_j=0` on its support for all large
`j`. If the target pair `(V,W_*)` is nondegenerate, its conditional
Gaussian tail at `W_*=0` is positive. The local-mask mass bounds and
soft-sign constants are therefore uniform for all sufficiently large
`j`. If the target cubic coefficient is nonzero, its lower bound is
uniform as well, by convergence of the Gaussian covariance tuple
`(V_j,W_j,G)` and uniform fixed moments.

The required order is:

1. Fix a sufficiently large operator cap `L` and the resulting `R,t`.
2. Choose `j` large enough that the Gaussian baseline error is less than
   a fixed fraction of the positive gain in (3).
3. Fix every coordinate, mask transition, cutoff, and perturbation step.
4. Only then let the matrix dimension tend to infinity.

Derivative and matrix-chaos constants are allowed to depend on `j,L,t,R`
in step 4. They multiply errors vanishing with dimension and do not alter
the universal constants in the limiting Gaussian calculation. This
resolves the otherwise invalid inference from a gain for each finite
approximation to a gain over the limiting baseline.

## 5. Transfer and remaining scope

Classical diagonal regularization gives a principal-submatrix operator
cap polynomial in the reciprocal deleted fraction. The loss in the
original `n^{3/2}` normalization is linear in that fraction. Equivalently,
the full-order random-refill regularization has loss `O(L^{-1/2})`.
Either loss is eventually smaller than (3). Thus, once the target
Gaussian baseline satisfies the explicit cubic and local-mass conditions,
the scheme gives a genuine universal strict improvement over that
baseline, not merely a strict improvement for each fixed bounded-op
class.

This implication supplies no numerical value for the strict gap without
additional quantitative choices. It also does not prove that the gain is
uniform over every maximizing mask in the entire old variational class,
nor does it decide convergence of `M_n/n^{3/2}`.

## 6. Final quantitative-normalization cross-check

I subsequently read the full final statement
`fresh_normalized_unmarked_gain_and_limit_order_2026_09_05.md` and checked
its constants independently. With its notation,
`epsilon=1/(2R)` is feasible, `t=a/(8CR)` gives the net pre-approximation
gain `a^2/(64CR^2)`, and allocating half of that gain to the finite
Gaussian approximation leaves `a^2/(128CR^2)`. Since
`R^2=9 log L`, the stated constant `a^2/(1152C)` is correct.
The principal-deletion normalization contributes exactly
`(r/n)^{3/2}`, not an additional square-root loss. Its polynomial operator
cap is enough for the logarithmic gain to dominate deletion.

The compact-neighborhood argument only requires nondegeneracy of the
pair `(V,W)`. The full triple `(V,W,G0)` need not be invertible: continuity
of the cubic expectation follows from Gaussian coupling and uniform
moments, since the response's `V,W` boundary events still have probability
zero. The instance-specific nonzero cubic check is supplied separately by
the exact finite-anchor cubic-coefficient certificate, not inferred here.
The complete final implication therefore passes this independent audit.
