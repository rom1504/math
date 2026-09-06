# Fresh adversarial audit of the strongest original lower theorem

Date: 2026-09-05. This is a fresh reconstruction of the final implication,
not a reliance on the number of earlier positive audits. The conclusion
passes with the exact scope stated below. In particular the headline
lower bound can avoid the unbounded-operator smooth-transport extension.

## 1. Conclusions checked

Let `c_minus=liminf M_n/n^(3/2)`. Let `J_anchor` be the exact Gaussian
value of the 21-anchor finite polynomial construction, and let
`C_scalar=sup_(V,alpha) E|U 1{|V|<=alpha}| 1{|V|>alpha}`, where `V`
ranges over unit first-chaos variables in the marked-tree Gaussian space.
The checked conclusions are

`c_minus > C_scalar >= J_anchor`,

with

`J_anchor >= .430658179405528602724053804634711026327173238336325190455581`.

The strict increment above `C_scalar` is not numerically quantified.
Neither convergence of the original minima nor completeness of an
enlarged Gaussian response family is asserted.

## 2. Minimal sufficient matrix dependencies

It suffices to use the following proved fixed-construction statements:

* The joint moment/CLT theorem for finite injective marked-tree fields.
* Their polynomial direct-energy identity.
* For each fixed odd `r>=3`, the polynomial weighted unmarked identity
  and uniform local joint Gaussian limit with variance
  `v_i=(B (B^2)^(circ r) B)_ii`.
* Fixed-fraction principal spectral deletion, proved by the diagonal
  Grothendieck majorant.

Every analytic extension needed for the final conclusion can be made
**after deletion**, when `||B||op<=L` is fixed. Thus the stronger
unbounded-operator Sobolev transport theorem and the raw even-input
covariance-operator theorem are not necessary dependencies of this
particular final lower bound.

Here is the explicit extension argument. At fixed finite old family,
for bounded Gaussian-a.e.-continuous responses `F,H`, approximate `F`
in Gaussian `L2` by a polynomial `P`, and `H` by a polynomial `Q`.
Preserve odd/even parity by symmetrizing the approximants. Uniform root
CLT and higher polynomial moments imply

`max_i E|F(X_i)-P(X_i)|^2 -> E|F(gamma)-P(gamma)|^2`,

and the analogous assertion for `H`. The products with polynomials are
uniformly integrable, since the responses are bounded and every fixed
polynomial has bounded higher moments. Matrix Cauchy--Schwarz bounds
the normalized energy error by

`L [ ||F-P||_(2,n) ||H||_(2,n) + ||P||_(2,n) ||H-Q||_(2,n) ]`,

where `||K||_(2,n)^2=n^-1 sum_i E K(X_i)^2`. The spin multiplier in
`B[S H]` preserves this norm. First send `n` to infinity at each fixed
polynomial, then improve the polynomial approximations. The Gaussian
right side is `E F U H`, continuous in these same `L2` norms by isometry.
No derivatives of the bounded responses are required for this argument.

For the weighted identity, let `Z=B h_r(BS)`. At fixed `L,r`, its
rootwise fourth moments are uniformly bounded. Replacing a bounded
local weight `M` by a polynomial `P` costs at most

`L ||F||_infinity (max_i E|M(X_i)-P(X_i)|^4)^(1/4)`
`                    * (max_i E Z_i^4)^(1/4)`.

Use Gaussian `L4` approximation for `M`, then `L2` approximation for
`F`; after `P` is fixed, the factor `P(X_i) Z_i` has bounded second
moments by the same fourth-moment estimate. Bounded deterministic
weights `d_i` only multiply these bounds by their fixed supremum.
This proves the needed bounded-response extension of the polynomial
weighted theorem. Polynomial approximants are analysis devices, never
unbounded Boolean means.

## 3. Finite-anchor existence and exact arithmetic

The anchor features are distinct normalized even Hermite monomials
and orthonormal. Ancestor closure makes each feature a function of the
fixed anchor coordinates alone. Every innovation in their first-chaos
orthogonal complement is independent of the joint anchor family.
For two innovations of correlation `q`, the proposed map has covariance
kernel `K(q)=sum_l w_l q^l`, with `sum w_l=1` and
`sum l w_l=D_Z<1`. For all `-1<=q<=1`,

`1-K(q) <= D_Z (1-q)`.

Thus the map is a strict contraction on the complete unit sphere of
the first-chaos orthogonal complement; convexity of that sphere is not
needed. The constructed fixed point is a genuine unit Gaussian variable
on the same space. The degree-200 polynomial is used only in this
Gaussian existence argument. It is not evaluated as a spin mean.

The exact 21-anchor script and both interval-helper modules were read
again. The Hermite grouping subtracts exactly the anchor monomials from
innovation degree zero. Starting total degree at two has already
removed the edge constant, so it is not subtracted twice. The common
factor `4 phi(alpha)^2` cancels from the derivative ratio. The saved
Gaussian integral series bounds allow initially increasing Taylor
terms: only the omitted alternating tail must be decreasing, and the
asserted `x^2/2<degree+2` ensures that condition.

An independent fresh execution of
`fresh_finite_anchor_fixed_point_certificate.main()` recovered the
entire recorded interval, with its output-file operation redirected
to an in-memory stream. In particular

`D_Z in [.990872910662428464868565727570052064530720110786980942487571,`
`         .990872910662428464868565727570052064530720110786980942487572]`,

and the lower-bound interval is exactly the one in Section 1 and the
certificate artifact. No changed result file or numerical solver was
used in this replay.

## 4. Gaussian closure and all parameter orders

For an infinite bounded even mask `H`, finite conditional expectations
preserve its range and parity and approximate it in `L2`. Smooth these
finite masks if desired. The isometry makes
`J(H)=E|UH|(1-H)` continuous, with Lipschitz bound `2||H-K||2`.
For a finite mask, project `UH` onto a finite ancestor-closed coordinate
set containing the response coordinates. Its omitted Gaussian tail is
independent of those coordinates and mean zero, so the finite old
energy is exactly `E F W_j`, not an approximate identity.

At fixed `L`, bounded feasible finite means and Section 2 therefore
realize every target value `J(H)` to arbitrary accuracy. The headline
anchor lower bound follows first for every fixed bounded-operator class
with the **same** value `J_anchor`, and then for all signings by the
deletion argument in Section 7. This route avoids any growing-depth
matrix iteration and any unbounded-operator smooth approximation.

For the strict gain, the order is more specific: fix `L`, choose
`R=3 sqrt(log L)` and `t=a/(8 C R)`, choose one finite smooth target
approximation accurate compared with `1/log L`, and only then send
matrix order to infinity. Analytic polynomial approximations inside
the fixed-response identity are taken after their respective matrix
limits, as in Section 2. No degree, depth, cutoff, or smoothing width
varies during a matrix-size limit.

## 5. Exact Boolean feasibility and cancellation

Write `K=S H+epsilon D`, and use `mu_plus=F+K`,
`mu_minus=-F+K`. Outside the local update support, `|F|+H<=1`.
On that support `H=0`, `|F|<=1/4`, and
`epsilon |D|<=1/2`, so both means lie in the cube pointwise.
For a hollow symmetric `B`, conditional independent rounding preserves
their expected quadratic energies exactly. Algebra gives

`[E H_B(mu_plus)-E H_B(mu_minus)]/(2n)`
` = n^-1 E F^T B[S H] + epsilon n^-1 E F^T B D`.

The two `K^T B K` terms cancel, including marked/unmarked cross terms
and the quadratic perturbation. There is no discarded remainder and
no assumption of favorable signed cancellation between child blocks.
Each rounded energy lies in `[-Q(B),Q(B)]`, so their half-difference is
at most `Q(B)`. This checks the factor of two and orientation precisely.

## 6. Mean standard deviation and uniform scalar escape

For odd `r=2k+1`, set `P=Q^(circ(k+1))` and let
`Psi(X)=Q^(circ k) circ X`, where `Q=B^2` is a correlation matrix.
The map is unital, completely positive, and trace-self-adjoint, with
`Psi(Q)=P` and `Psi(P)=Q^(circ r)`. Square-root Jensen yields

`Tr Q sqrt(Q^(circ r)) >= Tr P^(3/2) >= n`.

Unit column norms of the symmetric square root `B` then give
`sum_i sqrt((B Q^(circ r) B)_ii)>=n`. This step needs no operator cap
or flatness assumption. With `d_i=1/sqrt(max(v_i,1/4))`, local variance
is at most one and `mean(d_i v_i)>=1/2`. Consequently the principal
weighted term is at least `a t`; clipping loses at most
`L sqrt(C_M t tau_R)`, with
`tau_R<=C_0(R+1)exp(-R^2/2)`. For the chosen `R,t`, the clipping loss
is at most `a t/2` for all sufficiently large fixed `L`. The exact
paired identity gives gain of order `1/log L` after the old loss
`O(t^2)`. All positive Gaussian constants are independent of `L`.

The uniformization over scalar near-optimizers also passes. Their mask
masses lie in `[9/50,5/8]`; `a1=E G0 F>=sqrt(2p/pi)-p>=1/200` on this
interval. Boundedness of `E[F|G0]` prevents all its higher odd Hermite
coefficients from vanishing. Compactness of the covariance triples,
including singular full triples, then selects a finite degree set with
a uniformly positive coefficient threshold. The `(V,W)` pair remains
uniformly nondegenerate by the height barrier, so the local feasible
rectangle has uniformly positive density. Non-realizable triples added
by closing the PSD parameter set are harmless: the same bound on `a1`
holds throughout that larger set. The chosen finite degree set precedes
`L` and all matrix limits.

Here is an independently certified elementary substitute for the
qualitative height nondegeneracy argument. On the mass interval above,
`1/5<alpha<9/10`. Fix `epsilon=10^-6`, `q=1-epsilon`, and put
`G'=qG+sqrt(1-q^2)Z` with independent standard normals. On
`alpha-.001<=G<=alpha`, the threshold on `Z` for `G'>alpha` is less
than `1+.0009<2`. Therefore

`P(|G|<=alpha, |G'|>alpha) >= .001 phi(.9) P(Z>2)>5*10^-6`.

The exact Gaussian primitives verify `phi(.9)>1/4`, `P(Z>2)>1/50`,
and both threshold endpoint comparisons. Dividing by `p<=1` gives
`1-K_p(q)>5*10^-6`, hence `q-K_p(q)>4*10^-6`. The height projection
lemma consequently gives the explicit uniform residual lower bound
`||V-U g_p(V)||2 >= 1/500000`. Reproduction is in
`computations/fresh_uniform_height_barrier_certificate.py`.

## 7. Principal deletion and the final strict inequality

For every hollow signing, exact polarization gives `beta(A)<=4Q(A)`.
The diagonal Grothendieck SDP supplies `D>=+/-A` with
`Tr D<=K beta(A)`. Delete coordinates above
`K beta(A)/(epsilon n)`. Fewer than `epsilon n` coordinates are deleted,
and the remaining principal signing has normalized operator bound
`L<=C_1/epsilon` on a low-cap minimizing sequence, with `C_1` independent
of `epsilon<=1/2` and of `n`. No coherence assumption is introduced.

Principal monotonicity is exact for hollow matrices: average independent
omitted spins, leaving the chosen internal quadratic energy unchanged.
Thus `Q(A)>=Q(A_sub)`, with no refill or bridge penalty. The retained
order is at least `(1-epsilon)n`, yielding the precise factor
`(1-epsilon)^(3/2)` in the limit.

For the non-strict headline lower bound, fix `epsilon`, take the matrix
limit using the fixed-`L` bound `J_anchor`, and finally let
`epsilon` decrease to zero. For scalar escape the bounded-operator
bound is `C_scalar+c/log L`, with the same positive `c` for all large
`L`. Choose a sufficiently small **fixed** `epsilon` so that

`(1-epsilon)^(3/2) [C_scalar+c/log(C_1/epsilon)]>C_scalar`.

Such a choice exists because `epsilon log(C_1/epsilon)->0`. Only then
take the matrix-size limit along a minimizing subsequence. The retained
orders need not exhaust all integers: the fixed-operator theorem is
universal for every growing sequence. This establishes a genuine
positive strict gap and requires no exchange of an infimum with a
growing algorithmic depth.

## 8. Audit outcome and limitation

No gap was found in the displayed final implication, its signs,
normalization, compactness, or order of limits. The polynomial
tree-energy and all-odd weighted matrix theorems remain substantive
dependencies; this audit does not replace their separate diagram and
matrix-chaos proofs. It does show that the final lower theorem does
not additionally require the strongest unbounded-operator response
extensions developed during the campaign.

The independently proved cosquare cap separation is consistent with
all these conclusions: it refutes a proposed completeness statement,
not a lower-bound identity. Original convergence remains unresolved.
