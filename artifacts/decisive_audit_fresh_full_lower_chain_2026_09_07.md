# Fresh reconstruction of the .4333221116640807 original lower bound

Status: **PASS** at the fixed-construction, ordered-limit scope below.
This audit reconstructs the mathematical implications and reads the exact
certificate implementation; earlier audit verdicts are not proof steps.
The exact replay is byte-identical to the recorded certificate.  No actual
gap was found.  This is a lower bound, not a convergence result.

Write `H_A(x)=x^T A x/2` and `Q(A)=max_{x in {+-1}^n}|H_A(x)|` for a
hollow symmetric signing.  The conclusion is

    liminf_n M_n/n^(3/2)
      >= .433322111664080753415812928897579346558634648033693413106996
      >  .4333221116640807.

## 1. Removing the operator restriction does not remove edges by fiat

Exact hollow polarization gives `beta(A)=max_{x,y}|x^TAy|<=4Q(A)`:
use `(x+y)/2,(x-y)/2` in the cube, whose quadratic values are bounded
by the vertex maximum by independent rounding.  The diagonal
Grothendieck majorant has

    D >= A,-A,       Tr D <= K_G beta(A).

The simultaneous two-sided constraint is important.  Its SDP dual has
positive matrices `X_+,X_-` with `diag(X_++X_-)=1`.  Their Gram vectors
give unit vectors `(u_i^+,u_i^-)` and `(u_j^+,-u_j^-)`; the ordinary real
bilinear Grothendieck inequality bounds the dual objective by
`K_G beta(A)`.  Thus this is one simultaneous majorant, not two unrelated
one-sided deletions.

On a sequence with `Q(A)<=C n^(3/2)`, delete coordinates with
`D_ii>4K_G C sqrt(n)/epsilon`.  At most `epsilon n` coordinates are
deleted.  The retained principal signing of order `n'` has

    ||A_sub/sqrt(n'-1)||op <= C'(C)/epsilon

for fixed `epsilon<=1/2` and large `n`.  Its entries are still exactly
signs, and its row normalization is exactly `n'-1`, not the old order.
Principal cap monotonicity is exact: average independent omitted spins
while fixing the retained witness.  It incurs only the final factor
`(n'/n)^(3/2)>=(1-epsilon)^(3/2)`.  Sequences with unbounded normalized
cap already satisfy the desired lower bound.  Therefore it suffices to
prove the same Gaussian numerical value for every fixed operator cap
`||B||op<=L`, where `B=A/sqrt(n-1)`, and only then send `epsilon` to zero.

This reconstructs the deletion in
`fresh_final_lower_dependency_algebra_audit_2026_09_05.md`, Section 7.

## 2. The finite old Gaussian frame comes from actual injective trees

Use trees with one external root of degree one and every other vertex
of odd degree.  A tree with `d` nonroot vertices contributes the
injective multilinear field

    X_(T,i)=1/(sqrt(a(T))(n-1)^(d/2))
      sum_(f injective, f(root)=i) product_edges A_(f(u),f(v))
                                      product_(v not root) S_f(v).

The output is exactly independent of its own spin.  At a fixed moment
order, only pairings of all marked labels can lead.  Smaller label
counts lose `O(1/n)`.  A nonempty parity-edge graph has an edge whose
two free endpoints can be summed with the bilinear discrepancy bound,
costing `beta(A)/n^2`.  If the parity graph is empty, connectedness and
the exact edge count force a tree with every edge doubled.  Each marked
label has only two occurrences, so pairing propagates along the whole
tree: it cannot switch to a third copy.  Rooted automorphism factors
give the stated variance one.  The same argument at two selected roots
has cross-covariance `1_(T=U)(B^2)_ij`, retaining all its possible
coherence.  Gaussian moment determinacy justifies the uniform finite
local Gaussian law.

For a tree `T`, removing its external edge gives the even normalized
child Hermite monomial `h_T`.  Distinct trees give distinct such
monomials and every finite even Hermite monomial occurs.  Consequently

    U h_T = G_T

extends to an isometry from even Gaussian `L2` onto first Gaussian
chaos.  The direct marked energy identity is obtained by adding the
bridge edge and the explicit root spin to the same moment expansion.
Its only doubled-tree main term pairs that bridge with the top edge
of a selected response tree; the remaining children produce `h_T`.
This verifies the isometry's normalization rather than postulating a
state-evolution rule.

Read/reconstructed sources:
`fresh_tree_chaos_literature_audit_2026_09_05.md`, Sections 1--6;
`fresh_limit_hierarchical_tree_energy_2026_09_05.md`, Sections 1--3;
`fresh_limit_injective_input_gram_2026_09_05.md`, Sections 1--7.

## 3. Nonlinear return: trace norm, not an operator-norm inference

For a fixed odd polynomial `F` in a finite ancestor-closed old frame,
write `F=P1F+R`, `K_F=U^{-1}P1F`, and
`tau^2=||R||2^2`.  Every local Hermite monomial in `R` has odd local
degree at least three.  Its *original input-spin degree* can be much
larger; these two degrees are never identified.

The needed fixed-root tree tensor cuts are `O(n^(-1/2))`, while every
global cut including the output index is `O(1)`.  Both follow by
compressing the crossing-edge tensor product of `B` and charging each
crossing-isolated vertex to an internal tree edge.  With the root
fixed, its incident edge supplies the additional `n^(-1/2)`.  Deleting
collisions has respectively `O(n^(-1/2))` or `O(1)` Hilbert norm, as
required for these two different statements.

For a transported product of at least two branches, classify a proper
marked-slot cut by the number of split branches.  Two split branches
give two small factors and an `O(sqrt(n))` row absolute sum.  No split
branches gives `C_left^T diag(B_i) C_right`.  One split branch factors
through the remaining whole-branch root map.  Thus every proper cut
of the transported forest is small.  This does **not** apply to a
single old branch, whose return can retain non-Gaussian coherent atoms.

The normalized nuclear covariance statement is

    || Cov(R(X))-sum_(odd k>=3) w_k (B^2)^(circ k) ||_*/n -> 0,

where `w_k` are the squared local Hermite coefficients grouped by local
degree.  In a two-forest pairing, any component which is neither a
whole edge nor a star has a four-vertex path giving two disjoint proper
cut gains.  Two nontrivial stars likewise give two gains.  The sole
remaining partial case is one star plus whole edges; there must be at
least two whole edges since both local forest degrees are at least
three.  Keep one whole covariance factor, of Frobenius norm `O(sqrt(n))`,
against the star's entrywise `O(n^(-1/2))` bound.  Every partial term
has Frobenius norm `O(1)` and nuclear norm `O(sqrt(n))`.  Whole matches
give the exact Schur main term, with the Hermite factorials canceled.

Raw local Wick errors and squarefree slot deletion cost `o(n)` in
nuclear norm by

    ||E[UV^T]||_* <= sqrt(E||U||2^2 E||V||2^2).

This is why an operator-sized accumulation of tiny local errors is not
silently discarded.  Gaussian and Rademacher covariance of the exact
squarefree kernels agrees, but their whole distributions are not
identified by that fact.

The primary Gaussian-chaos implication was checked directly in
[Noreddine--Nourdin, Theorems 1.1 and 1.3](https://arxiv.org/pdf/1009.1310):
fixed chaos orders, covariance convergence, and vanishing component
fourth cumulants suffice, including singular covariance limits.
Degree-one Gaussian coordinates are handled directly by the same
derivative-covariance estimate.  Different original degrees remain
separate coordinates until this finite-vector step.  Sign-input
transfer is a separate argument below, not an assumption of this
Gaussian theorem.

Read/reconstructed sources:
`fresh_tree_matrix_derivative_and_weighted_third_chain_2026_09_05.md`,
Section 1; `fresh_limit_unmarked_cubic_projection_2026_09_05.md`,
Section 2; `resumed_bound_audit_nonlinear_forest_flattening_2026_09_06.md`,
Sections 1--4; and
`resumed_bound_audit_full_nonlinear_covariance_trace_2026_09_06.md`,
Sections 1--8.

## 4. The two exceptional mixed contractions really are addressed

Put `Z=B R(X)`, first with the exact finite squarefree forest main.
An old returned coordinate has the representation
`B X_T=Q V_T+small`, `Q=B^2`, where `V_(T,a)=N_a H_a` and `H_a` is
an exact pure Gaussian child forest excluding `N_a`.  Neither `B X_T`
nor `Q[S h_T(X)]` is declared Gaussian.

First, pairing an old own-spin input with a nonlinear forest has an
unmatched output of odd degree at least three.  Either a free-free
parity edge supplies the discrepancy gain, or at least two length-two
paths survive between the distinguished endpoints.  The latter yields
`n^(-1/2)(Q^(circ p)-I)` or `B circ Q^(circ p)`, `p>=2`.  Their
Frobenius bounds give total squared old/return covariance `O(1)`.

Most further mixed contractions are controlled by the return's proper
cuts.  For contraction of **all** return slots into the higher-degree
old return there are two separate cases.

* If `N_a` is not hit, use the positive covariance-of-squares identity
  for the exact pure child chaos `H_a`, not an arbitrary sum of
  chaoses.  Its joint old-local/return Gaussian approximation bounds
  the full contraction by the squared old/return covariances plus an
  error vanishing on average, uniformly in the other root.  The
  degree-bounded Gaussian creation inequality controls the sum
  `sum_a Q_ia N_a U_(a,i)`.
* If `N_a` is hit, the expression is a diagonal of `Q J B`.  The
  doubled diagram for `E||J||F^2` has nominal size `O(n)`.  Its two
  nonlinear forest copies have disjoint leading free labels.  At
  most their common edge to `a` can cancel at the other output root;
  each has at least three distinct branches.  A parity edge survives,
  giving `O(sqrt(n))`; extra label collisions contribute `O(1)`.
  Hence `E||J||F^2=o(n)` and bounded root transport gives vanishing
  averaged diagonal error.

For Rademacher inputs, fourth-order coordinate replacement of the
unmarked test is legitimate even if the old return has high influence:
it enters only linearly.  If `rho_i` is the maximum combined influence
size of the old local fields and the forest main, the total replacement
error is `O(rho_i^2)` and its root average tends to zero.  The direct
marked test instead retains the root spin exactly.  After own-input
removal it uses

    E[S_i S_j f_i h_j]=E[(D_j f_i)(D_i h_j)]  (i!=j),

because `f_i` is own-i-free and `h_j` own-j-free.  The old derivative
is `O(n^(-1/2))`; Cauchy--Schwarz with the row norm of `Q` and the
return's total influence makes the off-diagonal sum vanish.  The
diagonal `Q_ii=1` yields `K_F`.  This checks the actual marked return,
not merely a covariance-based surrogate for it.

These are the arguments in Sections 10--12 of the full nonlinear
covariance proof and Sections 4--5 of
`resumed_bound_audit_restricted_channel_center_update_2026_09_06.md`.
They are sufficient for the two parity-restricted, linearly returned
energy tests used next.  No claim of a full conditional law for the
coherent old return is needed.

## 5. From those tests to a feasible Boolean lower certificate

For the forest main variance `u_i`, the nuclear estimate implies

    average_i sqrt(u_i) >= tau-o(1).

Indeed for any correlation matrix `Q` and odd Schur mixture `R0`,
the completely positive Schur map `X -> Q^(circ k) circ X` and
square-root Jensen give `Tr Q sqrt(R0)>=n`.  With `B^2=Q`, unit column
norms then give `sum_i sqrt((B R0 B)_ii)>=n`.  A normalized nuclear
error changes the mean standard deviation by at most
`L sqrt(||error||_*/n)`.  This is scalar square-root continuity, not
an operator-Lipschitz claim.

The required zero-variance limit is also explicit.  Add an independent
Gaussian of fixed width `delta` to the return before approximating
tests, and use a soft sign of fixed width `eta`.  All Gaussian return
variances now lie in a compact interval bounded away from zero.  The
even/odd decomposition of the test in the root spin is exactly the two
marked/unmarked identities above.  The pointwise estimate

    y softsign_eta(y+e) >= |y|-eta-2|e|

removes the smoothing uniformly even when `u_i=0`.  Convexity and
monotonicity of `r -> E|k+rN|` then use the mean-standard-deviation
bound.  For fixed polynomial `F` this proves the analytic inequality

    liminf average_i E[H(X_i)|(BF(X))_i|]
        >= E H E_N |K_F+tau N|.

The polynomial need not be a feasible spin mean.  Approximate an actual
bounded odd `F` by polynomials only as an analysis device, paying
`L||F-P||2` on the matrix side.  Projection and residual norms are
continuous on the Gaussian side.  Take the matrix limit for each
fixed polynomial, then improve that approximation.

For bounded feasible `|F|+H<=1`, `H` even and nonnegative, the actual
means `mu_+=F+H sign(BF)` and `mu_-=-F+H sign(BF)` lie in the cube.
Independent coordinate rounding is exact for a hollow matrix.  Their
quadratic energy difference is

    H_A(mu_+)-H_A(mu_-)=2sqrt(n-1) sum_i H_i |(BF)_i|.

Each expected energy lies between `-Q(A)` and `Q(A)`.  Dividing by
two proves the needed cap lower bound, with no hidden extra factor.
Thus every such Gaussian pair supplies

    liminf Q(A_n)/n^(3/2) >= J(F,H):=E H E_N|K_F+tau N|.

On bounded feasible pairs,
`|J(F,H)-J(F',H')|<=(1+E|N|)||F-F'||2+||H-H'||2`.
Finite-coordinate conditional expectation and then common Gaussian
smoothing preserve `|F|+H<=1` and parity.  Therefore countable Gaussian
objects may be approximated by finite ancestor-closed constructions.
This does not run a depth depending on matrix order.

Read/reconstructed sources: the odd-Schur proof
`fresh_limit_odd_schur_standard_deviation_2026_09_05.md`; restricted
center Sections 7--9; full nonlinear covariance Sections 13--14.

## 6. Actual anchored core and the two-Gaussian policy

The 21 explicit anchor child multisets are distinct, even, and ancestor
closed.  Their inverse features are orthonormal.  With rational anchor
coefficients `rho`, let `s^2=1-|rho|^2>0`.  The degree-200 partial-OU
resolvent deletes exactly those anchor features and normalizes the
remaining even polynomial `h` to norm one.  Its innovation derivative
energy is bounded above by

    .999904932505320245853669704564681545270483582860128726229320 < 1.

For innovations of correlation `q`, the covariance is `sum_l w_l q^l`
and `1-sum_l w_lq^l <= (sum_l l w_l)(1-q)`, including negative `q`.
Hence `Z -> U h(anchors,Z)` is a strict contraction on the complete
unit sphere orthogonal to the anchors.  Convexity of that sphere is
not required.  Its fixed point gives genuine unit Gaussian `V=U g`,
where `g` is an even polynomial of degree 200 in the independent
22-coordinate core and **has norm exactly one**.

Set `m(V)=E[g|V]`, `alpha=91/125`, and

    f(V)=1_(|V|<=alpha) P(V)-(9/250)m(V),
    P(V)=.9969843-.0728676 V^2+.5715628 V^4-1.4557984 V^6.

The displayed decimals are exact rationals.  If
`A0=<g,f>`, `nu^2=||f||2^2-A0^2`, then `nu^2>.0348840132` and

    Z=(U f-A0 V)/nu

is a standard Gaussian independent of `V` on the *same* countable
space.  Its inverse is `(f-A0 g)/nu`.  For the frozen odd ternary
rectangle policy `F1(V,Z)`, take `H1=1-F1^2`.  This is the final pure
policy (`theta=1`), so no mixture/purification gate or old causal
inverse `u` is a mathematical dependency.  Its exact inverse and
residual variance are

    K=a g+b f(V),
    a=aV-aZ A0/nu,   b=aZ/nu,
    tau^2=E F1^2-aV^2-aZ^2,
    aV=E[V F1],      aZ=E[Z F1].

The last norm uses the orthonormal inverse pair, not `||m||2` in place
of `||g||2`.  The resulting residual variance is above `.0920313540`.

The new `Z` is not independent of the 22-dimensional old core.  Its
core covariance vector `b_core` is computed from the inverse isometry:
anchor entries use `<h_j,f>`, and the innovation entry is fixed by
`A0=sum rho_j<h_j,f>+s<h_*,f>`.  It is orthogonal to the coefficient
vector of `V`, and its squared norm is below `.3027511301`.

Consequently the exact conditional polynomial `E[g|V,Z]` must be used
inside the mask.  The source does so.  Conditional projection after
partial OU smoothing has coefficients

    sqrt(binomial(d,k)) (R^2+s^2r)^(d-k)
          [s(r-1)b_core,*]^k h_(d-k)(V) h_k(Z).

The resolvent integral is evaluated by the positive finite expansion

    integral_0^1 r^(a_res-1)(R^2+s^2r)^j(1-r)^k dr
      =sum_(ell=0)^j binomial(j,ell)(R^2)^(j-ell)s^(2ell)
             k!/product_(v=0)^k(a_res+ell+v).

Anchor projections are the corresponding finite multinomial
coefficients.  These formulas reconstruct the bivariate polynomial;
the independently computed `k=0` column agrees with `m` within the
outward intervals.

## 7. The one-sided numerical certificate has no uncharged tail

On each positive-V bin `[l,r]`, the policy is -1 below its rational
lower-Z threshold, zero between thresholds, and +1 above the upper
threshold.  It is extended by joint odd reflection; for `|V|>2` it is
`sign(V)`.  The 256 bins cover `[0,2]`, are adjacent, and contain the
old alpha boundary as a bin boundary.  The source checks all these
facts and ordered thresholds in `[-8,8]`.

The code's factor two on positive-V integrals is correct because both
`g` and the mask are jointly even.  Its first moments are exactly

    E_Z F1=1-Phi(upper)-Phi(lower),
    E_Z Z F1=phi(lower)+phi(upper),
    E_Z F1^2=1-[Phi(upper)-Phi(lower)].

The `|V|>2` support and V first moment are included as `2Phi(-2)`
and `2phi(2)`; their Z first moment is zero.  Their hole is zero, so
there is no missing nonzero objective contribution there.  Gaussian
Z tails are included by the full CDF formulas, not discarded at 8.

The Jensen lower bound uses the *hole-weighted* bin mean of `K`, not
the product of the hole mass with `m(V)`.  For bin hole mass `M` and
integral `Kbin`, its term is exactly

    Kbin[2Phi(Kbin/(M tau))-1]+2M tau phi(Kbin/(M tau)).

Keeping Z-Hermite degrees `0,...,19` has the rigorous principal-angle
error `||conditional tail||2 <= ||b_core||^20 ||g||2`.  To see this,
resolve the core into `V`, the unit direction parallel to `b_core`,
and the orthogonal remainder.  The new Gaussian is an OU observation
of that direction with correlation `||b_core||`; its degree-k
Hermite coefficient is attenuated by that correlation to power k.
The conditional projection norm is at most the full norm one.
Lipschitz continuity in the shift charges the entire bin sum by
`|a| ||b_core||^20`.  This does not assume mask/inverse independence.

Bins whose mass enclosure is too small are omitted only from a
nonnegative objective sum.  All their support and first-chaos moments
were included before that decision, and the global tail charge remains.

The rational interval implementation was read: `I` rounds every
operation outward to a `10^-60` grid; square roots use integer square
roots; pi is enclosed by Machin's identity and alternating arctangent
tails.  `phi_large` computes `exp(-x^2/64)` with a finite explicit
remainder and squares five times.  The CDF uses 129 integrated Taylor
terms: for `|x|<=8` the omitted alternating tail is decreasing and
bounded by `8*32^129/(129!*259)<10^-24`, even though early terms need
not decrease.  Interval arguments crossing zero are enclosed using
the density bound `phi(0)`.  All CDF argument ranges are asserted.
The degree-200 factorial normalizations are formed as exact rational
ratios before conversion when necessary, so no tiny reciprocal is
silently rounded to zero.  Central Gaussian power moments use a
separate decreasing alternating integrated series at `alpha<1`.

Read source modules, not only JSON: the birth certificate;
`resumed_response_full_center_certificate_2026_09_06.py`;
`resumed_response_conditional_v_certificate_2026_09_06.py`;
`resumed_response_feedforward_inverse_certificate_2026_09_06.py`;
`fresh_finite_anchor_fixed_point_certificate.py`;
`fresh_limit_rooted_lower_certificate.py`; and the `cdf_point`,
`phi_large` helpers in `fresh_limit_mask_ascent_certificate.py`.

Fresh replay command:

    .venv/bin/python computations/resumed_response_rich_core_birth_certificate_2026_09_06.py --policy computations/results/resumed_response_rich_core_optimized_rectangle_policy_2026_09_06.json --theta 1 --target 4333/10000 --output /home/math/quadra/tmp/decisive_audit_rich_lower_replay_2026_09_07.json

It reproduced

    pre-tail lower = .433325375702943007494944152418094290701027670940273741230068
    tail upper     = .000003264038862254079131223520514944142393022906580328123072
    final lower    = .433322111664080753415812928897579346558634648033693413106996.

The replay and canonical JSON have identical SHA256
`c6ab5199c366952bbfc831d26d1527149f2adaeb7120c4c5ec9b1d3ee7f86dcb`.
This is a fresh replay of the checked implementation, not a claim that
a second independently written full interval engine was used.

## 8. Limit order and deliberately excluded dependencies

For fixed deletion fraction, hence fixed `L`, first choose any desired
finite approximation to the countable Gaussian feasible pair.  For
each analytic polynomial approximation and fixed positive smoothing
widths, send the retained matrix order to infinity.  Improve those
polynomial approximations, remove the smoothing widths, and improve
the finite-coordinate approximation in their established L2 norms.
The resulting Gaussian value is independent of `L`.  Finally send
the deletion fraction to zero.  Equivalently each desired numerical
accuracy can be met by a finite list of fixed choices before its
matrix limit.  No degree, number of anchors, root variance cutoff,
or iteration depth grows during that limit.

The .433322 claim does not depend on the floating ascent being globally
optimal, the old partial-mixture optimum, an operator-norm estimate
for nonlinear raw covariance, Gaussianity of the coherent first-chaos
return, moment determinacy of arbitrary polynomials in a rich frame,
or the later high-degree full-feedback extensions.  The remaining
original convergence question is unaffected by this audit.
