# Fourth checkpoint: exact feasibility, dynamic escape, and nonlocal tests

Date: 2026-09-06, approximately 04:10 UTC. The six-hour campaign continues.
The original lower endpoint remains .4314603928237005; convergence is open.

## 1. Feasibility is now a theorem, not an assumed Gaussian closure

I derived the matrix-height recovery theorem and two researchers separately
reconstructed it. For normalized even vector features g, the Gaussian
correlation kernel K(P)=E[g(X)g(Y)^T] is continuous and Loewner monotone.
An actual canonical frame G=Ug(G) exists exactly when K^h(0) tends to I,
equivalently when there is no proper subfixed matrix. Compatible height
projections construct the unique solution and give its exact L2 tail.
The anchored form allows actual cyclic finite anchors.

The scalar critical case was already in yesterday's archive. The new
claim is the matrix/measurable/anchored recovery theorem, not rediscovery
of that scalar case. See
[the proof](resumed_director_creation_recovery_theorem_2026_09_06.md),
[first reconstruction](resumed_bound_audit_creation_height_converse_2026_09_06.md),
and [second reconstruction](resumed_response_matrix_recovery_independent_audit_2026_09_06.md).

A two-component rough feedback loop has a proper subfixed matrix for
every positive return strength r. I strengthened this to approximate
unrealizability: directional Minkowski traps the actual height covariance
below Q whenever both the kernel and residual fit that barrier. The
specific loop has an operator-L2 residual floor r^2/(40sqrt(2)). This is
not a lower bound on the original signing energy. Critical barriers need
not have positive approximation floors. The independent auditor checked
the inequality and constants in
[the quantitative theorem](resumed_director_quantitative_feedback_barrier_2026_09_06.md).

The finite-energy necessary condition in
[the feedback note](resumed_response_matrix_height_feedback_obstruction_2026_09_06.md)
is also reconstructed: diagonalize an input PSD direction, increase the
Gaussian covariance deficit quotient to its extended Dirichlet form,
and use compactness after removing its common zero kernel. If every
relevant output combination were strictly expansive, a finite covariance
deficit would already be a subfixed matrix. Equality remains admissible;
roughness in a feedforward direction is not by itself an obstruction.

## 2. Uniform improvement belongs to a fixed frame, not all depths

I read and reconstructed the continuous-core proof, the fractional-tie
escape, and the uniform two-stage argument. A continuous cyclic core plus
arbitrary finite causal appendages has no ternary stationary response
above phi(0). After causal coefficients are eliminated, the stationary
binary support is locally constant off two linear level planes. The
resulting scalar equation forces its continuous inverse response to take
only finitely many values, contradicting a nonzero gradient.

Fractional ties require a different argument: an independent even gate
purifies the support while preserving the conditional pair and the exact
objective. The gradient gains variance E[q(1-q)s(K)^2]. An arbitrarily
small flat purification followed by line search gives a quadratic
improvement; it does not make the original first derivative positive.

For a fixed continuous-core frame, full purification gives the exact
variance identity

    rho^2=E[(1-q)s(K)^2]-4B E[(1-q)K s(K)]
                            +4B^2|a|^2-|A|^2.

The q^2 terms cancel, making the post-purification gap continuous on the
weak-* compact high-value feasible set. It is everywhere positive, hence
has a positive minimum. This yields a uniform two-stage gain for that
fixed frame. The gate inverse must be orthogonal to both the old inverse
span and the old gradient. The frame changes after a step; no uniform
rate over successive frames follows.
[Proof](resumed_response_continuous_core_and_fractional_tie_escape_2026_09_06.md),
[independent audit](resumed_bound_audit_continuous_core_escape_2026_09_06.md).

## 3. Nonlocal arithmetic constructions: the exact licenses

I checked the primary Gauss-sum mapping, primitive and lower-level phases,
and the constant-channel correction. The order-49 exact replay uses
G_49=512(-1-i sqrt(7)chi) and G_7=832-448i sqrt(7)chi. An initially
incorrect factor eight in the latter was repaired before publication.
The director replay passed. The signed resonant carrier restores the
full seven-point rotation, including constants, in the limiting Boolean
test. Its level and field-extension limits are ordered explicitly.

Tensoring independent rotations and pulling back along a linear map gives
independent phases on projective quadratic-residue triples. I separately
derived C C^T=7^(r-1)I and the polar squared-defect identity. The phase
theorem does not give independently chosen phases at all frequencies or
arbitrary overlapping gates. Its exact remaining Boolean feasibility
condition is retained in
[the phase theorem](resumed_convergence_programmable_projective_phases_2026_09_06.md).

The vectorial Maiorana--McFarland map gives a genuinely different operator
F_q P_p F_q on balanced profiles. I reconstructed its complete finite
Fourier identity, including the exceptional delta and trace distribution,
and reran the independent integer verifier at q=8,16. No arbitrary
permutation or restoration of unbalanced constants is claimed.
[Proof](resumed_convergence_vectorial_mm_power_operators_2026_09_06.md).

My exact C++ enumeration of all balanced pairs at q=16, p=7 gives maximum
1792/512=7/2 for the weighted seed [-1,2;2,4], with 82,818,450 pairs after
simultaneous reversal. This is a complete finite test, not an upper bound
on the whole regularized norm. The program checks its finite-field trace,
orthogonality, and constant preservation in integer arithmetic.

Finally, the conditional-bent reduction is exact: if every large-block
conditional section is bent, its Walsh image is Boolean and all remaining
seed amplification is confined to the selector block. I checked the
pointwise Boolean upper reduction and common-bent-carrier equality. The
large bent block cannot amplify this restricted construction's norm.
[Proof](resumed_convergence_conditional_bent_reduction_2026_09_06.md).

## 4. Director judgment and next work

The response feasibility gap is materially narrower and has an exact
construction/obstruction criterion. The general high-value rough cyclic
stationary case and a rate over expanding frames remain unresolved. The
fixed-rule certificate ceiling sqrt(15)/8 still applies; it is not a
statement about actual signings.

The nonlocal track has new realizable tests but no upper-preserving map
from near minimizers to arbitrary orders. The original convergence gap
has not been removed. Continue with a general nonlocal synthesis or
falsifier and independent optimizer-specific analysis, rather than
tuning the weighted seed or adding further smoothness variants.
