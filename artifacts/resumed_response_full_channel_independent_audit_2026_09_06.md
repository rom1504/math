# Independent reconstruction of the full nonlinear center theorem

Date: 2026-09-06. Response-track audit of
`resumed_bound_audit_full_nonlinear_covariance_trace_2026_09_06.md`,
the original forest-flattening module, the exact injective-input Gram
module, and `resumed_full_forest_root_hit_director_2026_09_06.md`.

Status: the finite-polynomial modules and the bounded feasible response
handoff below pass an independent reconstruction. This is not an assertion
of a full conditional law for BF; only the two parity-tested identities
needed for the center comparison are proved. Numerical evaluation is a
separate step.

## 1. Precise conclusion reconstructed

On the countable old Gaussian tree space, let F be jointly odd and H
jointly even, with 0<=H<=1 and |F|+H<=1. Write

    K=U*F,  R=F-P_1F,  tau=||R||_2.

For finite Gaussian-a.e.-continuous responses, the full center theorem is

    liminf M_n/n^(3/2) >= E[H E_N |K+tau N|],                 (1)

where N is an independent standard Gaussian. Section 7 extends it to
arbitrary measurable responses on the countable old space. It remains
an original-problem lower bound after the same principal-deletion
argument used for the restricted theorem.

The entire norm of the higher LOCAL Gaussian chaos is used. Local
chaos degree and original input-spin chaos degree are different; the
proof always separates them at the required steps.

## 2. Covariance: why the nuclear-norm conclusion is sufficient and valid

For a fixed odd polynomial R with no first local chaos, local Wick
replacement gives a finite sum P of homogeneous old-branch forests,
with uniform local L2 error o(1). Every forest has at least three old
branches. The global root maps of the unsymmetrized branch products
are bounded by the Schur Gram inequality, and symmetrization is an
orthogonal projection, so their covariance operators are bounded.

In a covariance slot pairing, form the bipartite branch-matching graph.
The following exhaustive classification was checked directly:

* A one-edge component is a whole-branch covariance.
* A nontrivial star has an odd number p>=3 of leaves. Contracting one
  leaf into its split center yields O(n^-1/2) by a proper flattening.
* A connected component which is neither an edge nor a star contains
  a P4. Contract its disjoint outer edges first. The two inner vertices
  each have another neighbor, so both cuts are proper. The two gains
  multiply to O(n^-1), including when the graph has extra edges.

Contraction of all remaining shared slots cannot increase the product
of the Hilbert norms. This follows by repeatedly merging two tensors,
contracting ALL indices they share at that merge. Thus cycles in the
branch graph do not invalidate the bound.

A partial pairing with a nonstar component or two star components has
entry bound O(n^-1), hence Frobenius norm O(1). The only remaining case
is one star plus whole edges. There are at least two whole edges,
because the side containing the star center has local degree at least
three. One whole-branch Gram factor has Frobenius norm O(sqrt(n));
the star contributes O(n^-1/2) entrywise. This again gives O(1)
Frobenius norm. Hence all partial terms have nuclear norm O(sqrt(n)).

Whole-branch pairings give sum_k ||R_k||_2^2 Q^(circ k). In replacing
the exact old covariances by delta_TU Q, the other Schur factors have
bounded Schur multiplier norms because they have Gram factorizations
with bounded row vectors. A mere entrywise bound would not suffice.
The old operator covariance theorem therefore gives an o_op(1) error
in each finite whole-branch product. Normalization of repeated branches
cancels exactly against the local normalized Hermite multiplicities.

Local Wick errors and marked-slot collision deletion change covariance
by o(n) in nuclear norm using

    ||E[V W^T]||_* <= sqrt(E||V||^2 E||W||^2).

The collision-free homogeneous kernels have identical Gaussian and
Rademacher covariances. Direct local sign Wick replacement finishes
the transfer. Thus for actual raw R(X),

    ||Cov(R(X))-sum_k ||R_k||_2^2 Q^(circ k)||_*/n -> 0.

This is deliberately not promoted to operator-norm convergence.

For u_i=E[(BR(X))_i^2], the scalar square-root inequality gives an
o(1) error in their mean standard deviation after transport. Indeed
if E_n is the covariance error,

    avg_i |sqrt((B C B)_ii)-sqrt((B S B)_ii)|
       <= sqrt(Tr[Q |E_n|]/n)
       <= L sqrt(||E_n||_*/n).

The exact odd-Schur mean-standard-deviation theorem applies to the
convex mixture S/tau^2 and yields avg_i sqrt(u_i)>=tau-o(1).

## 3. The endpoint graph really eliminates old drift

For an exact own-spin input tree V_T,a and a fully injective forest
P_j with k>=3 branches, a nonzero covariance identifies the marked
sets bijectively. Root a is among those marked vertices but root j is
not. The prefactor for fixed a,j is n^-1/2. The parity graph has odd
vertices a,j, and its degree at j is EXACTLY k: the input tree has no
vertex j and cannot cancel any edge incident there.

If a free-free edge survives, bilinear discrepancy contributes another
n^-1/2 and the covariance matrix has Frobenius norm O(1). Otherwise
every nonisolated free vertex must join both endpoints, since its
parity degree is even. The graph is an optional a-j edge plus p length-
two paths, with p>=2. Its matrix is n^-1/2 Q^(circ p) or
B circ Q^(circ p), up to Frobenius-o(1) error. The correlation Schur
bound controls the first, and the row sum sum_j |Q_aj|^p<=L^2 controls
the second. Both have Frobenius norm O(1).

The old output collision has covariance operator O(1/n); the fully
injective forest has bounded covariance operator because collision
deletion is a common tensor-space projection. Covariance Cauchy--Schwarz
therefore preserves Frobenius O(1) when passing from V_T to X_T and
transporting through B. In particular,

    ||Cov(BP,X_T)||_F=O(1).

This is the all-root squared-covariance control needed below, not just
vanishing same-root covariance. Raw local errors are harmless at the
averaged scale where they are actually used.

## 4. Exceptional first-chaos transport: both full-contraction classes

The nonlinear forest vector is jointly Gaussianized by its small
proper flattenings; B X_T need not be Gaussian. All its mixed
contractions with a forest Z therefore have to be checked. Proper-Z
contractions follow from its small self-contractions. For a full Z
contraction use B X_T=Q[N H]+small error and split according to whether
the contraction hits the explicit root N_a.

Root not hit: the term is sum_a Q_ia N_a U_a,i. The U_a,i exclude N_a,
and the fixed-degree creation inequality bounds the squared norm by
C sum_a Q_ia^2 ||U_a,i||_2^2. The needed estimate follows in this order:

1. Treat H_a as its exact homogeneous child-forest input chaos, or
   replace a local Wick polynomial by that pure TOP input chaos first.
2. Apply the positive covariance-of-squares contraction identity to
   bound the full contraction by Cov(H_a^2,Z_i^2).
3. In the joint Gaussian approximation of the local old family and Z,
   Gaussian integration by parts expresses that covariance as a
   quadratic form in Cov(X_T,a,Z_i). Its coefficients are fixed
   polynomial Gaussian moments and uniformly bounded.

This gives C sum_T |Cov(X_T,a,Z_i)|^2+epsilon_i, with average epsilon_i
tending to zero uniformly in a. The error form follows from uniform
old flattenings and averaged Z collision/flattening errors. The
all-root Frobenius bound from Section 3 and bounded Q row norms now
make the averaged root-not-hit norm vanish. Applying the positivity
identity directly to a general sum of input chaoses would be invalid;
the first step above is essential.

Root hit: after fixing a finite slot-matching pattern the term is
diag(Q J B). Let the forest have D marked vertices and the even child
forest have d. In E||J||_F^2 the leading diagram has

    2+2(D-1)+(d-D+1)=D+d+1

free labels and normalization n^(-(D+d)), so absolute counting gives
O(n). Any extra identification costs a whole label, including a
remaining child slot equal to the external forest root j. In a leading
pattern the two contracted-label sets are disjoint. At root j the two
forest copies can share at most their edge to a; because each has at
least three incident edges, a parity edge survives. All labels are
summed, so bilinear discrepancy gives n^-1/2 and hence

    E||J||_F^2=O(sqrt(n)),
    avg_i E|[Q J B]_ii|^2=O(n^-1/2).

This reconstruction agrees with the director proof. It does not need
a general high-order random derivative-matrix operator theorem.

All mixed contractions thus vanish at the needed averaged scale.
The product formula gives the needed moment independence of Z from
the exceptional first-chaos block, jointly with the old local family.
That block may itself retain coherent nongeneric Q-spin atoms.

## 5. Rademacher input and the two restricted identities

Marked-slot projection preserves bounded row Hilbert norm and total
influence. Its averaged Hilbert error makes the average squared maximum
influence tend to zero. If c_i,a controls the first derivative of old
and forest inputs, put rho_i=max_a c_i,a. For the unmarked test
Y_i f_i with squarefree Y_i=(B X_T)_i entering linearly, fourth-order
Lindeberg has only Y_i D_a^4 f_i and 4 D_aY_i D_a^3 f_i. Its total
error is bounded by

    C [sum_a c_i,a^4 + (sum_a c_i,a^6)^(1/2)] <= C rho_i^2.

Uniform fixed-degree hybrid moments and bounded total influence were
checked in this estimate. The averaged error vanishes even if Y_i has
large individual influences. The direct Boolean marked test uses
own-root removal and a two-spin derivative; its off-diagonal channel
term costs at most

    (n sqrt(n))^(-1) sum_ia |Q_ia| ||D_a Z_i||_2=O(n^-1/2).

Thus for Z=B R(X), finite odd polynomial F, and even old polynomial M,
the required parity tests have Gaussian model

    E[M psi_odd(Z) BF] = EM E[Z psi_odd(Z)],
    E[S M psi_even(Z) BF] = E[MK] E[psi_even(Z)],

with Z modeled by N(0,u_i), independent of old X. Statements are
averaged with arbitrary bounded deterministic root weights. These are
not assertions about the full conditional law of BF.

## 6. Bounded feasible F: use polynomial channels only as test devices

Fix L and bounded feasible F,H on a finite old family. Approximate F
in Gaussian L2 by a fixed odd polynomial P, and set

    K_P=U*P,  R_P=P-P_1P,  Z_P=B R_P(X),  tau_P=||R_P||_2.

The actual feasible center test is

    C_i=H(X_i) tau_eta(S_i K_P(X_i)+Z_P,i+delta xi_i),

where tau_eta is the clipped linear soft sign and xi_i is independent
standard Gaussian. It satisfies |C|<=H regardless of how large P is.
The endpoint means remain the bounded ACTUAL ±F+C, not ±P+C.
Replacing BF by BP in the tested cross term costs at most

    L ||F-P||_2+o(1)

by averaged L2 transport and |C|<=1.

For each fixed P, adding delta xi bounds channel variance below by
delta^2. The common dominating Gaussian density then justifies the
bounded joint-test extension used in the restricted theorem. At root
variance u_i, the smoothed tested regression is

    S K_P + [u_i/(u_i+delta^2)](Z_P+delta xi).

Here the Gaussian variances used for the common-density approximation
are those of the bounded-row FOREST MAIN kernels. The raw channel
B R_P is only known to be close to that main vector in averaged L2;
its actual variance at every individual root need not have a uniform
cap, since normalized nuclear covariance control is not operator
control. The forest main variances are uniformly bounded by their
row Hilbert norms. For fixed eta the soft sign is Lipschitz, so the
raw/main averaged L2 error changes the bounded cross test by o(1).
The same error changes the mean standard deviation by o(1). Thus the
mean-standard-deviation lower bound still applies to the main variances
used here. No uniform root-variance claim for raw polynomial transport
is needed.

Equivalently this is the conditional expectation of S K_P+Z_P given
the smoothed channel, so the elementary inequality

    y tau_eta(y+e) >= |y|-eta-2|e|

gives a lower bound

    avg_i E H |K_P+sqrt(u_i)N| - eta-2delta E|N|.

This argument is valid even when some u_i vanish. It never approximates
sign(sqrt(u)N) uniformly in polynomial L2 as u decreases to zero.

The mean-standard-deviation theorem and convexity/monotonicity in the
Gaussian standard deviation give the bound E H|K_P+tau_P N|. Now
||K_P-K||_2<=||P-F||_2 and |tau_P-tau|<=||P-F||_2. First pass to the
matrix limit for fixed P,delta,eta and fixed L; then improve P and
remove eta,delta. This proves (1) at fixed L. The right side is
independent of L, so principal deletion removes the operator cap.

## 7. Arbitrary countable-field responses and precise limits

Condition F,H jointly on finite ancestor-closed old sigma-fields and
then apply the same Gaussian OU smoothing to both. Conditional Jensen
and positivity preserve |F_j|+H_j<=1, all bounds, and parity. The
resulting finite smooth responses converge to F,H in L2. Moreover

    ||K_j-K||_2<=||F_j-F||_2,
    |tau_j-tau|<=||F_j-F||_2.

The functional is continuous under this joint approximation:

    |E H_j E|K_j+tau_j N|-E H E|K+tau N||
      <= ||H_j-H||_2+(1+E|N|)||F_j-F||_2,

because ||K_j||_2^2+tau_j^2=||F_j||_2^2<=1. Apply the finite theorem
for each fixed j and then let j tend to infinity. This extends (1)
without compactness, a maximizing mask, or an iterative closure claim.

The full theorem is stronger than the finite edge-Hermite-channel
theorem, but remains a one-step variational mechanism. The exact
iteration/twin obstructions in the separate response-track artifact
remain applicable; the full norm does not automatically create fresh
independent center spins after the update.
