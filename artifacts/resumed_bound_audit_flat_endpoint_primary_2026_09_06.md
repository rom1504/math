# Independent primary/dependency audit: flat-involution terminal endpoint

Date: 2026-09-06. Audit checkpoint: 06:17 UTC.

Outcome: PASS in the explicitly fixed-rule, operator-close Hadamard
scope of `resumed_response_flat_involution_terminal_endpoint_2026_09_06.md`.
The two requested dependency checks close. The fresh-noise argument in
the first-order module should use the noncircular primary justification
in Section 3 below. No convergence theorem for the original optimum is
claimed here.

## 1. Haar old-tree moments: independent reconstruction

Let U be a balanced Haar orthogonal involution. Its normalized tracial
moments converge in L2 exactly to the balanced two-point law. The
orthogonally invariant injective-traffic theorem applies to every fixed
connected quotient multigraph: its limiting weight is the product of
cycle free cumulants for a bridgeless cactus, and zero otherwise.
This is Theorem 4.2 of
[Gorini--Jones--Kunisky--Pesenti](https://arxiv.org/pdf/2604.11729),
not the separate deterministic punctured theorem. The primary source and
the local source `tmp/paper_2604_11729/traffic-dist.tex` were checked.

For a moment of the canonical injective old trees at one common root,
write D for the total number of edge occurrences, also the total number
of marked nonroot positions. The common root is excluded from every
marked position. Nonzero centered input moments force marked blocks of
size at least two, so the quotient has v <= D/2+1 vertices. Individual
tree injectivity prevents loops. A loopless cactus with e edges has
e <= 2(v-1), with equality only for a doubled tree. Since e=D, every
survivor has equality throughout: all marked blocks have size two and
the quotient is a doubled tree. In particular no fourth or higher free
cumulant is present.

At a marked vertex exactly two tree occurrences meet. Each incident
edge must be doubled by the other occurrence, not by a second vertex of
the same injective copy. Propagation through the connected nonroot part
therefore pairs whole isomorphic rooted copies. Rooted-automorphism
normalization gives weight one per pair. These are exactly the Wick
moments of independent standard Gaussian tree coordinates. The argument
works for Gaussian or Rademacher inputs; blocks of size at least four
lose a quotient vertex and vanish. It also gives all fixed polynomial
moment bounds needed for cutoff approximation.

Using two independent colored input replicas gives the claimed averaged
conditional variance/covariance concentration. In particular the
fourth-cumulant statement must use both the annealed fourth moment and
the colored conditional-second-moment calculation, as the traffic note
does. The local averaged-root statement does not assert uniform-root
or operator-Gram universality.

## 2. Haar paired energy and the finite recursion

An energy term F(X_i) U_ij S_j H(X_j) adds one edge and one marked
position. Thus e=D+1 and v <= (D+1)/2+1. If root i is identified with an
H-side marked position, even this upper bound loses a vertex and no
leading doubled-tree term remains. Separately, the diagonal bridge is
o(1) by the vanishing maximum diagonal and the polynomial moment bounds.

In a surviving off-diagonal diagram, i is unmarked and the bridge must
double a top edge of one F-side tree occurrence. The explicit S_j pairs
with its child mark. Removing that occurrence exposes precisely its
normalized child-Hermite product on the H side. Its branches cannot
self-pair, because labels within that selected tree were injective.
Other F-side occurrences pair in whole copies. This gives the claimed
derivative/child-product paired identity, with its stated sign and no
extra factor. This reconstructs Sections 2--3 of
`resumed_convergence_haar_response_traffic_2026_09_06.md`.

Consequently, on both the original Hadamards and the Haar comparator,
the three normalized expectations

    ||X_T||^2/n,  ||B[S h_T(X)]||^2/n,
    <X_T,B[S h_T(X)]>/n

converge to 1. The second uses B^2=I+o_op(1). Their difference proves the
normalized-L2 tree recursion, not merely a formal tree expansion.
Cut off each even child polynomial first, then choose the already
constructed child approximations at the resulting Lipschitz tolerance.
The bounded update clip(S) h_cut(children) is globally Lipschitz in all
formal inputs. Induction constructs an actual fixed finite equivariant
first-order approximation simultaneously on both ensembles. The order
of cutoff and child approximation is essential and is correct in the
endpoint proof.

## 3. Primary state evolution and a noncircular nonsingularity repair

The relevant primary result is
[Wang--Zhong--Fan](https://arxiv.org/pdf/2206.13037), Proposition
2.7(b)(2) and Theorem 2.8. The local source
`tmp/paper_audit_2206/universalityAMP.tex` was checked directly.
Random signed-permutation conjugation of an exact incoherent involution
satisfies the proposition: even powers are I, odd powers are U, and
the diagonal trace error is at most twice the maximum entry. Its
spectral law is the balanced two-point law. The theorem gives almost
sure empirical Wasserstein-2 convergence to a deterministic Gaussian
state-evolution law, subject to nonsingular prescribed covariances.
It is substantially stronger than agreement of empirical moments.

Here is a precise way to verify that last hypothesis for an arbitrary
fixed Lipschitz first-order computation. Smooth the finite coordinate
maps first, preserving the signed parity by symmetrization. This has
arbitrarily small uniform map error and hence arbitrarily small
normalized-L2 recursion error at bounded operator norm. Serialize the
computation, and add a unique fresh Gaussian side-information variable
delta eta_t to the input of multiplication t. Earlier updates do not
use eta_t. Include an independent noisy initialization if necessary.

Use [Fan](https://arxiv.org/pdf/2008.11892), Assumption 4.2(d--e) and
Theorem 4.3. Smooth Lipschitz maps meet its regularity hypothesis. Its
nondegeneracy condition excludes an update that is a deterministic
linear combination of preceding Gaussian states and updates. The fresh
eta_t forces this exclusion: condition on all the other side information
and state variables, which leave a nonzero independent Gaussian term on
the new-update side only. The theorem establishes nonsingularity
inductively, so no next-step state evolution is assumed to prove its own
hypothesis. Its coefficient limits are the deterministic prescriptions
used by WZF; finite Lipschitz stability removes empirical-coefficient
errors if that formulation of Fan is used.

For an independent algebra check, Fan's equations (4.4)--(4.7) use
strictly lower triangular Phi and

    Sigma_t = sum_(j>=0) kappa_(j+2)
                    sum_(i=0)^j Phi^i Delta (Phi^(j-i))^T.

Adding fresh variance delta^2 to Delta_tt adds exactly
kappa_2 delta^2 e_t e_t^T, because Phi e_t=0. Here kappa_2=1 and
b_tt=kappa_1=0. This checks the innovation sign and normalization.
By itself that algebra would still need a PSD justification for the
unperturbed candidate; the preceding primary nondegeneracy argument
supplies the full induction without this gap.

The first-order-to-AMP change of coordinates is triangular: after
defining u_t from past coordinates, its prescribed coefficient row is
determined; then y_t=z_t+sum_(s<=t)b_ts u_s recovers the original
matrix-product register. There is no dependence on a not-yet-defined
update. All maps remain fixed and Lipschitz at every fixed delta.

## 4. Why the self-energy limit is deterministic and zero

Expose energy by one further multiplication. Empirical W2 convergence
controls the continuous quadratic coordinate test m_i (Bm)_i. Thus
each fixed smoothed/noisy rule has a deterministic energy limit, common
to U and -U, since their limiting diagonal distributions are identical.
Signed-permutation equivariance removes the random conjugation exactly
in law: the transformed jointly centrally symmetric iid seed rows
retain their original joint distribution.

Now remove the fresh noise and smoothing, keeping the rule fixed before
taking n to infinity. Finite Lipschitz stability bounds normalized
mean-square output error by a constant times the squared perturbation.
For bounded outputs and bounded operator norm this also bounds energy
error. The deterministic perturbed energy limits are therefore Cauchy;
the original rule has the same deterministic limiting property. This
step must precede the matrix-sign-charge argument, since the auxiliary
fresh noise need not preserve that charge.

For the actual charge-preserving tree approximants, X_T(-B)=-X_T(B),
F(-B)=-F(B), and the soft center C(-B)=C(B). Hence each self-energy
changes sign under B -> -B, exactly. Its common deterministic limit is
its own negative, and is zero. Feasibility bounds these normalized
energies uniformly, upgrading convergence in probability to L2.
The same reasoning applies to an ordered normalized-L2 limit of the
fixed approximants; it does not require uniformity over all rules.

The endpoint's hard-sign passage uses its separate flat-involution
local-channel theorem: nonzero bounded odd F has strictly positive
higher-chaos variance, giving the stated small-ball estimate. This
allows the soft center to approach the hard center in normalized L2.
The cross-energy limit is identified by that channel, not by the Haar
traffic argument alone. The dependency distinction is maintained.

## 5. Scope and verdict

The requested external dependencies and deterministic self-energy step
pass, with Section 3 above making the fresh-noise argument explicit.
All limits are fixed family/rule first, n second, approximation errors
last. There is no conclusion for growing depth, matrix-adapted seeds,
arbitrary non-equivariant rules, or a comparator obtained only from
Frobenius-near involution structure. The terminal-functional ceiling
therefore concerns the stated architecture on Hadamard sequences; it
is neither an upper bound on the Boolean optimum nor a proof of
convergence of M_n/n^(3/2).
