# Edge-transitive signed symmetry gives an exact radial ground-state law

Date: 2026-09-06. Elementary finite theorem, independently derived by the
director from the Paley radial identity. This is an explanation of that
identity's sufficient symmetry, not a claim about arbitrary minimizers.

Let n>=2, let A be a hollow symmetric order-n signing, and let
Q(A)=max_{x Boolean}|x^T A x|/2. Suppose a finite group consists of pairs
(P_g,epsilon_g), where P_g is a signed permutation matrix and
epsilon_g is a sign, with

    P_g^T A P_g = epsilon_g A.

Suppose its underlying permutation action is transitive on unordered pairs
of distinct coordinates. No transitivity on vertices or ground states is
additionally needed. The map epsilon is automatically a character.

## Exact theorem

For any oriented ground state (x,sigma), so sigma*x^T A x/2=Q(A), average

    X=P_g x,   Sigma=epsilon_g sigma

over the group. Every sample is again an oriented ground state. Its
off-diagonal covariance is exactly

    E[Sigma X_i X_j] = [2Q(A)/(n(n-1))] A_ij,   i != j.       (1)

The diagonal is not asserted zero: it equals sigma E epsilon_g. It is zero
if the group contains an element of negative epsilon, but hollowness makes
it irrelevant to all conclusions below.

Proof. Write M=E[Sigma XX^T]. Haar averaging on the finite group gives
P_h M P_h^T=epsilon_h M. The same identity holds for A. Therefore M_ij/A_ij
is invariant under the underlying permutation action on unordered pairs:
the two switching signs and the orientation character cancel. Pair
transitivity makes it a constant rho. Every sample is an oriented ground,
so summing A_ij M_ij over i<j gives Q(A)=rho*n(n-1)/2. This proves (1).

For any real symmetric hollow B, averaging its oriented energy gives

    Q(B) >= [2Q(A)/(n(n-1))] |sum_(i<j) A_ij B_ij|.          (2)

In particular, if B differs from A on d unordered signing entries,

    Q(B) >= Q(A) |1-4d/(n(n-1))|.                           (3)

These statements use exact ground states but no spectral hypothesis.
The distribution need not be computably available to prove its existence.

## Consequences and limitations

Cap values at a fixed order all have the parity of binom(n,2). Thus fewer
than n(n-1)/(2Q(A)) edge flips cannot lower Q(A). If a sequence in this
symmetry class has Q(A_n)/n^(3/2)->c, every signing at coefficient at most
c-delta must stay at least (delta/(4c)-o(1))*n^2 edge edits from A_n and
from -A_n. This last conclusion assumes c>0 and a fixed positive delta.

For Paley cores the affine group supplies exactly this signed action, with
the orientation character given by the quadratic character of the scale.
Its asymptotic c=1/2 is a separate theorem. The resulting edit lower bound
is (delta/2-o(1))*n^2, matching the explicitly calculated Paley law.

This principle does not prove that exact minimizing signings possess such
a symmetry, nor that their oriented ground-state covariance has a radial
representative. Averaging over all switchings/relabelings of an arbitrary
minimizer changes the matrix as well as the spin. It does NOT manufacture
(1) for one fixed matrix. In particular the new law cannot be substituted
for an optimizer-specific radial law without a separate proof.
