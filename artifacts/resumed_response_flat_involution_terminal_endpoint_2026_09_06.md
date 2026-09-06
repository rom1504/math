# Terminal self-energies vanish on flat involutions

Date: 2026-09-06. Status: the full local-channel specialization was
independently reconstructed by the bound-audit agent. The separate GFOM
parity argument and tree-recursion embedding were independently checked
by the convergence agent and the director. This note has a deliberately
restricted matrix class.

## 1. Matrix and algorithm scope

The primary original-signing application is a sequence of normalized
symmetric Hadamards U_n. They are deterministic exact involutions with

    U_n^2=I, max_(i,j)|U_(n,ij)|=n^(-1/2+o(1)).    (1)

We allow the operator-o(1) perturbations B_n obtained by hollowing and
normalizing such a sequence. More generally the perturbation conclusion applies
when an explicit operator-close comparator satisfying (1) is supplied
and the rules have the stated normalized-L2 stability AND the old-tree
Gaussian/marked identities used in Section 2. For arbitrary nonconstant
exact involutions satisfying (1), the same Lipschitz GFOM approximants
have the universal conclusion, but identifying arbitrary injective-tree
formulas with those approximants requires these old-tree identities
separately. They are verified here for the original Hadamard signings
and the Haar comparator. We do NOT infer
such a comparator from Frobenius-near involution structure, nor from
||B_n^2-I||op=o(1) alone.

Use finitely many of the canonical odd input-tree fields X_T, with their
child closure, and iid Rademacher seeds S. Let F be a fixed bounded odd
Gaussian-a.e.-continuous function of these fields and H a fixed bounded
even Gaussian-a.e.-continuous function, with H>=0 and |F|+H<=1 pointwise. Assume
||F||_(Gaussian L2)>0. Define

    C_n=H(X_n) sign(B_n F(X_n)), sign(0)=0,
    Q_B(v)=v^T Bv/2.

Then the conclusions below are

    Q_B(F)/n ->0 and Q_B(C)/n ->0 in L2,            (2)
    F^T B C/n ->J(F,H) in L2,                     (3)
    Q_B(C+F)/n -> J(F,H),
    Q_B(C-F)/n ->-J(F,H), both in L2,              (4)

where J is the full terminal Gaussian functional. Thus the asymptotic
expectation AND variance of the common normalized energy
[Q_B(F)+Q_B(C)]/n are zero. No quantitative finite-n rate is asserted.

The same conclusion applies to fixed finite approximations of the
actual rich cyclic core, and then to its Gaussian variational closure
with the usual ordered limits. It does not claim an infinite algorithm
or an arbitrary matrix-dependent signing satisfies (2)--(4).

## 2. An explicit GFOM approximation of every finite input-tree family

The applicability gap for injective tree formulas can be closed here
without claiming that graph polynomials are automatically first-order
algorithms. For a tree T put

    Y_T=B_n[S h_T(X)],

where h_T is its normalized even child-Hermite product. The old local
Gaussian law gives

    E||X_T||^2/n ->1, E||S h_T(X)||^2/n ->1.

Because B_n^2=I+o_op(1), also E||Y_T||^2/n ->1. The independently proved
marked hierarchical energy identity, with its first function X_T and
second function h_T, gives

    E<X_T,Y_T>/n ->E h_T^2=1.

Consequently

    E||X_T-B_n[S h_T(X)]||^2/n ->0.                (5)

The same argument holds on the balanced Haar-involution comparator:
the old local Gaussian and marked energy identities there are proved
in `resumed_convergence_haar_response_traffic_2026_09_06.md`, Sections
2--3. Thus (5) holds on BOTH ensembles needed for universality.

To turn (5) into a globally Lipschitz first-order approximation, proceed
by induction on tree size. First replace h_T by a bounded globally
Lipschitz EVEN cutoff h_(T,M). The old local CLT and its uniform fixed
polynomial moments make the resulting normalized-L2 error arbitrarily
small on both ensembles. Then substitute the already constructed
Lipschitz approximations to the smaller child fields; their errors
propagate through the fixed Lipschitz constant and ||B_n||op=1+o(1).
Choose child tolerances after choosing M. The base tree is B_n S.

The coordinate update can be written clip(S)*h_(T,M)(children), using
an odd clip to [-1,1]. It agrees with S on the Rademacher seed support
and is globally Lipschitz in ALL its formal arguments. This avoids a
hidden unbounded product in the GFOM definition. Every approximation
preserves signed-permutation equivariance. No polynomial substitution
or high-degree moment-determinacy assertion is needed.

Finally approximate bounded F,H by bounded Lipschitz functions preserving
their respective odd/even parities. Their local Gaussian L2 errors and
boundedness transfer the approximation to the tree fields. For a fixed
odd Lipschitz softsign psi_epsilon, the terminal vector
H psi_epsilon(BF) is therefore a fixed Lipschitz GFOM limit in normalized
L2, simultaneously on the deterministic and Haar ensembles.

## 3. Matrix-sign charge and vanishing self-energy

Every canonical input tree has an odd number of edges. Recursively its
edge count is one plus the sum of an EVEN number of odd child edge
counts. Therefore, exactly before any limit,

    X_T(-B,S)=-X_T(B,S).

This is a matrix-sign charge statement, different from ordinary global
spin parity. It implies

    F(-B)=-F(B), H(-B)=H(B),
    (-B)F(-B)=BF(B), C_epsilon(-B)=C_epsilon(B).

Thus both Q_B(F) and Q_B(C_epsilon) change sign under B->-B. The
bounded Lipschitz approximations in Section 2 can all be chosen to
preserve this charge exactly.

The fixed-rule universality module is
`resumed_convergence_involution_local_ceiling_2026_09_06.md`, Sections
1--3. Its primary input is Wang--Zhong--Fan,
[Universality of approximate message passing algorithms and tensor
networks](https://arxiv.org/pdf/2206.13037), Proposition 2.7(b)(2),
Theorem 2.8, and the finite-first-order reparameterization explained in
that module. These primary statements were checked directly. Under
(1), odd powers are U_n and even powers are I, so the required power
entry bounds hold; the limiting spectral law is balanced +/-1.

For every fixed equivariant Lipschitz rule, normalized empirical energies
therefore converge to deterministic limits identical for B_n and -B_n.
Exact charge makes the corresponding limits negatives of each other.
They must be zero. Signed-permutation equivariance removes the auxiliary
random conjugation; it is not being imposed on a fixed deterministic
matrix by wishful averaging. Since feasible outputs and bounded operator
norm make normalized energies uniformly bounded, this convergence is
in L2. It proves (2) first for fixed softsigns and all Lipschitz
approximations. The next section removes the softsign.

## 4. Flat involutions identify the full terminal local channel

For a fixed odd polynomial P write P=P1P+R, with

    K_P=sum_T E[P G_T]h_T, tau_P^2=||R||^2.

Equation (5), multiplied by B_n, gives

    B_n X_T-S h_T(X) ->0 in normalized L2.          (6)

Thus B_n(P1P)=S K_P(X)+o_L2. The coherent first-chaos part that prevented
a full conditional-law claim for GENERAL signings has disappeared here.

The full nonlinear forest theorem supplies the covariance identity

    Cov(R_main)=tau_P^2 I+E_n, ||E_n||_1/n ->0.

After multiplication by B_n, the output covariance is tau_P^2 I plus
normalized-nuclear-o(1). In particular its diagonal variances converge
to tau_P^2 in average absolute error. Proper flattenings, mixed old-field
contractions, and own-root removal give joint Gaussian convergence of
this channel with the old fields, and asymptotic independence of the
root spin. Raw local-Wick errors are removed in normalized L2 before
using this statement. At an averaged uniform root the resulting law is

    (S_i,X_i,(B_nP(X))_i)
         -> (S,G,S K_P(G)+tau_P N),                (7)

where S is a Rademacher, and S,G,N are independent. This conclusion is
SPECIAL to the flat involution reduction. It is not claimed by the
general restricted-response identities alone.

Approximate bounded odd F by fixed odd polynomials P in Gaussian L2.
The bounded operator norm controls B(F-P) in normalized L2, while
K_P->K_F and tau_P->tau_F. Passing first in matrix size and then in
polynomial approximation gives (7) with F. Since a nonzero bounded odd
Gaussian function cannot be a nonzero unbounded Gaussian first-chaos
linear form, tau_F>0.

Conditional on S,G, the last coordinate of (7) has density bounded by
1/(tau_F sqrt(2pi)). Hence the averaged probability of |BF|<=epsilon
tends to at most 2epsilon/(tau_F sqrt(2pi)). This proves normalized-L2
approximation of H sign(BF) by H psi_epsilon(BF), and removes the
softsign in Section 3. An alternate convention at exact finite-n zero
fields is asymptotically irrelevant under the same small-ball bound.

The joint law (7), uniform integrability from ||B_nF||^2/n<=O(1), and
bounded H also give

    (1/n) sum_i E[H(X_i)|(BF)_i|]
         ->E[H(G)|S K_F(G)+tau_F N|]=J(F,H).

The fixed-GFOM empirical convergence gives concentration of the cross
energy; the softsign and normalized-L2 approximations preserve it.
This proves (3). Combining (2)--(3) with the exact endpoint identity
proves (4).

For hollow B_n, independent Boolean rounding of a cube endpoint does
not change these normalized limits. Conditional on a mean vector m,
write the rounded signs as m+xi, with independent centered bounded xi.
The conditional variance of the linear energy error is O(||B_n m||^2)
=O(n), and that of its off-diagonal quadratic part is O(||B_n||_F^2)
=O(n). Hence the rounded energy divided by n differs from Q_B(m)/n
by o_L2(1). In the current ternary/binary policy the two endpoints are
already Boolean except on a negligible set of exact zero fields.

## 5. Consequence, and the exact limits of the conclusion

On the matrix class (1) and its stated stable hollow perturbations,
the actual outputs of THIS fixed terminal architecture have asymptotic
energy J(F,H), not an unidentified extra common-energy contribution.
The global terminal-functional ceiling therefore gives a <.45 ceiling
for this restricted randomized architecture on normalized Hadamards.
The actual Boolean optimum of those signings may nevertheless approach
1/2. The deterministic stationary pairs in the earlier iteration
falsifiers are not fixed equivariant odd-tree terminal rules from iid
seeds, so there is no contradiction.

The theorem does not cap arbitrary GFOMs, retained-energy iterations,
matrix-adapted initialization, nonlocal methods, or coordinate counts
and depths growing with n. It does not assert a rate uniform over
terminal rules. Mere near-involution Frobenius control is explicitly
outside its scope. For a Gaussian-L2-zero representative F, or a rule
with an unverified zero-field atom, only the fixed-softsign assertion
is automatic; the hard-sign passage must be treated separately.

## 6. Charge normal form for general fixed equivariant algorithms

This algebra explains precisely where an algorithm escaping the terminal
ceiling must differ. Let mu(B,S) be ANY fixed feasible signed-permutation-
equivariant Lipschitz GFOM output, with values in [-1,1]. Run the same
finite algorithm also at -B with the same seeds, and define

    F(B,S)=[mu(B,S)-mu(-B,S)]/2,
    C(B,S)=[mu(B,S)+mu(-B,S)]/2.

These remain fixed equivariant GFOM outputs, have respectively odd and
even matrix-sign charge, and satisfy the exact pointwise identity

    |F|+|C|=max(|mu(B,S)|,|mu(-B,S)|)<=1.

The parity/universality argument in Section 3 therefore gives

    Q_B(F)/n ->0, Q_B(C)/n ->0,
    Q_B(mu)/n-F^T B C/n ->0, all in L2.

This conclusion does not require an old-tree representation of F or C.
It DOES require the same fixed-rule universality hypotheses, and does
not apply to arbitrary matrix-dependent signings.

For each fixed epsilon>0 define the even-charge feasible soft best
response C_epsilon=(1-|F|) psi_epsilon(BF), where psi_epsilon is odd,
globally Lipschitz, bounded by 1, and satisfies
z psi_epsilon(z)>=|z|-epsilon. It is another fixed equivariant GFOM and

    F^T B C_epsilon/n >= |F^T B C|/n-epsilon.

Its self-energy still vanishes. Thus the soft best-response endpoint
F+C_epsilon asymptotically dominates the absolute energy of the original
output, up to epsilon. No zero-field assumption is needed for this
soft assertion. The normal form is therefore broad, but its charge-odd
F need NOT be an unmarked function of the canonical old Gaussian frame.
It can involve root seeds and registers of both charges. The restricted
Gaussian terminal-channel law, and consequently the .45 bound, do not
extend merely by this algebraic normal form. This identifies a precise
remaining mechanism: enlarge the charge-odd half itself beyond the old
unmarked creation family; adding more old-only masks does not do so.
