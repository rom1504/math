# Nonlinear old-response transport: a forest flattening lemma

Date: 2026-09-06, 01:35 UTC. This is follow-up work after freezing
`resumed_bound_audit_minimal_proof_2026_09_06.md`. It concerns a proposed
stronger conditional-channel theorem, not a dependency of that bound.

## What is proved here, and what is not

A deterministic tensor lemma proves that transporting a product of at
least two old tree kernels through one flat matrix row gives uniformly
small proper marked-slot flattenings. Applied to higher odd local Hermite
monomials, this supplies their Gaussian-input transported CLT after the
usual fixed-degree Wick replacement. The argument is elementary once
the old marked-tree global and fixed-root flattening bounds are available.

This note does NOT claim the full conditional formula for B F. In
particular, it does not silently infer independence from one-channel
covariances when the transported first-chaos part need not be Gaussian.

## 1. Deterministic hypotheses and conclusion

Let B be hollow symmetric with |B_ij|<=C/sqrt(n), bounded operator norm,
and uniformly bounded row Euclidean norms. For a finite list of branch
types t, let K_t,j be a tensor of fixed order d_t>=1. Suppose:

1. The global row map C_t:j -> K_t,j has operator norm O(1).
2. Each row tensor has Hilbert norm O(1).
3. Every proper slot flattening of each fixed-root K_t,j has operator
   norm O(n^(-1/2)), uniformly in j. For order one this is vacuous.

The normalized injective marked odd-tree kernels satisfy these
hypotheses under a fixed normalized operator cap. The global bound is
the all-vertex tree-cut estimate; the fixed-root bound has the additional
flat external edge. The arguments are given in
`fresh_tree_matrix_derivative_and_weighted_third_chain_2026_09_05.md`,
Section 1, and `fresh_limit_unmarked_cubic_projection_2026_09_05.md`,
Section 2; both are direct tensor-cut arguments, not CLT assertions.

For a fixed k>=2, put

    K_i = sum_j B_ij (K_t1,j tensor ... tensor K_tk,j).

Then ||K_i||_HS=O(1), and EVERY proper flattening of K_i has operator
norm O(n^(-1/2)). Both assertions are uniform in i. Symmetrization of
its finitely many marked slots preserves the conclusions.

## 2. Root maps for products of whole branches

For a nonempty collection A of branches define C_A(j)=tensor_{a in A}
K_ta,j. Its Gram matrix is the Schur product of the corresponding
positive semidefinite branch Gram matrices. Their operator norms and
diagonals are uniformly bounded. The Schur multiplier inequality

    ||P circ Q||op <= (max_j P_jj) ||Q||op,   P,Q>=0,

follows by representing P as a Gram matrix and compressing Q tensor I.
Iteration proves ||C_A||op=O(1). This includes collections with repeated
tree types. It is not based on pairwise decorrelation alone.

The Hilbert norm conclusion in Section 1 now follows immediately:
K_i is the row combination b_i C_{1,...,k}, and ||b_i||_2=O(1).

## 3. All proper cuts: classify how many branches are split

Fix a nontrivial bipartition of the marked slots. A branch is called
split if it has slots on both sides. Let s count split branches. Each
split branch has fixed-root matrix operator norm O(n^(-1/2)); an
unsplit branch is a vector of bounded norm.

If s>=2, the triangle inequality and tensor-product norm identity give

    ||flatten(K_i)||op
      <= sum_j |B_ij| O(n^(-s/2))
      = O(n^((1-s)/2)) = O(n^(-1/2)).

Here sum_j |B_ij|=O(sqrt(n)); no cancellation is assumed.

If s=0, some whole branches lie on each side because the cut is proper.
The flattening is

    C_left^T diag(B_i) C_right.

Both C maps have bounded operator norm by Section 2 and the diagonal
has norm O(n^(-1/2)).

If s=1, call the split-branch matrices L_j. Since k>=2, at least one
whole branch remains. Suppose first that all whole branches lie on the
left, with root map C. Factor the flattening as (C^T tensor I) T, where

    T x = (B_ij L_j x)_j.

The stacked operator has squared norm at most

    sum_j B_ij^2 ||L_j||op^2 = O(n^(-1)),

so the required bound follows. The all-right case is its transpose.
If whole branches lie on both sides, factor instead as

    (C_left^T tensor I)
       blockdiag_j(B_ij L_j)
    (C_right tensor I).

Its norm is O(n^(-1)), an even stronger estimate. These formulas use
the same j block on both sides and so do not replace a shared root by
independent summation indices.

All possibilities are covered. Symmetrization averages a fixed number
of permutations; a fixed cut after permutation is just another cut of
the original tensor, so the same bound applies.

The distinction k>=2 is real. For k=1, the single branch can be split
with NO whole branch left. The last factorization is unavailable, and
sum_j B_ij K_t,j can have an order-one proper flattening. This is exactly
the first-chaos transport exception, not an overlooked case of the lemma.

## 4. Gaussian-input consequence for higher local Hermite monomials

Let X_T,j=I_dT(K_T,j) be the fixed injective Gaussian old-tree chaoses,
normalized to limiting variance one. For a local normalized Hermite
monomial with multiplicities m_T and total local Hermite degree
k=sum_T m_T>=2, fixed-degree Gaussian product expansion gives

    product_T He_mT(X_T,j)/sqrt(m_T!)
      = I_D(sym tensor_T K_T,j^(tensor m_T))
          /sqrt(product_T m_T!) + error_j,

where D=sum_T m_T d_T, and max_j ||error_j||_2=o(1).
The tensor convention includes each branch's original normalization.
Whole contractions of equal branches are precisely the Hermite
subtractions. Full contractions of distinct tree types are o(1).
Every other contraction is bounded by a proper fixed-root flattening
times a bounded Hilbert norm. There are only finitely many contractions.
This establishes the displayed error for every fixed monomial; fixed
higher moments follow by hypercontractivity.

Transport preserves the error in AVERAGED L2:

    n^(-1) E||B error||_2^2
       <= ||B||op^2 n^(-1) E||error||_2^2 = o(1).

The new main kernel is precisely the kernel in Section 1, apart from
fixed normalization constants. For each proper contraction index r,

    ||K_i tensor_r K_i||_HS
       <= ||flatten_r K_i||op ||K_i||_HS = O(n^(-1/2)).

Thus its fourth cumulant tends uniformly to zero. Finitely many such
chaoses, together with finitely many old chaoses and unmarked-star
chaoses, are jointly asymptotically Gaussian along any subsequence on
which their covariance matrices converge. Different degrees may be
kept as separate vector coordinates and then summed. Bounded variances
ensure covariance subsequences exist; the conclusion permits singular
covariance limits.

The primary theorem used for this last step is Theorem 1.1 of
[Noreddine--Nourdin, *On the Gaussian approximation of vector-valued
multiple integrals*](https://arxiv.org/pdf/1009.1310), inspected directly
on 2026-09-06. Its assumptions are fixed chaos orders, convergence of
covariances, and vanishing component fourth cumulants. No covariance
nondegeneracy is required for that qualitative conclusion. This theorem
alone does not assert that any of the limiting covariances vanish.

For an odd polynomial R with zero first LOCAL Gaussian-chaos projection,
all its Hermite monomials have odd k>=3. Therefore this argument applies
to B R(X). The resulting Gaussian may a priori have an odd linear drift
correlated with the old Gaussian vector; that is a separate covariance
calculation, not needed to state the lemma.

## 5. What remains before a conditional-channel theorem is justified

The proposed target writes an odd old response as F=P1F+R and seeks,
for a nonlinear channel Z aligned with the edge-only Hermite projection,
a conditional regression of B F on (S_i,oldX_i,Z_i).

For B R, joint Gaussianity from Section 4 reduces the problem to
covariances, provided the Rademacher transfer and collision deletion are
also completed. The existing deterministic-weighted one-Z theorem can
identify the averaged covariance with Z: uniformity over bounded
deterministic weights even implies an averaged absolute covariance error
by choosing their signs. Old drift may be harmless for the intended
even-old tests, but it should not be declared absent without proof.

The first-chaos transport B(P1F) is different. Its representative
B X_T is close to Q[S h_T(X)] and need not be Gaussian. Its prominent
own-root term is S_i h_T(X_i). A Gaussian CLT for Z plus zero covariance
with B X_T does NOT imply independence or conditional linearity.
For example, in an abstract Gaussian model Y=Z^2-1 has zero covariance
with Z but nonconstant conditional mean. With parity matching, a product
Y=U(Z^2-1) creates the analogous mixed obstruction.

Most mixed tensor contractions between a Gaussianizing Z kernel and an
arbitrary bounded kernel vanish by Cauchy--Schwarz and Z's proper
self-contractions. The exceptional case contracts ALL slots of Z into
a higher-degree B X_T kernel. That leaves a non-scalar tensor, and
small self-contractions of Z do not by themselves control it. A kernel
of the form sym(K_Z tensor K_U) illustrates why. This full-contraction
case needs its own tree/flatness estimate before the desired tests
E[M_even psi_odd(Z) B F] and E[S M_even psi_even(Z) B F] are proved.

No counterexample in the actual signing class has been found here.
The forest lemma is affirmative progress on the nonlinear part, not
an endorsement of the still-unproved full conditional theorem.

## 6. Subsequent progress: the first-chaos root-spin tests

Added during the same resumed audit. Let L_i=(B X_T)_i. The one-level
representation and fixed operator cap give, in averaged L2,

    L_i = sum_j Q_ij S_j h_T(X_j) + o(1),  Q=B^2.

Let M be a fixed polynomial of old local fields and psi a fixed polynomial
of the unmarked odd channel Z_i. For this paragraph its fixed coordinate
derivatives have Lp norms O(n^(-1/2)), as proved for unmarked stars and
finite mixtures. Remove S_i from Z_i at Lp cost O(n^(-1/2)); the local
old fields already exclude it exactly. Put f_i=M(X_i)psi(Z_i), now
exactly own-i-free.

For j!=i, h_T(X_j) is own-j-free and the exact two-spin identity gives

    E S_i S_j f_i h_T(X_j)
      = E[(D_j f_i)(D_i h_T(X_j))].

Here D denotes the half-difference derivative on the Rademacher cube;
after differentiating, each factor is independent of both relevant own
coordinates. Each derivative has L2, and every needed higher Lp, norm
O(n^(-1/2)). Hence the expectation is O(1/n). Since
sum_j |Q_ij|<=sqrt(n)||Q_i||_2=O(sqrt(n)), the averaged off-diagonal
sum vanishes. The diagonal Q_ii=1 contributes

    E[M(X)h_T(X)] E[psi(sqrt(v_i)N)] + o(1),

by the uniform same-root old/Z joint Gaussian approximation. Thus the
desired root-spin channel test is proved, in particular for M even and
psi even. This argument does not falsely declare L_i Gaussian and is
unaffected by coherent finite-spin atoms in Q S.

## 7. A route through the exceptional full contraction

The following Gaussian calculation closes the exceptional contraction
identified in Section 5, conditional only on the already proved old-tree
directional-derivative matrix bound. It is recorded in detail so the
director can independently check and integrate it.

Work with a single odd unmarked-star order r>=3, with kernel

    Z_i=I_r(K_Z,i),  K_Z,i proportional to sum_l B_il b_l^(tensor r).

Fixed normalizing factorials are immaterial to the vanishing estimates.
Let H_a=h_T(X_a) and write the first-chaos transport as

    Y_i=sum_a Q_ia N_a H_a.

Every H_a excludes N_a exactly. Contracting all r slots of K_Z,i into
the Gaussian polynomial Y_i is, up to fixed constants,

    sum_l B_il D_bl^r Y_i.

Its product rule has two classes.

The root-not-hit class is

    sum_a Q_ia N_a U_a,i,
    U_a,i=sum_l B_il D_bl^r H_a.

One has max_{a,i}||U_a,i||_2=o(1). Here is a route to this fact which
does not assume generic asymptotic independence. For ANY two output
roots a,i, the finite old family at a and Z_i are jointly Gaussian in
the limit. Proper self-contractions of every component vanish. The
cross covariances vanish uniformly: unequal original chaos orders are
orthogonal; at equal order, the recursive contraction estimate for a
nonstar tree has global Frobenius norm O(n^(-1/2)), and the star case
is the matrix B[B circ Q^(circ(r-1))]B, whose middle factor has operator
norm O(n^(-1/2)). Collision deletion costs O(n^(-1/2)) in the covariance
because it pairs a fixed-root Hilbert error with the bounded Hilbert
norm of K_Z,i. Only uniform ENTRYWISE covariance control is needed
here; no claim about the exact collision-error covariance operator is
being made.

Consequently H_a, being a fixed polynomial of the old family, is
asymptotically independent of Z_i in all joint moments, uniformly in
(a,i). Replace H_a by its pure even-chaos Wick kernel with o(1) L2
error as in Section 4. For two multiple integrals, covariance of their
squares is a sum of nonnegative contraction-norm terms and in particular
dominates a fixed positive constant times the squared norm of the full
Z contraction. Its limit is zero by the joint moment independence.
This proves the assertion for U_a,i; finite-degree derivative contraction
is L2-continuous with a bounded K_Z,i kernel, so the Wick replacement
error is harmless. If H_a has degree below r, the contraction is zero.

The U_a,i remain own-a-free. The Gaussian root-creation inequality,
the analogue of the Boolean root-multiplication bound, is

    ||sum_a q_a N_a P_a||_2^2
       <= (d+1) sum_a q_a^2 ||P_a||_2^2

when all P_a exclude N_a and have degree at most d. It follows by
symmetrizing the added e_a slot in the Hermite coefficient tensor;
there is no annihilation term because of own-coordinate exclusion.
Applying it with q_a=Q_ia, whose Euclidean norm is uniformly bounded,
shows that the root-not-hit class tends uniformly to zero in L2.

The root-hit class, apart from its factor r, is

    sum_{a,l} Q_ia B_il B_la D_bl^(r-1) H_a
      = diag[Q (B circ K_(r-1)) B]_i,

where K_(r-1)(a,l)=D_bl^(r-1)h_T(X_a). The earlier matrix derivative
theorem gives || ||K_(r-1)||op ||_p<=polylog(n) for every fixed p.
The assumption r>=3 ensures at least one derivative, as required by
that theorem. Flatness gives

    ||B circ K_(r-1)||_F
       <= C n^(-1/2)||K_(r-1)||_F <= polylog(n).

Since ||Q||op and ||B||op are bounded, the squared L2 norm averaged
over i of the displayed diagonal is O(polylog(n)^2/n)=o(1).
This estimate uses averaged roots, not an unjustified uniform-root
operator-to-entry improvement.

Thus every full-Z mixed contraction with B X_T vanishes in averaged
squared Hilbert norm. Proper-Z mixed contractions already vanish by
the Cauchy--Schwarz contraction inequality and Z's proper
self-contractions. The one-level replacement is valid in averaged L2;
contracting against bounded K_Z,i preserves that error.

The standard Gaussian product formula now makes Z asymptotically
independent, in fixed polynomial moments and averaged roots, of the
joint collection consisting of B X_T and the old fields. This does NOT
say B X_T itself has a Gaussian limit. In particular, for odd psi,

    n^(-1) sum_i E[(B X_T)_i M(X_i) psi(Z_i)] -> 0.

Finite mixtures of odd r>=3 are handled by keeping their star-chaos
coordinates separate throughout the contraction argument. Their
joint Gaussian limit permits the final fixed polynomial of their sum.

## 8. Why high-influence first-chaos transport does not block input transfer

For the preceding expectation, Y_i=(B X_T)_i enters LINEARLY. It is a
squarefree input polynomial. The local test f_i=M(X_i)psi(Z_i), after
the standard squarefree unmarked-star replacement, is a polynomial of
low-influence squarefree fields. Its repeated derivative in one input
coordinate has fixed Lp bound O(n^(-k/2)) at order k, under Gaussian,
Rademacher, and mixed product inputs. Hypercontractivity supplies the
fixed higher moments; these are not pointwise derivative assertions.

In a fourth-order Lindeberg replacement of E[Y_i f_i], the first three
input moments match. Since D_j^2 Y_i=0 exactly, the fourth derivative
has only

    Y_i D_j^4 f_i + 4(D_j Y_i)D_j^3 f_i.

The first class sums to O(n^(-1)) times the relevant Y_i L2 bound.
For the second, use total influence
sum_j ||D_jY_i||_2^2<=degree(Y_i)||Y_i||_2^2, followed by
Cauchy--Schwarz over j. The sum is again O(n^(-1)) times ||Y_i||_2.
These bounds hold in averaged roots even if a particular root is
exceptional, because the averaged second moment of B X_T is bounded.
Taylor remainders under the mixed input laws are controlled by the
same fixed-degree moment estimates.

Therefore these TEST EXPECTATIONS transfer without asserting a
Gaussian law for high-influence Y_i or for Q S. The root-spin version
can instead use the direct Boolean argument of Section 6, retaining
S_i exactly rather than replacing its distribution.

For nonlinear B R, the forest main kernels have low influences. Delete
cross-branch repeated slots before transport: at each local root their
Hilbert norm cost is O(n^(-1/2)), by entry flatness and the O(n^(D-1))
forbidden tuples. B transports this error in averaged Hilbert norm.
The remaining squarefree main kernels have vanishing averaged maximum
influence by Section 3, and ordinary fixed-vector Lindeberg invariance
applies. Output-own-coordinate removal has vanishing averaged error
for the same reason. The local Rademacher Hermite-to-forest replacement
can also be proved directly from the marked-tree moment partitions:
the restricted leading pairings are precisely the Hermite/Wick ones.

These additions substantially narrow the original conditional-channel
gap. They are an explicit proof route for the two restricted test
families requested by the director, with the inherited directional
matrix bound and standard product-formula contraction identities made
visible. They have not yet been promoted here to a fully packaged
conditional-expectation theorem or a new numerical original-problem
lower bound.
