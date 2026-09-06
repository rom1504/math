# Arbitrary Fourier matching permutations compute the majorant norm

Date: 2026-09-06. Status: a general norm identity for an ENLARGED
operator model, and a precise sufficient realization target for R=T.
The matching permutations used below have NOT been realized by the
prescribed Hadamard outer family. Therefore this is not an R=T theorem
and does not change the original convergence interval.

## 1. Statement, including the restricted permutation support

For a fixed real symmetric k-by-k seed B let

    T(B)=(1/2) min_(D>=B,D>=-B) max_(x Boolean) x^T D x.

The previously established cut-covariance dual is

    2T(B)=max_(C in CUT_k) ||C^(1/2) B C^(1/2)||_* .        (1)

Let H_m be normalized Walsh on F_2^m. Write frequencies as (a,b)
in F_2 x F_2^(m-1), and let P range over permutations supported on
the affine hyperplane a=1, fixing EVERY frequency with a=0.
In particular P fixes the zero frequency. Put U_P=H_m P H_m.

Define the enlarged balanced-profile norm

    N_perm(B)=(1/2) sup_(m,P,F,G)
                              |E_x G(x)^T B (U_P F)(x)|,

where F,G have k Boolean columns, all with exact mean zero. Then

    N_perm(B)=T(B).                                         (2)

The upper bound holds for arbitrary orthogonal U, even without any
balance restriction. The substantive claim is that permutations on
one affine frequency hyperplane already suffice for the lower bound.

## 2. Upper bound for every orthogonal outer operator

Let D dominate B and -B. On its support factor

    B=D^(1/2) M D^(1/2),       ||M||_op<=1.

For arbitrary Boolean F,G and orthogonal U, Cauchy-Schwarz gives

    |E G^T B U F|
       <=sqrt(E F^T D F) sqrt(E G^T D G)
       <=max_(x Boolean) x^T D x.

Minimizing D proves N_perm(B)<=T(B). This argument does not depend
on symmetry of U or on the zero-frequency constraint.

## 3. The correct covariance-preserving polar map

Fix C in CUT_k and write S=C^(1/2) B C^(1/2). Let O be sign(S)
on its nonzero eigenspaces, extended to a symmetric orthogonal
operator on the zero eigenspaces. Choose the extension to commute
with the range projection Pi_C. Set

    A=C^(dagger 1/2) O C^(1/2).                             (3)

Here dagger denotes the Moore-Penrose inverse. Direct multiplication
on the range of C gives

    A^T C A=C,
    C A=C^(1/2) O C^(1/2),
    Tr(C A B)=||S||_* .                                    (4)

This is the ROW-vector convention: row covariance is transformed
by A^T C A. Reversing the two square-root factors is generally wrong.
If C is singular, every random vector of covariance C and mean zero
lies in its range almost surely, so no undefined kernel action is
used in the construction.

## 4. Empirical Walsh clouds converge in quadratic transport

Choose a symmetric Boolean law X in {+1,-1}^k with covariance C.
Such a law is obtained from any cut representation of C by assigning
equal mass to each atom and its negative. Let N=2^r and let X_z,
z in F_2^r, be independent copies. Define the unitary-normalized
Walsh rows

    Z_b=N^(-1/2) sum_z (-1)^(b dot z) X_z.                  (5)

The empirical measure N^(-1) sum_b delta_(Z_b) converges in
probability, in quadratic Wasserstein distance, to N(0,C).

Here is an elementary proof of the ingredients. Because X is
symmetric, the marginal law of Z_b is the same for every b and
is the usual normalized sum of iid bounded vectors. For b!=c,
multiply each X_z by its harmless sign (-1)^(b dot z). The pair
(Z_b,Z_c) then has the same law as the sum with coefficient pairs
(1,(-1)^((b+c) dot z)). Exactly half the signs are plus and half
minus. Its covariance is diag(C,C). The bounded-summand
characteristic-function expansion, with remainder O(N^(-1/2)),
therefore proves convergence to two independent N(0,C) vectors,
uniformly over distinct b,c.

For every bounded continuous test function, its empirical average
thus has the correct limiting mean and variance tending to zero.
This proves weak empirical convergence in probability. Parseval
gives the exact second-moment identity

    (1/N) sum_b ||Z_b||^2=(1/N) sum_z ||X_z||^2=k.          (6)

The limiting second moment is also Tr(C)=k. Weak convergence
together with convergence of these moments gives quadratic
transport convergence. This last implication can be proved by
truncating to a large ball, partitioning it into finitely many small
cubes, matching the mass in each cube, and bounding the remaining
cost by the second-moment tails. No independence of all Walsh rows
is claimed or needed.

The same conclusion holds for an independent Boolean sample Y_z,
and after multiplying its Walsh rows on the right by the fixed A
in (3), since A^T C A=C. Its empirical second moment converges to
Tr(A^T C A)=k by the ordinary bounded-vector law of large numbers.

## 5. Exact balance and matching only one frequency hyperplane

On F_2 x F_2^r define the exactly balanced Boolean profiles

    F(b,z)=(-1)^b X_z,       G(b,z)=(-1)^b Y_z.              (7)

Their normalized Walsh transforms vanish on the frequency
hyperplane a=0. On a=1 they are sqrt(2) times the Walsh clouds
of X and Y respectively. By Section 4, the two active point clouds

    (H F)(1,c),       (H G)(1,c) A

converge in quadratic transport to the same Gaussian law N(0,2C).
There consequently exist permutations pi of the N active rows
such that

    (1/(2N)) sum_c
       ||(H F)(1,pi(c))-(H G)(1,c)A||^2 ->0.              (8)

One can justify permutation matching without a transport theorem:
match as many points as possible inside each cell of a fixed fine
finite partition of a large ball, match all leftovers arbitrarily,
then use converging cell counts and second-moment tails. Equivalently,
optimal transport between equal-weight empirical measures is an
assignment problem and has a permutation optimizer.

Extend pi by the identity on a=0 and call its matrix P. Applying
Walsh orthogonality to (8) yields

    E ||U_P F-G A||^2 ->0.                                 (9)

Moreover E G^T G converges to C. Since B is fixed and G is Boolean,
the error (9) contributes at most
sqrt(k)||B||_op sqrt(E||U_P F-GA||^2) to the bilinear objective.
Equations (4) and (9) therefore imply

    E G^T B U_P F ->Tr(C A B)=||S||_* .                    (10)

Choose C maximizing (1). This proves the lower bound in (2).
The random construction only proves existence of deterministic
profiles and matching permutations along a subsequence; no random
outer matrix is substituted for the original problem.

## 6. Exact consequence and the unresolved realization gap

Suppose the prescribed regularized Hadamard family could realize
U_P=H_m P H_m as asymptotically isometric balanced spin-profile
tests for every permutation P supported on a=1. Then (2), together
with R<=T, would immediately prove

    R(B)=T(B) for every finite symmetric B.                 (11)

In particular the majorant floor k^(3/2)/2 for full sign seeds would
be unavoidable under this regularization. This would settle the
regularized-seed realization question, not automatically the
original convergence of M_n/n^(3/2).

What is actually proved so far is much narrower: independently
programmable projective Fourier PHASES, and conjugated PERMUTATIONS
arising as dual-index maps of certain vectorial dual-bent functions,
including field power permutations. Neither theorem realizes the
data-dependent arbitrary matching permutations used in (8).

Thus the current exact gap is an operator/Boolean realization gap:
can the available nonlocal constructions implement, or sufficiently
approximate on the selected finite Boolean profiles, the matching
permutations from (8)? It is not a remaining convex-duality gap,
Gaussian-covariance gap, or normalization issue.
