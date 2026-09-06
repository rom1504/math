# Exact same-space recovery from matrix height, including critical closure

Date: 2026-09-06. Status: director proof, independently reconstructed by
both the response researcher and the proof auditor.
This is a realization theorem for the proved Gaussian response dictionary,
not a realization theorem for arbitrary near-minimizing sign matrices.

## 1. Precisely specified objects

Use the countable canonical Gaussian tree space and its creation isometry
U from even L2 functions onto the first Gaussian chaos. The normalized
even Hermite monomial h_T of a tree's children maps to its Gaussian G_T.
Let F_h be the sigma-field of all tree coordinates of height at most h,
and Pi_h its first-chaos projection. F_0 is trivial. Fixed height can have
infinitely many coordinates. The elementary basis identity is

    Pi_(h+1) U f = U E[f | F_h].                         (1)

Let g:R^d -> R^d be measurable, jointly even, with

    E[g(Z)g(Z)^T]=I_d,     Z standard Gaussian in R^d.

No differentiability is assumed. For 0<=P<=I define

    K(P)=E[g(X)g(Y)^T],

where X,Y are jointly standard Gaussian with cross covariance P. These
are uncentered second moments: K(0)=Eg(Eg)^T. Equivalently K(P) is the
Gram matrix of conditional expectations when both copies share their
Gaussian component of covariance P. Consequently K is Loewner monotone,
K(I)=I, and it is continuous on this compact matrix interval. Continuity
follows by L2 approximation of g by polynomials, uniformly for every
coupling because both marginals remain standard Gaussian.

## 2. Recovery theorem

The following are equivalent:

1. There exists an orthonormal first-chaos frame G in THIS canonical
   space such that G=U g(G), component by component.
2. Starting P_0=0, the iterates P_(h+1)=K(P_h) tend to I.
3. There is no proper matrix 0<=P<=I, P!=I, with K(P)<=P.

If it exists, the frame in (1) is unique in the specified canonical
space. Its height-h projection is constructed explicitly below, and

    E ||G-Pi_h G||^2 = tr(I-P_h).                      (2)

The necessary direction and the subfixed obstruction were separately
proved in `resumed_response_matrix_height_feedback_obstruction_2026_09_06.md`.
Here is a constructive converse; no recursive-tree-process theorem is
being imported without a mapping.

Set Y_0=0. Inductively, with Cov(Y_h)=P_h, define the even function

    f_h = E_Z g(Y_h+(I-P_h)^(1/2) Z),
    Y_(h+1) = U f_h,                                   (3)

where Z is an external standard Gaussian used ONLY in this integral.
The resulting f_h is a function in the actual old Gaussian space, so
U f_h uses no external randomness. The joint law of two independent
integrations sharing Y_h proves Cov(Y_(h+1))=K(P_h).

More importantly, the variables in (3) have compatible projections:

    Y_h in F_h,       Pi_h Y_(h+1)=Y_h.                (4)

For the induction, (1) reduces (4) to E[f_h|F_(h-1)]=f_(h-1).
By the preceding instance of (4), Pi_(h-1)Y_h=Y_(h-1). Since these
are first-chaos Gaussian vectors, their difference is independent of
F_(h-1), with covariance P_h-P_(h-1). Integrating that difference
together with the external Gaussian in (3) produces exactly covariance

    (P_h-P_(h-1))+(I-P_h)=I-P_(h-1),

which proves the required conditional identity. The induction also
proves that P_h increases and remains between 0 and I.

If P_h tends to I, (4) shows that Y_h is L2 Cauchy, with first-chaos
limit G of covariance I. Its projection onto F_h is Y_h. Therefore
the conditional Gaussian law of G given F_h has mean Y_h and covariance
I-P_h. In particular f_h=E[g(G)|F_h]. Martingale convergence gives

    U g(G)=lim_h U f_h=lim_h Y_(h+1)=G.

This proves (2) implies (1), for merely L2 inverse features.

Conversely, for any such frame, applying (1) to its fixed-point equation
forces its successive projections to follow (3). This proves both
necessity and uniqueness. No independent standard Gaussian frame at
height zero was assumed.

Finally the monotone P_h converge to some P_infinity. Continuity gives
K(P_infinity)=P_infinity. Thus failure of (2) supplies a proper subfixed
matrix. Every proper subfixed matrix bounds every P_h by monotonicity,
so (2) and (3) are equivalent.

## 3. Anchored scalar critical closure

Suppose a finite ancestor-closed anchor bank is fixed, and a proposed
unit innovation h(A,z) is orthogonal to all anchor inverse features.
Expand h only in the innovation z:

    k(q)=E[h(A,X)h(A,Y)]=sum_l w_l q^l,
    w_l>=0, sum_l w_l=1,

where X,Y have correlation q and share the same anchors. If

    k(q)>q for every 0<=q<1,                            (5)

then a unique same-space Gaussian innovation exists. To justify the
anchored reduction in the sufficient direction, the full canonical
anchor coordinates have deterministic finite heights. In the full
iteration (3), their projections become EXACT after that finite
height. At that stage the anchor block of P_h is I and all cross
blocks are zero, so the remaining scalar recursion is k. Condition
(5) forces its limit to one. Apply Section 2.

In particular, assume

    sum_l l w_l=1,
    some w_l>0 with l>=2.

Strict convexity on (0,1), k(1)=1, and k'(1)=1 give k(q)>q there.
The mean-one identity and higher coefficient imply w_0>0, so (5)
also holds at zero. Thus CRITICAL derivative energy one is allowed;
a strict contraction margin is sufficient but not necessary.

If h has finite innovation degree and k''(1)=v>0, and the scalar
height deficit is e_h=1-q_h, Taylor expansion gives

    e_(h+1)=e_h-(v/2)e_h^2+O(e_h^3),
    e_h ~ 2/(v h).                                    (6)

The shift by the finite anchor height does not affect this asymptotic.
Indeed e_h tends to zero, and
1/e_(h+1)-1/e_h -> v/2; Cesaro summation proves (6).
For derivative energy greater than one, k(q)<q near one, giving a
proper subfixed matrix and ruling out this same scalar closure.

## 4. What has and has not been recovered

Archive comparison matters: the UNANCHORED scalar D<=1 theorem and its
uniqueness were already proved in
`fresh_tree_fixed_point_critical_independent_2026_09_05.md`, by a
Galton--Watson coefficient recursion. Critical scalar existence alone
is not new. The present addition is the full matrix necessary-and-
sufficient recovery theorem for measurable inverse features, its explicit
conditional-Gaussian construction, and the anchored multifeature scope.
The finite-polynomial critical tail in (6) is a quantitative corollary.

This theorem replaces an assumed contraction-based existence step by
an exact matrix criterion, with a constructive height approximation and
uniqueness. It admits critical cyclic inverse features while rejecting
specific unstable feedback. Finite height is not finite description:
after fixing h, ordinary Hermite approximation is still needed for a
finite matrix algorithm when g is not polynomial. All such choices
precede the matrix-order limit in the original response theorem.

It does not identify the optimum of the full response functional, remove
its fixed-rule ceiling, realize arbitrary action-limit sign matrices,
or establish convergence of M_n/n^(3/2).

For comparison, Mach--Sturm--Swart, *A new characterization of endogeny*,
[arXiv:1801.05253v4](https://arxiv.org/pdf/1801.05253), Theorem 1, relates
endogeny to bivariate uniqueness for genuine stochastic recursions on
independent child copies. The primary definitions and theorem were read.
Our U is not automatically such a recursion; (1)--(4) supply the actual
same-space construction instead of presuming that imported hypothesis.
