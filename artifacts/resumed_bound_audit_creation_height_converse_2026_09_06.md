# Exact converse to covariance-height exhaustion

Date: 2026-09-06. Director proposal, independently reconstructed here.
This removes a realization obligation for a specified Gaussian inverse
map. It does not prove the original signing minima converge.

Let U be the canonical even-Hermite to Gaussian-tree creation isometry,
F_h the Gaussian sigma-field of canonical tree coordinates of height at
most h, and Pi_h the corresponding first-chaos projection. Let
g:R^d->R^d be measurable, globally even, and satisfy E[g(N)g(N)^T]=I
for a standard Gaussian d-vector N. No continuity or Sobolev hypothesis
is imposed. Define

    K(P)=E[g(X)g(Y)^T],   0<=P<=I,

where X,Y are standard Gaussian d-vectors with cross covariance P.
This is an uncentered second-moment kernel. Set P_0=0 and
P_(h+1)=K(P_h).

## Theorem

The following are equivalent:

1. There exists an orthonormal Gaussian first-chaos frame G in this
   SAME canonical space satisfying G=U g(G), componentwise.
2. P_h tends to I.
3. There is no proper PSD matrix 0<=P<=I with K(P)<=P.

When such a frame exists, it is unique almost surely for the fixed map g
and the fixed canonical creation isometry U.

## Proof

Necessity is the audited identity

    Pi_(h+1) U f=U E[f|F_h].

It gives the exact height recursion for every endogenous solution, whose
height projections exhaust its first chaos. Monotonicity of K by Gaussian
conditional variance traps this recursion under every subfixed P, proving
1=>2=>3.

Conversely, K is continuous on the compact matrix interval [0,I]. To
check this even at singular Gaussian couplings, approximate g in its
standard Gaussian L2 norm by polynomials. Cauchy--Schwarz makes the
resulting kernel error uniform over P because the two marginals remain
standard; polynomial kernels are continuous. Thus the monotone sequence
P_h has a limit P_infinity satisfying K(P_infinity)=P_infinity. If (3)
holds, that fixed point must be I, proving 3=>2.

To construct the actual frame from (2), set Y_0=0 and recursively define

    f_h(y)=E_Z g(y+(I-P_h)^(1/2) Z),
    Y_(h+1)=U[f_h(Y_h)].                           (1)

The map f_h is even, so creation is applicable. Inductively Y_h is a
centered Gaussian first-chaos vector of covariance P_h, measurable in
F_h. Isometry and the shared-mean Gaussian representation give
Cov(Y_(h+1))=K(P_h)=P_(h+1).

The essential stronger induction is compatibility:

    Pi_h Y_(h+1)=Y_h.                              (2)

For the induction step, condition Y_h on F_(h-1). Its mean is Y_(h-1)
by the preceding compatibility, and its residual is independent Gaussian
with covariance P_h-P_(h-1). Combining that residual with the independent
Gaussian in f_h yields total covariance I-P_(h-1). The creation-height
identity now gives

    Pi_h Y_(h+1)
      =U E[f_h(Y_h)|F_(h-1)]
      =U f_(h-1)(Y_(h-1))=Y_h.

Hence these fields have orthogonal successive height increments. Since
P_h tends to I, Y_h converges in Gaussian L2 to a first-chaos frame Y
of covariance I. Compatibility implies Pi_h Y=Y_h. Conditional on F_h,
Y therefore has mean Y_h and deterministic residual covariance I-P_h.
Consequently

    f_h(Y_h)=E[g(Y)|F_h].

Martingale L2 convergence gives f_h(Y_h)->g(Y), and isometry in (1)
gives Y=U g(Y). This proves sufficiency without a continuity assumption
on g or an unproved limit interchange through a discontinuity.

Finally, any endogenous solution must have the recursively prescribed
height projections (1): start at zero and apply the exact conditional
Gaussian law at each height. Thus every such solution has Pi_h G=Y_h
for all h. Height exhaustion proves uniqueness. QED.

## Anchors and the scalar critical case

For a finite canonical ancestor-closed anchor bank, adjoin its inverse
features to the innovation inverse map and apply the full theorem. The
constructed anchor coordinates are exactly their canonical coordinates
once their finite heights are reached. Thereafter the covariance block
is diag(I_anchor,Q_h): the cross block vanishes because a PSD matrix
I-P_h with a zero anchor diagonal block has zero cross block. Therefore
the anchored innovation criterion is exactly the absence of a proper
subfixed matrix Q for K_A(Q), not merely a necessary test.

For one innovation K_A(q)=sum_l w_l q^l, w_l>=0 and sum w_l=1.
The criterion becomes K_A(q)>q for every q<1. In particular a genuinely
nonlinear critical map with sum l w_l=1 passes: convexity is strict
unless K_A(q)=q identically. The pure linear kernel K_A(q)=q fails,
as does derivative energy greater than one, including infinity.
This recovers the critical scalar extinction distinction; the new
converse also covers matrix-valued and merely measurable inverse maps.

For infinitely many frame coordinates, the construction and exhaustion
would require an appropriate operator-valued formulation. The finite-d
statement above does not silently assert that extension.

## Arbitrary actual cyclic anchors: relative-height reconstruction

The finite canonical ancestry restriction in the preceding anchor
paragraph is unnecessary. This extension independently checks Section 4
of `resumed_response_matrix_recovery_independent_audit_2026_09_06.md`.

Let A=U a(A) be any actual finite orthonormal Gaussian frame, including
a cyclic frame with infinitely many canonical coefficients. Inside
canonical first chaos define closed subspaces

    V_0=span(A),   V_(h+1)=U L2_even(sigma(V_h)).

Creation closure gives V_0 subset V_1; induction gives all subsequent
inclusions. The isometry makes each range closed. Since U1 is in V_1,
and each subsequent creation step includes the next canonical tree
height, the union is dense in the complete canonical first chaos.
For every jointly even f, the exact orthogonal projection formula is

    Proj_(V_(h+1)) U f = U E[f|sigma(V_h)].

Let h(A,z) be a finite vector of jointly even inverse functions with
E[hh^T]=I and E[a h^T]=0 under independent standard A,z. Define K_A(Q)
by sharing A and coupling the innovation vectors with covariance Q.
Set Z_0=0, Q_0=0 and

    Z_(h+1)=U E_N h(A,Z_h+(I-Q_h)^(1/2)N).

The joint Gaussian vector (A,Z_h) has independent blocks of covariance
I and Q_h. Inverse orthogonality therefore keeps Z_(h+1) perpendicular
to A. In particular Proj_(V_0)Z_1=0. The same Gaussian smoothing
calculation used above now proves, inductively,

    Proj_(V_h)Z_(h+1)=Z_h,   Cov(Z_h)=Q_h,
    Q_(h+1)=K_A(Q_h).

The constant anchors are measurable at every relative level, so the
conditional integration combines only the innovation covariances
Q_h-Q_(h-1) and I-Q_h. No finite canonical support of A is required.
Relative-height exhaustion and the L2 martingale argument prove that
there is an actual extension Z=U h(A,Z), independent of A, exactly
when Q_h tends to I; it is unique with A fixed. Equivalently, there is
no proper PSD innovation subfixed matrix K_A(Q)<=Q. In the successful
case the exact finite-level error covariance is I-Q_h.

## Quantitative stability of a strict subfixed obstruction

The director's quantitative extension also reconstructs directly.
Suppose Q is positive definite, Q<=I and Q!=I, and

    K_A(Q)<=a^2 Q,   0<=a<1.

For any actual standard first-chaos innovation frame Z perpendicular
to A, put E=Z-U h(A,Z). If Cov(E)<=b^2 Q with a+b<=1, the actual
relative-height projections Z_h=Proj_(V_h)Z have P_0=0 and obey

    sqrt(v^T P_(h+1)v)
      <=sqrt(v^T K_A(P_h)v)+sqrt(v^T Cov(E)v).

This is ordinary directional L2 Minkowski after the exact projection
identity. Monotonicity then traps P_h<=Q at every height, contradicting
P_h->I. Thus even an approximate realization has a positive residual
floor whenever a strict subfixed contraction is available.

For h_1=sqrt(r) A Z_2+sqrt(1-r)(A^2-1)/sqrt(2),
h_2=A sign(Z_1), with the single edge anchor A and 0<r<=1, take

    Q=diag(1-3r^2/8,1-r/2).

Its diagonal gaps Q-K_A(Q) are at least r^2/8 and r/20. For the latter,
arccos(1-3r^2/8)>=sqrt(3)r/2 and sqrt(3)/pi>11/20; the last inequality
follows already from pi<22/7 and 140 sqrt(3)>242. Hence
K_A(Q)<=(1-r^2/20)Q and Q>=I/2. Every actual frame consequently satisfies

    ||Cov(E)||_op^(1/2)
      >=(1-sqrt(1-r^2/20))/sqrt(2)
      >=r^2/(40 sqrt(2)).

This is a quantitative obstruction for this specified feedback map,
not a claim that all stationary response maps contain such a cycle.
