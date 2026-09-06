# A quantitative barrier for approximately realizable Gaussian feedback

Date: 2026-09-06. Status: director derivation independently reconstructed
by the bound-audit researcher. This is about actual Gaussian frames in the established creation
space, not a claim about all signings or their minimizing energy.

## 1. A stable version of the subfixed-matrix obstruction

Let A be an actual finite creation-closed orthonormal Gaussian anchor
frame. Its inverse features are measurable in A. Let Z be any actual
orthonormal Gaussian d-frame orthogonal to A in the canonical first
chaos. Let h(A,z) be jointly even, of identity second-moment matrix,
and orthogonal to the anchor inverse features. Define

    K_A(P)=E[h(A,X)h(A,Y)^T],

where the two standard Gaussian innovation copies share cross
covariance P and the same anchors. No equation Z=Uh(A,Z) is assumed.
Define its actual residual vector

    E=Z-U h(A,Z),       S=E[E E^T].                    (1)

Suppose there is a proper positive definite matrix 0<Q<=I, Q!=I,
and a number 0<=a<1 such that

    K_A(Q)<=a^2 Q.                                    (2)

Then for every b with 0<=b<=1-a, it is IMPOSSIBLE that

    S<=b^2 Q.                                        (3)

In particular

    ||S||op^(1/2) >= (1-a) sqrt(lambda_min(Q)).        (4)

This gives a nonzero distortion floor for the proposed equations,
not merely nonexistence of an exact fixed point.

## 2. Proof using actual relative-height projections

Let V_0=span(A), and recursively

    V_(j+1)=U L2_even(sigma(V_j)).

These are closed nested first-chaos subspaces. Nesting starts because
A is creation-closed and then follows by monotonicity. Their union
is dense in the canonical first chaos: U1 is in V_1, and induction
adds every finite canonical tree. The isometry gives the projection
identity onto V_(j+1), namely U followed by conditional expectation
onto sigma(V_j).

Let Z_j be the projection of Z onto V_j, and put P_j=Cov(Z_j).
Then P_0=0, P_j increases to I, and the Gaussian conditional law of
Z given sigma(V_j) has mean Z_j and covariance I-P_j. The projected
inverse response therefore has covariance K_A(P_j).

Assume inductively that P_j<=Q. Monotonicity and (2) give
K_A(P_j)<=a^2 Q. Projection cannot increase the residual covariance,
so (3) and the scalar L2 triangle inequality in EVERY direction v give

    sqrt(v^T P_(j+1) v)
       <=sqrt(v^T K_A(P_j)v)+sqrt(v^T S v)
       <=(a+b) sqrt(v^T Q v)
       <=sqrt(v^T Q v).

Thus P_(j+1)<=Q. Starting from zero traps all P_j below Q, in
contradiction to P_j tending to I. This proves the theorem. For (4),
an operator-norm residual strictly below its right side would imply
(3) with b<1-a since Q>=lambda_min(Q) I.

No covariance recursion was asserted for Z itself when the residual
is nonzero. The proof uses the projected actual equation (1) and
directional Minkowski, which is the needed quantitative replacement.

## 3. Explicit robust two-cycle obstruction

Use the actual anchor G_0=U1. For fixed 0<r<=1 take the proposed
inverse features

    h_1=sqrt(r) G_0 z_2+sqrt(1-r) h_2(G_0),
    h_2=G_0 sign(z_1),

where h_2(G_0) on the first line denotes the normalized second
Hermite polynomial, not the second component of this vector.
The notation can equivalently be replaced there by
(G_0^2-1)/sqrt(2). These features are orthonormal and anchor-orthogonal.
The exact diagonal kernel is

    K_A(diag(q_1,q_2))
       =diag(1-r+r q_2, (2/pi)arcsin(q_1)).

Choose

    Q=diag(1-3r^2/8, 1-r/2).

The first gap Q_11-K_A(Q)_11 is r^2/8. For the second, use
arccos(1-u)>=sqrt(2u), which follows from 1-cos(t)<=t^2/2:

    Q_22-K_A(Q)_22
       >=(sqrt(3)/pi-1/2)r > r/20 >=r^2/20.

The strict comparison sqrt(3)/pi>11/20 follows, for example, by
squaring the elementary rational bound pi<22/7. Hence

    K_A(Q)<=Q-(r^2/20)I <=(1-r^2/20)Q,
    Q>=I/2.

For ANY actual independent Gaussian innovation frame Z_1,Z_2,
orthogonal to G_0 in the canonical space, its residual in these
proposed equations must satisfy

    ||Cov(Z-Uh(G_0,Z))||op^(1/2)
       >= [1-sqrt(1-r^2/20)]/sqrt(2)
       >= r^2/(40sqrt(2)).                            (5)

The final inequality uses 1-sqrt(1-x)>=x/2. Thus even arbitrarily
large finite or countable canonical representations cannot approximate
this fixed positive-return feedback to zero error. At r=0 the floor
vanishes and the exact feedforward construction exists, as it should.

## 4. Scope of accumulation

An exact subfixed matrix on a critical boundary need not imply a
positive approximation floor. Hypothesis (2) is strictly stronger:
it demands contraction below a positive definite barrier in every
output direction. Critical shifts can have approximate solutions
without exact endogenous solutions. We do not erase that distinction.

This theorem can reject an approximate stationary-response proposal
when its displayed inverse equations have such a barrier. Applying
it to all high-value stationary responses would require deriving
their barrier uniformly; that remains open. It is not a recurrence
or a nonconvergence theorem for M_n/n^(3/2).
