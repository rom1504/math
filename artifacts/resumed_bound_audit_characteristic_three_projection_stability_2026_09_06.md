# Weighted characteristic-three bound and linear edit stability by projection

Date: 2026-09-06. The weighted extension and projection argument were
independently reconstructed by the convergence agent and director.
This is an actual-signing structural corollary. It imposes no spectral
flatness or operator-norm bound.

Let q=3^r>81 and define

    rho_q=(q-29)/(26(q-3)),
    gamma_q=(q-81)/(81(q-1)),
    c_(3,q)=(4/9)/sqrt(1-rho_q/10+rho_q^2/(400 gamma_q)).

Then c_(3,q) tends to

    c_3=2080/(9 sqrt(269441))
       =0.44523467985944279335138857059...

The input theorem and its full independent audit are
`resumed_convergence_characteristic_three_incidence_gain_2026_09_06.md`
and `resumed_bound_audit_characteristic_three_incidence_2026_09_06.md`.

## 1. The same incidence proof applies to weighted Cayley matrices

For every hollow real symmetric additive-Cayley matrix C on F_3^r,

    Q(C) >= c_(3,q) q/sqrt(q-1) ||C||_F.                     (1)

There is no requirement that the off-diagonal entries be signs or even
bounded by one. To verify (1), repeat the proved incidence argument
with K=2Q(C)/q^(3/2), t=lambda_C(0)/sqrt(q), and
h(a)=(8lambda_C(a)+lambda_C(0))/(9sqrt(q)). The three- and nine-point
Boolean profile inequalities, projective color counts, and spectral
deficit direction are unchanged for real weights.

Hollowness still gives sum lambda_C=0. The only changed identity is

    E h^2=(64/81) v^2+gamma_q t^2,
    v^2=||C||_F^2/[q(q-1)],

because sum lambda_C^2=||C||_F^2. All subsequent inequalities are
homogeneous. The same completion of the square therefore gives
Q(C)/q^(3/2)>=c_(3,q) v, exactly (1). For C=0 this is trivial.

## 2. Translation averaging is a norm-decreasing orthogonal projection

Let B be ANY hollow real symmetric matrix indexed by F_3^r. For each
translation z let T_z be its permutation matrix, and put

    P B=(1/q) sum_z T_z^T B T_z.

This is a hollow real symmetric Cayley matrix. P is the orthogonal
projection in the Frobenius inner product onto the translation-invariant
matrix subspace: it is self-adjoint and idempotent by the finite-group
average. The cap Q is convex and invariant under vertex permutations,
so

    Q(P B)<=Q(B).                                             (2)

Fix ANY hollow Cayley signing A on this same group. Then P A=A and
||A||_F=sqrt(q(q-1)). Orthogonality and Cauchy--Schwarz imply

    |<A,B>_F|=|<A,P B>_F|
                       <=sqrt(q(q-1)) ||P B||_F.

Combining with (1)--(2) proves the linear correlation bound

    Q(B)/q^(3/2)
        >=c_(3,q) |<A,B>_F|/[q(q-1)].                        (3)

Thus (3) holds for arbitrary hollow real B, not merely for a Cayley B
or for another signing. The reference A only needs to be a Cayley
signing. This is a convex-projection lower test, not an assertion that
the oriented ground-state covariance of A is radial.

## 3. Exact edit bound and its asymptotic consequences

If B is a signing differing from A on d unordered off-diagonal entries,
then every changed edge changes two Frobenius summands from +1 to -1.
Therefore

    <A,B>_F=q(q-1)-4d,

and (3) becomes

    Q(B)/q^(3/2)
      >=c_(3,q) |1-4d/[q(q-1)]|.                            (4)

In particular d=o(q^2) preserves the lower limit c_3. If, for a fixed
0<delta<c_3, a sequence instead obeys
Q(B_q)/q^(3/2)<=c_3-delta+o(1), then (4) implies distance at least

    (delta/(4c_3)-o(1)) q^2

from every characteristic-three Cayley signing A_q and from -A_q.
The statement also holds for all vertex switchings and relabelings of
such reference matrices, by conjugating B and A together; cap and edit
distance are invariant under this operation.

For comparison, the direct variance calculation gives the weaker
bound ||P B||_F^2>=q(q-1)-8d: if m_z directed entries of displacement z
are edited, the averaged kernel is A(z)(1-2m_z/q), and sum_z m_z=2d.
The correlation argument retains the stronger absolute linear factor
in (4), including its symmetry between A and -A.

The result does not say that arbitrary minimizing signings are close
to this class, nor that the minimum coefficient is below c_3. It
provides a uniform obstruction to low-cap sequences being so close.
