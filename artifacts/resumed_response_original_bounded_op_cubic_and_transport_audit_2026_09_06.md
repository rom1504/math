# Original bounded-operator signings: cubic self-energy and the exact new transport gap

Date: 2026-09-06. Status: the bound-audit agent independently reconstructed
the cubic self-energy, endpoint Frobenius bound, and own-spin argument.
This note concerns
arbitrary symmetric hollow signings with a fixed normalized operator cap,
not just involutions. It does not assert a universal mixed-charge gain.

## 1. Setup and a controlled old-frame self-energy

Let B=A/sqrt(n-1), with A symmetric hollow and all off-diagonal entries
in {+1,-1}, and ||B||op<=L fixed. Put Q=B^2 and

    m_3(B)=Tr(B^3)/n.

Exactly, Q_ii=1, |Q_ij|<=1, and Tr(Q^2)<=L^2 n. In particular
|m_3(B)|<=L. Let X be a fixed finite canonical old-tree family. For a
bounded odd Gaussian-a.e.-continuous old function F, write

    r^2=||P1F||_2^2, K_F=U^*F,
    Q_B(v)=v^T Bv/2.

Then

    E Q_B(F(X))/n = (r^2/2)m_3(B)+o(1).         (1)

The error tends to zero for each fixed F and L; no limit of m_3(B) is
assumed. Formula (1) is an EXPECTATION theorem. No new energy-variance
claim is inferred here.

### Proof for fixed polynomials

Write F=L_F+R, where L_F=P1F and R has no first local Hermite chaos.
The old covariance theorem gives

    Cov(L_F)=r^2 Q+o_op(1).

Consequently E Q_B(L_F)/n=(r^2/2)m_3(B)+o(1).

For the nonlinear piece, the full covariance theorem gives

    Cov(R)=sum_(odd k>=3) w_k Q^(circ k)+E_n,
    ||E_n||_*/n ->0.

The trace error is negligible because ||B||op<=L. For every k>=3,

    |Tr[B Q^(circ k)]|/n
      <=(1/sqrt(n-1)) (1/n) sum_(i!=j)|Q_ij|^k
      <=L^2/sqrt(n-1) ->0.                     (2)

Here B has zero diagonal. Thus E Q_B(R)/n=o(1).

The potentially troublesome first/nonlinear covariance does NOT require
an unproved normalized-nuclear bound on Cov(L_F,R). The endpoint lemma
already gives, for every old tree coordinate T,

    ||Cov(BR,X_T)||_F=O(1)

for the exact finite forest main. Therefore

    |E X_T^T BR|/n <= ||Cov(BR,X_T)||_F/sqrt(n)=o(1).

Raw local-Wick and collision errors are then removed by normalized-L2
Cauchy--Schwarz and the fixed operator bound. Summing the finitely many
linear coefficients proves E L_F^T BR/n=o(1), completing (1).

The invoked estimates are exactly Sections 1--7 and 10 of
`resumed_bound_audit_full_nonlinear_covariance_trace_2026_09_06.md` and
Sections 1--6 of `fresh_limit_injective_input_gram_2026_09_05.md`.
The forbidden one-branch/nonlinear star is not incorrectly claimed to
have a small covariance operator; the endpoint trace estimate avoids it.

### Bounded F

Approximate F in old Gaussian L2 by odd polynomials P. The energy
difference is bounded by

    |E Q_B(F)-E Q_B(P)|/n
      <= (L/2) (E||F-P||^2/n)^(1/2)
                    [(E||F||^2/n)^(1/2)+(E||P||^2/n)^(1/2)].

The old local law supplies the fixed-P limit. Also
||P1P||_2^2->||P1F||_2^2 and |m_3(B)|<=L. First take the matrix
limit, then the approximation. This proves (1) for bounded F without
claiming a covariance-operator limit for the raw nonlinear response.

## 2. The own-spin center and its creation field are already controlled

Let H be a bounded even Gaussian-a.e.-continuous old function, or an
explicitly ordered-L2 realization of such approximants. Then

    E Q_B(S H(X))/n ->0.                        (3)

For finite even Hermite polynomials, replace S H(X) in normalized L2
by the corresponding exact injective own-spin inputs. Their covariance
is ||H||_2^2 I+o_op(1), so the normalized energy trace vanishes since
Tr B=0. Bounded H follows by the same fixed-L approximation argument.
Only normalized L2 replacement is needed here; no illicit raw-product
operator-norm substitution is made. A completely arbitrary Borel
representative is not evaluated literally on finite-n field lattices:
Gaussian-null modifications there could invalidate the approximation.

Also the actual input/output relation gives

    B[S H(X)] = (U H)(X)+o_(normalized L2),      (4)

with finite Hermite approximation before the matrix limit and the
Gaussian L2 limit afterward when UH uses infinitely many coordinates.
Equation (4) holds in this general bounded-op signing class. It does
NOT require B^2=I. It follows from X_T=B V_T plus the controlled root
collision, after the own-input local-Wick approximation.

The marked cross identity is

    E[F^T B(SH)]/n -> E[H K_F]=:j.              (5)

Thus, if H>=0 and |F|+H<=1, the two initial feasible endpoints obey

    E Q_B(F+SH)/n = (r^2/2)m_3(B)+j+o(1),
    E Q_B(-F+SH)/n = (r^2/2)m_3(B)-j+o(1).      (6)

More usefully, the exact endpoint identity and Jensen give the
original-matrix conditional lower bound

    max_(x Boolean) |Q_B(x)|/n
      >= |j|+(r^2/2)|m_3(B)|-o(1).             (7)

The cube-to-Boolean step is valid because B is hollow. This is a
genuine extra original-class comparison term. It does not imply that
small m_3 forces involution structure, or improve the banked universal
constant when m_3 may vanish.

## 3. The first-chaos return is Q times the marked input, not that input

The correct decomposition of the next field is

    BF = Q[S K_F(X)] + Z_R +o_(normalized L2),
    Z_R=B(F-P1F)(X).                           (8)

For polynomial F this follows directly by multiplying the old creation
relation by B; bounded F is treated in the stated ordered limits. For
the elementary single-edge choice X=BS and P1F=bX, the coherent term
is exactly b Q S.

Only on an operator-near involution may Q[S K_F] be replaced by
S K_F. The existing full-noise identities allow old even tests times
odd functions of Z_R, and own-spin old even tests times even functions
of Z_R, with the exceptional transport entering LINEARLY. They do not
give a general conditional law of Q[S K_F] under Boolean input, nor
allow its nonlinear thresholding for free.

There is a concrete original-class witness to this distinction. Take
H_m a symmetric Sylvester Hadamard with trace zero, and hollow

    A_(2m)=J_2 tensor H_m - diag(J_2 tensor H_m).

Then B=A/sqrt(2m-1) is operator-close to
B_0=(J_2/sqrt2) tensor(H_m/sqrt m), so ||B||op->sqrt2 and
m_3(B)->0. But

    B_0^2=J_2 tensor I_m, ||B_0^2-I||op=1.

Its coherent return at paired coordinates is S_i+S_(i+m), not S_i.
This is a bounded-op hollow-signing example, not a proposed low-energy
minimizing sequence. It disproves the inference from bounded operator
norm, or from vanishing cubic moment, to the involution return law.

## 4. Exact next mixed-charge term not supplied by current identities

For the physical first endpoint u_1=F+SH, define the exact fields

    A_new=B(SH)+alpha F,
    T_new=BF+alpha SH.

Coupling the computations at B and -B with the same seeds, the next
threshold scores are A_new+T_new and -A_new+T_new. Away from ties their
charge halves are

    D_2=sign(A_new) 1{|A_new|>|T_new|},
    C_2=sign(T_new) 1{|T_new|>|A_new|}.           (9)

The same decomposition can always be made using fixed softsigns, so
no unverified tie assumption is needed to identify the transport issue.
Their two endpoint energies have common term Q_B(D_2)+Q_B(C_2) and
opposite cross term

    E[D_2^T B C_2]/n
      =(1/n) sum_ij B_ij E[
          sign(A_new,i) 1{|A_new,i|>|T_new,i|}
          sign(T_new,j) 1{|T_new,j|>|A_new,j|}].  (10)

Equation (4) identifies the old creation part of A_new. It does not
identify (10). The masks in (9) depend nonlinearly on the exceptional
return Q[S K_F] and on the already transported nonlinear noise Z_R.
They are not old-frame functions to which the existing forest covariance
or own-input Gram theorem may simply be reapplied.

Even a new proof of the full ONE-ROOT joint law of
(old fields, Q[S K_F], Z_R) would not alone identify (10): a cross-root
transport/covariance statement for the new mixed functions is needed.
The current full-noise restricted tests only license BF as a linear
tested multiplier; (10) instead contains B applied to a nonlinear
function of that field. Existing equations do not close that operation.

Equivalently, starting from the stronger terminal best response
C=H sign(BF), its cross energy with F has the banked full-J lower bound
and its F self-energy is controlled by (1), but E Q_B(C)/n is not
identified by these old-frame modules. A further damped update of
F+C requires the field

    B(F+C)=BF+B[H sign(BF)],                    (11)

and the joint transport of the second term in (11) is genuinely new.
The operator bound gives its norm bound, not an independent Gaussian
innovation or a positive uniform stability gap.

This does not establish impossibility of a universal mixed-charge gain.
One could seek new transport identities or inequalities avoiding the
full law. It states the precise closure that remains unproved. In
particular one need not identify every new self-energy to obtain a
lower bound: a uniform positive cross-energy bound in (10) alone would
suffice via the exact endpoint identity. What is not justified is to
import the involution paired-query law into (10) on arbitrary signings.
