# Independent audit of the discontinuous rank-deficient feedback obstruction

Date: 2026-09-06. I reconstructed every step of
`resumed_response_matrix_height_feedback_obstruction_2026_09_06.md`.
The stated theorem and the explicit two-innovation obstruction pass.
This audit is not a proof that every full-response stationary threshold
has the return-path property required by that example.

## 1. Creation and height must live in the same canonical space

The key identity is exact:

    Pi_(h+1) U f=U E[f|F_h].

Conditional expectation keeps exactly the normalized even Hermite
monomials involving only children of height at most h. Their created
parent has height at most h+1. Finite sums prove the identity, and
Gaussian L2 continuity proves it for all even L2 inputs. The edge
U1 is at height one. A fixed height contains infinitely many canonical
coordinates; the proof correctly does not replace height by finite
cardinality or degree.

For an orthonormal Gaussian FIRST-CHAOS frame G in that same space,
the conditional mean given F_h is Pi_h G and its conditional covariance
is the deterministic matrix I-P_h. Therefore the Gram of the conditional
expectations of g(G) is exactly K(P_h), and isometry gives
P_(h+1)=K(P_h). K is an uncentered second-moment kernel: replacing it
by centered covariance would incorrectly remove the root/constant
creation term K(0)=Eg(Eg)^T.

Every individual tree has finite height. Hence the increasing height
projections exhaust the entire canonical first-chaos Hilbert space,
and P_h tends to I for a finite orthonormal frame. A distributional
Gaussian solution in an enlarged probability space is not enough to
invoke this exhaustion; the source note explicitly requires endogeneity
in the original creation space.

## 2. Matrix monotonicity and anchored comparison

For arbitrary 0<=P<=Q<=I, including noncommuting P,Q, take independent
Gaussian pieces with covariance P,Q-P,I-Q. Conditioning g of their sum
on the first piece or on the first two pieces gives K(P),K(Q) as the
corresponding second-moment Grams. Conditional variance makes their
difference PSD. No simultaneous diagonalization is assumed.

Thus any proper P with K(P)<=P traps the exact recursion started at
zero, contradicting its limit I. A rank-deficient defect is sufficient;
one does not need a positive-definite gap in I-P or P-K(P).

At an anchored comparison matrix diag(I,Q), the two anchor copies are
identical. Their inverse features are functions of those anchors alone.
The cross blocks vanish by the stipulated orthogonality of the innovation
inverse features to every anchor inverse feature. The anchor block is
I and the remaining block is K_A(Q). This is a valid full-matrix
comparison from P_0=0, even before the actual height process has
finished revealing the anchors.

The local polynomial derivative L(D)=E[J_h D J_h^T] is also correctly
normalized. Gaussian covariance differentiation sums every ordered
matrix entry of D. A strict positive-definite expansion L(D)>=(1+eta)D
gives a subfixed I-epsilon D after controlling the polynomial Taylor
remainder. It would be invalid to use a single large derivative entry
as a replacement for this matrix inequality; the note does not do so.

## 3. Exact two-coordinate example, including the endpoint

For independent G0,Z1,Z2, the proposed inverse features

    h1=sqrt(r) G0 Z2+sqrt(1-r) h2(G0),
    h2=G0 sign(Z1)

are jointly even, norm one, mutually orthogonal, and orthogonal to 1.
With diagonal innovation cross covariance diag(q1,q2), direct expectation
gives diag(1-r+r q2,(2/pi)arcsin(q1)); the off-diagonal terms vanish
by independence of the two innovation pairs and by the odd G0 factors.

At Q=diag(1-r^2/2,1-r/2), the first entry of K_A(Q) equals Q exactly.
The second is strictly smaller for every 0<r<=1, since

    1-cos(pi r/4)<=pi^2 r^2/32<r^2/2.

The direction of the inverse-cosine inequality is correct, so
(2/pi)arcsin(1-r^2/2)<1-r/2. This produces the required subfixed matrix
even though only one entry has a strict defect. The exclusion is exact
for each positive r; it is not an interchange of r and matrix-size
limits or an inference from infinite Sobolev energy alone.

At r=0 there is an actual causal construction: first
Z1=U h2(G0), then Z2=U[G0 sign(Z1)]. Its second inverse is orthogonal
to 1 and h2(G0), and has norm one. Therefore Z2 is a standard Gaussian
independent of the previous frame. It lies in height three, despite
using infinitely many Hermite degrees. This confirms that discontinuity
alone is harmless in a feedforward addition; the smooth return edge
at r>0 is the crucial extra constraint.

## 4. Precise remaining obligation

The theorem supplies a test for a PROPOSED endogenous feedback frame.
It does not show that an arbitrary high-value stationary response
necessarily has a nonzero return path through every discontinuous
direction, nor does it cover every rank-deficient jump by fiat. Any use
of this obstruction to exclude a general stationary fixed point must
first construct its actual subfixed covariance matrix, or prove the
claimed return-path reduction. That necessity remains open in the
source artifact and in this independent audit.

## 5. Subsequent one-rough-own-coordinate theorem also passes

I separately reconstructed the added Section 7. For Q=I-epsilon e_d e_d^T,
partial OU expansion makes E=I-K(Q) exactly the Gram with multiplier
1-(1-epsilon)^ell. The smooth-complement generator-domain assumption
sum ell^2||h_i,ell||^2<infinity gives E_od=O(epsilon) by moving the
self-adjoint OU defect onto h_i. Merely finite first derivative energy
would not give this estimate; the stronger hypothesis is correctly
stated.

On the range of A=lim E_oo/epsilon, its inverse is O(1/epsilon), so
the cross Schur correction is O(epsilon). Every vector in ker A has
all positive own-coordinate Hermite coefficients exactly zero; its
rows in E, including the rough cross column, vanish exactly, not just
to first order. Meanwhile infinite own-coordinate Dirichlet energy
gives E_dd/epsilon tending to infinity. The Schur complement therefore
proves E>=epsilon e_d e_d^T and traps the height recursion below the
proper PSD subfixed matrix Q. This handles the singular smooth block.

The conclusion is valid with Q<=I and Q!=I; the gap I-Q is rank one,
not positive definite. Any displayed Q<I here should be understood as
proper inequality, not strict Loewner positive-definiteness. Feedforward
W=U H(V) has zero energy in its own W coordinate and is not excluded.
If several inverse features are rough in the same coordinate, the
smooth-complement hypothesis may fail and no conclusion follows.
