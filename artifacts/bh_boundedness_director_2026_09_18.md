# Director derivations: Boolean BH boundedness

Status: separate three-hour campaign, 2026-09-18. The signing problem is
paused. Results below concern the actual scalar Boolean supremum, not
a completely bounded or torus norm. External novelty is not asserted.

## 1. Uniform control for open and closed Hadamard chains

**Proved; independently reconstructed by the barrier researcher.**
Let d>=2 and let H_1,...,H_(d-1) be real or complex Hadamard matrices
of common order N. Thus U_j=H_j/sqrt(N) is unitary and every entry of
U_j has modulus N^(-1/2). Give each slot its own N Boolean variables.
The normalized open-chain polynomial is

    T(x_1,...,x_d)
      =N^(-1) x_1^T U_1 D(x_2) U_2 ... D(x_(d-1)) U_(d-1) x_d.

Its N^d distinct coefficients have modulus N^(-(d+1)/2). Consequently
its coefficient q_d norm is exactly one. Its scalar Boolean supremum is
at least 1/sqrt(2) over the real field, and at least sqrt(2)/pi over
the complex field. Thus its BH ratio is at most sqrt(2), respectively
pi/sqrt(2), uniformly in BOTH d and N.

Proof: for d>=3 fix the variables from x_3 onward and write the tail
as v, with ||v||_2=sqrt(N). Randomize x_2. Every coordinate of
U_1 D(x_2)v has squared Rademacher coefficient norm one. Khintchine's
L1 inequality gives expected coordinate modulus at least 1/sqrt(2).
In the real case choose x_1 to align all coordinates. In the complex
case use the elementary inequality

    max_(epsilon_i=+-1) |sum_i epsilon_i z_i|
       >=(2/pi) sum_i |z_i|.

Indeed maximize the real part after a global phase, and average that
phase. This proves the stated bound after division by N. For d=2
randomize x_d directly. The complex Rademacher L1 constant 1/sqrt(2)
follows from its real version: average real projections, then minimize
the concave average sqrt(lambda_1 cos^2 t+lambda_2 sin^2 t) at a
rank-one covariance with fixed trace. No tensor polarization is used.

The same conclusion holds for the cyclic polynomial

    C(x_1,...,x_d)=Tr(D(x_1)U_1 ... D(x_d)U_d),

where all d edges are normalized Hadamards. Its coefficient q_d norm
is sqrt(N). Fix x_3,...,x_d. The remaining bilinear coefficient matrix
has entries (U_1)_(ij) V_(ji), with V unitary. Every row has squared
norm 1/N. Random x_2 and the same sign/phase argument give Boolean cap
at least sqrt(N/2), respectively sqrt(2N)/pi. This includes the closed
cycle; closing a path does not create a degree-growing BH ratio here.

Scope: slots are disjoint, all edges are square of the same size, and
each has flat unitary entries. Identifying variables can create Walsh
collisions and degree collapse and is NOT covered by this proof.

Literature connection: the primary September 14 paper
[Hadamard Rigidity and Sharp Stability in the Completely Bounded
Bohnenblust--Hille Inequality](https://arxiv.org/abs/2609.16329)
identifies Hadamard chains as its exact extremizers. Our scalar estimate
above is elementary and does not depend on that classification proof.
Its fixed-degree stability theorem is in coefficient norm, with
degree-dependent constants. Neither fact supplies uniform scalar BH
control for arbitrary polynomials or even a dimension-free transfer of
coefficient perturbations to the scalar supremum. No such transfer is
claimed.

## 2. A useful exact diagnostic, not a solution

For nonzero f define pi_S=|fhat(S)|^2/||f||_2^2. For every positive
degree budget m the exact identity is

    log R_m(f)=log(||f||_2/||f||_infinity)
                  +H_(m/(m+1))(pi)/(2m).

Hence support at most exp(Cm), or ambient dimension at most Cm,
immediately gives a uniform ratio. A divergent family must escape these
simple support bounds and must compensate for its L2-to-supremum loss.
This identity alone is merely a reformulation, not a lower or upper
bound for the unrestricted constants.

The counterexample track now has a finite degree-four example with
ratio strictly above two. It escapes torsion-valued outputs through a
tangent perturbation and actual inter-block degree transfer. Ordinary
tensor amplification cannot turn this into divergence. The next test
is whether degree transfer can operate in a growing ambient/degree
geometry, such as balanced address lifts, without losing its cap or
spending the degree saved in normalization. See the counterexample
artifact for the exact coefficients and current evidentiary status.
