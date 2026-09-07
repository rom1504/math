# The same-positive-ground cosh kernel need not be positive semidefinite

2026-09-07. Exact finite counterexample on an ACTUAL minimizer, followed by a scalable leading-scale counterexample on ACTUAL bounded-cap full signings. The scalable construction is deliberately not claimed to be near-minimal or exactly minimal.

For a hollow symmetric signing A put P=max_x H_A(x)>0. The tested kernel is

    K_t(x,y)=cosh(t x^T A y/(2P)),

restricted to Boolean x,y with H_A(x)=H_A(y)=P. The proposed positive-semidefiniteness is false, despite the correct two-point inequality |x^T A y|<=2P.

## 1. Exact optimal order-five counterexample, for every t>0

Let A5 have coefficient -1 on the edges of the five-cycle and +1 on its chords, with diagonal zero. Its Boolean energies are exactly 0,+4,-4. To check this, negate the whole spin if necessary and classify its at most two negative coordinates: no negatives or one negative gives zero; two adjacent negatives give -4; two nonadjacent negatives give +4. Thus Q(A5)=P=4.

It is an actual exact minimizer: every order-five full signing has E_x H_A(x)^2=10, hence Q(A)>=sqrt(10), and its energies are even integers. Therefore M5>=4, attained by A5.

The three positive ground states

    x=(-1,1,-1,1,1),
    y=(-1,1,1,-1,1),
    z=(-1,1,-1,-1,1)

have exact Gram matrix

    [8 4 8]
    [4 8 8].
    [8 8 8]

Writing a=cosh(t) and b=cosh(t/2), their kernel matrix is

    [a b a]
    [b a a].
    [a a a]

The vector (1,1,-2) has quadratic value 2(b-a)<0 for every t>0. Equivalently the determinant is -a(a-b)^2<0. This is an exact three-point failure, not a numerical eigenvalue observation.

## 2. Uniform failure along the whole continuous ground set

On [-1,1]^5 the maximum of H_A5 remains 4. Every continuous maximizer lies on a Boolean positive-ground vertex or a positive-ground edge. Indeed multilinear interpolation expresses H(z) as the mean over the smallest Boolean face containing z. Equality at the maximum forces every vertex of that face to maximize. A face with two free coordinates is impossible, since its nonzero mixed quadratic coefficient A_ij would have to vanish.

The ten positive-ground vertices and ten positive-ground edges are obtained by cyclic rotations and global negation from

    z(s)=(-1,1,-1,s,1),  -1<=s<=1.

This small assertion is also checked exhaustively by the companion exact script. Let R be cyclic rotation. For this family the five cyclic Gram entries are

    z(s)^T A5 R^j z(s)
       =(8, -(7+s^2), 3+s^2, 3+s^2, -(7+s^2)).      (1)

At t=1, the corresponding five-point cosh kernel has a circulant eigenvalue

    lambda(s)=cosh(1)-phi cosh((7+s^2)/8)
                       +phi^(-1)cosh((3+s^2)/8),     (2)

where phi=(1+sqrt(5))/2. It decreases as s^2 increases: the negative derivative term has both the larger coefficient and larger positive sinh argument. Hence

    lambda(s)<=lambda(0)
       =cosh(1)-phi cosh(7/8)+phi^(-1)cosh(3/8)
       <-.0728<-.07.                                (3)

The last rational comparison is checked with outward interval arithmetic. In particular the failure is uniform over the ENTIRE continuous ground set, not only over the three Boolean corners from Section 1. Continuity gives a fixed neighborhood of that compact set in which a negative eigenvalue bounded away from zero persists.

## 3. Actual full-sign amplification at bounded normalized cap

There exist fixed constants C0<infinity and gamma>0 and actual full signings A_n, at every sufficiently large n divisible by five, such that

    Q(A_n)<=C0 n^(3/2),

and FIVE ACTUAL positive ground states have a normalized cosh kernel K_1 with smallest eigenvalue at most -gamma. Thus a vanishing-error positive-semidefinite repair of that normalized ground kernel is impossible uniformly over bounded-cap actual inputs.

Here is the construction. Write n=5m. Fix a sufficiently large constant L, independent of n, and define the hollow mean matrix

    M_n=(L/sqrt(n))(A5 tensor J_m).

Within each macro block the mean is zero, including its diagonal. For n>=L^2 all off-diagonal means lie in [-1,1]. The cyclic macro-block shift acts freely on unordered edges, with every orbit of size five. Independently choose one sign per edge orbit with the prescribed mean, and assign it to all five edges in that orbit. The resulting A_n is an ACTUAL symmetric hollow full signing, invariant under the cyclic macro-block shift.

There is a realization with

    Q(E_n)<=3n^(3/2),   E_n=A_n-M_n.                  (4)

For a fixed spin the error energy is a sum of independent centered orbit signs with coefficients b_O of absolute value at most five. The squared coefficient sum is at most 5 binom(n,2). Hoeffding followed by the union over all 2^n spins bounds the failure probability in (4) by

    2 exp[(log2-9/5)n]<1

for large n. Thus (4) is a genuine deterministic existence statement and preserves the cyclic invariance.

For a Boolean x let z_i be its average spin in macro block i. Then EXACTLY

    H_Mn(x)=(L/25)n^(3/2) H_A5(z),
    P(M_n)=Q(M_n)=(4L/25)n^(3/2).                    (5)

Choose any ACTUAL positive ground state x of A_n, with energy P_n. From (4)-(5),

    |P_n-(4L/25)n^(3/2)|<=3n^(3/2),
    H_A5(z)>=4-150/L.                               (6)

Compactness of [-1,1]^5 and the exact maximum-set description imply that every such z is uniformly close to the ground set of Section 2 when the FIXED constant L is sufficiently large. This choice does not depend on m.

The five cyclic shifts R^j x are actual positive ground states of A_n, because A_n was built invariant under that permutation. Their mean-matrix cross energies satisfy

    (R^r x)^T M_n(R^s x)
       =(L/25)n^(3/2)(R^r z)^T A5(R^s z).

The error cross energies are uniformly at most 4Q(E_n)<=12n^(3/2) by hollow polarization. Together with (6), division by 2P_n shows that their normalized Gram entries differ from

    (R^r z)^T A5(R^s z)/8

by O(1/L), uniformly in n and in the choice of the actual ground state. First use the uniform margin in (3) and continuity to choose a neighborhood, then choose the fixed L sufficiently large for (6) and the O(1/L) error to fit inside it. The actual five-point cosh kernel then has, for example, an eigenvalue at most -gamma for some fixed gamma>0. Also

    Q(A_n)<=(4L/25+3)n^(3/2),

so the inputs have a fixed bounded normalized cap, as claimed.

## 4. Exact scope

Section 1 disproves the universal finite-order PSD assertion even for exact minimizers. Section 3 disproves uniform asymptotic PSD with an o(1) normalized kernel repair over ALL bounded-cap actual full signings. It does not put that scalable family near the optimum: its cap constant is deliberately large because the mean signal dominates the orbit-rounding noise.

The more selective possibility that suitable asymptotic exact minimizers have approximately positive ground kernels is not disproved here. Neither is a seed-transfer construction that does not require this kernel's positive semidefiniteness. The result is an exact test of the proposed ground-set repair, not a no-go theorem for global flatification.
