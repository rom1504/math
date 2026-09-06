# Conditional bent blocks do not amplify the selector seed problem

Date: 2026-09-06. Status: exact scoped reduction, in response to the
weighted-CHSH prescribed-spectrum target. This is not an obstruction
to general Maiorana–McFarland superclasses, nonflat conditional blocks,
growing selector problems, or the original R norm.

## 1. The full conditional-bent theorem

For a symmetric k-by-k real seed B and d>=0 define

    C_d(B)=(1/2) max_(f_i Boolean on F_2^d)
                         E_c sum_j |sum_i B_ij (U_d f_i)(c)|,

where U_d is the normalized d-bit Walsh transform.

Let x range over F_2^(2m), z over F_2^d, and suppose each Boolean
input column F_i(x,z) is bent in x for every fixed z. In other words,

    b_i(a,z)=(U_(2m) F_i(.,z))(a) is Boolean for every a,z,i.

The tensor product Fourier identity is exact:

    (U_(2m+d) F_i)(a,c)=(U_d b_i(a,.))(c).                 (1)

For each fixed a the k functions b_i(a,.) are Boolean d-bit inputs.
Taking the optimal unrestricted Boolean output of the seed bilinear
form, then averaging over a, therefore proves

    (1/2) E_(a,c) sum_j |sum_i B_ij U_(2m+d)F_i(a,c)|
        <= C_d(B).                                        (2)

Conversely choose a common bent carrier h(x) and any optimal d-bit
inputs f_i(z), and set F_i(x,z)=h(x)f_i(z). Its transformed carrier
is also Boolean, so equality holds in (2). Hence the supremum over
ALL conditional-bent inputs, for every m, is exactly C_d(B).

This is not a mere bound on peak Fourier coefficients. It is an
exact identity for the full vector-valued, seed-weighted bilinear
optimization with unrestricted output signs. Making the conditional
bent block larger does not improve the selector norm at all.

## 2. Conditional Maiorana–McFarland permutations are included

For each seed coordinate i and selector z, let pi_(i,z) be any
permutation of F_2^m, and g_(i,z) any Boolean exponent. Set

    F_i(x,y,z)=(-1)^[x dot pi_(i,z)(y)+g_(i,z)(y)].

Directly summing over x gives

    U_(2m) F_i(a,b,z)
       =(-1)^[g_(i,z)(pi_(i,z)^(-1)(a))
                              +b dot pi_(i,z)^(-1)(a)].    (3)

This is Boolean for every a,b,z. Thus (2) includes arbitrary
selector-dependent permutations and signs, separately and
independently for every seed coordinate. No common permutation or
common dual function is assumed.

The full transform is a d-bit Walsh transform of the Boolean
z-profile in (3). A coupled prescribed-spectrum solution in this
architecture must already solve the corresponding d-bit coupled
problem; the large Maiorana–McFarland block does not remove it.

## 3. Stable version for approximately bent conditional blocks

Let Z_i(a,z)=U_x F_i(a,z), and suppose Boolean b_i(a,z) satisfy

    E_(a,z) sum_i |Z_i(a,z)-b_i(a,z)|^2 <= epsilon^2.

Walsh orthogonality preserves this error after transforming z.
Cauchy-Schwarz and the operator norm of B then give

    (1/2) E_(a,c) ||(U F)(a,c) B||_1
       <= C_d(B)+(sqrt(k)/2)||B||_op epsilon.              (4)

Thus asymptotically flat conditional carriers also only transport
the selector problem, up to a vanishing error. This includes using
an increasingly good scalar ultraflat carrier as a tensor factor.
It does not block correlated NONFLAT conditional blocks, which are
exactly what the remaining weighted target may require.

## 4. The concrete weighted target and fixed-selector consequences

For B=[-1,2;2,4], equality at T(B)=5/sqrt(2) requires Boolean input
and output profiles satisfying, in normalized L2,

    U F1=(-G1+2G2)/sqrt(2),
    U F2=(G1/2+G2)/sqrt(2).                                (5)

Both Boolean correlations are necessarily 3/4. Equivalently write
F1=h(1-2S), F2=h, G1=k(1-2T), G2=k, with S,T Boolean indicators
of density 1/8. Then (5) becomes

    U h=k(3/2-T)/sqrt(2),
    U(hS)=k(1/2-3T)/(2sqrt(2)).                             (6)

This is a coupled synthesis problem: constructing one scalar
two-level spectrum is not sufficient without the sparse partner
and the common output signs.

For fixed even d, every normalized d-bit Walsh value lies in a
fixed dyadic rational lattice. The target magnitudes for U F2,
1/(2sqrt(2)) and 3/(2sqrt(2)), are irrational, so their distance
from that finite alphabet is positive. Consequently the exact
conditional-bent class with this fixed selector dimension has a
positive polar defect independent of the size of the bent block.

For odd d in {1,3,5,7}, Schmidt's primary manuscript records the
known minimum peak normalized Walsh value sqrt(2). This exceeds
the larger target magnitude 3/(2sqrt(2)) by 1/(2sqrt(2)). Every
Boolean selector section therefore has at least one frequency
with squared amplitude error at least 1/8, and its mean squared
F2-target defect is at least 1/(8*2^d). Again this is a fixed-selector
statement, not a global R gap.

The source for that last small-d peak fact is Schmidt,
*Asymptotically optimal Boolean functions*,
[author manuscript](https://math.uni-paderborn.de/fileadmin-eim/mathematik/AG-Diskrete_Mathematik/Publications-schmidt/pw.pdf),
Introduction, immediately after the definition of mu_n. The main
conditional-bent identities (1)–(4) are elementary derivations here
and do not rely on a literature classification theorem.
