# Affine bent families give signed profile modules and a Toffoli gate

Date: 2026-09-06. Status: an exact extension of the dual-bent-index
realization lemma. Only an AFFINE family of components has to be
bent; the unused linear hyperplane need not consist of bent
components. This gives a concrete nonlinear three-bit index gate
and all unitriangular inversion gates as legitimate R tests. No
closure under overlapping gate composition is asserted.

## 1. Affine component identity

Let n be even. Suppose Boolean functions h,h* on F_2^n and maps
Phi,Psi:F_2^n -> F_2^r satisfy, for every a in F_2^r,

    H_n [h(x)(-1)^(a dot Phi(x))](u)
                     =h*(u)(-1)^(sigma(a) dot Psi(u)),    (1)

where sigma is a permutation of F_2^r. It may be nonlinear and
need not fix zero. Assume every nonzero character of Psi has
expectation of modulus less than one. (The same property for
Phi follows if desired from an injective component family, but
only output equidistribution is needed for the lower test.)

Use s independent copies, defining h_s as the product of their
h values and Phi_s as the sum of their Phi values, and similarly
on the output. Tensoring (1) preserves the same index sigma.
Every nontrivial character expectation of Psi_s tends to zero,
so Psi_s approaches the uniform law on F_2^r.

For an arbitrary scalar profile f, expand its probability Fourier
series. The exact intertwining is

    H_(ns)[h_s f(Phi_s)]
               =h*_s (T_sigma f)(Psi_s),
    T_sigma=H_r P_sigma H_r.                              (2)

There is NO constant-channel defect in (2). The profile's
constant coefficient uses the bent carrier h rather than the
constant physical function. In the language of the earlier
balanced-profile construction, the Fourier support has been
confined to the affine component hyperplane (1,a); characters
on the inactive hyperplane (0,a) never enter the calculation.

For every symmetric finite seed B and arbitrary Boolean profile
array F on F_2^r, equation (2) therefore gives

    R(B)>=(1/2) E_t ||(T_sigma F)(t) B||_1.               (3)

The physical input h_s F(Phi_s) is Boolean regardless of the
distribution of Phi_s. After transformation the common Boolean
carrier h*_s disappears inside the rowwise absolute value.
Uniform convergence of the output label law proves (3).
The physical outer order 2^(ns) is an even-dimensional Walsh
order, hence allowed through the regular R4 bilinear equivalence.

This is a SIGNED profile module. Products of two embedded
profiles lose the common carrier, so (2) does not license the
composition of arbitrary overlapping modules.

## 2. General affine inverse-matrix construction

Let M(a)=M_0+sum_(i=1)^r a_i M_i be an affine family of invertible
d-by-d binary matrices. Suppose another affine family N(b)
satisfies

    M(a)^(-1)=N(sigma(a))                                (4)

for a permutation sigma of F_2^r. Put

    h(x,y)=(-1)^(x^T M_0 y),
    Phi_i(x,y)=x^T M_i y,

and define h*,Psi from N_0,N_i with output coordinates reversed.
Summing first over x in the normalized Walsh transform forces
M(a)y=u and proves

    H_(2d) (-1)^(x^T M(a)y)(u,v)
                          =(-1)^(v^T M(a)^(-1)u).         (5)

Equations (4)--(5) are exactly (1). If the matrices N_1,...,N_r
are linearly independent, then every nonzero character of Psi
is a nonzero bilinear form. For a nonzero binary matrix L,

    E_(u,v) (-1)^(v^T L u)=2^(-rank L)<=1/2.             (6)

Thus the output law in s copies is uniform up to nontrivial
Fourier coefficients at most 2^(-s). This verifies every
equidistribution hypothesis used in Section 1.

## 3. An exact three-bit Toffoli index permutation

Take

             [1  a  c]
    M(a,b,c)=[0  1  b].
             [0  0  1]

Then

    M(a,b,c)^(-1)=M(a,b,c+ab).                            (7)

Hence the active affine dual-index map is the Toffoli gate

    sigma(a,b,c)=(a,b,c+ab).                              (8)

An explicit carrier and label map on F_2^3 x F_2^3 are

    h(x,y)=(-1)^(x dot y),
    Phi(x,y)=(x_1 y_2, x_2 y_3, x_1 y_3).

After swapping the two output coordinate triples, the same
carrier and label map appear. Thus the symplectic Walsh
involution has an invariant signed module with quotient
T_sigma=H_3 P_Toffoli H_3. Equations (2)--(3) apply to every
Boolean profile on the eight labels, including unbalanced ones.

There is no contradiction with the full vectorial quadratic
restriction: the inactive components sum_i a_i Phi_i have
singular polar matrices. They were never required to be bent.

## 4. All unitriangular inversion permutations

More generally take M=I+A, where A ranges over every strictly
upper-triangular d-by-d binary matrix. The r=d(d-1)/2 free
entries give an affine space of invertible matrices, and

    M^(-1)=I+A+A^2+...+A^(d-1)

belongs to the same affine space. Inversion is a permutation
of the r entry coordinates, indeed an involution. Section 2
therefore realizes H_r P_inv H_r as a full arbitrary Boolean
profile lower test for R.

The output entry (i,j) is the parity of all directed increasing
paths from i to j, each weighted by the product of its edge
entries. This supplies high-degree index permutations in one
legitimate construction, rather than by assuming products of
previously realized gates are available.

## 5. Exact remaining gap

Toffoli gates generate reversible circuits when compositions
are allowed, but that algebraic fact alone is irrelevant here:
we have proved individual affine inverse-matrix modules, not
arbitrary products on a shared profile space.

Possible catalytic use of the unitriangular path polynomials
would require a lossless ancilla construction: unwanted free
entry coordinates and computed intermediate values cannot be
postselected at a fixed positive probability cost. In particular,
fixing an arbitrary nonlinear graph of Fourier coordinates does
not automatically leave a Boolean profile after inverse Walsh.

The full index-programming problem, the equality R=T, and the
original convergence question all remain open after this result.

Exact integer replay is provided by
`computations/resumed_convergence_affine_bent_verify_2026_09_06.py`.
It passes for d=2,3,4, checks every component identity, every
profile on the eight-label Toffoli space, 66 profiles on the
64-label d=4 space, both carrier orientations, orthogonality,
and the nontrivial label-character contraction bound 1/2.
