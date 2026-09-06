# Resonant carriers restore the full seven-point phase action

Date: 2026-09-06. Status: an exact profile-intertwining calculation and
limiting lower theorem for the prescribed R norm. This strengthens the
mean-zero/contraction phase tests. It is a signed spin module, not a
unital invariant probability algebra, and does not prove composition
of overlapping gates or R=T.

The character/Gauss input and all normalizations are those in
`resumed_convergence_gauss_phase_lower_tests_2026_09_06.md`, which reads
Schmidt's primary Lemmas 4–7. No further number-theoretic theorem is
imported here, and no unknown higher-layer orientation is assumed.

## 1. Full seven-point family

On functions on Z_7 with uniform measure put

    A_ab=chi_7(a-b)/sqrt(7),
    R_theta=cos(theta)(I-E_7)+sin(theta)A,
    V_theta=E_7+R_theta.

Since A is real skew-symmetric and A^2=-(I-E_7), V_theta is a real
orthogonal operator, fixes constants, and obeys

    V_theta=exp(theta A).

Equivalently, after composing with inversion P_7, it is the symmetric
orthogonal involution

    V_theta P_7=E_7+cos(theta)(P_7-E_7)+sin(theta)A P_7.     (1)

The sign of A depends only on the harmless global phase convention.

**Theorem.** For every real symmetric finite seed B and every real
theta,

    R(B) >= (1/14) max_(F,G in {+/-1}^(7*k))
                         |Tr[G^T V_theta F B]|
          = (1/14) max_F ||V_theta F B||_1.                (2)

All finitely many vector columns are transported by the same
construction. Tensoring gives every finite product of the full
operators V_theta as a legitimate limiting bilinear test. These tests
can also be tensored with the independently proved probability-space
averaging reflections.

## 2. The exact digit split

Fix a resonant angle

    theta=2 pi l/7^r,       r>=1,

and then let e>r. Put N=7^e and M=7^(e-1). Write each spatial index
uniquely as

    a=b+tM,       0<=b<M,       t in Z_7.                  (3)

Let D_(e,theta) be the circulant version of the Gauss limiting
operator: compose the symmetric Hankel operator K_(e,theta) with
inversion on one side. Its nonconstant Fourier multiplier is

    exp(i eta_d chi_7(u) 7^(e-d) theta),
    j=7^(e-d)u, 7 does not divide u,

and its constant multiplier is zero. Choose the overall phase sign
so the primitive layer has the orientation defining A above. The
other eta_d remain arbitrary fixed signs.

Averaging a function over t in (3) retains exactly the Fourier
frequencies j divisible by 7. Thus the space of functions with zero
average on each seven-point fiber is exactly the primitive-frequency
space. On that space the multiplier depends only on j modulo 7, so
the action is fiberwise R_theta.

More explicitly, summing the primitive Fourier modes shows its
convolution kernel is supported on spatial differences divisible by
M; its seven values are the entries of R_theta. This proves the
fiberwise assertion without treating cyclic digits as independent
group coordinates.

Choose a Boolean function h of period L=7^(e-r) on Z_N. Because L
divides M, h is constant on each fiber in (3), and can be denoted
h(b). Its Fourier support consists of multiples of 7^r. Every
nonzero multiplier on that support equals one:

    exp(i eta_d chi_7(u) 7^(e-d) 2 pi l/7^r)=1

whenever e-d>=r. This is independent of all eta_d. Consequently

    D_(e,theta) h=h-mean(h).                               (4)

For any seven-point vector profile f, define the N-point profile

    (J_h f)(b+tM)=h(b) f(t).

Its fiber-centered part h(b)(f(t)-mean(f)) has primitive Fourier
support; its fiber mean is mean(f)h(b) and has Fourier support
divisible by 7^r. There are no intermediate frequency layers.
The two exact actions therefore give the full intertwining identity

    D_(e,theta) J_h f
       =J_h V_theta f-mean(h)mean(f) 1.                   (5)

The map J_h preserves L2 inner products and sends every Boolean
profile to a Boolean profile. It is not unital: J_h(1)=h.

## 3. Exact objective and vanishing defect

For any two finite vector profiles f,g and any fixed seed B, (5)
gives

    (1/N) Tr[(J_h g)^T D_(e,theta) J_h f B]
       =(1/7) Tr[g^T V_theta f B]
          -mean(h)^2 mean(g)^T B mean(f).                (6)

Choose h on one period of length L with (L+1)/2 plus signs and
(L-1)/2 minus signs. Then mean(h)=1/L. As e tends to infinity with
r fixed, the defect in (6) tends to zero quadratically in 1/L.

For each e, Schmidt's even-extension subsequences realize the
physical limiting Gauss action with its actual layer orientations.
Resonance makes those orientations irrelevant in (4), and the
primitive sign can be adjusted by theta's overall sign. Each fixed
e construction is a valid bilinear lower test by the finite-field
argument in the preceding artifact. Taking its extension limit
first and then e to infinity proves (2) for resonant angles.

Numbers 2 pi l/7^r are dense modulo 2 pi. The objective for each
fixed pair of finite profiles is continuous in theta. This proves
(2) for every real theta without asserting simultaneous independent
Gauss phases in different character layers.

## 4. Inversion and cyclic carries are only label permutations

The physical normalized Fourier operator gives K_(e,theta), whereas
the calculation used D_(e,theta)=K_(e,theta) P_N. To implement its
bilinear objective, use P_N J_h f as the input and J_h g as the
output. Inversion is a permutation of all N coset labels, so both
profiles remain Boolean.

In the coordinates (3) its exact carry rule is

    -(b+tM) = 0+(-t mod7)M                  when b=0;
    -(b+tM) = (M-b)+(-t-1 mod7)M            when b>0.

Thus there is no assumption that inversion acts independently on
the two displayed digits. Independent row/column label permutations
are permissible in the bilinear norm, and the previously proved R4
lift converts any such bilinear witness to the same-spin R norm.

## 5. Scope and the first weighted-seed check

The earlier seven-coset contraction ceiling 10 sqrt(6)/7 for the
weighted seed [-1,2;2,4] remains correct for that uncarried family.
It does not apply to (2): the carrier restores constants with a
vanishing defect. This is a strict expansion of the proved test
architecture, not a contradiction.

The bounded exploratory script
`tmp/resumed_convergence_full7_phase.py` exhausts all 2^13 Boolean
input profiles modulo global sign on a 501-angle grid, with ordinary
sign output optimization. It did not find a value above 7/2 for
[-1,2;2,4]. The angle grid is not an upper proof; the theorem (2)
does not assert that this single full seven-point family settles the
weighted seed.

Most importantly, J_h is a signed linear spin embedding:
`(J_h f)(J_h g)=fg` while `J_h(fg)=hfg`. Its image is not a unital
probability algebra. Tensoring independently realized modules is
valid, but multiplying arbitrary overlapping or noncommuting
actions has not been justified.
