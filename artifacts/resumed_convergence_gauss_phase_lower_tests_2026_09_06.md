# Unbalanced coset lower tests and an exact 49-coset phase audit

Date: 2026-09-06. This is a rigorous additional family of lower tests for
the prescribed regularized seed norm R. It does not prove R=T or improve
the original universal lower bound. Numerical searches below are lower
witness searches, not upper certificates.

## 1. Primary input and the phase coherence proviso

The primary source is Kai-Uwe Schmidt, *Asymptotically optimal Boolean
functions*, [author manuscript](https://math.uni-paderborn.de/fileadmin-eim/mathematik/AG-Diskrete_Mathematik/Publications-schmidt/pw.pdf),
Lemmas 4–7 and their proofs. We use its normalized binary-field Fourier
transform and index-two Gauss-sum evaluation. With v=7^e,
m=3*7^(e-1), q=2^(sm), and a character of order 7^d, its normalized
Gauss sum is

    -(-1)^s zeta^(sigma * 7^(e-d)*s),
    zeta=(-1+i sqrt(7))/sqrt(8),  sigma in {+1,-1}.

The proof of Lemma 7 establishes density of both even and odd powers of
zeta. Even s therefore allows a continuous phase limit while retaining
even Walsh dimension. The sign sigma is character-dependent; the source
alone does NOT state its coherence between different orders. Within a
fixed order its two Frobenius orbits are quadratic residues and
nonresidues, with conjugate signs. Relative layer orientations must be
proved or retained as fixed signs, not freely optimized.

Fix a compatible primitive character in the base field F_(2^m), and
lift it by the norm to each extension. For each d let eta_d be its fixed
orientation on characters whose unit index is a quadratic residue
modulo 7. Along an even-extension subsequence, and after an irrelevant
global sign, the nontrivial phase multiplier has the form

    lambda_j(theta)=exp(i eta_d chi(u) 7^(e-d) theta),
    j=7^(e-d)u, 7 does not divide u, 1<=d<=e.                (1)

The same theta occurs in every layer. Define lambda_0=0 and the real
symmetric v-by-v matrix

    K_theta(a,b)=(1/v) sum_(j=1)^(v-1)
                         lambda_j(theta) exp(2 pi i j(a+b)/v).
                                                                  (2)

Changing Fourier convention changes theta's sign or reverses coset
indices, with no effect on the optimized lower tests. Orthogonality of
characters gives

    K_theta 1=0,       K_theta^2=I-E_v.                    (3)

Thus K_theta is an orthogonal involution on mean-zero functions, but
only a contraction on all functions. No action on constants is being
adjoined.

## 2. Arbitrary Boolean labels are legitimate: a delta is harmless in L1

Let B be any fixed real symmetric k-by-k seed and F any v-by-k array
with entries in [-1,1]. On every nonzero field element y, assign the
row F_c indexed by y's multiplicative coset c. Assign zero at additive
zero. Write f_q for this array on F_q and W_q for the additive Walsh
matrix, with U_q=W_q/sqrt(q).

The exact character expansion used in Schmidt's Lemma 4, without
assuming balanced labels, gives for a nonzero element in coset c

    (U_q f_q)(a)=K_(q) F(c) - q^(-1/2) mean(F),             (4)

where K_(q) has the nontrivial normalized Gauss multipliers. At zero,

    (U_q f_q)(0)=(q-1)/sqrt(q) mean(F).                     (5)

Along the selected subsequence K_(q) converges, up to the chosen global
sign, to (2). For fixed v and B, (4), (5) imply

    lim_q (1/q) ||U_q f_q B||_1=(1/v)||K_theta F B||_1.     (6)

In particular a nonzero input mean creates a large value at the single
additive-zero point, but its contribution to this normalized L1
objective is only O(q^(-1/2)). There is no need to discard an entire
multiplicative coset, fill it by discrepancy, or impose mean-zero
output labels. Nor does this observation restore the missing constant
component in L2.

The bilinear Boolean norm is at least its value on any two cube-valued
vectors, by separate linear maximization over the cube. Choose the
output signs to be the signs of U_q f_q B. Order-four Walsh and the
prescribed R4=J4-2I4 are equivalent by independent signed row/column
permutations; tensoring gives this equivalence in every even Walsh
dimension. Independent equivalence suffices for the bilinear norm.
The previously proved equality between R and the regularized bilinear
norm consequently gives

    R(B) >= (1/(2v)) max_(F in {+/-1}^(v*k))
                                    ||K_theta F B||_1.    (7)

The same bound holds for cube F. In fact the convex objective is
maximized on Boolean F, so allowing cube labels does not strengthen
the maximum in (7). The output in (7) is unrestricted Boolean; requiring
its mean to vanish unnecessarily weakens the lower test.

Products of finitely many such limiting tests are also legitimate by
tensoring the finite Walsh constructions and taking sufficiently good
approximations in each factor. This does not authorize multiplying
overlapping, noncommuting quotient actions on a common feature space.

## 3. Seven points: an exact weighted-seed obstruction within this family

For e=1 one may choose orientation so that

    K_theta=cos(theta)(J_inv-E_7)+sin(theta) C,
    C_ab=chi(a+b)/sqrt(7).

Here J_inv sends a to -a, C^2=I-E_7, and C anticommutes with J_inv
on mean-zero functions. This is a continuous family of reflections,
not merely the averaging reflection obtained from the H144 quotient.

For the remaining order-three full-sign class, merging its duplicated
coordinates gives the weighted two-coordinate seed

    B=[-1 2; 2 4],        T(B)=5/sqrt(2),       q(B)=7/2.

Put P=sqrt(2) diag(1,4). Then P dominates both B and -B, and
M=P^(-1/2) B P^(-1/2) is an orthogonal involution. For any odd v,
Boolean input F and Boolean output G, Cauchy-Schwarz gives

    |(1/v) Tr[G^T K_theta F B]|
       <= sqrt( (Tr P - mean(F)^T P mean(F)) Tr P )
       <= Tr P sqrt(1-v^(-2)).                            (8)

Every Boolean column on odd v atoms has absolute mean at least 1/v.
Therefore the entire v-coset test in (7) is bounded above by

    (5/sqrt(2)) sqrt(1-v^(-2)).                            (9)

At v=7 this is 10 sqrt(6)/7, strictly below 7/2: squaring and clearing
denominators reduces the strict comparison to 2400<2401. Hence no
choice of seven-coset labels or continuous phase can improve this
seed's unamplified value. This is only a ceiling on this single coset
family, not on its tensor products or on R(B).

## 4. Exact relative orientation for v=49

Use F_2[x]/(x^21+x^2+1), generator g=x, and character
chi(g)=exp(2 pi i/49). The standalone verifier
`computations/resumed_convergence_gauss49_verify_2026_09_06.py` traverses
the full 2^21-1 multiplicative orbit before returning to one. This
proves all nonzero elements are units, hence the quotient is a field
and g is primitive. It computes the absolute-trace linear functional
by squaring every basis element and then accumulates integer counts

    b_a=sum_(t congruent a mod49) (-1)^Tr(g^t).

All 49 bins satisfy the exact simple pattern:

    a not divisible by7:  b_a=-81 for chi_7(a)=+1,
                                 47 for chi_7(a)=-1;
    a=7t:                 b_a=-337 for t=0,1,2,4,
                                687 for t=3,5,6.          (10)

The ordinary quadratic Gauss identity on Z_7, or direct seven-term
algebra, now gives

    G(chi^j)=512(-1-i sqrt(7) chi_7(j)),       7 does not divide j;
    G(chi^(7u))=832-448 i sqrt(7) chi_7(u),    u=1,...,6.   (11)

Because (-1-i sqrt(7))^7=832-448i sqrt(7), division by sqrt(2^21)
shows the primitive QR layer is zeta-bar and the order-seven QR layer
is zeta-bar^7. Thus the two layer orientations agree for this explicit
compatible character. After changing theta's overall sign, the exact
49-coset family in (1) has eta_1=eta_2=+1.

This is an exact integer certificate for the relative sign omitted by
the primary lemma's formulation. No corresponding orientation claim
for e>=3 is made here.

## 5. Bounded computation, with its scope explicit

The exploratory script `tmp/resumed_convergence_gauss_phase_search.py`
forms (2) and uses alternating Boolean best responses over a theta
grid, then local scalar refinement. For v=49 the two-coordinate
unweighted H2-equivalent seed gave a lower witness about 1.392744972.
For the weighted seed [-1,2;2,4], the initial bounded scan did not
exceed 3.498542274, its theta=0 mean-loss baseline. Neither observation
is an upper bound or evidence of a positive global phase defect.

The positive theorem is the continuous family of valid contraction
tests (7), its exact 49-layer coherence certificate, and the scoped
seven-point obstruction (9). The polar-action realization needed for
R=T remains open.
