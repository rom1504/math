# Independent audit of the seven-coset arithmetic-phase catalyst

Date: 2026-09-06. Status: the index-seven continuum and its prescribed-
outer bilinear lower tests are proved below. Higher-index relative phase
orientation and polar-coupling realization remain separate obligations.

## 1. Primary import and the even-extension quantifier

I read Lemmas 4--7 and their proofs in Kai-Uwe Schmidt, *Asymptotically
optimal Boolean functions*, [author manuscript](https://math.uni-paderborn.de/fileadmin-eim/mathematik/AG-Diskrete_Mathematik/Publications-schmidt/pw.pdf).
Lemma 6 permits every positive extension integer s; only the subsequent
Lemma 7 chooses odd s for its different application. Specializing e=d=1
in Lemma 6 gives, over F_(2^(3s)), normalized order-seven Gauss sums

    g_l=-(-1)^s zeta^(epsilon(l)s),
    zeta=(-1+i sqrt(7))/sqrt(8), epsilon(l) in {+1,-1}. (1)

Frobenius invariance of the canonical additive character gives g_(2l)=g_l,
and complex conjugation gives g_(-l)=conjugate(g_l). Thus epsilon is
constant on each of {1,2,4} and {3,5,6}, with opposite signs. Up to
inverting the chosen multiplicative character, epsilon(l) is the quadratic
character chi(l) modulo seven. This establishes simultaneous phase
coherence for all six nontrivial characters at index seven.

The square zeta^2=(-3-i sqrt(7))/4 lies in Q(sqrt(-7)) and is neither
of its two roots of unity, +1 and -1. Therefore zeta is not a root of
unity, and its argument divided by pi is irrational. Even powers of
zeta are dense on the circle. The extra minus sign in (1) for even s is
absorbed into theta by adding pi; it does not restrict the possible
limits. For every prescribed real theta there is an infinite subsequence
of EVEN positive s with

    g_l -> exp(i chi(l)theta),  l=1,...,6.           (2)

The field exponent 3s is then even, so the additive Fourier order is an
allowed power of four. No odd-dimensional Walsh factor is being smuggled
into the outer family.

## 2. Exact Fourier calculation and all normalization factors

Let N=2^(3s), fix a primitive multiplicative element, and index its seven
cosets by k in Z/7. For a vector f=(f_k), define its field lift by the
coset values and f(0)=0. Let fbar=(1/7)sum_k f_k. The normalized additive
Fourier transform is N^(-1/2)sum_y f(y)(-1)^Tr(ay).
At a nonzero output in coset j its value is EXACTLY

    (K_s f)_j-fbar/sqrt(N),
    (K_s)_(j,k)=(1/7)sum_(l=1)^6 g_l exp(-2pi i l(j+k)/7). (3)

At the additive zero the transform is (N-1)fbar/sqrt(N).
The trivial-character term in (3) is negative because its Gauss sum is
-1. Omitting the trivial character subtracts averaging; it does not
produce a constant-preserving probability-algebra action.

Let E be uniform averaging, P_(j,k)=1{j+k=0}, and
C_(j,k)=chi(j+k)/sqrt(7). The quadratic character Gauss identity gives

    K_theta=cos(theta)(P-E)+sin(theta)C.             (4)

In particular K_theta kills constants. On zero-mean functions it is the
director's U_theta=cos(theta)P+sin(theta)C. Direct integer character
sums give

    C1=0, C^2=I-E, PC=-CP.

Consequently K_theta is a real symmetric partial isometry and
K_theta^2=I-E. It is NOT an orthogonal constant-preserving extension on
all seven atoms. Every fixed finite collection of coset functions shares
the same subsequence (2), and convergence is uniform on its seven cosets.

The exact executable check
`computations/resumed_bound_audit_seven_phase_exact.py` verifies the integer
matrix identities and independently enumerates F8 and F64 trace-character
bins. The F64 bins are [5,-3,-3,1,-3,1,1], yielding the predicted even-
extension phase, up to the legitimate global character orientation.

## 3. Landing in precisely the allowed bilinear outer family

Use distinct notation R4(B) for the regularization with only the prescribed
R4=J4-2I4 tensor generator. For symmetric seed B the exact block identity

    z=(x,y,y,x),
    z^T(R4 tensor C)z=8 x^T C y

gives

    2R4(B)=sup_a beta(R4^(tensor a) tensor B)/(4^a)^(3/2). (5)

The converse inequality is q(C)<=beta(C)/2. The larger regularization
also allowing H144 satisfies R(B)>=R4(B); equality of that larger R
with the right side of (5) has not been proved here.

The additive Fourier sign matrix over F_(2^(3s)) is a row/column
permutation of the ordinary Sylvester Walsh matrix: a basis and its
trace-dual basis give the required independent coordinate identifications.
At even dimension, that matrix is independently row/column signed-
equivalent to R4^(tensor(3s/2)), by the elementary order-four Hadamard
equivalence and tensoring. Such independent equivalences preserve beta.
They need not preserve the same-spin quadratic maximum, which is why
the bilinear identity (5) is essential.

Lift fixed vector-valued coset functions f_k,g_i in [-1,1] and set ALL
their additive-zero values to zero. Then the outer normalized bilinear
objective is

    (N-1)/(7N) sum_j sum_(i,k) g_i(j) B_(i,k)
                  [(K_s f_k)(j)-fbar_k/sqrt(N)].

The finite bilinear sign maximum dominates this cube-valued test by
independent rounding. Taking the even subsequence gives the exact lower
test

    R(B)>=(1/2)|E_j sum_(i,k)g_i(j)B_(i,k)(K_theta f_k)(j)|. (6)

Neither f nor g actually needs zero coset mean in (6). The large
additive-zero spike is harmless here because the output lift was set to
zero at that one coordinate. Alternatively its normalized L1 contribution
is O(N^(-1/2)). What fails without zero mean is an L2 probability-algebra
identification, not the bilinear LOWER test. This agrees with the
convergence agent's stronger unbalanced-input observation.

With F the seven-by-seed-size matrix of input labels, optimizing each
output label in (6) yields

    R(B)>=(1/14)||K_theta F B||_(entrywise 1).        (7)

The zero-mean version in the director's note is therefore valid and is
a subfamily of (7). Constants are still annihilated by K_theta; (7)
does not repair that missing channel or identify a full orthogonal
seven-point model.

## 4. The remaining coherence and realization obligations

For index 7^e with e>=2, Lemma 6 gives phase magnitudes proportional to
s,7s,...,7^(e-1)s, but its displayed plus/minus sign depends on the
character. Frobenius fixes the sign on each orbit at a given character
order, and conjugation fixes its opposite orbit. Those facts alone do
NOT determine the relative orientation between different orders.
One must establish the cross-layer signs using compatible characters,
an additional Gauss-sum identity, or exact finite-field computation.
In particular a formula with fixed linked phases theta,7theta,... cannot
be imported from Lemma 6 while silently selecting each layer's sign.

Changing the primitive multiplicative character can reverse all layer
orientations together. It does not obviously permit reversing one layer
independently. The convergence agent is investigating this at index 49;
this audit does not pre-approve the result.

Tensor products of independently realized tests are allowed. Composing
two of these operators on an overlapping seven-point interface is not
automatically an allowed outer tensor construction. The continuum of
phases by itself does not prove arbitrary polar couplings, R=T for a
general seed, convergence of the original minima, or an improved universal
signing lower bound. Those are genuine further realization problems.

## 5. The index-49 integer certificate now passes independent reconstruction

The convergence audit subsequently supplied
`computations/resumed_convergence_gauss49_verify_2026_09_06.py`. I read
the entire verifier and replayed it. Its first version had a genuine
factor-eight error in its final power assertion:

    (-1-i sqrt(7))^7=832-448 i sqrt(7),

not 104-56 i sqrt(7). I reported this; the convergence agent corrected
both text and code. The corrected complete verifier PASSES.

The field proof does not rely on a primitive-polynomial table. In the
quotient by x^21+x^2+1, x is a unit because the polynomial has constant
term one. The code verifies that its powers first return to one after
exactly 2^21-1 steps. These unit powers are distinct; they exhaust every
nonzero residue, proving that the quotient is a field. The trace is
computed on every polynomial basis vector by 21 successive squarings;
linearity over F2 then justifies the bit-mask trace for every element.

Here is an independent algebraic reconstruction from the exact bins.
Write chi for the quadratic character modulo seven, with chi(0)=0.
Their whole 49-vector can be written as

    b(a)=-17-64 chi(a)
       +1{7|a}[192-512 chi(a/7)-512 1{a=0}].

For a unit frequency j, the constant term and the term periodic modulo
seven vanish after summing the other cyclic digit; the remaining
seven-term Gauss identity gives

    sum_a b(a)exp(2pi i ja/49)=-512-512 i sqrt(7)chi(j).

For j=7u, the periodic term instead contributes -448 i sqrt(7)chi(u)
and the exceptional terms contribute 832. Thus the primitive QR layer
is zeta-bar and the lower QR layer is zeta-bar^7 after normalization.
The corrected seventh-power identity proves their relative coherence.

I also reconstructed the weighted-seed ceiling in Section 3 of
`resumed_convergence_gauss_phase_lower_tests_2026_09_06.md`.
For B=[-1,2;2,4], P=sqrt(2)diag(1,4), the matrix
P^(-1/2) B P^(-1/2) is orthogonal. Since K^T K=I-E, weighted Frobenius
Cauchy--Schwarz gives

    |E G^T K F B|<=sqrt[(TrP-mean(F)^T P mean(F))TrP].

Every Boolean mean on odd v atoms has magnitude at least 1/v, yielding
the claimed half-energy bound (5/sqrt(2))sqrt(1-v^(-2)). At v=7 its
comparison with 7/2 is exactly 2400<2401 after squaring and clearing
denominators. This is only the contraction family's ceiling.

## 6. Resonant carrier restoration: independent proof of a stronger lower test

The convergence agent then proposed a way to restore the constant
channel. The following reconstruction checks both frequency support
and the nontrivial cyclic carry. The result is valid for bilinear tests;
it does not assert an invariant probability algebra of the original
symmetric outer matrix.

Fix a resonance theta=2pi l/7^r with r>=1. Choose e>r and let

    v=7^e, M=7^(e-1), m=7^(e-r).

Use the index-v Gauss phase limit with primitive-layer orientation
chosen positive. The common irrelevant global minus sign from even
extension is removed from the bilinear output. Its primitive frequency
multipliers are exp(i chi(j)theta). Every nonzero frequency divisible
by 7^r has multiplier ONE, regardless of the unknown orientation of
its lower character-order layer, because 7^r theta is a multiple of 2pi.
Middle-layer multipliers can remain unspecified.

Choose any +/-1 function h on Z/m with mean hbar, and extend it
periodically to Z/M and Z/v. Write a=b+tM uniquely with 0<=b<M and
0<=t<7. For arbitrary seven-point labels f, set

    x(b+tM)=h(b) f(t).                              (8)

The part h(b)[f(t)-mean(f)] has mean zero on every fine seven-point
fiber. Averaging over that fiber keeps exactly Fourier frequencies
divisible by seven, so this part uses ONLY primitive frequencies. The
remaining part h(b)mean(f) is periodic modulo m and therefore uses
only frequencies divisible by v/m=7^r. In particular no middle-layer
phase affects (8).

Factor the Hankel action as K_e=D_e P_e, where D_e is the circulant
Fourier-multiplier action and P_e is inversion. Inverting a cyclic index
has the exact carry

    -(b+tM)=(-b mod M)+(-t-epsilon_b)M,
    epsilon_b=1{b>0},                               (9)

with the fine coordinate understood modulo seven. On the primitive
subspace D_e acts independently on each fine fiber by the seven-point
circulant multiplier. On the low carrier part it acts as identity
minus global averaging. Therefore (8)--(9) give EXACTLY, in the
limiting coset operator,

    (K_e x)(b+tM)
       =h(-b mod M)[U_theta f](t+epsilon_b)
                       -hbar mean(f),              (10)
    U_theta=E_7+K_(7,theta).

Here U_theta is a full real symmetric orthogonal map preserving
constants. The seven-point inversion already appears in K_(7,theta);
the remaining carry is only the displayed cyclic output relabeling.

Choose the output labels

    y(b+tM)=h(-b mod M)g(t+epsilon_b).

Average (10). Since h^2=1 and every fine shift permutes seven atoms,

    E_a y(a)(K_e x)(a)
       =E_t g(t)(U_theta f)(t)
                       -hbar^2 mean(g)mean(f).      (11)

This identity also holds componentwise for vector labels and any fixed
seed bilinear form. Thus the error in this coordinated bilinear
compression is quadratic in the carrier mean.

Because m is odd, choose h with |hbar|=1/m. Let e-r tend to infinity
at fixed resonance. All fields and extension subsequences are chosen
after fixing the finite coset labels at each stage; a supremum over
allowed outer sizes permits these successive approximations. Resonant
angles with denominators 7^r are dense, so continuity in the finite
seven-point operator extends the resulting lower tests to every real
theta. In particular

    R(B)>=(1/14) max_(F Boolean) ||U_theta F B||_1,
    U_theta=E_7+cos(theta)(P-E_7)+sin(theta)C.        (12)

No coherence assertion for e>=3 was needed: resonance and the support
split entirely avoided the unknown middle layers. Formula (12) is
strictly broader than the earlier contraction tests and is not subject
to their odd-mean deficit bound.

The embeddings in (8) and (11) are SIGNED and differ on the two sides.
They are cube-preserving bilinear embeddings, not a common unital
probability-algebra embedding. Independent row/column sign changes and
the carry relabeling are allowed for beta; (5) supplies the subsequent
same-spin regularization. This distinction still forbids silently
composing overlapping U_theta actions as though all of them lived on
one fixed invariant algebra. A general polar-coupling theorem remains
unproved.

## 7. Independent audit of programmable projective phases

The complete theorem in
`resumed_convergence_programmable_projective_phases_2026_09_06.md`
passes a separate normalization and quotient reconstruction. For
G=F_7^r and one representative of each projective line, the matrix
C_(a,l)=chi_7(a dot l) satisfies CC^T=7^(r-1)I: independent row forms
give zero off diagonal, the diagonal full-vector sum is 6*7^(r-1),
and each projective line contributes six equal products. Consequently
every real phase obeying psi(ca)=chi_7(c)psi(a) is a linear combination
of the line characters. One tensor of the full seven-point rotations,
followed by the equal-fiber linear quotient, realizes its multiplier
exp(i psi(a)). Tensoring independently approximated signed physical
modules is legitimate for each fixed finite tensor size. This uses no
overlapping physical quotient composition.

Probability-normalized Fourier Parseval gives, for fixed Boolean
profiles F,G and real symmetric B, the exact optimized bilinear value

    |mean(G)^T B mean(F)|
       +2 sum_projective_l |sum_(c in {1,2,4})
                  Ghat(ca_l)^* B Fhat(ca_l)|.

Each projective line has just ONE adjustable phase for its three QR
frequencies; the three opposite frequencies are complex conjugates.
There is no freedom to phase the three QR modes independently.

For the weighted seed B=[-1,2;2,4], P=sqrt(2)diag(1,4), and
M=P^(-1/2)BP^(-1/2), the exact polar defect is

    Tr(P)-q^(-1)Tr(G^T U F B)
       =(2q)^(-1)||U F P^(1/2)-G P^(1/2)M||_F^2.

The placement of M is correct: expanding the inner product produces
P^(1/2)MP^(1/2)=B. The constant Fourier mode contributes half its
squared mismatch; every QR triple and its conjugate together contribute
the sum of the two triple squared norms minus twice the absolute
Hermitian inner product. This is the artifact's exact Equation (9).
Thus phase availability reduces this PARTICULAR family to Boolean
QR-triple phase-collinearity, not to an automatic equality R=T.

For odd r, chi_(7^r)(ca)=chi_7(c)chi_(7^r)(a), so the full Paley
rotation is included. The optional uniform-angle simplex formula also
has the claimed normalization: summing
chi_(7^r)(x)chi_7(Tr(ax)) over the six representatives on every line
and dividing by six gives chi_(7^r)(a)H(1), where
H(1)=sum_(Tr x=1)chi_(7^r)(x) has absolute value 7^((r-1)/2).
No factor six or square-root factor is missing.
