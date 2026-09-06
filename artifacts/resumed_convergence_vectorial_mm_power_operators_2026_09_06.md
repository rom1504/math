# Nonlinear vectorial-MM profiles realize conjugated power permutations

Date: 2026-09-06. Status: an exact finite-field Fourier calculation and
new limiting balanced-profile lower test for R. This bypasses the
conditional-BENT-input reduction because a nonlinear profile of a
vectorial bent map need not itself be bent. No full constant-channel
restoration, arbitrary permutation synthesis, or R=T is asserted.

## 1. The available balanced-profile operator

Fix q=2^m, and let p be any integer coprime to q-1. On F_q, let
P_p be the permutation of additive Fourier indices

    alpha -> alpha^p,       with 0->0.

Let F_q denote the normalized real additive Fourier/Walsh matrix
using the field trace pairing. The real orthogonal operator

    T_p=F_q P_p F_q                                          (1)

fixes constants. For every symmetric finite seed B and every
Boolean profile array F on F_q whose columns have uniform mean
zero, the construction below proves

    R(B)>=(1/(2q)) ||T_p F B||_1.                           (2)

The output signs in this lower test are unrestricted. Balanced
Boolean profiles are possible because q is even. For example the
weighted-CHSH target correlation 3/4 is compatible with balanced
profiles on every q divisible by 16.

## 2. Arbitrarily large compatible extensions

Choose a positive integer e congruent to -p modulo q-1, so
gcd(e,q-1)=1. There are infinitely many s with

    gcd(e,q^s-1)=1.

To see this, for each odd prime ell dividing e and not dividing q,
the order of q modulo ell exceeds one, since ell does not divide
q-1. Taking s congruent to one modulo the least common multiple of
these orders makes q^s-1 nonzero modulo every such ell. Primes
dividing q cannot divide q^s-1 in any case.

Put N=q^s and choose d inverse to e modulo N-1. Thus y->y^d
is a permutation of F_N. Define the vectorial map

    Phi(x,y)=Tr_(F_N/F_q)(x y^d),       x,y in F_N.         (3)

For a scalar profile f:F_q->{+1,-1}, use the Boolean input

    f_N(x,y)=f(Phi(x,y)).

Every nonzero additive component alpha of Phi is a Maiorana–McFarland
bent function. But f_N, for general nonlinear f, is not required to
be bent; that distinction is the entire point of this construction.

## 3. Exact transform for every finite N

Expand f using probability-normalized Fourier coefficients on F_q:

    f(t)=sum_alpha fhat(alpha) psi_q(alpha t).

Use the normalized Walsh transform on the pair (x,y), whose order
is N^2 and normalization is 1/N. For alpha!=0, summation over x
forces the unique equation alpha y^d=a, and gives

    (1/N) sum_(x,y) psi_N(alpha x y^d+a x+b y)
       =psi_N(b(a/alpha)^e)
       =psi_q(t alpha^(-e)),
    t=Tr_(F_N/F_q)(b a^e).                                 (4)

The alpha=0 term is N mean(f) 1_((a,b)=(0,0)). Consequently

    U_(N^2) f_N(a,b)
       =(T_p f)(t)-mean(f)
                         +N mean(f) 1_((a,b)=(0,0)).       (5)

Here -e is congruent to p modulo q-1, so the nonzero frequency
permutation in (4) is exactly P_p. Formula (5) includes a=0: then
t=0 and every nonzero component in (4) equals one.

Both Phi(x,y) and t(a,b) have the same elementary distribution:

    Pr(t=0)=1/N+(1-1/N)/q,
    Pr(t=u)=(1-1/N)/q for u!=0.                            (6)

Indeed, for nonzero y (respectively a), multiplication by y^d
(respectively a^e) is invertible and the field trace is uniform.

If mean(f)=0, equation (5) is the full T_p action with no removed
component. Applying (5) coordinatewise to the fixed finite Boolean
array F, optimizing output signs, and using (6), gives

    lim_(N->infinity) (1/(2N^2)) ||U_(N^2) F_N B||_1
       =(1/(2q)) ||T_p F B||_1.

Every outer dimension N^2=2^(2ms) is even-dimensional Walsh, hence
equivalent by independent signed row/column permutations to an
allowed R4 tensor power in the bilinear norm. This proves (2).

For unbalanced f the delta in (5) vanishes in normalized L1, but
the uniform constant component is still removed. The resulting
general test uses T_p-E_q, not the full T_p. No unbalanced full
extension is being claimed.

## 4. What this does and does not change

These are Fourier-CONJUGATED power permutations, in contrast to the
previously programmable Fourier-DIAGONAL projective phase operators.
The nonlinear index permutation can move feature Fourier directions
between different frequencies, so it is a genuinely different
structured realization mechanism. This is not obtained by simply
assuming two available overlapping gates can be composed.

It does not give F_q P F_q for every permutation P. In (4), the
homogeneity of the monomial inverse lets the a dependence be absorbed
into t. An arbitrary permutation in place of y^d generally loses
that property and produces a varying family of operators depending
on a; it cannot silently be replaced by one prescribed P.

The first bounded search in
`tmp/resumed_convergence_power_permutation_search.py` checked the
distinct Frobenius power classes on q=64 with balanced input best
responses for [-1,2;2,4]. It found no improvement over 7/2. This is
not a norm upper bound or a claim of a persistent defect.

The positive result is the exact realization (2), with the explicit
nonflat vectorial construction (3) and finite transform identity (5).

## 5. General dual-bent-index realization and a simpler amplification

The same argument applies to ANY finite binary vectorial dual-bent
pair, and direct sums remove the need for field extensions. Precisely,
suppose maps Phi,Psi:F_2^n -> F_2^m, n even, and a permutation sigma
of F_2^m fixing zero satisfy, for every nonzero alpha,

    U_n[(-1)^(alpha dot Phi)](u)
                  =(-1)^(sigma(alpha) dot Psi(u)).         (7)

For r independent blocks define Phi_r=sum_j Phi(x_j) and
Psi_r=sum_j Psi(u_j), with sums in F_2^m. Tensoring (7) preserves
the SAME sigma. Every nonzero character expectation of Phi_r and
Psi_r has magnitude 2^(-nr/2), since their components are bent.
Their value distributions therefore approach uniform on F_2^m.

For every scalar profile f, the exact character expansion gives

    U_(nr)[f(Phi_r)]
       =(F_m P_sigma F_m f)(Psi_r)-mean(f)
                                  +2^(nr/2)mean(f)delta_0.

Thus every fixed balanced Boolean profile gives the full conjugated
operator F_m P_sigma F_m as an R lower test. No limiting family of
different finite-field index permutations is needed: one finite
dual-bent pair with the desired sigma is enough.

For the power case, take the finite map Phi(x,y)=xy^d over F_q
itself, Psi(a,b)=ba^e, and sigma(alpha)=alpha^(-e). Direct sums
already provide arbitrarily large even Walsh outer dimensions.

This is a conditional realization lemma, not a claim that every
permutation sigma admits a dual-bent pair. Known families and their
exact index restrictions must be checked before using (7).

The independent integer verifier
`computations/resumed_convergence_vectorial_mm_verify_2026_09_06.py`
has been replayed successfully. It checks all 256 profiles on F_8
inside F_64, and 34 selected profiles on F_16 inside F_256, including
unbalanced profiles and the exact constant-delta term. It also
checks primitive field orbits, both trace conventions, label
distributions, and full finite Walsh transforms.
