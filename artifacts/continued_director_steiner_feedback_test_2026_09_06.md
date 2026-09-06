# An actual asymmetric-spectrum signing family for feedback tests

Date: 2026-09-06. The construction and spectral identity below are exact.
Its nonlinear-feedback consequence is now independently reconstructed using
the separate zero-first masked covariance theorem.
No minimizing or low-cap property is asserted.

## 1. Primary construction and exact mapping

Theorem 1 of Fickus, Mixon and Tremain, [Steiner equiangular tight
frames](https://arxiv.org/pdf/1009.5730), constructs a real frame from a
(2,k,v) Steiner system and a Hadamard matrix of order 1+(v-1)/(k-1).
The director read its full theorem and proof, pages 4–6 of the PDF, and
reconstructed the following binary-triple specialization directly.

Take the nonzero vectors of F_2^m as v=2^m-1 points. Its triples are
{a,b,a+b}; each distinct pair is in exactly one triple. Each point lies
in r=(v-1)/2=2^(m-1)-1 triples. There are b=vr/3 triples.
Take a Sylvester Hadamard H of order r+1 with first row all ones.

For each point, put the r nonconstant rows of H into its r incident
triple positions, in any fixed bijective order. Concatenate these v blocks
to obtain a b by N matrix E, N=v(r+1). Entries of E are zero or signs.
Two distinct columns belonging to the same point have inner product -1;
columns belonging to different points overlap in exactly one row and
have inner product plus or minus one. Every column has squared norm r.
Distinct rows of E are orthogonal, and each has squared norm 3(r+1).

Consequently

    S=E^T E-r I

is an EXACT hollow symmetric signing, and EE^T=3(r+1)I. Therefore

    S^2=(r+3)S+(N-1)I,          N-1=r(2r+3).                (1)

Its two eigenvalues are 2r+3, with multiplicity b, and -r, with
multiplicity N-b. The exact program
`computations/continued_director_steiner_signing_2026_09_06.py`
checks all incidence, Gram, sign, and polynomial identities; floating
eigenvalue computations are unnecessary.

## 2. What this actual family tests

Set B=S/sqrt(N-1) and gamma=(r+3)/sqrt(N-1). Then

    B^2=I+gamma B,
    ||B||op=(2r+3)/sqrt(N-1) ->sqrt(2),
    Tr(B^3)/N=gamma ->1/sqrt(2).                              (2)

These signings have a bounded normalized operator norm and a NONZERO
limiting cubic moment. They are not fixed trace-zero-Hadamard tensor
lifts, whose matrix-sign symmetry kills the very self-energy being tested.
They do not violate any result stated only for exact minimizers.

For Q=B^2, the cubic Schur transport has the exact simplification

    Q^(circ3)=I+[gamma^3/(N-1)]B,
    T=B Q^(circ3) B=Q+[gamma^3/(N-1)]B^3.                   (3)

In particular T=Q+O_op(N^-1), rather than I+o_op(1). The diagonal
of T is 1+gamma^4/(N-1). Its normalized off-diagonal correlations
are alpha_N B_ij, where

    alpha_N=[gamma+gamma^3(1+gamma^2)/(N-1)]
                /[1+gamma^4/(N-1)] ->1/sqrt(2).           (4)

## 3. Proved nonlinear-feedback consequence

The independently audited masked covariance theorem establishes, for a bounded
odd F with zero first Gaussian chaos and positive Gaussian variance,
the Gaussian covariance approximation for C=sign(BF(BS_0)), with fresh
independent Boolean seed S_0. On the family above, every fixed odd
Schur level beyond the first differs from I by O_op(N^-1). Thus its
normalized transported covariance tends to Q, and the arcsine formula
gives the actual expected normalized half-energy

    E[C^TBC]/(2N) -> 1/(pi sqrt(2)).                       (5)

For an old even mask H of mean p, the value is p^2/(pi sqrt(2)).
This is a positive self-energy, in contrast to the flat-involution
terminal theorem. Equation (5) does not follow merely from one-root
Gaussianity: its dependencies are the cross-root/nuclear covariance result
and hard-threshold approximation. The latter is justified here by the exact
positive variance floor T_ii>=||F||_2^2. See the complete independent proof in
`continued_audit_steiner_and_nonzero_first_falsifier_2026_09_06.md`, Section 2.
This proves a nonzero actual feedback self-energy on an explicit signing
family, not on arbitrary near-minimizers.

A concrete bounded zero-first F is a nonzero scalar multiple of
sin(z)-(exp(3/2)/2)sin(2z), scaled to have absolute value at most one.
The identity E[Z sin(tZ)]=t exp(-t^2/2) checks its zero first chaos.
The full response-cube feasibility |F|+H<=1 is a separate restriction
when claiming an energy for F plus the masked feedback; C alone is
a valid Boolean vector when H=1.
