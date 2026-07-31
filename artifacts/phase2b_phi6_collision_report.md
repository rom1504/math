# Phase 2B: a scalable `phi_6` collision

## Result and classification

Let `phi_6(A)` be the vector counting every switching/permutation/global-
negation class among all principal restrictions of orders 4, 5, and 6 of a
signing `A`.

**Verified exact finite result.** There are two order-10 signings `A` and `B`
such that

```math
\phi_6(A)=\phi_6(B),
\qquad \operatorname{cap}(A)=19,
\qquad \operatorname{cap}(B)=21.                    \tag{P2B.1}
```

Their root-gauged edge codes are respectively

```text
5850642905
28771662001
```

and their common 16-component profile is

```text
(45,165; 3,64,127,58; 0,2,11,2,11,78,39,8,55,4).
```

The semicolons separate orders 4, 5, and 6.  Exact enumeration of all 512
projective spins gives energy ranges `[-17,19]` and `[-17,21]`.

The search itself was **deterministic but nonexhaustive**.  It traversed an
odd affine permutation of the `2^35` root-gauged/global-negation
representatives and found the collision after 555 distinct candidates.  The
reported pair was then independently reconstructed and verified by
`phase2b_verify_phi6_collision.py`; no solver inference is needed for (P2B.1).

Reproduction:

```bash
g++ -O3 -std=c++17 -Wall -Wextra -pedantic \
  computations/phase2_profile_collision_n8.cpp \
  -o /home/math/quadra/tmp/phase2_profile_collision_search

/home/math/quadra/tmp/phase2_profile_collision_search --sample 10 1000000

.venv/bin/python computations/phase2b_verify_phi6_collision.py
```

The exact matrices, energy histograms, class profiles, hashes, and positivity
certificates are in
`computations/results/phase2b_phi6_collision_n10.json`.
The raw deterministic search record is
`computations/results/phase2b_profile_collision_search_n10.json`.

## Stronger oriented profile equality

The pair satisfies more than (P2B.1).  If global negation is *not* quotiented
out, the complete switching/permutation class histograms agree at every order
`r<=6`.  The numbers of oriented classes at orders 2 through 6 are

```text
1, 2, 3, 7, 16.
```

The verifier compares all these histograms exactly.  This stronger fact is
what makes the blowup below rigorous; equality only after global negation
would not automatically preserve profiles when vertices repeat.

## Balanced twin blowup preserves every profile through order 6

For an integer `L>=1`, define `T_L(A)` by replacing every base vertex `i` by
`L` twins, putting sign `+1` inside each twin class, and putting the base sign
`a_ij` on every edge between the classes `i` and `j`.  Define `T_L(B)` in the
same way.

**Proved profile preservation.** For every `k<=6`, the two blowups have the
same complete oriented switching/permutation profile, hence the same
`phi_6` profile.

To prove this, take a `k`-vertex restriction of a blowup.  Its support is an
`r`-vertex base restriction, `r<=k`, together with a positive occupancy vector

```math
(m_1,\ldots,m_r),\qquad 1\le m_i\le L,qquad\sum_i m_i=k. \tag{P2B.2}
```

The induced signing is obtained by replacing vertex `i` of the support by
`m_i` positive twins.  A switching/permutation equivalence between two
oriented base restrictions lifts by switching whole twin classes and
permuting their fibers.  For each oriented base class, summing over all
supports and all occupancy vectors (P2B.2) is therefore a fixed linear
transformation depending only on `k` and `L`.  The base histograms agree for
every `r<=6`, so the transformed histograms agree for every `k<=6`.

As an independent finite audit, the verifier constructs both twofold blowups
of order 20 and checks all `binom(20,k)` restrictions for `k=4,5,6`.

## Exact blowup caps

Write the ten spin sums inside the twin classes as

```math
z_i\in\{-L,-L+2,\ldots,L\}.
```

If `Q_A=A+I`, the blowup energy is exactly

```math
H_{T_L(A)}(z)=\frac12 z^{\mathsf T}Q_Az-5L.          \tag{P2B.3}
```

The diagonal entries of `Q_A` are positive.  Holding the other coordinates
fixed makes `z^T Q_A z` convex in each coordinate, so its maximum on the box
`[-L,L]^10` occurs at a vertex.  Since the positive one-sided caps of `A` and
`B` are 19 and 21,

```math
\max_z\frac12z^{\mathsf T}Q_Az=24L^2,
\qquad
\max_z\frac12z^{\mathsf T}Q_Bz=26L^2.               \tag{P2B.4}
```

For each matrix, the verifier gives all ten positive leading principal minors
of `Q+4I`.  Sylvester's criterion therefore proves `Q+4I` positive definite.
Consequently

```math
z^{\mathsf T}Qz>-4\lVert z\rVert_2^2\ge-40L^2.      \tag{P2B.5}
```

For `L>=3`, the positive extrema from (P2B.4), after the `-5L` correction in
(P2B.3), dominate the absolute value allowed by (P2B.5).  Hence

```math
\operatorname{cap}(T_L(A))=24L^2-5L,
\qquad
\operatorname{cap}(T_L(B))=26L^2-5L,                \tag{P2B.6}
```

and the exact cap gap is `2L^2`.

This is a **scalable falsifier of universal fixed-profile control**: equal
`phi_6` does not even determine cap to `o(N^2)`, where `N=10L`.

## Why this is not yet a low-cap landing falsifier

Both caps in (P2B.6) are `Theta(N^2)`, far above the project scale
`Theta(N^(3/2))`.  The construction therefore refutes a universal theorem
that `phi_6` controls cap, but it does not show that `phi_6` fails inside a
carefully restricted low-cap family.

There is a rigorous obstruction to the most direct embedding.  Suppose a
larger signing is divided into `K` modules of common size `s`, with every
cross block constant in sign so that unrooted module profiles can be
substituted safely.  Let `D` be the order-`K` signing of cross signs.  On spins
that are constant inside every module, all internal module energies contribute
a fixed constant, while the cross energy is `s^2 H_D(sigma)`.  Comparing its
maximum and minimum over `sigma` gives

```math
\operatorname{cap}(\text{parent})
\ge\frac{s^2}{2}\bigl(P_D+Q_D\bigr)
\ge\frac{s^2}{2}M_K.                                \tag{P2B.7}
```

The project's uniform lower bound `M_K>=c K^(3/2)` and `N=Ks` turn (P2B.7)
into

```math
\operatorname{cap}(\text{parent})
\ge (c/2)N^{3/2}\sqrt{s}.                            \tag{P2B.8}
```

Thus an `O(N^(3/2))` modular construction requires bounded `s`.  With bounded
modules, even replacing every module by a member of this collision changes
the sum of internal caps by only `O(N)`, below the `N^(3/2)` scale.  Growing
uniform modules preserve the profile but force excessive cross cap.

The constant-block obstruction does not apply to a common non-rank-one micro
pattern.  That escape in fact gives a scalable separation at the correct cap
scale.

## A low-cap-scale Hadamard separation

Let `H_k` be the symmetric Sylvester Hadamard matrix of order
`k=4^r`, `r>=1`, and put `D_k=diag(H_k)`.  Define

```math
S_A(k)=A\mathbin\otimes H_k
       +I_{10}\mathbin\otimes(H_k-D_k),              \tag{P2B.9}
```

and define `S_B(k)` analogously.  The diagonal is zero and every off-diagonal
entry is `+/-1`, so these are signings of order `N=10k`.

**Profile preservation.** For every restriction order at most six,
`S_A(k)` and `S_B(k)` have identical oriented profiles.  Indeed, a restricted
support on macro vertices sees the same fixed micro-coordinate pattern
multiplied by the corresponding base signs.  An oriented switching of a base
support lifts by switching every selected micro vertex in the corresponding
fiber.  Summing over all macro supports and all micro-coordinate occupancies
is again a class-dependent linear transformation of the equal base profiles.

This is the non-rank-one analogue of the twin-blowup proof; it does not require
rooted or externally colored equality.

Put `Q=A+I`.  Since

```math
S_A(k)=Q\mathbin\otimes H_k-I_{10}\mathbin\otimes D_k,
```

write `U=H_k/sqrt(k)` and decompose each Boolean fiber spin as
`x_i=p_i+q_i`, where `p_i` and `q_i` are its projections onto the `+1` and
`-1` eigenspaces of `U`.  With

```math
G^+_{ij}=\frac{\langle p_i,p_j\rangle}{k},\qquad
G^-_{ij}=\frac{\langle q_i,q_j\rangle}{k},
```

the matrices `G^+` and `G^-` are positive semidefinite and
`diag(G^+)+diag(G^-)=1`.  The leading quadratic term is

```math
\frac{k^{3/2}}2\operatorname{tr}\bigl(Q(G^+-G^-)\bigr). \tag{P2B.10}
```

The exact rational vector

```math
y=\frac1{200}(445,490,661,668,436,645,405,427,485,513) \tag{P2B.11}
```

satisfies

```math
\operatorname{Diag}(y)-Q/2\succ0,
\qquad
\operatorname{Diag}(y)+Q/2\succ0,
\qquad
\sum_i y_i=\frac{207}{8}.                            \tag{P2B.12}
```

This is an exact certificate, not a floating-point SDP claim.  After
multiplication by 200, the verifier checks every leading principal minor of
the two integer matrices `Diag(445,...,513) +/- 100Q`; all are positive, so
Sylvester's criterion proves (P2B.12).  Pairing the two PSD inequalities with
`G^+` and `G^-` in both signs bounds the absolute value in (P2B.10).
The separate program `phase2b_hadamard_lift_theorem_audit.py` independently
checks the same certificate by an exact rational `LDL^T` factorization and
checks the tensor witness at `k=4`.

Every nontrivial Sylvester matrix has trace zero, so the diagonal correction
in (P2B.9) contributes zero to every Boolean energy.  Therefore

```math
\operatorname{cap}(S_A(k))\le\frac{207}{8}k^{3/2}.  \tag{P2B.13}
```

For the second base, take a Boolean spin attaining its positive energy 21;
then `s^T(B+I)s=52`.  The order-4 Sylvester matrix has the Boolean `+2`
eigenvector `(-1,-1,-1,1)`, whose `r`-fold tensor power is a Boolean
`+sqrt(k)` eigenvector of `H_k`.  The product spin in (P2B.9) gives

```math
\operatorname{cap}(S_B(k))\ge26k^{3/2}.             \tag{P2B.14}
```

Combining (P2B.13)--(P2B.14) proves the **scalable correct-scale gap**

```math
\operatorname{cap}(S_B(k))-\operatorname{cap}(S_A(k))
\ge\frac18 k^{3/2}
=\frac{1}{8\,10^{3/2}}N^{3/2}.                      \tag{P2B.15}
```

Both families have cap `O(N^(3/2))`: (P2B.13) proves it for the first, while
the Hadamard spectral norm and the fixed order-10 base give it immediately for
the second.  The exact order-20 audit gives caps 54 and 52 for the first
Hadamard lift; heuristic larger-order values are preserved separately in
`phase2b_hadamard_profile_lifts.json` and are not used in the proof.

Equation (P2B.15) is a scalable obstruction at the project's cap scale.  It
proves that fixed `phi_6` cannot control cap even to `o(N^(3/2))` over all
low-scale signings.  It still does not prove a landing gap near the unknown
optimum: the certified constants here are substantially above the best known
upper constant.  Thus `phi_6` may remain a descriptive coordinate inside a
more restrictive near-optimal family, but it cannot be the quantitative state
that controls cap or composition by itself.
