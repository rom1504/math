# Independent audit: binary phase closure and controlled-index promotion

Date: 2026-09-06. Status: all stated signed-module identities pass.
No arbitrary-permutation or overlapping-composition theorem is proved.

Main audited artifacts:

- `resumed_convergence_h144_three_cell_2026_09_06.md`;
- Section 5 of `resumed_convergence_cubic_affine_mm_modules_2026_09_06.md`;
- `resumed_convergence_arbitrary_controlled_index_promotion_2026_09_06.md`.

## 1. The reflection input really supplies joint L2 modules

The exact H144 three-cell certificate gives an invariant probability
algebra with action Q3=2E3-I. The integer verifier was rerun and passed.
For the ternary simplex matrix with one column per projective line in
F3^r, each nonzero pulled-back character has 3^(r-1) nonzero coordinates.
This is odd. Consequently the tensor product of Q3 factors induces the
single reflection Q_(3^r), not merely product parity on all functions.

Partitioning its atoms gives exact weighted reflections

    Q_p f=2 E_p f-f,

inside physical allowed outers, for ternary-rational vectors p. Let
p tend to uniform measure on F2^m, with m FIXED. For a character chi_a,

    Q_p chi_a-Q_uniform chi_a
      =2(E_p chi_a-delta_(a,0)) 1.

Thus its L2(p) norm is at most 2||p-uniform||_1. Gauge by chi_b on both
sides and negate globally. The resulting reflection has the desired
single-frequency phase -1 at b and +1 elsewhere, with the SAME error
bound on every character. Physical Boolean gauges may be kept as input
and output carriers, so no unauthorized outer generator is introduced.

Tensor one such gauged factor for each required negative frequency and
pull back along the sum of their labels. The input character chi_a is
the tensor product of chi_a in all factors. Tensor telescoping gives
an L2 error bounded by the sum of the component errors: all exact
operators are isometries, and all Boolean target factors have norm one.
The finite sum-label distribution tends to uniform. This proves a
single pair of Boolean carriers and labels that works simultaneously
for ALL characters, not merely a separately optimized scalar norm test.

The number of factors can be exponential in m, but is finite and fixed
before approximation accuracy tends to zero. For a fixed Boolean
profile f, Fourier expansion bounds its output error by
sqrt(2^m) times the maximum character error. This is harmless in the
stated order of limits. No uniform-in-growing-m approximation is claimed.

If the negative-frequency set is empty, an identity signed module can
also be supplied with uniform labels: the normalized Walsh operator on
F2^m x F2^m, followed by the output swap, fixes the bent carrier
(-1)^(x dot y) times every character of x. Equivalently use any of the
explicit zero-phase modules in the audited construction.

## 2. The selector identity and normalization

Assume the preceding physical operator U has component identity

    U[h chi_(a,b)(Phi)] ~= epsilon(a,b) h* chi_(a,b)(Psi),

where Phi=(Phi1,Phi2), Psi=(Psi1,Psi2), with block dimensions r,k.
Take the arbitrary phase epsilon(a,b)=(-1)^(b dot f(a)), with
f:F2^r -> F2^k completely unrestricted.

Add physical selectors s,t in F2^k and apply normalized Walsh on these
2k bits, whose matrix normalization is 2^(-k). The input component is

    h(x)(-1)^[s dot(t+Phi2(x))+a dot Phi1(x)+c dot s].

At output selector (u,b), the t-sum is 2^k times the constraint s=b.
It cancels the normalization EXACTLY. The output is therefore

    (-1)^[(u+c) dot b] U[h chi_(a,b)(Phi)],

which, after the old identity, becomes

    h*(y)(-1)^[b dot(u+Psi2(y))]
           chi_(a,c+f(a))(Psi1(y),b).

This proves the common-carrier index permutation

    sigma_f(a,c)=(a,c+f(a)).

It is an involution regardless of f. Its squared L2 component error
is exactly the average over b of the old component squared errors at
(a,b), because the remaining selector u is uniform. There is no error
amplification. Both new label laws tend to uniform: their first blocks
are old uniform-limit marginals and their second blocks are independent
uniform selectors. Even Walsh dimension 2k remains in the allowed
R4-equivalence class.

The exact integer verifier was independently rerun. It passed all
16 scalar controls on two input bits, all 256 two-bit vector controls
on two input bits, and 17 three-bit vector examples. The proof above,
not the finite replay, establishes the general claim.

## 3. The resulting norm test and its actual scope

For every fixed finite symmetric seed B and finite Boolean profile F,
the physical signed module gives the legitimate bilinear lower test

    R(B)>=(1/2) E ||[H P_(sigma_f) H F] B||_1.

Uniform limiting labels and finite Fourier expansion justify the
handoff, and the existing R4 bilinear lift gives an allowed same-spin
witness. Unbalanced profile columns are permitted: the common carrier
also transports the constant character.

For a scalar target bit, the family has 2^(2^r) distinct index maps.
It therefore genuinely leaves the previously counted quadratic-pencil
families. Their subexponential counting obstruction must not be applied
to these arbitrary truth-table translations.

However sigma_f only translates inside each fixed control fiber. To
attempt an arbitrary permutation pi using f(a)=pi(a), restriction to
the input indices (a,0) yields outputs (a,pi(a)), not pi(a). The retained
a-block is a genuine extra Fourier label. Setting its spatial output
label to zero is postselection on a small set; it is not a Boolean
isometric embedding and cannot be silently normalized away.

Thus a lossless erasure identity or an actual simultaneous overlapping
composition construction remains necessary. Tensoring independently
realized modules is valid; invoking abstract reversible-circuit
universality is not a substitute for that physical construction.
