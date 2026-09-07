# Independent reconstruction of the balanced-width obstruction

2026-09-07. Source:
`flatify_construct_2026_09_07_balanced_overlap_obstruction.md`. PASS.
The width conclusion, rather than only an absolute-cap bound, is essential
to the source's arbitrary-child-polarity certificate obstruction.

## Elementary spectral estimate and uniform exceptional-pair control

For a unit eigenvector with eigenvalue magnitude lambda,
||v||infinity<=sqrt(n)/lambda, because every matrix entry has magnitude
at most one. Cube extension gives beta(A)>=lambda/||v||infinity^2, hence
lambda^3<=n beta(A)<=4nQ(A). Thus Q(A)<=C n^(3/2) implies
||A||op=O_C(n^(5/6)), without complex interpolation. Consequently

    Tr A^4=O_C(n^(11/3)), ||A1||^2=O_C(n^(8/3)).

For a uniform distinct pair u,v, both E[(A^2)_uv^2] and
E[(d_u-d_v)^2] are O_C(n^(5/3)). Outside O_C(n^(-1/6)) of the pairs,
|(A^2)_uv|<=n^(11/12) and (d_u-d_v)^2<=n^(11/6). These looser exponents
still give R2=2n+o(n), R1^2=o(n^2), uniformly in A and in the switched
center. The identities for R2 and R1 in the source are exact; the omitted
u,v terms in R1 cancel because A is symmetric.

## Pairing calculation and Gaussian-sized field

Conditional on a designated J pair, its remaining M=n-2 coordinates are
uniformly permuted. The first 2l are paired into I. Any unordered pair
occurs among these l pairs with probability 2l/[M(M-1)]. Summing
(r_i-r_j)^2 over all unordered pairs gives M R2-R1^2. This proves exactly

    E S=2l/(M-1)*(R2-R1^2/M).

Swapping two entries of the exposing permutation changes at most two
terms of S, each between zero and sixteen. Coupling conditional tails
by one transposition bounds every Doob increment by an absolute constant,
so Var(S)=O(n) uniformly. Thus S/n tends in probability to 2p, where
l=pn/2+O(1), after excluding the vanishing bad-pair fraction.

Conditional field coefficients have magnitude at most four. On S of
order n, their normalized maximum tends uniformly to zero. Expanding the
product of characteristic-function cosines proves the normal limit;
bounded normalized second moments give uniform integrability of the
absolute first moment. The mean field magnitude is therefore
(2sqrt(p/pi)+o(1))*sqrt(n). Linearity of expectation over the J pairs,
without requiring their fields to be independent, yields a pairing and
I signs with cross energy at least
[(1-p)sqrt(p)/sqrt(pi)-o(1)]n^(3/2).

## Width and precise certificate scope

Optimizing each J pair's sign gives cross energy L. Reversing all J spins
leaves the internal I and J energies unchanged and preserves balance,
so two actual full sign vectors have energies T+L and T-L. The width is
at least 2L. At p=1/3 this is
4/[3sqrt(3pi)]*n^(3/2)-o(n^(3/2)). Odd orders lose only O(n) on deleting
and restoring a vertex, giving the nearest balanced slice.

For two children, the width of the sum of their separately optimized
balanced energy ranges is the sum of their widths. Its maximum absolute
value is at least half that width, irrespective of global child polarity.
Adding the scalar dephased-bridge envelope at zero overlap therefore
costs at least [1+4/(3sqrt(3pi))-o(1)]n^(3/2), exceeding the relevant
2sqrt(2)c budget throughout c<1/2.

This is a failure of the SEPARATELY OPTIMIZED scalar-overlap upper
certificate. The bridge envelope may be very loose at the child-energy
extrema, so no actual parent lower bound follows. Nor does the result
restrict an operation which changes old internal edges. It nevertheless
rigorously identifies a nonvanishing energy-window exception on every
chosen center slice, including for actual optimal children.

## Fixed finite partitions

The source's subsequent Section 4 also passes reconstruction. Pairing
only inside macroscopic atoms changes the exact field expectation into
p times total R2 minus p times the atom-mean projection loss, with O_R(1)
denominator errors. The loss is ||P A(e_u-e_v)||^2 for the atom-constant
projection. Its mean is bounded by O(||PA||F^2/n) using uniform vertex
marginals; ||PA||F^2<=R||A||op^2. The squared column-correlation average
is controlled because within-atom pair density is at most O(1/epsilon)
times uniform pair density when atom sizes exceed epsilon*n. Thus the
same variance and CLT proof applies for fixed R and epsilon.

Omitting small atoms from the greedy construction does not license an
edge-deletion bound. The source instead extends BOTH witnesses by the
SAME random balanced spins on each omitted even atom. Omitted internal
energies cancel in their difference, and every added cross difference
has expectation zero. Some common extension preserves the constructed
width. Letting epsilon tend to zero after n recovers the full leading
constant uniformly over partitions with fixed R. At most R odd leftovers
cost only O_R(n). This legitimately defeats any bounded number of scalar
center overlaps, but asserts nothing for growing partition complexity.

## Arbitrary fixed-rank features

Section 4 of `flatify_construct_2026_09_07_arbitrary_subspace_bridge.md`
also passes independent audit. For an n-by-r orthonormal basis U, at most
rn/T^2 rows have ||sqrt(n)U_i||>T. Their contribution to ||U^T x|| is at
most sqrt(rn)/T, since restriction of U has operator norm at most one.
Quantizing the other scaled rows to coordinate mesh epsilon produces a
finite partition independent of n when r,T,epsilon are fixed. Its
unscaled matrix error has Frobenius norm at most sqrt(r)*epsilon, so
balanced spins in those classes have normalized projection at most
sqrt(r)*epsilon+sqrt(r)/T+o(1). Odd-class residuals are harmless because
the number and quantized values of the classes are fixed.

Applying the finite-partition width theorem gives width at least
(2L0-o(1))n^(3/2) on every fixed positive low-projection threshold, uniformly
over rank-r subspaces varying with n. Thus a uniform absolute-energy
threshold there must be at least L0. The source's three-case parent
certificate max{2a+1,c+a+sqrt(1-tau),2c+1-tau} is correctly derived; its
first term already fails for fixed-rank features throughout c<=1/2.
No assertion for polynomially growing rank follows from this quantization.
