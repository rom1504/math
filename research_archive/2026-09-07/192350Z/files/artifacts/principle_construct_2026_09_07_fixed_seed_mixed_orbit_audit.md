# Independent audit: fixed-seed mixed-orbit weave

2026-09-07. **PASS**, with two minor normalization/scope clarifications communicated to the author. I independently read the full candidate `principle_synthesis_2026_09_07_fixed_seed_mixed_orbit_weave.md` and the full original `continued_convergence_recursive_orbit_bound_2026_09_06.md`.

## Mixed-group estimate and the recursion

The input orbit covariance has norm exactly the signed-input Gaussian orbital average: its finite Gram matrix is nonnegative and has constant row sum. This statement does not require any symmetry of the OUTPUT group.

Ordinary-permutation invariant symmetric tensors of degree d have dimension at most p(d), rather than the signed group's p(d/2) in even degrees. Summing through degree D and applying the partition generating function gives exp(pi sqrt(2D/3)). The existing degree cutoff and Poisson tail therefore still cost only exp(O(sqrt m)). Orthogonal transformations preserve degree, so the covariance/rank argument proves the stated mixed-group estimate with no output signs inserted.

Ordinary output concatenation is projection onto a larger invariant subspace and has the stated direction. The deletion estimate improves to sqrt(m) because the raw Gaussian diagonal kernel is exactly one.

The recursion uses ordinary output norms at every node. Conditional on an input vector, its fresh signed input permutation is independent of both child bases. Hence the same signed pair tables and entropy rates as in the original recursion apply. At terminal nodes only, the mixed-group estimate controls the ordinary output norm by the signed INPUT orbital norm. Terminal dimensions and their energy bounds change constants only at fixed depth. There is no unsigned-input Bellman problem hiding inside the induction.

With H= sqrt(m) U^T, the last input signed permutation becomes a physical row gauge/permutation of H. It is not a final output-column sign matrix. This transpose convention is essential and checks algebraically.

## Fixed macro signs and the cap calculation

For raw Gaussian scalar features, replacing b by -b is a unitary reflection, acting by (-1)^d on degree d. Therefore every fixed macro edge sign is a norm-one twist on one endpoint of a tensor-network edge. Graph Cauchy--Schwarz bounds the absolute contraction by the product of vertex Hilbert norms even with these twists. After averaging ordinary column permutations, each squared vertex norm is the ordinary-permutation orbital average of its signed spectrum.

The deficit identity sums ordered macro pairs; each off-diagonal macro edge occurs twice. Exponentiating -tD/(2k) gives exactly the product exp[-t(v_i(j)-sigma S_ij v_j(i))^2] over unordered edges. The Markov factor is exp(t gamma m^2). Conditioning on each diagonal coordinate leaves independent uniform permutations of the other coordinates, and dropping the nonnegative diagonal deficit terms has the correct direction. The deletion bound then pays only a polynomial factor per row. Physical hollowness changes the cap by at most N/2, not by a leading macroblock deletion.

Independent fibre bases produce a product of EXPECTED row sums. The old fixed-depth row certificate consequently applies uniformly to every prescribed macro S. It does not depend on the numerical cap of S.

## Exact non-erasure statistic

A fixed physical row is a concatenation of independently signed uniform terminal rows. At the first split this follows because each child's fresh input gauge makes the law of a queried physical row independent of the queried index. The common root sign and the possible relative minus sign do not correlate the children: both child row laws are already symmetric. Induction proves the claim for all fixed-depth leaves.

The terminal physical Hadamard should be normalized to have an all-positive row; equivalently the author may dephase both its rows and columns. With this convention, averaging products over four distinct coordinates of a signed-uniform terminal row gives 1/(q-3). The second distinct-coordinate moment is zero after averaging the positive and balanced rows, and odd moments vanish by row-sign symmetry.

Under the final ordinary permutation, four queried macro coordinates either fall in a single leaf block or have zero expectation: odd occupancies vanish, and a two-plus-two split factors into two zero second moments. Thus

    mu4=b q(q-1)(q-2)/(m)_4

is exact. Five independent fibres give the displayed mu4^5 times macro K5 product. At fixed depth this moment is of order m^-5, so conditioning on a cap event of failure probability exp(-c m^2) preserves its sign for all large m.

The K5 reconstruction argument correctly recovers every four-cycle product by symmetric differences of four five-cliques on seven vertices. These products determine the OFF-DIAGONAL complete signing modulo vertex switching and global sign. Arbitrary diagonal signs are not determined by K5 products; this is a wording clarification, not a construction failure.

The result is therefore genuinely seed-sensitive in law and is not equivalent to absorbing S into independent final column signs. Its current cap bound is nevertheless seed-blind numerically. The polynomially small cycle statistic is not a leading pressure gain, and the independent-edge tilted stability penalty cannot be imported without a new argument.

## Exact order-eight sector response: independent replay PASS

I also read `principle_synthesis_2026_09_07_holonomy_gap_and_exact_test.md` and replayed its exact integer/rational checker. It enumerates 112 flat spectra among 1120 retained four-coordinate input words, and all 23,551 degree-zero-or-four subgraphs. The reported counts, matrices, cap values, moments, and polynomial coefficients reproduce.

The raw sign kernel is c_t(1+u_t S ab). Since the seven-port row moments vanish except at degrees zero and four, a surviving subgraph with s active vertices has exactly 2s edges. Its term is therefore A_s(S)(-u_t^2/35)^s, giving the claimed polynomial and equality between S and -S. The two-replica row norm is independently

    rho^2=c_t^7[1+binom(7,4)u_t^4/35^2]
         =c_t^7(1+u_t^4/35).

The opposite response signs can be checked without the Bernstein implementation. For z in [-1/35,0], the two cubic factors obey

    q_580(z)<=-8+72/1225+5/42875<0,
    q_722(z)>=24-68/35-24/1225>0.

Multiplication by z^5<0 proves the strict increase for atlas580 and strict decrease for atlas722 at every positive t. This is an exact favorable SELECTABLE child response in the specified actual finite row sector, not a statement uniform over all sectors or an asymptotic cap improvement.

For the general tensor criterion, averaging the common actual row law first gives the exact spin-orientation factor 2^(mk+1) times the averaged contraction. Substitution of Gamma and a_m into Markov's bound has the stated normalization. The seed target in the active cap interval gives a positive deficit parameter gamma, as required by the displayed event bound.
