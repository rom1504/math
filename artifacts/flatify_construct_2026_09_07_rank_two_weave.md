# Rank-two Hadamard-tile construction and the joint-profile obstacle

Date: 2026-09-07, begun approximately 08:10 UTC. Construction proposed by
the director; exact algebra and numerical diagnostics reconstructed here.
This is an actual dense correlated sign construction. No new asymptotic
upper bound or seed-transfer inequality has yet been proved.

## Exact construction

Let k=2m admit a Hadamard matrix. For each fibre i choose a full k-by-k
Hadamard F_i and label its columns by (j,b), j in [m], b in {0,1}.
Let S be an arbitrary symmetric full sign matrix of order m. Put

    C_(i,a),(j,b) = S_ij/2 * f_i(a,j)^T H2 f_j(b,i),
    f_i(a,j) = (F_i(a,(j,0)), F_i(a,(j,1))).

For any two Boolean two-vectors u,v, u^T H2 v is +2 or -2.
Thus C has exactly sign entries and is symmetric. In the orthonormal
feature basis diag(F_i/sqrt(k)), C is k/sqrt(2) times the symmetric
orthogonal operator swapping i,j and applying H2/sqrt(2), with factor S_ij.
Therefore C^2=N I for N=mk=2m^2. Also trace C=0: each internal diagonal
is S_ii F_i(a,(i,0))F_i(a,(i,1)), whose sum vanishes by orthogonality.
Hollowing gives Q<=N^(3/2)/2. The objective remains the full Boolean cube,
not its spherical relaxation.

The reciprocal-column matching lower no longer gives one half. For one
matched pair, the two selected columns at each endpoint are orthogonal
sign vectors. Their attainable two-feature vectors obey ||h||_1<=k,
and individual columns attain (k,0). Thus the cross-block bilinear cap is
exactly k^2/2. Independent matching-component reversals give

    Q >= m k^2/4 = N^(3/2)/(2 sqrt(2))

for even m, with arbitrary internal blocks. This does not close the gap
between .353553 and .5.

If every F_i is the canonically paired Sylvester H_m tensor H2, the
construction factors as the reciprocal self-weave of S tensor H2.
Independent signed column permutations genuinely change feature routing;
row permutations and row switches are merely parent cube gauges.

IMPORTANT: under fully independent column signs, every outer S_ij can be
absorbed into the two columns at one endpoint of edge ij. Consequently
the random cross-block ensemble loses the input seed entirely. It could
produce a universal upper construction, but is not by itself a transfer
of Q(S).

For asymptotic analysis it is cleaner to DELETE all within-fibre rank-two
tiles and complete those fibres independently with O(k^(3/2))-cap actual
signings. In feature coordinates deletion simply sets the 2m loop
coordinates to zero; the remaining cross operator still has norm1.
Completion costs O(m k^(3/2))=O(N^(5/4)), a genuine power-saving error.
This avoids a loop-concentration issue: the original two loop columns can
carry leading Fourier energy despite using only O(m) feature slots.

## Numerical diagnostics

`computations/flatify_construct_2026_09_07_rank_two_weave.py` verifies exact
signs, symmetry, trace zero and the Hadamard identity for each instance.
It uses independent signed column permutations of a Sylvester frame,
with exact stored minimizing outer children when available. Coordinate
ascent with 1024 starts in each polarity gives ONLY lower bounds:

* m=4, N=32: canonical witness cap80 (.441941738); nineteen randomized
  instances produce cap84 or86 (.464038825 or .475087369).
* m=8, N=128: canonical witness cap680 (.469563097); randomized instances
  give witnesses676 through684 (.466800961 through .472325233).

Full actual matrices and verified spin witnesses are saved in
`computations/results/flatify_construct_2026_09_07_rank_two_weave_m{4,8}.json`.
Subsequent independent audit upgrades the witnessed86 order-32 records
to EXACT caps: every symmetric trace-zero Hadamard of order32 has Q<=86.
The proof is in `flatify_adversary_2026_09_07_rank_two_tile_audit.md`.
The 80/84 records and every order-128 record remain lower bounds only.

## Why a joint profile bound is nontrivial

For fixed fibre spins, write their two-column feature coefficients as a,b
at one end and c,d at the other. Independent feature phases give exact
cross-edge exponential moment

 .5[cosh(t a(c+d)/2)cosh(t b(c-d)/2)
    +cosh(t a(c-d)/2)cosh(t b(c+d)/2)].

Distinct cross edges use disjoint phases, so these moments factor conditional
on column magnitudes. Column permutations couple the allocation of magnitudes
to edges. A typed-row joint empirical-distribution argument is plausible;
neither positive-semidefinite graph Cauchy--Schwarz nor separate channel
caps are justified.

The following exact relaxed source must be handled. With probability1/2
take a uniformly from {+/-1}^2 and b=H2 a/sqrt(2), and with probability1/2
reverse the endpoint roles. Then a=H2 b/sqrt(2) exactly. The common scalar
marginal nu has masses

    nu(0)=1/4, nu(+/-1)=1/4 each, nu(+/-sqrt(2))=1/8 each.

H(nu)=2.25 log2, while joint four-coordinate entropy is3 log2. Consequently
D(pi||nu^4)=6 log2. This source alone is sufficiently expensive for an
unrestricted cube union bound at the usual D/4 normalization. But dilute
the whole four-vector by an all-zero atom of mass1-p and rescale the active
part by1/sqrt(p). Variance remains1 and alignment remains exact, whereas
D tends to0 as p tends to0. Therefore a variational bound allowing arbitrary
scalar sources while charging the full log2 cube entropy cannot prove a gap.
It must retain the entropy/count of ACTUAL Walsh-compatible row profiles.

This is not an actual Boolean counterexample. The displayed two nonzero
amplitude levels have irrational ratio sqrt(2), whereas all coefficients
of one finite Walsh transform share a rational lattice. Approximation by
actual Boolean profiles is unproved here. The source diagnoses only an
unrestricted distributional relaxation.

For clarity, the diluted source's exact divergence is

 D(p)=4 h(3p/4)-h(p)+3p log3-2p log2,

where h is binary entropy in natural units. Its scalar support fraction
is r=3p/4 and D(p)/4 is asymptotic to (2/3)h(r), not (3/4)h(r).

One elementary actual-cube counting bound can be proved without any
distributional assumption. For ANY invertible k-by-k linear transform,
the number L_s of Boolean inputs whose output has support at most s obeys,
for every integer q>=s,

    L_s <= binom(k,q) 2^q / binom(k-s,q-s).

Indeed a fixed q-coordinate output subspace has dimension q and contains
at most2^q Boolean input points: select q input coordinates whose projection
is injective on that subspace. Each support of size<=s lies in at least
binom(k-s,q-s) such q-subspaces. Double count. Taking s~rk and q~2rk,
for fixed r<=1/2, gives log L_s <= k h(r)+o(k). This improves the direct
bound k[h(r)+r log2], but still exceeds the relaxed alignment cost
(2/3)k h(r). It does NOT certify a bad actual profile: the subspace count
ignores Fourier idempotence and support structure. Sparse affine/tensor
Boolean constructions have substantially smaller row entropy.

## Finite profile transport tests and their precise scope

The exact Walsh magnitude census at k=4,8,16 has respectively2,3,8
histograms. Every count sums to2^k; each histogram records |H_k x|, not its
normalization. Scripts and results use the prefix `walsh_profile_census`.
The number of possible histograms in general is exp(O(k^(2/3))):
nonzero magnitude counts n_j satisfy sum j^2 n_j=k^2, and

 log prod_j (1-exp(-t j^2))^(-1)
 <= C t^(-1/2) sum_l l^(-3/2).

Optimizing t of order k^(-4/3) bounds the number of partitions into squares;
the zero count is determined. Thus histogram choices themselves cost
o(N) entropy across m rows when k=2m. This does not determine the much
more important counts WITHIN each histogram.

Here is the precise heterogeneous proxy used in the finite calculations.
Let nu_i be the magnitude law |H_k x|/sqrt(k) for histogram type i,
s_i=log(number of rows of type i)/k, and p_i type proportions. For
a,b,c,d>=0 define

 K_t(a,b,c,d)=E_signs exp[t(AC+AD+BC-BD)/sqrt(2)],

where A,B,C,D independently attach uniform signs to a,b,c,d. For scalar
potentials lambda_i on the support of nu_i, put

 G_p(t)=inf_lambda {
   sum_ij p_i p_j log E_(nu_i^2 nu_j^2)
      [K_t(a,b,c,d) exp(lambda_i(a)+lambda_i(b)
                       +lambda_j(c)+lambda_j(d))]
   -4 sum_i p_i E_nu_i lambda_i }.

The associated normalized cap expression is

    [sum_i p_i s_i + G_p(t)/4]/t.

This dual imposes only ROW-AVERAGED marginals: a type-i row may route
different coefficient distributions toward different neighboring types.
It does NOT impose nu_i separately on every i-j edge. Magnitude signs
are correctly integrated in K_t, not optimized for free.

For homogeneous k4/8/16 profiles, common-temperature worst expressions
are respectively .473227236, .472692704, .480315977 on the saved grid.
For k16 all28 two-type pairs with weights1/4,1/2,3/4 were tested at
t=4,8,16,32. Their common-temperature worst values are respectively
.504912263, .481819175, .480455021, .485034743. The maximum numerical
marginal residual in the homogeneous tests is below3.6e-8.

These are finite-profile transport proxies, NOT certified cap bounds:
the relevant row order must grow with m, higher-order mixtures have not
been optimized, loop/truncation bookkeeping must be included in a finite
partition argument, and no asymptotic bound on actual histogram counts
has been proved. The scripts are `rank_two_profile_transport.py` and
`rank_two_mixed_transport.py`, with corresponding JSON results.

The director has now supplied a rigorous finite conditioning certificate
in `flatify_director_rank_two_typed_certificate_2026_09_07.md`. This track
independently reconstructed its exact conditioning, unused-coordinate
factors, polarity symmetry, type overhead and actual completion: PASS.
Its inverse temperature is half the one used in the proxy scripts because
it normalizes x^T C x rather than H_C. No growing-order inequality has yet
been verified. Mixed-profile numerical marginal residuals are below9.2e-7.

## Primary literature scope

Carlet, Danielsen, Parker and Sole, “Self-dual bent functions,”
IJICOT1(4),384--399(2010), DOI10.1504/IJICOT.2010.032864, studies odd-variable
Sylvester Rayleigh quotients. Its Theorem3.3 supplies a LOWER bent-concatenation
bound; Theorem3.4 recursively improves witnesses. Neither supplies the
strict-subhalf uniform upper needed here. The reported finite normalized
Rayleigh value .883883 corresponds to our cap normalization .4419417.
Primary author-uploaded text was inspected at
https://www.researchgate.net/publication/210230166_Self-dual_bent_functions ;
the former author PDF https://www.ii.uib.no/~matthew/sdbentFinal.pdf now
redirects to a university department page, although indexed text remains.

A relevant exact-profile entropy input is Potapov,
https://arxiv.org/html/2303.16547 (v3,18November2024): the binary logarithm
of the bent-function count is at most (11/32+o(1))2^n, and the paper gives
separate s-plateaued count estimates. These statements concern exact
flat/plateaued spectra only. They do not justify entropy bounds for
approximately aligned or arbitrary profiles used in the present problem.
