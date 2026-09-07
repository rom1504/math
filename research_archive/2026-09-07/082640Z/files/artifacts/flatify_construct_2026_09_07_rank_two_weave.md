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
No numerical value above is an upper certificate.

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
