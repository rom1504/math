# Dense correlated construction trials: self-weaves and closed-quad trades

## Later diagnostic: free active diagonals (08:02–08:07 UTC)

The hollow objective permits resetting diagonal signs before every trade,
at exactly zero cap cost. For a four-set Q whose row product on all OUTSIDE
columns is constant, first switch one active vertex if that constant is -1.
Set each active diagonal to the product of its other three principal entries.
All four rows now have even product at every column. Test principal total
divisibility by 8 and conjugate as above; hollow the output. This is a valid
actual-sign operation and need not preserve the previous full completion's
spectrum. Among the 64 principal off-diagonal sign patterns, 32 pass: four
with all-positive forced diagonal, four all-negative, and 24 with two of each.

The exact script `computations/flatify_construct_2026_09_07_free_diagonal_trades.py`
tests this operation. The supplied adversary order-16 parent has cap 30 and
50 extremal rays. A first trade reduces this to 34 extremal rays while keeping
cap 30; all 24 first-step candidates have cap 30 or 32. A 501-state cap-aware
walk finds no cap 28. The Sylvester cap-32 input similarly retains cap 32
through 501 states. Results with both constant-parity signs allowed have
`free_diagonal_trades_signed_` filenames; initial positive-only 121-state
runs are preserved separately. These are heuristic orbit explorations,
not exhaustive orbit or lower-bound claims. This flexibility has not yet
produced a scalable cap inequality.

Date: 2026-09-07. Work begun after the 07:22 checkpoint. The exact algebra
below has been reconstructed jointly with the director. Numerical screens
are finite diagnostics, not asymptotic cap upper bounds.

## 1. Self-weave: an actual construction with a universal half-floor

Given any symmetric full sign seed S of order m, define the symmetric full
sign matrix of order m² by

    C_((i,a),(j,b))=S_ij S_aj S_bi.

Remove its diagonal when evaluating Q. For a spin matrix X with row i
equal to x_i, put h_i(j)=sum_a S_aj X_ia. The full quadratic is
.5 sum_ij S_ij h_i(j)h_j(i), before the diagonal correction.

The director found the decisive exact lower bound. Pair fibres whose seed
diagonal signs are opposite. In a pair i,j, choose x_i=alpha_i S_:j and
x_j=alpha_j S_:i, with alpha_i alpha_j selected to make the cross energy
sigma m². Their two internal energies cancel because they equal
.5(S_ii+S_jj)(<S_:i,S_:j>²-m)=0. Leave majority-diagonal fibres unpaired
and use their own seed columns, giving internal energy
sigma(m²-m)/2. Independent sign reversal of each whole component cancels
all cross-component energies in expectation. With p positive and q negative
seed diagonal entries, choose sigma to be the majority sign. Then

    Q(C-hollow diagonal)>=m³/2-|p-q|m/2.

The same argument with one COMMON retained row selector of size k gives
mk²/2-|p-q|k/2, a normalized half-floor .5sqrt(k/m)+O(1/sqrt(mk)).
No orthogonality or cap assumption on S is needed. Full self-weave therefore
cannot transfer sub-half seed caps asymptotically.

The director subsequently strengthened this to the exact universal matching
floor for ANY reciprocal-column weave with full-sign row bases F_i and
arbitrary internal fibre signings D_i. For a matching of fibres, choose
x_i=column_j(F_i), x_j=sigma S_ij column_i(F_j). Every internal energy is
unchanged by sigma; call their sum T. Independent whole-component reversals
cancel all nonmatching cross terms in expectation, leaving
T+sigma sum_matched k_i k_j. Choose sigma=sign(T). Thus

    Q>=max_matching sum_matched k_i k_j.

For even m and equal fibre sizes k this is exactly mk²/2, with NO diagonal
error and NO common-selector restriction. The original special proof above
is retained only as the exploration record; the stronger argument supersedes
it. It was independently reconstructed by this track.

## 2. A literal symmetric full-sign-preserving orthogonal trade

Let H be ANY symmetric full sign matrix. Suppose four rows indexed by Q
have entrywise product+1. On those coordinates let

    R_Q=J4/2-I4,

and let R be identity elsewhere. This is a symmetric orthogonal matrix.
Every outside column restricted to Q has even sign parity. R maps any such
vector to another sign vector: balanced vectors are negated, while constant
vectors are fixed. Hence RH has full sign entries on those rows.

To preserve signs after conjugation, check the four-by-four principal block
B=H_QQ. It has even parity in each row and column. The exact criterion is

    R_Q B R_Q has only signs  <=>  sum_ij B_ij is divisible by8.

For completeness, row sums of B are0 or+/-4. Write r_i=(row sum_i)/4 and
s=sum_i r_i. Then (R_QBR_Q)_ij=B_ij-2r_i-2r_j+s. If s is odd all entries
are even and feasibility fails. If s is even, constant rows cannot have
opposite signs; the cases zero, two, or four constant rows give signs
directly. This proves the criterion without invoking a classification.

For an admissible Q, H'=R H R is a symmetric full signing, isospectral to H.
If H is Hadamard then H' is Hadamard. This is a genuine old-edge operation,
not a signed permutation. Repeated overlapping trades can change a dense
edge set; no cap monotonicity is assumed.

Related primary switching literature is Orrick,
https://arxiv.org/html/math/0507515. Our exact symmetric conjugation and
its cap tests are proved directly above rather than imported from that paper.

## 3. It can escape obvious exact half-floor structures

Exact enumeration at order16 gives the Sylvester matrix20 extremal rays,
spanning all16 dimensions (8 in each eigenspace). The first admissible trade
in lexicographic order leaves12 extremal rays spanning only12 dimensions
(6 in each eigenspace). Thus no complete Boolean eigenbasis remains.

An exhaustive455-candidate check also shows this first switched matrix has
no simultaneous partition into four groups of four with every4-by-4 block
rank one. Such a partition WOULD be detected by taking its group containing
vertex0: projective column signatures of its four rows must form exactly
four groups of four (orthogonality excludes fewer); symmetry then forces
the same row partition. All candidates fail for the switched matrix, while
the initial Sylvester positive control passes.

Nonetheless, a401-state exact random-walk screen at order16 always retained
cap32. An81-state adaptive order64 experiment repeatedly selected trades
that destroyed known Boolean eigenvectors, but an exact verified +/-8
Boolean eigenvector was found at EVERY state by integer constraint solving.
Thus these experiments yield no cap gain or asymptotic upper construction.

The distinction is important: Shi et al.,
https://arxiv.org/html/2203.16439v2, Conjecture1 concerns existence for SOME
Hadamard matrix at each even square order; it does not assert every symmetric
Hadamard has a Boolean eigenvector. We did not use such a universal premise.

## 4. Directly applying trades to doubled actual optimal seeds

Let S=A_m+diag(d), and start C0=H2 tensor S. It is a symmetric full signing
with trace zero; its spectrum agrees with sqrt(2)diag(S,-S). For ANY seed,
the four rows top i,top j,bottom i,bottom j form a closed quadruple. Their
principal-block criterion holds when d_i=d_j=d.

The director's exact child calculation is:

* Outside the active pair, the new top child undergoes transposition i,j;
  the bottom child undergoes switching at i,j.
* The top ij edge becomes-d and the bottom ij edge becomes+d.
* Therefore if A_ij=-d, both children remain signed-permutation copies of
  the actual optimal children. Restoring their gauges gives bridge T S T,
  where T is the negative transposition at i,j.
* If A_ij=d, one edge in EACH child is additionally flipped. This costs at
  most2 per child, but exact child optimality must not be silently claimed.

There is a concrete repeatability obstruction. For simultaneous disjoint
paired trades, a cross tile is H2 tensor B2, where B2 is the seed's2-by-2
cross block between the two vertex pairs. Conjugating it by both R_Q maps
it to full signs iff product(entries B2)=+1, equivalently rank(B2)=1.
All16 patterns were independently enumerated. Thus a dense batch of such
trades is not automatically feasible for a generic optimal seed.

## 5. Exact finite seed screens

Every reported cap below is computed by full spin enumeration. In each
case diagonal completion was optimized exhaustively, except the explicitly
marked all-positive completion. Trade walks are NOT exhaustive orbit searches.

    seed m    seed Q    doubled initial cap    best screened cap
       5         4              13                   13
       6         5              18                   18   (all-positive diag)
       7         9              25                   25
       8        10              32                   32
       9        12              39                   39
      10        13              44                   44

For m8 with the deliberately unoptimized all-positive diagonal, trades
actually decreased cap40 to36, demonstrating nontrivial cap change, but not
an improvement over the best diagonal completion. At m7, cap25 is below
the aligned-child finite target2sqrt(13/6)*9, whereas m8's cap32 is above
its target2sqrt(15/7)*10. No uniform flatification inequality follows.

The independent adversarial track's different experiment preserves actual
children and varies a Hadamard bridge; its order16 cap30 construction is
not in this Hadamard-isospectral whole-matrix orbit and is not contradicted
by our cap32 screen.

Scripts and complete result matrices are under the prefix
`computations/flatify_construct_2026_09_07_` and corresponding results:
`symmetric_hadamard_switch`, `hadamard_eigen_escape_n64`, and
`doubled_seed_trades_m*`. The scripts preserve failed runs/initial choices
as separate result files. No fixed-step or random-walk cap theorem has been
proved, and no target-child optimality hypothesis has been replaced by the
finite experimental evidence.
