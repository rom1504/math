# Independent reconstruction: joint rare packets give actual biased response

2026-09-07. **Director proof reconstruction and dependency audit.** The
positive construction is in `principle_invent_2026_09_07_stratified_marked_selector.md`.
This note records independent reasoning, not a previous audit verdict used
as a proof step. Convergence of the original minimax sequence is NOT proved.

## 1. What is actually constructed

Use a fixed Sylvester Hadamard node of order L=32, retain exactly31 of
the32 physical coordinates in every node group, and use independent
deep recursive child bases. One child carries a marked constant row;
all other children have the ordinary law. Dephasing cancels input signs
exactly, but leaves an unsigned permutation. Its groups are independent
of the displayed child bases. Independent non-DC output signed
permutations and fair outer edge signs remain available.

The resulting physical transform has a DC column and selected non-DC
column excess O_L(sqrt(m log m)). Exact balanced-column sign repair
has operator error O_L(m^(1/4) log^(5/4)m). In the rank-one block weave
this gives uniform energy error o(N^(3/2)), N=mq, q/m->31/32.
Thus the bulk has actual sign entries between fibres and exactly zero
block row and column sums. It is not a fractional projection passed
off as a signing. Its diagonal-fibre zeros are an explicitly separated
bulk operator, not a full minimax witness by themselves.

## 2. The new information inequality, independently derived

Let X_r lie within O(r) of a fixed finite alphabet Z in {0,b_1,...,b_d},
with P(Z!=0)=O(r). Nonzero coarse atoms are separated, so Z is a function
of X_r for small r. Write ell=log(1/r), t=tau ell, tau>0 fixed.
Then, uniformly over all atom probabilities,

    E_t(X_r)=-ell sum_j P(Z=b_j) min(tau b_j^2,1)+o(r ell).       (A)

Here E is the ORIGINAL envelope g_t(E Var(X|label))-I(X;label).
Since Var(X_r)=O(r), the Gaussian formula implies
`-tV<=g_t(V)<=-tV+2t^2V^2`, uniformly over all conditional channels.
Hence E_t=-J_t+O(r^2 ell^2), where J_t=inf[I+t E Var].

For an independent lower proof of J_t, put Y=E[X_r|label], choose
delta=ell^(-1/4), and set A={|Y|>delta}. Jensen gives E Y^2=O(r),
so P(A)<=O(r/delta^2). For every rare coarse atom b_j let
u_j=P(A|Z=b_j). Categorical data processing, followed by binary
data processing on A and dropping the nonnegative baseline term, gives

    I(X_r;label)>=ell sum_j p_j u_j-O(r log ell).

On A-complement the squared prediction error, restricted to rare atoms,
is at least sum_j p_j(1-u_j)(|b_j|-delta-O(r))_+^2. Add t times
this inequality and minimize each u_j in [0,1]. This proves the lower
bound J_t>=ell sum_j p_j min(1,tau b_j^2)-o(r ell), uniformly.
For the reverse, reveal only the coarse atoms with tau b_j^2>1 and
pool all others. The information cost is their total mass times ell
plus O(r); residual prediction error is the unrevealed second moment
plus O(r^2). This proves (A) without assuming a preferred latent channel.

## 3. The exact finite packet inequality

Let H be ANY dephased Hadamard of order L, S a nonempty proper subset
of columns, s=|S|, and C_j=sum_(l in S) H_jl. Then

    sum_j min(C_j^2/(L-1),1)>=1+s/(L-1).                     (B)

The DC coefficient is s. If it is not clipped, no coefficient is,
and Parseval gives Ls/(L-1)>=1+s/(L-1). Otherwise it pays1.
The non-DC energy is s(L-s), and their squared coefficients are at
most min(s,L-s)^2. Their clipped sum is consequently at least

    s(L-s)/max(L-1,min(s,L-s)^2)>=s/(L-1).

Both entries of the denominator are at most (L-1)(L-s). This proves
(B) with no classification of subsets or typical-spectrum assumption.
The proof is linear in s and therefore survives arbitrary packet mixtures.

For a node retaining k=L-1 coordinates, a rare minority packet S has
coarse outputs -2C_j/sqrt(k). Its physical entropy is paid ONCE, while
(A) charges all its child outputs. With the ONE common temperature
t=(1/4)log(1/r), (B) yields

    Gamma_t(1-2r)<=-[k/(L(L-1))] r log(1/r)+o(r log(1/r)).    (C)

The empirical group law is optimized before this bound, not assumed
product. Its entropy excess is at most P(S!=empty)log(1/r)+O_L(r),
and E|S|=kr. These identities justify every factor in (C).

## 4. Exact cap normalization and coverage

For the centered weave, the row variance is v=1-a^2 and the defect
identity is D_sigma/(2q)=m sum_i v_i-2sigma H_bulk/q.
Thus a common-temperature sufficient row condition for coefficient b is

    Gamma_t(a)+t(1-2b sqrt(p))(1-a^2)<0,  p=k/L.

Equation (C) gives the limiting rare coefficient sqrt(p)/2. This is
an UPPER bound for actual biased slices; it does not assert an actual
matching lower witness on those same slices.

For every fixed density-ratio bound K, sufficiently small fixed r>0,
and every positive coefficient margin, one may choose a single finite
child depth so that the actual construction controls ALL words whose
fibres are constant or have minority fraction in [r,Kr], with

    |H_bulk(x)| <= [sqrt(p)/2+o_(r->0)(1)] sqrt(N)
                         sum_i ||P_i x_i||^2 + o_N(N^(3/2)).

Signs of fibre means, active subsets and physical spin counts are all
included. For positive macroscopic variance the row margin beats the
subleading type/depth errors. For arbitrarily small total variance,
the deterministic operator bound pays an arbitrarily small additive
N^(3/2) error. There is no microscopic relative-error assertion.
The stronger log-density band in the construction artifact uses an
interval of clipping thresholds; it still does NOT cover arbitrary
widely separated density scales.

The limit order is: fix L, desired coefficient/additive margins and
bias band; choose sufficiently large FINITE child depth; then take
physical order to infinity. Shrinking margins/bias bands uses a diagonal
choice afterward. No growing-depth type enumeration is silently assumed.

## 5. Same-law balanced certificate: a subtle symmetry requirement

For Sylvester H32, translating input columns multiplies EACH FIXED
output row by a sign. The DC output is unchanged. Averaging an
unstratified input tuple law under translations therefore preserves
the raw DC law and each symmetrized non-DC law exactly, increases
entropy, and makes every input marginal equal to the prescribed source.
The general precision-Schur inequality now proves the old E bound at
this root. No global-reversal averaging of the raw marked law is needed.

Stratification costs epsilon=h(31/32)-log32/32 in the one-row exponent.
The independent directed interval replay, at60 decimal digits, proves
E_(97/20)(nu_(31/32))<=-4/5:29,895 boxes,14,948 leaves. Outward rounding
at65 digits then gives the SAME-law balanced-face coefficient

    [97/20+(31/32)log2-4/5+epsilon]
                 /[(97/10)sqrt(31/32)] < .497761196,

while sqrt(31/32)/2 < .492125493 for the rare bands. These are not
new global minimax upper constants: intermediate mixed profiles remain.

An attempted H24 shortcut is deliberately NOT promoted. Its exact
unsigned column automorphisms can permute output rows. Although this
preserves a sum before averaging, it may MIX different child laws;
E is not known to have the concavity needed to preserve that sum.
The H24 numerical constants and automorphism witness are preserved,
but they do not prove the same-law balanced inheritance. The rare
packet theorem for H24 is independently valid and does not need that step.

Reproducible checks: `computations/principle_director_stratified_selector_check_2026_09_07.py`
verifies69,896 proper subsets exhaustively at orders4,8,12,16, plus20,048
order24 test subsets, and the outward normalization intervals. These
checks supplement (A)--(C); they do not replace their proofs.

## 6. What this removes, and what it does not

This is a positive joint-response principle: physical packet entropy is
charged before its algebraically linked channels are separated. The
result beats the rigorously falsified unstratified scalar variance-ratio
certificate, whose optimized near-constant coefficient tends to
1/(2sqrt(p)). It controls actual sign bulk slices, not just a new state.

Still unproved are the complete intermediate/multiscale mixed response,
its compatibility with an arbitrary favorable seed, and any liminf-to-
all-order comparison for M_n. Neither this sub-half bulk theorem nor
the new necessary near-minimizer complexity law makes the original
convergence claim a proved consequence.
# Closing audit: stratified sparse-port realization

Root independently read Rudelson--Vershynin Lemma3.6 and its proof in
the author PDF, pages7--11 (printed numbering), on2026-09-07. Its input
is an arbitrary deterministic collection of bounded-coordinate vectors,
not an iid sampling law. For one omitted Hadamard row per L-group,
aggregate expected Gram is M I, M=m/L. Independent nonidentical
symmetrization and the deterministic lemma give

    D<=2a sqrt(D+M),  a<=C sqrt(s)log^2(m),
    D<=4a^2+2a sqrt(M).

Thus D/m=o(1) for s=o(m/log^4 m). Conditioning EACH physical selector
given its full frame on the restricted-Gram/excess good event costs a
uniform density factor at most2. This preserves prior nonnegative row
partition bounds up to factor2, independent fibres, and column-orbit
invariance; it does NOT claim conditioned child independence. Exact
balanced repair has global operator error o(sqrt(m)); orthogonal
projection can only decrease every restricted upper singular norm.
The signed reciprocal-swap bound then yields relative variance
coefficient sqrt(p)/2+o(1) for all words with at most s+1 active fibres,
simultaneously every outer seed and symmetric mask. This is a modified
actual law, not an unproved claim that the unconditioned law satisfies
all simultaneous events. Canonical proof and independent audit are the
synthesis sparse-active artifact Section6 and construct stratified-RIP
audit. No positive-density active theorem follows from this lemma.
