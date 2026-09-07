# Stratified marked selectors and a strict near-constant bulk response

2026-09-07. **Positive actual ensemble operation and finite block exponent;
Sections1--7 independently reconstructed PASS by both director and
construction agent.** Sections8--9 sharpen the response and state the
bounded-ratio heterogeneous-band consequence. This does not claim the
mixed response at EVERY bias, nor original minimax convergence.

## 1. Exact L-way marked Hadamard node

Fix a Hadamard order L, dephase its first row to all ones, and write its
sign matrix as H_L. At ambient order m=LM form the normalized Hadamard

    U=diag(U_0,...,U_(L-1)) (H_L/sqrt(L) tensor I_M) g,

where g is a uniform signed input permutation and the independent child
bases have order M. Mark an output row in child0 and dephase all physical
input coordinates by its signs. If r is that child's marked row sign
vector, the marked row signs before g are L copies of r. Exactly as in
the binary marked-column proof, the dephased matrix is

    Ubar=diag(Ubar_0,U_1 D_r,...,U_(L-1) D_r)
                        (H_L/sqrt(L) tensor I_M) Pi,      (1)

where Pi is a uniform UNSIGNED permutation. Conditional on child0, all
other displayed children retain their independent ordinary laws, because
their independent fresh signed input permutations absorb D_r. Thus they
are independent of the marked child and of Pi.

The canonical input groups are the L coordinates with a common index
in [M]. Pull them back by Pi to obtain a random partition of the physical
input coordinates into groups of size L. This partition depends only on
Pi, not on ANY of the displayed child matrices. All dependence statements
therefore survive conditioning the selector on this partition.

Choose exactly k coordinates uniformly in each group, independently
between groups, where 1<=k<L and p=k/L. The total physical selector size
is q=kM=p m. If P_V projects onto the marked child's preterminal input
subspace (the span of the normalized group indicators), then EXACTLY

    P_V D_T P_V=p P_V.

For every vector v supported on the retained physical coordinates,

    ||P_V v||^2<=p||v||^2.                         (2)

This is a deterministic all-vector statement, not a typical-word claim.

## 2. Entropy cost and invariances

Relative to a uniform q-subset of [m], the stratification event has
probability

    binom(L,k)^M / binom(m,q),

so its negative logarithm per ambient coordinate tends to

    epsilon_(L,k)=h(p)-L^(-1)log binom(L,k).        (3)

For p bounded away from0 and1 this is O(log L/L); the same tends to zero
for k=L-1. This cost is fixed before m tends to infinity. Conditioning a
nonnegative one-row partition function costs at most this exponent.
Thus, WHEN the unconditioned root inherits the old certificate as proved
under the symmetry hypothesis in Section10, a sufficiently large fixed L
preserves that strict certificate whenever its margin exceeds
epsilon_(L,k), with the usual small continuity payment if p is adjusted
to k/L. Arbitrary L-way roots are not assumed to have this inheritance.

Independent output signed permutations of the non-DC columns remain
available: the stratification event depends on the physical input
partition, not those output permutations or any edge sign. Different
fibres use independent node/selector samples. The marked DC column can
still be aligned with its self-port and omitted, at the standard
polynomial orbital-deletion cost.

Column excesses remain small enough for exact balanced-column repair.
Each non-DC full Hadamard column has total sum zero. Its selected sum
is a sum over independent groups of uniformly sampled k-subsets, with
mean p times that total, hence zero, and group summands bounded by L.
For fixed L, Hoeffding and a union bound give O_L(sqrt(m log m)) excess
uniformly over all relevant columns and fibres. This is o(m/log^2 m),
so the repair theorem applies. Conditioning and repair preserve the
old balanced-face certificate up to the stated entropy payment and an
o(N^(3/2)) energy error.

## 3. An explicit constrained finite block exponent

Let W be the FINITE set of pairs w=(T,x), where T is a k-subset of [L]
and x is a sign assignment on T. For a desired physical fibre mean a,
define the L scalar block outputs

    y_j(w)=sum_l H_L(j,l)1_T(l)(x_l-a)/sqrt(k),
                              0<=j<L.             (4)

Let mu range over probability laws on W with

    E_mu sum_(l in T) x_l=k a.

The marked output j=0 has mean zero. Write nu_0=law_mu(y_0), and for
j>0 write nu_j=sym law_mu(y_j). Define the explicit block functional

    Gamma_(t,L,k)(a)=L^(-1) sup_mu
       {H(mu)-log binom(L,k)+sum_(j=0)^(L-1) E_t(nu_j)}.  (5)

The E_t here is the ORIGINAL envelope g_t(E Var)-I, not E[g_t(Var)]-I.
The law mu is a finite alphabet law on specified physical group words,
not an arbitrary feedback state.

Equation (5) is an actual one-row upper exponent. Indeed each sequence
of M group words has selector weight binom(L,k)^(-M). Its empirical
type has count exp(M H(mu)+O(log M)). The L child coordinate arrays
are exactly the linear images (4). The ordinary children have the old
deep exponent E_t; the marked child has the proved marked-spine upper
exponent E_t. Children are independent by (1), and Gaussian orbital
norms obey the concatenation product bound. Choosing all child depths
sufficiently large but FIXED and then M tending to infinity yields (5)
up to any prescribed exponent tolerance. Source alphabets are finite
and uniformly bounded for fixed L, so these estimates are uniform in a.

## 4. Uniform rare-atom information lemma

Fix a finite set of nonzero real numbers {b_1,...,b_d}, a constant K,
and tau>0. Let r decrease to zero and suppose Z takes values in
{0,b_1,...,b_d}, with P(Z!=0)<=K r. Suppose also |X-Z|<=K r, with the
coarse label Z determined by X (disjoint O(r) neighborhoods suffice).
Put ell=log(1/r) and t=tau ell. Then uniformly over all such laws,

    E_t(X)=-ell sum_(b!=0) P(Z=b) min(tau b^2,1)
                         +o(r ell).               (6)

Here the possible O(r) within-cluster fluctuations may have arbitrary
laws. In the present finite group application they are finite.

Proof of the needed upper bound. Let J_t=inf_L[I(X;L)+t E Var(X|L)].
Since E X^2=O(r), the Gaussian formula and 0<=rho<=2tV give uniformly

    -tV<=g_t(V)<=-tV+2t^2 V^2,
    E_t(X)=-J_t(X)+o(r ell).                       (7)

For ANY channel L let q(L)=P(Z!=0|L), delta=ell^(-1/2), and call a
label high if q>delta. The pointwise binary relative entropy is
nonnegative. On high labels,

    D(Ber(q)||Ber(P(Z!=0)))
       >=q ell-O(q log ell+q),

uniformly, since h(q)<=q(log(1/delta)+1). Dropping the nonnegative
low-label contribution gives

    I(X;L)>=P(Z!=0,high) ell-O(r log ell).          (8)

For low labels, |E[Z|L]|^2<=max_b b^2 q^2. Thus their contribution to
E Var(Z|L) is at least

    sum_b b^2 P(Z=b,low)-O(delta r).

The inequality Var(X|L)>=(1-delta)Var(Z|L)
-delta^(-1)E[(X-Z)^2|L] transfers this to X with an additional
O(delta r+r^2/delta) error. Discard the nonnegative high-label variance.
Adding t times this bound to (8) proves

    I+t E Var(X|L)
       >=ell sum_b P(Z=b) min(1,tau b^2)-o(r ell).

The error is uniform over channels and laws, proving the required
direction of (6). For the reverse direction reveal exactly those coarse
atoms with tau b^2>1 and pool the others. Their entropy cost is
ell times their total mass plus O(r), uniformly over finite atom
probabilities, and the residual mean squared prediction error is the
sum of b^2 masses of unrevealed atoms plus O(r^2). This proves (6).

## 5. Rare packet mixture expansion

Put a=1-2r. For w=(T,x), let S be the subset of T carrying minus signs,
and s=|S|. The mean constraint is E_mu s=k r. The outputs (4) are

    y_j=2[r sum_(l in T)H_L(j,l)-C_j(S)]/sqrt(k),
    C_j(S)=sum_(l in S)H_L(j,l).                   (9)

They differ by O_L(r) from the finite coarse values -2 C_j(S)/sqrt(k).
Moreover P(S!=empty)<=k r. Symmetrization of an ordinary child only
splits the signs of its nonzero coarse atoms; formula (6), which depends
on their squares, is unchanged. Therefore for t=tau log(1/r), (6) is
uniform in EVERY admissible mu.

The selector/spin entropy satisfies

    H(mu)-log binom(L,k)<=H(S).

Indeed conditional on S there are at most binom(L-s,k-s) choices of T,
never more than binom(L,k). If theta=P(S!=empty), then
r<=theta<=k r and

    H(S)<=theta log(1/r)+O_L(r).

Substituting the rare lemma into (5) proves the UNIFORM mixture bound

    Gamma_(t,L,k)(1-2r)
       <=[log(1/r)/L] E_mu {1_(S!=empty)
               -sum_j min(4tau C_j(S)^2/k,1)}+o(r log(1/r)), (10)

where the right side is maximized over admissible mu. There is ONE
common tau, not a packet-dependent temperature or an exchanged minimax.

## 6. A universal Hadamard clipping inequality

For every nonempty proper subset S of [L], put s=|S|. Then

    sum_j min((L+1)C_j(S)^2/L^2,1)>=1+s/L.        (11)

To prove this, C_0=s and |C_j|<=s. If the DC term is not clipped, none
is clipped, and Parseval gives a sum (1+1/L)s>=1+s/L. Otherwise the DC
term contributes1. For the remaining terms,

    sum_(j>0) C_j^2=s(L-s),
    max_(j>0) C_j^2<=min(s,L-s)^2,

the latter using zero sums of the non-DC Hadamard rows. Hence their
clipped sum is at least

    s(L-s)/max {L^2/(L+1), min(s,L-s)^2}>=s/L.

The last inequality follows because L(L-s) dominates BOTH entries of
the maximum when1<=s<=L-1. This proves (11) without a classification
of packet supports or any randomness assumption.

Take

    tau=k(L+1)/(4L^2).

Then (10)--(11) and E s=k r give

    Gamma_(tau log(1/r),L,k)(1-2r)
       <=-(k/L^2) r log(1/r)+o(r log(1/r)).        (12)

This is a strictly negative entropy coefficient uniform over all packet
mixtures and all selector/spin correlations admitted by the construction.

## 7. Actual strict near-constant variance response

The exact centered-weave profile normalization is the one in the marked
column artifact: a row with variance v=1-a^2 contributes the Markov term
t(1-2b sqrt(p))v. Since v=4r(1-r), (12) implies

    limsup_(r->0) [1+Gamma_(tau log(1/r),L,k)(1-2r)
                         /(tau log(1/r) 4r(1-r))]
                    /(2 sqrt(p))
       <= L/[2(L+1) sqrt(p)].                    (13)

For k=L-1 this is STRICTLY below1/2. More generally strictness holds
when p>(L/(L+1))^2. This differs from the unstratified scalar-envelope
floor 1/(2 sqrt(p)): the physical rare packet pays its entropy once,
but the Hadamard node spreads it into multiple separately penalized
child channels. Equation (11) makes that gain uniform under mixtures.

The statement (13) first gives a homogeneous biased-slice coefficient.
Section9 supplies the common-temperature consequence for bounded-ratio
near-constant bands. Arbitrarily many logarithmically separated density
scales still require a further argument.

For each fixed small positive r and any b strictly above the right side
of (13), the finite-depth row certificate has a strict margin for all
fibres with that bias (and its sign reversal). The standard product
Finner bound then gives an ACTUAL ensemble cap on those slices, and
operator-norm balanced repair transfers it with additive o(N^(3/2)).
All depths and L are fixed before the ambient order. This is a positive
new part of the missing mixed response, not a complete all-profile seed
recovery theorem.

## 8. Sharper clipping reaches the matching-scale coefficient

The following strengthening improves the response in (13):

    sum_j min(C_j(S)^2/(L-1),1)>=1+s/(L-1),       (14)

for EVERY nonempty proper S. If the DC term is not clipped, none is,
and the sum is Ls/(L-1)>=1+s/(L-1). Otherwise the non-DC sum is at least

    s(L-s)/max {L-1,min(s,L-s)^2}>=s/(L-1).

Indeed (L-1)(L-s) dominates L-1, and dominates min(s,L-s)^2 because
L-s<=L-1 and min(s,L-s)<=L-s. Adding the clipped DC term proves (14).

Take now the common temperature scale

    tau_*=k/[4(L-1)].

Equations (10) and (14) improve (12) to

    Gamma_(tau_* log(1/r),L,k)(1-2r)
       <=-[k/(L(L-1))] r log(1/r)+o(r log(1/r)).   (15)

Its variance-response coefficient is

    b_*(L,k)=(1-1/L)/(2 sqrt(p)).                  (16)

At k=L-1 this is sqrt(p)/2<1/2, with tau_*=1/4. This is an UPPER
coefficient at the familiar matched-column scale, not an assertion that
a matching lower witness has been constructed on every homogeneous
near-constant slice. The weaker bound (13) is retained as its independently
checked precursor.

For the concrete choice L=24,k=23,

    p=23/24, tau_*=1/4,
    b_*(24,23)=sqrt(23/24)/2 approximately .48947250518628047.

The director has also completed an independent directed-interval
certificate E_(4.85)(nu_(23/24,0))<=-41/50. Paying the exact selector fee
(3) would give a SAME-ENSEMBLE balanced-face coefficient below .498649,
provided the L24 source-envelope inheritance is separately verified.
Its complete interval record is
`../computations/results/principle_director_stratified_balanced_certificate_2026_09_07.json`.
The scalar interval certificate alone must not be confused with that
ensemble-symmetry obligation. The analogous inheritance for Sylvester
L is proved in Section10. Once both ingredients are present, the actual
stratified and repaired ensemble has strict-subhalf certificates on
BOTH the all-balanced face and sufficiently near-constant bands below.
This does not cover mixed collections of balanced and near-constant
fibres, or the intervening biases.

## 9. Uniform heterogeneous bounded-ratio bands

Fix K>=1. For a small fixed r>0, allow each physical fibre independently
to be either EXACTLY constant, or to have minority fraction

    r_i=(1-|a_i|)/2 in [r,K r].

There is no restriction on the sign of its mean or on the set of active
fibres. Choose ONE common t=tau_* log(1/r). The rare-information proof
is uniform with this reference r: each active row has rare packet mass
at most kK r, atom displacement O_L(K r), and E s=k r_i. The entropy
bound is

    H(S)<=P(S!=empty) log(1/r)+O_(L,K)(r).

Thus (14) gives, uniformly for every active row,

    Gamma_t(a_i)<=-[k/(L(L-1))] r_i log(1/r)
                                      +o_(L,K)(r log(1/r)).

Since r_i>=r, the error is also o(r_i log(1/r)). On constant rows the
exact partition contribution has zero leading exponent and zero variance;
one may count their two sign choices separately, at total exp(O(m)) cost.
Summing the row bounds therefore gives the same coefficient (16) for
ALL such heterogeneous profiles at once. Their total variance is
sum_i q*4r_i(1-r_i).

More precisely, for every b>b_*(L,k), all sufficiently small fixed r
admit a sufficiently deep but FIXED child ensemble whose repaired
balanced bulk satisfies, with asymptotically high construction probability,

    |H_bulk(x)|<=b sqrt(N) sum_i||P_i x_i||^2+o(N^(3/2))   (17)

uniformly over this entire band. To see the uniform union bound, first
restrict total variance to at least delta N for fixed delta>0. The strict
row-exponent margin is then order m^2, while the choices of integer fibre
means and active sets have log count O(m log m). The usual product Finner
estimate applies with the SAME temperature on every edge. For total
variance below delta N, the deterministic operator bound
||B_bulk||op<=m+o(m) gives an absolute O(delta N^(3/2)) bound. Send the
order to infinity first, then delta to zero. Balanced repair costs
o(N^(3/2)) uniformly in all physical words. This proves the additive
uniform formulation, without claiming relative control at microscopic
variance.

In particular, whenever the same-law balanced certificate is validated
as in Section10, that construction controls the balanced face and any
sufficiently small fixed bounded-ratio band with a common strict-subhalf
upper coefficient. Simultaneous control of
arbitrarily many widely separated minority scales, or of arbitrary
intermediate fibre means, remains OPEN. No convergence conclusion is
being inferred from this restricted but actual positive response theorem.

## 10. Same-law balanced inheritance requires an explicit symmetry

The finite constrained exponent (5) is valid for ANY dephased H_L.
However the old scalar E_t(nu_p) upper bound for the UNCONDITIONED L-way
root is not automatically valid for every H_L. Its input tuple law has
only an average marginal constraint, while the general precision-Schur
inequality has one envelope per individual marginal. It would be invalid
to replace those by E_t of the average law without further reasoning.

A sufficient condition is a transitive group of UNSIGNED input-column
permutations P such that H_L P=D H_L, where D is DIAGONAL signs, fixing
the marked DC row with positive sign. For the Sylvester translation
group, every non-DC character takes both signs equally often. Average a
tuple policy under this group ALONE. The raw DC law is unchanged; every
other output becomes its already-used symmetrized law. Relative entropy
decreases, and transitivity makes EVERY input marginal exactly nu_p.
There is no preliminary global-reversal average of the raw DC law.
The arbitrary orthogonal precision-Schur theorem now gives

    sum_j E_t(output_j)-D(tuple||nu_p tensor...tensor nu_p)
                                <=L E_t(nu_p).

Thus the unconditioned one-row exponent is at most p log2+E_t(nu_p),
and conditioning pays exactly (3). This proves the needed inheritance.

For a Sylvester H_L, L=2^d, the hypothesis is elementary: translating
input columns by a vector in F_2^d multiplies each output character row
by its sign. The DC row is fixed and the translation group is transitive.
This gives a fully structural safe family without any automorphism
classification. For L24=H12 tensor H2, the construction agent has verified
an explicit transitive unsigned-column automorphism action. However its
OUTPUT action also permutes rows. This does not by itself suffice:
averaging the tuple then mixes different child laws, and convexity of E
allows their summed envelope value to decrease. Merely saying that every
individual group element preserves the summed value would miss this
mixture issue. The exact automorphism is useful data, but the current
safe same-law inheritance theorem is the Sylvester case. Neither a full
monomial input-sign action nor an output-permutation-only symmetry is
being silently substituted for the proved diagonal-output condition.

## 11. Polynomially broad minority-density bands at ONE temperature

There is a useful strengthening of Section9. For EVERY

    L/2<=c<=L-1,

the clipping inequality holds in the form

    sum_j min(C_j(S)^2/c,1)>=1+s/c.               (18)

If the DC term is not clipped, the sum is Ls/c and (L-1)s>=c. If it
is clipped, the non-DC lower bound is
s(L-s)/max(c,min(s,L-s)^2). This is at least s/c because

    min(s,L-s)^2/(L-s)<=L/2<=c.

Consequently the same coefficient b_*(L,k) holds uniformly for any
temperature-to-log-density ratio in the entire interval

    k/[4(L-1)]<=t/log(1/r_i)<=k/(2L).             (19)

The rare-information lemma and its packet entropy estimates are uniform
when that ratio ranges over this compact positive interval.

Fix any alpha with

    1<=alpha<=2(L-1)/L.

For small fixed r, allow every fibre to be either exactly constant or
to have minority density anywhere in

    r^alpha<=r_i<=r.

Choose the SINGLE common temperature t=[k/(2L)]log(1/r). Then (19) is
satisfied for every active row. With c_i=k log(1/r_i)/(4t), (18) gives

    Gamma_t(a_i)<=-(4t/L)r_i+o(t r_i)

uniformly over all these densities and both signs of the means. The
same product Finner, fixed-depth, and low-total-variance argument from
Section9 proves (17) at every b>b_*(L,k) for this MUCH BROADER band.

For L32, one can take any alpha up to31/16. Thus a single actual ensemble
controls minority scales from r to nearly r^2 at a common strict-subhalf
variance coefficient. This is not an interchange of temperatures inside
an edge kernel: the temperature is literally the same on all edges.
Arbitrarily large ratios of log minority densities remain outside this
proved statement.

## 12. A concrete audited same-matrix theorem at L32

The safe Sylvester choice is now fully certified, with no H24 symmetry
assumption. Set L=32,k=31,p=31/32 and take t_0=97/20. The director's
60-digit directed-interval certificate proves E_(t_0)(nu_p)<=-4/5.
After the exact stratification fee

    epsilon=h(31/32)-log32/32,

65-digit outward evaluation proves

    [t_0+p log2-4/5+epsilon]/[2t_0 sqrt(p)]
             < .497761196,
    sqrt(p)/2 < .492125493.                       (20)

The source certificate is
`../computations/results/principle_director_sylvester_stratified_certificate_2026_09_07.json`
(29,895 checked boxes,14,948 accepted leaves). The independently evaluated
normalization intervals are in
`../computations/results/principle_director_stratified_selector_check_2026_09_07.json`.
The director's complete mathematical reconstruction is
`principle_director_joint_packet_audit_2026_09_07.md`.

The construction agent sharpened (18) to the exact one-hole tradeoff

    min_S [sum_j min(C_j^2/c,1)-1]/|S|
                    =min(1/c,2/L,L/c-1),  1<=c<=L.

Its proof and attaining singleton/complement/half-column packets are in
`principle_construct_2026_09_07_rare_packet_temperature_curve.md`.
For a chosen b between sqrt(p)/2 and1/2, put delta=1-2b sqrt(p).
The actual band theorem extends to any

    alpha<2(1-delta)/(L delta).                   (21)

This uses ONE common temperature over the band, with strict slack in
the indicated c interval. Its proof is the same uniform rare lemma and
finite-depth union argument, now using the exact clipping curve.

In particular b=499/1000 permits alpha=3 at L32. Indeed (21) reduces
to delta<1/49, equivalent to

    sqrt(31/32)>24000/24451,

which follows by squaring positive rationals. For every sufficiently
small fixed r>0 there are arbitrarily large actual repaired sign bulks
B from the SAME ensemble with BOTH properties:

1. Every word balanced in every fibre has |H_B(x)|<(.499+o(1))N^(3/2).
2. Every word whose fibres are constant or have minority fractions in
   [r^3,r] has

       |H_B(x)|<=(.499) sqrt(N) sum_i||P_i x_i||^2
                                             +o(N^(3/2)).

Signs of the fibre means and all choices of active fibres are included.
These are simultaneous properties of one actual matrix sequence: choose
a common sufficiently large finite child depth for the balanced and band
temperatures, apply both moment bounds, and intersect their high-probability
events with the uniform repair event. Exact within-fibre balance and full
cross-fibre sign entries are preserved.

The word sets in properties1 and2 are separate regions. A single word
mixing balanced fibres with rare fibres, or arbitrary intermediate means,
is NOT asserted to satisfy the displayed variance bound. The bulk may
be filled by constant-mode blocks to obtain a full signing, but no cap
for that full completion is inferred here. This concrete two-region
theorem is the positive checkpoint; its remaining mixed-profile gap is
not hidden by the numerical constants.

## 13. An explicit mixed cone tolerating arbitrary exceptional fibres

The director independently reconstructed the following exponent budget.

There is a small but genuine extension beyond separate word regions.
Use L32,k31, b=.499, and the sharp-coefficient band

    r^(3/2)<=r_i<=r,
    t=[k/(2L)]log(1/r).

Put delta=1-2b sqrt(p) and d=1/L-delta>0. For sufficiently small fixed r,
the sharp clipping plateau and the uniform rare remainder give every
good active row the strict weighted exponent bound

    Gamma_t(a_i)+t delta(1-a_i^2)<=-2t d r_i.      (22)

Exactly constant rows have zero variance and zero leading exponent.
Now allow a set B of exceptional fibres whose words are COMPLETELY
arbitrary, including balanced or intermediate-bias words. Since every
orbital norm is at most1, an exceptional row's entire spin partition has
exponent at most p log2. Its Markov variance term is at most t delta.
For small r, t>=p log2, so its total contribution is at most t(1+delta).

Consequently the same common-temperature proof still works whenever

    |B|<=[d/(1+delta)] sum_(good active i) r_i.     (23)

Indeed the total exceptional exponent is then at most t d sum r_i,
while (22) contributes at most -2t d sum r_i. The remaining strict
margin pays all fixed-depth and type errors. Summing over the possible
exceptional sets costs only exp(O(m)), and the low-total-variance case
is handled by the same uniform additive operator estimate.

Thus the actual variance bound with b=.499 holds on this explicitly
defined MIXED cone, not only when every active fibre is in the same
near-constant band. Condition (23) is quantitative: the exceptional
count is controlled by the total minority mass of the good fibres.
There is no assertion that an arbitrary full-sign seed forces its
maximizing words into this cone, and unrestricted mixed-profile recovery
remains open.
