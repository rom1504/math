# Quantitative hierarchy compression and actual seed-preserving dilution

2026-09-07. **PASS, including Section7's sharper constant; independently
reconstructed in `principle_construct_2026_09_07_hierarchy_complexity_audit.md`.** This strengthens
the earlier zero-complexity-center and fixed power-barrier criteria. It
uses the FULL absolute near-extreme superlevels, permits positive
normalized factorization complexity at every fixed level, and allocates
the balancing budget quantitatively across their hierarchy.

The operation preserves every old edge and fills every new edge by a
sign. No convergence assertion is made. It supplies a new necessary
small-energy complexity condition for sequences realizing the liminf.

## 1. Statement and normalizations

For a matrix of sign rows indexed by a nonempty code E subset {+/-1}^n,
write gamma(E) for its gamma2 factorization norm: E=LV with each column
of V of Euclidean norm at most one and each row of L of norm at most
gamma(E). Always1<=gamma(E)^2<=n.

Let A_n be any sequence of hollow full signings with
Q(A_n)/n^(3/2)->c<infinity. Define

    E_n(eta)={x:Q(A_n)-|H_A(x)|<=eta n^(3/2)},
    e(eta)=limsup_n gamma(E_n(eta))^2/n.

Let h be binary entropy with natural logarithms, and use its monotone
cap

    H(g)=h(min(g,1/2)),  0<=g<=1.

The cap is essential: h(g) itself decreases above1/2 and would give a
spurious zero at g=1. Set

    K=limsup_(eta down to0) e(eta) H(e(eta))/eta.     (1)

If K<infinity, then for EVERY tau>4096K and every sufficiently small
fixed epsilon>0 there are actual full-sign extensions B_(n+q), with
q=floor(epsilon n), whose old principal block is exactly A_n and

    limsup_n [Q(B_(n+q))-Q(A_n)]/n^(3/2)
        <=tau epsilon+epsilon^(3/2).                (2)

In particular, when K=0 the incremental cost can be made o(epsilon)
in the ordered limits n first, then epsilon down to0. The constant4096
is deliberately safe, not optimized.

Let c_* =liminf_n M_n/n^(3/2)>0. For ANY sequence of actual signings
whose normalized caps tend to c_*, necessarily

    limsup_(eta down to0)
       e(eta)H(e(eta))/eta >=3c_*/8192.              (3)

Thus, along arbitrarily small energy windows, their full-superlevel
gamma2-square density must reach the scale

    sqrt(eta/log(e/eta)),

up to a positive constant depending on c_*. This is substantially
stronger than the generic linear-in-eta lower bound coming from a
Boolean face near any bounded-cap maximum.

## 2. The exact finite hierarchy resource

Fix levels eta_j=eta_0 2^(-j), j=0,...,J. Suppose for a particular n
we have factorizations of E_n(eta_j) with

    gamma(E_n(eta_j))^2<=g_j n,  0<g_j<=1/2.

For fixed epsilon,tau>0 define weights

    a_j^2 =160 epsilon g_j[H(g_j)+epsilon log2]
            /(eta_j/2+tau epsilon)^2,  0<=j<J,

    a_J^2 =160 epsilon g_J[H(g_J)+epsilon log2]
            /(tau epsilon)^2.                     (4)

If sum_j a_j^2<=1/2, concatenate all the factorization columns with
these weights and an identity component of weight1/sqrt2:

    v_i=(e_i/sqrt2, a_0 V_(0),i,...,a_J V_(J),i).

Every v_i has norm at most one. Apply the Gram--Schmidt Walk, started
at zero, to these n vectors. Its output is an ACTUAL sign vector z,
and the primary subgaussian statement gives, for every test vector theta,

    E exp(t<theta,sum_i z_i v_i>)
       <=exp(40t^2||theta||^2/2).                  (5)

The imported normalization is Theorem1.4 and the MGF definition in
[Bansal--Dadush--Garg--Lovett](https://theoryofcomputing.org/articles/v015a021/v015a021.pdf),
read and independently audited in the balanced-augmentation artifacts.
The identity feature is retained, not traded away to balance the codes.

For x in E_n(eta_j), testing only the j-th factorization block gives
a subgaussian variance proxy40g_j n/a_j^2 for z dot x. Testing only
the identity component gives the universal proxy80n for EVERY x.
Take q independent walk outputs as the q columns of a new sign matrix C.
For fixed old and new words x,y, the bridge x^T C y consequently has
variance proxy40qg_j n/a_j^2 on level j and80qn universally.

The factorization-to-entropy estimate used here is exact in scope:

    VCdim(E)<=gamma(E)^2,
    log|E|<=n H(g)+o(n) if gamma(E)^2<=gn.           (6)

The first follows by testing all labels on a shattered set against its
factorization and averaging their squared signed-column sum. The second
is Sauer's recursion and the elementary binomial entropy bound. These
were proved in the balanced-augmentation artifact; no entropy estimate
for a continuum of fields is substituted for the finite code.

Use the union thresholds

    T_j^2=(80q g_j n/a_j^2)
             log[8(J+1)|E_n(eta_j)|2^q],

    U^2=160qn log[8*2^(n+q)].                      (7)

The two-sided subgaussian tails and a union bound show that, with
probability at least1/2, ALL level bounds |x^T C y|<=T_j and the
universal bound |x^T C y|<=U hold simultaneously for all relevant x,y.
Indeed each level costs at most1/[4(J+1)] probability and the universal
bound costs at most1/4. Fix any successful realization.

For fixed parameters as n tends to infinity, (4), (6), and (7) yield

    limsup T_j/n^(3/2)
       <=(eta_j/2+tau epsilon)/sqrt2  (j<J),
    limsup T_J/n^(3/2)<=tau epsilon/sqrt2,
    U/n^(3/2)->sqrt[160epsilon(1+epsilon)log2].       (8)

The factor160 in (4) leaves slack for every finite-n entropy and union
error. In particular all the same limits hold with the larger right
sides lacking1/sqrt2.

If x lies in E_n(eta_j) but not E_n(eta_(j+1)), its old energy deficit
is greater than eta_j/2. The corresponding bridge in (8) therefore
raises its absolute energy to at most Q(A_n)+tau epsilon n^(3/2)+o(n^(3/2)).
The innermost level has bridge at most tau epsilon n^(3/2)+o(n^(3/2))
without needing a positive old deficit. Outside E_n(eta_0), the
universal bound cannot win once epsilon is small enough that the last
quantity in (8) is less than eta_0.

Finally put any sign matrix of order q with cap at most q^(3/2) in the
new principal block. Such a signing exists for all large q by the
elementary independent-edge sign union bound. For every parent word,

    |H_parent|<=|H_old|+|bridge|+Q(new block).

This pays the full absolute bridge, respecting independent reversal of
all new spins. It proves (2), provided the finite resource condition
sum a_j^2<=1/2 can be met. No old edge was altered.

## 3. The hierarchy costs O(K/tau), not a power of its depth

Choose kappa with K<kappa<tau/4096. By (1), for all sufficiently small
eta one can choose g(eta)>e(eta), still at most1/2, with

    g(eta) H(g(eta))<=kappa eta.                    (9)

For each fixed finite hierarchy, the limsup defining e(eta) then makes
its factorizations valid for all sufficiently large n. This is the only
order of quantifiers used; no uniform convergence at n-dependent tiny
energy levels is assumed.

For g<=1/2 the entropy chord gives H(g)>=2(log2)g. Hence (9) also gives

    g(eta)<=sqrt[kappa eta/(2log2)].                (10)

First sum the part of (4) containing g_j H(g_j). Put
u_j=eta_j/(2tau epsilon). Then

    sum_(j<J) 160epsilon kappa eta_j
                    /(eta_j/2+tau epsilon)^2
       =(320kappa/tau)sum_(j<J) u_j/(1+u_j)^2
       <=1280kappa/tau.                            (11)

The last bound is elementary and uniform in the depth: split the
dyadic sequence at u=1. On the large side bound by1/u and sum a
geometric series; on the small side bound by u and do the same.
Their total is at most four. This adaptive allocation is what avoids
an artificial logarithmic-depth loss.

The remaining epsilon*g_j part is negligible. By (10), after the same
substitution it is bounded by a constant times

    sqrt(kappa epsilon)/tau^(3/2)
       sum_j sqrt(u_j)/(1+u_j)^2
       =O(sqrt(kappa epsilon)/tau^(3/2)).           (12)

The dyadic sum is at most five: its two geometric ratios are2^(-3/2)
and2^(-1/2). This tends to zero for fixed kappa,tau.

Choose J so that eta_J is between epsilon^2/2 and epsilon^2, after
shrinking epsilon below eta_0. The special innermost weight in (4)
is, by (9)--(10),

    O(kappa eta_J/(tau^2 epsilon)+g_J/tau^2)=o(1).  (13)

Combining (11)--(13),

    sum_j a_j^2<=1280kappa/tau+o(1)<5/16+o(1)<1/2.

This proves the finite resource premise for all sufficiently small
fixed epsilon, and completes the extension theorem. It explicitly
accounts for the accumulated hierarchy, terminal, entropy, and universal
tail errors.

## 4. Consequence for the actual liminf

Suppose the sequence realizes c_*, but the left side of (3) is smaller
than3c_*/8192. Choose tau with

    4096K<tau<3c_*/2.

The actual parents from (2) then have, for every sufficiently small
fixed epsilon,

    limsup_n Q(B_(n+q))/(n+q)^(3/2)
       <=[c_*+tau epsilon+epsilon^(3/2)]/(1+epsilon)^(3/2)
        <c_*.

This contradicts the defining global liminf, along the parent orders.
Thus (3) holds. It concerns every liminf-realizing sequence, not every
exact-minimizer order in the absence of a convergence theorem.

For small g, gH(g) is comparable to g^2 log(e/g). Inverting this
monotone function yields the announced square-root-over-log scale
along arbitrarily small eta. The limsup and this sequence-of-scales
qualification are important; the theorem does not assert the bound
at every sufficiently small energy window.

## 5. Explicit center covers and compression

The theorem also gives a quantitative cover criterion, without requiring
the center complexity to be o(n) at each fixed scale. If E_n(eta) is
covered at Hamming radius r_n n by a center code F_n with
gamma(F_n)^2<=g_n n, then

    gamma(E_n(eta))/sqrt(n)<=sqrt(g_n)+2sqrt(r_n).   (14)

To see this, choose a center for each row x and write x=f(x)+(x-f(x)).
The repeated-center matrix has gamma at most gamma(F_n). The residual
matrix factors through the identity with row norm at most2sqrt(r_n n).
The direct-sum factorization, with optimized column weights, proves
subadditivity and hence (14).

Thus one may replace e(eta) in the sufficient hypothesis by any upper
bound of the form

    (sqrt(g_center(eta))+2sqrt(r_cover(eta)))^2.

This permits nonzero complexity densities that decay with the energy
window, and accumulates the hierarchy with the exact resource (4).
It includes the earlier power and critical-cover constructions but is
not limited to a barrier about actual grounds or a preselected fixed
center family.

For example, if the full superlevels themselves satisfy
e(eta)=O(eta^s) for any fixed s>1/2, then K=0 and the actual extension
has o(epsilon) cost. More sharply the same follows from
e(eta)=o(sqrt(eta/log(e/eta))). Such behavior is impossible along a
liminf-realizing sequence by (3).

## 6. A concrete check on the near-half critical families

In `principle_invent_2026_09_07_near_half_critical_landscape.md`, all
fibre-constant words form a center code of gamma2 square at most m,
while n=mq and q tends to infinity. The linear transverse barrier there
places E_n(eta) within O_c(eta)n of that code. Equation (14) gives
e(eta)=O_c(eta), so K=0. The present theorem really does dilute those
actual critical two-ground families after enlarging and compressing the
center code. Their two-ground critical barrier is not mistaken for an
intrinsic obstruction to adaptive balancing.

## 7. Sharper exact resource constant via a telescoping hierarchy

The baseline4096 above is independently audited and can be improved
without a new probabilistic input. Let phi=(1+sqrt5)/2 and

    A_*=20phi^5=221.803398874989... .

Then the extension theorem holds for every tau>A_* K. Consequently
the sharper necessary condition for liminf-realizing signings is

    limsup_(eta down to0) e(eta)H(e(eta))/eta
        >=3c_*/(40phi^5).                          (15)

Here are the full changes to the certified baseline argument. Use a
geometric hierarchy eta_j=eta_0 R^(-j), R>1, whose shell deficit is at
least eta_j/R. Retain an identity feature of squared weight theta>0,
and replace160 in (4) by80(1+zeta), with zeta>0. The code weights are

    a_j^2=80(1+zeta)epsilon g_j[H(g_j)+epsilon log2]
                  /(eta_j/R+tau epsilon)^2,

with the same terminal denominator tau epsilon. The asymptotic union
threshold has slack1/sqrt(1+zeta), while the universal identity test
has variance proxy40n/theta. Since theta is fixed before epsilon tends
to zero, the universal bridge is still below the fixed eta_0 for small
epsilon. Thus it suffices that sum a_j^2<1-theta.

Put r=sqrt(R), u_j=eta_j/(R tau epsilon), and f(u)=u/(1+u)^2. Direct
algebra gives the pointwise telescoping majorant

    f(u) <= (r+1)/(4(r-1))
               [1/(1+u/r)-1/(1+ru)].              (16)

For clarity, the ratio of the left side to the bracket equals

    [1+(r+1/r)u+u^2]/[(r-1/r)(1+u)^2],

whose maximum occurs at u=1 and is (r+1)/(4(r-1)). Since
ru_(j+1)=u_j/r, the brackets in (16) telescope and their sum is at
most one. The entropy part of the total resource is therefore at most

    80(1+zeta)(kappa/tau)
           R(r+1)/(4(r-1)).                       (17)

The function r^2(r+1)/(r-1) is minimized for r>1 at the golden ratio:
its logarithmic derivative is2/r+1/(r+1)-1/(r-1), which vanishes exactly
when r^2-r-1=0. At R=phi^2, the constant in (17), before1+zeta, is

    80 phi^2(phi+1)/(4(phi-1))=20phi^5=A_*.

The epsilon*g part still tends to zero by the same square-root estimate
as (12), now using geometric ratios R^(-1/2) and R^(-3/2); the terminal
level eta_J of order epsilon^2 still has vanishing resource. If
tau>A_* K, choose kappa>K, then small zeta,theta>0, so that

    (1+zeta) A_* kappa/tau <1-theta.

For all sufficiently small fixed epsilon the full finite resource fits.
All the remaining limit-order and extension steps are unchanged. This
proves the sharper statement and (15). The displayed decimal is only
the elementary value20phi^5, not a numerical change to the cap interval.

## 8. Actual separated near-extremal phases

The factorization conclusion has a concrete metric consequence, not a
rank assertion. Fix an energy width eta and write e=e(eta)>0. For ANY
sequence of Hamming-radius rn covers F_n of E_n(eta), (14) and the
trivial gamma(F_n)^2<=|F_n| give

    limsup_n |F_n|/n >=(sqrt(e)-2sqrt(r))_+^2.       (18)

This follows by taking a subsequence on which the superlevel gamma2
square divided by n tends to e; no lower bound at every order is
asserted. In particular, with r=e/16, every such cover obeys

    limsup_n |F_n|/n >=e/4.                         (19)

Choose F_n to be a maximal subset of E_n(eta) with pairwise distances
STRICTLY greater than rn. Maximality makes it an rn cover, so along an
n-subsequence it has at least(e/8)n words. Half of them have the same
energy polarity, once eta is small compared with c_*; retaining that
half gives at least(e/16)n same-polarity near-extremal words, still
separated by more than(e/16)n coordinates. The equality constants in
the limsup statement (19) were deliberately reduced for this finite
subsequence formulation.

Combining this with (15), there are arbitrarily small fixed eta for
which e is at least a positive multiple of
sqrt(eta/log(e/eta)). Hence liminf-realizing signings have, along the
corresponding order subsequences, a linear-size family

    Omega(n sqrt(eta/log(e/eta)))

of same-polarity nearly extremal phases, with pairwise Hamming distance

    Omega(n sqrt(eta/log(e/eta))).

The order subsequence may depend on the chosen fixed eta. The statement
does not exchange the n and eta limits or claim every level and every
large order have this packing. Unlike single-spin flips, these phases
are separated at a fixed positive fraction of n for each selected eta,
while their normalized energy gaps are only eta.
