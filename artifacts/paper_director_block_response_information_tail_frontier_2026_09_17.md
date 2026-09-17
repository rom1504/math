# Information cost of approaching the sharp compositional tail floor

2026-09-17. **Proved; independent discrepancy and localization reconstructions PASS.** This combines
the sharp large-rank block-code response floor, finite binomial information
bounds, and exact physical hot/cold realization. It is an independent
benchmark theorem; no conclusion about original minimizing signings follows.
External novelty is unestablished.

## 1. Precise two-resource asymptotic quantity

Fix K>1 and write kappa=sqrt(2/pi). For n=br let C_(b,r) be all Boolean
words constant on each of r blocks of size b. Let J_K(delta) be the infimum
of

    liminf_j D(nu_j || Uniform_{b_j r_j})/r_j

over sequences b_j->infinity, r_j->infinity, of centered, EXACTLY isotropic,
full-support physical sign laws with linear subGaussian proxy at most K,
and

    limsup_j max_(x in C_(b_j,r_j))
       E_(nu_j)|h.x|/sqrt(b_j r_j) <= kappa/sqrt(K)+delta.

The infimum over an empty collection is infinity. Then

    lim_(delta down to0) J_K(delta)/log(1/delta)=1-1/K.     (1)

Thus at a FIXED tail budget K, approaching its best possible large-rank
response forces divergent information PER FEATURE. This is distinct from
the simpler r log(1/epsilon) response-price theorem, in which the allowed
tail budget can grow with the desired response. Neither information cost
nor covariance alone captures this joint constraint.

All limits are explicit: K is fixed first; for each fixed positive delta
the physical dimensions tend to infinity; delta tends to zero last.

## 2. Lower bound: near equality forces two separated norm regions

Set S_a=b^(-1/2) sum_(i in block a) h_i and Y=||S||/sqrt(r).
Exact isotropy gives E Y^2=1. The assumed physical K-proxy passes to S.
The sharp block-code proof gives, uniformly in b and the law,

    E Y <= 1/sqrt(K)+delta/kappa+o_r(1),                  (2)

and exponentially small upper tails of Y above sqrt(K(1+eta)) for every
fixed eta>0, including their first and second tail moments.
For completeness, the first assertion averages all query signs z and
uses

    E_z|z.S| >= kappa||S||-3||S||infinity,
    E||S||infinity <=sqrt(2K log(2r)).

The second follows from Gaussian integration of the linear MGF:
Eexp(t||S||^2)<=(1-2Kt)^(-r/2).

The function g(y)=y-y^2/sqrt(K) is nonnegative on [0,sqrt(K)].
Equations (2) and E Y^2=1 imply
E g(Y)<=delta/kappa+o_r(1). Its negative part has expectation tending
to zero: first truncate at sqrt(K)+eta, use the vanishing far tail, and
then send eta to zero. Consequently, for each fixed small a>0,

    limsup_j P(a<Y<sqrt(K)-a) <=C_K delta/a.

Also Markov's inequality gives

    P(Y>=sqrt(K)-a)<=1/(sqrt(K)-a)^2=1/K+O_K(a).

Thus the event E_a={Y<=a} satisfies

    liminf_j nu_j(E_a)>=1-1/K-C_K(a+delta/a).              (3)

This is a quantitative near-equality consequence, not an assumed
two-phase description of all admissible laws.

## 3. Small-ball reference cost and entropy

Under the independent physical cube, the S_a are independent normalized
binomial sums. The finite estimate

    E exp(-t|S_a|)<=C(1/t+b^(-1/2))

is proved in the block-code information theorem by a lattice geometric
series and the central binomial maximal-atom bound. On E_a,
sum|S_a|<=ar. Markov's inequality with t=1/a yields, when b>=a^(-2),

    Uniform(E_a) <= (C1 a)^r.                             (4)

Binary data processing of relative entropy gives for any event E

    D(nu||Uniform)>=nu(E)log(1/Uniform(E))-log2.

Take a small enough that C1a<1. Combining (3)--(4) gives

    J_K(delta)
      >=[1-1/K-C_K(a+delta/a)] log(1/(C1a)).              (5)

The inequality is uniform over all admissible sequences, so the infimum
in the definition introduces no additional limit exchange. Set
a=delta^alpha for any fixed 0<alpha<1, take delta down to zero, and
then alpha up to one. This proves the lower limit in (1).
More quantitatively, choosing a=delta log(1/delta) in (5) gives

    J_K(delta)>=(1-1/K)[log(1/delta)-loglog(1/delta)]-O_K(1).

An O(1) additive remainder relative to the leading logarithm alone is
not proved. The coefficient in (1) does not rely on such a remainder.

## 4. Upper bound: exactly isotropic physical hot/cold blocks

Fix a small v in (0,1). In ONE block of size b use the actual laws

    dP_(v,b)/dUniform proportional
       exp[(1-1/v)(sum h_i)^2/(2b)],
    dP_(K,b)/dUniform proportional
       exp[(1-1/K)(sum h_i)^2/(2b)].

They are centered, permutation invariant and have full support. The
reconstructed rank-one feature theorems show, as b->infinity,

    v_b=E_(P_v) S^2 ->v,    k_b=E_(P_K) S^2 ->K,
    D(P_(u,b)||Uniform) ->(u-1-log u)/2, u in{v,K}.

The hot law is EXACTLY K-subGaussian at every b. The cold law is
(1+o_b(1))-subGaussian; hence it also has proxy at most K eventually.
These are physical MGF statements, not inferences from covariance.

For sufficiently large b, v_b<1<k_b. Choose one common hot label with
probability p_b=(1-v_b)/(k_b-v_b). Conditional on this label, sample
the r blocks independently from the appropriate block law. Within each
block the averaged second magnetization moment is exactly one, so all
off-diagonal physical covariances vanish. Different blocks are centered
and conditionally independent. The joint physical law is therefore
centered and EXACTLY isotropic. It has full support and proxy at most K.

Every block-constant query has the same distribution, by within-block
sign symmetry. Conditional scalar CLT gives, as b,r->infinity,

    response ->kappa[(K-1)sqrt(v)+(1-v)sqrt(K)]/(K-v)
       =kappa/sqrt(K)+kappa(1-1/K)sqrt(v)+O_K(v).          (6)

Uniform moment bounds from the common K-proxy justify this CLT along
joint sequences; alternatively a diagonal sequence suffices for the
infimum in J_K. No prescribed relation between b and r is needed.

Entropy convexity and tensorization within each phase bound the
information per block by

    [(K-1)/(K-v)](v-1-log v)/2
       +[(1-v)/(K-v)](K-1-log K)/2 + o_b(1)
      =(1-1/K)log(1/sqrt(v))+O_K(1).                     (7)

Choose sqrt(v)=delta/[2kappa(1-1/K)] for sufficiently small delta.
The strict factor-two margin in (6) permits the finite limits. Equation
(7) proves J_K(delta)<=(1-1/K)log(1/delta)+O_K(1), completing (1).

## 5. Interpretation and scope

The limiting response floor is achieved by allocating variance to a
shared hot phase, but exact attainment forces the remaining positive
probability mass of the global normalized feature norm toward zero. The binomial small-ball
cost charges precisely this fraction, 1-1/K, in the leading information
divergence. The proof therefore combines a tail constraint, an extremal
query requirement, physical realizability, and composition.

It does not identify low-dimensional order parameters for original
minimizers, give their optimal bridge law, or prove a cross-order defect.
The value of this benchmark is its sharp quantitative interaction between
two resources that the separate covariance and entropy statements do
not capture.
