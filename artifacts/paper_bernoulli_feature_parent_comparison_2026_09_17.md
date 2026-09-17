# Feature laws compare complete child maxima, on a declared old code

2026-09-17. Final combination of the all-offset scalar feature theorem
with the director's codewise mean-branch comparison. Independently
reconstructed by the director. This is a comparison of actual sign
bridges to a Gaussian VARIANCE MIXTURE, not a numerical parent bound.

## 1. A finite all-offset comparison theorem

Let nu be a centered physical law on{+-1}^n, and let G be a centered
law on R^n. Assume both have linear subGaussian proxy K and that

    sup_(x Boolean,s real)
      |E_nu|s+h.x/sqrt(n)|-E|s+G.x/sqrt(n)||<=epsilon. (1)

Take q independent physical columns h_j and q independent columns G_j
of the respective laws. Fix an old signing A, a child signing D of
order q, and a nonempty declared code C subset{+-1}^n, all chosen
BEFORE these columns. Set

    Q_C(h)=max_(x in C,y Boolean)
          |H_A(x)+H_D(y)+sum_j(h_j.x)y_j|,

and define Q_C(G) by the same expression. Then

    |E Q_C(h)-E Q_C(G)|
      <=q sqrt(n) epsilon+sqrt[4Kqn log(2|C|)].      (2)

All q new spins and both absolute polarities are retained. No
independence between responses at different old words is assumed.
The old code IS restricted; control of its complement is not part of
this theorem. The child may depend on the old data but not on the newly
sampled bridge. The proof equally permits arbitrary deterministic old
offsets and arbitrary functions of the new spins at each declared index.

## 2. Why the scalar comparison retains the child exactly

Index I by (x,sigma) in C times{+-1}, write x_i=sigma x,
b_i=sigma H_A(x), and d_i(y)=sigma H_D(y). Then

    Q_C(h)=max_i F_i(h),
    F_i(h)=b_i+max_y[d_i(y)+sum_j(h_j.x_i)y_j].

For a fixed index and all but one normalized scalar field frozen,
the remaining dependence is EXACTLY

    max{sqrt(n)u+B_+,-sqrt(n)u+B_-}
      =(B_++B_-)/2+sqrt(n)|u-(B_--B_+)/(2sqrt(n))|.

Successive replacement of the q independent scalar fields using(1)
therefore gives

    |E F_i(h)-E F_i(G)|<=q sqrt(n) epsilon.        (3)

This is why a uniform affine-offset comparison is used rather than
merely comparing the unshifted absolute field mean or its covariance.

## 3. One fully paid old-code fluctuation term

For a centered K-subGaussian scalar V and an independent copy V',
every sqrt(n)-Lipschitz function phi satisfies

    E exp[t(phi(V)-E phi(V))]
      <=E cosh[t(phi(V)-phi(V'))]
      <=E cosh[t sqrt(n)(V-V')]<=exp(Kn t^2).       (4)

The first inequality is symmetrization, the second uses the Lipschitz
bound and monotonicity of cosh in absolute value, and the last follows
from the two scalar MGF bounds. This is valid without assuming V is
symmetric. Doob martingale tensorization proves that a separately
sqrt(n)-Lipschitz function of q independent such scalars has centered
MGF proxy2Kqn. Each F_i has those separate Lipschitz constants.

Consequently, for either ensemble, its expected maximum lies in

    [max_i E F_i, max_i E F_i+sqrt(4Kqn log|I|)].

Equation(3) compares the two maximum means. The two one-sided intervals
then give(2) with ONE fluctuation term, not two. At a single signed
index the corresponding log|I| term is correctly zero; one old word
still has two signed indices when both absolute polarities are included.

## 4. Application to the actual growing-rank feature realization

Fix finitely many v_j>0 and weights pi_j with sum pi_j v_j=1. Let P
be a rank-r projection with leverage max P_ii<=Lr/n. The independently
proved cold/hot construction and exact physical covariance repair give
one centered, exactly isotropic, full-support physical law nu for
r=o(sqrt(n)), with a fixed common linear proxy K. Its scalar comparator
is the genuine column law

    choose V=v_j with probability pi_j,
    G | V ~ N(0, I+(V-1)P).                       (5)

Every matrix in(5) is positive definite, with eigenvalues1 and V. The
mixture has covariance I, but its distribution is NOT replaced by
N(0,I). The variance mixture is the mechanism retaining the cheaper
absolute profile at projected queries.

The [hot proof](paper_bernoulli_hot_feature_laws_2026_09_17.md), (22a),
is uniform over all affine offsets. For the
[cold proof](paper_localization_growing_rank_cold_tilt_2026_09_17.md),
symmetry writes the difference of shifted-absolute expectations as
the same characteristic-function difference multiplied by cos(ts);
its absolute value is bounded by the already paid integral, uniformly
in s. The exact-isotropy repair also changes shifted-absolute means
by O(repair mass), uniformly s, after subtracting the common constant
|s| and using ||s+z|-|s||<=|z|. Thus(1) holds with epsilon=o(1).

If q=O(n) and log|C|=o(n), (2) gives an o(n^(3/2)) expected restricted-
parent comparison, including the fixed child's FULL energy and every
new-spin configuration. This is a leading-order sign realization of the
declared Gaussian-mixture parent problem, not its solution.

One may choose C to be a subexponential nearcode supplied by the actual
[cloned-block preparation](paper_director_cloned_block_regularization_2026_09_17.md).
That composition retains the explicit preparation cap cost. It does
NOT control far old words, establish that this nearcode has favorable
feature capture, or bound the value of the mixture-Gaussian parent.
No improved original extremal constant or convergence conclusion is
claimed by this restricted-code theorem.

The scalar transform, covariance repair, information budget and
concentration hypotheses each have separate proofs. Matching a Gaussian
profile without those physical-law inputs would not imply(2).

The localization track independently read and reconstructed the full
finite theorem and feature application: PASS, including the exact child
max identity, Doob MGF, single fluctuation interval, and retained mixture
distribution. The director independently reconstructed the same argument.
The discrepancy track then independently read the full file and returned
PASS as well, including the bridge-independent child assumption and
the single fluctuation term.

## 5. Quantitative composition with actual cloned-block preparation

Fix0<gamma<2/3 and0<=rho<1/2 before taking n to infinity. The actual
cloned-block preparation theorem, at the SAME old order n, replaces
any supplied full signing A by a full signing W with

    Q(W)<=Q(A)+O(n^(3/2-gamma/2)),
    C={x: |H_W(x)|>=Q(W)-T},
    T=n^(3/2-gamma),
    log|C|=O(n^(1-gamma/4)log n).                 (9)

All absolute polarities and ties are included in C. This is the stated
fixed-gamma preparation theorem, not a cost-free code hypothesis.

Choose feature frames of rank r=O(n^rho), leverage at most Lr/n, and
fixed variance-mixture parameters. The all-offset scalar realization
and the conservative Frobenius covariance repair bound give

    epsilon=O(n^(rho/2-1/4)+n^(rho-1/2)).         (10)

The newer trace-sensitive repair estimate is stronger, but is not
needed for the following rate. Let q=O(n), and choose any fixed actual
child D independently of the subsequent bridge. Apply(2) with old
child W and its COMPLETE declared nearcode C. Since rho<1/2, the
second summand in(10) is smaller than the first. Therefore

    |E Q_C^W(h)-E Q_C^W(G)|
       =O(n^(3/2)[n^(rho/2-1/4)
                           +n^(-gamma/8)sqrt(log n)]).  (11)

The first term is q sqrt(n)epsilon. The second is the exact old-code
fluctuation term, since sqrt(qn log|C|)
=O(n^(3/2-gamma/8)sqrt(log n)). Constants may depend on the fixed
gamma,rho,leverage,variance parameters and bound for q/n, not on n.

For the explicit choice rho=1/3 and gamma=1/2,

    old preparation cap cost: O(n^(5/4)),
    declared absolute nearcode window: T=n,
    restricted parent comparison error: O(n^(23/16)sqrt(log n)). (12)

Indeed the feature replacement contributes O(n^(17/12)), which is
smaller than the displayed fluctuation contribution n^(23/16)sqrt(log n).
The parent comparison retains ALL q new spins and both polarities.

The superscript W in(11) matters: BOTH parent ensembles have the
prepared old signing W. The preparation gives the cap inequality in(9),
not a uniform pointwise bound on H_W-H_A that could be inserted freely
inside an arbitrary parent optimization. Nor does(11) assert that C
has favorable feature capture, bound the mixture-Gaussian parent VALUE,
or exclude old words outside C. Consequently these quantitative rates
are not an original extremal upper bound or a convergence theorem.

The localization track independently audited all of Section5: PASS.
Its exact rational exponent replay returned5/4,1,17/12,23/16 for
preparation,window,scalar replacement and old-code selection respectively.
At fixed-density q proportional to n, the displayed selection-error
certificate divided by the protected window is of order
n^(7gamma/8)sqrt(log n). This comparison bound therefore does not itself
pay escape from the declared nearcode. It is a limitation of the bound,
not a lower bound on the actual transfer error.
