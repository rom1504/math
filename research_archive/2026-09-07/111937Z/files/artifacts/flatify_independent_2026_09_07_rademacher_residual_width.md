# Rademacher greedy half-width outside arbitrary small subspaces

2026-09-07. This strengthens the constant in the director's Gaussian
greedy obstruction. It is an obstruction to a separately paid residual
operator envelope, not to the actual correlated bridge operation.

## Theorem

Let A_n be hollow symmetric sign matrices with Q(A_n)=O(n^(3/2)). Let
X_n be ANY prescribed real subspaces with rank r_n=o(n). There
are Boolean x_n^+,x_n^- with

    ||P_X x_n^+||^2+||P_X x_n^-||^2=o(n),

and

    [H_A(x_n^+)-H_A(x_n^-)]/2
      >=[c_R-o(1)] n^(3/2),
    c_R=(2/3)sqrt(2/(3*pi)) = 0.3071059... .             (1)

Thus the half-width of an asymptotically near-orthogonal spin slice is
at least c_R n^(3/2). The projection tolerance can be chosen uniformly
over subspaces of the stated rank and matrices with a fixed cap bound.

## 1. Two matrix estimates paid by the actual cap

Write beta(A)=max_{u,v in signs}|u^T A v|. Cube polarization gives
beta(A)<=4Q(A). Every row of A belongs to the real cube; consequently

    sum_j |(A^2)_ij| <= beta(A).                        (2)

Since |(A^2)_ij|<=n, summing squares in (2) yields

    tr(A^4)<=n^2 beta(A)=O(n^(7/2)).                    (3)

Also ||A||op^2<=beta(A), though the Frobenius argument below only needs
(3). No unproved spectral regularity of optimizing children is used.

The subsequently proved Grothendieck factorization bound strengthens
(3) to

    tr(A^4)<=K_G^2 beta(A)^2=O(n^3).                   (3a)

See `flatify_construct_2026_09_07_spectral_fourth_moment.md`, independently
reconstructed here: factor A=D_p^(1/2) T D_q^(1/2) with probability
weights and ||T||op<=G(A), the vector Grothendieck relaxation. Factoring
the middle A^T in AA^T A gives a vector representation of AA^T A with
maximum norm product at most G(A), since |A_ij|<=1. Hence
<A,AA^T A><=G(A)^2<=K_G^2 beta(A)^2. The factorization follows from
finite SDP duality as spelled out in that artifact. Using (3a) is what
upgrades the original rank-o(sqrt(n)) theorem to ALL sublinear ranks.

## 2. Actual Rademacher greedy spins

Partition the vertices deterministically into P,Q, with |P|/n->2/3,
q=|Q| and q/n->1/3. Let B=A[P,Q]. Draw independent uniform signs
eta_j for j in Q. For i in P put

    s_i=sign(sum_j B_ij eta_j),

using independent fair tie signs when the sum is zero. Define

    x^+=(s,eta),           x^-=(-s,eta),
    Z=sum_i |sum_j B_ij eta_j|.

Every row has exactly q sign coefficients, so, for a sum S_q of q iid
uniform signs,

    E Z=|P| E|S_q|=[c_R+o(1)]n^(3/2).                 (4)

Here E|S_q|/sqrt(q)->sqrt(2/pi), by the scalar central limit theorem and
uniform integrability from E(S_q/sqrt(q))^2=1. Exactly

    H_A(x^+)-H_A(x^-)=2Z,       0<=Z<=Q(A).            (5)

The internal P and Q terms cancel in the difference. This yields
half-width, with no assumption that the two energies have opposite signs.

## 3. Covariance estimate for the greedy signs

Let K=E[x^+(x^+)^T]. All coordinates have zero mean. We claim

    ||K||F=O(sqrt(n)).                                 (6)

For two distinct rows i,l of B let

    rho=(BB^T)_il/q.

Their greedy-sign correlation satisfies the uniform estimate

    |E s_i s_l|<=C(|rho|+q^(-1/2)).                    (7)

For |rho|>1/2 this follows trivially after enlarging C. Otherwise gauge
the first row's coefficients to all +1. The second row has a agreements
and b disagreements, a+b=q, a,b>=q/4. Their two fields are U+V,U-V,
where U,V are INDEPENDENT Rademacher sums of lengths a,b.

The one-dimensional Berry--Esseen inequality gives Kolmogorov errors
O(a^(-1/2)) and O(b^(-1/2)) for the standardized U,V. Applying it to
their absolute values and then integrating the two independent marginal
CDFs shows that P(|U|>|V|) differs by O(q^(-1/2)) from the corresponding
Gaussian probability. Boundary ties have probability O(q^(-1/2)), since
U+V and U-V each have the S_q law up to signs. Therefore

    E s_i s_l=(2/pi)arcsin(rho)+O(q^(-1/2))

uniformly when |rho|<=1/2. The Gaussian identity follows from planar
rotation symmetry; arcsin is Lipschitz on this interval. This proves (7).

The sole quantitative CLT input here is scalar, iid, bounded, mean-zero,
unit-variance Berry--Esseen. For a primary reference one may take d=1
in V. Bentkus, *On the dependence of the Berry--Esseen bound on dimension*,
J. Statist. Plann. Inference 113 (2003), 385--402,
https://doi.org/10.1016/S0378-3758(02)00094-0 . Both nonempty sums above
have length at least q/4, so no degenerating-variance version is needed.

Summing the squared bound (7), including the diagonal separately, gives

    ||E ss^T||F^2
      <=C[q^(-2)||BB^T||F^2+|P|^2/q+|P|].             (8)

Schatten norm contraction under row/column projections implies
||BB^T||F^2=||B||S4^4<=||A||S4^4=tr(A^4). By (3a), the right side of
(8) is O(n). The earlier elementary (3) alone gave O(n^(3/2)), sufficient
for the former rank-o(sqrt(n)) version.

The Q,Q block of K is the identity. For the mixed block, permutation
symmetry of a Rademacher sum gives the EXACT formula

    E[eta_j s_i]=B_ij E|S_q|/q.

Indeed summing q identically distributed gauged influences recovers
E[S_q sign(S_q)]=E|S_q|; random ties contribute zero. Its squared
Frobenius norm is at most |P|, by E|S_q|<=sqrt(q). Together with (8)
this proves (6). The covariance of x^- is a diagonal-sign conjugate of
K, so it has the same Frobenius norm.

## 4. Removing large projections without losing the energy gain

For any orthogonal projection P_X of rank r,

    E||P_X x^+||^2 <=||P_X||F ||K||F <=C sqrt(rn),

and the same holds for x^-. Since r=o(n), both expectations are
o(n), uniformly as claimed. Choose a deterministic epsilon_n->0 slowly
enough that both projection squares are at most epsilon_n n except on
an event of probability o(1). Markov's inequality suffices.

By (5), deleting that event removes at most Q(A)o(1)=o(n^(3/2)) from
E Z. Thus some remaining outcome satisfies Z>=E Z-o(n^(3/2)), while
both projection constraints hold. Equations (4)--(5) prove (1).

## 5. Consequence for the proposed scalar residual envelope

For two low-cap actual children A,D and arbitrary rank-o(n)
subspaces X,Y, define the separately paid envelope

    E(A,D;X,Y)=max_{x,y} [ |H_A(x)+H_D(y)|
      +sqrt(n)||P_Xperp x||||P_Yperp y|| ].

Apply (1) separately to the children. Between the plus/plus and
minus/minus choices, at least one internal absolute sum is at least
the sum of their half-widths. Each chosen spin has residual norm
sqrt(n)(1-o(1)). Consequently

    E(A,D;X,Y)>=[1+2c_R-o(1)]n^(3/2)
       =[1.6142118...-o(1)]n^(3/2).                   (9)

This is strictly above the required balanced parent budget 2sqrt(2)M_n,
even if only the elementary limsup M_n/n^(3/2)<=1/2 is used.

The envelope is an UPPER certificate for an actual selected bridge;
a lower bound on this certificate does not lower-bound the actual
parent cap. The new low-rank localized full-sign operation remains
valid. Any successful use must exploit actual bridge/child correlations
on the residual states rather than paying its full operator norm on
every such state.
