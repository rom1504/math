# Every fixed-size iid frame has an actual super-half floor

2026-09-07. **Proof complete; independent audit requested.** This extends
the independently audited two-row theorem to every fixed k>=2 in the
iid directed-feature law. It does not assert a general-k orthogonal-frame
repair; the separate two-row exact-balance coupling already covers k=2.

## Actual construction and statement

Fix k>=2 and a prescribed hollow symmetric signing S of order n. At each
directed macro port i->j, choose a vector

    a_ij=(1,eta_ij,1,...,eta_ij,k-1),

whose last k-1 entries are independent fair signs, independently over all
directed ports. For i<j put the actual cross-fibre block

    W_ij=S_ij a_ij a_ji^T,     W_ji=W_ij^T.

Fill each physical diagonal fibre by any fixed hollow signing. The full
parent has order N=kn, and its first coordinate in each fibre induces
EXACTLY the prescribed child S. Every physical off-diagonal entry is a
sign. Uniformly over every fixed S, without even a cap assumption,

    Q(W)/(kn)^(3/2) >= c_k-o_P(1),

    c_k = [(2/pi)(k-1)^(3/2)
                  +sqrt(2/pi)sqrt(k-1)]/k^(3/2).

The constants c_k increase strictly in k>=2, so their minimum is

    c_2=0.5071738708131547... >1/2.

No event uniform over an adaptively chosen S is claimed. Correlated
directed frames, rare choices, and globally rewritten parent operations
are outside this theorem.

## 1. Independent BLOCK Gaussian replacement

Regroup physical coordinates so the first coordinates form the prescribed
child. At each unordered macro edge the k^2 physical entries constitute
one independent random block. Its (1,1) entry is the deterministic S_ij;
all other entries have mean zero. For indices r,s,r',s' in {0,...,k-1},

    E[a_ij,r a_ij,r']=1_(r=r'),

and likewise at the other endpoint. Thus the centered block vector has
identity covariance on its k^2-1 random entries. This assertion is only
about covariance: the block entries are NOT independent.

Replace each centered block by independent standard Gaussian entries,
retaining the deterministic S corner. Define the soft absolute maximum
in n-normalized energy units by

    F_beta(V)=(1/beta) log sum_(x in {+/-1}^(kn),sigma=+/-1)
                                      exp[beta sigma H_V(x)/sqrt(n)].

It differs from Q(V)/sqrt(n) by between zero and
`(kn+1)log(2)/beta`. Along any block direction z, the third derivative
has magnitude at most

    C beta^2 n^(-3/2) ||z||_1^3.

Indeed it is a third Gibbs cumulant of a linear combination of at most
k^2 spin products; each spin product has magnitude one. This bound is
independent of the values of all deterministic offsets, including S.
Taylor replacement of a block cancels its constant, linear, and quadratic
terms because the means and covariance agree. For fixed k both its sign
and Gaussian third absolute l1 moments are bounded by a constant C_k.
Summing over fewer than n^2/2 blocks gives

    |E F_beta(W)-E F_beta(W_G)| <= C_k beta^2 sqrt(n).

The fixed physical diagonal fibres contribute only O_k(n) in unnormalized
energy and can be omitted throughout this comparison. Taking beta=n^(1/6)
and restoring ground-energy units yields

    |E Q(W)-E Q(W_G)| <= O_k(n^(4/3)).                 (1)

The O_k(n) filling correction is included. This is a quantitative actual
cap universality statement for this independent-block law, not merely
an entrywise covariance heuristic.

## 2. Removing the prescribed deterministic child for a LOWER bound

Write W_G=A_S+G, where A_S is S supported on its first-coordinate child
and G is the centered Gaussian matrix on the other edge positions. The
function

    Phi(A)=E Q(A+G)

is convex. Since G and -G have the same law and Q(-V)=Q(V), Phi is even.
Therefore `Phi(A_S)>=Phi(0)`. This step only LOWER-bounds the expected
parent cap; it does not claim that deleting a child preserves its cap.

Let d=(k-1)n. Up to O_k(n) expected absolute edge cost, the Gaussian
matrix G is the model

    [[0,C],[C^T,D]],

with C a full n by d iid standard Gaussian bridge and D an independent
hollow iid Gaussian symmetric matrix. The discrepancy consists only of
the originally absent within-fibre edges, O_k(n) many. Adding them is
legitimate for this expectation lower bound since the cap changes by at
most their total absolute coefficients.

Choose y by any positive-energy witness for D, then choose
`x_i=sign((Cy)_i)`. Conditional on D and y, the n bridge fields are
independent N(0,d). Hence

    E[x^T C y]=n sqrt(2d/pi).

There is no first-child energy after the convexity reduction. It follows
that

    E Q(G) >= E max_y H_D(y)+n sqrt(2d/pi)-O_k(n).       (2)

## 3. A self-contained Gaussian child witness

For a d by d GOE matrix Z with off-diagonal variance one and diagonal
variance two, let lambda_1 be its largest eigenvalue and v a unit top
eigenvector. Orthogonal invariance of the Gaussian density makes v uniform
on the sphere and independent of the eigenvalues; conditional on v and
the eigenvalues the other eigendirections are uniform in v-perp. Put
`y=sign(v)`. Sphere integration gives

    E||v||_1^2=1+2(d-1)/pi.

Conditional averaging of the remaining eigendirections therefore gives

    E[y^T Z y]
       = E lambda_1 * E||v||_1^2
         + E[(Tr Z-lambda_1)/(d-1)]
                                 * [d-E||v||_1^2]
       = (2d/pi) E lambda_1.

The top eigenvalue is simple almost surely; its global sign does not
affect the quadratic witness. Removing the Gaussian diagonal subtracts
Tr Z/2 from every Boolean energy, with expectation zero. Consequently

    E max_y H_D(y) >= (d/pi) E lambda_1.

For completeness the needed lower bound on E lambda_1 is only
`liminf E lambda_1/sqrt(d)>=2`, not an edge fluctuation theorem. Fixed
closed-walk moments give convergence in probability of the empirical
measure of Z/sqrt(d) to the semicircle distribution: even tree walks
give the Catalan moments, odd leading moments vanish, and connected-pair
counts give vanishing variance. The second-moment bound ensures tightness,
and the compactly supported moment-determinate semicircle law identifies
every subsequential limit. Thus for every fixed epsilon>0 a positive
fraction of eigenvalues lies above 2-epsilon with probability tending to
one. Also lambda_1>=Tr Z/d, so its negative part has expectation at most
E|Tr Z|/d=O(d^(-1/2)). This proves the required expectation lower bound.
Hence

    E max_y H_D(y) >= (2/pi-o(1))d^(3/2).              (3)

This proof uses only orthogonal Gaussian invariance, sphere integration,
and fixed-walk moments. The auxiliary Gaussian diagonal is part of the
witness randomization, not an unallowed physical diagonal in D.

## 4. Constant and concentration

Equations (1)--(3) give

    liminf E Q(W)/n^(3/2)
       >= (2/pi)(k-1)^(3/2)+sqrt(2/pi)sqrt(k-1).

For a=2/pi and b=sqrt(2/pi), differentiate
`(ar+b)sqrt(r)/(r+1)^(3/2)`, r=k-1. The derivative has the sign of
`r(3a-2b)+b`, which is positive because 3a>2b. Therefore c_k increases.

Finally each independent unordered macro block changes at most k^2
physical edges, with cap sensitivity at most 2k^2. Bounded differences
over O(n^2) independent blocks concentrates Q on scale O_k(n), with
probability at most exp(-c_(k,delta)n) of a fixed delta n^(3/2)
deviation below its mean. Fixed within-fibre fillings alter cap by only
O_k(n). This proves the claimed actual high-probability floor.

The k=2 balanced-frame theorem is separately audited. Extending the
coupling to mutually orthogonal k-row physical frames for every k would
require a new finite-type repair argument and is not included here.
