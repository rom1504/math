# Actual balanced augmentation and a necessary complexity of near-minimizers

2026-09-07. **Proved construction argument; independent audit requested.**
This is an actual full-sign extension preserving every old edge. It gives
a structural exclusion for asymptotically optimal signings, not convergence.
The elementary version uses signature pairings. A primary vector-balancing
theorem strengthens its range to sublinear factorization complexity.

## 1. Statement

Write H_A(x)=sum_(i<j) A_ij x_i x_j and Q(A)=max_x |H_A(x)|. Suppose
A_n is a sequence of symmetric hollow full signings and there are sets
G_n of Boolean reference patterns and fixed constants kappa>0, 1<alpha<2
such that, for EVERY Boolean x,

    Q(A_n)-|H_(A_n)(x)|
       >=kappa n^(3/2) [d_H(x,G_n)/n]^alpha.               (1)

Assume the sign matrix whose rows are G_n admits a factorization

    G_n=L_n V_n,
    max_i ||(V_n)_:i||_2<=1,
    max_(g in G_n) ||(L_n)_g:||_2<=gamma_n,
    gamma_n^2=o(n).                                     (2)

Then for every sufficiently small fixed epsilon>0 and q=floor(epsilon n)
there exists a symmetric hollow full signing A_n^+ of order n+q, whose
old n by n principal block is EXACTLY A_n, and

    limsup_n [Q(A_n^+)-Q(A_n)]/n^(3/2) <=D(epsilon),
    D(epsilon)=o(epsilon) as epsilon decreases to zero.   (3)

One valid rate, with constants depending only on kappa and alpha, is

    D(epsilon)=O((epsilon log(1/epsilon))^(alpha/(2alpha-2))
          +epsilon^(2alpha/(2alpha-1))+epsilon^(3/2)).     (4)

In particular (2) holds whenever |G_n|=o(n): take V=G/sqrt(|G|)
and L=sqrt(|G|)I. Thus sublinearly many uniformly robust valleys are
enough for the construction. Factorization permits more centers when
their coordinate patterns have low-dimensional structure.

The absolute energy in (1) includes both energy polarities. It automatically
requires G_n to contain every exact maximizer of |H_A|, including its
global reversal. There is no silent assumption about the cross-block
sign after reversing a child.

## 2. Elementary signature-pairing precursor

If G has R reference patterns, partition old vertices by their R-bit
signatures across G. Pair vertices within every signature cell. At most
ell<=2^R vertices remain unpaired. For each new vertex independently,
give the two old vertices in every pair opposite independent random
cross signs; give every unpaired vertex its own independent fair sign.
This is a genuine full +/-1 cross column v.

For any g in G each paired contribution to v dot g vanishes. If x has
Hamming distance d from g, at most d pairs contribute nonzero coefficients
to v dot x, each of magnitude two; the unpaired coefficients have magnitude
one. Thus

    E exp(lambda v dot x)
       <=exp[lambda^2(4d+ell)/2].                         (5)

The leftovers are RANDOM, not paid by a deterministic q ell bound.
Consequently this elementary construction already proves (3) whenever
2^R=o(n), and more generally when the actual number of odd signature
cells is o(n) and log R=o(n). At the centers, its cross contribution
is of scale q sqrt(ell), hence normalized scale epsilon sqrt(ell/n).

## 3. Stronger cross-column distribution from Gram--Schmidt balancing

The imported result is Theorem 1.4 of
[Bansal, Dadush, Garg and Lovett, The Gram--Schmidt Walk (2019)](https://theoryofcomputing.org/articles/v015a021/v015a021.pdf).
For vectors of norm at most one and initial coloring zero, it supplies
a random full sign coloring whose vector discrepancy has MGF bounded
by exp(40||theta||_2^2/2). This is an all-directions MGF statement, not
only coordinate tail bounds. The theorem and its exact convention were
checked in the primary publication.

Apply it to the vectors

    w_i=((V_n)_:i,e_i)/sqrt(2),                          (6)

whose norms are at most one. Let v be the resulting full sign vector.
For x and a nearest g in G, use the test vector

    theta=(sqrt(2)(L_n)_g:, sqrt(2)(x-g)).

Its inner product with sum_i v_i w_i is exactly v dot x, while
||theta||_2^2<=2gamma_n^2+8d_H(x,G). Therefore

    E exp(lambda v dot x)
       <=exp[lambda^2(80gamma_n^2+320d_H(x,G))/2].         (7)

The identity coordinates in (6) are essential: they preserve the small
Hamming-increment control while the other coordinates balance the reference
family. Independent samples of v form the q actual cross columns.
No real-valued cross matrix, scalar covariance rounding, or matrix-sign
replacement is being used.

## 4. The center count is automatically subexponential

For completeness, (2) implies log|G_n|=o(n); this need not be an extra
assumption. The sign family G_n has VC dimension at most gamma_n^2.
Indeed if d coordinate positions are shattered, for every sign vector
epsilon on them there is a row vector l_epsilon of norm at most gamma_n
with <l_epsilon,V_i>=epsilon_i on those positions. Hence

    d <=gamma_n ||sum_i epsilon_i V_i||_2.

Average over independent fair epsilon_i and use Cauchy--Schwarz and
||V_i||<=1 to obtain d<=gamma_n sqrt(d), so d<=gamma_n^2.

The elementary Sauer recursion bounds a binary family of VC dimension
D by sum_(j<=D) binom(n,j): splitting on the last coordinate gives the
recursion f(n,D)<=f(n-1,D)+f(n-1,D-1), with the usual boundary cases.
Thus for gamma_n^2=o(n),

    log|G_n|<=n h2(floor(gamma_n^2)/n)+o(n)=o(n).          (8)

When gamma_n^2 is below one, nonempty sign families already force gamma_n
at least one by Cauchy--Schwarz, so no exceptional denominator is hidden.

## 5. One simultaneous bridge bound for all old and new spins

Let C be the n by q cross matrix with independent columns from (7).
For fixed x,y, with x at distance d from G and y in {+-1}^q,
x^T C y is subgaussian with variance proxy

    q(80gamma_n^2+320d).

At distance d there are at most |G| binom(n,d) old configurations. Union
over every shell d=0,...,n, all old configurations, both tails and all
2^q new spin words gives, with probability at least 1/2, simultaneously,

    |x^T C y| <=T_d,
    T_d=sqrt(2q(80gamma_n^2+320d)
          log[4|G|(n+1)2^q binom(n,d)]).                 (9)

There is no unnecessary confidence term of order n in the logarithm:
constant positive existence probability suffices. By (2), (8), and
log binom(n,d)<=n h2(d/n), the normalized envelope converges uniformly
in r=d/n to at most

    F_epsilon(r)=sqrt(640epsilon r[h2(r)+epsilon log2]).  (10)

Choose any new q by q hollow signing B_q with Q(B_q)<=q^(3/2). Such a
signing exists by the elementary independent-edge tail and union over
all Boolean words; no existing strict upper theorem is required here.
Set A_n^+ to have blocks A_n,C,C^T,B_q. For every x,y,

    |H_(A_n^+)(x,y)|
       <=Q(A_n)-kappa n^(3/2)r^alpha+T_d+Q(B_q).          (11)

This inequality pays the genuine absolute cross contribution. In particular
it respects the child-reversal identity; no cancellation between separately
bounded channels is assumed.

## 6. Why the added cost is smaller than the dilution gain

Since h2(r)<=r log(e/r), (10) is at most a universal constant times

    sqrt(epsilon) r sqrt(log(e/r))+epsilon sqrt(r).       (12)

Split the negative term kappa r^alpha equally between the two summands.
For the first term, take a cutoff

    r_0=C_(alpha,kappa)
          (epsilon log(1/epsilon))^(1/(2alpha-2)).

For r>=r_0, the first summand is at most (kappa/2)r^alpha
once the constant is large and epsilon is small. For r<=r_0, the
function r sqrt(log(e/r)) is increasing, so its value at r_0 bounds
the whole positive excess by
O((epsilon log(1/epsilon))^(alpha/(2alpha-2))).

The ordinary one-variable maximum of C epsilon sqrt(r)-(kappa/2)r^alpha
is O(epsilon^(2alpha/(2alpha-1))). Together with Q(B_q), this proves
(4). Both powers before the child term exceed one precisely as needed:
alpha/(2alpha-2)>1 when alpha<2, and 2alpha/(2alpha-1)>1 always.
All constants are independent of n and of the center family.

First fix epsilon, then let n grow so gamma_n^2/n and log|G_n|/n vanish;
only afterward decrease epsilon. There is no assumption of a growing
center family being uniform in an unproved finite-dimensional limit.

## 7. A necessary geometry of asymptotically optimal signings

Let c_inf=liminf_n M_n/n^(3/2)>0. Suppose a subsequence of signings has
Q(A_n)/n^(3/2)->c_inf and satisfies (1)--(2) with fixed kappa,alpha.
Choose a sufficiently small fixed epsilon so

    D(epsilon)<c_inf[(1+epsilon)^(3/2)-1].

Then (3) yields actual signings along the larger orders with

    limsup Q(A_n^+)/(n+floor(epsilon n))^(3/2)
       <=[c_inf+D(epsilon)]/(1+epsilon)^(3/2)<c_inf,

a contradiction to the definition of c_inf. Thus no liminf-realizing
sequence can have a uniformly subquadratic robust-valley landscape with
sublinear factorization complexity as in (1)--(2).

In particular asymptotic minimizers cannot have only o(n) reference
valleys with a common kappa>0 and exponent alpha<2. They must violate
at least one of: small center complexity, uniform coercivity, or the
subquadratic exponent. This is a structural theorem about the ORIGINAL
full-sign problem obtained from an explicit successful augmentation.

It is not a proof that arbitrary minimizers admit the favorable geometry;
in fact it proves they cannot. The theorem explains a genuine obstruction
to simple seed extension: near-optimal signs must retain enough competing
near-extremal structure to prevent balanced cross columns from making
dilution profitable. No convergence conclusion is asserted.
