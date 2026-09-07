# A center-conditioned joint-Gaussian sector theorem for actual sign bridges

Date: 2026-09-07. Proved sector theorem. This is an actual full rectangular
sign bridge, not a Haar or Gaussian replacement. It does not control all
physical spin pairs, and does not prove original convergence.

The scalar Gaussian-profile coupling is extended from
`flatify_adversary_2026_09_07_gaussian_profile_sector.md`. Root proposed the
bipartite rank-two realization and the center-conditioned two-state target.
The antipodal transport coupling, finite Gram-tilted determinant bound,
and rational certificate are derived below.

## 1. Exact sign ensemble and two-state profiles

Let k=2m be a Hadamard order and n=mk. There are m left and m right
fibres, each of size k. In each fibre independently use a fixed k-order
Hadamard followed by a uniform signed column permutation. Assign its k
columns in pairs to the m opposite fibres. On each fibre pair put

    C_ij=(1/2) F_i[:,pair j] H_2 G_j[:,pair i]^T.

Every entry is a sign: for u,v in {+1,-1}^2, u^T H_2 v is +2 or -2.
Writing U,V for the block diagonal normalized Hadamards, exactly

    C/sqrt(n)=U K V^T,                                  (1)

where K is orthogonal, pairing reciprocal two-coordinate slots by
H_2/sqrt(2). In particular ||C||op=sqrt(n). There are no unused coordinates
and no missing diagonal blocks in this rectangular construction.

Fix arbitrary physical centers x0,y0 and candidate spins x,y, all of
length n. Set alpha=4/5, and first suppose

    <x0,x>/n=<y0,y>/n=alpha.                            (2)

In a left fibre the unrandomized two-state word consists of k vectors

    v_l=((F^T x0_i)_l,(F^T x_i)_l)/sqrt(k) in R^2.

Use its empirical law on the antipodal quotient R^2/{+1,-1}, whose metric
is d([v],[w])=min(||v-w||,||v+w||). Let gamma_G be the quotient of the
centered Gaussian law with covariance

    G=[[1,4/5],[4/5,1]].

Define delta_L^2 as the average over left fibres of the squared quotient
W_2 distance to gamma_G, and delta_R similarly. The joint-Gaussian sector
consists of pairs with delta_L,delta_R<=delta. The two replica coordinates
must be treated jointly; two separate Gaussian magnitude profiles are
not sufficient for this definition or theorem.

## 2. Quotient matching gives a finite, exact construction-law coupling

For k iid N(0,G) vectors, choose epsilon=epsilon_k tending to zero such
that their empirical quotient law is within epsilon of gamma_G with
probability p>=1/2. This follows from empirical W_2 convergence for a
fixed finite-second-moment law in R^2, also proved by truncation into a
finite grid followed by the law of large numbers. No growing-dimensional
normal approximation is used.

Independently in every fibre sample such a Gaussian word conditioned on
this event. Optimal k-to-k quotient matching and the triangle inequality
couple it to a signed permutation of the physical word with row squared
error at most k(delta_i+epsilon)^2. The physical marginal is exactly a
uniform common signed permutation: after matching, apply an independent
uniform signed permutation simultaneously to both words. The Gaussian
conditional law is invariant under this action, and the physical orbit
becomes uniform. Quotient matching chooses the common sign of both
replicas, not their signs independently.

Equivalently, quotient W_2 equals W_2 of the symmetrized lifts: a quotient
coupling can be lifted with common fair orientation and the cheaper sign,
while projecting any lifted coupling cannot increase its cost. This
also explains why exact two-dimensional empirical-type counting, with
its potentially linear entropy overhead, is unnecessary.

Write A,B for the n-by-2 randomized physical feature arrays, and a,b for
their coupled Gaussian arrays. Minkowski yields, with eta=delta+epsilon,

    ||A-a||F,||B-b||F <=eta sqrt(n),
    ||A||F=||B||F=sqrt(2n),
    ||a||F,||b||F<=(sqrt(2)+eta)sqrt(n).

Therefore, with e=eta(2sqrt(2)+eta),

    ||A^T K B-a^T K b||F <=e n,
    ||A^T A-a^T a||F,||B^T B-b^T b||F <=e n.             (3)

By (2), A^T A=B^T B=nG exactly. Dropping all the Gaussian conditioning
events in a nonnegative expectation costs at most 2^(2m)=exp(o(n)).
All large coordinates are paid by the W_2 error and (3).

## 3. General finite Gram-tilted determinant estimate

Let T=diag(t0,t1) with t1>0, and let L,R be symmetric 2-by-2 matrices
such that

    P=[[G^(-1)+L,-T],[-T,G^(-1)+R]] >0.

Define

    F(T,L,R)=(1/2)Tr((L+R)G)-(1/2)logdet(P)-logdet(G),
    D(T,L,R)=||T||F+(||L||F+||R||F)/2.

Because K is orthogonal, replacing the Gaussian array b by Kb leaves
its law and its Gram matrix unchanged. The n coordinate pairs then are
independent N(0,G) pairs. Exact Gaussian integration gives

    E exp[Tr(T a^T K b)
       -(1/2)Tr L(a^T a-nG)-(1/2)Tr R(b^T b-nG)]
       =exp[n F(T,L,R)].                                (4)

Put h0=x0^T C y0/n^(3/2), h1=x^T C y/n^(3/2). On the physical event
|h0|<=c and h1>=b, the compensated physical exponent is at least
n(t1 b-|t0|c). Equations (3)--(4), Markov, and the conditioning cost prove

    P(|h0|<=c, h1>=b)
      <=exp[2m log2-n(t1 b-|t0|c-F-eD)].                 (5)

For -h1>=b use -T, with the same L,R and determinant. Both polarities
therefore cost only a factor two. This is a fixed-physical-pair estimate
uniform over the stated joint-profile sector.

## 4. A fully rational useful point

Take

    T=diag(-8/5,10/3),
    L=R=[[64/45,-16/9],[-16/9,35/9]],
    b0=37/50.

The inverse precision has the exact form

    P^(-1)=[[G,J],[J,G]],
    J=[[0,2/5],[2/5,37/50]],
    det(P^(-1))/det(G)^2=25/108.

It is positive definite, and substitution in (4) gives

    t1 b0-F=I0:=(1/2)log(108/25)
                    =.731627701128009... .             (6)

The matrix identities and positive leading minors are checked exactly
in `computations/flatify_independent_2026_09_07_joint_gaussian_bridge_certificate.py`.
Its directed interval calculation certifies

    I0-2h(1/10)>.08.                                   (7)

For intuition only, in orthonormal replica coordinates the covariance
cross-block here is [[0,2/3],[2/3,5/18]]. The optimized finite-rank Haar
rate is the same determinant expression. The Haar block density appears
in Jiang, *Annals of Probability* 34 (2006), Lemma 2.5:
[primary paper](https://arxiv.org/pdf/math/0601457).
Our actual-sign proof above uses Gaussian integration and transport,
not a substitution of Haar matrices for signs.

## 5. Uniform control of an actual physical cloud sector

Consider the exact Hamming spheres at distance n/10 from x0 and y0,
along compatible n divisible by 10. They have at most exp[2n h(1/10)]
pairs. Combining (5)--(7), the probability that a pair in their joint-
Gaussian sector violates |h1|<=b while |h0|<=c is at most

    2 exp[2m log2+n{2h(1/10)-I0
                    -(10/3)(b-b0)+(8/5)c+eD}].         (8)

For delta_k->0, epsilon_k->0 and c_k=n^(-1/4), take b=b0. The right
side tends to zero exponentially. Independently, the center has

    E h0^2=1/n,

since each randomized physical feature vector has covariance I_n,
left and right are independent, and K is orthogonal. Thus
P(|h0|>n^(-1/4))<=n^(-1/2). For all sufficiently large compatible orders
there consequently exist actual sign bridges satisfying both

    |x0^T C y0|<=n^(5/4),
    |x^T C y|<=(37/50)n^(3/2)

for EVERY physical pair in the two Hamming spheres whose joint-Gaussian
profile error is at most delta_k. This is a simultaneous sector bound,
not merely an expected-energy assertion.

## 6. Exact original-value relevance and remaining obligations

For actual children A,D, restrict further to the typical-energy windows
around their chosen same-polarity noisy extrema. If both child caps are
at most (c+o(1))n^(3/2), their combined absolute energy in these windows
is at most (32c/25+o(1))n^(3/2). The actual parent identity therefore gives
on the above joint-Gaussian cloud sector

    |H_A(x)+H_D(y)|+|x^T C y|
       <=[37/50+32c/25+o(1)] n^(3/2).                  (9)

At c=.493608094 this is more than .02 below 2sqrt(2)c; the directed
certificate checks that margin. Thus the center-conditioned actual-sign
estimate genuinely escapes the unconditional uniform-one-state union
obstruction identified in `flatify_independent_2026_09_07_correlated_bridge_shell_barrier.md`.
More generally (9) is below 2sqrt(2)c only if
c>(37/50)/(2sqrt(2)-32/25); no smaller-c seed closure is claimed.

Two substantial gaps remain before a whole-parent theorem: physical
pairs outside the joint-Gaussian sector are uncontrolled, and the other
energy shells and their centers have not been covered at a paid entropy
cost. In particular Gaussianity of a typical row does not license a union
bound over all exceptional rows. The theorem is along the stated
compatible orders; no all-order child-preserving interpolation is asserted.
