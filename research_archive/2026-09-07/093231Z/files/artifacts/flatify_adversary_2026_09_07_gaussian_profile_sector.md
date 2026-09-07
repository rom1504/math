# A uniform actual Gaussian-profile sector bound for rank-two cross signs

2026-09-07. Proved sector theorem, not a bound on all spins and not a
selected-child transfer theorem.

Write c_G=sqrt(15)/8=.4841229182759271... . Let k=2m be a Hadamard order,
N=mk, and use any fixed Hadamard F of order k with independent signed
column permutations in its m fibres. Delete all within-fibre tiles for
the moment. In microcoordinates the cross operator K is symmetric with
norm one: it swaps each reciprocal pair and acts by H_2/sqrt(2), and is
zero on the two loop coordinates per fibre. Thus

    Z_x=x^T C_cross x/sqrt(N)=u_x^T K u_x,
    ||u_x||^2=N.

Let gamma_+ be the law of |G| for standard normal G. A physical spin x has
row magnitude laws rho_i=empirical(|F^T x_i|/sqrt(k)). Define

    delta(x)^2=(1/m) sum_i W_2(rho_i,gamma_+)^2.

For every sequence delta_k->0, there are actual cross-sign outputs for
which ALL spins with delta(x)<=delta_k have absolute energy at most

    [c_G+o(1)] N^(3/2).

More generally the bound is c_G+delta_k+o(1), uniformly for bounded delta_k.
This uses actual physical spin rows and the count at most 2^N, not an
assumed source-entropy replacement. Filling fibre blocks by low-cap signs
adds only O(N^(5/4)). It does not bound the remaining nongeneric spins.

## 1. A finite coupling with an explicit uniform error

Choose epsilon>0 such that k iid half-normal samples have empirical law
within W_2 distance epsilon of gamma_+ with probability p>=1/2. Such
epsilon=epsilon_k can tend to zero. For a fixed physical x, independently
in each row sample a Gaussian vector g_i conditioned on this empirical
event. Match its sorted absolute coordinates to the sorted word rho_i;
give the latter the Gaussian coordinate signs and ordering. The resulting
u_i is exactly a uniformly signed permutation of the physical word.
Consequently it has exactly the required construction law for this x.
Couplings need not agree for different x: only the fixed-spin probability
bound is used before the union bound.

One-dimensional sorted matching and Minkowski give

    ||u-g|| <= (delta+epsilon) sqrt(N),
    ||u||=sqrt(N),       ||g||<=(1+epsilon)sqrt(N).

Therefore, uniformly over all physical profiles with delta(x)<=delta,

    |u^T K u-g^T K g| <= e N,
    e=(delta+epsilon)(2+epsilon).                         (1)

The conditioning probability over all rows is at least 2^(-m), so dropping
the conditioning event in a nonnegative Gaussian expectation costs at most
2^m, subexponential in N. Large coordinates are paid in W_2 and in (1);
small mass alone is not a justification for discarding their energy.

## 2. Exact Gaussian determinant and a finite probability estimate

The eigenvalues of K are +1 and -1, each with multiplicity m(m-1), and
zero with multiplicity 2m. For iid standard Gaussian g and a=1+2lambda>2t,

    E exp(t g^T K g-lambda||g||^2)
      = (a^2-4t^2)^[-m(m-1)/2] a^(-m).                  (2)

Take

    a=16, lambda=15/2, t=2sqrt(15), c_G=sqrt(15)/8.

Then a^2-4t^2=16 and 2t c_G=lambda. The logarithm of (2) is exactly
-(N+2m)log2. If |u^T K u|>=2qN, the coupling implies a Gaussian quadratic
value of magnitude at least (2q-e)N, while ||g||^2<=(1+epsilon)^2N.
Markov, (2), both polarities, and the at-most-2^N physical spins therefore
bound the probability of ANY bad spin in the sector by

    exp[(1-m)log2
         -2t(q-c_G)N +t e N +lambda(2epsilon+epsilon^2)N]. (3)

In particular, at

    q=c_G(1+epsilon)^2 +(delta+epsilon)(2+epsilon)/2,       (4)

this probability is at most 2^(1-m)<1 for m>=2. The unused loop coordinates
are retained in (2); their determinant contribution even pays the row
conditioning loss at this parameter. Replacing q by q+eta improves (3)
by exp(-2t eta N).

Equation (4) is a finite, explicit, uniform error bound. It is an upper
bound on the restricted original absolute quadratic cap of actual cross
signs. No positive-semidefinite property of the four-variable kernel is
required, and no random real matrix is substituted for the output signs.

## 3. Why epsilon can vanish, with an elementary rate

For a half-normal variable Y, truncate at R=sqrt(2 log k). On [0,R],
W_2^2<=R W_1 and

    E W_1(empirical_k(Y clipped at R),law(Y clipped at R))
       <= R/(2sqrt(k)),

by integrating the elementary variance bound for the empirical CDF.
The two truncation errors have expected squared transport cost bounded by
E[Y^2 1_(Y>R)]=O(sqrt(log k)/k). Squaring the three-term triangle inequality
thus gives

    E W_2(empirical_k(Y),gamma_+)^2=O(log k/sqrt(k)).

Markov supplies p>=1/2 with epsilon_k=O(k^(-1/4)sqrt(log k)). Therefore
the error in (4) is delta_k+O(k^(-1/4)sqrt(log k)), uniformly for bounded
delta_k. The additional internal-sign completion is O(k^(-1/2)) after
normalizing by N^(3/2). No optimized Gaussian empirical-process rate is
being invoked.

## 4. Direct Gaussian typed dual and the rate function

For rho=gamma_+, take the even potential phi(z)=-lambda z^2. The exact
four-coordinate edge expectation is (a^2-4t^2)^(-1). Its two unused loop
coordinates contribute a^(-1), while the row potential correction is
k lambda. This reproduces (2) with the norm correction, directly inside
the finite typed inequality. At the asymptotic optimum a satisfies

    a^2-a-4t^2=0.

The used coordinates of the tilted joint Gaussian have variance one.
The two unused loop coordinates instead have variance 1/a, so a finite
row's average variance is 1-(2/k)(1-1/a), not exactly one. The difference
is O(1/k) and is already included in (2).

For a requested normalized cap q<1/2, choose

    a=1/(1-4q^2), lambda=(a-1)/2, t=q a.

The limiting Gaussian rate is

    I(q)=-(1/4)log(1-4q^2).

The all-spin entropy threshold I(q)=log2 gives q=sqrt(15)/8. This is the
Gaussian-profile sector only. No claim that all Hadamard row profiles are
W_2-near Gaussian is made.

## 5. Independent audit of root's mixed Gaussian/flat extension

Root subsequently proposed assigning each physical row, before the random
permutations, either to gamma_+ or to the flat magnitude law delta_1. Let
b be the flat fraction, g=1-b, and let delta^2 be the average squared W_2
distance to these assigned reference laws. The same coupling works, using
independent Rademachers in flat rows and conditioned Gaussians elsewhere.
Its quadratic error is still at most (delta+epsilon)(2+epsilon)N.

At a=16, t=2sqrt(15), lambda=15/2, let
C=log cosh(4sqrt(30)) and D=a^2-4t^2=16. The exact edge factors, with the
negative quadratic tilt on Gaussian coordinates only, are

    GG: D^(-1),
    GF: a^(-1) exp(4t^2/a),
    FF: cosh(2sqrt(2)t).

For GF, H_2 applied to a two-bit flat sign vector has only one nonzero
coordinate, of magnitude two. Direct Gaussian integration gives the
displayed factor. Thus the leading log moment per N, including Gaussian
norm compensation, is

    F(b)=lambda g-(g^2/4)log16
         +(gb/2)(-log16+15)+(b^2/4)C
        =F(0)+(b^2/4)(log16-30+C).

The linear term cancels exactly. Moreover log16-30+C<-5: use
log16<3, sqrt(30)<11/2 and logcosh(v)<v for v>0.

If s=gm is the number of Gaussian rows, the exact finite edge/loop
correction, including conditioning by probability at least 2^(-s) and
both polarities, is

    -(s/2)log16-((m-s)/2)C+s log2+log2 <=(1-m)log2.

Consequently root's mixed-sector conclusion passes independently:

    Q_sector/N^(3/2)
      <= c_G-5b^2/(8t)+delta+o(1),

uniformly in the assignment and in b. The Gaussian norm compensation
error is only lambda g(2epsilon+epsilon^2); replacing g by one gives a
uniform bound. No counting of bent functions is needed. This is still a
sector theorem: it does not claim every actual row is close to one of
these two reference laws.
