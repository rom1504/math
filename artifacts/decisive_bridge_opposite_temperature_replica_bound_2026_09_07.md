# Exact opposite-temperature replica bound at a width optimizer

Let A minimize L_A(lambda)=(log Z_A(lambda)+log Z_A(-lambda))/2 over
actual hollow signings of order N, with lambda>0 the raw temperature.
Let X,Y be independent Gibbs spins at temperatures +lambda,-lambda,
and write C^+_ij=E XiXj, C^-_ij=E YiYj. Both C matrices are positive
semidefinite correlation matrices. This note does not identify width
optimization with the original absolute-pressure optimization.

## 1. Exact mask representation

Set u=(X+Y)/2 and v=(X-Y)/2. Their supports partition the vertices into
the match mask S and its complement, with signs on the occupied sites.
Directly,

    H_A(X)-H_A(Y)=2 sum_(i in S,j notin S) A_ij u_i v_j.

Hence Z_A(lambda)Z_A(-lambda) is exactly the sum, over all masks S,
of the bipartite partition on A_(S,S^c) at raw temperature 2lambda.
In this product Gibbs law, the overlap is

    (1/N) X dot Y=(2|S|-N)/N.

Its squared expectation equals Tr(C^+ C^-)/N^2, including diagonals.

## 2. Edge optimality retains the inactive-mask probability

For an edge e=(i,j) put W_e=XiXj-YiYj, taking values 0,+2,-2. Then

    p_e=P(W_e!=0)=(1-C^+_e C^-_e)/2.

The ratio of the paired partitions after flipping A_e is exactly

    E exp(-2lambda A_e W_e)
      =1+p_e(cosh(4lambda)-1)
                  -(A_e E W_e/2)sinh(4lambda).

It is at least one by actual width-sign optimality. Therefore

    A_e(C^+_e-C^-_e)
        <=tanh(2lambda)(1-C^+_e C^-_e).              (1)

Summing and keeping the diagonal contribution gives

    Tr(C^+ C^-)/N^2
       <=1-4 L'_A(lambda)/(N^2 tanh(2lambda)).        (2)

This is stronger than discarding the probability that the two replicas
agree across an edge. It uses one common actual signing for both
temperatures, not separately optimized branch signings.

## 3. Normalization and integrated consequence

Let psi_N(beta)=min_A L_A(beta/sqrt(N))/N. At almost every beta an
active minimizing branch has derivative psi'_N(beta), so (2) becomes

    Qopp_N(beta):=Tr(C^+ C^-)/N^2
       <=1-4 psi'_N(beta)/(sqrt(N) tanh(2beta/sqrt(N))).

Equivalently,

    psi'_N(beta)
      <=[sqrt(N)/4]tanh(2beta/sqrt(N))(1-Qopp_N(beta)). (3)

The optimized envelope is Lipschitz on compact beta intervals; branch
switches do not invalidate integrating this a.e. inequality. Since
psi_N(0)=log2, and tanh x<=x, it follows that

    integral_0^beta s Qopp_N(s) ds
       <= beta^2/2-2(psi_N(beta)-log2).              (4)

The optimizer and its two Gibbs matrices may vary measurably with s.
This bound is not a proof that the opposite-temperature overlap vanishes.
It measures its integrated size by the deficit from the Gaussian/annealed
quadratic pressure. Unless that deficit is o(1), (4) allows a nonzero
macroscopic common covariance.

## 4. What the signed-response l1 estimate does and does not remove

The exact full-response estimate separately gives

    sum_(i<j)|C^+_ij-C^-_ij|=O_beta(N^(3/2)).

Since every entry difference is at most two, this implies
||C^+-C^-||_F^2=o(N^2). It does NOT imply that either matrix has small
Frobenius norm: both can share a common positive-semidefinite component.
The exact identity

    ||C^+||_F^2+||C^-||_F^2
       =||C^+-C^-||_F^2+2 Tr(C^+ C^-)

shows that (2), rather than the signed-response estimate alone, is the
relevant additional control. Its present scale is an annealed-energy
deficit, not a sublinear interpolation error.

## 5. Pinning does not automatically upgrade the optimizer quantifier

With quenched labels z, one fixed minimizing A only satisfies an
AVERAGED flip inequality. It does not permit choosing A_e after seeing
z. The Jensen argument controls |E_z R_e(z)| after a radial summation,
not E_z |R_e(z)| or E_z R_e(z)^2. Pinning the global orientation at
sublinear pressure cost therefore does not by itself repair this gap.
If A is reoptimized for the pinned objective it is still a single
signing before z is observed; if A is optimized separately for each
pinned orientation, the problem has changed to separate one-sided
minima. Those minima are trivial at leading order in this campaign.

For endogenous Gibbs pinning labels L, the exact relation is

    log Z_(A^e)-log Z_A
        =log E_(L under A) exp(log Z_(A^e|L)-log Z_(A|L)).

Thus unconditional optimality gives a log-mgf inequality, not
conditional or even averaged conditional optimality. Any proposed
pinning upgrade must account for this likelihood/entropy term rather
than silently applying the edge inequality inside each pinned law.

There is already an exact homogeneous triangle counterexample to the
conditional upgrade. All triangle signings have the same width pressure:
switching leaves only the triangle product, and changing that product
interchanges the two branches. Thus the all-positive triangle is an
actual global width minimizer for every raw lambda>0. Its positive-branch
edge correlation is

    C^+_e=(exp(4lambda)-1)/(exp(4lambda)+3)
             > tanh(lambda).

The inequality is strict because, on writing q=exp(2lambda)>1, the
difference is 2(q-1)^2/[(q^2+3)(q+1)]. Hence conditioning on the positive
orientation does not preserve the optimizer edge inequality, even at
an exact globally minimizing homogeneous width signing. With a symmetric
quenched global orientation field, all triangle classes remain tied by
the same switching/global-sign symmetry. This is a finite quantifier
counterexample, not an asymptotic dense-sign obstruction to a weaker
pinning theorem with a separately proved error budget.
