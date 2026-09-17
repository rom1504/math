# Linear block sharpness with bounded normalized full-sign cap

2026-09-17. New deduction, pending independent audit. This strengthens
critical_full_cap_embedding.md: the deterministic child bias is now
only of order m^(-1/2), so the compared caps are O(N^(3/2)), not N^2.
The constants are deliberately nonoptimal and need not be close to the
minimizing constant.

## 1. Input row law and a spectral-moment lemma

The director's tilted row law in dimension m has exact mean zero,
covariance I_m, an absolute Euclidean L-subGaussian constant, and

    E|sum_j epsilon_j|-E|sum_j Z_j| >= eta sqrt(m)          (1)

for one fixed eta>0 and all sufficiently large m. The ordinary sign
law also has these properties with L=1. Use a common L>=1 for both.

If B consists of n independent rows of either law, then

    E ||B||op^2 <= 32 L^2 (n+m).                         (2)

Indeed, for fixed unit x,y, x^T B y is L-subGaussian. Two 1/4-nets
have combined size at most 9^(n+m), and their bilinear maximum is at
least ||B||op/2. Hence

    P(||B||op^2 > 8 L^2[(n+m)log9+log2+s]) <= exp(-s).

Integrating gives (2). No independence within a row is used.

## 2. Deterministic weakly planted right child

Fix constants K>=1 and lambda>0, set m=Kn (integer K), and let
rho=lambda/sqrt(m), which lies in [0,1] for all large m. There exists
a symmetric hollow full-sign D such that

    D = rho(J-I)+R,       ||R||op <= C sqrt(m), C=8.       (3)

To prove existence, take independent upper-triangle signs of mean rho.
For each unit u, the centered quadratic u^T R u has MGF <=exp(t^2).
A 1/4-net of size 9^m and the symmetric-net inequality give positive
probability of ||R||op<=8sqrt(m). This is a direct bounded-entry
argument, not an imported random-matrix premise.

Write beta=H_D(1). From (3),

    |beta-rho m(m-1)/2| <= (C/2)m^(3/2).                 (4)

Let A be ANY n-vertex full-sign child with Q(A)<=C_A n^(3/2).
It may in particular be an exact minimizer. Set

    W(B) = [ A  B ; B^T D ],     N=n+m.

Every entry of W(B) is a sign. D and A are the SAME for the two row
laws being compared.

## 3. Positive maximum: quantitative approximate pinning

For a right word y, use simultaneous reversal of x,y to arrange that
its negative set S has size k<=m/2; write s=1_S. Then

    H_D(y)-beta
       = -2rho k(m-k)-2s^T R 1+2s^T R s
       <= -lambda sqrt(m)k+2Cm sqrt(k)+2C sqrt(m)k.

For lambda>=4C this is at most

    -(lambda/2)sqrt(m)k+2Cm sqrt(k).

The bridge change from y=1 is at most
2sqrt(n)||B||op sqrt(k). Optimizing the resulting quadratic in sqrt(k)
therefore proves, for P(B)=max_(x,y) H_W(x,y) and r=B1,

    beta+sum_i|r_i|-Q(A) <= P(B)
      <= beta+sum_i|r_i|+Q(A)
         +4C^2 m^(3/2)/lambda
         +4n||B||op^2/(lambda sqrt(m)).                 (5)

The lower bound simply chooses y=1 and x_i=sign(r_i).

## 4. Negative polarity costs no leading uncontrolled term

Let N_-(B)=max_(x,y) -H_W(x,y). Since
-rho[(sum y)^2-m]/2<=rho m/2,

    N_-(B) <= Q(A)+rho m/2+(C/2)m^(3/2)
                    +sqrt(nm)||B||op.

For m>=4 and lambda>=16(C+C_A+1), (4) implies

    beta-2Q(A)-rho m/2-(C/2)m^(3/2)
                                      >= (lambda/16)m^(3/2).

Using P(B)>=beta-Q(A) and (u-t)_+<=u^2/(4t), we obtain

    (N_-(B)-P(B))_+ <= 4n||B||op^2/(lambda sqrt(m)).     (6)

Consequently the ABSOLUTE, unrestricted full-sign cap satisfies

    beta+sum_i|r_i|-Q(A) <= Q(W(B))
      <= beta+sum_i|r_i|+Q(A)
          +4C^2m^(3/2)/lambda
          +8n||B||op^2/(lambda sqrt(m)).                (7)

All bounds before taking expectation hold pointwise in the sign bridge.

## 5. Leading-order separation and its scope

Let W_0 use iid fair sign rows, and W_1 use the tilted rows. Combining
(1), (2), and (7), with m=Kn, gives

    E Q(W_0)-E Q(W_1)
      >= eta n sqrt(m)-2C_A n^(3/2)
         -[4C^2+256L^2(K+1)/K^2]m^(3/2)/lambda.         (8)

The last displayed coefficient follows exactly from
8n/(lambda sqrt(m))*32L^2(n+m).

First choose fixed integer K so that 2C_A/sqrt(K)<=eta/4.
Then choose fixed lambda sufficiently large that it meets Section 4
and the last term in (8) is at most (eta/4)n sqrt(m). Thus

    E Q(W_0)-E Q(W_1)
        >= (eta/2)n sqrt(m)
        = [eta sqrt(K)/(2(1+K)^(3/2))] N^(3/2).         (9)

Both models use n independent random blocks of size m=Theta(N),
with EXACT matching mean and covariance I_m and a common uniform
subGaussian constant. All child entries are common deterministic
full signs. Therefore (9) disproves an o(N^(3/2)) all-offset cap
universality statement at linear block size under those assumptions.
Comparison of both models to their common Gaussian reference would
contradict (9), so no direct Gaussian pinning analysis is required.

In contrast to the earlier hard-ferromagnet embedding,

    E Q(W_a)=O_(K,lambda,L,C_A)(N^(3/2)),  a=0,1.

The same order bound holds with probability 1-exp(-cN), by the net
tail for ||B||op and the deterministic inequalities above. This is a
bounded-normalized-cap example in expectation and with high probability;
the full support can contain exceptional high-cap bridges.

The fixed constants K and lambda are large. The construction is NOT
near-minimizing and does not refute universality confined to optimizing
or quantitatively near-optimizing parent matrices. The left child can
be an actual minimizer, but the right child is deliberately planted.
