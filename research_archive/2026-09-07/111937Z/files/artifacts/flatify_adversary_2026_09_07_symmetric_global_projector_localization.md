# Symmetric global full-sign projector operation

2026-09-07. Actual same-order symmetric sign recovery around a symmetric Hadamard seed. This may change ALL old edges. It is not a recovery theorem for the original weighted two-child target.

## 1. Exact compatible-order statement

Let H be a symmetric full sign Hadamard of order n, so H²=nI, and let P be ANY prescribed orthogonal projection of rank r>=1. If a=r log n/n tends to zero, there exists a hollow symmetric full signing C with, simultaneously for every Boolean x,

    |H_C(x)|<= (sqrt(n)/2)||(I-P)x||²
                  +O(n^(3/2)(r log n/n)^(1/4)+n),     (1)

where H_C(x)=x^T C x/2. Constants are absolute. The construction selects one diagonal sign gauge of H, projects the real target on both sides, contracts amplitudes, and independently rounds symmetric entry pairs. Thus r=o(n/log n) is allowed. Rank zero is the hollowing of H itself.

## 2. A one-gauge relative entry bound

Let L=log(16n²), s be a uniform sign vector, H'=D_s H D_s, and T=(I-P)H'(I-P). Then H' is symmetric Hadamard, T is symmetric, and ||T||op<=sqrt(n).

For every i,j, the absolute value of (PH')_ij is the absolute value of the Rademacher sum sum_u P_iu H_uj s_u. Its squared coefficient norm is P_ii. Hoeffding therefore gives simultaneously

    |(PH')_ij|<=sqrt(2L P_ii)                          (2)

except with probability at most 2n²exp(-L). The corresponding H'P bound follows by symmetry, without another failure event.

For the double term write (PH'P)_ij=s^T M s, where M_uv=P_iu H_uv P_vj. Let v=sqrt(P_ii P_jj). The Frobenius norm of M is exactly v; symmetrizing cannot increase it. Its diagonal trace has absolute value at most v. Removing its diagonal gives a hollow symmetric matrix N with ||N||F,||N||op<=v.

Here is an elementary quadratic-chaos bound with the constants needed. For Z=s^T N s/2, random-bipartition Jensen gives Z as the mean of twice a cross bilinear sum. Hoeffding in one side and Gaussian linearization of the positive quadratic exponential in the other give

    log E exp(tZ)<=2t²v²,  |t|<=1/(sqrt(8)v).

The cross-matrix squared Frobenius norm is at most ||N||F²/2; this proves the displayed determinant bound exactly as in `flatify_adversary_2026_09_07_bent_basis_gauge_screen.md`. At t=1/(sqrt(8)v), Chernoff implies

    P(|s^T N s|>12Lv)<=2exp(-L), L>=1.

The bound is trivial when v=0. Adding the diagonal trace and taking a union bound yields

    |(PH'P)_ij|<=13L sqrt(P_ii P_jj)                  (3)

simultaneously, except with probability at most 2n²exp(-L). Thus (2)--(3) hold with probability at least 3/4.

Fix such a gauge and put a_i=sqrt(13L P_ii). Since every H'_ij is a sign,

    |T_ij|<=1+sqrt(2L)(sqrt(P_ii)+sqrt(P_jj))
                       +13L sqrt(P_ii P_jj)
             <=(1+a_i)(1+a_j).                       (4)

## 3. Symmetric contraction and sign recovery

Set D_ii=(1+a_i)^(-1), and B=DTD. This is symmetric and entrywise in [-1,1]. The operator bound on T and sum_i(1-D_ii)²<=13Lr give

    beta(T-B)<=2n sqrt(13Lr).                         (5)

The exact projection deficit is n²-||T||F²=2nr-||PH'P||F²<=2nr. Every row and column of T has squared norm at most n. Therefore the full-matrix variance budget is

    V=n²-||B||F²
      <=2nr+2n sum_i(1-D_ii²)
      <=2nr+4sqrt(13)n sqrt(Lnr).                    (6)

Independently for i<=j, choose a sign Cfull_ij with mean B_ij and copy it to the symmetric position. For a Boolean bilinear pair x,y, the coefficient of a random off-diagonal entry in x^T(Cfull-B)y is x_i y_j+x_j y_i, of magnitude at most two. Thus its variance is at most 2V and each centered summand has magnitude at most four. Union Bernstein over all 2^(2n) pairs, with u=(2n+2)log 2, gives a realization satisfying

    beta(Cfull-B)<=2sqrt(Vu)+(8/3)u.                 (7)

This explicitly handles the dependence of the two symmetric matrix entries. Equations (5)--(7) give beta(Cfull-T)=O(n^(3/2)(rL/n)^(1/4)+n) when rL/n tends to zero.

Finally remove the diagonal of Cfull to obtain the hollow signing C. Since those diagonal entries are signs, this changes every quadratic value by at most n/2. Also

    |x^T T x|<=sqrt(n)||(I-P)x||².

Combining these facts proves (1). The bilinear error estimate holds on the entire real cube by separate linearity; this will be used for zero-padding below.

## 4. All-order availability is a separate input

The compatible-order theorem only assumes a symmetric Hadamard seed of the stated order. For an all-order extension, it suffices to have symmetric Hadamard orders N>=n with N/n tending to one. Embed P by zeros into N dimensions, apply the theorem there, and restrict the recovered hollow signing to the original n coordinates. Zero-padded spins are allowed by the cube bound above and have unchanged projection norms. The extra leading-factor cost is at most n(sqrt(N)-sqrt(n))/2.

A fully explicit power-saving supply uses ANY ordinary Hadamard H_h, not a prime in a prescribed congruence class. Define a symmetric sign matrix of order N=h² by

    W_[(i,a),(j,b)]=(H_h)_ib (H_h)_ja.

Directly summing over j,b gives W²=h² I=N I. The already independently audited unrestricted-prime Paley/Baker--Harman--Pintz construction supplies an ordinary Hadamard order

    h>=sqrt(n), h-sqrt(n)=O(n^(21/80)).

Consequently N-n=O(n^(61/80)), and the leading-factor cost is O(n^(101/80)). Thus at EVERY sufficiently large order n, with arbitrary prescribed rank r=o(n/log n), an actual hollow symmetric full signing satisfies

    |H_C(x)|<= (sqrt(n)/2)||(I-P)x||²
       +O(n^(3/2)(r log n/n)^(1/4)+n^(101/80)),        (8)

uniformly in Boolean x. For r<=n^(1-delta), the first error is O(n^(3/2-delta/4)(log n)^(1/4)); the second is also power-saving. This all-order argument does NOT assume symmetric Hadamards at the prime-derived ordinary Hadamard order itself. It first squares that order by the explicit W construction.

## Scope

This is a genuine global symmetric sign construction: it need not preserve any internal edges of a seed. Its reference target is a gauged, projected symmetric Hadamard, not the row-regular weighted sum of actual minimizing children. Choosing P=0 gives only the familiar half-scale spectral bound. No favorable original-value recurrence, strict sub-half upper bound, or convergence conclusion follows from this operation alone.
