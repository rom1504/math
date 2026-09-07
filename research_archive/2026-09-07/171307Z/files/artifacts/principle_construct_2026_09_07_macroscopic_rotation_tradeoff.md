# Exact cap cost of continuously flattening opposite weighted blocks

2026-09-07. A seed-retaining global orthogonal operation with exact row-energy and amplitude accounting. Its worst-case cap loss is sharp and can be leading on actual bounded-cap inputs. This does not exclude a different operation or a favorable selection of exact minimizers.

## 1. An exact amplitude-homogenizing rotation

Let A be a hollow symmetric full signing of order n, and define

    k=sqrt((2n-1)/(n-1)),
    W0=k diag(A,-A).

Every row of W0 has squared norm 2n-1, the full-sign target at order 2n. For 0<=theta<=pi/8, orthogonal block mixing (and an immaterial whole-block sign switch) produces

    Wtheta=[ alpha A   beta A ]
           [ beta A   -alpha A],
    alpha=k cos(2theta),  beta=k sin(2theta).

This is an exact orthogonal conjugate of W0. More importantly, its squared matrix has the same diagonal as W0^2, because both blocks involve the SAME A. Thus it remains exactly row-square-regular, not merely Frobenius-normalized. Every nonmatching entry has amplitude alpha or beta; the n matching entries remain zero. Its maximum amplitude decreases from k to alpha. At theta=pi/8 both amplitudes are k/sqrt(2)=1+O(1/n), so contracting that negligible excess and filling the matching costs only O(n) in cap.

## 2. Exact seed-sensitive energy identity

Write P(A)=max H_A, R(A)=-min H_A and width(A)=(P+R)/2. Then

    Q(W0)=k(P+R)=2k width(A).

For parent spins x,y put u=(x+y)/2 and v=(x-y)/2. Set

    U=H_A(u), V=H_A(v), Z=u^TAv.

The supports of u,v are disjoint and partition the old vertices. Direct expansion gives

    H_Wtheta(x,y)=2[beta(U-V)+alpha Z].                (1)

Consequently the actual excess is governed by the joint induced-energy/cross-energy profile of the ORIGINAL A, not its spectrum or a chosen signed cancellation. At the terminal angle this specializes to the archived clique-flip orbit identity; it does not remove that known payment.

## 3. Sharp universal relative-cap tradeoff

For every hollow symmetric A, without a full-sign assumption,

    cos(2theta) Q(W0) <= Q(Wtheta)
       <=[1+sin(2theta)] Q(W0).                       (2)

The lower bound follows by averaging over whole-child reversal to remove the cross block and then maximizing the remaining difference of child energies.

For the upper bound set

    r=alpha/(k+beta),
    f(a,b)=U a^2+V b^2+Zab.

Both f(1,r) and f(r,-1) are energies of vectors in the original cube. Hollow multilinearity implies that their difference has absolute value at most P(A)+R(A)=2 width(A). Exact algebra gives

    beta(U-V)+alpha Z
       =(k+beta)/2 [f(1,r)-f(r,-1)].

Apply (1) and maximize to obtain Q(Wtheta)<=2(k+beta)width(A), which is the upper bound in (2). This improves the cruder formulation using Q(A): the natural denominator is exactly the cap of W0, even when A is unbalanced.

In terms of the new maximum amplitude K'=alpha, the upper factor is

    1+sqrt(1-(K'/k)^2).

Thus the universal loss permitted by this rotation is of square-root order in a small amplitude decrease, not linear order.

The constant is sharp for the two-dimensional cube-energy constraint. The extremal profile is

    U=-V=1/(1+r^2),   Z=2r/(1+r^2),

whose quadratic has range [-1,1] on the square. It attains beta(U-V)+alpha Z=k+beta. The following construction realizes the same obstruction asymptotically with actual full sign A of bounded normalized cap.

## 4. Macroscopic actual-sign sharpness and leading-loss examples

Fix theta>0 and let r be as above (it depends only on theta). Split n=2m vertices into two equal macro blocks. For a constant L>0 choose a hollow mean matrix

    M=(L/sqrt(n)) [ J_m-I_m    r J_m ]
                   [ r J_m    -J_m+I_m ].

For n sufficiently large every mean entry belongs to [-1,1]. Independently round its upper-triangle entries to signs, symmetrically, obtaining actual A. Hoeffding and a union over 2^n spins show that a realization exists with

    Q(A-M)<=2 n^(3/2).                                (3)

This constant is independent of L. The matrix M is negated by the signed block swap (x,y)->(y,-x), so its cap equals its width. Optimizing its two magnetizations gives, with g=(1+r^2)/2,

    Q(M)=(L g/4)n^(3/2)+O(L sqrt(n)).

The profile u=1 on the first block and v=1 on the second yields

    Q(Wtheta(M))>=2(k+beta)Q(M)-O_theta(L sqrt(n)).

Using (2) on the error E=A-M and (3),

    Q(Wtheta(A))-Q(W0(A))
       >=2beta Q(M)-(4k+2beta)Q(E)-O_theta(L sqrt(n)). (4)

For every fixed theta>0, choose L sufficiently large, still fixed before n grows. Equation (4) is then a positive constant times n^(3/2). Yet A is an actual full signing with Q(A)=O_theta,L(n^(3/2)), W0 and Wtheta are exactly row-square-regular, and the maximum amplitude strictly decreases by a fixed amount.

Moreover, sending the fixed constant L to infinity after fixing theta makes the actual ratio Q(Wtheta(A))/Q(W0(A)) approach the upper factor 1+sin(2theta). Thus the tradeoff is sharp over bounded-normalized-cap actual-input families; the bound on that normalized cap is allowed to depend on the desired proximity to sharpness.

This is a genuinely macroscopic mean structure, not deletion of a small exceptional set. The starting W0 already has the two comparable weighted blocks of interest, and the orthogonal transformation retains all information about A.

## 5. Scope and archive boundary

The terminal theta=pi/8 clique-flip formulation and its stronger counterexamples on actual asymptotic minimizers are already in `decisive_independent_h2_exact_profiles_and_algebra_2026_09_07.md` and `decisive_independent_h2_subhalf_planting_counterexample_2026_09_07.md`. They are not new claims here.

The new calculation is the complete continuous angle/amplitude tradeoff (2), with an exact relative-cap denominator and a macroscopic actual-sign sharpness construction. It shows that this natural global mixing cannot have a universal power-saving cap error throughout the bounded-cap weighted class. It does NOT show that every angle is bad for selected exact minimizers, and it does not rule out other correlated block exchanges or genuinely different global operations.

## 6. Exact seed-sensitive initial cost

For any fixed finite nonzero A let P=max H_A and -R=min H_A. Directly in the two original child spins,

    H_Wtheta(x,y)=k cos(2theta)[H_A(x)-H_A(y)]
                  +k sin(2theta) x^T A y.

The maximum of finitely many differentiable functions has a right derivative equal to the maximum derivative among the active functions. At theta=0, the active positive-energy pairs satisfy H_A(x)=P and H_A(y)=-R; swapping the children treats the negative-energy active pairs. Reversing either whole child leaves its quadratic energy unchanged and reverses x^T A y. Therefore

    d/dtheta Q(Wtheta)|_(0+)
       =2k max { |x^T A y| : H_A(x)=P, H_A(y)=-R }.  (5)

In particular this derivative is nonnegative. The maximum amplitude decreases only quadratically, k-alpha=2k theta^2+O(theta^4). When the opposite-extremizer coupling in (5) is positive, the initial cap payment is linear in the angle and hence square-root in the amplitude decrease. When it vanishes, (5) alone gives no strict improvement or loss: near-active configurations may become active later. The statement is an exact finite identity, not a uniform-in-n derivative interchange, and no lower bound on this coupling for selectable exact minimizers is asserted.
