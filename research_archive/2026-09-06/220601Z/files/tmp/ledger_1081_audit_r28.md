# Audit of ledger section 10.81.1 and equations (10.852)--(10.855)

## Verdict

The mathematical content and all asymptotic exponents are correct. The
random-hash constants agree with the independent planted-partition proof.
Equation (10.855) is also correct once its implicit orientation is made
explicit. There is one optional factor-two sharpening and two wording
clarifications; none changes the frontier.

## Section 10.81.1

1. Equation (10.848) is exact. The saturation identity is
   Pr(hit)=|E_T K_P|/2^n. Since E_T=-E_T and global sign lies in K_P,
   there are only 2^(k-1) effective translates, so the amplification factor
   is correct. Randomizing P cannot change this pointwise-in-P bound.

2. The conditioning inequalities are correct. If h=Pr(hit) and the row-good
   event has mass g at least 1/2, then Pr(hit|G) is at most 2h.
   Inclusion--exclusion, minimized over g in [1/2,1], gives the lower bound
   max(0,2h-1).

   Claim-strength clarification: the zero lower bound shows that the
   half-mass row theorem cannot prevent complete erasure. It does not itself
   construct an A,S,P for which actual complete erasure occurs. Thus
   “conditioning can therefore erase” is safest as “the theorem cannot
   exclude conditioning from erasing.”

3. Equation (10.849) is exact in the stated normalization. The Frobenius norm
   counts ordered off-diagonal entries:

   m(m-1)(1-p_2)^2 +
   [n(n-1)-m(m-1)]p_2^2 = n(n-1)p_2(1-p_2).

   Principal compression gives ||J_S A J_S||op <= ||A||op, so the operator
   bound is correct.

4. Equation (10.850) is valid but has an optional factor-two sharpening. As
   written, it unions the two-sided Hanson--Wright bound (prefactor 2) over
   all 2^k block words, producing 2^(k+1). Because z and -z give the same
   quadratic form, one may union over 2^(k-1) representatives and replace
   2^(k+1) by 2^k. The displayed bound remains true and the difference is
   asymptotically immaterial. For literal completeness, declare T >= 0 before
   invoking Hanson--Wright.

5. Equation (10.851) and its exponents are correct. With
   L_n=n^(3/4-c), k=o(L_n), and T=omega(n^(3/2-c)), the linear exponent is
   omega(L_n). The quadratic exponent is at least omega(n^(1-2c)), hence
   omega(L_n), because c<1/4. At T >= eta n^(3/2), the linear exponent is
   Omega(n^(3/4)).

   Wording clarification: “exactly scale-matched” is most precise at
   T=Theta(n^(3/2-c)). For the larger class T=O(n^(3/2-c)), the coarse bound
   simply does not rule out success.

6. The conclusion is properly scoped to uniform random diagonals. It does not
   falsify an adapted common row-good law in (10.838).

## Random hashing, (10.852)--(10.855)

1. Equation (10.852) is exact. With bar_e=k^-1 1, the mean is
   (Ax)bar_e^T. The summand norm is sqrt((n-1)(1-1/k)), and the two variance
   matrices are exactly those displayed.

2. The existence constants in (10.853) check. Rectangular Bernstein with
   u=log(8(n+k)) fails with probability at most 1/8. Binomial Chernoff gives
   block-balance failure at most 2k exp[-n/(12k)] <= 1/8 under
   n/k >= 12 log(16k), so the events intersect. From
   ||(Ax)bar_e^T||op^2=R_2(x)/k and (a+b)^2 <= 2a^2+2b^2, the constants are
   precisely 2, 8, and 16/9.

3. Equation (10.854) has the correct powers. For
   k asymp n^(3/4-c)/log n and c<1/4,

   n^2/k=O(n^(5/4+c)log n)=o(n^(3/2)),

   so nu=O(n^(3/2)), k nu u=O(n^(9/4-c)), and
   k L^2 u^2=O(n^(7/4-c)log n), which is lower order.

4. Equation (10.855) is correct, with an implicit orientation. If x_S=y
   underlies an exact oriented child optimizer, choose sigma so that
   sigma y^T A[S]y=Q(A[S]). Then

   |x^T H_S x| >= sigma x^T H_S x
   = Q(A[S])-p_2 sigma x^T A x
   >= Q(A[S])-p_2 q_n.

   This verifies the sign and projective/orientation factor. It would be
   clearer to include sigma in the prose or this one-line derivation, but the
   displayed inequality is valid.

5. The quantifier warning is accurate: constructed P,D_x depend on (S,x).
   The result is local and proves neither one common coset in (10.837) nor the
   common-law statement (10.838).
