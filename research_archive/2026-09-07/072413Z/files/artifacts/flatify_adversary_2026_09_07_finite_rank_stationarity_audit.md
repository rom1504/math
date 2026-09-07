# Independent audit: signed finite-rank surgery and stationarity

Date: 2026-09-07. Verdict: PASS for
`flatify_construct_2026_09_07_finite_rank_surgery_stationarity.md`.
This is a valid actual-sign, same-order operation and a genuine consequence
of global cap optimality, not an assumed covariance law.

Subsequent stronger implementation: see
`flatify_adversary_2026_09_07_target_contraction_surgery.md`. Contracting the
target before rounding eliminates both masking and Schur bias, removes
Grothendieck from the operation, and improves fixed-budget error to
O(n^(-1/4)). The proof audited below remains valid as a historical route.

## 1. Mask and signed majorant

Write G=UTU^T, Pi=UU^T, ||T||op<=1. If Z removes b leverage scores above
mu, then b<=r/mu. For a Boolean x,

    |x^TGx-(Zx)^TG(Zx)|
      <=||x-Zx|| (||x||+||Zx||)<=2sqrt(nb).

The diagonal discrepancy is separately bounded by
sum_(removed i)|G_ii|<=sum_i Pi_ii=r. Thus the normalized masking cost
is at most theta sqrt(b/n)+theta r/(2n). This works for signed T; no
PSD assumption on G is hidden in the estimate.

Diagonalize T, set z_l=Zv_l, and K=sum_l |lambda_l| |z_l||z_l|^T.
Then K is PSD, |(ZGZ)_ij|<=K_ij, and K_ii<=mu, hence K_ij<=mu.
The proposed flip probability p_ij=t(K_ij+A_ij(ZGZ)_ij)/2 lies in
[0,tmu] and yields mean A-t offdiag(ZGZ)-t(A circ K). Both the
nonnegative majorant and the Schur bias are necessary and correctly paid.

## 2. Reconstructed diagonal-majorant argument

For symmetric hollow A define beta(A)=max_(x,y Boolean)|x^TAy|. Hollow
polarization gives beta(A)<=4Q(A), by writing the bilinear form as the
difference of two quadratic forms at points of the cube.

The SDP minimizing Tr D over diagonal D>=A and D>=-A has dual

    max <A,X-Y>,  X,Y PSD, diag(X+Y)=1.

Strict primal feasibility gives equality. Write X_ij=<a_i,a_j> and
Y_ij=<b_i,b_j>. The unit vectors u_i=(a_i,b_i), v_i=(a_i,-b_i)
satisfy <u_i,v_j>=X_ij-Y_ij. The real bilinear Grothendieck inequality
therefore bounds the dual objective by K_G beta(A). Its finite matrix
form is Theorem 1.1 in [Pisier's survey](https://webusers.imj-prg.fr/~gilles.pisier/grothendieck.UNCUT.pdf),
independently checked during this audit. No symmetric-Grothendieck variant
with a different constant is being substituted.

Consequently some nonnegative diagonal D has D>=+/-A and
Tr D<=K_G beta(A). Schur multiplying these two PSD inequalities by K
gives D circ K>=+/-(A circ K). Since D circ K is diagonal,

    Q(A circ K)<=1/2 sum_i D_ii K_ii<=2K_G mu Q(A).

This validates the bias term 2K_G theta C mu sqrt(n) after normalization.

## 3. Rounding, edit count, and the rank threshold

The expected number of changed edges is at most tnr/2, because
sum_l |lambda_l| ||z_l||_1^2<=nr. Centered coefficients are bounded by
2 and have total variance at most 4sum p_ij<=2tnr. Bernstein with
a=(n+2)log2 and threshold sqrt(4tnra)+4a/3 controls every Boolean
quadratic response simultaneously with probability at least 1/2.
Markov bounds the probability of more than 2tnr edits by 1/4, so both
properties hold for an actual full signing with positive probability.

Combining these terms gives exactly the displayed e_n. Put
alpha=r/sqrt(n), mu=alpha^(1/3)/sqrt(n). The masking and Schur terms
are respectively theta alpha^(1/3) and 2K_G theta C alpha^(1/3),
while the rounding term is 2sqrt(theta(a/n)alpha). The probabilities
are valid when theta alpha^(1/3)<=1. Thus r=o(sqrt(n)) and fixed theta
give a uniform o(1) normalized error without any incoherence assumption.
For fixed r the leading explicit rate is O(n^(-1/6)). The trivial r=0
case is handled separately rather than inserting alpha=0 into the formula.

## 4. Global optimality and minimax: the quantifiers are sound

For every T the rounding construction gives some actual A'(T). If A is
globally optimal, Q(A'(T))>=M_n. A common A' for all T is not needed:
the error bound uniform in T establishes the lower bound on the minimum
over T of the perturbed cap.

Let p(x)=U^Tx/sqrt(n), h(x)=H_A(x)/n^(3/2), q=Q(A)/n^(3/2).
Replacing the off-diagonal form by the full form changes each response
by at most theta r/(2n), since Tr G=Tr T and |Tr T|<=r. Hence

    min_(||T||op<=1) max_(sigma,x)
       [sigma h(x)-theta sigma p(x)^T T p(x)/2]>=q-epsilon.

Here epsilon includes near-optimality excess, surgery error, and this
additional diagonal payment. Maximizing over probability laws on the
finite signed spin set makes the objective bilinear on compact convex
sets. Minimax and operator/nuclear norm duality then give a law mu with

    E[q-sigma h]+(theta/2)||E[sigma p p^T]||_*<=epsilon.

The first term is nonnegative by the definition of q. This proves both
near-ground concentration and signed covariance balance simultaneously.
For the conditioned law, Markov gives bad mass at most epsilon/eta;
||p p^T||_*=||p||^2<=1 bounds the removed covariance contribution.
Dividing by the retained mass gives precisely the stated conditional bound.

This law depends on the chosen subspace and is not claimed to work
simultaneously for every subspace. It need not have equal polarity masses,
need not be isotropic, and gives no upper control of added-row absolute
response. Those explicit limitations are essential; no cross-order
recurrence follows from the stationarity statement alone.
