# Rank-one majorant removes Grothendieck from finite-rank surgery

Date: 2026-09-07. Independently checked by the constructive agent after
proposal. This strengthens/simplifies the audited finite-rank operation.

Subsequently superseded by target contraction before rounding; see
`flatify_adversary_2026_09_07_target_contraction_surgery.md`. That stronger
operation removes this masking step and improves the rate. This note
preserves the intermediate simplification and its correct proof.

Let G be any symmetric operator contraction. Put ell_i=|G|_ii, where
|G| is the spectral absolute value, and tau=Tr|G|. Then

    |G_ij|<=sqrt(ell_i ell_j).

This follows by writing G=|G|^(1/2) sign(G)|G|^(1/2) and applying
Cauchy-Schwarz. For G=UTU^T one may instead use ell_i=(UU^T)_ii and
tau=r, uniformly over all ||T||op<=1.

Mask coordinates with ell_i>mu. On retained coordinates set
d_i=sqrt(ell_i), otherwise d_i=0, and choose the RANK-ONE majorant

    K=d d^T.

It dominates |ZGZ| entrywise, K_ii<=mu, and
sum_(i<j)K_ij<=1/2(sum_i d_i)^2<=n tau/2. Thus every probability,
variance, and edit-count estimate of the original construction survives
with tau in place of r.

But now C circ K=diag(d) C diag(d), for any hollow real coefficient
matrix C. Since max_i d_i<=sqrt(mu), separate-affinity cube contraction
gives the sharper elementary bound

    Q(C circ K)<=mu Q(C).

No Grothendieck theorem or diagonal SDP majorant is required. In the
uniform-sign error formula replace 2K_G C0 by C0; in the weighted
extension replace 2K_G C0/b by C0/b. All other terms stay unchanged.

The leverage mask deletes at most tau/mu coordinates, its diagonal
payment is at most tau, and the operator contraction still bounds the
quadratic masking error. Consequently the construction is controlled by
nuclear budget tau=o(sqrt(n)), NOT necessarily rank. High-rank signed
operators with sufficiently small nuclear norm are allowed.

## A global operator-norm covariance consequence

Taking the entire symmetric nuclear-norm unit ball as perturbations,
rather than an operator-norm ball on a chosen subspace, gives a uniform
error O(n^(-1/6)) for fixed theta at an exact cap minimizer. The mask may
depend on G; uniformity of its bound is all minimax needs.

Nuclear/operator duality therefore gives ONE law mu on signed spin states,
not a separate law for each direction, with

    E[q-sigma h(x)]
       +(theta/2)||E[sigma xx^T/n]||op <= O(n^(-1/6)).

Here h=H_A/n^(3/2), q=M_n/n^(3/2), and theta>0 is fixed. The diagonal
correction is at most theta/(2n), already smaller than the error. This is
a full-space signed covariance statement in OPERATOR norm. It does not
give a small nuclear norm or a small trace: the dimension factor prevents
deducing equal polarity masses or midpoint balance from it.

For the weighted extension, perturbations still must vanish on forbidden
zero edges. The full nuclear unit ball is available only when the profile
has no forbidden edges; supported nuclear balls remain available otherwise.
