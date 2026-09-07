# Fifth checkpoint: an actual all-order subspace bridge and its certificate limit

2026-09-07, continuing six-hour campaign. Not a stopping point.

## Construction actually achieved

For arbitrary prescribed subspaces X,Y of total dimension r>=1 in R^n,
there is an actual full sign bridge C with

    |x^T C y| <= sqrt(n)||P_Xperp x||||P_Yperp y||
                      +O(n^(5/4)*sqrt(r))

uniformly for Boolean x,y, at EVERY sufficiently large order. Thus
r<=n^(1/2-2delta) gives a power-saving O(n^(3/2-delta)) remainder.
This is not merely an assumed profile or a Gaussian model.

Proof dependencies: low-rank orthogonal localization into sparse coordinate
spaces; nuclear-magnitude diagonal contraction; independent sign rounding
whose total variance is O(n*||Delta||_*); uniform cube discrepancy rather
than only a cap bound; zero embedding into nearby Paley Hadamard orders.
Baker--Harman--Pintz Theorem 1 supplies h-n=O(n^.525), with its actual
primary statement independently read by the director. Both prime congruence
classes are handled algebraically, so no prime-in-progressions input occurs.

Canonical proofs:

- `flatify_construct_2026_09_07_arbitrary_subspace_bridge.md`;
- `flatify_adversary_2026_09_07_rectangular_orthogonal_rounding_audit.md`;
- `flatify_adversary_2026_09_07_power_saving_all_order_localization.md`.

## Sharp limitation proved on actual children

For EVERY Q(A)=O(n^(3/2)) signing and EVERY rank-o(sqrt(n)) subspace,
the near-orthogonal spin set retains energy half-width at least

    [4/(3*pi*sqrt(3))-o(1)] n^(3/2).

This uses an explicit correlated Gaussian-sign witness, the elementary
bound ||A||op^2<=beta(A)<=4Q(A), a positive Schur-series covariance
bound, and conditioning paid by |cross energy|<=Q(A).
See `flatify_director_subspace_residual_width_2026_09_07.md`.

Consequently ANY certificate that separately maximizes child energy
envelopes and the above residual-norm bridge bound costs at least
1+8/(3*pi*sqrt(3))>1.4900 in n^(3/2) units. The equal-child target is
below 2sqrt(2)*.493608094<1.397. This is a certificate obstruction, not
a lower bound on the actual bridge-completed signing. Joint alignment
or global old-edge changes remain allowed and unresolved.

## Correction to the previous local assessment

Recent Gaussian/paired-profile sector results are mathematically valid but
unnecessarily restrictive for their fixed-center energy windows. A plain
dephased Hadamard bridge already gives

    |xCy|<=n^(3/2)sqrt((1-alpha^2)(1-beta^2))+O(n),

for ALL spins, with alpha,beta their overlaps with arbitrary chosen centers.
For energies bounded by c*alpha^2 and c*beta^2 with c<=1/2, the total
coefficient is at most 1. Fourier regularity is not the obstruction for
these windows. Whole-child high-energy coverage is unpaid. The older
conditional calculations and failed experiments are preserved, not erased.

Further verified work: a bounded-degree Boolean bilinear tile classification;
balanced-slice width lower bounds; fixed-partition localization; precise
scopes of seed-spread and tensor floors; actual optimizer noisy-shell
variance and an exact optimal-order-5 PSD-slack counterexample. These do
not improve the original recurrence and are not counted as doing so.

## Updated frontier and immediate continuation

The original interval remains [.4333221116640807,.493608094).
No convergence or nonconvergence proof exists. The selected-child weighted
block flatification with O(N^(3/2-delta)) total cap loss remains open.

Next: attack actual JOINT child/bridge values or a genuinely dense global
operation, not finer scalar projection envelopes. Preserve all research
at this checkpoint and continue substantive work through the campaign.
