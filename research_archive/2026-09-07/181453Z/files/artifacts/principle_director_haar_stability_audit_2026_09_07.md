# Independent reconstruction of the Haar stability falsifier

2026-09-07. **Verified**, independently reconstructed from the stated model,
not inferred from a prior audit. This is about Haar weighted involutions, not
an original full-sign construction or the optimal signing constant.

For even n let U=O diag(I,-I) O^T with equal multiplicities and Haar O.
Write r=x^T Ux/n for a fixed Boolean x. Orthogonal invariance gives
(1+r)/2~Beta(n/4,n/4), hence the upper-tail rate
I(r)=-log(1-r^2)/4 at speed n. Conditional on r, the component of
Ux/sqrt(n) perpendicular to x is exactly a uniform spherical vector of
length sqrt(1-r^2). This follows from the full stabilizer of x, and does
not require independent entries or an asymptotic Gaussian replacement.

After switching x to 1, represent that sphere by
(g_i-gbar)/sqrt(n q), q=n^-1 sum(g_i-gbar)^2. A maximizing Boolean x
for the hollow quadratic must obey x_i(Ux)_i>=U_ii. The event
max|U_ii|<=1/1000 has probability tending to one by the same Beta tails.
For 0<=r<=31/32, stability therefore implies all standardized projected
Gaussian coordinates exceed -(31+32/1000)/sqrt(63).

On |gbar|<=.01 and q<=1.01^2 this implies g_i>=-4 for every i. The
probability is at most Phi(4)^n plus the two exceptional probabilities.
The elementary Mills bound, normal mean tail and chi-square Chernoff bound
in the companion proof give a uniform upper bound 4 exp(-n/50000).
No conditioning or independence of the diagonal event and Ux is used:
the diagonal event merely supplies a necessary inequality on Ux.

Set delta=1/50000 and r0=sqrt(1-exp(2delta)/16). Then
I(r0)=log2-delta/2. Union-bounding ONLY stable candidates with
r in [r0,31/32] gives exp(-delta n/2+O(log n)). Larger r already has
I(r)>log2 and is excluded by the ordinary tail bound. Applying the same
argument to -U handles the negative extreme. Since Tr U=0, hollowing
leaves every Boolean energy EXACTLY unchanged. Since ||U||op=1,
boundedness also gives the claimed expected-cap limsup.

Thus the normalized half-energy cap is asymptotically at most
sqrt(1-exp(1/25000)/16)/2, strictly below sqrt(15)/8. All constants are
fixed before n tends to infinity. There is no interchange of an optimizing
threshold with n, nor a replica or metastability-counting premise.

The identity relating the old Gaussian terminal functional to the fixed-
vector Haar log moment-generating function explains the shared old number
sqrt(15)/8. The strict improvement proves that coincidence is annealed
slack, not an exact ground-state primal/dual theorem. Retaining joint local
stability is a substantive missing condition; it is not yet an operation
that transfers arbitrary optimal signings between orders.

Full proof and exact elementary constants:
`principle_invent_2026_09_07_haar_stability_strict_gap.md`.
