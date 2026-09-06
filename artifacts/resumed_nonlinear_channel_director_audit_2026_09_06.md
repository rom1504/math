# Director reconstruction: nonlinear channel and feasible center update

Date: 2026-09-06. This is an independent reconstruction of the proof in
`resumed_bound_audit_restricted_channel_center_update_2026_09_06.md` and its
deterministic forest module, not an appeal to the author's positive audit.
The original convergence question is not settled by this theorem.

## Structural checks

1. Old kernel global flattenings are bounded at fixed operator cap. For a
   proper fixed-root cut, assign the external root to the side containing its
   child: its flat edge is internal and the fixed root removes the otherwise
   present sqrt(n) isolated-index factor. This gives the extra n^(-1/2).
   Collision deletion has fixed-root Hilbert norm O(n^(-1/2)).
2. The nonlinear forest argument really needs at least two branches. If no
   branch is split by a cut, its two nonempty whole-branch root maps factor
   through diag(B_i). If one branch is split, at least one whole branch
   supplies the bounded root map; its stacked split tensors have norm
   O(n^(-1/2)). If two or more are split, triangle inequality already suffices.
   Branch product root maps are bounded by the positive Schur Gram inequality.
   This verifies all proper cuts, not only a subset of natural cuts.
3. Local Gaussian Hermite order and physical input-spin chaos order must be
   distinguished. The local Wick polynomial has a pure top input-chaos kernel
   plus a rootwise vanishing L2 error: whole-copy contractions give its Hermite
   subtractions and every proper contraction is small. The same leading-spin
   pairings prove the Rademacher local replacement; one must not merely evaluate
   a repeated-index Gaussian Wick polynomial on signs.
4. Fixed-L averaged L2 transport preserves those errors. Repeated marked-slot
   deletion commutes with the final matrix-row summation, and costs averaged
   Hilbert norm o(1). Proper singleton flattenings bound influences of the
   remaining squarefree nonlinear kernels. Own-output-spin deletion is harmless
   for these kernels, but not for first-chaos transport.
5. Multivariate chaos normality follows only after covariance subsequences are
   extracted. I inspected Noreddine--Nourdin Theorem 1.1 in the primary PDF:
   it requires fixed chaos orders, covariance convergence, and fourth-cumulant
   convergence, and does not require nonsingular limiting covariance. The
   forest bounds give its hypotheses. This does not Gaussianize B X_T.
6. For the latter exceptional term, B X_T is approximated by Q[S h_T(X)].
   Its full contraction with a star channel splits into a root-not-hit and a
   root-hit term. Before the covariance-of-squares argument, replace h_T by its
   pure top chaos as in item 3; its small fixed-degree error stays small after
   contraction. Joint old-at-any-root/star independence then kills every full
   contraction into this child polynomial. Multiplication by the omitted root
   input is controlled by the degree-bounded creation inequality, not a sum
   of n absolute errors.
7. The creation inequality follows by Hermite coefficients: because P_a is
   own-a-free, multiplying by N_a only creates an index of multiplicity one.
   Each output multiindex has at most d+1 possible origins a. Cauchy--Schwarz
   therefore gives the factor d+1 times sum_a q_a^2 ||P_a||_2^2.
8. The root-hit term is exactly the diagonal of Q(B circ K_(r-1))B, up to a
   fixed normalization. The inherited directional-derivative matrix theorem
   gives polylogarithmic operator moments. The flat Schur factor converts
   its Frobenius norm to polylog(n), hence its normalized diagonal L2 norm
   vanishes. No uniform-root conclusion is needed.
9. Every other mixed contraction is proper in the Gaussianizing star kernel
   and vanishes by Cauchy--Schwarz. Wick expansion then factors polynomial
   moments from the whole background block (old fields, B X_T). Zero linear
   covariance alone would not establish this step.
10. Gaussian-to-sign replacement uses squarefree Y=B X_T only linearly.
    In its fourth derivative there are just Y D^4 f and 4(DY)D^3 f.
    The test fields have coordinate derivative bounds n^(-k/2); the total
    influence of Y is degree times its L2 norm. Summing the replacement
    errors costs O(n^(-1))||Y||_2, with no low-influence assumption on Y.
    Thus coherent Q S atoms do not invalidate this expectation comparison.
11. For the marked test, after removing own S_i from the star channel,
    the exact two-spin difference identity makes each off-diagonal term
    O(n^(-1)); the Q row has l1 norm O(sqrt(n)). The own-root term is kept
    exactly with S_i^2=1. The separate Gaussian-IBP module gives a consistent
    independent derivation, including the exceptional cubic third derivative.

## Bounded response and limit-order checks

The deterministic weighted linear theorem identifies the averaged absolute
covariance error by choosing its sign as a deterministic bounded weight.
Together with the newly proved joint laws and mixed contractions this gives
the two restricted nonlinear identities; no unrestricted conditional-law
claim is necessary.

For bounded feasible F,H, polynomial approximants are used only for analysis.
The matrix norm controls their energy error at each fixed L. Keep the finite
channel fixed while approximating F. Adding delta times an independent Gaussian
to the channel bounds all its comparison variances away from zero; a single
Gaussian density dominates this compact variance family. This justifies uniform
bounded-test approximation even when some original v_i vanish.

The actual common-center means +/-F+H sign(BF) are exactly cube-feasible.
Their half energy difference is sum_i H_i |(BF)_i|, so there is no quadratic
remainder. The softsign test loses at most eta+2s delta E|N| uniformly in
the local variance. After the fixed construction's matrix limit, send analysis
approximation errors, eta, and delta to zero. Jensen in sqrt(v_i), and the
proved Schur mean-standard-deviation inequality, remove the row dependence.

Finally fix a principal deletion fraction, apply the theorem at its resulting
fixed L, and then send the fraction to zero. The Gaussian lower functional
is independent of L. This order is valid even though finite approximation
sizes and error constants depend on L.

## Verdict

This reconstruction finds the restricted nonlinear-channel proof and feasible
center-update inequality sound, subject to the explicit pure-top-chaos
clarification in item 6 being retained in the source. The additional independent
response-agent reconstruction is complete and agrees; its full checks are in
`resumed_response_center_theorem_independent_audit_2026_09_06.md`.
The result supplies one further feasible universal lower-bound mechanism; it
does not license repeating the update with its now-dependent center spin, does
not construct an all-order upper bound, and does not prove convergence.

I also inspected the primary [Nourdin--Rosinski Theorem 3.4 and equations
(3.5)--(3.6)](https://arxiv.org/pdf/1112.5070). They require all mixed
unsymmetrized contractions, not merely symmetrized ones. The proper-star
estimates control those norms. In the exceptional full-star contraction,
the remaining slots are symmetric already, so the derivative calculation
controls the required full norm. The conclusion needed is block polynomial
moment independence; no moment-determinacy assumption on the exceptional
first-chaos transport is being silently added.
