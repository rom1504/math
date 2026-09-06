# Independent response-track audit of nonlinear center transport

Date: 2026-09-06. The response agent independently read and reconstructed
`resumed_bound_audit_restricted_channel_center_update_2026_09_06.md` and
`resumed_bound_audit_nonlinear_forest_flattening_2026_09_06.md` in full.
Result: the finite-response center theorem passes this adversarial audit.
No prior audit was treated as a substitute for the following checks.

## 1. Forest transport and the first-chaos exception

I reconstructed all three deterministic cut cases. With zero split branches,
the two nonempty whole-branch root maps flank `diag(B_i)`. With at least
two split branches, the triangle bound gains two factors `n^-1/2` against
only `||B_i||_1=O(sqrt(n))`. With one split branch, at least one other
whole branch exists; the stacked or block-diagonal factorization uses its
bounded root map. The shared root label is retained in these factorizations.
Repeated branch types are controlled by their positive Schur Gram products.

The resulting small proper flattenings apply to every transported higher
LOCAL Hermite monomial. I checked that local Hermite degree and original
input chaos order are kept distinct. Local Wick replacement cancels whole
equal-branch contractions, while proper and unequal-branch contractions
have vanishing fixed-root errors. The errors are transported only in
averaged L2 with fixed `||B||op`; a false uniform-root transport estimate
is not used.

The first local chaos is correctly separated. Its transport can contain
coherent own-spin atoms through `Q[S h_T(X)]`; the proof never declares
this entire vector Gaussian.

## 2. Full contractions, including the essential pure-chaos ordering

For the first-chaos part, proper-star contractions vanish using the star's
own small self-contractions and Cauchy--Schwarz. The exceptional full-star
contraction is genuinely treated separately.

The root-not-hit argument requires the following order:

1. Replace `h_T(X_a)` by its single pure TOP input-chaos Wick kernel.
2. Obtain old-family/star joint moment independence at arbitrary roots
   `(a,i)` from their joint Gaussian limit and vanishing cross covariances.
3. Apply positivity of the covariance-of-squares contraction formula to
   these two PURE multiple integrals.
4. Restore the `o(1)` replacement error using bounded fixed-degree
   derivative/contraction operators.

Applying the positive formula directly to a general sum of chaoses would
not be justified. The detailed proof does use the valid order above.
Because the whole star is contracted, the remaining tensor is symmetric
in the other kernel's surviving slots, so the derivative representation
controls the required full contraction norm.

The resulting `U_(a,i)` remains own-`a`-free. I independently checked the
Gaussian creation bound by normalized Hermite coefficients: an output
multiindex can receive at most its total degree many own-root creations.
The root-hit term has the exact diagonal matrix expression
`diag[Q(B circ K_(r-1))B]`; flatness gives only an averaged square bound,
which is exactly what the theorem uses.

The any-root covariance estimates also check: the star/tree case has
the middle factor `B circ Q^(circ(r-1))`, whose absolute row sums are
`O(n^-1/2)`, and the nonstar case has a vanishing global Frobenius
contraction matrix. The collision covariance errors require only entrywise
control and are paired against an already bounded star kernel.

The primary checks were made directly, not merely by recalling a CLT name.
[Nourdin--Rosinski, Theorem 3.4 and equations (3.5)--(3.6)](https://arxiv.org/pdf/1112.5070)
give block moment-independence from vanishing unsymmetrized mixed
contractions and bounded second moments, with covariance of squares
dominating each contraction norm. This is the needed statement even when
the other block is not Gaussian. [Noreddine--Nourdin, Theorem 1.1](https://arxiv.org/pdf/1009.1310)
allows the vector Gaussian approximation from covariance convergence and
component fourth-cumulant convergence without a nonsingular limiting
covariance assumption.

## 3. Input replacement and exact marked tests

The nonlinear forest kernels admit low-influence replacement after
collision deletion. The repeated-slot count is performed locally before
transport; its error is then propagated in averaged L2. A per-column
error is not summed without its weights.

The first-chaos quantity `Y_i=(B X_T)_i` is exactly square-free and occurs
linearly in its unmarked test. Thus a fourth coordinate derivative has
only `Y_i partial^4 f_i` and `4 partial Y_i partial^3 f_i`. The repeated
derivatives of the square-free low-influence local test have the asserted
scales. Total influence, followed by Cauchy--Schwarz over the replaced
coordinate, handles the second term without assuming small individual
influences of `Y_i`. Its averaged L2 bound suffices.

For the marked test I checked the exact Boolean two-spin identity after
removing the own spin from the unmarked channel. Each differentiated
factor excludes both relevant own spins. The off-diagonal `Q_ij` sum is
bounded by `sqrt(n)||Q_i||_2`, leaving `O(n^-1/2)` after the two derivative
factors. The diagonal `Q_ii=1` yields precisely `E[M h_T] E psi`.
No Gaussian replacement of the distinguished root spin is needed.

## 4. Zero variance, bounded responses, and the final means

The auxiliary Gaussian `delta xi` makes all channel variances bounded
below before polynomial approximation. Its density domination by the
largest-variance Gaussian is in the correct direction and has a constant
depending only on fixed `L,delta`. Finite polynomial separable tests span
the required jointly bounded tests while preserving their parities.

The tested regression after adding this noise has coefficient
`s v/(v+delta^2)`, not `s`; the proof retains this coefficient. Conditioning
back on the unsmoothed channel and using
`y tau_eta(y+e)>=|y|-eta-2|e|` removes the auxiliary noise uniformly,
even at `v=0`. This avoids the invalid uniform polynomial approximation
of `sign(sqrt(v)N)` as `v` tends to zero.

Finally, `mu_plus/minus=+/-F+H sign(BF)` are exactly feasible, and their
half energy difference is `sum H|BF|` with the correct factor. Jensen
in the root standard deviations, rather than their variances, uses the
proved Schur inequality in the correct direction. The resulting gain
does not depend on the fixed operator cap, so principal deletion can
be removed after the fixed-response matrix limit.

## 5. Verdict and scope

No unresolved step was found in the finite-response center theorem.
The response track also derived the marked polynomial identity independently
through Gaussian own-input integration by parts and whole-functional
hybrid replacement, in `resumed_response_marked_nonlinear_transport_2026_09_06.md`;
that alternative is not needed to close this audit.

This verifies a successful feasible marked/unmarked response mechanism.
It does not supply an iteration theorem for the enlarged response family,
a sharp constant, matching upper recovery, or convergence of the original
normalized minima.
