# Independent audit: adaptive spectral regularization at disorder-entropy scale

Date: 2026-09-07. Full reconstruction PASS for
`decisive_independent_entropy_compatible_spectral_regularization_2026_09_07.md`.

The real-cube core argument uses the simultaneous diagonal majorant with
trace <=4K_G Q(B). Deleting diagonal entries STRICTLY above the displayed
threshold gives fewer than epsilon n vertices, including the integer
boundary case. Selecting the first acceptable deletion set from a fixed
finite list is measurable; no measurable SDP selector is needed.

For a fixed mask, the cap and operator-norm refill tails have the stated
normalizations. The cap test uses at most 2^n projectively reduced signed
spin states, h<=epsilon n^2 independent uniform entries, and threshold
2sqrt(epsilon)n^1.5. The operator test uses sum(2x_i x_j)^2<=2, a 1/4-net
of size <=9^n, and the factor-two net comparison. Their failure probabilities
sum to less than 1/2 for n>=2. The good event depends only on the mask and
the refill, not on the retained matrix values. Conditioning therefore costs
at most log2 relative entropy, uniformly over all masks actually selected.

In the joint reference law the mask is independent uniform among 2^n sets.
The deterministic adaptive mask costs exactly n log2. Conditional on any
reference mask, original retained coordinates and newly sampled coordinates
are independent uniforms; consequently the OUTPUT reference remains the
full uniform coefficient cube. Data processing thus gives the claimed
output KL bound despite complete dependence of the actual selected core
on the input matrix. The unused original coordinates may remain in the
joint reference before applying the output map; no extra entropy is paid.

The deterministic output costs are q(B')<=q(B)+2sqrt(epsilon) and
delta(B')<=delta(B)+h/d<=delta(B)+4epsilon. Its normalized operator bound
is exactly 4K_G C/epsilon+8. The outer Gibbs mass with q>C is at most
exp(-alpha d), since F>=q and Jensen gives L_n<=C-1. Conditioning on the
complement has variational cost at most
-log(1-exp(-alpha d))/(alpha d).

Finally the restricted integral uses an UNNORMALIZED reference restriction.
Its variational formula is therefore correctly expressed with D(P||U_d),
not entropy against a normalized restricted law. Substituting the output
law proves the stated complete inequality. Fixed alpha,tau, then large n,
then epsilon decreasing to zero is a legitimate uniform approximation.

The result preserves the n^2-scale disorder entropy, not merely one good
matrix. It does not supply bounded-operator all-order recovery, continuity
of variance occupancy from another profile, or a microstate-volume limit.
