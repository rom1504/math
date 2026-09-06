# Independent audit of bounded-response nuclear dual splitting

2026-09-06. PASS under the explicit hypotheses of
`transfer_adversary_bounded_response_nuclear_dual_split_2026_09_06.md`,
which was read in full. This is not an audit of a rich full-return cap
theorem: the local query BC remains a separate obligation.

## 1. The two different local comparisons are sufficient

For the sparse dual part the tested field is the pair (X_i,X_j), both
of fixed original degree P. The literal R_i has no degree-P part.
Consequently lower-degree mixed contractions are proper on X, and full
contractions into a higher R degree are proper on the RIGHT R tensor.
The listed small cuts make the pair jointly Gaussian independently of
the ACTUAL (W_i,R_i). A singular two-by-two covariance is harmless for
the characteristic-function argument. Fixed-degree moment bounds license
the polynomial h_k(X_j) factor. This proves the required uniform entry
formula without assuming polynomial density in the law of R_i.

For the dense dual part the field is Y_i=sum_j Msmall_ij U_j. Its
degree is p=kP, so the only non-proper contraction against R is now its
degree-p covariance. The row-Euclidean alias hypothesis and the row l2
bound of Msmall control precisely that covariance. Higher R degrees
have small RIGHT cuts, while lower ones use the Y proper cut. This is
a different local comparison from the sparse part; neither alone would
prove the nuclear theorem.

## 2. Norms and order of limits

For ||M||op<=1 and fixed tau>0, the large-entry part has row and column
l1 norms at most 1/tau. Therefore its total l1 norm is at most n/tau,
its operator norm is at most 1/tau, and
`||Msmall||op<=1+1/tau`. The mask does not preserve operator norm.
The sparse error is at most n eta_n/tau.

The exact covariance comparison of U is used only after tau is fixed.
Its dense variance error divided by n is at most
`(1+1/tau)^2 ||Cov(U)-K^{circ k}||_*/n`, tending to zero. The Schur
main has bounded operator norm and Msmall has Frobenius norm at most
sqrt(n), so its averaged variance bound is independent of tau.

After a fixed row variance cap, stable Boolean Stein supplies a modulus
going to zero as tau goes to zero. The response is bounded and depends
on the entire retained local query list, whose JOINT law is preserved.
The discarded roots cost O(R_0^-1/2), from averaged second moments.
Taking n first, then tau, then R_0, uniformly over the original dual
matrix proves the nuclear statement. No unproved quantitative covariance
rate is required in this fixed-tau version.

The dense target matrix costs O(tau n), because k>=2 and
`sum_ij |K_ij|^k <= const tr(K^2)=O(n)`. This completes the dual test;
the entrywise target is not being discarded without a norm estimate.

## 3. The same-source alias bound and uniform diagram cuts

For a k-branch probe of original degree P greater than every old
primitive degree M, pulling the one higher old-source transport back
gives k disjoint Hall matches. Each is a proper or retained-label
influence merge. The source entry is O(n^-k/2); its one flat transport
loses at most sqrt(n). Thus the alias entry is O(n^((1-k)/2)), and the
row squared sum is O(n^(2-k)), tending to zero for k>=3. This does not
apply to a later source built out of the probe itself.

For the intended fixed injective odd-tree and forest diagrams, the
UNIFORM global-cut premise can be checked directly by the banked
boundary-graph theorem, rather than inferred from generic polylogarithmic
circuit estimates. Declare the output root an extra boundary port.
At every nonroot vertex the degree parity agrees with the number of
seed marks; internal even-mark vertices have even degree, free marked
vertices odd degree. Thus a bridge-forest leaf not containing a free
mark would have odd exiting degree but only even internal vertex degrees,
which is impossible. The root port covers its own exceptional leaf.
All edge matrices have fixed bounded operator norm, and finite equality
patterns preserve these properties. Hence every selected diagram has
every global cut bounded by a fixed constant.

This checks the premise only for those literal finite tree/forest
diagrams. An arbitrary bounded-function projection or arbitrary later
history is not such a diagram by definition. Its uniform global cuts
and alias property would require their own proof.

## 4. Honest endpoint of the result

The theorem makes bounded response tests possible without an exponential
moment of a high-degree literal coherent variable. It therefore removes
the mixed-nuclear bottleneck in the stated finite source setting.
It does not make BC-beta E independent of the probe. That separate
full-return query is exactly what is needed before inferring an actual
rich-history spin-flip gain or improving the original lower bound.
