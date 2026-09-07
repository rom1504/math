# Unfinished: weak energy reward does not yet imply iid bridge universality

2026-09-07. Historical attempted route, followed by the now-proved resolution in Section 2. The new quenched universality theorem closes the transfer gap in the regime p/n tending to infinity; it does not settle the bounded-p/n regime or rare favorable selection.

For low-cap children and p=||A||op||D||op, the Gaussian bridge law R=I-rho A tensor D/p has rank-one-index process covariance differing from iid by at most beta(A)beta(D)/p=O(n³/p), uniformly over two Boolean pair indices. At the natural bridge variance n² this is o(n²) if p/n tends to infinity. Standard Gaussian comparison could therefore transfer the Gaussian bridge maximum up to o(n^(3/2)) in that regime.

This observation does not establish the corresponding result for the ACTUAL SIGN bridge. The exact-covariance small-source MGF theorem bounds every deterministic source, including increments, but does not by itself transfer a supremum or a Gibbs source chosen after sampling. Treating it as full process universality would be an unjustified step.

An attempted coupling is also insufficient with current bounds. Couple G=R^(1/2)z to iid z, and add small independent Gaussian noise to both components to make the joint covariance nonsingular. The new MGF theorem then applies to their sign difference, but its covariance contains a quadratic coupling term controlled by

    (x^T A² x)(y^T D² y)/p².

The cap only supplies x^T A²x<=n beta(A)=O(n^(5/2)), from ||Ax||infinity<=n and ||Ax||1<=beta(A). Thus the available uniform bound for this product is O(n^5/p²), which need not be o(n²) in the allowed p=O(n^(3/2)) range. This coupling does not pay the supremum error.

No actual-law iid obstruction, Gaussian universality theorem, or original recurrence is claimed here. A sharper optimizer-specific response bound or a genuine correlated-process comparison would be additional work.

## 2. Resolution using the new quenched comparison

The canonical full law obstruction, independently read and checked, is `flatify_independent_2026_09_07_opposite_spectral_nearmin_law_obstruction.md`. It additionally supplies explicit planted-clique additive near-minimizers and the paid moving gap epsilon_n=n^(-1/7). The derivation below records this agent's independent reconstruction and allows unchanged deterministic parent-energy offsets throughout the Gaussian comparison.

The subsequently proved theorem `flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md`, independently audited in `flatify_adversary_2026_09_07_quenched_gaussian_sign_audit.md`, supplies exactly the missing quenched transfer. Here is the resulting scoped actual-law obstruction.

Assume Q(A),Q(D)<=C n^(3/2), p=||A||op||D||op, and p/n tends to infinity. Fix 0<=rho<1 and sample the ACTUAL sign bridge C=sign G with covariance of G equal to I-rho A tensor D/p. Its covariance-matched Gaussian bridge Y has vector covariance

    I-kappa A tensor D, kappa=(2/pi)arcsin(rho/p)=O(1/p).

Let J be an iid standard Gaussian bridge. For two Boolean pair indices, including an optional common absolute-polarity bit, the covariance difference of their bridge energies has magnitude at most

    delta_n=kappa beta(A)beta(D)=O(n³/p).

This uses |x^T A x'|<=beta(A)<=4Q(A), and similarly for D. The number of pair indices is at most 2^(2n+1). A finite Gaussian softmax interpolation bounds the expected maximum difference by O(sqrt(delta_n log(2^(2n+1))))=O(n²/sqrt(p)). One can verify this without an equal-variance premise: the derivative of the expected log softmax is bounded by the inverse-temperature squared times delta_n, and the softmax/maximum error is at most log(number of indices) divided by that inverse temperature. Optimizing gives the displayed bound. Arbitrary identical deterministic child-energy offsets are allowed in this comparison.

The new quenched sign theorem therefore gives

    E Q(parent with C)=E Q(parent with J)+o(n^(3/2)).

The iid Gaussian parent's cap is at least beta(J), by independent reversal of each whole child. The previously independently audited heat-martingale/SK comparison gives

    E beta(J)>=[2 c_heat-o(1)]n^(3/2),
    c_heat=Gamma(3/4)/(sqrt(pi)Gamma(5/4))
          =0.7627597635... .

Thus the actual opposite-spectral sign parent has the same lower floor in expectation. On the full parent order 2n, its normalized coefficient is at least c_heat/sqrt(2)=0.5393526012..., above the current favorable target.

### Typicality is paid

This expectation statement also holds with probability tending to one, at any fixed smaller leading coefficient. Conditional on W in G=W+sqrt(epsilon)z with epsilon=1-rho, the entries are independent signs. Changing one bridge entry changes the parent cap by at most two, so conditional variance is O(d), d=n², by bounded differences.

Let F(W) be the conditional expected cap. Its partial derivative in W_i is f'(W_i)/2 times the difference between the expectations with that sign fixed to plus or minus one. That difference has magnitude at most two. Hence |partial_i F|<=|f'|<=C epsilon^(-1/2), uniformly and regardless of the deterministic children. Gaussian Poincare, with Cov(W) operator norm at most two, gives Var F(W)=O_epsilon(d). Total variance is therefore O_epsilon(n²)=o(n³), and Chebyshev proves the claimed concentration at the n^(3/2) scale.

This rules out a typical successful output from this particular fixed-rho random law in the regime p/n tending to infinity. It does NOT rule out an exponentially rare favorable draw, a selection/conditioning mechanism, bounded-p/n children, or unrestricted changes of the old internal edges.
