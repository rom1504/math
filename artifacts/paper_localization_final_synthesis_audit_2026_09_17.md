# Localization final synthesis audit

2026-09-17, closing campaign audit. Scope: full independent read of
`paper_portfolio_final_synthesis_2026_09_17.md`, with rank, information,
uniformity, full-parent and original-problem implications checked against
the frozen canonical proofs. Also independently read the complete
quantitative parent corollary in
`paper_bernoulli_feature_parent_comparison_2026_09_17.md`, Section5.

## 1. Mathematical conclusions: PASS

The growing-rank physical feature statement correctly uses fixed variances
and leverage constant, r=o(sqrt(n)), exact isotropy/full support, and
separately proved every-direction subGaussian control. The scalar estimate
is uniform in ALL Boolean queries and ALL deterministic real offsets.
The conservative stated error
(r^2/n)^(1/4)+r/sqrt(n) is valid; the later sharper Frobenius repair bound
only improves its second term. The information bound
-(r/2)sum pi_j log(v_j)+o(r) is valid and uniform in the rank budgets
used by the covariance-capture minimax theorem.

The fixed-K AppendixB and scalar-offset AppendixC of the combined
realization theorem were independently read in full: PASS. In particular
the Gaussian repair uses Frobenius entrywise-sine control, not an invalid
operator-norm sine contraction; its constant can be chosen so that the
common proxy is EXACTLY any prescribed K>1. Fixed accuracy precedes the
order limit. No simultaneous shrinking-accuracy rate is asserted.

The parent comparison retains the entire fixed child and all new spins,
but only on a code declared before the new columns. The Gaussian
comparator is the actual VARIANCE MIXTURE, not an ordinary covariance-
matched Gaussian. Conditional scalar replacement compares each signed
index's mean, and the two one-sided expected-maximum intervals pay ONE
old-code fluctuation term. No joint-query coupling is inferred from the
pointwise scalar theorem.

The block-code information theorem's final upper construction is finite
for every r when b>=4096 epsilon^(-2). The shared-phase law is exactly
isotropic and full-support; its information and subGaussian costs are
both paid. The fixed-K floor and joint information/tail frontier have
the correct order of limits: fixed K and accuracy, then physical block
dimensions, then accuracy tending to zero. The coefficient1-1/K in the
frontier is independently reconstructed, including the loglog correction
in the lower estimate. No O(1) two-sided additive frontier error is claimed.

The stated operator-free low-cap response theorem has the correct
uniformity over all full signings and every word in the indicated
macroscopic energy code. Its explicit2^(-67) discount was checked against
the separate rational certificate, not confused with this track's looser
sixth-power formula. The exact rational replay passed again in this audit.
The inherited half-range theorem and scalar comparison remain its explicit
dependencies; a practical n0 is not supplied.

## 2. Synthesis wording corrections communicated during this audit

1. The physical scalar error is O(e_n), so the displayed parent inequality
   needs C_(L,v)q sqrt(n)e_n, not an unsupported coefficient one. Corrected
   in the refreshed synthesis.
2. The comparison of the old-code selection scale with its protected
   window is stated for q PROPORTIONAL to n. It is not generally true
   for tiny q. Corrected in the refreshed synthesis.
3. The convex-query information lower bound uses an UNNORMALIZED average
   reduction Delta from the independent mean. The normalized discount
   costs order n information on bounded-covariance query laws. Corrected
   in the refreshed synthesis.
4. The smooth-energy likelihood falsifier requires cap-bounded actual
   energies H_A/n AND a normalization bounded below. Both are essential;
   the final synthesis now includes them.
5. The .786393873897 floor requires EXACT ISOTROPY of the bounded-band
   Gaussian-sign mixture, not just a band condition on arbitrary
   correlated components. Corrected in the refreshed synthesis. Its
   formulation for all exactly isotropic mixtures (not only paired ones)
   is valid: the Schur/arcsine operator bound places every normalized
   component query variance in[2/3,4/3]; exact mixture isotropy makes
   their mean one. The endpoint chord lower bound for the concave square
   root gives kappa[sqrt(2/3)+sqrt(4/3)]/2, with the uniform scalar-CLT
   error averaged over components. Removing isotropy would make the
   sentence false.

The added equal-variance affine-offset tradeoff was also read in full:
PASS. Equal variance makes the integrated response difference zero;
the negative triangle has area at least delta^2/2, the paid subGaussian
tail area is at most delta^2/4, and a positive shift cost delta^2/(8R)
follows. This does not assert that an actual child creates that shift.

The final positive complement for independent uniform-scale-mixture
offsets was also reconstructed: conditional on A=a, differentiating
the averaged absolute Gaussian response gives
F_a''(v)=-phi(a/sqrt(v))/(2v^(3/2)), including its a=0 limit.
The stated curvature lower bound on [l,u] and strong-Jensen gain
follow. This is a conditional cavity-law criterion, not an assertion
about actual optimizing children. The synthesis now explicitly takes
both block size and block rank to infinity for the sharp fixed-K floor.

## 3. Quantitative restricted-parent corollary: PASS

For fixed0<gamma<2/3 and0<=rho<1/2, cloned preparation gives

    cap cost O(n^(3/2-gamma/2)),
    window T=n^(3/2-gamma),
    log|C|=O(n^(1-gamma/4)log n).

Feature rank r=O(n^rho) gives scalar error
O(n^(rho/2-1/4)+n^(rho-1/2)). The latter summand is smaller. For q=O(n),
the restricted full-child comparison is therefore

    O(n^(3/2)[n^(rho/2-1/4)+n^(-gamma/8)sqrt(log n)]).

At (rho,gamma)=(1/3,1/2), exact rational arithmetic gives respectively
5/4 for preparation,1 for the window,17/12 for scalar parent replacement,
and23/16 for the selection contribution. The total bound is
O(n^(23/16)sqrt(log n)). BOTH ensembles use the prepared old signing W;
the cap bound for W does not provide a pointwise replacement of H_A by
H_W in an arbitrary parent. Favorable Gaussian-mixture value, useful
feature capture, and escape outside C remain unproved.

## 4. Final original-problem implication and replay boundary

The synthesis correctly makes NO new extremal-constant, convergence,
nonconvergence, or favorable cross-order-defect claim. The full inherited
strict upper chain is not represented as independently recertified in
this campaign. The next-step recommendation asks for a full parent VALUE
with all-energy escape paid, not merely a smaller scalar-response constant.

Final localization reruns all PASS: quartic isotropic laws, cold Fourier
feature identities, finite block-code information, low-cap constants,
Gaussian-edge stability, actual-sign noise stability, and nonlocal sign
response. Their JSON outputs remain in
`tmp/paper_portfolio_2026_09_17/localization/`. They validate finite identities
and replay arithmetic; no analytic asymptotic theorem is inferred from
floating-point tests. Canonical proofs retain their independent audit
records and the original-versus-new-versus-unproved scope distinctions.
