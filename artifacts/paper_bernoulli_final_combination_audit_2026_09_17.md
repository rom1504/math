# Final composition and offset audits

2026-09-17, closing campaign audit by the Bernoulli track. These are
independent reconstructions of named canonical results, not additional
claims of original-problem convergence or external priority.

## 1. Quantified prepared-code transfer

The complete finite theorem and application in
[feature-parent comparison](paper_bernoulli_feature_parent_comparison_2026_09_17.md)
received full independent reconstruction by the director, localization
and discrepancy tracks. Section5 was separately checked by localization:

    preparation cost: n^(3/2-gamma/2),
    protected window: n^(3/2-gamma),
    code entropy: O(n^(1-gamma/4)log n),
    scalar feature error: O(n^(rho/2-1/4)+n^(rho-1/2)),
    parent error: O(n^(3/2)[n^(rho/2-1/4)
                               +n^(-gamma/8)sqrt(log n)]).

Here0<gamma<2/3,0<=rho<1/2 are fixed and q=O(n). At
(rho,gamma)=(1/3,1/2), the four relevant exponents are5/4,1,
17/12 and23/16, with the final selection term carrying sqrt(log n).
Both ensembles use the prepared old signing W. The preparation cap
inequality is not a pointwise bound on H_W-H_A. At fixed-density q,
the displayed selection certificate exceeds the protected window by
a factor n^(7gamma/8)sqrt(log n). This does not prove an actual-error
lower bound; it prevents using this certificate alone to remove the
old-code restriction.

## 2. Equal variance and affine offsets

The complete director
[offset tradeoff](paper_director_equal_variance_offset_tradeoff_2026_09_17.md)
was read and independently reconstructed, including its positive
complement. All constants pass.

For centered Z, D_Z(s)=E|s+Z|-|s| is nonnegative and

    integral D_Z=E Z^2,
    integral_(|s|>R) D_Z=E(|Z|-R)_+^2.

These follow by Tonelli from the two stop-loss formulas on the positive
and negative half-lines. No density or symmetry assumption is needed.
If X,Y have variance1 and common subGaussian proxy K, their difference
d is2-Lipschitz and integrates to zero. The condition d(0)<=-delta
forces negative area at least delta^2/2. Each tail integral is at most
4K exp(-R^2/(2K)), so R=sqrt(2K log(32K/delta^2)) leaves positive
area at least delta^2/4 inside[-R,R]. Hence some such offset has
d(s)>=delta^2/(8R). Atomic laws cause no endpoint problem.

For the positive complement let S=A U, with U uniform[-1,1],
A>=0 independent and E A finite. At fixed A=a, Gaussian heat
differentiation gives

    F'_a(v)=[2Phi(a/sqrt(v))-1]/(2a),
    F''_a(v)=-phi(a/sqrt(v))/(2v^(3/2)).

The a=0 limits and differentiation under the A-average are valid:
both derivatives are uniformly bounded on every positive compact
v-interval. For independent V in[l,u], E V=1, strong concavity yields

    F(1)-E F(V)
      >=E exp(-A^2/(2l))*Var(V)/(4sqrt(2pi)u^(3/2)).

The finite-integrability hypothesis is present in the canonical theorem.
Alternatively the same calculation can be formulated for the finite
excess E[|S+sqrt(v)G|-|S|] without E A finite.

The child interpretation is correctly conditional: the cavity offset
depends only on other columns and is independent of the replaced field,
but its uniform-scale-mixture distribution is an EXTRA hypothesis at
every hybrid stage. Neither ordinary symmetry nor Gaussian bridge
inputs imply it. The classical scalar inequality is not a full-parent
value theorem.

## 3. Final synthesis scope

The full [final synthesis](paper_portfolio_final_synthesis_2026_09_17.md)
was read. Its actual low-cap response theorem, growing-rank physical
realization, sharp block information/tail benchmark and remaining
value/escape obligations agree with the independently audited sources.
The one wording correction sent to the director was to specify EXACT
ISOTROPY in the spectral-band mixture floor: the nonpaired mixture
extension uses isotropy to make the mean component query variance one.
A single cold component need not obey that floor. This qualification
is explicit in the canonical spectral-band proof.

No improved original extremal constant, cross-order defect,
convergence or nonconvergence conclusion has been supplied.

## 4. Updated feature dependencies

The updated
[growing-rank realization](paper_director_growing_rank_variance_realization_2026_09_17.md)
Appendices B and C were then read in full. Both pass independent
reconstruction. In Appendix B the c_K-dependent Gaussian-sign repair
uses Frobenius control of ENTRYWISE sine, not an invalid operator
contraction. It retains the exact prescribed K-proxy. The eigenvalue
threshold makes every adaptive retained frame uniformly diffuse, and
the finite-dimensional minimax set preserves the entropy and MGF
budgets. Its endpoint information coefficient is
(1-1/K)log(1/theta)+O_K(1), with theta fixed before the order limit.

In Appendix C, symmetry gives the shifted Fourier difference with
factor cos(ts); the hot shifted Stein test remains1-Lipschitz; and
subtracting |s| before the small physical repair removes dependence
on offset size. Thus the all-offset rate is genuinely uniform. It
does not allow selecting the offset after observing the same column,
and it does not itself compare a joint maximum over query words.
