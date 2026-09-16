# Polynomial Boolean BH: focused transfer audit

2026-09-16. Final synthesis after independent audits; authorized session
05:13:45--07:13:45 UTC. This is a bounded paper investigation, not a
restart of the autonomous original-problem loop.

## Verdict for the original problem

No convergence, nonconvergence, or improvement of the original numerical
interval was proved. The paper is mathematically relevant and its proof
withstands independent reconstruction, but its polynomial degree loss does
not remove the signed high-moment obligation. Several proposed alternative
uses were tested beyond that observation and have precise scoped barriers.

The reported original frontier remains

```math
.4333221116640807\le\liminf_n\frac{M_n}{n^{3/2}}
\le\limsup_n\frac{M_n}{n^{3/2}}
\le U_0-\zeta<.493608094,
```

with the same positive tiny zeta and dependencies as the
[September 7 synthesis](principle_director_final_synthesis_2026_09_07.md).
This campaign reconstructed the dependencies needed for its new claims,
not the entire earlier lower/upper proof chain. Stored finite cap witnesses
through n=14 were replayed; their global optimality is not a proof input here.

## 1. Primary theorem and a positive analytic deduction

The complete primary source is Ivanisvili,
[*Polynomial growth of Bohnenblust--Hille constants on the Hamming cube*](https://arxiv.org/html/2609.12427v1).
For degree at most m, its coefficient exponent is q_m=2m/(m+1); its
dimension-free BH bound is O(m^27), and its weighted square-function bound
is O(m^22). Coordinate restrictions, mixed norms, complex scalars,
noise damping, the fixed-dimensional contraction supremum, and the
full-level-flat influence hypothesis were independently checked.

A checked refinement of the SAME proof gives

```math
\left(\sum_{r=1}^m r^{-5}
 \|\widehat f^{=r}\|_{2r/(r+1)}^2\right)^{1/2}
 \le C m^3\|f\|_\infty,
\qquad
\|\widehat f\|_{q_m}\le C' m^{11/2}\|f\|_\infty.
```

The contraction is at most (27/32)^(1/8)<1, with exact finite checks and
an analytic all-degree tail. The improved low-level extraction uses a
self-contained Joukowski/Cauchy estimate. This is a secondary analytic
deduction, not a claim of external novelty or of progress on M_n.
See the [primary reconstruction and refinement](bh_2026_09_16_independent_paper_audit.md)
and the [fresh, separately initialized reconstruction](bh_2026_09_16_fresh_refinement_check.md).
The latter also tightens the harmless closure multiplier by using the
exact sum of squared low/high levels; neither exponent nor M_n changes.

## 2. Powers: the exact remaining coefficient problem

Let c be the reduced Walsh coefficients of H_A^k, let
L=||H_A||_(2k), and F=||c||_(4k/(2k+1))^(1/k). If pi is the normalized
squared coefficient law and alpha=2k/(2k+1), then EXACTLY

```math
\log(F/L)=\frac{H_\alpha(\pi)}{4k^2}.
```

Thus at k proportional to n the coefficient norm and the ordinary
moment differ by a factor 1+O(1/n). For every K=o(n), one can choose a
single signing at each order whose entire collection of power certificates
2<=k<=K is o(n^(3/2)), even if the power is selected adaptively afterward.
At superlinear k the moment uniformly recovers the original cap
multiplicatively; additive recovery uses a bounded normalized-cap class.

A stronger ACTUAL-CLASS estimate also holds. If Q(A)<=C n^(3/2), then
uniformly over every K>=2 and every such A,

```math
\frac{\max_{2\le k\le K}
 \max\{F_k(A),W_s(H_A^k)^{1/k}\}}{n^{3/2}}
\le C_C\left(n^{-1/4}+\sqrt{K/n}+K/n^{3/4}\right).
```

It remains true with the actual reduced degree in the coefficient norm.
Thus all exact minimizers have vanishing adaptive power certificates for
K=o(n^(3/4)). The planted full-sign family shows that this exponent cannot
be increased for the whole bounded-cap class; it is not minimizer sharpness.

The top coefficient is k! times a signed hafnian, but lower boundaries
contain genuine signed Eulerian cancellation. Keeping only favorable
diagrams is invalid: the exact order-14 regression exhibits a huge excess
from a selected positive diagram class over the complete moment.

Both the older exp(O(sqrt(m log m))) BH bound and the new polynomial bound
have root loss 1+o(1) in the linear-degree regime. Even an elementary
dimension-dependent coefficient/L2 comparison is already enough there.
Both old and new root losses are also summable on geometric degree scales,
with vanishing total overhead as the starting degree grows; the new bound
does not repair a hidden geometric-accumulation problem in this comparison.
The missing LINEAR-degree optimized signed-moment estimate did not get easier.
See [powers and flattening](bh_2026_09_16_power_analysis.md).

## 3. Strongest target-class filter theorem

Put Y_A=H_A/n^(3/2). For every fixed C,c>0 there is eta>0 such that,
uniformly for ALL full signings with Q(A)<=C n^(3/2) and all scalar
polynomials P of growing degree d<=eta n^(3/4), bounded by one on [-c,c],

```math
\frac{W_s(P(Y_A))}{d^a}\longrightarrow0\quad(a>1),
\qquad
\frac{\|\widehat{P(Y_A)}\|_{4d/(2d+1)}}{d^b}
 \longrightarrow0\quad(b>1/2).
```

Here s>=0 is fixed and W_s^2=sum_(r>=1)r^(-2s)||fhat^{=r}||_(q_r)^2.
The same statements hold for mean-square uniform frozen-coordinate tests
on free sets I with |I|>=2d, chosen independently of the frozen signs;
random partitions into such eligible sets are included. Frozen coefficient
sums are retained exactly. This does not cover maximization over frozen
assignments, adaptive selection I(y), or every smaller-fibre degree test.

This includes actual exact minimizers without spectral-flatness assumptions:
the elementary random-sign bound already supplies C=1. Dependencies are
the safe cap-to-operator estimate ||A||op<=sqrt(8Q(A)), a reconstructed
quadratic tail estimate, Bernstein--Walsh, exact Walsh degree via odd
hafnians, and a low-degree Taylor calculation. Under the stronger
||A||op=O(sqrt(n)) hypothesis the barrier covers EVERY growing sublinear d.

In particular, neither the paper's exponents nor the refined exponents
can distinguish a false fixed cap by these tests in the stated window.
This is a certificate theorem, not a statement that all minimizing energy
laws coincide. The endpoint scale is sharp for the larger bounded-cap class:
a planted clique of n^(3/4) vertices in a Sylvester signing produces
constant-energy tail mass exp[-O(n^(3/4))], detectable by Chebyshev/BH
at degree Theta(n^(3/4)). These examples are NOT certified minimizers.
See [filter theorem](bh_2026_09_16_nonlinear_filter_audit.md) and the
[independent audit](bh_2026_09_16_ports_filters_independent_audit.md).

## 4. Other mechanisms and exact scope of their failure

| Tested mechanism | Verified outcome | What remains unexcluded |
| --- | --- | --- |
| Unrestricted degree-two BH constant | A six-spin cap-five witness limits this direct asymptotic certificate to 5/30^(3/4)<.391 | A genuinely stronger dense/optimizer-specific inequality |
| Scalar signed permutation flattening | A fully flattened output level retains at most eta=||c||_1^2/(D||c||_2^2) of the input variance; a stored square has eta=27/35 | New coefficient dispersion structure or noncontractive transformations with paid norm cost |
| Hilbert permutation orbit | Coefficient lengths flatten, but scalar BH/influence fails even on full row-balanced sign orbits; actual Hilbert cap is Theta(n), not Theta(n^(3/2)) | A different scalar mechanism preserving the actual extremum |
| Generic joint linear/quadratic ellipse | Without a stronger zero-field input, no such norm-only ellipse bootstrap improves a seed w>=.4088703384 | Field-law-sensitive or optimizer-specific joint inequalities |
| Weighted reciprocal ports | Valid Schur inequality; in the unmasked repaired p=31/32 model, every core either gives only o(N^(3/2)) gain or belongs to a sector already bounded by (.016387114+o(1))N^(3/2) | Boolean-specific joint bounds, additional masks, or other packet mechanisms |
| Multi-moment Chebyshev filter | The archive ALREADY proves O((n/d)^2) cap recovery; L2 detection is stronger than the BH-conditioned test | A theorem evaluating the required optimized moment data |

The [mechanism analysis](bh_2026_09_16_mechanism_analysis.md) and
[Bohr/ellipse audit](bh_2026_09_16_bohr_ellipse_independent_audit.md)
give exact hypotheses. Generic positive quadrature supplies a matching
moment-oracle lower example, not an actual quadratic-signing example;
see [Chebyshev comparison](bh_2026_09_16_chebyshev_moment_recovery.md).
The [planted-tail theorem](bh_2026_09_16_planted_tail_sharpness.md)
and [independent Chebyshev/sharpness audit](bh_2026_09_16_chebyshev_independent_audit.md)
make the filter barrier's precise bounded-cap scope explicit.
The port dichotomy uses full nonself sign columns and the archived uniform
operator-controlled repair hypothesis. It is not a Boolean eigenvector or
energy-attainment claim; the small fixed-cap sector need not be o(N^(3/2)).

## 5. A genuine secondary application: Boolean Bohr radii

Combining polynomial BH with the checked real-Boolean homogeneous Bohr
theorem gives, for BOTH arbitrary coefficients and the complete flat-sign
class minimized over degree,

```math
R_n,\ R_n^{\rm flat}
=\sqrt{\frac{\log n}{n}}
 \left(1+O\left(\frac{\log\log n}{\log n}\right)\right).
```

Every exact minimizing degree is log n+O(sqrt(log n log log n)).
This improves the rate provided by the earlier subexponential proof.
However, removing all fixed degrees, including degree two, leaves that
minimum exactly unchanged eventually. The aggregate radius therefore
does not recover the original fixed-quadratic minimum through this identity.
Real Boolean versus complex torus norms and flat versus arbitrary
coefficients have exact finite separations. See the
[primary-source mapping](bh_2026_09_16_bohr_flat_separation.md).

## 6. Exact original gap and research decision

For k=ceil(alpha n), define
mu_(n,alpha)=min_A ||H_A/n^(3/2)||_(2k). Uniformly in A,

```math
2^{-(n-1)/(2k)}\frac{M_n}{n^{3/2}}
 \le\mu_{n,\alpha}\le\frac{M_n}{n^{3/2}}.
```

Existence of lim_n mu_(n,alpha_j) along any sequence of FIXED
alpha_j tending to infinity would imply convergence, by sending j to
infinity LAST. Neither
uniform convergence in alpha nor a common value of those limits is needed.
The limsup-minus-liminf gap is at most
(2^(1/(2alpha_j))-1) lim_n mu_(n,alpha_j), which vanishes by the
uniform cap bound. This remains an old
obligation, not a strict reduction furnished by the paper. The stronger
archived Chebyshev moment proxy changes its recovery rate but not its
unproved optimization/realization step.

**Director judgment:** no longer BH-to-signing campaign is justified by the
evidence obtained here. Reopen only with a concrete new signed-coefficient
estimate or a joint mechanism that preserves the actual scalar extremum
and supplies a quantitative original-value implication. The analytic BH
refinement and Boolean Bohr corollary may warrant separate mathematical
communication after further external novelty checking; that is a different
objective. Do not resume the general research loop automatically.

All finite diagnostic norms/quadratures are distinguished from exact
integer/rational certificates and analytic proofs. Proofs, counterexamples,
primary downloads, failed calculations and replay outputs are preserved;
preservation itself is not verification.

Publication: substantive checkpoint `3bdaf9e` and closing proof/archive
commit `430aa44` were pushed normally. The closing snapshot has1,888
hash-verified payloads; the post-publication audit matched all1,863
remaining research working files. See the
[preservation and dependency audit](../research_archive/publication_audit_bh_2026_09_16.md).
The only intentionally untracked ordinary working originals are a reviewed
compiled build and older raw research already preserved losslessly.
