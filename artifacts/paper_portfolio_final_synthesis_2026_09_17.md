# Six-hour paper-combination campaign: final mathematical synthesis

2026-09-17,16:42:21--22:42:21UTC. Completed; paused for assessment. The user authorized
this bounded campaign, not an automatic extension. Original convergence
and nonconvergence both remain unproved. No original extremal constant
improves in this campaign. Results below mean proved here with independent
reconstruction, not a claim of external priority or formal verification.

## 1. Reconstructed tools and what survived transfer

| Primary mechanism | Reconstructed input/output | Actual transfer and limit |
|---|---|---|
| [Stochastic localization](https://arxiv.org/abs/2109.00709), with [Eldan's decomposition](https://arxiv.org/abs/1811.11530) | Gaussian observation of a finite law; conditional covariance, information and integrated precision estimates; support remains physical | Conditional precision can be combined with rounding, but information and mean response alone do not bound the full parent maximum. |
| [Gram--Schmidt walk](https://arxiv.org/abs/1708.01079), sharpened biased form in [Harshaw et al.](https://arxiv.org/abs/1911.03071) | Prescribed vector means and norm constraints give physical sign rounding with a subGaussian proxy; conditional precision costs add | The exact sharp constant is usable. It does not make an arbitrary dense residual lower order. |
| [Bernoulli-process decomposition](https://arxiv.org/abs/1305.4292) | An index set is decomposed into an l1-controlled part and a Gaussian-chain-controlled part, with universal comparison constants | It compares physical rounded processes with Bernoulli widths. Preparation additionally uses read-two information and Sauer counting; stability additionally uses mixed-tail chaining. Constants are not silently replaced by one. |

Additional primary arguments were reconstructed only as needed:
[Krivine/Grothendieck rounding](https://arxiv.org/abs/1711.10595),
[mixed-tail chaining](https://arxiv.org/abs/1309.3522),
[convex-distance concentration](https://arxiv.org/abs/math/9406212),
and [Brascamp--Lieb covariance bounds for log-concave measures](https://arxiv.org/abs/1106.0709).
Canonical chronological track files retain the decisive calculations,
failed transfers, references, and independent audit records:
[localization](paper_localization_2026_09_17.md),
[rounding](paper_discrepancy_2026_09_17.md),
[Bernoulli](paper_bernoulli_2026_09_17.md).

## 2. Strongest theorem about actual signings

For every sufficiently large hollow full signing A with
Q(A)<=.5 n^(3/2), there is ONE centered, exactly isotropic physical sign
law nu, with linear subGaussian proxy3/2, such that

    sup_(|H_A(x)|>=.30 n^(3/2)) E_nu|h.x|/sqrt(n)
       <= sqrt(2/pi)-2^(-67).

The law protects the entire indicated energy set simultaneously. There
is no operator-norm, low-rank, query-count or query-covariance assumption.
It need not have an efficient sampler or a short description.

The proof combines a signed-covariance alternative, Grothendieck
localization to a small COORDINATE set, and the half-range budget paid
by its full-sign complement. It is not ordinary same-map rounding.
The prior universal half-range proof was reread and its rational
certificate replayed. Three independent reconstructions check the new
argument. Its fully explicit discount is tiny; n0 is not made effective.
[Proof and dependencies](paper_director_low_cap_uniform_response_2026_09_17.md).

This removes a real structural hypothesis. It neither reaches the roughly
3c/2 parent slope nor controls maximization over all old and new spins.
Thus it is not a convergence theorem or a better bound on M_n.

## 3. Main general combination: non-Gaussian sign realization

Let U have orthonormal columns, rank r=o(sqrt(n)), and leverage
max_i (UU^T)_ii<=Lr/n. Fix finitely many positive variances v_j and probabilities pi_j
with sum pi_j v_j=1 BEFORE taking the order limit. There are centered,
full-support, EXACTLY isotropic physical sign laws nu with fixed linear
subGaussian control and relative entropy O(r), such that uniformly over
ALL Boolean x and ALL real offsets s,

    E_nu|s+h.x/sqrt(n)|
      =sum_j pi_j E|s+sqrt(1+(v_j-1)t_x) G|+O(e_n),
    t_x=||U^T x/sqrt(n)||^2,
    e_n=(r^2/n)^(1/4)+r/sqrt(n),       G~N(0,1).

The information bound is at most -(r/2)sum pi_j log(v_j)+o(r).
The theorem combines physical hot/cold quadratic tilts, positive convex
auxiliary measures, operator covariance control, a separate every-direction
MGF estimate, and a paid exact covariance repair. No covariance-to-tail
inference is used. The cold result is zero-field; it does NOT import an
all-external-field covariance theorem or efficient Glauber mixing.
[Complete statement](paper_director_growing_rank_variance_realization_2026_09_17.md).

The uniform OFFSET is consequential. For any fixed child signing D,
q independent bridge columns, and any old-spin code C chosen before
the bridge, replace the physical columns by the actual Gaussian variance
mixture with covariance I+(V-1)UU^T conditional on V. Then

    |E Q_C(physical)-E Q_C(mixture)|
      <=C_(L,v) q sqrt(n) e_n+sqrt(4Kqn log(2|C|)).

ALL new spins and both absolute polarities remain optimized. A single
new field enters its fully optimized branch as a shifted absolute value,
which proves the replacement without a union bound over new spins.
For q=O(n), log|C|=o(n), the error is o(n^(3/2)).
[Finite comparison](paper_bernoulli_feature_parent_comparison_2026_09_17.md).

Combining this with same-order physical cloned-block preparation supplies
subexponential nearcodes for any prescribed subleading window, at a
subleading cap cost. It DOES NOT show their favorable feature capture,
bound the Gaussian-mixture parent value, or exclude far old words.
For q proportional to n, the certified old-code fluctuation exceeds the
protected window in the available power-scale construction. This is a limitation of the bound,
not a proved fluctuation lower bound.

## 4. Sharp independent benchmark: information AND tail budgets

For all block-constant Boolean queries on r blocks, physical response
at most epsilon costs

    r log(1/epsilon)-O(r) <= minimum KL
                          <=r log(1/epsilon)+O(r).

Here 0<epsilon<1/4 and block size b>=4096 epsilon^(-2) suffice for
the finite two-sided statement at every rank. With fixed linear subGaussian
proxy K>1, as block size b and rank r both grow the optimal response floor
is kappa/sqrt(K), kappa=sqrt(2/pi). Write
J_K(delta) for the minimum limiting information PER BLOCK needed to
approach this floor within delta. Then

    J_K(delta)=(1-1/K+o(1))log(1/delta).

The physical order limits precede delta down to zero. The proof converts
near-equality in a norm/tail inequality into forced small-ball mass, pays
the binomial entropy cost, and attains it with exactly isotropic common-phase
hot/cold sign laws. This is more than covariance matching or a rate-distortion
definition. All three researchers independently reconstructed the argument.
[Sharp theorem](paper_director_block_response_information_tail_frontier_2026_09_17.md).

Replacing the common phase by independent phases preserves every block
marginal and exact covariance, but restores response kappa. The shared
label may have bounded entropy while the induced total correlation is
extensive. This is an explicit compositional distinction, not a model of
original minimizing ground-code geometry.
[Exact counterexample](paper_bernoulli_shared_phase_block_composition_2026_09_17.md).

## 5. Falsifiers, information limits and formalism changes

- EXACTLY ISOTROPIC Gaussian-sign mixtures whose latent spectra stay in[1/2,3/2] have
  response floor .786393873897..., above the useful slope below .75.
  Non-Gaussian laws escape that class. This is not a universal impossibility.
- For query covariance norm<=L, any UNNORMALIZED average physical
  response discount Delta from the independent mean costs
  KL>=Delta^2/(64L)-2. Thus ALL o(n)-degree positive Walsh densities
  fail to give fixed normalized discounts on bounded-covariance duals.
  Such a dual is not proved for every actual minimizing code.
- Uniformly bounded-Lipschitz likelihoods of finitely many actual energies
  H_A/n, with Q(A)<=C n^(3/2) and normalizer bounded below by a fixed
  positive constant, leave scalar response at kappa+O(n^(-1/4)). Extensive tilts,
  growing descriptions and rare conditioning are outside this statement.
- Independent sign-noise stability is proved after cloned preparation,
  but its certified error exceeds the mean contraction gain.
- Equal variance forces an offset tradeoff: a normalized unshifted
  discount delta under common subGaussian proxy K causes a shifted
  absolute-response loss at least delta^2/(8R) for some |s|<=R,
  R=sqrt(2K log(32K/delta^2)). This follows from the exact integrated
  stop-loss identity; it does not assert that actual children create
  that damaging offset. [Proof](paper_director_equal_variance_offset_tradeoff_2026_09_17.md).

The final check also supplies a positive complement. If an independent
offset has law S=A U, with A>=0 independent of U uniform[-1,1]
and E A finite, then
F(v)=E|S+sqrt(v)G| is concave in v. On[l,u] its negative curvature is
at least E exp(-A^2/(2l))/(2 sqrt(2pi) u^(3/2)). A mean-one variance
mixture therefore improves this AVERAGED shifted response by at least
that constant times Var(V)/2. This gives an explicit child-field shape
criterion under every hybrid replacement. It is independently verified;
actual optimizing-child cavity fields are NOT proved to have that shape.
Indeed one actual sign edge with unequal Gaussian field variances already
has a nonconcave SIGNED-branch value. Its full absolute cap is concave,
so the counterexample is not promoted to an absolute-parent obstruction.
A separate final proof does falsify unrestricted FULL absolute variance
concavity: for EVERY order-four signing, including exact minimizers,
independent field variances(v,10000,10000,10000) give F''(1/4)>.3024.
An independent proof and exhaustive rational enumeration verify the claim.
It remains an unequal-field finite example, not an isotropic bridge or
an asymptotic minimizing-sequence obstruction.
It can moreover be embedded in TWO centered exactly isotropic Gaussian
mixture field laws with common proxy10000: split one variance only in a
rare hot phase and compensate covariance by an unchanged cold phase.
The full absolute value increases. This rules out universal monotonicity
from global isotropy alone, but the fields are not physical sign columns
or unconditionally independent coordinates. The embedding was independently
reconstructed; it is still not an asymptotic obstruction.

The useful formal change is concrete: retain the joint variance phase,
its tail budget and information cost, and compare affine-offset responses
when a child is to remain optimized. A covariance matrix or unshifted
mean alone is insufficient. No new theory name is needed.

The positive and negative theorems are consistent quantitatively. If a
query law has covariance norm<=L, its top-r mass is at most rL/n;
bounded-rank capture cannot protect it with a fixed discount. The general
entropy inequality then charges a linear information budget for such a
discount. Low-rank physical realization solves the correlated-query regime,
not this extensive-information regime. The actual low-cap response theorem
is not claimed to have a subextensive-information representation.

## 6. Verification, original consequence and next decision

The original reported interval is unchanged:

    .4333221116640807 <= liminf M_n/n^(3/2)
       <=limsup M_n/n^(3/2) <=U0-zeta <.493608094.

The lower dependency used here was reaudited; the entire inherited strict
upper chain was preserved, not recertified from scratch this campaign.
Exact replay checks all42 physical-column primal/dual games on stored
orders4--14; their global signing-optimality labels remain imported.
Analytic theorems are not inferred from floating-point quadrature/MGF grids.
All promising combined claims received independent full-proof reconstruction.

Original gaps: no theorem gives useful low-dimensional capture of actual
minimizing high-energy sets; no favorable full parent value is proved;
and no all-energy estimate excludes escaping old states. Consequently no
sublinear full cross-order defect, convergence or nonconvergence follows.

Strongest justified next campaign: test the full parent VALUE for one
non-Gaussian physical law selected using actual low-cap energy geometry,
with its all-energy escape cost paid from the outset. Require a power-saving
defect or a scoped scalable falsifier; do not spend a campaign optimizing
only the tiny response constant. The present paper method merits further
targeted use, not automatic continuation or a claim that convergence is near.

The reconstructed Gaussian transforms, entropy duality and convexity tools
are classical. Focused primary-literature checks, including the distinct
[single-negative-outlier Ising theorem](https://arxiv.org/abs/2512.22803),
did not establish external novelty of our extensions; none is claimed.
Research preservation and replay provenance are in the
[publication audit](../research_archive/publication_audit_papers_2026_09_17.md).
