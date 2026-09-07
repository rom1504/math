# Independent audit: sparse stability of balanced bulk

2026-09-07. **PASS**, Sections2--4a of
`principle_synthesis_2026_09_07_balanced_bulk_sparse_stability.md`.

The single-spin change of the ferromagnetic ratio completion is exactly
-2 sigma x f -4b sqrt(N) x s/q+4b sqrt(N)/q. Thus the minority threshold
includes the positive physical 2b/q correction. At minority fraction at
most1/4 its absolute value is at least b.

For any common light threshold A rho, condition on the heavy edge signs.
The absolute aggregate conditional mean on the minority coordinates is
bounded by L sqrt(epsilon/p0)|R| <= b|R|/4. Jensen therefore gives the
convex hinge expectation at least b|R|/2. The hinge Lipschitz constant is
at most L sqrt(A rho |R|). The previously proved product-convex zero-event
inequality gives exp[-c q r_i/rho], c=b^2/(64L^2 A), uniformly in all
positive allowed minority fractions. No central-limit approximation is
used; zero Lipschitz makes the zero event empty.

The heavy-tail classification is deterministic given the word and frame
array, independent of outer edge signs. Although the per-row proof
conditions on its own heavy signs, it gives an unconditional event bound,
so graph Finner applies legitimately to the shared independent edge
variables with exponent1/2.

The active-fibre global-max count is correct: reset constant fibres to
positive, leaving both bulk energy and squared magnetization unchanged.
This preserves global maximality and removes the spurious 2^m charge in
the sparse-active count. Good rows cost exp(q h(r0)+O(log q)); bad rows
cost 2^q. The stated strict entropy inequality pays their union.

The stronger minority-mass version needs no such reset. Take rho=rbar.
If mass-good rows carry at least1-delta of total minority mass, their
Finner penalty is exp[-c(1-delta)N/2]. All words with rbar<=r0 are counted
by 2^m sum_{s<=r0N} binom(N,s), hence exp[N h(r0)+o(N)]. For sufficiently
small fixed r0 this proves exclusion of the stated nonconstant LOCAL
maxima with failure exp(-c'N).

The premise on incoming heavy energy is substantive. This result neither
excludes the heavy-tail/intermediate-bias branch nor supplies an all-bias
cap theorem. The actual fair edge signs must retain their conditional
independence after choosing frames and masks.

The later entropy-adaptive Section4b also passes. At exact total
minority count s and r=s/N, its threshold
eta=(4/c)[h(r)+log2/q+2logN/N] gives a fixed-count union probability
at most exp[-N h(r)-m log2-4logN] <=2^(-m)N^(-4).
Summing at most N counts and both polarities gives2^(1-m)N^(-3).
Hence the mass-good fraction tends to zero along nearconstant local
maximizers with r tending to zero. The threshold is asserted only when
eta<=1; the incoming-tail condition remains essential.
