# The stronger physical-response slope: finite failures and asymptotic scope

2026-09-17. Bernoulli-track independent test of the proposed stronger
target beta(C_0(A))/sqrt(n)<=3c(A)/2, where c(A)=Q(A)/n^(3/2).
Equivalently beta(C_0(A))<=3Q(A)/(2n). This is NOT implied by the
new fixed-discount low-cap theorem.

## 1. Exact finite comparison on the complete ground codes

The [full physical-column games](paper_discrepancy_full_column_games_2026_09_17.md)
include every Boolean column and every absolute ground word. Their
stored matrices' global-minimum labels are imported archive provenance;
this audit independently checks the caps and complete rational primal
and dual inequalities, not a new global classification.

| Order/case | Q | unrestricted beta | isotropic beta | 3Q/(2n) |
|---|---:|---:|---:|---:|
| 6 | 5 | 5/3 | 5/3 | 5/4 |
| 8, class 0 | 10 | 3/2 | 7/4 | 15/8 |
| 8, class 1 | 10 | 3/2 | 3/2 | 15/8 |
| 10 | 13 | 53/25 | 144/65 | 39/20 |
| 12 | 18 | 9/5 | 9/5 | 9/4 |
| 14 | 21 | 98/39 | 98/39 | 9/4 |

Thus the literal all-order target fails already at orders 6,10,14,
even WITHOUT requiring isotropy. At order14 its exact ratio is
392/351>1, with normalized response about .67157953 versus target
.60133779. At order12 the ratio is instead exactly4/5: the unrestricted
isotropic law improves on the best maximal-radial subclass, whose
response is exactly9/4. A radial-only test would miss this gain.

The order12 gap is a GROUND-only fact. At deficit window2, its
unrestricted value is452/185 and isotropic value1323/536, both above
9/4; the latter target ratio is147/134. Order14's window2 value is
4837/1768 in both games, also larger than its ground value.

## 2. Transparent exact witnesses at orders12 and14

At order12, the uniform law on the20 projective ground words is
exactly isotropic. Every absolute Gram row has13 zeros,6 entries4,
and one entry12. Its response is therefore9/5 at every ground word.
Enumeration against ALL 2^11 physical projective columns gives
minimum9/5, so this law is simultaneously a primal and a dual.
Its signed covariance is an involution but is not proportional to A.

At order14, the uniform law on all156 projective ground words is
exactly isotropic. Every absolute Gram row has values
0,2,4,6,10,14 with respective multiplicities57,49,21,21,7,1.
The row sum392 gives98/39. Enumeration against ALL 2^13 physical
projective columns has the same minimum, attained precisely at the
ground columns. This is again a matching primal/dual certificate,
not merely a sampled-column numerical optimum.

These identities are preserved with full exact geometry in
[the discrepancy ground-code analysis](paper_discrepancy_actual_ground_dual_geometry_2026_09_17.md).

## 3. What these failures do and do not settle

They refute a finite universal statement with no error term. They
DO NOT refute

    beta(C_0(A_n))/sqrt(n)<=3c(A_n)/2+o(1)

for every sequence of exact minimizers as n tends to infinity, or a
weaker selectable-family statement. The order14 normalized cap is
about .40089186, below the certified asymptotic lower constant, so
it cannot simply be read as an asymptotic example at that cap value.

No operation supplied here turns these fixed finite response games
into arbitrarily large exact minimizing full signings while preserving
both the normalized cap and the response-game ratio. Coherent blowup
changes the cap scaling; a tensor or Hadamard completion changes the
complete ground code and its physical-column game. Declaring one of
those transformations harmless would be an unproved essential step.

The other asymptotic stress tests in this campaign do not fill that
gap either. Polynomial-size arbitrary codes can have response kappa
but need not be actual high-energy sets. The full-sign hard high-energy
construction is explicitly separated from the ground by a macroscopic
energy gap. The whole block-rank-one dictionaries have a sharp response
transition but cannot all be near-ground words of a positive-cap full
signing. None is an asymptotic exact-minimizer falsifier.

## 4. The available principled positive bound is weaker in quantifiers

The [selectable deletion baseline](paper_bernoulli_selectable_deletion_baseline_2026_09_17.md)
gives a subsequence of asymptotically minimizing principal cores B_n
and ONE physical column h_n with the all-query inequality

    |h_n.x|<= (3c_inf/2+o(1))sqrt(n)
                 +Q(B_n)-|H_(B_n)(x)|.

Consequently their ground-code unrestricted game satisfies the desired
leading slope along that selectable subsequence. But those cores need
not be exact minimizers at their own order, the one-column law is not
isotropic, and repeating it across many new vertices multiplies the
deficit allowance. It does not establish the all-order or reusable-law
claim currently under investigation.

The new low-cap theorem instead supplies one exactly isotropic,
(3/2)-subGaussian law for a complete macroscopic high-energy code,
but only with a fixed discount from kappa. Its certified2^(-67)
discount is far from the required3c/2 slope. These results address
different missing pieces and should not be conflated.

## 5. Solver-free independent replay

`computations/paper_bernoulli_2026_09_17_ground_slope_audit.py`
reads the frozen rational certificates, recomputes each actual cap and
complete code, verifies every primal query inequality and every dual
physical-column inequality, and checks isotropy with exact integers.
It uses no LP solver. It passes42,576 all-column lower checks and
1,716 all-query upper checks, for all stored ground cases and the
order12/order14 deficit2 cases. Python compilation passes.

Output: `tmp/paper_portfolio_2026_09_17/bernoulli/ground_slope_audit.json`.
The independent replay verifies the finite response games, not a new
proof of the imported global-minimality labels or an asymptotic claim.
