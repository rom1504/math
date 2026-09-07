# Independent audit: balanced prescribed-partition obstruction

Date: 2026-09-07. Verdict: PASS for
`decisive_independent_balanced_power_nearmin_obstruction_2026_09_07.md`.

The reverse Minkowski calculation is correct: for r=2/3,
`u -> (u^(1/r)+d)^r` is convex, so its two-point Jensen inequality
gives aggregate energy gain at least `2^(1/r)d=2 sqrt(2)d`.

The simultaneous diagonal majorant has trace at most `4 K_G Q(B)`.
Applied separately on each N/2 shore, threshold `16 K_G C sqrt(N)`
leaves at least N/4 usable indices. Restriction to a k-subset has cap
at most half the sum of its diagonal entries, hence the stated
`8 K_G C k sqrt(N)` bound. This does not require either shore to
intersect a preselected global spectral core.

For the bounded dyadic sequence `u_j=m_(2^j)^(2/3)`, the positive
increment has liminf zero; otherwise u_j would eventually increase
by a fixed positive amount. Along this subsequence the two child
minimum bounds give exactly `a^r+b^r >= P^r-e_j`, without missing the
factor from their N/2 orders.

The chosen `delta=max(1/log N,sqrt(e_j))` tends to zero and satisfies
`e_j<=delta^2`. Integer clique order
`k=floor(sqrt(2 delta) N^(3/4))` gives normalized clique energy
`d=delta+o(delta)`. Each old block costs at most
`eta=O(sqrt(delta) N^(-1/4))=o(delta)`. Opposite child orientations
cause no difficulty: the parent upper estimate is a triangle
inequality, while each child's lower estimate tests its own old
oriented ground state.

All cap variables stay in a fixed positive compact interval, so the
power errors are Lipschitz. Consequently the powered defect is
`N[c0 delta-O(e_j+eta)]`, with the positive constant stated in the
proof, and eventually exceeds `c N/log N`.

Crucial scope: only the STARTING parent is exact. The MODIFIED parent
is an original normalized near-minimizer with excess
`O(delta N^(3/2))`; that excess is not negligible at the much smaller
scale of a reverse-Fekete error. The theorem says nothing against a
selected exact parent, an existential partition, or a minimum-value
inequality. It supplies a balanced pointwise obstruction only.
