# Marked feedback energy: exact source estimates and the remaining star

Date: 2026-09-06. Status: working proof boundary, not a completed energy
theorem. Local separation is proved separately in
`continued_feedback_marked_local_noise_separation_2026_09_06.md`.

Subsequent resolution: the complete energy proof now appears in
`continued_feedback_first_marked_history_energy_projection_2026_09_06.md`.
The source estimate in Section 2 remains a proof module used there.
Section 4 below records the precise historical gap, which was resolved
by all-global cuts, two-factor flat transport, and a separate source-level
split/no-split argument; it was not bypassed by the local law alone.

## 1. Frozen intended energy statement

Use `G=BS`, `D=S h2(G)`, `Y=BD`, and the literal returns `QS,QD`.
For bounded odd regular `f(G,Y)`, write its first local Gaussian
projection as `b0 G+b1 Y`, and put

```math
Z=B[f(G,Y)-b_0G-b_1Y],\qquad V=b_0QS+b_1QD,
\qquad C=H(G,Y)\psi(V+Z).
```

The desired energy projection has the same shape as the earlier scalar
theorem: retain `c0=H E_N psi(V+sigma N)` literally, retain its exact
cross with a linearly tested residual, and identify any remaining noise
self-energy. This note does not yet assert that the coefficient of the
linear residual can be replaced by its deterministic mean.

## 2. Summable covariance of exact old nonlinear source forests

Let `P,P'` be fixed finite linear combinations of old local Wick
forests, each with odd local branch count at least three. Work with
their EXACT squarefree source kernels before the next root transport.
Then

```math
\max_i\sum_j|\mathbb E P_iP'_j|+
\max_j\sum_i|\mathbb E P_iP'_j|=O(1).                       (1)
```

This strengthens the previous normalized-nuclear covariance statement
for the forest main only. It does not assert the same operator-level
bound for every raw bounded response.

First consider unrestricted Gaussian product kernels. The complete
branch-pairing classification in
`resumed_bound_audit_full_nonlinear_covariance_trace_2026_09_06.md`
gives (1) directly, as follows.

- Whole-branch pairings have at least three bounded-op Gram factors.
  Two such factors already have bounded absolute row and column sums
  after entrywise multiplication, by Euclidean Cauchy--Schwarz.
- A non-star component or two nontrivial star components give an
  entrywise `O(1/n)` term. Its absolute row and column sums are bounded.
- A sole nontrivial star has an entrywise `O(n^{-1/2})` bound and at
  least TWO additional whole-branch factors, since both local branch
  counts are at least three. The latter factors have bounded absolute
  row and column sums after entrywise multiplication. Thus this case
  even has absolute row and column sums `O(n^{-1/2})`.

There are finitely many fixed-degree pairings. Crosses of different
original input degrees vanish exactly.

Now let `K,K'` be two such global forest root maps, and let `Pi` be
the common marked-label distinctness projection. Both root maps have
bounded operator norm. The deleted row Hilbert norms are `O(n^{-1/2})`
by the established local forest collision bound. The covariance
correction is exactly

```math
E=K(I-\Pi)(K')^{\mathsf T}.
```

Every row of `E` has Euclidean norm at most
`||K_i(I−Pi)|| ||K'(I−Pi)||op=O(n^{-1/2})`, and similarly for columns.
Cauchy--Schwarz in the row or column index then bounds its absolute
sums by `O(1)`. Equivalently, for a self-covariance the correction is
positive semidefinite with diagonal `O(1/n)` and operator norm `O(1)`.
This proves (1) with exact squarefree sources.

Consequently, for any transported residual-source channel `Z'=BP'`,

```math
\|\mathbb E[P(Z')^{\mathsf T}]\|_{\max}
 =\|\operatorname{Cov}(P,P')B\|_{\max}=O(n^{-1/2}).          (2)
```

All root maps involved remain bounded. Equation (2) is an actual
entrywise-small return factor for fixed exact source forests.

## 3. What (2) settles in a star expansion

Expand a transported noise branch into its outer flat root edge and
its at least three old `G` or `Y` branches. If a transported-noise leaf
consumes a group of WHOLE old branches at the center, the group's local
branch count is odd.

If that count is at least three, (2) supplies an entrywise
`O(n^{-1/2})` factor. If the count is one, the only equal-original-degree
possibility in this canonical scope is an old `Y` branch paired with
the transported `h3(G)` residual. Its covariance has operator norm
`O(n^{-1/2})` by the exact cross formula in the local-separation note.

A whole coherent linear leaf gives a bounded-op linear covariance
factor. A whole degree-three marked leaf can match a whole old `Y`
branch, giving `B` followed by its bounded-op coherent transport, or
three old `G` branches, giving the explicitly small `h3(G)--D` cross.
These cases reproduce the small noise-leaf factors needed by the
earlier scalar star proof.

## 4. The exact unresolved configuration

The previous paragraph does not cover an old `Y` branch whose three
marked slots split among several coherent marked `D` branches, while
those `D` branches also meet other old branches or a cross-root
coherent coefficient. The old `Y` proper cut is small, but a marked
`D` proper cut need not be: its own-spin slot can have influence of
order one. A single small cut alone does not establish a negligible
energy trace.

In particular, before replacing a random local noise-linear coefficient
`A_i(W_i)` by `a_i=E A_i`, one must establish

```math
\frac1n\mathbb E\sum_{i,j}B_{ij}
  [A_i(W_i)-a_i]Z_i c_j^0\longrightarrow0.                  (3)
```

Local independence at the SAME root does not imply (3). Nor do the
polylogarithmic coherent Walsh root maps alone imply it. The canonical
identity `D_l=e_l` times two flat `B_l` row factors suggests an
additional own-slot factor in the troublesome split networks, but
that factor has not yet been completely classified.

The computation
`continued_feedback_marked_old_tree_projection_2026_09_06.py`
now tests (3) directly using `A=cos(QD)`, `c0=sin(QD)`, and
`Z=B[sin(Y)−exp(−1/2)Y]`, as well as testing the full candidate energy.
It includes actual hollow Steiner signings and their positive twin
lifts. Finite agreement or shrinking discrepancies do not settle (3).
