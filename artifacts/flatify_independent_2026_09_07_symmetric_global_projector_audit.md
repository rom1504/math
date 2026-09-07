# Independent audit: symmetric global projector localization

2026-09-07, final bounded audit. Full read of
`flatify_adversary_2026_09_07_symmetric_global_projector_localization.md`:
PASS. This verifies the operation and its stated scope, not an original
recurrence or convergence claim.

- One-gauge linear term: removing the outside sign s_j inside an absolute
  value leaves a genuine Rademacher sum with variance P_ii, including
  the u=j summand; no false independence from s_j is needed.
- Double term: M_uv=P_iu H_uv P_vj has Frobenius norm
  sqrt(P_ii P_jj). Symmetrizing and removing the diagonal cannot increase
  that Frobenius norm. The diagonal trace is bounded by the same number.
  The hollow symmetric N therefore has both Frobenius and operator
  norm at most v=sqrt(P_ii P_jj).
- For Z=s^TNs/2, cut Jensen, Gaussian linearization, and
  ||N[cross]||F^2<=||N||F^2/2 give log E exp(tZ)<=2t^2v^2 for
  |t|<=1/(sqrt(8)v). At |s^TNs|>12Lv, Chernoff has exponent
  -(3/sqrt(2))L+1/4<=-L for L>=1. Thus the 13L bound after restoring
  the trace is valid, with the displayed simultaneous probability.
- The relative majorant (1+sqrt(13LP_ii))(1+sqrt(13LP_jj)) is exact
  as an upper bound. The contraction and variance budgets are
  2n sqrt(13Lr) and V<=2nr+4sqrt(13)n sqrt(Lnr), respectively.
- Symmetric pair rounding gives variance at most2V, not V, because
  an off-diagonal random sign has bilinear coefficient x_i y_j+x_j y_i.
  Its centered summand is at most4. Therefore
  2sqrt(Vu)+(8/3)u is a valid uniform Bernstein threshold. Hollowing
  costs at most n/2 in each quadratic value.
- For any ordinary Hadamard H_h, W_(i,a),(j,b)=H_ib H_ja is symmetric,
  and its square is h^2 I by two independent Hadamard orthogonality
  sums. No symmetric Paley matrix at order h is assumed.
- The banked ordinary-Hadamard gap h-sqrt(n)=O(n^(21/80)) gives
  h^2-n=O(n^(61/80)). Multiplying sqrt(h^2)-sqrt(n) by n/2 yields
  exactly the stated O(n^(101/80)) leading-factor padding cost.

The recovered matrix is genuinely hollow symmetric and changes global
edges if needed. Its reference target remains a projected symmetric
Hadamard, not an unpaid weighted composition of optimizing children.
