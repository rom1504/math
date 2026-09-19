# Twisted chiral campaign: final mathematical assessment

Bounded authorization: three ACTIVE hours, begun September18 at17:09:11UTC,
interrupted after approximately50 minutes, resumed September19 at18:58:29UTC,
with the closing endpoint21:08:29UTC. This file records the closing evidence;
it does not authorize an automatic extension. Original convergence is OPEN.
The previously reported asymptotic interval is unchanged and its older upper
proof chain is not claimed to have been reconstructed in full here.

## 1. Exact object and finite claims

Write Q(A)=max_x |x^T A x|/2 and beta(A)=max_(x,y)|x^T A y|. Let F(A)
be the minimum cap of the allowed order2n twisted double with child A.
For B=G^T A G, with G a signed permutation, the exact identity is

```math
Q\begin{pmatrix}A&B+\operatorname{diag}d\\B+\operatorname{diag}d&-A\end{pmatrix}
=\max_{I,z}\left\{2|z_I^TA_{I,I^c}z_{I^c}|
 +|2H_{B[I]}(z_I)-2H_{B[I^c]}(z_{I^c})+d(I)-d(I^c)|\right\}.
```

Set y=z on I and y=-z on its complement. This gives the two summands
before absolute values. Reversing z on one part changes only the cut term,
so maximizing its sign gives the displayed sum of absolute values. This
retains joint maximization; it does not claim cancellation after separate
channel bounds. Chirality centers the range but does not certify optimality.

Complete child/twist enumerations yield the following. Equivalence means
switching, permutation and global sign reversal of the child.

| child order | number of optimal-child classes | parent family minima |
| ---: | ---: | --- |
| 3 | 1 | 5 |
| 4 | 1 | 10 |
| 5 | 1 | 13 |
| 6 | 1 | 18 |
| 7 | 3 | 21,21,25 |
| 8 | 2 | 30,30 |
| 9 | 9 | five classes33; four classes37 |
| 10 | 2 | 44,44 |

The small cases at parent2 and4 give1 and4. Recorded global optima through14
are reproduced by favorable child selection. The entries30 at16 and33 at18
are upper witnesses, not new global optimality certificates.

At child10, TWO independent complete root-gauged censuses recover the same
725,760 minimizers in two classes. Both complete twist families have minimum44.
A nonoptimal conference child has Q=15 rather than13 and parent cap40.
Thus exact optimal-child inheritance fails at this finite order. This does
not falsify asymptotic near-minimizer transfer. In particular M20=40 is NOT
proved. The external claim of global optimality through20 remains unverified.

The conference parent40 has a standalone arithmetic proof using a symmetric
Hadamard completion of order20. An explicit spin attains40. The arithmetic
lemma reconstructs a case of the classical maximum-excess result, not a
claimed new external discovery. A second inequivalent cap40 parent has no
Hadamard sign-diagonal completion. Two bounded joint searches made4,139,221
and1,415,920 exact candidate evaluations respectively without finding38;
these heuristic failures give no global or family lower bound beyond the
separately completed exhaustive searches.

Sources: `twisted_chiral_computation_2026_09_18.md`,
`twisted_chiral_switching_table_2026_09_19.md`, and
`twisted_chiral_hadamard20_2026_09_19.md`. All are in this artifacts directory.

## 2. Strongest uniform constructive theorem: two-sided near-order transfer

For every full-sign A of order n with Q(A)<=K n^(3/2), and 1<=r<=n/4,
there is ONE principal restriction C of size n-r such that Q, beta and F
all change by at most

```math
300[(r/n)\beta(A)+n\sqrt r]+r=O_K(n\sqrt r).
```

The essential extra step is cycle-skipping repair of the chosen signed
permutation. This preserves the SAME-child twisted-family constraint; an
ordinary deletion of a parent need not do so. Random-row symmetrization,
factor-one scalar contraction and the image-containment proof were each
independently reconstructed, with exact finite checks.

Conversely, given ANY child and ANY chosen allowed parent D, one may insert
r child vertices, extend the twist by the identity, and retain both old
matrices as exact induced submatrices. With N=n+r, L=nr+binom(r,2), and
V=8nr+4binom(r,2), one completion simultaneously satisfies

```math
Q(A')\le Q(A)+\sqrt{2L(N+2)\log2},
\quad\beta(A')\le\beta(A)+\sqrt{8L(2N+1)\log2},
```

```math
Q(D')\le Q(D)+r+\sqrt{2V(2N+1)\log2}.
```

Only new CHILD edges are independently sampled; correlated parent entries
retain the orbit relation. The exact old/new coefficient matrix has TT^T=4I,
giving the variance8n per new vertex. A union bound leaves positive
probability for all three estimates at once. This is an actual extension,
not an assumed favorable bridge. Its concentration ingredient is classical.

Consequently T_n=min_A F(A) satisfies |T_(n+r)-T_n|<=C n sqrt(r), uniformly
for r<=n/4. Its normalized values have a relative-order1/2-Hölder modulus.
This removes the orbit-compatibility obstacle to filling relative-o(1) order
gaps. It does NOT construct a dense good sequence or cross a fixed ratio
with vanishing error. The n sqrt(r) discarded-row norm scale is unavoidable
for that proof method, not a lower bound on the actual cap difference.

Proofs: `twisted_chiral_near_order_2026_09_19.md` and Section16 of
`twisted_chiral_symmetry_followup_2026_09_19.md`. Independent reconstruction:
`twisted_chiral_uniform_final_audit_2026_09_19.md`.

## 3. Joint cancellation: a sharp theorem, but a leading realization gap

For hollow symmetric A, the coordinate-diamond resource

```math
\Delta(A)=\sup_{|a_i|+|b_i|\le1,\ |c_i|+|d_i|\le1}
 |H_A(a)-H_A(c)+b^TAd|
```

satisfies Delta(A)<=4beta(A)/3. The factor is sharp, proved by an exact
14-term rational certificate and a weighted six-type witness. Hollow
full-sign sharpness is at QUADRATIC cap scale, not the minimizing scale.
The known implementation uses balanced clones, whose child cap is at least
4Q(A)-n. These clones are outside the near-minimizer regime. Conference
examples and the inherited lower bound show a leading resource-to-native
gap. Therefore the certificate does not prove F(A)<=4beta(A)/3+o(n^(3/2)).
Its coarse block approximation has error O(n^2/sqrt(log n)), which is too big.

An exact native identity Q([[A,C],[C,-A]])=max_J beta(C_J) retains the whole
cut optimization and is not a strict reduction. Switch-only twists are
untwisted doubles of partial seed reversals; actual permutations strictly
improve several small optimal-child families. Neither observation supplies
the required uniform cancellation estimate.

Sources: `twisted_chiral_uniform_2026_09_18.md`,
`twisted_chiral_director_native_gap_2026_09_19.md`, and
`twisted_chiral_switch_only_reduction_2026_09_19.md`.

## 4. A tested new selection mechanism and its exact limits

At fixed parent spins and a conditional pair of twisted child spins, all
single-spin stability inequalities can be kept jointly. Averaging the
remaining sector permutations and matching signs gives an EXACT product
of two generating permanents. Brégman bounds give a coefficientwise,
polynomial-time conditional upper bound from two row-degree arrays. A
scaled-capacity refinement also retains column compatibility. These are
valid selection certificates, not a proved estimate of the remaining
exponential child-pair average.

The exact average is already insufficient on the smallest test cases: its
expected positive stable count above the true family minimum is5/4 at
parent6 and4/3 at parent8. It first certifies7 and12, respectively, just
as the corresponding raw projective first moment does. Independent full
group enumerations agree with the conditional permanent formula exactly.
The projective cutoff1 cannot be increased on symmetry grounds: actual
allowed parents can have only one violating positive stable projective state.

The larger wind tunnel also finds substantial filtering but no useful
sampled count below one. Samples are not population certificates. A
high-cap Hall example shows why row counts alone miss simultaneous matching
compatibility; a proved O(n log n) dominance-matching algorithm repairs
support feasibility after the fields are computed, but does not count matches.
Separately, explicit cap-(1/2)n^(3/2) Hadamard-derived children possess
high-energy conditional pairs for which stability probability is exactly1.
Thus even low cap does not give a uniform per-pair stability penalty. Their
frequency under the actual twist distribution is the unresolved issue.
Neither example rules out an averaged theorem for selectable near-minimizers.

Sources: `twisted_chiral_stable_permanent_2026_09_19.md`,
`twisted_chiral_stability_row_bound_audit_2026_09_19.md`, and
`twisted_chiral_stability_wind_tunnel_2026_09_19.md`.

## 5. Scalable failures, with scope retained

For a symmetric full sign matrix H, vec(H) is an explicit Boolean vector
giving vec(H)^T(H tensor H)vec(H)=tr(H^4). Thus

```math
Q(\operatorname{hollow}H^{\otimes2k})/q^{3k}
\ge\tfrac12[(\operatorname{tr}H^4/q^3)^k-((\operatorname{tr}H)^2/q^3)^k].
```

This shortens an archived even-tensor obstruction: Hadamard completions tend
to1/2, non-Hadamard fixed completions diverge. Heterogeneous P=product_tensor
(H_j tensor H_j) also has normalized cap at least1/2-o(1). For a chiral D,
the optimal diagonal Gram defect is exactly

```math
\Omega_*=\sum_{i\ne j}(D^2)_{ij}^2+2q(q-2)-4\sum_i|(D^3)_{ii}|.
```

Every sign diagonal completion has normalized even-tensor cap at least
(1/2)(1+Omega_*/q^3)^k. The two order20 cap40 witnesses have Omega_*=0 and960;
their best possible self-square caps are respectively exactly4,000 and at
least4,480 at order400. The latter is a lower witness, not an exact cap.
This excludes literal paired self-tensor amplification, NOT restrictions,
new twists or global changes of internal edges. See
`twisted_chiral_tensor_witness_2026_09_19.md` and its integer certificate.

The fractional-shell suggestion has a different precise obstruction:
at width4n, one- and two-spin flips of a maximizer already span the entire
edge-feature space. Fixed-delta shells of width delta n^(3/2) therefore
eventually leave no exact null direction. Narrow-shell refinements must
also control the previously unprotected states becoming new maximizers.
This is not an impossibility theorem for every fractional-rounding scheme.

## 6. Exact remaining convergence obligations and research judgment

A sufficient native step would select an asymptotically minimizing child
and an allowed twist with

```math
Q(D)\le2\sqrt2\,M_n+O(n^{3/2-\delta}),\qquad\delta>0,
```

uniformly in all sufficiently large child orders (or on an appropriate
closed seed class with the same liminf). No such theorem was proved.
Power saving gives a summable normalized error under iteration, but one
doubling subsequence alone does not prove convergence. A second independent
integer multiplier with summable errors, or uniform favorable proportional
thinning of the constructed descendants, would suffice. The proved
near-order modulus is weaker than that coverage requirement. Exact
countermodels to the doubling-only inference are preserved in the uniform
artifact; they are not counterexamples to convergence of M_n itself.

The campaign therefore delivers genuine uniform NEAR-order construction,
complete finite-family evidence, a sharp joint resource theorem and scoped
scalable failures. It does not improve the asymptotic interval or resolve
convergence. Neither finite cap40 nor the joint4/3 certificate is a new
asymptotic signing upper bound.

Director judgment: stop at the bounded endpoint; do not automatically extend
the same search. One concrete possible next campaign is a bounded attempt
to estimate the stability-permanent tail averaged over ACTUAL low-cap child
spin pairs, starting with the exact small failures and demanding a strict
subleading-error native selection theorem. That campaign is justified only
by a new way to control this average; more finite caps, row-only relaxations
or Hadamard tensor examples do not supply it. No current evidence supports
promising that a longer run of the same mechanisms will prove convergence.

Canonical proofs, exact certificates, searches, parameters, negative results
and independent reports are committed; scratch research is hash-preserved
in the September19 dated research archive. Verification verdicts never
replace the displayed proof dependencies.
