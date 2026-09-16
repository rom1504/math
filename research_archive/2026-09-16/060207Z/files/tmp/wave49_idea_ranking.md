# Wave 49 independent idea ranking

Evidence used after generation: Wave 48 `Updated frontier` (§10.101.4) and
`STEERING.md` at the Wave 48 cutoff. This scratch ranking makes no assertion
that a proposed lemma is true.

1. **Project-row ground-lift Pareto theorem.** Seek a maximal or near-parent
   selector `S` and child ground `y` satisfying simultaneously
   `||A[S]y||_2^2=O(n^(9/4-c))` and
   `q_n-Q(A[S])=O(n^(3/4-c))`. By (10.1223)--(10.1224), this gives a
   project-row completion and hence `D_C>0`. Test exact minimizers for the
   best joint Pareto value and derive it from block replacement or averaged
   local margins. This attacks the newly isolated mass obligation directly.

2. **Cubic high-harmonic coarea inequality.** Expand the non-strict boundary
   target exactly. With
   `Xi_z=sum_(j>=2)(lambda_1-lambda_j)W_j(z)`, it is equivalent to
   `E[1_C a_z Xi_z] <= delta E[1_C a_z^3]`. Seek a minimizer-specific proof
   from port regrets/full completion profiles, or a scalable family violating
   it. This is sharper and more testable than asking vaguely for negative
   degree--escape correlation.

3. **Full-slack cross-block Bellman pressure.** Combine the projective slack
   identity (10.1225) over a weighted overlap cover and prove an all-state
   margin of the form `s_[n](d)+4 sum_(e in F)M_(d,e)>0`. Any averaging must
   charge the positive slack layers that migrate in `A8`. A useful theorem
   would quantitatively dominate migration by selector-cover pressure; an
   exact jointly realizable equality is its falsifier.

4. **Profile-conditioned arbitrary-ground agreement.** Partition selectors
   by principal cap, child-ground local-margin profile, and external Gram;
   within one affordable cell seek a common or decodable project-row full cut
   whose event satisfies (10.795). The construction must use the retained
   deficit/tolerance and must not impose `L_S>=q_n` or another complement
   threshold. This is the most faithful direct fallback but currently has no
   structural theorem forcing a large cell.

5. **Partial-completion vector balancing under minimizer flatness.** Given a
   child ground `y`, minimize `||A[:,S]y+A[:,T]w||_2` over outside signs `w`.
   Prove a signing-specific partial discrepancy bound at the project scale,
   perhaps from a quantitative near-conference estimate for exact minimizers.
   The all-positive signing is the law-free obstruction; a positive result
   must explicitly use minimizer flatness not presently known.

6. **Threshold-integrated row-supported coarea.** Apply layer cake jointly to
   `a_z(h)`, `B_z(h)`, and row truncation, seeking a threshold where both
   `D_C(h)>0` and the non-strict cubic inequality hold. An exact integral
   identity that makes the signed excess nonpositive would prove the whole
   leading package. Generic coarea and uniform-center versions are already
   false, so the center weighting must remain incidence- and row-sensitive.

7. **Finite exact coarea phase diagram.** Exhaust all exact minimizers through
   order eight, the order-nine example, and order-ten samples over thresholds
   and core scales. Record whether any project-row class has `D_C>0` and
   normalized boundary at most one, plus the first failing harmonic level.
   This is diagnostic rather than asymptotic, but a universal finite failure
   would rapidly redirect the leading route.

8. **Slack-aware fractional-to-Boolean rounding.** Use the lower bound
   `gamma_frac(S)>=(q_n-2K)_+` near the diagonal and round a fractional block
   flip while reserving margin for states of slack `t`. A concrete target is
   an additive integrality gap strictly below `gamma_frac` once all states
   with `s<=gamma_frac` are coupled through Bellman profiles. The full-edge
   fractional wall and A8 slack migration are mandatory tests.

9. **Dictator classification of favorable selector families.** Since
   non-strict coarea plus FKN forces a favorable family close to a coordinate
   dictator, derive the exact signing consequences of such a family: one
   vertex must be almost always included or omitted when a fixed row-good
   center is near-optimal. Seek a replacement contradiction or a direct
   recurrence from that vertex structure. This could convert the inverse
   theorem from an endpoint into a usable proof mechanism.

10. **Higher-order exact-minimizer obstruction search.** Use MILP/SAT to seek
    an order-ten or order-eleven failure of tight decomposition, a positive
    full-slack cross-block integrality gap, or absence of project-row ground
    lifts. Exact samples can falsify strengthenings but cannot by themselves
    prove convergence, so this remains a supporting diagnostic.

Selected attacks: combine 1, 2, 6, and 7 into one coherent project-row
coarea package; assign 4 as the direct bare-tail attack; and assign 3 with 8
as the full-slack cross-block attack. Ideas 5 and 9 are focused successors if
the leading package exposes the corresponding structure.
