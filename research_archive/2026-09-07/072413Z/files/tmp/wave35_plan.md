# Wave 35 plan

Evidence read: README, Wave 34 Updated frontier (§10.87.5), and STEERING.md
through Wave 34. The only user-stated research objective remains convergence
of `M_n/n^(3/2)`; the ideas below are research judgments, not directives.

## Ranked ideas

1. **Rare-center free-energy optimizer.** Penalize both the uniform soft
   degree and row square,
   `Psi(z)=-log E_(U_m) exp(-lambda a_z(S))+gamma R_2(z)`, choose a global
   cube minimizer, and derive its exact one-bit Euler inequalities. The
   distance increment is supported on selectors containing the bit, while the
   row-square increment is explicit. Try to sum these inequalities into the
   strict soft margin leading to (10.967), or construct an actual-minimizer
   obstruction. This targets a maximum center rather than a circular center
   average.
2. **Projective list-agreement from two-selector replacement.** Treat the
   favorable child labels on each uniform `m`-set as a local projective list.
   Derive a two-block global-replacement inequality forcing two random lists
   to have approximately consistent restrictions on their overlap. Then use
   a complete-complex low-soundness agreement theorem, or an independently
   proved majority decoder, to obtain one global label on an
   `e^{-O(TL_0)}` selector subfamily. Explicitly retain orientation, list
   size, Hamming error, and the decoded cut's row square.
3. **Harmonic posterior-information charge.** Rewrite every endpoint cost
   `C_j` as a conditional KL and seek an exact chain rule charging
   `sum_j C_j` and adverse context migration to selector-posterior information
   or entropy production. Test whether marginalizing the global orientation
   before interpolation removes its unanchored edge without increasing the
   parent entropy. Prove (10.983), or give a finite exact counterexample to
   the proposed charge.
4. **Direct arbitrary-cut free energy.** Optimize over cuts the uniform
   lower-tail functional of `widehat ell(S,d)` with an `R_2(d)` penalty.
   Derive flip and block-exchange conditions for an optimizer, avoiding both
   center distance and adversarial selector laws, and try to prove (10.795)
   directly.
5. **Soft multi-selector block replacement.** Put several independently
   uniform selectors in one penalized response game. Ask whether shared
   outside edge responses preserve a project-scale penalized value even though
   the one-selector zero-slack game has the `A_6/A_9` walls. The target is a
   quantitative partial collision, not favorable mass for every selector.
6. **Johnson lower-tail plateau.** For a fixed candidate cut or center,
   compute exact transposition increments and Dirichlet energy of the selector
   deficit on the Johnson slice. Combine a minimizer-specific variance bound
   with reverse hypercontractivity or small-set expansion to turn a genuinely
   large seed family—not one Johnson ball—into the required subexponential
   lower tail.
7. **Random-chain stability of principal grounds.** Along a uniformly random
   deletion permutation, couple exact child ground lists at neighboring
   sizes. Bound the total number of large projective ground-state changes by
   a telescoping global-replacement budget, then decode a center stable for
   many `m`-subsets.
8. **Free-coordinate terminal flow in expectation.** Average the exact
   terminal min-cut (10.923) over uniform selector batches and use shared free
   coordinates to bound conflict paths without paying Johnson entry. Seek a
   direct `O(k_0)` total-flow event with probability `e^{-O(rL_0)}`.
9. **Quadratic susceptibility at an optimized cut.** Apply the coarsest-row
   closure identity (10.914) not at a parent ground but at the cut optimizing
   idea 4. Test whether its second variation supplies the missing nonlinear
   selector tail while the first moment remains circular.
10. **Pressure diagnostic under randomized deletion.** Express the best
    one-step pressure reward in terms of the same posterior information used
    in idea 3. Look only for a checkable square-root mechanism or an
    unbounded-family falsifier; do not promote constant shortfall without it.

## Selected independent attacks

- Agent A: idea 1, exact rare-center Euler inequalities, asymptotic scaling,
  and small exact-minimizer falsification tests.
- Agent B: idea 2, a quantified projective/list agreement interface plus the
  precise two-selector replacement lemma and row-control input it needs.
- Agent C: idea 3, endpoint conditional-KL chain rules, orientation
  marginalization, and migration control or counterexample.
- Main agent: idea 4, the direct arbitrary-cut lower-tail free energy and its
  exact cube/block optimality conditions.

The agreement-theorem lead came from primary-source searches, especially
Dinur--Steurer, *Direct Product Testing* (ECCC TR13-179), and
Dikstein--Dinur, *Agreement theorems for high dimensional expanders in the
small soundness regime* (arXiv:2308.09582). No theorem is imported until its
parameter regime, list/projective variant, and quantitative conclusion are
reconstructed.

All scratch work stays under `/home/math/quadra/tmp/`; Python uses
`/home/math/quadra/.venv/bin/python`. Agents must not edit tracked files and
must distinguish proofs, proposed lemmas, finite numerics, and falsifiers.
