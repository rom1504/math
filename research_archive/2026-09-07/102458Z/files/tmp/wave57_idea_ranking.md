# Wave 57 idea ranking

Read against README.md, the Wave 56 Updated frontier, and STEERING.md.
The objective below is the user-stated convergence objective; the individual
routes are research hypotheses of the main agent.

1. **Sector-split active-face/child-face dichotomy.**  Split the common
   near-ground law by orientation and retain the PSD Gram matrices
   `E[xx^T|tau]`.  Combine the exact global energy identity
   `sum_e s_e=(q-Delta)/2` with the child ground-face normal cone.  Seek a
   quantitative dichotomy: common block escape either lands in another
   locally near-ground child word (hence the bare tail) or pays enough
   exposed-face loss to support a second exact-minimality perturbation.
   Falsify this with an exact signing if no such scale-sensitive inequality
   survives.  Promise: highest; directly targets the orientation obstruction.

2. **Mesoscopic core supersaturation from a favorable incidence graph.**  At
   `ell~sqrt(nH)` and, for `c>1/12`, below `n^(5/6)`, formulate and try to prove
   the exact distinct-selector/bounded-multiplicity lemma that upgrades
   `e^{-O(H)}` incidence to `P_ell-lambda_2>=e^{-O(H)}`.  Use dependent random
   choice or slice codegree inequalities, and test whether random/matching
   incidences give a sharp counterexample.  Promise: high and independently
   testable; capacity was established in Wave 56 but correlation is open.

3. **Row-penalized selected-prior dual.**  Add captured project row to the
   favorable-incidence LP/KL variational formula (10.1356), derive its exact
   obstruction certificate, and ask whether high row on every low-information
   favorable fibre yields a field-excess perturbation contradicting exact
   minimality.  Test on stored exact minimizers.  Promise: high; attacks the
   second exact bridge without returning to the falsified uniform prior.

4. **Randomized field-weighted excess certificates.**  Randomize the choices
   of the `r_i` child-positive edges in (10.1324).  Relate expected reversed
   certificate mass to the full child energy loss, so cancellation by
   negative-to-positive edges is either quantified or exposed as a sharp
   counterexample.

5. **Iterated migration packing.**  If every near-ground state has excessive
   row, extract many disjoint high-field blocks and migrate along each.  Prove
   that bounded overlap forces exponentially many separated near-ground
   states, enough for a low-information collision, or find the exact packing
   bottleneck.

6. **Global marginal mass plus top-block control.**  The common law satisfies
   both `sum_e E s_e=(q-E Delta)/2=Theta(n^(3/2))` and
   `h_r(E s)=o(r)`.  Use cut-norm discrepancy or a Lorentz estimate to find a
   large selector class on which the internal marginal sum has controlled
   sign and variance; determine whether this improves bare local deficit.

7. **Orientation-split Gram rigidity.**  Study the two PSD correlation
   matrices from the `tau=+/-` portions of the common law.  Seek a trace or
   principal-minor inequality preventing the abstract one-state all-negative
   response on `Theta(n^(3/2))` child certificates.

8. **Child ground-face Hoffman bound.**  Quotient by all child co-maximizers
   and seek an error bound between signed disagreement with the ground face
   and local deficit.  Exact finite examples can determine whether a uniform
   constant, a scale-dependent constant, or no useful bound is possible.

9. **Fix `c=1/8` and optimize the threshold window.**  Specializing removes
   unnecessary exponent freedom: `ell~n^(13/16)` lies polynomially below
   `n^(5/6)`.  Optimize threshold location and required partner multiplicity,
   then state the weakest signing-specific cluster lemma sufficient there.

10. **Near-ground Gibbs prior with row tilt.**  Analyze
    `nu(d) proportional exp(-beta Delta(d)-gamma R(d))`; use the exact
    selected-prior chain rule to identify a free-energy inequality that would
    yield both `Z_nu>=e^{-O(H)}` and captured row.  Retain only if its partition
    function can be controlled without the retired uniform-completion tail.

Selected independent assignments: 1 (with 4,7,8 as tools), 2 (with 9), and
3 (with 5,10 as possible continuations).
