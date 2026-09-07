# Wave 16 idea ranking

The current frontier is §10.68.4.  Any proposal must respect the exact range
wall (10.570), the response identity (10.575)--(10.576), the order-nine wall
(10.579)--(10.581), and the A5/A7/A8 timing walls.

1. **Polyhedral response dual -> localized service (selected).**  Regard the
   common-mosaic response as a support function of its finite signed-state
   profile.  Dualize the minimal response increase over local minimizers and
   ask whether a separating probability distribution on active global states
   decomposes the compulsory increase into block deficits that can be placed
   on endpoint/successor buckets.  This is the most direct route from (10.576)
   to the temporal telescope and is exactly checkable on `A_9`.

2. **Exact active-state audit of the `A_9` capture wall (selected).**  Enumerate
   the 40 optimal order-six completions for the `6+1+1+1` partition, their
   active global states, local deficit marginals, endpoint pairs, and induced
   deletion chains.  Formulate the weakest plausible response-to-temporal
   inequality and test it simultaneously against A5, A7, and A8.  A finite
   counterexample would prevent another false asymptotic target; a surviving
   identity would specify the correct theorem.

3. **Active-profile seminorm instead of global `Q`-distance (selected).**
   Define a one-sided perturbation seminorm only on states active or nearly
   active for the cross mosaic.  Prove the exact deterministic inequalities
   needed for hybrid response, then seek a partition estimate for that
   seminorm which can remain small while internal excess is leading.  This
   explicitly avoids the self-defeating centered norm in (10.570).

4. **Selected-child response potential.**  Add a response/margin term to the
   normalized child excess in (10.531), then choose a principal child rather
   than average over all children.  Search for a potential whose one-step
   drop pays the response increase and whose order-only remainder is
   tail-summable under (10.529).

5. **Mixed `A_9`--conference/Hadamard amplification.**  Analyze completed
   mixed tensors and nonuniform substitutions not covered by (10.582)--
   (10.585).  Either produce a globally minimizing growing wall or prove a
   product-ground/spectral no-go.  This is an independent falsification route.

6. **Sparse order-window contraction.**  Contract the causal tree into a
   logarithmic or adaptively chosen family of order windows.  Try to bound
   (10.543)'s maximum-antichain deficiency by boundary response plus the exact
   sibling range profile (10.561), allowing timing failures inside a window
   but requiring tail-summable costs across windows.

7. **Approximate boundary-state compression.**  Replace the exact
   `2^(n-2)` min-plus boundary table by an optimizer-restricted or additive
   `o(n^(3/2))` approximation.  Determine whether its reciprocal partition
   function admits a negative-moment interpolation with `o(n^2)` logarithmic
   error.  Exact compression is already ruled out by (10.508).

8. **Conditioned partitioning on active cross states.**  Choose a partition
   after first sampling a small certificate family of active states, and
   control only their response margins.  The theorem must not imply
   `Q(E_P)=o(n^(3/2))`, which would close it by (10.570); its output should be
   a one-sided response bound plus leading internal excess.

9. **Inverse MaxCut / factorization-norm structure.**  Audit whether the
   inverse theorem of Balla--Hambardzumyan--Tomon
   (arXiv:2506.23989) or the graph Grothendieck framework of
   Alon--Makarychev--Makarychev--Naor can turn unusually small one-sided
   response into a large homogeneous rectangle, hence a recursive block with
   controlled cross response.  First translate all hypotheses and scales;
   discard the route if they only give the already known lower bound.

10. **Regularized sequence convergence.**  Define a windowed lower envelope
    or infimal convolution of `q_n/n^(3/2)` for which a weak gluing/deletion
    inequality is genuinely subadditive.  Then seek a principal-restriction
    stability estimate transferring convergence from the regularization back
    to every `q_n`.  The transfer, not formal Fekete subadditivity, is the
    concrete target.

Selected tasks 1--3 are the most direct positive routes and are mutually
independent at the technical level: finite-dimensional response duality,
exact obstruction discovery, and a partition-sensitive analytic norm.  Task
5 is the first reserve if one of them closes quickly.
