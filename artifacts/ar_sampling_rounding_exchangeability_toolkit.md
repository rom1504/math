# AR literature toolkit: sampling, rounding, exchangeability, and prescribed statistics

Date: 2026-08-16.

This packet was assembled before comparing the mechanisms with the archived AR
attempts. It records theorem hypotheses at the scale relevant to
​\(T_A=A/\sqrt n\), rather than collecting graph-limit analogies by keyword.

## 1. Dense sampling and ordinary graphons

1. **Borgs--Chayes--Lovász--Sós--Vesztergombi, _Convergent
   sequences of dense graphs I_**
   ([arXiv:math/0702004](https://arxiv.org/abs/math/0702004)). Their sampling
   theorem (Theorem 4.7 in the source numbering) approximates a bounded
   graphon in cut distance by an order-​\(n\) sample, at every (n), with an
   error tending to zero; the general worst-case bound is logarithmic rather
   than a fluctuation-scale estimate. This is an all-order theorem for a
   uniformly bounded kernel.

2. **Lovász--Szegedy, _Limits of dense graph sequences_**
   ([DOI 10.1007/s00039-006-0273-5](https://doi.org/10.1007/s00039-006-0273-5)).
   A (W)-random graph converges almost surely to its bounded graphon, and
   finite simple graphs are dense in the graphon cut topology.

3. **Backhausz--Szegedy, _Action convergence of operators and graphs_**
   ([arXiv:1811.00626](https://arxiv.org/abs/1811.00626)). Theorem 1.2 gives
   subsequential compactness under a uniform ​\(\infty\to1\) bound and a
   (P)-operator limit under a uniform ​\(p\to q\) bound with finite
   (p,q). It is a compactness theorem, not a finite-model recovery theorem.
   Section 11 is especially diagnostic: normalized i.i.d. sign matrices
   concentrate in action distance around deterministic representatives at
   every order (Lemma 11.2), but Proposition 11.1 still passes to a
   subsequence. The paper explicitly leaves convergence along all natural
   orders open.

4. **Braides--Cermelli--Dovetta, _Gamma-limit of the cut functional on
   dense graph sequences_**
   ([arXiv:1806.03436](https://arxiv.org/abs/1806.03436)). The recovery
   sequence rounds partition/Young-measure variables while the dense graph
   kernels are already assumed to converge as graphons. It does not recover
   the edge kernel itself from a fluctuation-scale operator object.

**Scale translation.** On an (n)-point uniform space, the integral kernel
representing the standard matrix operator (A/\sqrt n) has entries
​\(\sqrt n,a_{ij}\). Its ​\(L^p\) norms diverge like ​\(\sqrt n\). Hence
bounded-graphon sampling theorems do not apply. Centering a sign matrix gives
ordinary graphon limit zero while erasing exactly the ​\(n^{3/2}\) Boolean
fluctuation retained by action convergence.

## 2. Exchangeability and projective realization

5. **Diaconis--Janson, _Graph limits and exchangeable random graphs_**
   ([arXiv:0712.2749](https://arxiv.org/abs/0712.2749)). Theorem 5.3 and
   Sections 5--6 identify infinite jointly exchangeable graph arrays with
   random graphon mixtures. Conditional on the graphon and i.i.d. vertex
   labels, edges are independent.

6. **Aldous, _Representations for partially exchangeable arrays of random
   variables_**
   ([J. Multivar. Anal. 11 (1981)](https://doi.org/10.1016/0047-259X(81)90099-3)),
   together with Hoover's array representation, supplies the general
   representation behind item 5.

The exact consequence for AR is proved in
`exchangeable_recovery_obstruction.md`: projective exchangeability plus
tight ​\(\lVert A_n\rVert_{op}/\sqrt n\) forces the graphon mean kernel to
vanish, hence forces i.i.d. Rademacher edges. A direct online-greedy proof
then gives cap coefficient at least
​\((2/3)\sqrt{2/\pi}=0.531923\ldots>1/2\). Thus projective sampling is not
an extremal recovery mechanism. Vertex-exchangeability separately at each
order imposes no comparable restriction: one may take the uniform orbit of
any deterministic signing, so it supplies no construction by itself.

## 3. Microcanonical graphons and fluctuation theory

7. **Chatterjee--Varadhan, _The large deviation principle for the
   Erdős--Rényi random graph_**
   ([arXiv:1008.1946](https://arxiv.org/abs/1008.1946)). The LDP is in graphon
   cut topology at speed (n^2), with a graphon entropy rate function.

8. **den Hollander--Mandjes--Roccaverde--Starreveld, _Ensemble equivalence
   for dense graphs_**
   ([arXiv:1703.08058](https://arxiv.org/abs/1703.08058)). Microcanonical
   constraints are exact per realization, but equivalence is measured by
   relative entropy divided by (n^2). Their main examples also show genuine
   breaking of ensemble equivalence under frustrated finite constraints.

9. **den Hollander--Markering, _Breaking of ensemble equivalence for dense
   random graphs under a single constraint_**
   ([arXiv:2107.04351](https://arxiv.org/abs/2107.04351)). Even one fixed
   subgraph-density constraint can have a positive speed-​\(n^2\) entropy
   gap and a spectral signature.

10. **Chatterjee--Dan--Bhattacharya, _Higher-Order Graphon Theory:
    Fluctuations, Degeneracies, and Inference_**
    ([arXiv:2404.13822](https://arxiv.org/abs/2404.13822)), and
    **Bhattacharya--Chatterjee--Janson, _Fluctuations of Subgraph Counts in
    Graphon Based Random Graphs_**
    ([arXiv:2104.07259](https://arxiv.org/abs/2104.07259)), describe Gaussian
    and higher-order degenerate fluctuations for any fixed collection of
    motif counts.

These results constrain finitely many ordinary dense statistics. Balanced
signings with normalized caps (0.5) and i.i.d. signings with cap above
​\(0.5319\) have the same zero ordinary graphon limit. Fixed motif
fluctuations therefore do not control the universal spin maximum or the
directed action one-profile. No cited microcanonical theorem works at the
required local-global extreme scale.

## 4. Matrix discrepancy and sign rounding

11. **Marcus--Spielman--Srivastava, _Interlacing Families II_**
    ([Annals 182 (2015)](https://doi.org/10.4007/annals.2015.182.1.8)). Their
    Kadison--Singer/Weaver theorem signs rank-one positive semidefinite
    summands with operator discrepancy
    ​\(O(\sqrt\alpha+\alpha)\) under isotropy and individual norm bound
    ​\(\alpha\).

12. **Kyng--Luh--Song, _Four deviations suffice for rank 1 matrices_**
    ([arXiv:1901.06731](https://arxiv.org/abs/1901.06731)). For independent
    finite-support scalar choices on rank-one outer products, some outcome
    differs from its mean by operator norm at most four times the intrinsic
    matrix standard deviation.

13. **Dadush--Jiang--Reis, _A new framework for matrix discrepancy_**
    ([arXiv:2111.03171](https://arxiv.org/abs/2111.03171)). Their mirror-descent
    partial-coloring framework proves Matrix Spencer for low-rank and
    block-diagonal families, with discrepancy on the natural square-root
    scale (and logarithmic factors in the general parameter ranges).

14. **Bandeira--Bölcskei, _Matrix Discrepancy for Representations of Finite
    Groups_**
    ([arXiv:2606.12181](https://arxiv.org/abs/2606.12181)). For every finite
    group (G), signs can be chosen so that the signed sum in the left
    regular representation has operator norm ​\(O(\sqrt{|G|})\).

15. **Reis, _An Algebraic Matrix Spencer Theorem_**
    ([arXiv:2606.16005](https://arxiv.org/abs/2606.16005)). Theorem 1.3 proves
    an ​\(O(\sqrt n)\) discrepancy bound for (n) operator-norm-bounded
    matrices lying in a finite-dimensional (C^*)-algebra whose intrinsic
    dimension is ​\(O(n)\); the paper also treats low-rank perturbations.

All these are constant-times-square-root bounds. For rounding a dense
weighted order-​\(n\) matrix, an operator residual (c\sqrt n) changes
​\(Q\) by as much as ​\(c n^{3/2}/2\), a fixed leading amount. AR needs
​\(o(\sqrt n)\) operator error, or a genuinely joint same-spin cancellation
theorem that is not a norm bound. Moreover, if a positive fraction of a
weighted matrix's entries stay a fixed distance from ​\(\{\pm1\}\), every
sign rounding has Frobenius residual ​\(\Omega(n)\) and hence operator
residual ​\(\Omega(\sqrt n)\). Thus no uniform little-​\(o\) operator
rounding theorem is possible for genuinely fractional inputs.

## 5. Fixed constraints and absorption

16. **Kuperberg--Lovett--Peled, _Probabilistic existence of regular
    combinatorial structures_**
    ([arXiv:1302.4295](https://arxiv.org/abs/1302.4295)). A lattice local
    central limit theorem yields exact balanced subsets for a fixed
    finite-dimensional function space under symmetry, divisibility, and
    bounded integer-basis hypotheses.

17. **Keevash, _The existence of designs_**
    ([arXiv:1401.3665](https://arxiv.org/abs/1401.3665)). Natural divisibility
    is sufficient, for all sufficiently large admissible orders, for fixed
    design parameters; the general decomposition theorem assumes
    pseudorandomness/extendability and a robust fractional decomposition.

These are genuine all-order (subject to explicit divisibility) absorption
mechanisms for fixed local templates. At fixed ​\((C,\epsilon)\), the AR
one-profile has a finite abstract net, as proved in
`minimal_all_order_action_recovery.md`. The unresolved mismatch is the
quantifier: a realizing signing must satisfy the outer-profile condition for
all (q^n) fixed-alphabet vertex colorings. The cited design theorems enforce
fixed local linear statistics, not a universal local-global coloring profile.

## 6. Frozen toolkit verdict

The literature supplies no all-order recovery theorem for the required
fluctuation object. It gives four exact boundaries:

1. bounded graphon sampling is at the wrong normalization;
2. action convergence itself remains subsequential even for random signs;
3. projective exchangeability collapses to a nonextremal i.i.d. object; and
4. modern rounding/absorption theorems stop at a fixed leading discrepancy or
   at finitely many local constraints.

The only plausible theorem shape left within these domains is a **joint
directed-profile rounding/absorption theorem** whose forced microscopic sign
residual is absorbed before taking the supremum. No theorem in this packet
provides such absorption.
