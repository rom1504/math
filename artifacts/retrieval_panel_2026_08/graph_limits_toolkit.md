# Dense graph limits, weighted graphons, action limits, and fluctuation-scale extrema

Status: retrieval packet checked through 2026-08-15; not a proposed solution.

Retrieval date: 2026-08-15. Scope: theorem-level mechanisms for obtaining limits, reconstructing finite objects at every fixed order, and deciding whether optimization over exponentially many labelings survives the chosen topology. The decisive normalization distinction is stated first.

## 1. Normalization dictionary

Let \(A_n=(a_{ij})_{i,j\le n}\) be symmetric and, unless stated otherwise, hollow. Let \(W_{A_n}\) be its step kernel on the \(n^2\) equal cells of \([0,1]^2\). All sums below are over ordered pairs \(i,j\); papers summing undirected edges once differ by a factor of two.

\[
\|W_{A_n}\|_\square
=\frac1{n^2}\max_{S,T\subseteq[n]}
\left|\sum_{i\in S,j\in T}a_{ij}\right|,
\qquad
\int W_{A_n}f_\sigma\otimes f_\sigma
=\frac1{n^2}\sigma^\top A_n\sigma .
\]

Thus ordinary bounded graphon/cut convergence is an \(n^2\)-scale theory. For a fixed spin alphabet \([q]\) and fixed coupling matrix \(J\), its natural global observable is

\[
\frac1{n^2}\max_{\phi:[n]\to[q]}
\sum_{i,j}a_{ij}J_{\phi(i),\phi(j)} .
\]

This maximizes over \(q^n\) assignments, but is nevertheless cut-continuous when \(q\) and \(J\) are fixed.

The fluctuation-scale Boolean observable has a different representation:

\[
\frac1{n^{3/2}}\max_{\sigma\in\{\pm1\}^n}
\left|\sigma^\top A_n\sigma\right|
=
\max_{\sigma\in\{\pm1\}^n}
\left|\frac1n\langle \sigma,(A_n/\sqrt n)\sigma\rangle_{\ell_2,\mathrm{sum}}\right|
=
\max_\sigma\left|\mathbb E_I[\sigma_I(\sigma A_n/\sqrt n)_I]\right|.
\]

Here \(A_n/\sqrt n\) is a matrix \(P\)-operator on the uniform probability space \([n]\), with ordinary matrix multiplication. In graphon-operator convention, \(W_{A_n}\) represents \(A_n/n\); therefore the kernel representing the same action as \(A_n/\sqrt n\) is \(\sqrt n\,W_{A_n}\), which is unbounded and whose weight laws are not tight for \(\pm1\) entries. This is exactly why bounded graphon compactness and probability-graphon tightness do not automatically capture the \(n^{3/2}\) coefficient, while action convergence can.

### What the main topologies see

| Framework | Compactness hypothesis | Stable observables | What happens to centered \(\pm1\) matrices |
|---|---|---|---|
| Bounded signed graphons, \(\delta_\square\) | uniform \(L^\infty\) bound | fixed motif densities; fixed-\(q\) ground states/free energies at \(n^2\) | quasirandom realizations converge to the zero kernel; the \(n^{3/2}\) coefficient is a convergence-rate/tangent datum, not part of the limit |
| Probability-graphons | tight edge-weight laws; automatic if weight space compact | every fixed sampled weighted subgraph; decorated homomorphism densities; fixed-\(q\) overlays/quotients | retains the one-edge law, but a linear signed ground state still vanishes at \(n^2\); scaling the kernel by \(\sqrt n\) destroys tightness |
| \(P\)-variables of Zucal | tightness of values | the same dense weighted data as probability-graphons | despite similar profile language, this is an entry-law theory, not the operator-action theory |
| Action convergence of \(P\)-operators | subsequential compactness from bounded \(\|\cdot\|_{\infty\to1}\); represented limits under a uniform \(p\to q\) bound with \(p<\infty\) | joint laws of bounded test vectors and their operator images; operator norms; rescaled matrix actions | \(A_n/\sqrt n\) can have a nontrivial limit and can retain the \(n^{3/2}\) quadratic optimum under a uniform \(2\to2\) bound |
| Higher-order graphon fluctuations | a fixed graphon \(W\), then sampling \(G(n,W)\) | second-order laws of finitely many fixed motif counts | detects fixed-motif fluctuations, not maxima over \(2^n\) configurations |
| Spin-glass/Parisi theory | random Gaussian disorder with a specified covariance | typical Boolean ground-state energy at \(n^{3/2}\) | gives a sharp probabilistic coefficient, but not a deterministic compact topology for arbitrary matrices |

## 2. A useful action-limit bridge for the Boolean extremum

The following is a short deduction from the definitions in Backhausz--Szegedy, not a theorem stated verbatim in that paper.

Suppose \(H_n\to H\) in action distance and \(\sup_n\|H_n\|_{2\to2}<\infty\). Define

\[
\Phi(H)=\sup_{\|v\|_\infty\le1}\left|\mathbb E[v(vH)]\right|.
\]

The one-profiles converge in Hausdorff--Lévy--Prokhorov distance. The uniform \(2\to2\) bound gives uniform \(L^2\) bounds on the output coordinate, hence uniform integrability of the product \(v(vH)\). It follows by taking suprema over the converging profile sets that \(\Phi(H_n)\to\Phi(H)\).

If \(H_n\) is a finite symmetric hollow matrix, its quadratic form is affine in each coordinate separately, so maxima and minima over \([-1,1]^n\) are attained at vertices. Consequently,

\[
\Phi(A_n/\sqrt n)
=\frac1{n^{3/2}}\max_{\sigma\in\{\pm1\}^n}
|\sigma^\top A_n\sigma|.
\]

This gives a plausible deterministic limit mechanism along every action-convergent, uniformly spectrally bounded subsequence. Important caveats:

- Bare \(\|\cdot\|_{\infty\to1}\) compactness supplies only first-moment tightness, not the uniform integrability used above; the \(2\to2\) hypothesis matters.
- Hollow symmetry matters for the exact cube-to-Boolean reduction. A diagonal can be separated because it contributes a fixed amount on Boolean spins.
- Backhausz--Szegedy prove only subsequential almost-sure action convergence for their iid \(\pm n^{-1/2}\) nonsymmetric matrices and explicitly leave full-sequence convergence open.
- This bridge identifies what action convergence would preserve; it does not identify the limiting operator or the spin-glass constant.

## 3. Compactness, full sequences, recovery, and extrema

These implications should not be collapsed:

1. **Compactness** gives a convergent subsequence. It does not imply that the original sequence converges or that all subsequential limits give the same extremal value.
2. **Actual full-sequence convergence plus continuity** gives convergence of continuous observables. BCLSV ground-state energies and Zucal overlays are examples at \(n^2\).
3. **\(\Gamma\)-convergence plus compactness/equicoercivity** gives convergence of minima and identifies cluster points of minimizers. Braides--Cermelli--Dovetta also construct recovery sequences at every sufficiently large index, rather than merely on a subsequence.
4. **Sampling from a limit object** gives finite approximants at every sample size and, under a coupled infinite sample, almost-sure convergence of the whole sampled sequence. It does not force exact deterministic constraints such as a prescribed number of \(+1\) spins unless a rounding/correction argument is added.
5. **All fixed orders** means “for every fixed motif size \(k\), spin number \(q\), or profile size \(k\).” It is not uniform control when that parameter grows with \(n\).

## 4. Source cards

### 4.1 Lovász--Szegedy, Limits of dense graph sequences

Primary source: [JCTB 2006, DOI](https://doi.org/10.1016/j.jctb.2006.05.002); [arXiv version](https://arxiv.org/abs/math/0408173). Peer reviewed.

- **Exact result.** For a simple graph parameter \(f\), the following are equivalent: \(f\) is a pointwise limit of homomorphism densities \(t(\cdot,G_n)\); \(f=t(\cdot,W)\) for a symmetric measurable \(W:[0,1]^2\to[0,1]\); and \(f\) satisfies normalized multiplicativity plus any of several reflection-positivity/connection-matrix conditions (Theorem MAIN).
- **Finite reconstruction.** For every graphon \(W\), the coupled \(W\)-random graphs \(G(n,W)\) converge to \(W\) almost surely for all finite \(F\). For \(F\) on \(k\) vertices,
  \[
  \Pr(|t(F,G(n,W))-t(F,W)|>\varepsilon)
  \le 2\exp[-\varepsilon^2n/(18k^2)].
  \]
- **Proof mechanism.** Nested weakly regular partitions yield a martingale of step kernels and hence a measurable limiting kernel; the converse uses \(W\)-sampling, Azuma concentration, and Borel--Cantelli.
- **Boundary.** The invariants are fixed homomorphism densities and the scale is \(n^2\). No second-order rate or \(n^{3/2}\) extremal coefficient is encoded.

### 4.2 Borgs--Chayes--Lovász--Sós--Vesztergombi I

Primary source: [Advances in Mathematics 2008, DOI](https://doi.org/10.1016/j.aim.2008.07.008); [arXiv](https://arxiv.org/abs/math/0702004). Peer reviewed.

- **Exact metric theorem.** A sequence of weighted graphs with uniformly bounded edge weights is left-convergent iff it is Cauchy in \(\delta_\square\). For edge weights in \([-1,1]\),
  \[
  |t(F,G_1)-t(F,G_2)|\le4|E(F)|\delta_\square(G_1,G_2),
  \]
  while sufficiently close densities of all \(k\)-vertex graphs imply
  \(\delta_\square(G_1,G_2)\le22/\sqrt{\log_2 k}\).
- **Compactness.** Graphons with values in a fixed finite interval, modulo cut-distance zero, form a compact metric space.
- **Testing theorem.** A bounded simple-graph parameter is testable iff it has a cut-continuous graphon extension up to an \(o(1)\) finite-size error; equivalently, it converges along every convergent graph sequence whose orders tend to infinity.
- **Proof mechanism.** Weak regularity, counting and inverse-counting lemmas, sampling, and compactness of the graphon quotient.
- **Boundary.** Compactness alone is subsequential. The metric is normalized by \(n^2\), so an \(O(n^{3/2})\) discrepancy tends to zero.

### 4.3 Borgs--Chayes--Lovász--Sós--Vesztergombi II

Primary source: [Annals of Mathematics 2012, DOI](https://doi.org/10.4007/annals.2012.176.1.2); [author-hosted PDF](https://renyi.hu/~sos/2011_Convergent_Sequences_of_Dense_Graphs_II_Multiway_Cuts_and_Statistical_Physics.pdf). Peer reviewed.

- **Exact equivalence.** For simple \(G_n\) with \(|V(G_n)|\to\infty\), left/cut convergence, Hausdorff convergence of all finite quotient sets, right convergence, convergence of every microcanonical finite-spin ground-state energy, and convergence of every microcanonical free energy are equivalent (Theorem 2.8).
- **Normalization and exponential labelings.** For fixed \(q,J,h\),
  \[
  -\widehat E(G,J,h)=
  \max_{\phi:V(G)\to[q]}
  \left\{\frac1{|V|}\sum_u h_{\phi(u)}
  +\frac2{|V|^2}\sum_{uv\in E(G)}J_{\phi(u)\phi(v)}\right\}.
  \]
  Thus the theorem controls a maximum over \(q^{|V|}\) assignments. The microcanonical version restricts color proportions to a prescribed vector \(a\).
- **Important non-equivalence.** Unconstrained ground-state/free-energy convergence is implied by left convergence but is strictly weaker; naive right convergence and spectral convergence are also weaker.
- **Proof mechanism.** Weak regular step approximations; quotient-set compactness; support-function/convex-geometric recovery of quotients from energies; rounding fractional partitions to vertex partitions.
- **Boundary.** \(q\) and \(J\) are fixed and the energy is divided by \(n^2\). For centered sign matrices converging to zero in cut distance, the theorem correctly returns a zero graphon-scale Boolean energy and says nothing about its \(n^{3/2}\) coefficient.

### 4.4 Diaconis--Janson, Graph limits and exchangeable random graphs

Primary source: [arXiv:0712.2749](https://arxiv.org/abs/0712.2749). Published in Rendiconti di Matematica, 2008.

- **Exact correspondence.** Distributions of random proper graph limits are in bijection with distributions of exchangeable infinite random graphs \(H\), via
  \[
  \mathbb E\,t_{\rm ind}(F,\Gamma)=\Pr(H|_{[k]}=F).
  \]
  Moreover \(H|_{[n]}\to\Gamma\) almost surely in graph-limit space.
- **Extreme points.** Deterministic graph limits correspond exactly to extreme exchangeable laws. Extremality is equivalent to dissociation of restrictions on disjoint vertex sets and to triviality of the relevant tail sigma-field.
- **Representation.** Every exchangeable graph is a mixture of \(G(\infty,W)\) laws; deterministic \(W\) gives an extreme law, modulo graphon weak isomorphism.
- **Proof mechanism.** Projective consistency/Kolmogorov extension, reverse martingales, and the Aldous--Hoover representation.
- **Boundary.** This is a representation and distributional correspondence, not a quantitative topology for fluctuation-scale optimization. A global mixing variable is essential for non-dissociated arrays.

### 4.5 Abraham--Delmas--Weibel, Probability-graphons

Primary source: [Innovations in Graph Theory 2025, DOI](https://doi.org/10.5802/igt.7); [arXiv](https://arxiv.org/abs/2312.15935). Peer reviewed.

- **Objects/topology.** A probability-graphon is a measurable probability kernel \(W(x,y;dz)\) on a Polish weight space. The weak-isomorphism quotient equipped with the Lévy--Prokhorov cut metric is Polish.
- **Compactness theorem.** A family is relatively cut-compact iff its averaged edge-weight laws are tight. If the weight space is compact, the entire quotient is compact.
- **Characterization.** Cut convergence is equivalent both to convergence of every homomorphism density decorated by bounded continuous edge functions and to convergence in distribution of \(G(k,W_n)\) for every fixed \(k\ge2\).
- **All-size reconstruction.** Under one coupled sample, \(G(k,W)\to W\) almost surely. The second sampling lemma gives, with probability at least \(1-\exp[-k/(2\log k)]\),
  \[
  \delta_{\square,\mathcal F}(G(k,W),W)\le22/\sqrt{\log k};
  \]
  the analogous latent weighted graph has constant \(21\).
- **Proof mechanism.** Measure-valued weak regularity, Prokhorov tightness, decorated counting/inverse-counting, sampling concentration, and Borel--Cantelli.
- **Boundary.** The object records conditional edge-weight laws at dense scale. It does not by itself retain realization-specific global correlations or a subleading linear signed ground state. The kernel rescaling needed to promote \(n^{3/2}\) to \(n^2\) is non-tight.

### 4.6 Zucal, Probability graphons: the right convergence point of view

Primary source: [arXiv:2407.05998](https://arxiv.org/abs/2407.05998). ArXiv preprint; no journal publication was verified by the retrieval date.

- **Exact theorem.** For \(W_n\to W\), unlabelled probability-graphon cut convergence is equivalent to convergence of every overlay functional \(\mathcal C(W_n,G^\beta)\), and to Hausdorff convergence of every fixed-\(k\) quotient set \(\mathcal Q_k(W_n)\). The Cauchy formulation assumes tightness.
- **Why it handles exponentially many assignments.** For a decorated graph \(G^\beta\) on \([k]\),
  \[
  \mathcal C(W,G^\beta)
  =\sup_{(S_1,\ldots,S_k)}
  \sum_{i,j}\int_{S_i\times S_j}\int\beta_{ij}(z)\,W(x,y;dz)\,dx\,dy,
  \]
  with prescribed part masses when desired. Finite step kernels turn this into a maximum over \(k^n\) labelings.
- **Proof mechanism.** Cut-continuity of overlays, finite decorated step approximation, compact quotient sets, and separation by support functionals. The paper explicitly notes that the implication from quotient/overlay data back to cut convergence is non-effective.
- **Boundary.** Fixed \(k\), bounded continuous decorations, and dense \(n^2\) normalization. This strengthens the global interpretation of probability-graphons but does not recover subleading extrema.

### 4.7 Zucal, Probability graphons and P-variables

Primary source: [arXiv:2408.07572](https://arxiv.org/abs/2408.07572). ArXiv preprint; no journal publication was verified by the retrieval date.

- **Objects.** A \(P\)-variable is a measurable \(W:\Omega_1^2\times\Omega_2\to\mathbb R\). Its \(k\)-profile consists of joint laws of
  \[
  (f_1(x),f_1(y),\ldots,f_k(x),f_k(y),W(x,y,z))
  \]
  for bounded vertex test functions. This should not be confused with the operator profiles \((v,vA)\) of action convergence.
- **Exact equivalence.** On tight sequences, convergence in the \(P\)-variable profile metric is equivalent to cut convergence of the associated probability-graphons.
- **Compactness/reconstruction.** Tightness is equivalent to relative compactness. Aldous--Hoover-style matrices \(M_{ij}^{(n)}=W(X_i,X_j,Y_{ij})\) converge to \(W\) almost surely in the profile metric; deterministic finite weighted graphs are dense.
- **Proof mechanism.** Reduce profiles to function partitions and quotient data, then invoke probability-graphon right convergence and sampling.
- **Boundary.** The model has vertex variables and independent pair variables but no global \(U_0\), so it does not encode arbitrary nonlocal random dependence. For raw \(\pm1\) weights it retains the law; for the graphon rescaling \(\sqrt n\,A_n\) needed by the \(n^{3/2}\) energy, tightness fails.

### 4.8 Backhausz--Szegedy, Action convergence of operators and graphs

Primary source: [Canadian Journal of Mathematics 2022, DOI](https://doi.org/10.4153/S0008414X2000070X); [arXiv](https://arxiv.org/abs/1811.00626). Peer reviewed.

- **Objects.** A \(P\)-operator is a bounded linear map \(A:L^\infty(\Omega)\to L^1(\Omega)\). Its \(k\)-profile is the set of laws of
  \[
  (v_1,\ldots,v_k,v_1A,\ldots,v_kA),\qquad |v_i|\le1.
  \]
  Action distance is the weighted sum of Hausdorff--Lévy--Prokhorov distances between these profile sets.
- **Compactness.** Uniform \(\infty\to1\) boundedness gives an action-convergent subsequence. If an action-Cauchy sequence is uniformly \(p\to q\) bounded with \(p<\infty\), it has a representing limit \(P\)-operator with the same bound. For \(q>1\), the bounded weak-equivalence classes form a compact metric space.
- **Unification.** Restricted to graphons, action convergence equals graphon convergence; restricted to graphings, it equals local-global convergence.
- **Random matrices.** For iid mean-zero \(\pm n^{-1/2}\) nonsymmetric matrices \(H_n\), every infinite set of sizes contains a further infinite set along which \(H_n\) converges almost surely in action distance to a \(P\)-operator. The full sequence is explicitly left open.
- **Proof mechanism.** Profile tightness and diagonal compactness; construction of a probability-space/operator representative from a countable algebra; for random matrices, a column-exposure Azuma argument plus spectral-norm bounds.
- **Boundary.** The paper does not state a Parisi formula, identify the Wigner action limit, or state the Boolean-energy continuity deduction in Section 2. It does provide the topology and compactness in which \(A_n/\sqrt n\), rather than \(A_n\), is the natural nontrivial object.

### 4.9 Hatami--Lovász--Szegedy, Limits of locally-globally convergent graph sequences

Primary source: [GAFA 2014, DOI](https://doi.org/10.1007/s00039-014-0258-7); [arXiv](https://arxiv.org/abs/1205.4356). Peer reviewed.

- **Exact framework.** For graphs of maximum degree at most \(d\), local-global convergence requires Hausdorff convergence, for every radius \(r\) and color number \(k\), of the distributions of rooted colored \(r\)-balls obtainable from all \(k\)-colorings.
- **Limit theorem.** Every local-global convergent bounded-degree sequence has a graphing representative, with all coloring profiles reproduced in the limit.
- **Proof mechanism.** Ultraproducts, separable invariant sub-sigma-fields, and approximation of measurable colorings.
- **Boundary.** The degree bound is essential in this formulation. Dense signed matrices scaled by \(1/\sqrt n\) are outside it; the relevance here is conceptual and through the action-convergence theorem that recovers local-global convergence on graphings.

### 4.10 Braides--Cermelli--Dovetta, \(\Gamma\)-limit of the cut functional

Primary source: [ESAIM COCV 2020, DOI](https://doi.org/10.1051/cocv/2019029); [arXiv](https://arxiv.org/abs/1806.03436). Peer reviewed.

- **Exact setup.** Let \(W_n\to W\) in cut norm after relabeling, let \(L=\{\ell_1,\ldots,\ell_N\}\) be a fixed finite label set, and let \(f:L^2\to\mathbb R\). For vertex-label functions,
  \[
  F_n(u)=\int W_n(x,y)f(u(x),u(y))\,dx\,dy
  =\frac1{n^2}\sum_{i,j}A^{(n)}_{ij}f(u_i,u_j).
  \]
- **\(\Gamma\)-limit.** In the narrow Young-measure topology,
  \[
  F_n\ \Gamma\hbox{-converges to}
  I(\nu)=\int W(x,y)\!\int f(\lambda,\mu)
  \,d\nu_x(\lambda)d\nu_y(\mu)\,dx\,dy.
  \]
- **Recovery at every order.** The proof approximates the Young-measure mass functions by rational step functions and interlaces pure labels in the \(n\)-mesh to produce a full recovery sequence.
- **Constrained minima.** If prescribed class sizes satisfy \(j_k^{(n)}/n\to j_k\), the constrained minima converge; Young-measure cluster points of minimizers minimize the constrained limit functional (Proposition 14).
- **Proof mechanism.** A continuous-convergence lemma combines cut-norm uniformity over all rectangles with weak-star convergence of label mass functions; then explicit chattering/rounding supplies limsup recovery.
- **Boundary.** Full cut convergence is assumed rather than obtained from compactness. At the zero graphon, the \(\Gamma\)-limit is zero, so different \(n^{3/2}\) minima remain indistinguishable.

### 4.11 Chatterjee--Varadhan, LDP for \(G(n,p)\)

Primary source: [European Journal of Combinatorics 2011, arXiv](https://arxiv.org/abs/1008.1946). Peer reviewed.

- **Exact theorem.** The laws of \(G(n,p)\) in the graphon quotient \((\widetilde{\mathcal W},\delta_\square)\) satisfy an LDP with speed \(n^2\) and good rate
  \[
  I_p(W)=\frac12\int\!\left[
  W\log\frac Wp+(1-W)\log\frac{1-W}{1-p}
  \right]dx\,dy.
  \]
  Closed sets get the upper bound and open sets the lower bound.
- **Conditioning consequence.** On a closed rare-event set satisfying matching interior/closure infima, conditional samples concentrate in cut distance near the compact set of rate-function minimizers.
- **Proof mechanism.** Product Bernoulli large deviations in a weak topology, weak regularity and finite covers for exponential approximation, lower semicontinuity, and compactness of graphon quotient.
- **Boundary.** Speed \(n^2\) describes dense large deviations. Typical \(n^{3/2}\) optimization is a fluctuation/moderate-scale phenomenon that collapses at this topology and speed.

### 4.12 Dionigi--Zucal, Large deviations for probability graphons

Primary source: [arXiv:2509.14204](https://arxiv.org/abs/2509.14204). Submitted preprint as of the retrieval date; treat as unrefereed.

- **Exact claimed theorem.** Let the edge-weight space \(Z\) be compact Polish and let all undirected edge weights be iid with reference law \(\nu\). Their probability-graphon laws satisfy an LDP in unlabelled cut distance with speed \(n^2/2\) and good rate
  \[
  I_\nu(W)=\int_{[0,1]^2}D(W(x,y)\|\nu)\,dx\,dy.
  \]
  Equivalently, the paper writes the bounds using \((2/n^2)\log\Pr(\cdot)\).
- **Proof mechanism.** Prove the finite weight-space theorem using weak regularity, finite-dimensional LDP bounds, and exponential tilting; pass to compact Polish \(Z\) using finite discretizations and Dawson--Gärtner; identify the projective rate with integrated relative entropy.
- **Boundary.** Independent edges, compact weight space, and dense large-deviation speed. It is not a theorem about correlated disorders or the \(n^{3/2}\) typical ground state.

### 4.13 Chatterjee--Dan--Bhattacharya, Higher-Order Graphon Theory

Primary source: [arXiv:2404.13822](https://arxiv.org/abs/2404.13822). Listed by the Institute of Mathematical Statistics as a forthcoming Annals of Statistics paper at the retrieval date; arXiv version consulted.

- **Exact scaling theorem.** For \(G_n\sim G(n,W)\), a fixed nonempty graph \(H\) with \(h\) vertices, and \(X(H,G_n)\) the number of unlabeled copies, center by
  \[
  (n)_h\,t(H,W)/|\operatorname{Aut}(H)|.
  \]
  If \(W\) is \(H\)-irregular, divide by \(n^{h-1/2}\) and obtain a Gaussian linear Wiener integral. If \(W\) is \(H\)-regular, divide by \(n^{h-1}\) and obtain the sum of an independent Gaussian edge-noise term and a second Wiener integral.
- **Marginal regular case.** The second Wiener integral can be written as a possibly infinite weighted sum \(\sum_\lambda\lambda(Z_\lambda^2-1)\), with weights from the spectrum of the two-point conditional kernel \(W_H\).
- **Joint theorem.** The paper gives the joint limit for any finite collection containing both regular and irregular motifs, preserving their dependence through common Wiener integrals.
- **Proof mechanism.** Generalized \(U\)-statistic/Hoeffding projections; the first nonzero projection determines the normalization; multiple Wiener--Itô integrals encode the limit.
- **Boundary.** These are fluctuations of finitely many fixed motifs. Higher degeneracy can make even the stated regular-scale limit vanish. The theorem does not cover a supremum over \(2^n\) spin assignments.

### 4.14 Auffinger--Chen, zero-temperature Parisi formula

Primary source: [Annals of Probability 2017, DOI](https://doi.org/10.1214/16-AOP1173); [arXiv](https://arxiv.org/abs/1606.05335). Peer reviewed.

- **Model/normalization.** On \(\Sigma_N=\{\pm1\}^N\),
  \[
  H_N(\sigma)=
  \sum_{p\ge2}\frac{c_p}{N^{(p-1)/2}}
  \sum_{i_1,\ldots,i_p}g_{i_1\ldots i_p}
  \sigma_{i_1}\cdots\sigma_{i_p}
  +h\sum_i\sigma_i,
  \]
  with \(\sum2^pc_p^2<\infty\) and covariance \(N\xi(R)\). The SK case is \(\xi(s)=s^2/2\).
- **Exact theorem.** The ground state per spin
  \(L_N=N^{-1}\max_\sigma H_N(\sigma)\) converges almost surely and
  \[
  \operatorname{GSE}
  =\inf_{\gamma\in\mathcal U}
  \left\{\Psi_\gamma(0,h)
  -\frac12\int_0^1t\xi''(t)\gamma(t)\,dt\right\},
  \]
  where \(\Psi_\gamma\) solves the zero-temperature Parisi PDE with terminal value \(|x|\).
- **Relevance to \(n^{3/2}\).** In the \(p=2\) term, the Hamiltonian contains a raw Gaussian quadratic sum divided by \(\sqrt N\). A finite nonzero limit of \(N^{-1}\max H_N\) therefore means that the corresponding raw quadratic maximum is of order \(N^{3/2}\), with convention-dependent constant supplied by the Parisi variational problem.
- **Proof mechanism.** Take the zero-temperature limit of the finite-temperature Parisi formula; use a stochastic-control representation of the PDE, compactness of scaled order parameters away from \(t=1\), and cancellation of the possible terminal atom/singularity.
- **Boundary.** Gaussian mean-field disorder with a prescribed covariance. Extension to Bernoulli signs requires a separate universality result; this paper is not a deterministic graph-limit theorem.

### 4.15 Alon--Naor, cut norm via Grothendieck

Primary source: [SIAM Journal on Computing 2006, DOI](https://doi.org/10.1137/S0097539704441629); [author-hosted PDF](https://web.math.princeton.edu/~naor/homepage%20files/cutnorm.pdf). Peer reviewed.

- **Exact algorithmic result.** For a real matrix \(A\), computing
  \(\|A\|_C=\max_{I,J}|\sum_{i\in I,j\in J}a_{ij}|\) is MAX-SNP hard, but a polynomial-time semidefinite relaxation and rounding algorithm finds \(I,J\) with value at least \(\rho\|A\|_C\), for an absolute \(\rho>0.56\).
- **Mechanism.** Relate the rectangle norm to the bilinear sign norm
  \(\max_{x_i,y_j\in\{\pm1\}}\sum a_{ij}x_iy_j\), apply Grothendieck's inequality to its vector SDP relaxation, and use Krivine-type randomized rounding.
- **Relevance.** This is a robust computational way to approximate global extrema over exponentially many rectangles/sign pairs, and it underlies algorithmic weak regularity.
- **Boundary.** The bilinear signs \(x,y\) are independent. A same-spin quadratic restriction \(x=y\) is a different optimization problem; the paper supplies comparison/approximation machinery, not a graph-limit compactness theorem or an asymptotic coefficient.

## 5. Practical synthesis

The strongest literature-backed route depends on the desired scale:

- For any fixed finite spin model normalized by \(n^2\), use cut convergence plus BCLSV II or probability-graphon overlays. These results already handle exponentially many assignments and, with actual full-sequence convergence, give full-sequence convergence of the optima.
- For recovery of nearly optimal finite labelings from a continuum minimizer, use the explicit Young-measure recovery construction of Braides--Cermelli--Dovetta; it works at every sufficiently large index and supports asymptotic class-size constraints.
- For weighted random edges and all fixed finite sample orders, use probability-graphon compactness/sampling; add a global Aldous--Hoover variable when the random array has nonlocal dependence.
- For \(\max_\sigma\sigma^\top A_n\sigma/n^{3/2}\), the scale-correct deterministic candidate is action convergence of \(A_n/\sqrt n\) under a uniform spectral bound. Ordinary graphon or probability-graphon limits do not carry this coefficient.
- For iid mean-field disorder, action compactness supplies possible limit objects, while the Parisi formula supplies the actual typical ground-state value. Identifying a universal full action limit and proving that it reproduces the Parisi coefficient remain separate tasks.

## 6. Avoidable overclaims

- Cut convergence of \(A_n\) to zero does not imply that the \(n^{3/2}\)-normalized Boolean maximum tends to zero.
- Probability-graphon convergence of raw \(\pm1\) entries to the constant law \((\delta_{-1}+\delta_{+1})/2\) does not determine realization-specific subleading extrema without additional independence/universality hypotheses.
- Compactness does not upgrade subsequences to a full sequence; a uniqueness or Cauchy argument is still needed.
- Convergence for every fixed \(q\) does not imply uniform control for \(q=q(n)\).
- The random-matrix proposition in Backhausz--Szegedy is subsequential and nonsymmetric as stated.
- The Auffinger--Chen theorem is Gaussian; Bernoulli universality should be cited separately before transferring its constant.
- The 2024 Zucal papers and the 2025 probability-graphon LDP were treated here according to their preprint status; the latter was explicitly still submitted/unrefereed in the sources checked.
