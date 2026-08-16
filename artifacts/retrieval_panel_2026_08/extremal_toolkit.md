# Extremal/probabilistic-combinatorics retrieval toolkit

Status: retrieval packet checked through 2026-08-15; not a proposed solution.

Retrieved through 2026-08-15.  This is a theorem-and-mechanism packet, not a proposed solution.  The ambient normalization used to compare papers is a symmetric matrix (A=(a_{ij})), (a_{ii}=0), (a_{ij}\in\{\pm1\}), with

\[
 Q_A(S):=\sum_{\{i,j\}\subset S}a_{ij}=\tfrac12\mathbf1_S^TA\mathbf1_S,
 \qquad D_{\rm ind}(A):=\max_{S\subseteq[n]}|Q_A(S)|.
\]

For disjoint (S,T), (A(S,T)=Q_A(S\cup T)-Q_A(S)-Q_A(T)).  Hence every disjoint rectangle is at most (3D_{\rm ind}(A)); decomposing arbitrary (X,Y) into (X\cap Y,X\setminus Y,Y\setminus X) gives

\[
 2D_{\rm ind}(A)\le \|A\|_{\square}:=\max_{X,Y}|\mathbf1_X^TA\mathbf1_Y|\le 11D_{\rm ind}(A).
\]

These are only constant-factor bridges: they do not preserve an exact leading constant.  Also (A(S,S^c)=Q_A([n])-Q_A(S)-Q_A(S^c)), so a cut-only assertion should not silently be identified with an induced-sum assertion.

## Direct discrepancy and signed-graph results

### 1. Random unconditional convergence of Rademacher chaos in (L_\infty) and sharp estimates for discrepancy of weighted graphs and hypergraphs

**Sergey V. Astashkin and Konstantin V. Lykov (arXiv 2024; Math. Ann. 2025).** [arXiv:2412.20107](https://arxiv.org/abs/2412.20107); [DOI 10.1007/s00208-025-03257-9](https://doi.org/10.1007/s00208-025-03257-9).

- **Theorem/normalization.** For an edge-weighted graph (G=(V,E,w)), define
  \[
  \operatorname{disc}(G)=\min_{\theta\in\{\pm1\}^{E}}\max_{U\subseteq V}
  \left|\sum_{e\in E(G[U])}\theta_e w_e\right|.
  \]
  Their Theorem 7 gives, with universal implicit constants,
  \[
  \operatorname{disc}(G)\asymp
  \mathbb E_\theta\max_{U\subseteq V}\left|\sum_{e\in E(G[U])}\theta_e w_e\right|
  \asymp \sum_{v\in V}\left(\sum_{e\ni v}w_e^2\right)^{1/2}.
  \]
  Thus (w\equiv1) on (K_n) gives (\Theta(n\sqrt{n-1})=\Theta(n^{3/2})).  Theorem 8 gives the analogous estimate for a complete weighted (d)-uniform hypergraph, with constants depending only on fixed (d); the unweighted scale is (n^{(d+1)/2}).
- **Mechanism.** Identify the induced-set maximum with an (L_\infty)-norm of second-order Rademacher chaos and a modified cut norm; decouple to the multiple Rademacher system; prove random unconditional convergence; reduce the matrix norm to sums of row/column ℓ2 norms through Khintchine-type estimates.
- **Quantifiers/order.** Every finite weighted graph, hence every order.  The lower estimate applies to **every signing**; the upper estimate is existence of a signing, and even the expected performance of a uniformly random signing.  It is not merely a special construction family.
- **Boundary.** The constants hidden by (\asymp) are not identified, so this settles exact order but neither a leading constant nor convergence of the normalized optimum.

### 2. Imbalances in (k)-colorations

**Paul Erdős and Joel Spencer (1971/72).** *Networks* 1, 379–385. [DOI 10.1002/net.3230010407](https://doi.org/10.1002/net.3230010407).

- **Theorem/normalization.** For every fixed uniformity (r\ge2), the minimum, over two-colourings (\theta:E(K_n^{(r)})\to\{\pm1\}), of
  \(
  \max_{U\subseteq[n]}|\sum_{e\in\binom U r}\theta_e|
  \)
  is between (c_r n^{(r+1)/2}) and (C_r n^{(r+1)/2}).  For (r=2) this is the foundational (\Theta(n^{3/2})) result.
- **Mechanism.** Moment/averaging estimates force a large imbalance for every colouring; a random colouring plus tail estimates/union bounds supplies the matching-order upper bound over all (2^n) vertex subsets.
- **Quantifiers/order.** The lower bound is uniform over all edge signings; the upper bound selects a random/existential signing.  All sufficiently large (n), not an arithmetic subsequence.
- **Boundary.** Fixed-factor order only.  It does not produce an exact constant or a finite-size recurrence.

### 3. Discrepancy in graphs and hypergraphs

**Béla Bollobás and Alexander Scott (2006).** In *More Sets, Graphs and Numbers*. [Author PDF](https://people.maths.ox.ac.uk/~scott/Papers/disc.pdf); [DOI 10.1007/978-3-540-32439-3_3](https://doi.org/10.1007/978-3-540-32439-3_3).

- **Theorem/normalization.** If (G) has order (n), density (p=e(G)/\binom n2), and (p(1-p)\ge1/n), put
  \[
  \operatorname{disc}_p^+(G)=\max_U\{e(G[U])-p\tbinom{|U|}{2}\},\quad
  \operatorname{disc}_p^-(G)=\max_U\{p\tbinom{|U|}{2}-e(G[U])\}.
  \]
  Their Theorem 1 proves
  \[
  \operatorname{disc}_p^+(G)\operatorname{disc}_p^-(G)
  \ge \frac{p(1-p)n^3}{6400},
  \]
  hence (\operatorname{disc}_p(G)\ge \sqrt{p(1-p)},n^{3/2}/80).
- **Mechanism.** Weighted discrepancy and averaging force both a positively and a negatively deviating induced subgraph; the product estimate prevents one side from being hidden by global density bias.
- **Quantifiers/order.** Uniform for every graph satisfying the density hypothesis and every order.
- **Boundary.** The centring is by the graph's own density (p), whereas a signed matrix is centred at zero.  Translating between them requires retaining the total edge bias.  The numerical constant is a lower-bound constant, not an optimum constant.

### 4. Positive discrepancy, MaxCut, and eigenvalues of graphs

**Eero Räty, Benny Sudakov and István Tomon (arXiv 2023; Trans. AMS 2026).** [arXiv:2311.02070](https://arxiv.org/abs/2311.02070); [DOI 10.1090/TRAN/9551](https://doi.org/10.1090/TRAN/9551).

- **Theorem/normalization.** With (p=e(G)/\binom n2) and (d=2e(G)/n), their one-sided discrepancy is the (\operatorname{disc}_p^+) above.  For fixed (\varepsilon>0), they prove
  \[
  \operatorname{disc}_p^+(G)\ge
  \begin{cases}
  \Omega(\sqrt d\,n),&1\le d\le n^{2/3},\\
  \Omega(n^2/d),&n^{2/3}\le d\le n^{4/5},\\
  \Omega(d^{1/4}n/\log n),&n^{4/5}\le d\le(1/2-\varepsilon)n.
  \end{cases}
  \]
  The first two regimes are best possible for (d\ll n^{3/4}).  For a (d)-regular graph, (\operatorname{disc}_p^+(G)\le (\lambda_2/2)n+d), and regular-graph discrepancy is related within constants to minimum-bisection deficit and MaxCut surplus.
- **Mechanism.** Semidefinite programming, PSD/Schur-product arguments, spectral moments, and extremal strongly regular examples.
- **Quantifiers/order.** Lower bounds are uniform over every graph in the stated density ranges; sharpness uses selected graph families.  All sufficiently large orders in the lower statements.
- **Boundary.** One-sided, density-centred discrepancy is not the same minimax functional as absolute signed induced discrepancy.  The dense range close to density (1/2) has genuine obstructions (for example complete bipartite graphs), so sparse estimates should not be extrapolated.

### 5. Subgraph discrepancies in the complete graph

**Micha Christoph, Lior Gishboliner and Michael Krivelevich (2026 preprint).** [arXiv:2602.04069](https://arxiv.org/abs/2602.04069).

- **Theorem/normalization.** Given any (f:E(K_n)\to\{\pm1\}), the discrepancy of a copy (F'\cong F) is ( |\sum_{e\in E(F')}f(e)|).  For every fixed (\varepsilon>0), if an (n)-vertex (F) has no isolated vertices and (\Delta(F)\le(1-\varepsilon)n), every colouring of (K_n) contains a copy with discrepancy (\Omega(\varepsilon n)).  If (F) is (d)-regular with (d\le(1-\varepsilon)n), the bound improves to (\Omega(\sqrt{\varepsilon d}\,n)), best possible in order.  They also determine the asymptotically optimal linear constants for (K_r)-factors and for 2-factors.
- **Mechanism.** A biased-bisection dichotomy; random embedding into a biased host bisection; otherwise many switchable vertex pairs with divergent neighbourhoods give a hypergeometric (\sqrt d) gain.  Factor constants use finite coloured-Ramsey block decompositions and local extremal optimization.
- **Quantifiers/order.** Uniform over **every host signing**, with a maximum over labelled copies of a prescribed (F).  Factor statements are on their natural divisibility/admissibility orders.
- **Boundary/status.** This is copy discrepancy, not induced discrepancy on a freely chosen vertex set, and its scale is generally linear in (n).  As of the retrieval date it is an arXiv preprint; no machine-only claim was identified.

### 6. Factorization norms and an inverse theorem for MaxCut

**Igor Balla, Lianna Hambardzumyan and István Tomon (arXiv 2025; Math. Ann. 2026).** [arXiv:2506.23989](https://arxiv.org/abs/2506.23989); [DOI 10.1007/s00208-026-03355-2](https://doi.org/10.1007/s00208-026-03355-2).

- **Theorem/normalization.** An (m\times n) Boolean matrix with bounded (\gamma_2)-norm (or bounded normalized trace norm) contains an all-zero or all-one submatrix with a constant fraction of the rows and columns, the fraction depending only on the norm bound.  Their graph application says: for every fixed (K), if a graph with (m) edges satisfies (\operatorname{MaxCut}(G)\le m/2+K\sqrt m), then it contains a clique of order at least (c_K\sqrt m).
- **Mechanism.** Factorization/trace-norm estimates yield a large homogeneous rectangle; an iterative structural extraction and extremal graph arguments convert it to a clique.
- **Quantifiers/order.** Uniform inverse theorem for every Boolean matrix/graph satisfying the norm or MaxCut hypothesis; all finite orders.  Examples determining necessity are selected families.
- **Boundary.** It diagnoses structure when cut surplus is exceptionally small; it does not evaluate absolute induced discrepancy.  Passing from a signed centred matrix to a (0/1) Boolean matrix changes both the all-ones direction and the normalization.

## Cut norms, spectra, algebraic and random constructions

### 7. Approximating the cut-norm via Grothendieck's inequality

**Noga Alon and Assaf Naor (2006).** *SIAM J. Comput.* 35, 787–803. [Author PDF](https://www.cs.tau.ac.il/~nogaa/PDFS/cutnorm3.pdf); [DOI 10.1137/S0097539704441629](https://doi.org/10.1137/S0097539704441629).

- **Theorem/normalization.** For an (m\times n) real matrix,
  \[
  \|A\|_C=\max_{S,T}\left|\sum_{i\in S,j\in T}a_{ij}\right|,
  \quad
  \|A\|_{\infty\to1}=\max_{x_i,y_j\in\{\pm1\}}
       \sum_{ij}a_{ij}x_i y_j.
  \]
  Lemma 2.1 gives (\|A\|_C\le\|A\|_{\infty\to1}\le4\|A\|_C), with equality (\|A\|_{\infty\to1}=4\|A\|_C) when every row and column sum is zero.  Grothendieck's inequality bounds the vector SDP relaxation of (\|A\|_{\infty\to1}) within (K_G), and randomized rounding gives a polynomial-time constant-factor approximation.
- **Mechanism.** Sign/indicator expansion followed by Grothendieck's inequality and randomized hyperplane rounding.
- **Quantifiers/order.** Every real matrix, every finite dimension; this is an algorithmic/norm theorem, not an existence result for a special graph family.
- **Boundary.** The factors (4) and (K_G), plus polarization from bilinear to quadratic forms, destroy an exact induced-discrepancy constant.  The paper also proves hardness barriers for exact cut-norm computation.

### 8. On orthogonal matrices

**R. E. A. C. Paley (1933).** [DOI 10.1002/sapm1933121311](https://doi.org/10.1002/sapm1933121311).

- **Theorem/normalization.** If (q\equiv1\pmod4) is a prime power, the quadratic character of (\mathbb F_q) gives a symmetric conference matrix (C) of order (q+1): (C_{ii}=0), (C_{ij}\in\{\pm1\}), and (CC^T=qI).  Thus (\|C\|_{\rm op}=\sqrt q), and for every (S),
  \[
  |Q_C(S)|=\tfrac12|\mathbf1_S^TC\mathbf1_S|
       \le\tfrac12\sqrt q,|S|
       \le\tfrac12(q+1)\sqrt q.
  \]
- **Mechanism.** Finite-field quadratic characters and their exact orthogonality/Gauss-sum identities.
- **Quantifiers/order.** A **selected algebraic signing**, initially only at arithmetic orders (q+1).  Principal submatrices preserve the displayed bound.  Using the prime number theorem in the progression (1\bmod4), one may choose a prime (q\ge n-1) with (q/n\to1), restrict (C) to (n) vertices, and obtain an all-order selected-family bound ((1/2+o(1))n^{3/2}) under the normalization at the top.
- **Boundary/status.** This is an upper construction only, not a lower bound for every signing and not proof that (1/2) is optimal.  General conference/Hadamard existence is not known at all admissible orders; the all-order statement above uses nearby prime orders and restriction, not the Hadamard conjecture.  No conjectural existence assertion is used.

### 9. The eigenvalues of random symmetric matrices

**Zoltán Füredi and János Komlós (1981).** *Combinatorica* 1, 233–241. [DOI 10.1007/BF02579329](https://doi.org/10.1007/BF02579329).

- **Theorem/normalization.** In the mean-zero specialization, a symmetric matrix with independent, uniformly bounded upper-triangular entries of variance (\sigma^2) has all eigenvalues in ([-(2\sigma+o(1))\sqrt n,(2\sigma+o(1))\sqrt n]) with probability tending to one.  For a random Seidel matrix, (\|A\|_{\rm op}\le(2+o(1))\sqrt n), so simultaneously for every (S\subseteq[n]),
  \[
  |Q_A(S)|\le\tfrac12\|A\|_{\rm op}|S|
      \le(1+o(1))n^{3/2}.
  \]
- **Mechanism.** High trace moments and enumeration of closed walks in which every independent entry is repeated.
- **Quantifiers/order.** A random **selected construction ensemble**, available at every (n); the subset bound is uniform over all (2^n) subsets once the matrix is sampled.  It is not a statement that every signing has small operator norm.
- **Boundary.** Spectral domination is generally lossy for a maximum over (0/1) vectors and therefore does not identify the optimum induced-discrepancy constant.

## Limits, interpolation and universality

### 10. The thermodynamic limit in mean field spin glass models

**Francesco Guerra and Fabio L. Toninelli (2002).** *Commun. Math. Phys.* 230, 71–79. [arXiv:cond-mat/0204280](https://arxiv.org/abs/cond-mat/0204280); [DOI 10.1007/s00220-002-0699-y](https://doi.org/10.1007/s00220-002-0699-y).

- **Theorem/normalization.** For the Gaussian Sherrington–Kirkpatrick Hamiltonian (H_N(\sigma)=N^{-1/2}\sum_{i<j}g_{ij}\sigma_i\sigma_j), the quenched free energy per spin has a unique (N\to\infty) limit, in expectation and almost surely.  The same framework yields convergence of the ground-state energy (N^{-1}\max_\sigma H_N(\sigma)), equivalently an (N^{-3/2})-normalized maximum for the unnormalized quadratic sum.  The paper treats a broader family including (p)-spin models.
- **Mechanism.** Smooth interpolation between one (N)-spin system and two independent systems of sizes (N_1,N_2); the covariance/overlap calculation gives the needed sub/superadditivity, then Fekete's lemma and concentration give expectation and a.s. limits.
- **Quantifiers/order.** All integer sizes, but for a **random Gaussian ensemble** and a maximum over all spin vectors.
- **Boundary.** It proves existence of a random ground-state constant, not an outer minimum over deterministic signings.  A (\{\pm1\}^N) spin maximum is related to signed cuts/quadratic forms, but a (\{0,1\}^N) induced-set process has additional linear/total-sum terms.

### 11. Universality in Sherrington–Kirkpatrick's spin glass model

**Philippe Carmona and Yueyun Hu (2006).** *Ann. Inst. H. Poincaré Probab. Statist.* 42, 215–222. [DOI 10.1016/j.anihpb.2005.04.001](https://doi.org/10.1016/j.anihpb.2005.04.001).

- **Theorem/normalization.** Let (\xi_{ij}) be independent with (\mathbb E\xi=0), (\mathbb E\xi^2=1), and (\mathbb E|\xi|^3<\infty).  Their Theorem 1 proves the normalized free energy converges to the same limit as in the Gaussian environment, with comparison error at most (9\mathbb E|\xi|^3\beta^3/\sqrt N).  If
  \[
  S_N(\xi)=\max_{\sigma\in\{\pm1\}^N}\sum_{1\le i,j\le N}\xi_{ij}\sigma_i\sigma_j,
  \]
  Theorem 2 gives (N^{-3/2}S_N(\xi)\to e_\infty) a.s. and in mean, with normalized expectation comparison (O((1+\mathbb E|\xi|^3)N^{-1/6})).  Bernoulli disorder is included.
- **Mechanism.** Approximate integration by parts and Lindeberg-style replacement compare one disorder entry at a time with a Gaussian; martingale concentration and a zero-temperature passage handle the ground state.
- **Quantifiers/order.** All orders, for selected random disorder ensembles; after sampling, the maximum is uniform over exponentially many spin configurations.
- **Boundary.** Universality transfers the **random max** constant across entry laws.  It does not imply universality of a deterministic minimum over signing matrices or equality with induced-subset maxima.

### 12. The Parisi formula

**Michel Talagrand (2006).** *Annals of Mathematics* 163, 221–263. [DOI 10.4007/annals.2006.163.221](https://doi.org/10.4007/annals.2006.163.221).

- **Theorem/normalization.** For the Gaussian SK model at every inverse temperature (\beta), the limiting quenched free energy equals Parisi's variational functional (an infimum over distribution functions/order parameters).  Together with the zero-temperature limit, this identifies the Gaussian ground-state constant variationally rather than merely proving it exists.
- **Mechanism.** Guerra's interpolation gives the variational upper bound; Talagrand's cavity/induction and concentration analysis proves the matching lower bound.
- **Quantifiers/order.** An all-order limiting statement for the selected Gaussian random model, not for arbitrary deterministic matrices.
- **Boundary.** The variational constant belongs to the SK spin process.  Importing it into induced (0/1) discrepancy would require a proved process-level reduction, including magnetization and total-sum terms; none of the cited theorem supplies that reduction.

### 13. Nearly subadditive sequences

**Zoltán Füredi and Imre Z. Ruzsa (arXiv 2018; revised 2026).** [arXiv:1810.11723](https://arxiv.org/abs/1810.11723).

- **Theorem/normalization.** If (f\ge0) is nondecreasing, (\sum_{n\ge1}f(n)/n^2<\infty), and
  \[
  a(n+m)\le a(n)+a(m)+f(n+m)
  \]
  holds for all sufficiently large comparable pairs (n\le m\le\mu n) for some (\mu>1), then (a(n)/n) converges.  Even with (f=0), it suffices to have subadditivity only on such comparable pairs.  Conversely, if (\sum f(n)/n^2=\infty), they construct a nearly (f)-subadditive rational sequence whose slopes (a(n)/n) attain every rational number.
- **Mechanism.** Balanced binary aggregation controls accumulated errors when the series is summable; a correction transform reduces near-subadditivity to the exact form.  A tailored convex/slope construction gives maximal nonconvergence in the divergent regime.
- **Quantifiers/order.** Abstract all-integer sequence theorem, with no graph assumptions.
- **Boundary/status.** It applies only after the combinatorial parameter has been linearized and a recurrence with the stated error has actually been proved.  Monotonicity or restriction by itself is not near-subadditivity.  The cited source is an arXiv manuscript; the results are proved, not labelled conjectural or machine-only.

### 14. Limits of dense graph sequences

**László Lovász and Balázs Szegedy (2006).** *J. Combin. Theory B* 96, 933–957. [arXiv:math/0408173](https://arxiv.org/abs/math/0408173); [DOI 10.1016/j.jctb.2006.05.002](https://doi.org/10.1016/j.jctb.2006.05.002).

- **Theorem/normalization.** A dense graph sequence for which every fixed homomorphism density (t(F,G_n)) converges is represented by a symmetric measurable graphon (W:[0,1]^2\to[0,1]), and conversely every graphon occurs as such a limit.  The associated cut-metric compactness packages simultaneous control of dense rectangular statistics.
- **Mechanism.** Weak regularity, martingale/measure representation, and sampling from a graphon.
- **Quantifiers/order.** Arbitrary convergent dense sequences and all fixed test graphs; not a selected arithmetic family.
- **Boundary.** Cut distance is normalized by (n^2).  Every (O(n^{3/2})) signed fluctuation vanishes after graphon normalization, so ordinary graphon compactness cannot distinguish candidate second-order constants.  A second-order refinement would be extra input, not a consequence of this theorem.

## Decomposition, completion and compression of exponentially many witnesses

### 15. The existence of designs via iterative absorption: hypergraph (F)-designs for arbitrary (F)

**Stefan Glock, Daniela Kühn, Allan Lo and Deryk Osthus (preprint 2016; Memoirs AMS 2023).** [arXiv:1611.06827](https://arxiv.org/abs/1611.06827); [DOI 10.1090/MEMO/1406](https://doi.org/10.1090/MEMO/1406).

- **Theorem/normalization.** Fix an (r)-uniform hypergraph (F).  Let (g_i(F)) be the gcd of the (F)-degrees of its (i)-sets.  For all sufficiently large (n), the necessary divisibility conditions
  \[
  g_i(F)\mid\binom{n-i}{r-i}\qquad(0\le i<r)
  \]
  suffice for an exact decomposition of (K_n^{(r)}) into edge-disjoint copies of (F).  Their stronger host theorem covers clique-distribution-regular hypergraphs and yields quasirandom, resilience, and high-minimum-degree variants.
- **Mechanism.** Iterative absorption along a vortex, approximate/fractional decomposition, regularity boosting, and a prebuilt absorber that consumes the final divisible leave.
- **Quantifiers/order.** Uniform existence for every sufficiently large **admissible** (n); divisibility can leave infinitely many excluded orders.  The theorem selects a decomposition, not a signing.
- **Boundary.** Edge decomposition does not make discrepancies of blocks cancel for every vertex subset.  If a completion leaves (L) uncovered edges, their worst absolute contribution is at most (L); retaining a leading (n^{3/2}) constant therefore needs (L=o(n^{3/2})), stronger than the usual (o(n^2)) approximate-decomposition guarantee.  Exact absorption removes the leave only on admissible orders.

### 16. Hypergraph containers

**David Saxton and Andrew Thomason (arXiv 2012; Invent. Math. 2015).** [arXiv:1204.6595](https://arxiv.org/abs/1204.6595); [DOI 10.1007/s00222-014-0562-8](https://doi.org/10.1007/s00222-014-0562-8).

- **Theorem/normalization.** In a standard corollary of the container theorem, for fixed (r) and (\varepsilon>0), if an (r)-uniform hypergraph (H) on (N) vertices has average degree (d) and codegrees satisfying (\Delta_j(H)\le c_{r,\varepsilon}\tau^{j-1}d) for (2\le j\le r), then there is a family (\mathcal C) such that every independent set lies in some (C\in\mathcal C), every (H[C]) has at most (\varepsilon e(H)) edges, and
  \(
  \log|\mathcal C|=O_{r,\varepsilon}(N\tau\log(1/\tau)).
  \)
  The paper's intrinsic statement uses its codegree function (\delta(H,\tau)) and also supplies small fingerprints determining the containers.
- **Mechanism.** A deterministic scythe algorithm records a sparse fingerprint while high-degree vertices force the remainder into a lower-density container; iteration converts supersaturation into entropy bounds.
- **Quantifiers/order.** Uniform for every hypergraph meeting the codegree condition; all finite orders.
- **Boundary.** Containers compress exponentially many **independent or sparse sets**, not arbitrary vertex subsets with a large signed quadratic sum.  They become relevant only after a suitable witness hypergraph and supersaturation/codegree theorem are supplied.  They do not by themselves union-bound the (2^n) induced-set process.

### 17. Quick approximation to matrices and applications

**Alan Frieze and Ravi Kannan (1999).** *Combinatorica* 19, 175–220. [DOI 10.1007/s004930050052](https://doi.org/10.1007/s004930050052).

- **Theorem/normalization.** For an (m\times n) matrix (A\in[-1,1]^{m\times n}) and (\varepsilon>0), one can construct
  (D=\sum_{t=1}^{O(\varepsilon^{-2})}d_t\mathbf1_{S_t}\mathbf1_{T_t}^T)
  such that
  \[
  \|A-D\|_\square\le\varepsilon mn.
  \]
  The decomposition and approximate cut norm can be found by randomized sampling in polynomial time.
- **Mechanism.** Greedy extraction of a large discrepant rectangle; each cut-matrix subtraction decreases Frobenius energy, so only (O(\varepsilon^{-2})) iterations occur.
- **Quantifiers/order.** Every bounded matrix, every dimension; uniform control of all (2^{m+n}) rectangles.
- **Boundary.** Fixed-(\varepsilon) weak regularity operates at the dense (n^2) scale.  To resolve (n^{3/2}) one needs (\varepsilon\asymp n^{-1/2}), for which the bound permits (O(n)) pieces and no finite-dimensional compact template.  Thus the theorem explains both how exponential families can be compressed and why its standard form loses the target fluctuation scale.

## Order-transfer and nonconvergence audit

The following deductions are elementary consequences of the cited constructions/sequence theorems and make the order quantifiers explicit.

1. **Restriction from above preserves induced discrepancy exactly.** If an (N\times N) signing (A_N) is available and (I\subseteq[N]) has size (n\), then
   \[
   D_{\rm ind}(A_N[I])\le D_{\rm ind}(A_N).
   \]
   Therefore a construction on orders (N(n)\ge n) with (N(n)/n\to1) transfers an (n^{3/2}) leading upper constant to every order.  Paley conference matrices plus nearby primes (1\bmod4) are a concrete cited instance.  This is a limsup/construction statement, not convergence of the optimum.

2. **Extension from below has a sharper gap requirement.** Extending an (m)-vertex signing arbitrarily to (n) vertices gives only
   \[
   D_n\le D_m+\bigl[m(n-m)+\tbinom{n-m}{2}\bigr]
       \le D_m+n(n-m).
   \]
   To make the added term (o(n^{3/2})), one needs (n-m=o(\sqrt n)).  Merely having relative gaps (n-m=o(n)) is insufficient in this direction.

3. **Design completion has the same scale test.** An uncovered/retouched set of (L) edges can change any induced sum by at most (L).  Approximate decomposition with an (o(n^2)) leave does not automatically preserve a leading (n^{3/2}) constant; the needed estimate is (L=o(n^{3/2})), unless extra cancellation is proved.

4. **Near-subadditivity can force a limit, but only with a summable error.** After rescaling a parameter so its main term is linear, Füredi–Ruzsa gives convergence from comparable-pair near-subadditivity when (\sum f(n)/n^2<\infty).  Their converse shows that allowing a divergent error budget can permit maximal oscillation.  Bounds (c\le D_n/n^{3/2}\le C), monotonicity under restriction, or an arithmetic-subsequence limit alone do not imply convergence.

5. **Dense compactness and random universality address different quantifiers.** Graphon compactness discards every (n^{3/2}) fluctuation.  SK interpolation/universality controls a maximum over exponentially many configurations for a sampled random environment, but does not exchange that random maximum with a minimum over deterministic signings.

## Status and strongest takeaways

- The sharpest directly matched source is Astashkin–Lykov: for every weighted graph, the minimax induced discrepancy and the expected random-sign discrepancy are both equivalent to the explicit incident-ℓ2 sum.  For (K_n), this is a fully all-order (\Theta(n^{3/2})) theorem, with untracked constants.
- Paley gives a selected algebraic family with spectral constant (1/2) in the displayed induced-sum upper bound, initially at conference orders; principal restriction plus nearby primes transfers that limsup bound to every order.
- Guerra–Toninelli, Talagrand, and Carmona–Hu are the clean model for “interpolation proves a limit, a variational formula identifies it, and Lindeberg replacement proves universality,” but their quantifier is a random spin-glass maximum, not a deterministic signing minimum.
- Iterative absorption is the robust exact-completion mechanism, while Füredi–Ruzsa gives the exact summability threshold for a near-subadditive limit argument.  Both require error estimates at the (o(n^{3/2})) scale to retain a leading constant here.
- No primary source located in this retrieval proves an exact asymptotic constant, or even convergence of the normalized optimum, for the precise minimax induced-signing parameter.  None of the positive results above is flagged by its authors as machine-only.  The recent Christoph–Gishboliner–Krivelevich item and the revised Füredi–Ruzsa manuscript are preprints; general Hadamard/conference existence remains conjectural but is not used in the all-order restriction statement.
