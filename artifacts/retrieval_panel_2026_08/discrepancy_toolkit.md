# Discrepancy, vector balancing, and extreme-process toolkit

Checked through **2026-08-15**. This is a retrieval packet, not a proposed solution.

## 1. Translation dictionary and the main structural warning

Write

\[
Q_A(x)=\sum_{1\le i<j\le n}a_{ij}x_ix_j,
\qquad a_{ij},x_i\in\{-1,1\},
\]

and let (m=\binom n2). The naked problem minimizes over (A) the two-sided quantity

\[
\max_{x\in\{\pm1\}^n}|Q_A(x)|.
\]

Three very close literatures use genuinely different native objects.

| Native object | Objective | Exact relation / mismatch |
|---|---|---|
| Signed-graph switching | Minimum number (l(A)) of negative edges after switching vertices | (\max_x Q_A(x)=m-2l(A)), but (\max_x|Q_A(x)|=m-2\min\{l(A),l(-A)\}). Thus the absolute-value problem requires **simultaneous** control of a signing and its antipode; the usual covering radius is one-sided. |
| Induced-subgraph discrepancy | (g_A(S)=\sum_{\{i,j\}\subset S}a_{ij}) | If (S_x=\{i:x_i=1\}), (T=\sum_{i<j}a_{ij}), and (d_i=\sum_{j\ne i}a_{ij}), then (4g_A(S_x)=T+\sum_i d_ix_i+Q_A(x)). The linear degree term prevents an automatic sharp-constant transfer. |
| Gale--Berlekamp / bilinear switching | (r^TAc), with independent (r,c\in\{\pm1\}^n) | Independent signs on two vertex classes are more flexible than the symmetric diagonal restriction (r=c=x). |
| Matrix discrepancy | (\min_x\|Bx\|_\infty) for a fixed list of linear rows | It colors variables for finitely many prescribed linear constraints; (Q_A(x)) is quadratic and its outer minimization chooses the coefficients rather than the coloring. |
| Bernoulli/Gaussian processes | (\mathbb E\sup_{t\in T}\langle\varepsilon,t\rangle) or (\mathbb E\sup_{t\in T}\langle g,t\rangle) | These handle exponentially or infinitely many correlated tests, but for random coefficients. They do not commute the random/process expectation with the deterministic outer minimization over (A). |

The first identity is particularly important: the standard maximum-frustration/covering-radius theorem for (K_n) can be very sharp while saying little about the two-sided absolute objective.

## 2. Closest primary theorems

### 2.1 Erdős--Spencer: induced complete subgraphs have exactly the (n^{3/2}) scale

**P. Erdős and J. Spencer, “Imbalances in k-Colorations,” Networks 1 (1972), 379--385.** [DOI](https://doi.org/10.1002/net.3230010407); [primary scan](https://users.renyi.hu/~p_erdos/1971-05.pdf).

**Theorem.** For fixed (k\ge1), color every (k)-subset of ([n]) by (\pm1), and set

\[
H_k(n)=\min_g\max_{B\subseteq[n]}
\left|\sum_{W\in\binom Bk}g(W)\right|.
\]

For all sufficiently large (n), (c_kn^{(k+1)/2}\le H_k(n)\le C_kn^{(k+1)/2}). The case (k=2) is precisely the minimum, over edge signings of (K_n), of the maximum signed edge sum in an induced clique.

**Proof mechanism.** The upper bound uses an independent random signing, a binomial/normal tail, and a union bound over all (2^n) subsets. For (k=2), the lower bound first produces disjoint (B_1,B_2) with a large bipartite sum and uses

\[
g(B_1)+g(B_2)+g(B_1,B_2)=g(B_1\cup B_2).
\]

For general fixed (k), an iterated Littlewood--Offord argument produces a large block-product coefficient. The final transfer to one induced set is a finite-dimensional compactness statement: on the space of bounded-degree polynomials, coefficient norm is bounded by a constant times the supremum norm on ([0,1]^k).

**Boundary.** This proves the right order and supplies a rare-event upper bound plus a compactness transfer, but no limiting leading constant. The displayed translation formula shows why its induced-set statistic is not identical to (Q_A).

### 2.2 Brown--Spencer: the bipartite switching analogue

**Thomas A. Brown and Joel H. Spencer, “Minimization of (\pm1) Matrices Under Line Shifts,” Colloquium Mathematicum 23 (1971), 165--171.** [DOI](https://doi.org/10.4064/cm-23-1-165-171); [publisher record/free text](https://www.impan.pl/en/publishing-house/journals-and-series/colloquium-mathematicum/all/23/1).

Let

\[
G_n=\min_{A\in\{\pm1\}^{n\times n}}
\max_{r,c\in\{\pm1\}^n}|r^TAc|.
\]

**Theorem.** The foundational line-shift result gives the (n^{3/2}) order; in the now-standard normalization the asymptotic bracket is

\[
\left(\sqrt{2/\pi}+o(1)\right)n^{3/2}
\le G_n\le (1+o(1))n^{3/2}.
\]

**Proof mechanism.** For the lower bound, choose row signs randomly and then choose each column sign to make its column sum nonnegative. The expected total is (n\mathbb E|S_n|=(\sqrt{2/\pi}+o(1))n^{3/2}). For the upper bound, take an (n\times n) submatrix of a Hadamard matrix of order (N=(1+o(1))n); orthogonality and Cauchy--Schwarz give (n\sqrt N), while the availability of a near-size order is a number-theoretic size transfer.

**Boundary.** The paper gives independent row and column switches. It neither proves convergence of (G_n/n^{3/2}) nor transfers its constants to a symmetric zero-diagonal form evaluated only at (r=c).

### 2.3 Pellegrino--Raposo: current finite constants for the same board

**Daniel Pellegrino and Anselmo Raposo Jr., “Upper Bounds for the Constants of Bennett's Inequality and the Gale--Berlekamp Switching Game,” Mathematika 70 (2024), e12229.** [DOI](https://doi.org/10.1112/mtk.12229); [arXiv](https://arxiv.org/abs/2111.00445).

**Theorem 5.2.** For every positive integer (n),

\[
\frac1{\sqrt2}\le \frac{G_n}{n^{3/2}}
\le \frac{75\sqrt{17}}{289}<1.08.
\]

The paper also records the asymptotic bracket (\sqrt{2/\pi}+o(1)\le G_n/n^{3/2}\le1+o(1)).

**Proof mechanism.** Hadamard blocks give bilinear forms of norm at most (n\sqrt N); the authors combine explicit small Hadamard orders with tensor/product constructions and norm interpolation for Bennett/Kahane--Salem--Zygmund constants.

**Boundary.** This is a useful audit of constants and dimension padding, but still no limit and still a bipartite rather than symmetric quadratic object.

### 2.4 Solé--Zaslavsky: exact coding formulation, and why it is one-sided

**Patrick Solé and Thomas Zaslavsky, “A Coding Approach to Signed Graphs,” SIAM Journal on Discrete Mathematics 7 (1994), 544--553.** [DOI](https://doi.org/10.1137/S0895480189174374); [author PDF](https://people.math.binghamton.edu/zaslav/Tpapers/cas.sidma1994.pdf).

For a graph (\Gamma) with (n) vertices, (m) edges, and (c) components, let (D(\Gamma)) be the maximum, over signings, of the minimum number of negative edges left after switching.

**Lemma 2 / main identification.** (D(\Gamma)) is exactly the covering radius of the binary cutset (cocycle) code (C^*(\Gamma)).

**Theorem 1.**

\[
D(\Gamma)\ge \frac m2-\sqrt{\frac{\ln2}{2}\,m(n-c)}.
\]

For simple bipartite (\Gamma), Theorem 2 gives (D(\Gamma)\le(m-\sqrt m)/2). The paper also records the exact complete-graph value (D(K_n)=\lfloor (n-1)^2/4\rfloor).

**Proof mechanism.** Switching classes are cosets of (C^*(\Gamma)), frustration is coset minimum weight, and the lower bound is the binary sphere-covering bound with code dimension (n-c). The bipartite upper bound uses that the code contains the all-one word and its dual has minimum distance at least three.

**Boundary / no-go.** The exact relation for the naked statistic is

\[
\max_x|Q_A(x)|=m-2\min\{l(A),l(-A)\}.
\]

Covering radius maximizes (l(A)) alone. Indeed (D(K_n)=m/2-O(n)) is compatible with a small **one-sided** maximum for an extremal signing while its antipode can have a large maximum. Coding radius by itself does not control the absolute objective.

### 2.5 Christoph--Gishboliner--Krivelevich: 2026 graph-copy discrepancy

**Micha Christoph, Lior Gishboliner, and Michael Krivelevich, “Subgraph Discrepancies in the Complete Graph,” 2026.** [arXiv:2602.04069](https://arxiv.org/abs/2602.04069). **Status:** February 2026 v1 preprint; unrefereed as of the check date.

**Theorems 1.1--1.2.** There is an absolute (c>0) such that, for large (n):

- every (n)-vertex graph (F) with no isolated vertices and (\Delta(F)\le(1-\varepsilon)n) has, in every signing of (K_n), a copy of discrepancy at least (c\varepsilon n);
- if (F) is (d)-regular and (d\le(1-\varepsilon)n), some copy has discrepancy at least (c\sqrt{\varepsilon d}\,n), best possible in (d,n).

For (K_k)-factors they identify an extremal constant (\lambda_k) through a one-parameter bipartite construction and prove a ((\lambda_k-o(1))n) theorem. Every prescribed 2-factor has a copy with at least ((2/3-o(1))n) edges of one color, and (2/3) is sharp.

**Proof mechanism.** A biased bisection of the guest is randomly embedded into a biased host bisection. If all host bisections are nearly unbiased, the proof pairs vertices with large symmetric-difference neighborhoods and selectively switches pairs; hypergeometric anti-concentration accumulates (\Theta(\sqrt d)) gain per useful pair. The exact factor constants come from stability around the extremal bipartite construction.

**Boundary.** This is a genuine modern mechanism for (\lambda+o(1)) graph-copy constants, but it maximizes over embeddings of a fixed guest, not over a switching orbit of one edge signing.

## 3. Foundational linear, hereditary, and vector discrepancy

### 3.1 Spencer: entropy partial coloring

**Joel Spencer, “Six Standard Deviations Suffice,” Transactions of the AMS 289 (1985), 679--706.** [DOI](https://doi.org/10.1090/S0002-9947-1985-0784009-0).

**Core theorem.** Every set system of (n) sets on (n) elements has a (\pm1) coloring of discrepancy at most (6\sqrt n). More generally, the paper's partial-coloring framework yields the standard (O(\sqrt{n\log(2m/n)})) scale for (m\ge n) bounded linear constraints.

**Proof mechanism.** Quantize the values of many linear forms. Entropy bounds show that among exponentially many sign vectors two have the same rounded image; their difference colors a constant fraction of the coordinates while controlling every row. Repeating on geometrically shrinking uncolored sets makes the discrepancy increments summable.

**Boundary.** The entropy argument beats a union bound for a prescribed family of linear forms. It does not directly encode a quadratic supremum whose (2^n) tests depend on the coefficient signing being optimized.

### 3.2 Lovász--Spencer--Vesztergombi: hereditary and determinant transfers

**László Lovász, Joel Spencer, and Katalin Vesztergombi, “Discrepancy of Set-Systems and Matrices,” European Journal of Combinatorics 7 (1986), 151--160.** [DOI](https://doi.org/10.1016/S0195-6698(86)80041-5).

For (A\in\mathbb R^{m\times n}), define

\[
\operatorname{detlb}(A)=
\max_{k,B}|\det B|^{1/k},
\]

where (B) ranges over square submatrices.

**Theorems.** Under the usual (\pm1) discrepancy normalization,

\[
\operatorname{herdisc}(A)\ge\tfrac12\operatorname{detlb}(A),
\qquad
\operatorname{lindisc}(A)\le2\operatorname{herdisc}(A).
\]

The second inequality rounds every fractional point in a cube; equivalent (0/1) conventions move factors of two.

**Proof mechanism.** The determinant bound is a volume obstruction: too-small discrepancy boxes cannot cover the relevant parallelepiped/lattice quotient. Linear discrepancy is reduced to hereditary discrepancy by iterative binary rounding on the remaining fractional coordinates.

**Boundary.** These are robust size/restriction transfers, but determinant lower bounds can lose the sharp polylogarithmic factors exhibited below.

### 3.3 Banaszczyk: convex Gaussian vector balancing

**Wojciech Banaszczyk, “Balancing Vectors and Gaussian Measures of n-Dimensional Convex Bodies,” Random Structures & Algorithms 12 (1998), 351--360.** [DOI](https://doi.org/10.1002/(SICI)1098-2418(199807)12:4%3C351::AID-RSA3%3E3.0.CO;2-S).

**Theorem.** If (v_1,\ldots,v_N\in\mathbb R^m) satisfy (\|v_i\|_2\le1/5) and a symmetric convex body (K\subset\mathbb R^m) has standard Gaussian measure at least (1/2), then signs exist with (\sum_i\varepsilon_iv_i\in K). Equivalently, with (\|v_i\|_2\le1), the sum lies in (5K).

**Proof mechanism.** Gaussian translation/measure inequalities let the signs be selected while maintaining a positive-measure feasible translate; the final geometric argument produces a full coloring rather than accumulating partial-coloring errors.

**Consequences and boundary.** Taking (K=tB_\infty^m) gives (t=\Theta(\sqrt{\log m})), hence the Komlós/Beck--Fiala (O(\sqrt{\log m})) type bounds. The cube's Gaussian-measure threshold is itself a (\sqrt{\log m}) barrier for this black-box use. The theorem balances a given vector list; it does not choose quadratic coefficients.

### 3.4 Lovett--Meka: quantitative Edge-Walk partial coloring

**Shachar Lovett and Raghu Meka, “Constructive Discrepancy Minimization by Walking on the Edges,” SIAM Journal on Computing 44 (2015), 1573--1582.** [DOI](https://doi.org/10.1137/130929400); [arXiv](https://arxiv.org/abs/1203.5747).

**Partial-coloring theorem.** Given (v_1,\ldots,v_m\in\mathbb R^n), (x_0\in[-1,1]^n), and thresholds (c_j\ge0) with

\[
\sum_{j=1}^m e^{-c_j^2/16}\le n/16,
\]

their randomized Edge-Walk finds (x\in[-1,1]^n) with

\[
|\langle x-x_0,v_j\rangle|\le c_j\|v_j\|_2
\]

for every (j), while at least (n/2) coordinates are arbitrarily close to (\pm1). The success probability is bounded below by a constant and the algorithm is polynomial time.

**Proof mechanism.** A Brownian/Gaussian walk is constrained to the face determined by nearly tight coordinates and discrepancy facets. Martingale tails control every row, and an (\ell_2) potential shows that many coordinates hit the boundary. Geometric recursion produces a full coloring.

**Boundary.** The entropy condition prices a fixed finite list of linear constraints. Recursion can accumulate logarithms when row energies fail to shrink.

### 3.5 Bansal--Dadush--Garg--Lovett: a full subgaussian coloring distribution

**Nikhil Bansal, Daniel Dadush, Shashwat Garg, and Shachar Lovett, “The Gram--Schmidt Walk: A Cure for the Banaszczyk Blues,” Theory of Computing 15 (2019), Article 21.** [DOI](https://doi.org/10.4086/toc.2019.v015a021); [arXiv](https://arxiv.org/abs/1708.01079).

**Published Theorem 1.4.** For (\|v_i\|_2\le1) and (x_0\in[-1,1]^N), a polynomial-time randomized algorithm outputs (x\in\{\pm1\}^N), preserving already integral coordinates, such that

\[
Y=\sum_i(x_i-x_{0,i})v_i
\]

is (\sqrt{40})-subgaussian in every direction:

\[
\mathbb E e^{\langle\theta,Y\rangle}
\le e^{40\|\theta\|_2^2/2}
\quad\text{for all }\theta.
\]

**Proof mechanism.** A pivot variable is updated along a Gram--Schmidt-orthogonal direction until another variable freezes. A phase decomposition and an exponential supermartingale control good and bad update times.

**Boundary.** This gives distributional control in all linear directions and a constructive constant-factor Banaszczyk theorem. A union bound over many unrelated directions can still be expensive, and the theorem does not characterize a Boolean quadratic process.

## 4. Factorization, tensorization, and sharp surrogate gaps

### 4.1 Matoušek--Nikolov--Talwar: (\gamma_2) factorization and exact tensorization

**Jiří Matoušek, Aleksandar Nikolov, and Kunal Talwar, “Factorization Norms and Hereditary Discrepancy,” International Mathematics Research Notices 2020(3), 751--780.** [DOI](https://doi.org/10.1093/imrn/rny033); [arXiv](https://arxiv.org/abs/1408.1376).

Define

\[
\gamma_2(A)=\min_{A=BC}
\|B\|_{2\to\infty}\|C\|_{1\to2};
\]

geometrically, it is the least (\ell_\infty)-radius of a centered ellipsoid containing all columns of (A).

**Theorems.** For an (m)-row matrix,

\[
\frac{\gamma_2(A)}{C\log m}
\le \operatorname{herdisc}(A)
\le C\sqrt{\log m}\,\gamma_2(A).
\]

Most importantly for product arguments,

\[
\gamma_2(A\otimes B)=\gamma_2(A)\gamma_2(B)
\]

exactly. The norm is also transpose invariant, monotone under submatrices, satisfies a triangle inequality, has an SDP formulation, and has a nuclear-norm dual.

**Proof mechanism.** The upper bound factorizes through Euclidean space and applies Banaszczyk. The lower bound uses the dual/nuclear norm and determinant lower bounds. Exact multiplicativity turns product set systems into Kronecker products.

**Boundary.** Tensorization can establish exact exponential growth for (\gamma_2), but conversion to actual hereditary discrepancy loses logarithmic factors. It therefore does not preserve a sharp leading constant automatically. This matrix factorization norm is not Talagrand's generic-chaining functional, despite the shared symbol (\gamma_2).

### 4.2 Jiang--Reis: the near-sharp universal detLB upper bound

**Haotian Jiang and Victor Reis, “A Tighter Relation Between Hereditary Discrepancy and Determinant Lower Bound,” SOSA 2022, 308--313.** [DOI](https://doi.org/10.1137/1.9781611977066.24); [full version](https://arxiv.org/abs/2108.07945).

**Theorem.** For every (A\in\mathbb R^{m\times n}), one can algorithmically find a coloring with

\[
\|Ax\|_\infty
\le C\sqrt{\log m\,\log n}\;\operatorname{detLB}(A),
\]

and hence the same bound holds for (\operatorname{herdisc}(A)).

**Proof mechanism.** A hereditary partial-vector-discrepancy SDP is bounded by detLB. A random-walk partial coloring, analyzed with Freedman-type martingale tails, rounds a constant fraction; geometric recursion gives the extra (\sqrt{\log n}).

**Boundary.** This is the strongest general transfer in this chain, but the next result shows its logarithms are real rather than merely proof artifacts.

### 4.3 Li--Nikolov: sharp determinant and vector-relaxation obstructions

**Lily Li and Aleksandar Nikolov, “On the Gap Between Hereditary Discrepancy and the Determinant Lower Bound,” SIAM Journal on Discrete Mathematics 38 (2024), 1222--1238.** [DOI](https://doi.org/10.1137/23M1566790); [arXiv](https://arxiv.org/abs/2303.08167).

**Theorems 1--3.** For each fixed (\varepsilon\in(0,1)), (n\ge2), and

\[
n\le m\le 2^{n^{1-\varepsilon}},
\]

there is a (0/1) (m\times n) matrix with

\[
\frac{\operatorname{herdisc}(A)}{\operatorname{detLB}(A)}
\gtrsim_\varepsilon\sqrt{\log m\,\log n}.
\]

For every real (m\times n) matrix the ratio is (O(\sqrt n)). Their discrete Haar matrix (A_k), with (n=2^k), also satisfies

\[
\operatorname{hervecdisc}(A_k)
\gtrsim\sqrt{\log n}\;\operatorname{detLB}(A_k).
\]

**Proof mechanism.** Let (P_N) be the incidence matrix of the power set. Kronecker amplification gives

\[
\operatorname{detLB}(P_N\otimes A)
\le\sqrt{eN}\operatorname{detLB}(A),
\qquad
\operatorname{disc}(P_N\otimes A)
\ge(N/2)\operatorname{disc}_1(A).
\]

The Haar/tree matrix has bounded detLB but normalized (\ell_1) discrepancy (\Omega(\sqrt k)). The universal (O(\sqrt n)) bound instead goes through a volume lower bound and geometrically shrinking partial colorings.

**Boundary / no-go.** A determinant-only certificate can miss actual hereditary discrepancy by the full (\sqrt{\log m\log n}) factor over nearly the whole range of (m); it can already miss hereditary vector discrepancy by (\sqrt{\log n}). Exact tensor operations on a surrogate do not remove this gap.

## 5. Suprema over exponentially many tests

### 5.1 Bednorz--Latała: the Bernoulli theorem

**Witold Bednorz and Rafał Latała, “On the Boundedness of Bernoulli Processes,” Annals of Mathematics 180 (2014), 1167--1203.** [DOI](https://doi.org/10.4007/annals.2014.180.3.8); [arXiv](https://arxiv.org/abs/1305.4292).

For bounded (T\subset\ell_2), put

\[
b(T)=\mathbb E\sup_{t\in T}\sum_i t_i\varepsilon_i,
\qquad
g(T)=\mathbb E\sup_{t\in T}\sum_i t_ig_i.
\]

**Bernoulli theorem.** There are (T_1,T_2) with (T\subset T_1+T_2) such that

\[
\sup_{t\in T_1}\|t\|_1\le Cb(T),
\qquad
g(T_2)\le Cb(T).
\]

Together with the elementary reverse inequality, this characterizes (b(T)), up to universal constants, by an (\ell_1) part plus a Gaussian/majorizing-measure part.

**Proof mechanism.** A multiscale sequence of chopping maps isolates coordinates that pay in (\ell_1); adaptive partitions, concentration, and Bernoulli minoration leave a Gaussian-chaining remainder.

**Boundary.** This is the right theorem when a union bound over exponentially many correlated linear tests is wasteful. It controls an expectation under random Rademacher coefficients, not a deterministic outer minimum over coefficient signings, and it is linear rather than quadratic.

### 5.2 Cai--Chen--Shu--Wang--Zou: many good colorings and prefix constraints

**Dongrun Cai, Xue Chen, Wenxuan Shu, Haoyu Wang, and Guangyi Zou, “Revisit the Partial Coloring Method: Prefix Spencer and Sampling,” 2024.** [arXiv:2408.13756](https://arxiv.org/abs/2408.13756). **Status:** v1 preprint; no peer-reviewed version verified as of the check date.

**Theorem 1.1.** For (A\in\{0,1\}^{m\times n}), (n\ge m), an efficient Gaussian-measure algorithm finds a partial coloring with (\Omega(n)) integral entries and

\[
\max_{t\le n}\left\|\sum_{i\le t}A(\cdot,i)x_i\right\|_\infty=O(\sqrt m).
\]

Recursion gives a full prefix coloring with discrepancy (O(\sqrt m\log(O(n)/m))).

**Theorem 1.2 / 5.1.** There is an efficient sampler whose output always satisfies (\|Ax\|_\infty=O(\sqrt m)), while for every fixed (\epsilon\in\{\pm1\}^n),

\[
\Pr[x=\epsilon]=O(1.9^{-0.9n}).
\]

Thus its min-entropy is (\Omega(n)), not merely its Shannon entropy.

**Proof mechanism.** Highly correlated prefix constraints are handled by a small-deviation estimate for Gaussian partial-sum processes and the Šidák--Khatri inequality. For sampling, leverage scores identify a coordinate that can be forced to either sign while remaining in a safe subspace; a martingale potential proves that this can be repeated for (0.9n) coordinates.

**Boundary.** This directly addresses abundance and point probabilities of good signings, but only for a prescribed linear system. It does not turn the target's outer minimizer into a high-entropy family.

### 5.3 Bansal--Jiang: affine spectral independence beyond the union bound

**Nikhil Bansal and Haotian Jiang, “Decoupling via Affine Spectral-Independence: Beck--Fiala and Komlós Bounds Beyond Banaszczyk,” STOC 2026.** [arXiv:2508.03961](https://arxiv.org/abs/2508.03961); [official STOC 2026 accepted-papers list](https://acm-stoc.org/stoc2026/accepted-papers.html).

**Main theorems.** For an (m\times n) set-system matrix of column degree (k):

- a basic bound is (\widetilde O(\sqrt k\log^{1/4}n));
- if (k=\Omega(\log^2n)), discrepancy is (O(\sqrt k)), resolving Beck--Fiala in that regime;
- generally, discrepancy is (\widetilde O(\sqrt k+\sqrt{\log n})), where the tilde hides polylogarithms in (\log n).

For real matrices whose columns have (\ell_2)-norm at most one, the Komlós bound improves from (O(\sqrt{\log n})) to (\widetilde O(\log^{1/4}n)). All results are polynomial-time.

**Proof mechanism.** An SDP-guided discrete Brownian motion includes affine spectral-independence constraints on combinations of rows, not only ordinary covariance constraints. This decouples their joint evolution and controls how many rows have large discrepancy at every time, avoiding the rowwise union-bound loss.

**Boundary.** This is the strongest verified 2025--26 linear-discrepancy advance in the packet. Its state remains a fractional coloring and its observables remain linear rows; no quadratic indexed supremum theorem is stated.

### 5.4 De--Nadimpalli--O'Donnell--Servedio: dimension-free Gaussian supremum sparsification

**Anindya De, Shivam Nadimpalli, Ryan O'Donnell, and Rocco A. Servedio, “Sparsifying Suprema of Gaussian Processes,” STOC 2026.** [arXiv:2411.14664](https://arxiv.org/abs/2411.14664); [official STOC 2026 accepted-papers list](https://acm-stoc.org/stoc2026/accepted-papers.html).

**Theorem 1.** Let (T\subset\mathbb R^n) be bounded and (X_t=\langle g,t\rangle) the canonical Gaussian process. For every (\varepsilon>0), there are

\[
S\subset T,\qquad |S|=2^{2^{O(1/\varepsilon)}},
\]

and shifts (c_s\in\mathbb R) such that

\[
\mathbb E\left|
\sup_{t\in T}X_t-
\sup_{s\in S}(X_s+c_s)
\right|
\le\varepsilon\,
\mathbb E\sup_{t\in T}X_t.
\]

The size is independent of both (n) and (|T|).

**Proof mechanism.** Talagrand's majorizing-measure hierarchy is pruned into (2^{2^{O(1/\varepsilon)}}) clusters. One representative per cluster and an offset for its expected within-cluster supremum approximate the entire random supremum in (L_1).

**Boundary / no-go.** This is Gaussian, permits non-centered shifts, and is distributional rather than a pointwise sparsifier. The paper's Example 32 proves that if one insists on a centered proper sub-process (S\subset T), dimension-independent sparsification is impossible: even normalized coordinate maxima require (n^{1-o(1)}) representatives for small fixed error. It therefore supplies both a compression mechanism and a precise warning against a naive finite-witness reduction for Boolean/quadratic extrema.

## 6. What the literature says about limits, transfers, and extreme events

### Asymptotic limits and leading constants

- The closest classical results usually prove only an order: Erdős--Spencer gives (\Theta(n^{3/2})) for induced cliques; Brown--Spencer gives a nonmatching constant bracket for the bipartite board. Neither establishes convergence of its normalized optimum.
- Exact multiplicativity of (\gamma_2) proves product exponents and exponential growth rates for the surrogate. The polylogarithmic comparison to hereditary discrepancy is too lossy to return a leading constant.
- Erdős--Spencer's compactness step converts a nonzero block coefficient into some induced-set imbalance with a constant depending on fixed degree. Compactness gives positivity, not a computable asymptotic constant.
- The clearest recent (\lambda+o(1)) constants occur in graph-copy discrepancy. Christoph--Gishboliner--Krivelevich reduce the constant to an extremal one-parameter bipartite construction and prove stability around it. Their indexing family is different from a switching orbit.
- No source in this packet proves that the normalized two-sided symmetric quadratic optimum has a limit.

### Transfers between sizes and dimensions

- **Restriction/heredity:** taking column submatrices is built into hereditary discrepancy and keeps detLB/(\gamma_2) monotone.
- **Fractional-to-integral:** LSV binary rounding and partial-coloring recursion transfer a fractional point to signs, with a constant or a geometric sum of scale-dependent errors.
- **Near-size embeddings:** Brown--Spencer/Pellegrino--Raposo embed an (n\times n) board in a Hadamard order (N\ge n), paying (\sqrt{N/n}). This works because the bilinear norm is stable under taking a rectangular block.
- **Tensor products:** (\gamma_2) is exactly multiplicative, while detLB and actual discrepancy require inequalities and can separate by sharp logarithms.
- **No generic quadratic padding theorem:** adding or deleting vertices changes every Boolean state and introduces new degree, constant, and cross terms. None of the retrieved hereditary or Hadamard transfers proves an (o(n^{3/2})) error for the symmetric quadratic objective across arbitrary nearby sizes.

### Rare and extreme events over exponentially many signs

- A direct union bound can already be scale-correct when the coefficient signing is random: a fixed Boolean state has variance (\Theta(n^2)), and (2^n) tests force deviations of order (n^{3/2}), exactly as in Erdős--Spencer's induced-set upper bound.
- Spencer's entropy method exploits collisions among many sign vectors instead of bounding each bad event. Cai et al. strengthen the abundance viewpoint by constructing a sampler with exponentially small point masses.
- Bednorz--Latała replaces cardinality by geometry: an (\ell_1) exceptional part plus a Gaussian/chaining part characterizes a Bernoulli supremum.
- Gram--Schmidt Walk supplies a full subgaussian coloring distribution in all linear directions. Affine spectral independence goes further by controlling the joint evolution and the **number** of large rows at every time.
- Gaussian supremum sparsification compresses a random maximum only after allowing shifts; the paper proves that centered proper compression can be dimension dependent.
- All five mechanisms concern a fixed process/constraint family. None justifies exchanging their expectation or high-probability bound with the outer minimization over edge signings in the naked problem.

## 7. Closest mechanisms and hard boundaries, ranked

1. **Exact native identity:** signed-graph switching/cutset codes give the cleanest algebraic model, but their covering radius is one-sided. The need to control both (A) and (-A) is an exact, not cosmetic, obstruction.
2. **Correct-scale direct analogue:** Erdős--Spencer proves (\Theta(n^{3/2})) for induced cliques using a union-bound upper and Littlewood--Offord/compactness lower. Its degree term blocks an exact identification with (Q_A).
3. **Best constant-bearing cousin:** Brown--Spencer and Pellegrino--Raposo give explicit constants and a near-dimension Hadamard transfer for independent bipartite switches. Symmetry and the zero diagonal are outside their theorem.
4. **Best tool for exponentially many correlated tests:** the Bernoulli theorem/majorizing measures. It is an expected random linear-process theorem, not a deterministic min--max quadratic theorem.
5. **Best current beyond-union-bound linear method:** affine spectral independence (STOC 2026). Its SDP controls linear row discrepancies; no theorem in the paper covers quadratic Boolean observables.
6. **Best product surrogate:** matrix (\gamma_2) tensorizes exactly. Li--Nikolov show why sharp constants cannot be inferred from detLB/vector/factorization surrogates without a separate lossless comparison.
7. **Best compactness/compression warning:** Gaussian supremum sparsification works with offsets, while proper centered sparsification can require (n^{1-o(1)}) witnesses. Thus a dimension-free reduction of a Boolean extreme to finitely many original states needs additional structure not supplied by generic process theory.

## 8. Reliability notes

- The classical and 2014--24 journal/conference results above are peer reviewed.
- Bansal--Jiang and De--Nadimpalli--O'Donnell--Servedio are listed by the official STOC 2026 program and have public full versions.
- Cai--Chen--Shu--Wang--Zou and Christoph--Gishboliner--Krivelevich are treated only as preprints; their statements are not used as settled foundations.
- No conjectural equality or unverified leading constant is asserted for the naked quadratic problem.
