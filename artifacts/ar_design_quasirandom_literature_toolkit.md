# Quasirandom signings and dense-design machinery: theorem packet

Literature checked through 16 August 2026.  The question audited here is deliberately narrow: which published or current primary theorems can supply exact symmetric hollow \(\{\pm1\}\)-matrices at dense orders, which can repair dense combinatorial structures below the \(n^{3/2}\) energy scale, and which actually control the universal action profile of \(T=A/\sqrt n\)?

## Executive verdict

The strongest immediately applicable construction fact is an elementary corollary of the Paley conference construction and the prime number theorem in arithmetic progressions: **for every sufficiently large integer \(n\)** there is an exact symmetric hollow \(\{\pm1\}\)-matrix \(A_n\) with
\[
\|A_n\|_{2\to2}\le (1+o(1))\sqrt n.
\]
It is obtained by taking an \(n\)-vertex principal compression of the Paley conference matrix associated with the least prime \(p\equiv1\pmod4\) just above \(n\).  This is asymptotically optimal in operator norm and gives the optimal-order \(O(n^{-1/2})\) dense cut discrepancy.  It is, however, only a low-norm comparison family: no cited theorem says that these compressions recover an arbitrarily selected action limit or one-profile.

The strongest repair-scale design fact is the typical-host packing theorem of Glock--Kuehn--Lo--Osthus: for a fixed graph \(F\), an appropriate dense typical host has an \(F\)-packing whose leave has bounded maximum degree, hence only \(O(n)\) uncovered edges.  This comfortably clears an \(o(n^{3/2})\) **energy-edit** budget.  Exact design theorems give zero leave at every sufficiently large *admissible* order.  Yet all of these results impose a fixed template and fixed local/lattice constraints; none controls all \(2^n\) spins, or realizes a prescribed law of \((f,Tf)\), uniformly over \(f\).

The closest primary result formulated in the right topology is Backhausz--Szegedy's random-matrix section.  It concerns full, nonsymmetric matrices with independent \(\pm n^{-1/2}\) entries.  It proves concentration in action distance and almost-sure convergence only after passing to a subsequence of any infinite order set; the paper explicitly leaves convergence along all natural numbers open.  It neither provides the symmetric/hollow model nor a target-dependent all-order recovery theorem.

| Theorem family | Orders | Quantitative strength | Constraint class | Universal action-profile verdict |
|---|---:|---:|---|---|
| Paley symmetric conference | \(q+1\), \(q\equiv1\pmod4\) prime power | exact \(C^2=qI\) | exact signs, strong spectrum | no target realization |
| Paley principal compression + PNT(AP) | every large \(n\) | \(\|A_n\|\le(1+o(1))\sqrt n\) | exact signs, optimal cut scale | no target realization |
| random / derandomized signing of \(K_n\) | every \(n\) | \((2+o(1))\sqrt n\) | exact signs | spectral only |
| Chung--Graham--Wilson | arbitrary large orders | qualitative \(o(n^2)\) discrepancy | every fixed graph density | too coarse at \(A/\sqrt n\) scale |
| Backhausz--Szegedy random matrices | subsequences | action-distance concentration | all fixed \(k\)-profiles | nonsymmetric/full; all-order left open |
| Keevash / iterative absorption | every large admissible \(n\) | exact, zero leave | fixed design/template | no spin/profile control |
| typical-host near decomposition | every large host satisfying typicality | graph leave \(O(n)\) | fixed \(F\) | energy-scale repair only |
| high-girth / conflict-free designs | every large admissible \(n\), or partial systems | exact for fixed girth; polynomial leave in partial result | bounded-size conflicts | not global spin constraints |
| spread design distributions | every large admissible \(n\) | simultaneous inclusion-probability bounds | all fixed packings \(\mathcal S\) | “spread” is not action-profile universality |

## 1. Conference and signed-matrix constructions

### 1.1 Exact conference matrices and their order obstruction

A conference matrix of order \(m\) is a zero-diagonal \(\{\pm1\}\)-off-diagonal matrix \(C\) satisfying
\[
CC^{\mathsf T}=(m-1)I.
\]
If \(C\) is symmetric, its eigenvalues are \(\pm\sqrt{m-1}\), each with multiplicity \(m/2\).  The Paley construction gives a symmetric conference matrix of order \(q+1\) whenever \(q\) is a prime power with \(q\equiv1\pmod4\).  See [Paley, *On Orthogonal Matrices* (1933)](https://doi.org/10.1002/sapm1933121311) and the modern spectral treatment [Haemers--Parsaei Majd, *Spectral symmetry in conference matrices*](https://arxiv.org/abs/2004.05829).

Exact symmetric conference existence is not an all-order theorem.  Necessarily \(m\equiv2\pmod4\), and \(m-1\) must be a sum of two squares; see [van Lint--Seidel, *Equilateral point sets in elliptic geometry*](https://pure.tue.nl/ws/files/1655945/593474.pdf).  Thus an argument requiring a conference identity at every integer order cannot follow from the known Paley family.

### 1.2 A derived all-order compression corollary

Let \(N\to\infty\), choose the least prime \(p\equiv1\pmod4\) with \(p+1\ge N\), and take an \(N\times N\) principal submatrix \(A_N\) of the Paley conference matrix of order \(p+1\).  The prime number theorem in the fixed progression \(1\pmod4\) gives \(p/N\to1\).  Principal compression gives
\[
A_N=A_N^{\mathsf T},\qquad (A_N)_{ii}=0,\qquad (A_N)_{ij}\in\{\pm1\},
\qquad \|A_N\|\le\sqrt p=(1+o(1))\sqrt N.
\]
This corollary is not stated as such in the conference-matrix papers; it is a direct combination of Paley, PNT(AP), and the variational characterization of operator norm.  It genuinely holds at every sufficiently large order, not merely on a density-one or ratio-dense subsequence.

The scale is optimal: every such signing satisfies
\[
\|A_N\|\ge \|A_N\|_{\mathrm F}/\sqrt N=\sqrt{N-1}.
\]
Consequently, for all \(S,T\subseteq[N]\),
\[
|\mathbf1_S^{\mathsf T}A_N\mathbf1_T|
 \le \|A_N\|\sqrt{|S||T|}=O(N^{3/2}),
\]
or \(O(N^{-1/2})\) after dense \(N^2\)-normalization.  Also, for \(H_A(x)=\tfrac12x^{\mathsf T}Ax\),
\[
\max_{x\in\{\pm1\}^N}H_{A_N}(x)
 \le \tfrac12N\|A_N\|=(\tfrac12+o(1))N^{3/2}.
\]
These are uniform scalar inequalities over spins.  They do **not** identify the set of empirical laws \(\{\operatorname{Law}(f,A_Nf/\sqrt N):|f|\le1\}\), and therefore are not a universal one-profile realization theorem.

### 1.3 All-order random and deterministic low-norm signings

For a hollow symmetric matrix with independent Rademacher entries above the diagonal, the classical random-matrix theorem of [Füredi--Komlós, *The eigenvalues of random symmetric matrices*](https://doi.org/10.1007/BF02579329) gives, for each fixed \(c>2\),
\[
\Pr\{\|A_n\|\le c\sqrt n\}=1-o(1).
\]
Thus exact signings with \(\|A_n\|\le(2+o(1))\sqrt n\) can be selected separately at every order.  Modern sharp nonasymptotic bounds include [Bandeira--van Handel](https://arxiv.org/abs/1408.6185) and the 2024 universality refinements of [Brailovskaya--van Handel](https://doi.org/10.1007/s00039-024-00692-9).

There is now also an algorithmic all-order statement.  Theorem 7.5 of Wang--Lau--Zhou, [*Derandomizing Matrix Concentration Inequalities from Free Probability* (2026 preprint)](https://arxiv.org/abs/2601.08111), says that if \(G\) is a \(k\)-regular graph on \(d\) vertices and \(k\gtrsim\log^4d\), then a polynomial-time algorithm finds an edge signing with
\[
\|A(x)\|\le2\sqrt k\left(1+O\left(
 \frac{\log^{2/3}d}{k^{1/6}}+
 \frac{\log^{3/4}d}{k^{1/4}}+
 \frac{\log d}{k^{1/2}}
\right)\right).
\]
Applied to \(K_d\), this is an exact symmetric hollow signing at every large \(d\) with norm \((2+o(1))\sqrt d\).  Neither the probabilistic nor deterministic theorem permits a prescribed action-limit target: their output guarantee is operator norm.

## 2. What quasirandom graph equivalences do—and do not—control

[Chung--Graham--Wilson, *Quasi-random graphs*](https://mathweb.ucsd.edu/~ronspubs/89_05_quasi_graphs.pdf), prove, for dense graphs of fixed edge density (classically \(1/2\)), the asymptotic equivalence of such properties as:

- correct counts for every fixed graph \(F\);
- the correct \(C_4\) count together with the edge count;
- uniform subset/cut discrepancy \(o(n^2)\);
- one principal adjacency eigenvalue of order \(n\) and all remaining eigenvalues \(o(n)\).

This is a qualitative graphon-scale equivalence.  It does not turn \(o(n^2)\) into the fluctuation-scale \(O(n^{3/2})\), and it deliberately forgets information visible after centering and multiplying by \(n^{-1/2}\).  Paley graphs happen to have nontrivial eigenvalues of order \(\sqrt n\), so they satisfy the stronger cut scale, but their fixed-pattern quasirandomness still does not prescribe an action profile.

The topology gap is explicit in [Backhausz--Szegedy, *Action convergence of operators and graphs*](https://arxiv.org/abs/1811.00626).  Dense graph limits send unnormalized independent signed matrices to the zero graphon, whereas the \(n^{-1/2}\)-normalized matrices have nontrivial profiles.  In their Section 11, \(H_n\) has **all entries independent** and equal to \(\pm n^{-1/2}\), so it is neither symmetric nor hollow.  Lemma 11.2 chooses deterministic \(M_n\) at every size with \(\|M_n\|\le3\) and \(d_M(M_n,H_n)\to0\) in probability.  Proposition 11.1 then says that every infinite order set contains an infinite subsequence along which \(H_n\) converges almost surely in action distance.  The authors explicitly leave open whether all natural numbers form a good sequence.  Hence this result supplies action compactness/concentration, not all-order recovery of a selected target and not the exact matrix class at issue.

## 3. Exact designs, approximate decompositions, and repair scale

### 3.1 Exact existence at every sufficiently large admissible order

For fixed integers \(q>r\ge1\) and \(\lambda\ge1\), the necessary divisibility conditions for an \((n,q,r,\lambda)\)-design are
\[
\binom{q-i}{r-i}\mid \lambda\binom{n-i}{r-i}
\qquad(0\le i<r).
\]
[Keevash, *The existence of designs*](https://arxiv.org/abs/1401.3665), proves that these conditions are sufficient for all sufficiently large \(n\).  [Keevash's 2024 short proof](https://arxiv.org/abs/2411.18291) states in particular: for every fixed \(q>r\), all sufficiently large \(K_q^r\)-divisible \(K_n^r\) have an exact \(K_q^r\)-decomposition.

[Glock--Kuehn--Lo--Osthus, *The existence of designs via iterative absorption*](https://arxiv.org/abs/1611.06827), proves the arbitrary fixed-template version.  For a fixed \(r\)-graph \(F\), put
\[
d_i(F)=\gcd\{|F(S)|:S\subseteq V(F),\ |S|=i\},
\]
where \(F(S)\) is the set of edges containing \(S\).  A host \(G\) is \((F,\lambda)\)-divisible when \(d_i(F)\mid\lambda|G(S)|\) for every \(i<r\) and every \(i\)-set \(S\).  Their theorem gives an exact \((F,\lambda)\)-design of every sufficiently large divisible complete host (and substantially more general hosts).

Thus the order conclusion is “all sufficiently large admissible orders,” not all integers.  For fixed parameters the admissibility constraints are fixed congruence conditions, so the admissible orders form bounded-gap residue classes when nonempty and are ratio-dense.  The theorem cannot remove an inadmissibility obstruction without changing the requested exact design.

### 3.2 Typical-host theorem and the bounded-leave result

The same paper gives a particularly sharp approximate tool.  An \(r\)-graph \(G\) on \(n\) vertices is \((c,h,p)\)-typical if for every family \(\mathcal A\) of at most \(h\) distinct \((r-1)\)-sets,
\[
\left|\bigcap_{S\in\mathcal A}G(S)\right|
 =(1\pm c)p^{|\mathcal A|}n.
\]
For a fixed \(r\)-graph \(F\) on \(f\) vertices, set \(q_0=2ff!\) and \(h=2^r\binom{q_0+r}{r}\).  Their exact typical-host theorem fixes \(p,c>0\) with
\[
c\le \frac{0.9(p/2)^h}{q_0^r4^{q_0}},
\]
and then supplies \(n_0,\gamma>0\) such that every \(n\ge n_0\), every \(\lambda\le\gamma n\), and every \((c,h,p)\)-typical \((F,\lambda)\)-divisible host has an exact \((F,\lambda)\)-design.

Their near-optimal packing theorem, under the analogous explicit bound
\[
c\le \frac{0.9p^h}{q_0^r4^{q_0}},
\]
gives an \(F\)-packing whose uncovered \(r\)-graph \(L\) has \(\Delta(L)\le C(F,p)\).  For \(r=2\), this means
\[
e(L)\le Cn/2=O(n).
\]
This is much stronger than the generic \(o(n^2)\) leave furnished by a nibble and is the cleanest theorem here for an \(o(n^{3/2})\) edit budget.  Its hypotheses nevertheless concern a fixed \(F\) and finitely many common-neighborhood statistics, not a target action profile.

For comparison, the foundational [Pippenger--Spencer theorem](https://doi.org/10.1016/0097-3165(89)90074-5) assumes an almost-regular fixed-uniformity hypergraph with maximum codegree \(o(D)\) relative to degree \(D\), and yields asymptotically optimal matchings/edge colorings.  In design applications this gives an \(o(n^r)\) leave.  For graph edges, \(o(n^2)\) alone need not be \(o(n^{3/2})\); a bare nibble has insufficient rate for the stated energy normalization.

### 3.3 Recent bounded-conflict and spread refinements (2024–2026)

- [Delcourt--Postle, *Proof of the High Girth Existence Conjecture via Refined Absorption*](https://arxiv.org/abs/2402.17856): for every fixed \(q>r\ge2\) and fixed girth threshold \(g\), every sufficiently large \(n\) satisfying \(\binom{q-i}{r-i}\mid\binom{n-i}{r-i}\) for all \(i<r\) has an exact \((n,q,r)\)-Steiner system of girth at least \(g\).  This is zero leave and all admissible orders, but \(g\) and the forbidden configuration sizes are fixed before \(n\).

- [Glock--Joos--Kim--Kuehn--Lichev, *Conflict-free hypergraph matchings*](https://arxiv.org/abs/2205.05564) (JLMS 2024): for fixed \(\ell,s,t\), there is \(\varepsilon>0\) and, for every sufficiently large \(m\), a partial \((m,s,t)\)-Steiner system of size
  \[
  (1-m^{-\varepsilon})\binom mt/\binom st.
  \]
  Every \(j\) chosen blocks, \(2\le j\le\ell\), span more than \((s-t)j+t\) points.  The leave is \(O(m^{t-\varepsilon})\).  When \(t=2\), the theorem does not assert \(\varepsilon>1/2\), so its printed bound alone does not certify \(o(m^{3/2})\).  Its general conflict theorem can forbid polynomially many conflict instances, but only under explicit bounded-size and bounded-codegree/dependency hypotheses.

- [Delcourt--Kelly--Postle, *Thresholds for \((n,q,2)\)-Steiner Systems via Refined Absorption*](https://arxiv.org/abs/2402.17858) (published 2026): for fixed \(q>2\), all sufficiently large \(K_q\)-divisible \(n\) admit a probability distribution on exact \(K_q\)-decompositions satisfying
  \[
  \Pr[\mathcal S\subseteq\mathcal H]
   \le \left(n^{-(q-6)/2-\beta}\right)^{|\mathcal S|}
  \]
  for every \(K_q\)-packing \(\mathcal S\), for some \(\beta(q)>0\).  Equivalently, it is \(n^{-(q-6)/2-\beta}\)-spread.  Their threshold theorem then embeds an exact design in \(\mathcal G^{(q)}(n,p)\) a.a.s. when \(p\ge n^{-(q-6)/2}\) (nontrivial as a sparse-host assertion for the relevant larger \(q\)).  The simultaneous quantifier over all packings is strong, but it bounds inclusion probabilities; it is not a statement about all spin functions or empirical action laws.

Keevash's [*The existence of designs II*](https://arxiv.org/abs/1802.05900) and [*Coloured and directed designs*](https://arxiv.org/abs/1807.05770) allow labeled complexes, colors, directions, resolutions, and other extra data under regularity and lattice/divisibility hypotheses.  These frameworks can encode a fixed finite menu of block statistics.  The template, labels, and constraint types remain fixed as \(n\to\infty\); no theorem there permits a growing family of \(2^n\) global spin constraints.

## 4. Finite-statistic realization is weaker than profile realization

[Lovasz--Szegedy, *Limits of dense graph sequences*](https://arxiv.org/abs/math/0408173), show that every graphon is realized by finite graphs; the sampled graph \(G(n,W)\) exists at every order and, for each fixed graph \(F\), satisfies \(t(F,G(n,W))\to t(F,W)\) almost surely.  For a fixed finite list of \(F\)'s, bounded-difference/U-statistic estimates give the natural \(O_{\mathbb P}(n^{-1/2})\) error (and hence deterministic realizations with constant-times \(n^{-1/2}\) error for each fixed list).

The quantifier is the limitation.  The finite list is fixed before \(n\); it does not range over all spins or all measurable test functions.  Worst-case graphon sampling also has only logarithmic general cut-norm guarantees—see the sampling theorem in [Borgs--Chayes--Lovasz--Sos--Vesztergombi](https://arxiv.org/abs/math/0702004)—whereas the \(n^{-1/2}\) cut scale for centered signs relies on special cancellation.  Most importantly, ordinary graphon convergence sends centered dense sign noise to zero, precisely while \(A/\sqrt n\) can retain a nontrivial action profile.

## 5. Error and quantifier audit

1. **Energy edits.** Flipping \(m\) undirected signs changes \(H_A(x)=\tfrac12x^{\mathsf T}Ax\) by at most \(2m\), uniformly in \(x\).  Thus \(m=o(n^{3/2})\) is enough to preserve the normalized extremal energy.  A bounded-degree design leave has \(m=O(n)\), so it clears this test.

2. **Operator/action edits.** The same \(m\) flips give \(\|\Delta A\|_{\mathrm F}=\sqrt{8m}\), hence only the generic bound \(\|\Delta A\|\le\sqrt{8m}\).  To infer \(\|\Delta A\|=o(\sqrt n)\) from support size alone requires \(m=o(n)\).  An \(O(n)\) leave gives merely \(O(\sqrt n)\), so bounded leave by itself does **not** prove vanishing action perturbation after division by \(\sqrt n\).  Cancellation or a stronger spectral repair theorem would be needed.

3. **Independent rounding.** Centered independent rounding errors in a dense matrix normally have operator norm \(\Theta(\sqrt n)\), as reflected in the random-matrix bounds above.  After the \(1/\sqrt n\) normalization this is order one, not \(o(1)\).  Entrywise unbiasedness is therefore not an action-stability theorem.

4. **Fixed versus universal tests.** Designs, high-girth systems, colored designs, and graphon sampling control fixed block types or bounded-size configurations.  Conflict-free matching may impose many instances of those fixed local types; spread distributions quantify over all packings but only through inclusion events.  None supplies the required target-dependent statement
\[
\text{for every bounded spin/test function }f,
\quad \operatorname{Law}(f,A_nf/\sqrt n)
\text{ approximates a prescribed profile element,}
\]
nor the converse Hausdorff inclusion preventing spurious profile elements.

5. **All-order versus admissible/subsequence.** Paley conference matrices themselves occupy an arithmetic subsequence; Paley compression supplies all integer orders but forfeits the conference identity.  Design existence supplies every sufficiently large admissible order, which is ratio-dense for fixed congruence data but not every integer.  Backhausz--Szegedy supplies action-convergent subsequences and explicitly does not prove all-order convergence.

## Bottom line

The literature provides three strong but separate modules:

- exact all-order low-norm signings, including asymptotically optimal Paley compressions;
- exact or \(O(n)\)-leave decomposition machinery at all large admissible orders;
- compactness and concentration for action profiles of an independent nonsymmetric random-matrix model.

No primary theorem located through 2026 composes these modules into a target-dependent recovery result for exact symmetric hollow signings at every sufficiently large order.  In particular, no result found upgrades fixed local densities, bounded conflict avoidance, spectral/cut discrepancy, or spread block probabilities to realization of a prescribed **universal directed one-profile** at the \(A/\sqrt n\) scale.  That missing quantifier-and-scale bridge is the substantive gap, not merely an order-divisibility issue.
