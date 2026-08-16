# Energy entropy, overlap entropy, and what they cannot remember

**Status.** Definitions frozen before archive comparison; primary-literature audit completed; exact counterexample and boundary theorems proved below. This report does **not** pursue convergence of the project quantity \(M_n\).

## 1. Answer to the naked question

There are two different answers, depending on what “identical entropy” is allowed to forget.

1. **Different normalized maxima: no, under any edge-resolving interpretation relevant here.** Exact energy counts already determine the maximum. More significantly, for a homogeneous quadratic Boolean landscape, even the strictly positive upper-tail entropy determines the limiting normalized maximum: Boolean noise around a maximizer creates exponentially many states above every strict sub-edge level. Thus the maximum cannot be hidden in a lone zero-entropy spike. Full energy-resolved pair entropy is stronger still and can expose an isolated edge even for a nonquadratic landscape.

2. **Different labeled block-coupling response: yes, even with exact finite-\(n\) equality of the energy histogram and the full global energy-energy-overlap histogram.** Put a Curie--Weiss quadratic on the left half of the coordinates in one landscape and on the right half in the other. Swapping the halves is an exact overlap-preserving bijection, so all global one- and two-replica counts agree. Yet pinning, penalizing, or constraining overlap on the fixed left block gives an order-one different zero-temperature response.

The qualifier “labeled” is essential. A global pair-overlap count determines every response that is a function only of the two energies and the **global** overlap, provided zero-entropy support is retained. It does not determine a response to a fixed coordinate block. Panchenko-style synchronization can close this gap in special random Gibbs models, but synchronization is an additional theorem, not a consequence of pair entropy for an adversarial deterministic landscape.

## 2. Primary-literature audit and the relevant lessons

The terminology is not uniform across the literature, so three objects must be kept separate: the density of configurations at an energy, the complexity (number) of metastable states or critical points, and the overlap-constrained free energy.

- Derrida's [Random Energy Model](https://doi.org/10.1103/PhysRevLett.45.79) makes the microcanonical lesson explicit: the density of levels has an extensive entropy in the bulk, and the extremal energy sits where that entropy reaches zero. Derrida's [GREM](https://doi.org/10.1051/jphyslet:01985004609040100) adds correlations and progressive hierarchical freezing. The REM/GREM comparison is the canonical warning that a one-point energy law does not record organization of states.

- Monasson's [metastable-state entropy](https://doi.org/10.1103/PhysRevLett.75.2847) counts localized metastable states, not raw Boolean configurations. Franz and Parisi's [overlap potential](https://arxiv.org/abs/cond-mat/9503167) is a free energy constrained to have prescribed overlap with a reference equilibrium configuration. Auffinger, Ben Arous, and Cerny's [Kac--Rice complexity](https://arxiv.org/abs/1003.1129) counts critical points of spherical \(p\)-spin landscapes. These are related compressions, but they are not interchangeable.

- Auffinger and Chen's [zero-temperature Parisi formula](https://arxiv.org/abs/1606.05335) identifies the mixed \(p\)-spin ground-state limit by a zero-temperature variational problem. The basic finite-volume compression behind this passage is lossless only after the inverse-temperature limit is retained: for an extensive Hamiltonian,

  \[
  m_n\le \frac{1}{\beta n}\log\sum_x e^{\beta H_n(x)}
  \le m_n+\frac{\log 2}{\beta},
  \qquad m_n=\frac1n\max_x H_n(x).
  \]

  Thus equality of the entire limiting pressure for arbitrarily large \(\beta\) forces equality of ground-state limits. A single temperature or an entropy curve with its zero-support edge deleted need not.

- Ghirlanda and Guerra's [identities](https://arxiv.org/abs/cond-mat/9807333) constrain the law of overlap arrays under asymptotic Gibbs measures. Panchenko proved that the identities force [ultrametricity](https://doi.org/10.4007/annals.2013.177.1.8). In the [multi-species SK model](https://arxiv.org/abs/1310.6679), the corresponding identities yield synchronization: every species overlap is a nondecreasing Lipschitz function of the total overlap. This is exactly the sort of hypothesis that would defeat the block counterexample below. It is not available for an arbitrary deterministic quadratic form.

- Barahona's [complexity theorem for Ising spin glasses](https://doi.org/10.1088/0305-4470/15/10/028) and the explicit study of [co-Ising graphs](https://pmc.ncbi.nlm.nih.gov/articles/PMC4103701/) are useful deterministic counterweights to the mean-field probabilistic theory. The latter exhibits nonisomorphic systems with identical classical Ising spectra but distinguishable refined responses. This is the same inverse-problem phenomenon isolated here, although the theorem below uses an even simpler anchored-coordinate construction.

## 3. Definitions, including the zero-entropy resolution

Let

\[
\Sigma_n=\{-1,+1\}^n,
\qquad
H_n(x)=\sum_{1\le i<j\le n}J^{(n)}_{ij}x_ix_j,
\]

and let \(a_n>0\) be the declared energy normalization. Put

\[
e_{H,n}(x)=\frac{H_n(x)}{a_n},
\qquad
R_n(x,y)=\frac1n\sum_{i=1}^n x_i y_i.
\]

Nothing below requires \(J_{ij}\in\{\pm1\}\). When comparison with the dense-signing project scale is desired, take \(a_n=n^{3/2}\).

For intervals \(I,I_1,I_2,Q\), define the one- and two-replica counts

\[
N^1_{H,n}(I)=\#\{x:e_{H,n}(x)\in I\},
\]

\[
N^2_{H,n}(I_1,I_2,Q)
=\#\{(x,y):e_{H,n}(x)\in I_1,\ e_{H,n}(y)\in I_2,\ R_n(x,y)\in Q\}.
\]

The upper-tail count is

\[
N^\uparrow_{H,n}(t)=\#\{x:e_{H,n}(x)\ge t\}.
\]

All logarithms in this report are natural.

### 3.1 Four resolutions that must not be conflated

**Exact finite-\(n\) data.** The counting measures \(I\mapsto N^1_{H,n}(I)\) and \((I_1,I_2,Q)\mapsto N^2_{H,n}(I_1,I_2,Q)\) are retained exactly. Exact energy equality trivially fixes the maximum.

**Support-sensitive exponential data.** For a fixed open box \(B\), use

\[
\mathsf s_H(B)=\limsup_{n\to\infty}\frac1n\log N_{H,n}(B),
\qquad \log0=-\infty.
\]

Here a singleton has rate \(0\), while an empty box has rate \(-\infty\). This is the weakest conventional exponential object that still distinguishes absence from zero entropy.

**Clipped exponential data.** Replacing the preceding expression by

\[
\overline{\mathsf s}_H(B)
=\limsup_{n\to\infty}\frac1n\log\bigl(1+N_{H,n}(B)\bigr)
\]

identifies the empty set, one state, and every \(e^{o(n)}\)-state cloud. It is legitimate only if the requested observable is known to be the closure of the positive-rate region. The quadratic noise theorem below proves that the normalized maximum is such an observable. Exact degeneracy and the last energy gap are not.

**Under-resolved summaries.** A scalar Shannon entropy of the random variable \(H_n(X)\), an entropy sampled at finitely many energy values, or the unconditional overlap entropy is not an energy-resolved microcanonical profile. For example,

\[
\#\{(x,y):R_n(x,y)=q\}
=2^n\binom{n}{(1-q)n/2}
\]

is independent of \(H_n\). Likewise, multiplying \(H_n\) by a nonzero scalar is a bijective relabeling of its exact energy values and hence preserves the scalar Shannon entropy of \(H_n(X)\), while changing the maximum. Such summaries are too coarse for the naked question.

### 3.2 Definition of “identical” used in the results

The statements below use three explicit notions.

- **Edge-identical energy entropy:** the support-sensitive capacities agree for every rational open energy interval.
- **Positive-tail-identical entropy:** the limiting rates of \(N^\uparrow_{H,n}(t)\), clipped or unclipped, agree for every continuity point \(t\), and only the locus where the rate is strictly positive is used.
- **Exact global pair identity:** \(N^2_{H,n}(I_1,I_2,Q)\) agrees for every finite-\(n\) choice of sets. The counterexample satisfies this strongest notion.

Any assertion of “identical entropy” that does not specify which of these it means is incomplete at the edge.

## 4. Why a different normalized maximum cannot hide

### 4.1 General support-sensitive no-go theorem

**Theorem 4.1 (support remembers the maximum).** Let \(H_n,G_n\) be arbitrary landscapes whose normalized maxima converge to \(m_H,m_G\). If their support-sensitive energy capacities agree on every rational open interval, then \(m_H=m_G\).

**Proof.** Suppose \(m_H>m_G\). Choose rational \(r,s\) with

\[
m_G<r<m_H<s.
\]

For all sufficiently large \(n\), an \(H_n\)-maximizer has normalized energy in \((r,s)\), so the capacity of that interval is at least \(0\). For all sufficiently large \(n\), \(G_n\) has no state in \((r,s)\), so its capacity is \(-\infty\). This contradicts equality. \(\square\)

Pair data is redundant for this theorem.

### 4.2 Quadratic noise thickening: positive entropy already reaches the edge

The preceding theorem uses the distinction \(0\ne-\infty\). Homogeneous quadratics allow a stronger statement that survives clipping.

**Theorem 4.2 (universal quadratic high-energy cloud).** Let \(f_n:\Sigma_n\to\mathbb R\) be any sequence of homogeneous quadratic forms and put \(M_n=\max_x f_n(x)>0\). For \(0<\theta<1\),

\[
\liminf_{n\to\infty}\frac1n
\log\#\{x:f_n(x)\ge\theta M_n\}
\ge
h\!\left(\frac{1-\sqrt\theta}{2}\right),
\tag{4.1}
\]

where \(h(p)=-p\log p-(1-p)\log(1-p)\). More precisely, for every

\[
0<\delta<\frac{1-\sqrt\theta}{2},
\]

a positive fraction of the Hamming sphere of radius \(\lfloor\delta n\rfloor\) around a maximizer lies above \(\theta M_n\), for all sufficiently large \(n\).

**Proof.** At each \(n\), fix a maximizer \(x^\star\), choose a uniformly random \(r\)-subset \(S\subset[n]\), and flip \(x^\star\) on \(S\). Write the resulting state as \(X=x^\star z\), where \(z_i=-1\) on \(S\) and \(z_i=1\) elsewhere. For \(i\ne j\),

\[
\mathbb E z_i z_j
=\lambda_{n,r}
:=\frac{(n-2r)^2-n}{n(n-1)}.
\]

Homogeneity of degree two gives

\[
\mathbb E f_n(X)=\lambda_{n,r}M_n.
\]

Let \(p=\mathbb P\{f_n(X)\ge\theta M_n\}\). Since \(f_n(X)\le M_n\),

\[
\lambda_{n,r}M_n
\le pM_n+(1-p)\theta M_n,
\]

and therefore

\[
p\ge\frac{\lambda_{n,r}-\theta}{1-\theta}.
\]

If \(r/n\to\delta<(1-\sqrt\theta)/2\), then \(\lambda_{n,r}\to(1-2\delta)^2>\theta\), so \(p\) is bounded below by a positive constant. The sphere has size

\[
\binom nr=\exp\{(h(\delta)+o(1))n\}.
\]

Letting \(\delta\uparrow(1-\sqrt\theta)/2\) proves (4.1). \(\square\)

The absolute-cap version is identical: set \(K_n=\max_x|f_n(x)|\), choose \(\varepsilon_n\in\{\pm1\}\) so that \(\max_x\varepsilon_nf_n(x)=K_n\), and apply the theorem to \(\varepsilon_nf_n\). Thus \(\#\{x:|f_n(x)|\ge\theta K_n\}\) obeys the same lower bound.

This theorem was independently frozen for the present question and then found already proved in the archive's entropy_energy_dichotomy.md; its use here is a new information-theoretic corollary, not a claim of a new noise lemma.

**Corollary 4.3 (positive tail entropy determines any subsequential maximum).** Let \(H_n\) be homogeneous quadratic, let \(a_n>0\), and pass to a subsequence on which

\[
\frac1{a_n}\max_xH_n(x)\longrightarrow m\ge0.
\]

Define

\[
s_H^\uparrow(t)=\liminf_n\frac1n\log N^\uparrow_{H,n}(t).
\]

Then

\[
\boxed{
m=\sup\{t:s_H^\uparrow(t)>0\}.}
\tag{4.2}
\]

The same endpoint is obtained if \(\log N\) is replaced by \(\log(1+N)\).

**Proof.** If \(t>m\), the tail is eventually empty. If \(t<m\), choose \(0<\theta<1\) so that \(\theta m>t\) when \(t\ge0\); when \(t<0\), any fixed \(\theta\in(0,1)\) works. Theorem 4.2 gives exponentially many states with energy at least \(\theta\max H_n>t a_n\) eventually. \(\square\)

Consequently, two quadratic sequences with genuinely identical positive upper-tail entropy cannot have two different limiting normalized maxima. This conclusion does not assume or prove that either maximum sequence converges; it applies on any proposed pair of convergent subsequences.

### 4.3 Full cross-energy pair entropy resolves even an isolated edge

The phrase “pair-overlap entropy” is strongest when it includes both energy coordinates, rather than only overlaps among states at the same threshold.

**Proposition 4.4 (an edge-resolving full pair profile detects a rare high
state).** Consider arbitrary landscapes with normalized energies in a common
compact interval `[-C,C]`.  Suppose the retained family of energy bins
contains a fixed bin `I` whose closure lies strictly above the smaller limiting
maximum and whose interior contains the larger limiting maximum.  For any
fixed finite overlap partition and fixed finite second-energy partition, the
two clipped energy--energy--overlap profiles cannot agree on all product bins.

**Proof.** Pair a maximizer `x*` of the larger landscape with all `2^n`
states `y`.  The fixed second-energy and overlap partitions place at least
`2^n/K` pairs in one recurring product cell, where `K` is independent of
`n`.  Its clipped pair entropy has exponent at least `log 2`.  In the smaller
landscape the first-energy bin `I` is eventually empty, so every corresponding
product cell has count zero. `square`

The edge-resolution hypothesis is necessary.  An arbitrarily coarse fixed
energy partition may put both distinct maxima in the same bin and need not
detect their separation.

This proposition does not apply if “pair entropy” means only the diagonal profile with both replicas constrained to the same rare band, a Gibbs overlap distribution at one temperature, or unconditional overlap counts. Theorem 4.2 handles the first of those omissions for quadratic upper tails.

## 5. What a full global pair profile does determine

Let

\[
\mathcal T_{H,n}
=\{(e_{H,n}(x),e_{H,n}(y),R_n(x,y)):x,y\in\Sigma_n\}
\]

be the support of the exact pair counting measure.

**Proposition 5.1 (global two-replica response is already in the support).** For constants \(b_1,b_2\) and any function \(\phi\) on the finite overlap grid,

\[
\mathcal R_{H,n}(b_1,b_2,\phi)
:=\max_{x,y}
\{b_1e_{H,n}(x)+b_2e_{H,n}(y)+\phi(R_n(x,y))\}
\]

satisfies

\[
\mathcal R_{H,n}(b_1,b_2,\phi)
=\max_{(u,v,q)\in\mathcal T_{H,n}}
\{b_1u+b_2v+\phi(q)\}.
\]

Hence exact global pair identity forces exact equality of every such response. Hausdorff convergence of the support gives the analogous asymptotic statement for continuous \(\phi\).

**Proof.** This is a change of indexing: the objective depends on \((x,y)\) only through its image in \(\mathcal T_{H,n}\). \(\square\)

This proposition draws a hard boundary. A proposed counterexample with different response to a function only of the two energies and global overlap must be exploiting a zero-support resolution loss. A block coupling depends on a different observable and is not covered.

## 6. Exact counterexample for block-coupling response

Fix \(n=2m\) with \(m\) even, and split the labeled coordinates into

\[
L=\{1,\ldots,m\},
\qquad
R=\{m+1,\ldots,2m\}.
\]

Define homogeneous quadratic forms

\[
Q_L(x)
=\frac{2}{m}\sum_{\substack{i<j\\ i,j\in L}}x_ix_j
=\frac{(\sum_{i\in L}x_i)^2-m}{m},
\]

\[
Q_R(x)
=\frac{2}{m}\sum_{\substack{i<j\\ i,j\in R}}x_ix_j
=\frac{(\sum_{i\in R}x_i)^2-m}{m}.
\]

Their maxima are \(m-1\), and their minima are \(-1\). Let \(a_n>0\) be any requested normalization and set

\[
H_n^{L}=\frac{a_n}{m-1}Q_L,
\qquad
H_n^{R}=\frac{a_n}{m-1}Q_R.
\tag{6.1}
\]

Thus both normalized maxima are exactly one. In particular, \(a_n=n^{3/2}\) puts the example at the dense-signing project's energy scale. Its nonzero coefficients are \(2a_n/[m(m-1)]=\Theta(n^{-1/2})\) and its coefficient Frobenius norm is \(\Theta(\sqrt n)\), so this is a naturally normalized SK-scale quadratic, although it is weighted and has zero cross-block coefficients rather than being a complete \(\{\pm1\}\)-signing.

### 6.1 Exact equality of global energy and overlap data

Let \(\pi\) swap \(L\) and \(R\) coordinatewise. Then

\[
H_n^R(\pi x)=H_n^L(x),
\qquad
R_n(\pi x,\pi y)=R_n(x,y).
\]

Therefore \(x\mapsto\pi x\) is a bijection of every energy fiber and \((x,y)\mapsto(\pi x,\pi y)\) is a bijection of every energy-energy-global-overlap fiber. Explicitly, for all sets \(I,I_1,I_2,Q\),

\[
N^1_{H^L,n}(I)=N^1_{H^R,n}(I),
\]

\[
\boxed{
N^2_{H^L,n}(I_1,I_2,Q)
=N^2_{H^R,n}(I_1,I_2,Q).}
\tag{6.2}
\]

This is exact finite-\(n\) equality, so it survives every entropy normalization and every near-zero convention.

### 6.2 A one-replica fixed-block profile separates them

Write

\[
u_L(x)=\frac1m\sum_{i\in L}x_i
\]

and, for feasible \(u\), define

\[
\Gamma^L_{H,n}(u)
=\frac1{a_n}\max\{H_n(x):u_L(x)=u\}.
\]

Then

\[
\Gamma^L_{H^L,n}(u)
=\frac{mu^2-1}{m-1},
\qquad
\Gamma^L_{H^R,n}(u)=1.
\tag{6.3}
\]

At \(u=0\), the limits are \(0\) and \(1\), respectively.

### 6.3 A genuine two-replica block-overlap response separates them

Define the left-block overlap

\[
R_L(x,y)=\frac1m\sum_{i\in L}x_i y_i
\]

and the constrained zero-temperature response

\[
\Theta^L_{H,n}(0)
=\frac1{2a_n}
\max\{H_n(x)+H_n(y):R_L(x,y)=0\}.
\tag{6.4}
\]

For \(H_n^R\), choose both right blocks uniform, and choose the left blocks with zero mutual overlap. This gives

\[
\Theta^L_{H^R,n}(0)=1.
\tag{6.5}
\]

For \(H_n^L\), put \(u=u_L(x)\), \(v=u_L(y)\). Global sign flips of either replica preserve its energy and preserve the constraint \(R_L=0\), so assume \(u,v\ge0\). If \(p_{--}\) is the fraction of left-block coordinates on which \(x_i=y_i=-1\), then

\[
p_{--}=\frac{1-u-v+R_L(x,y)}4
=\frac{1-u-v}{4}\ge0.
\]

Hence \(u+v\le1\), and \(u^2+v^2\le1\). Consequently

\[
Q_L(x)+Q_L(y)
=m(u^2+v^2)-2\le m-2.
\]

Equality is attained by taking \(x_L\) uniform and \(y_L\) balanced. Therefore

\[
\boxed{
\Theta^L_{H^L,n}(0)
=\frac{m-2}{2(m-1)}\longrightarrow\frac12,
\qquad
\Theta^L_{H^R,n}(0)=1.}
\tag{6.6}
\]

Thus exact equality of the full global pair histogram coexists with an asymptotic block-overlap response gap of \(1/2\).

### 6.4 A soft quadratic block perturbation also separates them

Let

\[
P_L(x)=\frac{(\sum_{i\in L}x_i)^2}{m}=Q_L(x)+1,
\qquad
\widetilde P_L=\frac{a_n}{m-1}P_L.
\]

At coupling strength one, define the absolute block-perturbed cap

\[
\mathcal B^L_{H,n}
=\frac1{a_n}\max_x|H_n(x)-\widetilde P_L(x)|.
\]

For \(H_n^L\), the expression inside the absolute value is the constant \(-a_n/(m-1)\), so \(\mathcal B^L_{H^L,n}=1/(m-1)\to0\). For \(H_n^R\), the two blocks optimize independently and

\[
\mathcal B^L_{H^R,n}=\frac{m+1}{m-1}\longrightarrow1.
\tag{6.7}
\]

The separation therefore exists for both a constrained block overlap and a soft quadratic block perturbation.

### 6.5 What the construction does and does not claim

The two landscapes are isomorphic after forgetting labels. That is not a defect for a block-response question: the fixed block \(L\) is part of the external experiment, and an isomorphism is not allowed to move the apparatus while claiming the same response. If the query is declared invariant under simultaneous relabeling of the landscape and every external block, this example is intentionally identified.

The construction does not yet give a pair of nonisomorphic dense hollow \(\{\pm1\}\)-signings with identical exact global pair histograms. That stricter realization problem is left open in Section 10.

## 7. Composition stress test

The counterexample is not destroyed by adding an independent common module.

Let \(K_p\) be any quadratic landscape on a disjoint coordinate set, with maximum \(b_p\), and form

\[
\widehat H^{L}=H_n^L\oplus K_p,
\qquad
\widehat H^{R}=H_n^R\oplus K_p.
\]

The permutation \(\pi\) swaps \(L,R\) and fixes the \(K_p\) coordinates. It preserves the total overlap

\[
R_{n+p}=\frac{n}{n+p}R_n+\frac{p}{n+p}R_p,
\]

so exact one- and two-replica equality persists under direct-sum composition.

For the two-replica constraint \(R_L=0\), the common module contributes \(2b_p\) to both responses. With normalization \(2(a_n+b_p)\), the response gap is

\[
\frac{a_n}{2(a_n+b_p)}+o(1).
\]

It remains order one whenever \(a_n/(a_n+b_p)\) stays bounded away from zero. The construction is therefore stable under composition with any comparably scaled independent environment.

There are two complementary closure facts.

1. **Equivariant couplings preserve indistinguishability.** If a common added coupling \(C(x,z)\) satisfies \(C(\pi x,z)=C(x,z)\), the same permutation proof continues to work.

2. **A labeled bridge is exactly where the missing information enters.** A coupling that sees \(L\) but not \(R\), or that uses \(R_L\) and \(R_R\) separately, is not a function of the global overlap

   \[
   R_n=\tfrac12(R_L+R_R).
   \]

   Global pair entropy cannot close the composition because it has forgotten which summand carries the energetic rigidity. In a synchronized multi-species Gibbs model one may have \(R_L=L_L(R_n)\) and \(R_R=L_R(R_n)\); absent that theorem, the decomposition of \(R_n\) is genuine hidden state.

The maximum no-go theorem also passes this stress test. Direct sums add maxima, and their upper-tail entropies combine by max-plus convolution. The quadratic noise cloud prevents a comparably scaled quadratic module from contributing an isolated maximum without also contributing positive entropy at every strict sub-edge level.

This is consistent with the archive's exact microcanonical composition results: exponential counting can compose in the bulk while still missing an endpoint unless support or an appropriate refined scale is retained. No composition argument in this report is used to infer convergence of \(M_n\).

## 8. Consequences outside the \(M_n\) problem

1. **Pinned and locally perturbed Ising optimization.** A density of states plus global overlap histogram is not a sufficient statistic for the response to a local field, a frozen subsystem, an antiferromagnetic penalty on a region, or a block-specific replica constraint. The relevant descriptor must include the joint energy/block-magnetization profile or block-resolved overlaps.

2. **Multi-species spin glasses.** The counterexample is a deterministic demonstration of why synchronization matters. If species overlaps are not proved to be functions of total overlap, a scalar Parisi-style overlap coordinate is not compositionally closed under species-specific perturbations.

3. **Graph and Hamiltonian identification.** Classical co-Ising data can fail to identify response to refined fields. For an anchored experiment, even an isomorphic relabeling can be physically distinguishable because the probe is tied to a specified vertex set. This is the labeled inverse-problem version of the co-Ising phenomenon.

4. **Coding theory.** A weight enumerator or pair-distance enumerator is invariant under coordinate permutation, but shortening, puncturing, or imposing reliability penalties on a fixed coordinate block is not. Block-resolved enumerators play the same role as the \(R_L\)-resolved entropy above.

5. **Modular robust optimization.** A summary intended to support later composition must be sufficient for the future query class, not merely for the unperturbed optimum. Global energy/overlap summaries support global two-replica couplings by Proposition 5.1; they do not support arbitrary labeled bridges.

6. **Landscape sampling and learning.** Learning the bulk microcanonical curve can identify a quadratic normalized edge through Corollary 4.3, but it cannot recover exact edge multiplicity, gaps, optimizer labels, or local susceptibility. Those require moderate-deviation or support-level data.

## 9. Comparison with the repository after freezing the candidates

The archive was consulted only after the three candidates in Section 1 were fixed.

- artifacts/entropy_energy_dichotomy.md already contains the universal Hamming noise-cloud theorem, with the sharper entropy bound used in Theorem 4.2. The present report imports it as the decisive resolution theorem for energy entropy.

- artifacts/good_signing_entropy_threshold.md makes the parallel warning on **disorder space**: \(\log(1+Z)\) at speed \(n^2\) identifies an empty class with a single switching orbit, while a refined \(n^{3/2}\log n\) scale detects Hamming thickening. That is distinct from the spin-configuration entropy studied here but confirms the need to state the zero-entropy scale.

- artifacts/microcanonical_disorder_counting_composition.md proves an exact lower-tail product theorem yet records an endpoint/changing-temperature obstruction. This matches the composition boundary in Section 7: bulk count composition is not automatically edge or query sufficiency.

- artifacts/entropic_franz_parisi_bernoulli.md derives exact two-constraint kernels depending on \(q^2\) and explains why fixed-replica information does not decide a zero-violation event. The current Proposition 5.1 says precisely what full support would decide, while the block example identifies an observable omitted by global \(q\).

- artifacts/soft_cap_overlap_obstruction.md shows that bridge reveal costs depend on the full matrix of rank-one block observables, not on a lone scalar overlap. The deterministic construction in Section 6 is the static, exact version of the same loss of block location.

- extremal_information/drafts/rate_distortion_report.md distinguishes sufficiency for a full perturbation-query experiment from retention of an unperturbed optimum. The block counterexample is an explicit zero-distortion collision for the global entropy sketch but an order-one error for a fixed-block query.

No main state or ledger file is modified by this report.

## 10. Open steps and final verdict

The main logical question is settled at the stated level:

\[
\boxed{
\begin{array}{l}
\text{edge-resolved energy entropy (or positive tail entropy for quadratics)}\\
\text{cannot coexist with different limiting normalized maxima;}\\[2mm]
\text{exact energy plus exact global pair-overlap counts can coexist}\\
\text{with different labeled block-coupling response.}
\end{array}}
\]

Three stricter construction problems remain worthwhile.

1. **Dense-sign realization.** Find two hollow \(\{\pm1\}\)-coefficient quadratic sequences, preferably nonisomorphic modulo switching and permutation, with the same full energy-energy-global-overlap counts (exactly or at a declared scale) but a separated anchored block response. Co-Ising/homometric graph constructions and finite lifts are plausible starting points. This is an open realization step, not needed for the weighted-quadratic counterexample.

2. **Intrinsic higher-replica separation.** If all coordinate labels are quotiented out, Proposition 5.1 rules out separation by a global two-replica objective when pair support is retained. The next honest target is a three-replica response for two landscapes with the same two-point distance-colored data, or an intrinsic composition that needs a species decomposition unavailable from total \(q\).

3. **Minimal augmentation.** Determine whether the collection of block-resolved pair entropies for a generating family of blocks is query-sufficient, or whether higher multi-overlaps are unavoidable. Under Ghirlanda--Guerra identities and synchronization the augmentation may collapse to the total-overlap law; for deterministic adversarial landscapes there is no such collapse theorem.

Finally, even when the normalized maximum is fixed by the closure of the positive entropy region, zero entropy remains essential for finer questions. Exact ground-state degeneracy, the last gap, Poisson extremal statistics, optimizer localization, and response to an anchored perturbation all live below the resolution of a bare \(n^{-1}\log\) count.
