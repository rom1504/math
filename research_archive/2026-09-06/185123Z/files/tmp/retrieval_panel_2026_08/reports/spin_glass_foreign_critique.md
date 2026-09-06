# Phase-5 foreign critique: spin glass versus terminal coset drift

## Bottom line

The genuinely foreign mechanism is an **all-ambient-root Terwilliger localizer**: rooted triple moments with coverage inequalities imposed at every ambient word, not merely at codewords. It can encode the top-two-energy-layer incidence defining terminal drift and therefore suggests a new certificate/disproof program.

\(L_{\mathrm{drift}}\) retains an exact zero-temperature observable, but its unique-zero clause is a uniform no-metastable-disorder-traps hypothesis. Moreover, infinite square-field Paley conference roots force that zero to be \(1/2\). Thus \(L_{\mathrm{drift}}\) would prove the sharp constant \(1/2\), not merely convergence, and its tractability is correspondingly much lower than its original B-level presentation suggests.

No all-temperature pressure continuation survives: the state \((r,b)\) provably omits coset-histogram information needed by every fixed-temperature Laplace transform.

## 1. Foreign mechanism and possible new route

Gijswijt--Polak's covering SDP inserts all-root inequalities
\[
  \sum_i\lambda_i A_i(u)\geq\gamma
\]
as scalar and matrix localizers inside Terwilliger triple-orbit moment matrices, for **every ambient root** \(u\). Spin-glass replicas normally root at Gibbs configurations, average over quenched disorder, and retain overlaps; they do not impose localizers simultaneously at every deterministic disorder coset.

For a coset \(U\) represented by a signing \(a\), write
\[
 Q=Q(a),\qquad g_v=(Q-a\!\cdot\!v)/2,\qquad
 N_v=\{e:a_ev_e=-1\}.
\]
The verified one-edge formula is
\[
 b(U)=N-\left|\bigcup_{g_v\leq1}N_v\right|.                 \tag{1}
\]
Thus \(b\) is an incidence statistic among an arbitrary root, its top two augmented-cut energy layers, and its coordinate neighbors. All-root triple localizers are precisely capable of retaining such data.

The narrower target they suggest is
\[
 b(U)=0,\quad z(U)\in I
 \Longrightarrow |z(U)-c|\leq\eta_n,\qquad\eta_n\to0,       \tag{2}
\]
uniformly over roots. Since every deepest coset has \(b=0\), (2) alone yields convergence. A dual certificate for (2), or a primal pseudo-distribution violating it, is genuinely different from Guerra--Toninelli pressure gluing.

This is only a route specification. The covering SDP optimizes unrestricted cover size at prescribed radius, not the radius of the fixed augmented cut code. Schrijver and the 2026 moving-projection theorem root at code points and are packing results. Growing harmonic degree \(k=\Theta(n)\) signals the needed scale but supplies no arbitrary-root localizer.

Normal-code amalgamation also does not transfer: it needs acceptable coordinates and compatible amalgams, while cut-code dimension and projective column types grow with order. The exact Zetterberg all-syndrome proof uses special finite-field equations and constant radius.

## 2. Paley conference audit: the zero is forced to \(1/2\)

Let \(q=s^2\) run through odd square prime powers and \(n=q+1\). The symmetric Paley conference switching class contains a regular representative \(C_q\) with
\[
 C_q^2=qI,\qquad C_qx=sx
\]
for some \(x\in\{\pm1\}^n\). Hence the spectral bound is attained:
\[
 Q(C_q)=\frac12x^{\mathsf T}C_qx=\frac{ns}{2},\qquad
 z(C_q)=\frac12\sqrt{\frac q{q+1}}\longrightarrow\frac12.  \tag{3}
\]

These roots have \(b=0\). The projective Paley switching group is transitive on coordinate pairs. For a maximizing \(x\), the switched conference signing \(c_{ij}x_ix_j\) has a negative edge (it cannot be the all-positive matrix). Flipping that edge raises its displayed energy by \(2\); pair transitivity modulo switching and global complement gives the same flip behavior at every edge. Thus no edge is outward.

The first square case \(q=9,n=10\) checks this exactly: \(Q=15\), all \(45\) edge flips have \(Q=17\), and \(b=0\). This is the terminal-band dead end already found by archive verification.

Apply \(L_{\mathrm{drift}}\) to the infinite roots (3). Uniformity gives \(\beta(z(C_q))\to0\). Continuity gives \(\beta(1/2)=0\), and the unique-zero clause forces \(c=1/2\). Since deepest roots also have \(b=0\), the lemma then implies
\[
 \frac{M_n}{n^{3/2}}\longrightarrow\frac12.
\]
Therefore the numerical content of \(L_{\mathrm{drift}}\) includes the sharp universal lower theorem. It remains a strictly compressed observable, but it is not a neutral convergence reduction.

## 3. Zero-temperature preservation versus metastability

At the observable level the extreme is preserved exactly:
\[
 r(U)=\frac{N-Q(a)}2,\qquad z(U)=\frac{Q(a)}{n^{3/2}},
\]
and (1) uses maximizers and states only two energy units below them. It keeps the joint \(A,-A\) channel and introduces no annealed or positive-temperature approximation.

However, \(b=0\) is exactly one-edge stability of the disorder objective: no single coupling flip lowers \(Q\). A nondeep root with \(b=0\) is a metastable local minimum in coupling space. The unique-zero assertion says every terminal one-edge-stable disorder is asymptotically at the same energy. It therefore renames metastability unless supported by an exclusion theorem.

Equation (1) explains the obstruction: \(Q\) records only the top height, while \(b\) records the union geometry of all top-two-layer sign patterns. That union can vary macroscopically at fixed \(Q\).

## 4. Strongest failure reason and fastest test

The strongest failure is a scalable non-Paley family of terminal dead ends whose normalized energy stays separated from \(1/2\), or two same-energy root families with separated normalized drift. The exact order-ten warning already has
\[
 Q=15,\ b=0,\ z=0.474341\ldots,\qquad
 M_{10}=13,\ z_*=0.411096\ldots.
\]
One order is not an asymptotic contradiction. An infinite family with
\[
 b(U_n)=0,\qquad z(U_n)\leq\tfrac12-\delta                 \tag{4}
\]
would contradict the Paley-forced unique zero immediately. Alternatively,
\[
 z(U_n)-z(V_n)\to0,\qquad
 |b(U_n)-b(V_n)|/N\nrightarrow0                            \tag{5}
\]
would kill scalar lumpability. Exact orders six through eight already show finite same-layer drift variation.

The all-root SDP should be used disproof-first to amplify (4) or certify (5), not to fit a continuous drift curve. Because \(L_{\mathrm{drift}}\) now contains the sharp \(1/2\) lower bound plus uniform no-metastability, I assign truth confidence \(0.03\) and tractability \(0.01\).

## 5. Strongest spin-glass rescue

The closest structural theorem is Huang--Sellke's exponentially branching near-ground-state construction and quadratic lower-tail suppression for spherical Gaussian disorder. It provides substantially more landscape geometry than a Parisi value and is qualitatively capable of addressing the union in (1).

A rescue would require a deterministic Ising/Rademacher incidence analogue uniform over every terminal signing: outside \(z=1/2+o(1)\), the top-two layers must leave a positive fraction of edges uncovered. The known theorem is spherical, Gaussian, probabilistic, and not uniform under outer deterministic disorder selection. Worse, many diverse near-ground states can cover more edges in (1), reducing \(b\), so its mechanism has no automatic favorable sign. It inspires no current proof of \(L_{\mathrm{drift}}\).

## 6. All-temperature continuation after archive collision

None survives in the \((r,b)\) state. The antipodal pressure is
\[
 \widehat Z_n(A;\beta)=2^{-n}\sum_{c\in\mathcal C_n^+}
 \exp\!\left(\frac{\beta}{\sqrt n}(N-2d(a,c))\right),
\]
which depends on the full outer distribution. Verification gives actual cosets with identical \((r,b)\) and different outer distributions, so their pressures differ for generic \(\beta\).

A growing-degree arbitrary-root Terwilliger hierarchy might approximate this transform, but at zero-temperature accuracy the coding audit shows that such closure can restore the full signed coset histogram. Fixed-degree and codeword-root versions are coset-blind. Thus pressure continuation either loses the extreme or collides with full-information reconstruction.

The frozen \(L_{\mathrm{Lap}}\) implication remains logically correct but separate; the foreign packet supplies none of its uniform deterministic binary-completion hypothesis. The only surviving foreign route is a zero-temperature, Paley-aware, disproof-first terminal certificate.

## Verdict

- **Foreign mechanism:** all-ambient-root Terwilliger coverage localizers.
- **Novelty:** a real certificate/falsifier route, not an established proof.
- **Paley collision:** forces the unique drift zero and any limit to equal \(1/2\).
- **Failure:** scalable terminal traps below \(1/2\), or nonvanishing same-layer drift spread.
- **Best rescue:** an unavailable deterministic incidence analogue of Huang--Sellke.
- **All-temperature continuation:** none.
