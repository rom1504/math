# Phase-5 extremal critique of the coding packet

## Verdict

The moving-projection idea is genuinely foreign to ordinary extremal-combinatorics practice, but it does **not** presently supply an absorption, design, random-algebraic, or other cross-order theorem for the augmented cut codes.  Its verified content is same-order packing.  Moving the projection from a codeword to an arbitrary ambient coset changes the scalar two-point kernel into a weighted outer/coset enumerator; at the degree needed here this restores the full rooted state.

Fixed-\(C\) spectral regularization gives at most a selected all-order upper construction with coefficient \(C/2\), while terminal \(L_{\rm drift}\) is obstructed by scalable algebraic local maxima: square-field Paley roots have \(b=0\) and \(z\to1/2\).  Thus any unique-zero drift law either fails or already contains the much stronger conclusion that the optimum tends to \(1/2\).

## 1. Foreign experiment: move a nontrivial harmonic projection with the deep-hole root

Write \(E=\binom n2\) for the Hamming length and \(C_n^+\subseteq\mathbb F_2^E\) for the augmented cut code.  The covering window is

\[
 r(U)=\frac E2-\Theta(n^{3/2})
      =\frac E2-\Theta(E^{3/4}). \tag{1}
\]

The imported packing mechanism chooses a nontrivial Boolean harmonic module \(\mathcal E_k\), embeds it across Fourier levels \(k,\ldots,L\), moves that copy with each code point, and uses a Perron vector of the resulting block-Jacobi path.  Strong-Gelfand branching makes overlaps scalar kernels of the distance between **code points**.

The natural covering experiment is to put the moving projection at an arbitrary ambient root \(U\), localize it against \(C_n^+\), and hope that its projected mass or leading Jacobi eigenvalue depends on a small root state.  The root quantity has the form

\[
 \sum_{c\in C_n^+}K_{k,L}(d(U,c)), \tag{2}
\]

so it is a weighted outer distribution of the coset \(U+C_n^+\).  It is not determined by the codeword pair distribution.

### Test result

1. **Wrong root in the proved theorem.** Schrijver and the 2026 moving-projection theorem root at code points and assume pairwise separation.  Covering needs every ambient \(U\).  The Gijswijt--Polak localizers retain every root, but optimize the size of an unrestricted cover at prescribed radius, not the radius of this prescribed growing linear code.
2. **Growing degree restores the hidden state.** Fixed degree cannot resolve the \(E^{3/4}\) displacement.  The archive audit requires degree \(k=\Theta(n)=\Theta(\sqrt E)\).  Then \(\log\binom E k=\Theta(n\log n)\); after arbitrary-root localization, the hidden harmonic support is not \(\exp(o(n))\).  Exact closure can recover the signed coset/energy histogram.
3. **No graph-order branching theorem.** Passing from \(n\) to \(n+1\) adds an entire star of \(n\) Hamming coordinates and changes the code dimension and stabilizer.  No cited Jacobi or Terwilliger theorem intertwines the arbitrary-root modules across this change with \(o(n^{3/2})\) error.
4. **The asymptotic regime is not imported.** The strict moving-projection packing exponent is proved for fixed relative distance below \(1/2\), whereas (1) approaches \(1/2\).  No uniform edge-window limit with the needed coefficient is supplied.

Therefore the mechanism produces no actual cross-order theorem.  It is a potentially strong finite same-order SDP, not a length-transfer construction.

### Strongest possible rescue

Graham--Sloane amalgamation is the strongest genuine all-syndrome length-transfer theorem in the coding packet, but it requires normal codes with compatible acceptable coordinates and an amalgamated code identity.  No such identity is known for \(C_n^+\to C_{n+1}^+\).

From extremal combinatorics, Glock--Kühn--Lo--Osthus iterative absorption is the strongest order-realization mechanism: for a **fixed** template it can remove the leave exactly, and bounded-order repairs cost \(O(n)=o(n^{3/2})\).  It does not handle the required growing harmonic degree/template, nor does an edge design impose (2) simultaneously for every ambient root.  Astashkin--Lykov controls exponentially many induced witnesses at the correct order of magnitude, but with the wrong functional and untracked constants.  An actual rescue would need both ingredients in one new theorem: arbitrary-root harmonic localizers closed under star extension, plus a summable \(o(n^{3/2})\) boundary error.

The decisive experiment is to construct the \(n\to n+1\) intertwiner explicitly.  If its boundary-sector rank is \(\exp(\Omega(n))\), or two roots with identical proposed projected state have separated next-order values, the route stops.

## 2. Fixed-\(C\) spectral regularization and all-order realization

Interpret the spectral proposal as producing signings with

\[
 \|A_n\|_{\rm op}\le (C+o(1))\sqrt n.
\]

It gives only

\[
 Q(A_n)\le\frac12n\|A_n\|_{\rm op}
          \le\left(\frac C2+o(1)\right)n^{3/2}. \tag{3}
\]

Any fixed \(C>1\) pays a fixed leading loss.  Even \(C=1+o(1)\) proves only the selected upper frontier \(1/2\), not a comparison between arbitrary minimizers or different orders.  A rank-bounded or bounded-profile correction is no cure: an \(\Theta(\sqrt n)\) operator perturbation changes the right side of (3) by \(\Theta(n^{3/2})\).

Paley plus nearby primes does give a legitimate selected all-order realization.  Restricting an order-\(N\) conference matrix to \(n\) coordinates yields

\[
 Q(B)\le\frac n2\sqrt{N-1}; \tag{4}
\]

if \(N/n\to1\), this is \((1/2+o(1))n^{3/2}\).  But arbitrary principal restriction is not monotone for the spin functional; (4) works only because operator norm is inherited.  Extension from a lower order without cancellation needs an order gap \(o(\sqrt n)\), not merely \(o(n)\).

The strongest extremal rescue would be a **uniform spectral stability theorem** saying that every near-minimizer is within \(o(n^{3/2})\) switching norm (or \(o(n^{3/2})\) edge edits) of an all-order conference-like model, followed by an \(O(n)\)-repair absorption theorem.  No cited inverse-MaxCut, discrepancy, or design theorem has these hypotheses.  The fastest falsifier is a low-\(Q\) family with no \(o(n^{3/2})\)-close conference-like representative, or a fixed positive gap between its cap and the best fixed-\(C\) spectral bound.

## 3. Terminal \(L_{\rm drift}\) and the Paley obstruction

For a coset/root \(U\), let

\[
 z(U)=\frac{E-2r(U)}{n^{3/2}}=\frac{Q(U)}{n^{3/2}},
 \qquad b(U)=|\{e:r(U+e)=r(U)+1\}|.
\]

The exact rooted identity is

\[
 b(U)=E-\left|\bigcup_{g_v\le1}N_v\right|, \tag{5}
\]

so \(b\) is controlled by the incidence union of the top two energy layers, not by distance alone.  Finite exact censuses already show nondeep dead ends and same-layer variation; (5) explains why complete-regularity intuition is unavailable.

The scalable obstruction is stronger.  Along square-field Paley orders there are roots \(P_n\) with

\[
 b(P_n)=0,
 \qquad z(P_n)\longrightarrow\frac12. \tag{6}
\]

If \(L_{\rm drift}\) holds uniformly with \(b/E=\beta(z)+o(1)\), (6) forces \(\beta(1/2)=0\).  If \(\beta\) has a unique zero \(c\), then \(c=1/2\).  Every deepest coset also has \(b=0\), so the lemma then forces

\[
 \frac{M_n}{n^{3/2}}\longrightarrow\frac12. \tag{7}
\]

Thus \(L_{\rm drift}\) is no longer a neutral quotient proving convergence to an unknown constant.  In the presence of (6), it either fails or contains the sharp universal \(1/2\) lower theorem.  Calling \((z,b/E)\) noninjective remains correct, but noninjectivity does not make the asserted terminal law easier.

The strongest conceivable extremal rescue is a local-to-global stability theorem for (5): every terminal-band union covering all edges must either be a Paley/conference configuration with \(z=1/2+o(1)\), or be globally deepest; combined with a sharp universal discrepancy lower bound this would prove (7).  Existing Bollobás--Scott or Astashkin--Lykov bounds give only unspecified constants and do not control the top-two-layer incidence union.  Hypergraph containers could help only after a new \(n^{3/2}\)-scale supersaturation/codegree theorem; their standard fixed-\(\varepsilon\) form is too coarse.

The decisive falsifier is now exact: find any second infinite \(b=0\) family with \(z\to c\ne1/2\).  In particular, an upper construction with \(M_n\le(1/2-\delta)n^{3/2}\) along the same square-field orders would contradict the unique-zero law immediately.  Conversely, proving that no such family exists is essentially the missing sharp stability/lower-bound theorem.

## 4. Cross-domain answers

| Mechanism | Strongest failure | Strongest plausible extremal rescue | Preserves zero temperature, all orders, joint cancellation, no fixed loss? | Fastest test |
|---|---|---|---|---|
| Moving Terwilliger projection | Arbitrary-root localization is the weighted outer histogram, and no \(n\to n+1\) module intertwiner exists | Fixed-template iterative absorption plus a new arbitrary-root harmonic closure theorem | **No**: only same-order packing is proved | Compute boundary-sector rank and next-order separation for equal projected states |
| Fixed-\(C\) spectrum | \(C/2\) is a fixed leading coefficient and selected spectral examples give no universal comparison | Paley restriction for realization, plus a new uniform near-minimizer stability theorem | **No** unless \(C=1+o(1)\) and stability has \(o(n^{3/2})\) error | Exhibit a low-cap family spectrally far from every conference-like model |
| \(L_{\rm drift}\) | Paley terminal traps force the unique zero to \(1/2\), making the lemma sharp-limit-equivalent or false | New top-two-layer local-to-global stability/supersaturation theorem | **Formally yes**, but only by assuming the sharp missing theorem | Construct a second scalable \(b=0\) family with a different limiting \(z\) |

Overall: the foreign packet supplies a useful falsifier and a precise rooted observable, but no actual cross-order theorem.  The moving-projection route is packing-blind to arbitrary deep holes; fixed spectral regularization is selected and coefficient-losing; and square-field Paley terminal traps remove the claimed tractability advantage of \(L_{\rm drift}\).
