# Primary-theorem audit for the continuous outer-limit route

Evidence cutoff: 2026-09-07. This is a scope audit, not a theorem import.
The soft-flatness and entropy-compatible spectral-regularization results
are proved in the adjacent `decisive_independent_*` artifacts.

## 1. The normalization that a limit theorem must handle

The continuous model has d=n(n-1)/2 disorder variables and free energy

    -(alpha d)^(-1) log E_(uniform cube)
      exp[-alpha d (Q(B)/n^(3/2)+tau delta(B))].

Thus the disorder entropy and the penalty are at speed n^2, while the
inner spin Hamiltonian has n spins and ground-state scale n^(3/2).
Replacing Q/n^(3/2) by the beta-spin pressure produces a negative
partition moment whose exponent is

    -alpha d/(beta n)=-alpha(n-1)/(2beta),

not a fixed negative replica exponent. Ordinary quenched SK limits
therefore do not imply this outer limit.

## 2. Latest chronological-entropy result does not remove subsequences

Primary source: David Jekel, *Free information geometry and the model
theory of noncommutative stochastic processes*, arXiv:2604.12212v2,
6 July 2026: https://arxiv.org/html/2604.12212v2 .

This is a genuinely relevant enlargement of the variational language:
chronological formulas are closed under partial suprema/infima and
appropriate heat evolutions. The paper establishes entropy chain rules,
optimal transport properties, and a stochastic-control representation.
However, its limiting statements retain an ultrafilter:

- Section 1.2 defines the entropy by normalized log-volumes along an
  ultrafilter on matrix orders.
- Section 5.5, around equations (5.45)--(5.46), states a Varadhan-type
  principle and LDP along that ultrafilter.
- Immediately after (5.46), the author explicitly identifies independence
  of the pressure from the ultrafilter as necessary for a true LDP.
- Section 1.3 explains that convergence of the theories of matrix
  algebras themselves is not known in the relevant generality.

Consequently this result cannot be invoked to prove that our outer
pressure has an all-order limit. Also, our Boolean-vector supremum and
entrywise cube constraint would require an additional encoding with a
distinguished diagonal structure; no such encoding is supplied here.

## 3. Convex matrix-potential results have different hypotheses

Primary source: Jekel, *An elementary approach to free entropy theory
for convex potentials*, Analysis & PDE 13 (2020), 2289--2374:
https://msp.org/apde/2020/13-8/apde-v13-n8-p02-p.pdf .

Its matrix-potential hypotheses include uniform convexity/semiconcavity
and approximation of the gradients by trace-polynomial functions.
Our potential includes the concave term -tau d^(-1)||B||_2^2, and the
Boolean cap is not an ordinary spectral trace polynomial. Neither
hypothesis has been verified. The related convex Gibbs-law result
arXiv:1906.10051 has the same broad restriction and is not imported.

## 4. Orthogonally invariant spin-glass theorems do not justify a rare tilt

Primary source: Zhou Fan and Yihong Wu, *The replica-symmetric free
energy for Ising spin glasses with orthogonally invariant couplings*,
arXiv:2105.02797: https://arxiv.org/abs/2105.02797 .

The inspected theorem is a sufficiently-high-temperature limit for
couplings orthogonally invariant in law. The entrywise cube model is
not orthogonally invariant, and the n^2-speed outer tilt is not a
typical-disorder quenched expectation. The result therefore supplies
neither the needed distributional replacement nor the needed all-beta
outer free-energy limit.

## 5. Action compactness has the correct cap topology but not recovery

Primary source: Backhausz and Szegedy, *Action convergence of operators
and graphs*, Canadian Journal of Mathematics 74 (2022), 72--121:
https://arxiv.org/abs/1811.00626 . The exact earlier project audit is
`fresh_action_limit_convergence_recovery_2026_09_05.md`.

Bounded-operator action compactness is useful after our entropy-preserving
regularization. But it does not supply all-order realizability or
exponential-volume recovery. The original source only proves subsequential
action limits for normalized iid sign matrices, explicitly leaving the
all-order assertion open. The present search found no primary theorem
that closes that specific gap. Variance-penalty continuity is a separate
obligation unless occupancy is explicitly included as a mark.

## Conclusion

The continuous outer model is now rigorously reduced to bounded operator
norm with only a vanishing asymptotic variational error. A fixed-L
all-order Laplace principle remains unproved. No inspected dense-LDP,
matrix-potential, or action-compactness theorem supplies it automatically.
The 2026 chronological-entropy paper is valuable precisely because it
states the surviving ultrafilter-independence problem explicitly.
