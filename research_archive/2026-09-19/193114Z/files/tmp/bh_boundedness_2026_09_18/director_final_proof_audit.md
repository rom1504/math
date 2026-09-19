# Director closing proof audit, 2026-09-18

Status: mathematical verification notes, not a new theorem or priority claim.
The three-hour endpoint remains 07:27:50 UTC.

The director reread the lower-bound proof in counterexamples Sections 7--8,
including the explicit five-bit coefficient table, the degree-(1,3) transfer
pair, pointwise tangent identity, old-support convexity, rational denominator
margin and exact degree-eight normalization term. Both independent integer
programs certify that last term's degree; no nominal degree twelve is used.
The all-order passage fixes epsilon=2^-104 FIRST, then lets k tend to infinity;
the exponentially slow but positive normalization rate is enough. Picking
the largest k with M_k<=m is legitimate because M_k has bounded gaps. The
error comparison uses all 2^(10k) coefficients and remains bounded at q_m.

The director independently reread the symmetric-polynomial proof, the
block-affine odd-map/parity-lattice argument, the general sphere-valued junta
induction and the one-block stability proof. The junta induction uses a
nonzero degree-(m-1) outside coefficient to lower-bound each variable's
survival probability; it does not assume independent survival events. The
atom obstruction fixes m before sending Hadamard order N to infinity, and
only afterward chooses m to defeat a uniform residual contraction.

The final flat-diagonal certificate was checked directly: the two error
terms use different Cauchy--Schwarz factorizations. In the second one,
A=E^*U has A^*A=|E|^2 because U and E are commuting diagonal operators;
rho is NOT assumed to commute with the adjacency operator. The n delta
term is explicit. For fixed nonconstant seeds, V>=I>0 controls the full
counting-measure coefficient error, giving the actual-ratio conclusion.
Extra variables and simultaneously growing seeds are not covered silently.

The quadratic entropy witness's Y signs depend on the whole auxiliary
configuration, but density diagonals depend only on the independently
chosen Z coordinates. Averaging therefore gives the exact uniform diagonal.
Linear terms are handled by the fictitious-sign symmetry of the imaginary
product; they do not require assuming that an affine field is centered.

Primary-source cross-check: the director read ABKRT arXiv:2010.12629,
Section 4.3, Theorem 18 and Lemmas 20--24. Their Schur-multiplier factorization
provides the same commutator/degree ingredient; it applies to complex
multipliers as well. The flat-diagonal L2 error conclusion is separately
proved in the repository. The source's pointwise Boolean approximation
result is not imported as an L2 theorem. No novelty inference follows.

The final exact address-tangent program was additionally replayed after
the thirteen-program closing suite. PASS: address widths 1,2,3; F support
64,256,1024; H_j support 114,588,2616 respectively; H_j squared norm 6;
outside squared masses 309/128,687/256,1443/512; actual degrees 7,8,9.
The newly added exact parity restriction also passes. This is a finite
check of the scoped address obstruction, not an asymptotic BH lower bound.

Manual dependency scan: the 26 BH source programs use explicit seeds and
constants, the standard library, numpy/scipy/sympy/mpmath, or named tracked
sibling modules. They do not need hidden scratch research inputs. Runtime
package versions are preserved by the archive utility. Search programs
and floating-point diagnostics remain distinct from exact certificates.

All three researchers independently approved the closing synthesis after
scope corrections. The original signing problem remains paused. Neither
uniform BH boundedness nor divergence has been proved.
