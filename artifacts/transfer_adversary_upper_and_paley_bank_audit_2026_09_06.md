# Theorem-bank audit: threshold upper integration and Paley ground laws

Date: 2026-09-06, approximately 19:29 UTC. Both requested theorem
audits **PASS**, with two minor wording clarifications reported to their
owners. This audit supplies no new Paley variation or convergence claim.

## 1. Threshold-credit integration

Audited `transfer_director_threshold_cap_integration_2026_09_06.md`
and its checker against the standalone construction's defect identity
and fixed-depth certificate.

For every fixed strict margin `a'<a_*`, a sufficient fixed recursion
depth gives
`p log2 + B^r Phi_4 +4(1-sqrt(p)) <= -a'`. Passing from the defect
exponential back to the full quadratic partition multiplies by
`exp(4m^2)`. Thus the annealed logarithmic coefficient is
`4sqrt(p)-a'`, not `4-a'` or a half of either expression.
Hollowing adds at most `4|tr W_T|/k<=4m` to the partition logarithm,
which is negligible on the `m^2` scale.

With `N=mk`, `k/m->p`, physical inverse temperature `8/k`, and
`Q(A)=c N^(3/2)`, the threshold pressure lower bound divided by
`m^2` is

```math
8q^2\sqrt p\,c+p\{h(\delta)+g_*\}+o(1),
\qquad q=1-2\delta.
```

This verifies both the factor `8sqrt(p)` and the numerator `p*g_*`.
Comparison with the annealed budget, followed by its legitimate
existence extraction, yields the displayed cap formula with `a'`.
Letting `a'` increase to `a_*` gives the asserted final coefficient.

The source's fixed-depth/all-order limit order is valid: fix a strict
margin, fix sufficient recursion depth, let the construction order
increase, perform the same asymptotically dense-order principal deletion,
and finally let the margin vanish. The threshold parameters remain
fixed. Its extra `O(sqrt(N))+log(N+1)+O(1)` errors vanish after division
by `m^2` and introduce no new residue or minimizer dependence.

One sentence was flagged for clarification: the endpoint budget
`4sqrt(p)-a_*+o(1)` should not be read as available at a single fixed
depth. State it first for every `a'<a_*`, then take the already specified
ordered limit. This is a wording repair, not a gap in the final theorem.

The parameter margins are strict:
`(33/32)^2 p>1` and `64/p<(33/4)^2`. Hence normalized operator norm
at most `1/sqrt(p)+o(1)` and `beta=8/sqrt(p)+o(1)` eventually lie in
the threshold theorem's ranges; `beta>8` in the limit.

The predecessor credit is exactly `v^2/(2p)`, since
`beta^2/128=1/(2p)` at the limiting weave temperature. Subtracting
the two cap formulas gives
`sqrt(p)(g_*-g_clip)/(8q^2)>0`, as claimed.

The interval checker passed. Its entropy and square-root bounds are
outward rational bounds, and every subtraction uses the correct endpoint.
It verifies `new_upper < clipped_lower < clipped_upper < 0.499432211`
in exact arithmetic. The diagnostic improvement is approximately
`3.868301613445815e-16`; the safe displayed upper endpoint is unchanged.
Running the requested checker rewrote its canonical result JSON with
the deterministic replay result, as the checker is designed to do.

## 2. General signed-symmetry ground-law theorem

Audited Section 3 of
`transfer_fresh_cavity_and_isotropic_ground_law_2026_09_06.md`
independently of its square-field eigenvector construction.

For a hollow SIGNING `A`, average one positive absolute ground state
over signed-permutation symmetries preserving `A`. The averaged law
stays on positive exact ground states. Pair transitivity makes
`A_ij E[X_i X_j]` constant over unordered edges, because the permutation
signs cancel between these two factors. Its sum is the unordered-edge
energy `Q(A)`. Consequently, with `D=binom(n,2)`,

```math
K_+=E[XX^T]=I+\frac{Q(A)}D A.
```

There is no missing factor two: `D`, not `n(n-1)`, is the number of
terms in that energy. An anti-symmetry `T^T A T=-A` makes both oriented
extrema equal and sends this law to negative exact grounds. Orthogonality
also gives `T A T^T=-A`, so its image has second-moment matrix
`K_-=I-[Q(A)/D]A`. The equal mixture is exactly isotropic.

These matrices are raw second moments, as explicitly defined in the
source; no zero-mean hypothesis or independence of spins is required.
The assertion is for signings, not arbitrary real edge weights.

## 3. All-prime-power Paley symmetries

Let the field order be any odd prime power congruent to one modulo
four. Translations preserve the finite difference character and fix
the border. For inversion, `chi(-1)=1` gives the stated character
identity for two nonzero finite points. The remaining pairs also check:

- zero/nonzero uses `chi(-z)=chi(z)`;
- infinity/nonzero uses `chi(-1/z)chi(z)=1`;
- zero/infinity keeps the border entry one.

Thus the source's signed inversion preserves the full conference matrix,
not just its finite block. Multiplication by a nonsquare flips every
finite off-diagonal entry; the extra switch at infinity flips the border
entries too, yielding the required anti-symmetry.

Translations and inversion are pair-transitive: translate one finite
point to zero, invert it to infinity, and translate the other finite
point to zero. Pairs already containing infinity require just the last
translation. All translations are available over the prime-power field;
the proof is not limited to the prime-field numerical cases.

This establishes the isotropic exact-ground law for every stated Paley
conference without requiring a Boolean spectral eigenvector or knowledge
of its exact cap.

## 4. Cavity cost and square-field law

For an appended sign row `b`, maximizing the new spin gives the exact
identity `max_z |H_A(x)+z b.x|=|H_A(x)|+|b.x|`. Under an isotropic
ground law, `E(b.x)^2=n`, so some exact ground pays at least `sqrt(n)`.
Since Paley conference order `n` is even, all dot products are even
integers; the improvement is `2 ceil(sqrt(n)/2)`. Combining this with
`Q(C)<=n sqrt(n-1)/2` gives the source's uniform derivative-scale gap
`(1/4-o(1))sqrt(n)`.

The more general slack/covariance criterion also checks: square the
pointwise bound `|b.x|<=Delta+g(x)`, average, and use
`b^T K b>=n lambda_min(K)`. Nonnegativity of `Delta` follows by testing
an exact ground, so taking the square root introduces no sign issue.

The independent square-field law is valid. Elements of the subfield's
multiplicative group are squares in the quadratic extension. Character
sums are `r-1` within a fibre and `-1` between fibres, giving the claimed
Boolean eigenvectors. Balanced fibre signs have pair expectation `-1/r`.
A square finite difference enters the subfield under a uniform square
multiplier with probability `2/(r+1)`; a nonsquare difference never does.
These facts give `I+C/r` exactly, with construction multiplicities
included. The nonsquare anti-map gives `I-C/r` and hence isotropy.

For `n=r^2+1`, odd `r`, parity rounds the cavity payment to `r+1`,
giving the displayed excess `r/4+1` over the leading derivative.
The finite replay passed all integer identities for `r=3,5,7` and the
complete prime-field ground laws for field orders `5,13,17`. It correctly
enumerates both old-spin and appended-row cubes modulo their separate
global signs. The exact finite extension values are `E(C_6)=9` and
`E(C_10)=19`; the larger square-field outputs are lower bounds, not
claimed optima.

One wording clarification was sent to the owner: the block statement
in Section 6 concerns an arbitrary SIGN bridge. For a general real
bridge its isotropic second-moment lower bound is the Frobenius norm,
not automatically `sqrt(nh)`. No block theorem is inferred here.

## 5. Original convergence scope

The source's proposed power-saving insertion bound for suitably chosen
EXACT minimizers would indeed imply convergence by a summable-error
recurrence. The newly proved Paley obstruction does not refute that
target: there is no large-order exact-minimizer hypothesis for the Paley
parents, and the square-field family's limiting coefficient is one half.
Finite exact-minimizer coincidences cannot replace this missing
asymptotic hypothesis. Both audited sources keep that distinction.

No temporary research files were produced. Apart from the integration
checker's prescribed canonical result replay, this combined audit is
the only new file from the task.
