# Sixth checkpoint: actual correlated response, minimizer geometry, and all-order Paley saturation

Date: 2026-09-06. This records continuing research, not the campaign's end.
The director reconstructed the arguments below rather than treating earlier
positive audits as proof steps.

## 1. Improved original lower bound

```math
\liminf_{n\to\infty} M_n/n^{3/2}\ge 0.4333221116640807.
```

This is an analytic asymptotic theorem with an exact rational interval
evaluation. It is not inferred from finite signings. The full nonlinear
response theorem, finite Gaussian realization, and spectral deletion are
the already reconstructed dependencies. Finite approximations precede the
matrix-order limit; the operator cutoff is removed last.

The new feasible policy is in the SAME two-dimensional actual Gaussian
frame `(V,Z)`, where `V=Ug` is the existing rich cyclic core and
`Z=(Uf(V)-A V)/nu`. The inverse `g` has norm one. Its scalar projection
`E[g|V]` is not substituted for it inside the inverse norm or a correlated
masked expectation.

The first successful partial step gave `.4320510727290484`. A subsequent
conditional-cell ascent selected a pure policy in this frame. Its 256
rational center rectangles are the only exploratory data read by the
exact certificate. The optimization routine is not a proof dependency.
The pure policy's residual variance is rigorously positive, approximately
`.09203135406363138`. Its inverse is a known linear combination of `g`
and the even scalar function `f(V)`.

The director read the full exact certificate and independently replayed:

```sh
.venv/bin/python computations/resumed_response_rich_core_birth_certificate_2026_09_06.py --policy computations/results/resumed_response_rich_core_optimized_rectangle_policy_2026_09_06.json --theta 1 --target 4333221116640807/10000000000000000 --output /home/math/quadra/tmp/resumed_director_optimized_policy_certificate_replay.json
```

The pre-tail lower sum is `.43332537570294300749494415241809429`.
The full omitted conditional Hermite tail is charged by at most
`.000003264038862254079131223520514945`. The final exact lower endpoint is

```text
.433322111664080753415812928897579346558634648033693413106996
```

The separate bound auditor reconstructed the mechanism and obtained an
identical certificate. Root checks included all support and tail moments,
actual old/new Gram products, parity of mirrored rectangles, the positive
beta-sum resolvent formula, outward interval arithmetic, and omission only
of nonnegative tiny-mass objective terms. In particular the certificate
uses `E[H g]`, not the generally false `E[H] E[g|V]`.

Proof and evidence:

- [Actual response birth](resumed_response_rich_core_response_birth_2026_09_06.md).
- [Same-frame ascent and independent audit](resumed_bound_audit_rich_core_cell_ascent_2026_09_06.md).
- [Exact certificate](../computations/results/resumed_response_rich_core_optimized_policy_certificate_2026_09_06.json).

## 2. A general certificate for correlated births

The [principal-angle theorem](resumed_director_conditional_birth_tail_theorem_2026_09_06.md)
was independently reconstructed by two researchers. If the projection of a
new orthogonal Gaussian block into the old first chaos has norm at most
`theta<1`, the conditional old inverse's total new-degree tail above `K`
has L2 norm at most `theta^(K+1)||g||_2`. Gaussian covariance SVD and
Parseval prove this uniformly in the old inverse's degree and dimension.

Weighted Jensen and the first-argument Lipschitz property of the Gaussian
absolute-value response transfer this to the objective with error at most
`|a|sqrt(EH)theta^(K+1)||g||_2`. The payment is made before dividing by
conditional mask mass. Thus the theorem has no rare-mask denominator.
There is no claim of a uniform angle or gain under indefinitely many births.

## 3. Exact-minimizer critical-window geometry

For a signing locally minimizing under single coefficient flips, switch
any oriented ground to the positive all-one state, with cap `Q` and
nonnegative local fields `ell_i`. If every oriented gap-two state had
projective distance at most `r`, the opposite orientation is excluded by
`Q>r(n-r)+1`. In a positive near-state, each of its own local fields is
at least `-1`. A flipped vertex therefore satisfies `ell_i<=2r-1`.

Consequently every near-state fixes `T={i:ell_i>2r-1}`. If
`n(2r-1)<Q`, then `sum_T ell_i>Q`, whereas a negative clique on `T`
would give `sum_T ell_i<=cut(T)<=Q`. A positive edge inside `T` is thus
present. Flipping that coefficient lowers every gap-two energy by two
and cannot lift any other oriented energy to `Q`: contradiction.

Hence every ground has another gap-two state at distance at least
`(c/2-o(1))sqrt(n)` when `Q>=(c-o(1))n^(3/2)`. No operator hypothesis is
needed. The director's first, weaker `n^(1/4)` argument is superseded by
the near-state local-field step.

Independent further geometry was checked directly. If all gap-two states
have the same orientation, vertices with identical near-state response
codes form negative cliques. Random independent class switches imply
`sum_C |C|^2<=2Q+n`, hence at least `1+log_2(n^2/(2Q+n))` near-states.
With mixed orientations the correct quantity is instead the overlap of
the code law with its orientation-twisted involution. It can vanish on
a small transversal, so the single-orientation count does not extend
automatically. These are actual-minimizer statements, not Gaussian-model
or asymptotic-near-minimizer statements.

- [Gap-two proof](resumed_bound_audit_gap_two_geometry_2026_09_06.md).
- [Orientation dichotomy](resumed_convergence_gap_two_orientation_geometry_2026_09_06.md).
- [Unconditional sparse-flip count](resumed_bound_audit_critical_window_state_count_2026_09_06.md).

The last count uses a bounded integer positive combination of abundant
edge-alignment patterns, with denominator at most
`D_k=ceil((k+1)^((k+1)/2))`. If `M_n>(2^k+2)D_k`, at least `k+1`
oriented states have gap at most `2D_k`. Its LP determinant bound,
availability of distinct edges, and protection of all unlisted states
were checked. It gives a logarithm-over-logarithm count in every fixed
`C sqrt(n)` window, without an orientation assumption. None of these
results yet controls joint discrepancy against a new incident row.

## 4. All-order Paley saturation

The director proposed and independently rechecked the proof in
[the Paley artifact](resumed_convergence_paley_squarefree_audit_2026_09_06.md).
For every odd prime power `q=1 mod 4`, the Paley core has BOTH normalized
Boolean Rayleigh endpoints tending to `+1` and `-1`, along ALL such `q`.
The core and bordered-conference half-energy caps therefore tend to `1/2`.

The essential new ingredient is a squarefree Hermite approximation to the
sign of a sum of independent cosine coordinates, uniform over characteristic.
Its frequencies have only distinct coefficients `0,+1,-1`. Select these
frequencies afresh in every field with one prescribed quadratic character,
while excluding a finite list of short additive aliases. A generic-line
reduction to the squarefree quadratic Weil bound proves existence uniformly.
Only finitely many moments are transferred; full trace-label independence
is neither true in general nor assumed. Oddness of the approximating
polynomial makes its squared error against sign continuous at zero.

The import was checked against the exact statement in Kim--Yip--Yoo,
[*Explicit constructions of Diophantine tuples over finite fields*,
Lemma 2.1](https://link.springer.com/article/10.1007/s11139-024-00888-5#Sec2).
All line-polynomial and characteristic hypotheses are verified in the proof.
The limit order is accuracy, polynomial degree, number of cosine variables,
moment cutoff, then every sufficiently large field. This is not a claim
about arbitrary signings or positive-density Paley restrictions.

## 5. Nonlocal positive construction and its exact limitation

[Arbitrary controlled Fourier-index translations](resumed_convergence_arbitrary_controlled_index_promotion_2026_09_06.md)
are genuinely realizable with common Boolean carriers and arbitrary
truth-table control functions. The root replayed all integer selector
checks. However the resulting primal operator is a direct sum of previously
available Fourier-sign phase operators. Its optimized seed norm is exactly
the average of those block norms. Its exponentially large formal index
family therefore adds no norm power by itself. Overlapping control networks
are not licensed by this construction.

## 6. A ceiling for the entire final response functional

The director independently reconstructed and replayed
[the uniform upper certificate](resumed_bound_audit_full_response_uniform_upper_2026_09_06.md).
For arbitrary feasible `F,H`, not merely fixed frames or scalar masks,
set `p=E|F|`, `h=1-p`, and `r=||P1F||_2`. Increasing the mask and variance,
then applying concavity of `x -> Psi(sqrt(x),t)` under the mask weight,
gives

```text
J(F,H) <= h Psi(r/sqrt(h),p-r^2)
        <= h Psi(m(p)/sqrt(h),p-m(p)^2),
m(p)=2 phi(Phi^(-1)(1-p/2)).
```

The second inequality uses the exact Gaussian moment body and monotonicity
in `r^2`. A 109-interval rational monotone envelope covers all `p`, including
the endpoint tails, and proves `sup J < .45`. Its exact upper enclosure is
`.449998927491931736970926031607336183`. This is a ceiling on the stated
lower-bound functional, NOT on actual signings or even all endpoint energies
of the algorithm that produced it. Those energies also contain a common
quadratic term that this certificate discards.

The [fixed-frame variational theorem](resumed_response_fixed_frame_variational_limit_2026_09_06.md)
also passes director reconstruction. Its global optimum exists; finite
symmetric cell optima increase to it with uniform explicit approximation
error. The global modulus handles zero variance; a positive-value superlevel
has a quantitative variance margin. Its KKT conditions are only restricted
stationarity, and the numerical policy iteration is not claimed globally
optimal. A fixed analytic core plus finite causal appendages still admits
an unrestricted strict escape. These two facts are compatible: every fixed
frame can be non-sharp while the entire terminal objective stays below .45.

## 7. A general nonlocal Ramsey saturation mechanism

The director and the independent bound auditor reconstructed
[the Fourier-involution theorem](resumed_convergence_fourier_involution_ramsey_saturation_2026_09_06.md).
Uniformly over growing finite abelian groups and all even Fourier
multipliers equal to +/-1 off the trivial character, there are Boolean
functions of asymptotically zero mean whose absolute Rayleigh quotient
tends to one. The zero-frequency multiplier is merely bounded by one.
Centering is important: otherwise the constant Boolean eigenvector makes
the unqualified norm statement trivial.

For large cyclic order, the finite Deuber theorem supplies a monochromatic
positive `(r,K,1)`-set. Positivity of ALL its signed expressions forces
its generators to be `K`-superincreasing. They are therefore short-alias
free when embedded in a sufficiently large cyclic group. For bounded
exponent, a large elementary abelian subgroup and the finite projective
Ramsey theorem supply an echelon basis. Every squarefree signed sum has
first nonzero coordinate +/-1, so evenness of the multiplier gives a
common sign. Even-order groups have an exact order-two Boolean character.
All finite thresholds are chosen before group order grows.

Both imports were read in the primary research paper
[Frankl--Graham--Rödl, *Iterated Combinatorial Density Theorems*,
Sections 3 and 6](https://fanchung.ucsd.edu/ron/papers/90_04_iterated.pdf).
The squarefree sign approximation is the separately audited analytic
component; no assertion of full prime-field coordinate independence is
made. Only one Rayleigh sign is guaranteed for a general multiplier.

Over additive finite fields the construction is robust to a vanishing
fraction of bad spectral frequencies: a random nonzero scalar dilation
avoids the bad set on an entire FIXED Ramsey host. It follows that EVERY
additive-Cayley hollow signing with operator norm at most
`(1+o(1))sqrt(q)` satisfies `Q(A_q)/q^(3/2)->1/2`. Parseval supplies the
required spectral concentration. This is an actual-signing theorem,
broader than the Paley family, but not a theorem for non-Cayley matrices
or Cayley matrices with a fixed spectral excess.

## 8. Continuing judgment

The original lower endpoint improves; neither the upper nor a comparable
cross-order recurrence improves. Convergence and nonconvergence remain open.
Work continues on the discarded common endpoint energy and the exact reach
of the Ramsey mechanism. The new .45 ceiling makes the former a more
discriminating target than further numerical tuning of this terminal
functional. The six-hour authorization remains active.
