# Independent audit of `exact_sign_locking_boundary.md`

## Verdict

**PASS, with one scope clarification rather than a mathematical repair.**
The sparse identity, scale accounting, averaged locking defect, and
pointwise failure of prescribed locking are correct.  The finite verifier
passes.  The result is exactly as scoped: it kills a universal one-layer
duplicate-and-lock architecture, not interacting multilayer compilation.

## 1. Sparse identity

For fixed `x`, every edge variable has optimum

```math
\max_(y_ij=+-1)y_ij(x_i-t_ijx_j)
=|x_i-t_ijx_j|=1-t_ijx_ix_j.
```

Therefore EL.2 is pointwise exact for arbitrary `F`; no quadratic property
of `F` is used.  Its `N=k+binom(k,2)` order and `N^(3/4)` signal accounting
are also correct.

The phrase “two orientations recover the absolute value” should continue to
be read literally: the second orientation changes `(F,T)` to `(-F,-T)`.
Thus this is not a single-future compilation of an absolute response when
the old child cannot itself be sign-switched.

## 2. Averaged locking bound

For fixed `x`, free optimization of the duplicate block gives
`||R^Tx||_1`, while `z=Dx` gives `x^TRDx`; hence EL.6 is nonnegative.  Under
uniform `x`, each column field is a sum `S_k` of `k` Rademachers, and

```math
E_x ||R^Tx||_1=k E|S_k|,
\qquad
E_x x^TRDx=tr(RD).
```

Since `|tr(RD)|<=k`, the maximum is at least
`k E|S_k|-k`, and the stated Khintchine constant follows.  This verifies
EL.7--EL.10.  The pointwise adversary is also correct: for a selected column
`j`, choose every `x_i`, `i!=j`, so that
`d_jx_jR_ijx_i=-1`; then the prescribed-sign field is at most
`1-(k-1)<0` for `k>=3`.

The pin in EL.3 is essential.  The audited argument does not show that this
bad `x` remains selected after optimizing a particular flat child and query.
The draft already states this limitation correctly.

## 3. A stronger architecture-specific companion

The independent note
[`exact_disjoint_star_compiler_barrier.md`](exact_disjoint_star_compiler_barrier.md)
freezes a separate lower bound for all independent-star compilers, allowing
arbitrary cancellation among their higher Fourier levels.  If star `a` has
old-spin fan-in `d_a`, exact compilation of the complete cut shell forces

```math
\sum_a d_a^(3/2)>=k(k-1).
```

Thus `m=O(k)` independent auxiliaries require maximum fan-in
`Omega(k^(2/3))`, every such compiler uses `Omega(k^(3/2))` old--new
incidences, and the all-positive shell requires at least `k-1` auxiliaries
regardless of fan-in by a PSD Gram-completion rank argument.  Endpoint-local
stars require exactly one auxiliary per old edge.

This is strictly stronger than checking only the particular one-layer
locking map, but it remains architecture-specific: auxiliary--auxiliary
interactions can correlate selector choices and escape the pair-Fourier Gram
decomposition.  That correlated-selector class is the correct remaining
target.
