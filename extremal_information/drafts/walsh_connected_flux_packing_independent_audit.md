# Independent audit: connected Walsh flux packing

Audit target: [`walsh_connected_flux_packing.md`](walsh_connected_flux_packing.md).

Verdict: **PASS**, with the normalization and information-model limitations
already stated in Section 6 of the draft.  I found no hidden dependence of a
query on the encoded flux word and no reconstruction of the full Boolean
landscape.

## Algebraic state cube

In each four-coordinate chunk the marked triple is `(u,v,u+v)`, with `u,v`
independent and even.  Its kernel is therefore exactly `{000,111}`.  Disjoint
chunk supports make the full kernel the direct sum of these `h` kernels and
make every cross-chunk pairing zero.  The unused `(4h+1)`st coordinate makes
the characteristic-root fibre empty for every flux word.  Direct calculation
also gives zero diagonal Gram entries and, within chunk `i`, either `(0,0,0)`
or `(1,1,1)` on its three off-diagonal entries.  Thus the construction really
has `h` independent relation-cycle flux bits after all the advertised other
data are frozen.

The local triples are coordinate permutations of those in Theorem 21.21.
Such a permutation preserves the binary dot product, the characteristic
vector, the normalized Walsh matrix, and the Boolean cube.  Hence the
transport of the local bound and of the common favorable witness is valid.
With `n=q^2`, every bridge ceiling is `qn=n^(3/2)`, and the target terms have
the stated normalization

```math
M_0={9\over2},\qquad M_1={3(1+\sqrt{17})\over4}.
```

## Padding and query independence

For every Boolean pair, orthogonality gives

```math
|q x^T F_Ey|\le q\lVert x\rVert_2\lVert y\rVert_2=n^{3/2}.
```

When the target flux is zero, putting its one common `F_E`-fixed witness in
all blocks saturates every nonnegative connector, even when connectors meet
the target or contain cycles.  In the unit-flux state, the target and
connector bounds are valid pointwise and may therefore be added.  This proves
the one-sided padding lemma without assuming that the connected optimization
decomposes.  Nonnegativity is essential and is correctly included in the
hypothesis.

For query `i`, the public graph and all public coefficients are determined by
`i,h,H,gamma_h`; none depends on `sigma`.  The matrices `C_(a_ir^sigma)` are
the data at the marked ports, not query coefficients.  Thus the query is a
fixed port template evaluated on the hidden marked state.  In the positive
activation version all `3h` ports have nonzero child weight and the same
connected support in every query.

The Hamiltonian-path construction has no within-gadget connector, has exactly
`3h-1` connectors, and its union with the local triangles has maximum degree
four.  The dense construction supplies precisely all cross-gadget edges, so
adding the local triangles gives `K_(3h)`.  These claims hold also across the
three transitions in the displayed path order when `h>=2`; `h=1` is handled
separately.

## Gap and information conclusion

One non-target triangle is uniformly bounded in absolute value by
`(3/2+3)n^(3/2)=(9/2)n^(3/2)`.  Therefore the total activation perturbation is
exactly bounded by `epsilon_0 n^(3/2)` with `epsilon_0=9/200`.  Two maxima are
being compared, so the certified separation is

```math
\Delta_*
={3(5-\sqrt{17})\over4}-{9\over100}
=0.567670780786\ldots>0.
```

For any two distinct marked flux words, a fixed public query indexed by one
differing coordinate separates their scalar responses by this amount.  Thus
a deterministic summary that uniformly answers the declared family of `h`
queries to error below `Delta_* n^(3/2)/2` must be injective on the `2^h`
states.  This proves the claimed `h`-bit lower bound.  It does **not** prove an
average-case or randomized mutual-information bound without an additional
error model.

The scale qualification is genuine.  Here `m=4h+1`, one Walsh block has
`n=2^(2m)` Boolean variables, and the whole query has `N=3hn` variables.  The
separation becomes `Delta_*/(3h)^(3/2)` in units of `N^(3/2)`.  Accordingly,
the result establishes scalar visibility and additive response memory for
independent flux cycles, not an extensive total-density gap or entrywise
recovery of an arbitrary Gram matrix.

## Verification audit

I reran

```text
./.venv/bin/python extremal_information/experiments/verify_walsh_connected_flux_packing.py
```

and obtained `connected Walsh flux packing checks passed: 152`.  Inspection
confirms that it enumerates the complete relation/root/Gram state through
`h=5`, checks both graph topologies and all weights, and verifies the symbolic
sector and gap arithmetic.  It intentionally imports, rather than
re-enumerates, the Boolean local theorem 21.21; the scalable proof uses that
separately audited theorem plus orthogonal norm bounds, so this is not an
unsupported extrapolation from the finite checks.
