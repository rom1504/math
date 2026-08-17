# Audit: finite-speed lexicographic count carrier

Scope: independent proof audit of
`drafts/multispeed_lexicographic_hypograph_scout.md`.  I checked the
valuation algebra, recovery construction, finite composition statement,
finite-type benchmark, and both growing-fibre examples.

## Verdict

**PASS AFTER SCOPE REPAIRS.**  The fixed-alphabet theorem is correct and the
two asymptotic calculations are correct.  The principal repair is conceptual:
the two-grid example falsifies ordinary Hausdorff convergence of the
*unscaled coefficient graph* as a sufficient state.  It does not falsify an
unspecified genuinely multiscale topology, since a topology resolving the
leading coefficient to `o(a_L/a_1)` already distinguishes the two meshes.
The Vandermonde example is the cleaner intrinsic demonstration that a
pointwise coefficient roof omits saddle multiplicity.

## 1. Theorem MS.1

### Valuation algebra

Correct for nonnegative sequences for which the displayed valuations exist
(and hence which are eventually positive unless their valuation is
`-infinity`).  If `u>_lex v` and `j` is the first unequal coordinate, then

```math
\log w_n-\log z_n
=a_{j,n}(u_j-v_j)+o(a_{j,n})\longrightarrow+\infty.
```

The later coordinates are `o(a_{j,n})`, and the original remainder
`o(a_{L,n})` is also `o(a_{j,n})`.  Equal valuations incur only `log 2`, which
is `o(a_{L,n})`.  Products add every coefficient.  The proof also handles an
eventually-zero factor in the evident absorbing way.  No cancellation is
allowed here; positivity is essential.

### Finite composition and queries

Correct.  A fibre of a map between fixed finite alphabets contains only
finitely many positive products, so repeated use of the two-term law gives
the lexicographic maximum.  The same applies to the finite query sum after
adding the valuation of `exp(V_n(q))`.  Calling this “exact composition” is
acceptable only as **exactness of the asymptotic valuation**, not equality of
finite-`n` counts.

### Recovery

Correct as an **abstract count recovery up to a common leading
normalization**.  After choosing `C+sigma_1(q)>0`, dominance of `a_1` makes
every finite exponent tend to `+infinity`; flooring changes its logarithm by
`o(1)`.  This construction does not preserve a prescribed total mass, code,
graph, or other model constraints.  The draft already acknowledges that
limitation, and it should remain explicit wherever the theorem is promoted.

### Growing fibres

The largest-summand bound is correct.  With uniform `o(a_L)` term errors,
`log |I_n|=o(a_L)` is a uniform sufficient condition.  The equal-summand
example proves **worst-case sharpness**: `exp(ca_L)` co-maximizers add `c` to
the last coordinate.  It is not an iff condition for a particular structured
family—large fibres can still be harmless when all but subexponentially many
terms are sufficiently suppressed.  “Sharp branching boundary” should
therefore be read or relabelled as “sharp uniform worst-case boundary.”

## 2. Positive benchmark

The calculation is correct for a fixed finite **acyclic** rooted type
skeleton (fixed depth and a fixed finite set of path types).  Every path
contributes the vector

```math
\left(\sum_{e\in p}s_e,\sum_{e\in p}t_e\right),
```

and the finite path sum selects its lexicographic maximum.  This is a genuine
two-speed distinction invisible to the speed-`n` carrier.  It does not cover
a type graph iterated for a depth growing with `n`, because then the path
fibre is growing and MS.1(2) no longer applies without checking its entropy at
the last speed.  The draft's explicit “deterministic finite-level skeleton”
qualification is therefore necessary.

## 3. Two-grid calculation

The numerical asymptotics are correct.  Writing `q=k/n` on the fine grid,

```math
\sum_{q\in G_n}\lfloor e^{n(1-q^2)}\rfloor
=e^n\sum_{k=-n}^n e^{-k^2/n}+O(n)
=e^n(\sqrt{\pi n}+o(\sqrt n)).
```

On the coarse grid `q=k/sqrt(n)`, the corresponding Gaussian sum converges
to `sum_(k in Z)e^{-k^2}`.  Thus the collapsed valuations at speeds
`(n,log n)` really are `(1,1/2)` and `(1,0)`.

The necessary scope repair concerns what “same coefficient hypograph” means.
The graphs

```math
\{(q,1-q^2,0):q\in G_n\},\qquad
\{(q,1-q^2,0):q\in \widetilde G_n\}
```

have the same limit in ordinary product/Hausdorff topology.  But the coarse
mesh cannot recover a general fixed `q` to second-speed precision: its
descriptor error is `O(n^{-1/2})`, hence its leading exponent error is
generically `O(sqrt(n))`, not `o(log n)`.  A topology explicitly requiring
leading-coordinate recovery to `o(log n/n)` would distinguish the grids.
Accordingly the example proves:

> ordinary limiting coefficient graphs (or a bare pointwise limiting roof)
> do not retain tangent counting mass.

It does **not** prove that every possible “multiscale hypograph topology”
forgets that mass.  The draft should name the ordinary topology rather than
leave “multiscale coefficient hypograph convergence” undefined.

## 4. Vandermonde calculation

Correct after adding the standard hypotheses `p in (0,1)` and choosing an
integer subsequence with `pn` integral.  Uniform Stirling asymptotics near the
interior saddle give

```math
\log {n\choose pn}=nh(p)-\tfrac12\log n+O(1).
```

The largest term in Vandermonde's sum occurs at `k=pn` and has vector
`(2h(p),-1)`, whereas

```math
\log {2n\choose 2pn}=2nh(p)-\tfrac12\log n+O(1),
```

with vector `(2h(p),-1/2)`.  The missing half-coordinate is exactly the
`Theta(sqrt(n))` saddle window.  This example avoids the topology ambiguity:
even exact pointwise coefficient maximization fails when the decomposition
fibre has entropy at the last retained speed.

It falsifies the **bare pointwise coefficient carrier**, not every decorated
finite-speed carrier.  A tangent density/amplitude plus nondegenerate Hessian
is a plausible Morse-class repair, but the draft has not proved its
sufficiency, minimality, or closure.  The sentence about such a repair should
be labelled as a proposed next theorem, not as an established conclusion.

## 5. Required edits before canonical promotion

1. Specify that the two-grid graphs converge in ordinary
   descriptor/coefficient Hausdorff topology; avoid an undefined claim of
   convergence in a scale-aware multiscale topology.
2. State worst-case sharpness, rather than necessity, of the branching
   condition.
3. Add `p in (0,1)` and the integrality subsequence to Vandermonde.
4. Keep the positive benchmark restricted to a fixed finite acyclic
   skeleton.
5. Mark the tangent-mass/Hessian repair as conjectural.

With those repairs, the director judgment is justified: MS.1 is a rigorous
subexponential refinement of the existing response algebra, while the
growing-fibre examples explain why it is not yet an autonomous compactness
theory.
