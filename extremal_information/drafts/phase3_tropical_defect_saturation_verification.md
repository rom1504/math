# Independent verification: tropical defect saturation

**Status.** Independent adversarial audit of
`phase3_tropical_defect_saturation.md`.  The central defect identity, its
arbitrary-tree interpretation, the convex contrast, and the fixed-chart
syndrome theorem are correct after the scope corrections below.  I found no
counterexample to the substantive claims.

## 1. The abstract theorem behind TDS.8

TDS.8 is a specialization of a metric-monoid identity.  Stating this first
would both shorten the proof and clarify exactly what is new.

### Proposition V-TDS.1 (power-orbit defect law)

Let `(M, product, 1, d)` be a commutative monoid with an extended metric such
that right translation is nonexpansive:

```math
d(xz,yz)\le d(x,y).
```

For fixed `b in M`, put `P_b(x)=xb`.  Then, for every `m>=1`,

```math
\sup_{x_1,\ldots,x_m\in M}
d\!\left(\prod_{i=1}^mP_b(x_i),
P_b\!\left(\prod_{i=1}^m x_i\right)\right)
=d(b^m,b).                                      \tag{V-TDS.1}
```

If every leaf is blurred once and `s` further blurs occur at arbitrary
internal nodes of a composition tree, the final state is

```math
\left(\prod_i x_i\right)b^{m+s}.                \tag{V-TDS.2}
```

Consequently its worst-case defect relative to one final blur is exactly
`d(b^(m+s),b)`.  Uniform stability under arbitrary size, depth, and internal
reblurring is equivalent to

```math
\sup_{k\ge1}d(b^k,b)<\infty.                    \tag{V-TDS.3}
```

#### Proof

Commutativity gives

```math
\prod_i P_b(x_i)=\left(\prod_i x_i\right)b^m.
```

Translation nonexpansiveness bounds the distance from the once-blurred
product by `d(b^m,b)`.  Setting every `x_i=1` gives equality.  Each additional
internal `P_b` contributes one further factor of `b`, proving (V-TDS.2) and
the remaining assertions. `square`

For min-plus profiles, `M` is the convolution monoid and the finite kernel
power orbit decreases to `b_*`; hence (V-TDS.3) equals
`||b-b_*||_infinity`.  For compact convex response bodies under Minkowski
addition, it becomes the orbit `K,2K,3K,...`, whose radius grows linearly.
Thus TDS.2 and TDS.2b are opposite specializations of the same exact law.

This abstraction is elementary rather than a replacement for the
model-specific content.  The informative step in a given model remains to
show that a lossy element has a bounded power orbit and that one blur is
already accurate on the relevant landscape class.

## 2. Audit of the min-plus statements

Let profiles be **proper** when they have at least one finite value.  Proper
extended profiles are closed under convolution, and convolution with a
finite kernel makes them finite-valued.  On this domain the sup-norm proof of
TDS.8 is valid verbatim:

```math
(P_bf_1)\star\cdots\star(P_bf_m)
=(f_1\star\cdots\star f_m)\star b^{\star m}.
```

Translation nonexpansiveness supplies the upper bound, while
`f_i=delta_0` attains it.  This also proves that equality in TDS.8 is not an
artifact of a loose triangle inequality.

The kernel powers decrease because a zero-cost step can be appended.  A
minimum-cost walk on the finite Cayley graph may have all repeated-vertex
cycles deleted, since costs are nonnegative.  It then has at most `|G|-1`
nonzero steps, and zero steps pad it to every larger length.  Thus the powers
stabilize finitely at the shortest-path closure.  The claims that `b_*` is the
greatest subadditive minorant and that `b_* star b_*=b_*` follow.

The internal-reblur claim is exact, not merely an upper bound.  With `m`
blurred leaves and `s` internal reblurs, the total exponent is `m+s`; no tree
shape remains after associativity and commutativity.  Comparing this state to
the one-blur exact product costs at most `Delta(b)`.

There is an important interpretation boundary.  TDS.2 controls the
**algebra defect relative to one blur**.  It does not by itself control error
relative to the unblurred exact product, and it does not by itself prove any
state-space compression.  The former requires the extra one-shot hypothesis
`||P_bf-f||<=alpha` in TDS.2a; the latter requires a separate quotient or
description-size argument.  With this wording, the claim of a genuinely
nonexact yet arbitrary-depth-stable approximate homomorphism is correct.

The `Z/3Z` example is correct: for `b=(0,1,2+eta)`, `0<eta<=1`, the closure is
`(0,1,2)` and the exact orbit defect is `eta`.

## 3. Necessary corrections and scope edits

1. **Exclude the identically-infinite profile.**  The sentence saying that
   convolution with a finite kernel produces finite output is false for
   `f identically +infinity`.  TDS.2 should quantify over proper extended
   profiles, or more simply over finite-valued profiles.  The norm comparison
   is then well-defined.

2. **Separate indicator kernels from the finite-kernel hypotheses.**
   `iota_K` takes `+infinity`, so Corollary TDS.3 is not literally an instance
   of Section 1 as stated.  Its convolution identities remain correct in the
   extended reals, and the complete-profile defect is infinite whenever a
   later sumset strictly exceeds `K`; say this directly.

3. **State the relevant Hamming error explicitly.**  For `f(x)=|x|` and
   `s=min(mr,w)`,

   ```math
   ||f-P_{B_s}f||_infinity=s,
   \qquad
   ||P_{B_r}f-P_{B_s}f||_infinity=s-r.           \tag{V-TDS.4}
   ```

   The draft states the first quantity correctly, but the second is the
   finite Lipschitz-class analogue of the one-blur-versus-many-blur algebra
   defect.  Recording both prevents the raw smoothing error from being read
   as the additional composition error.

4. **Clarify the convex-body domain.**  If all response bodies are compact
   convex, support functions prove the displayed Hausdorff equality directly.
   If the supremum is meant over arbitrary compact bodies, Minkowski
   nonexpansiveness gives the upper bound and choosing all exact bodies
   `{0}` gives equality.  The present proof silently uses convexity of the
   bodies, although the proposition only explicitly assigns it to `K`.

5. **Give the range of the fixed-chart parameter.**  TDS.4 requires
   `0<=r<=w`; the nontrivial asymptotic following TDS.19 should state
   `0<epsilon<1/2`.  The endpoint `epsilon=1/2` may be handled separately by
   taking `H=G`.

6. **Justify the word “strict.”**  The quotient is genuinely noninjective,
   but the draft does not exhibit a collision.  For example, when
   `1<=r<=w-2`, choose nonzero `h in H` and
   `v=e_(r+1)+e_(r+2)`.  The distinct supports

   ```math
   E union {v},\qquad E union {v+h}
   ```

   have the same image support in `G/H`, hence the same quotient word
   profile.  A one-line witness would establish strictness rather than leave
   it implicit in the state count.

7. Remove the duplicated sentence `outside it. Then` at the start of
   Section 3.

None of these edits changes a theorem constant.

## 4. Audit of the fixed-chart syndrome quotient

The identity

```math
\lambda_{S\cup T}=\lambda_S\star_G\lambda_T
```

is correct even when `S` and `T` overlap: a shortest union representation can
be assigned to the two supports, while concatenating two representations can
only become shorter after duplicate atoms cancel.  Subgroup min-filtering is
an exact homomorphism from convolution on `G` to convolution on `G/H`, proving
TDS.16.

Since the common basis lies in every support, `lambda_S` is one-Lipschitz in
the coordinate Hamming metric.  Every `H`-coset has diameter `r`, so its
minimum is within `r` of every value in that coset.  This gives exactly the
interval in TDS.17 and midpoint error `r/2`.  The quotient has `2^(w-r)`
integer entries in `{0,...,w}`, yielding the asserted state-count upper bound.
These constants are correct.

The qualifier **fixed chart** is essential.  This is an exact algebra only
for the union-closed class of supports containing the same basis `E`; it does
not supply compatible quotients for arbitrary support-selected bases.  Within
that stated class it is both genuinely lossy and stable under compositions of
arbitrary size and depth.

## 5. Reproducible finite checks

`experiments/verify_phase3_tropical_defect_saturation.py` checks:

- every kernel with costs in `{0,1,2,3}` on cyclic groups of orders two
  through four;
- 171,624 proper-profile nonexpansiveness comparisons, including exact
  attainment by `delta_0` and finite power stabilization;
- the two distinct Hamming errors in (V-TDS.4) through width eight;
- TDS.16 and TDS.17 for every fixed-chart support pair through width three,
  and 512 deterministic support pairs for every chart at width four; and
- explicit noninjectivity of a nontrivial fixed-chart quotient.

The saved result is
`experiments/phase3_tropical_defect_saturation_results.json`; all checks pass.

## Verdict

Promote the power-orbit defect law as the abstract theorem, then present
tropical shortest-path saturation and convex linear growth as its two sharp
specializations.  After the seven scope edits above, TDS.1--TDS.4 are
rigorous.  The result is a real general law about arbitrary-depth response
composition, while its limitations are equally clear: bounded algebra defect
does not alone imply bounded distortion from the true response or a smaller
feature state.
