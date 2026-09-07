# Wave 37: internal replacement controls one endpoint, not completion width

## Status

The interval identities, deletion/replacement inequality, and two-sign
self-exposure dichotomy below are **proved**.  The `A_6,A_8,A_9` tables are
**exact finite evidence**, checked by
`tmp/width_replacement_r37_check.py`.  No project-scale upper bound on the
completion-width deficit is obtained.

The main negative conclusion is quantitative: replacing one internal block
by either sign/switching of an exact order-`m` minimizer can force only one
outside-completion endpoint to height `q_n-q_m`.  Even if both replacements
are exposed by the prescribed child label, the resulting bound on `b_S` has
bulk size `(q_n-q_m)/2=Theta(n^(3/2))` at fixed density.  This is much larger
than the required `O(n(d+1))`.  Witness migration makes the statement still
weaker, and cannot be silently discarded.

## 1. Exact fiber interval bookkeeping

Fix `S`, write `q=q_n`, and canonically represent a projective word `y` on
`S`.  With the notation of (10.1034), put

```math
Z_+(y)=\max f_{S,y},\qquad Z_-(y)=-\min f_{S,y},
```

and define its half-width and interval center by

```math
w(y)=\frac{Z_+(y)+Z_-(y)}2=q-b(y),
\qquad
a(y)=\frac{Z_+(y)-Z_-(y)}2.
\tag{R37.W1}
```

Uniform outside spins have mean zero in `f_(S,y)`, so `Z_+,Z_-` are both
nonnegative.  Hence

```math
\boxed{|a(y)|\le q-b(y).}
\tag{R37.W2}
```

If `e_A(y)=y^TA[S]y`, the two parent response orientations are bounded by
`q`, and adding/subtracting those bounds gives the exact fiber constraint

```math
\boxed{|a(y)+e_A(y)|\le b(y).}
\tag{R37.W3}
```

Thus small completion-width deficit already requires the center of the
outside interval to cancel the signed child energy very accurately.

## 2. What exact block replacement really proves

Let `A^(0,S)` be the parent with the internal `S`-block set to zero, and let
`C` be any exact order-`m` minimizer, `Q(C)=q_m`.  Replacing the zero block by
any switching or either global sign of `C` gives a complete signing.  Exact
order-`n` minimality and the triangle inequality therefore give

```math
q\le Q(A^{(0,S)}+C)\le Q(A^{(0,S)})+q_m.
```

On the other hand, the zero-block norm is exactly

```math
Q(A^{(0,S)})
=\max_{[y]}\max\{Z_+(y),Z_-(y)\}
=\max_{[y]}\{q-b(y)+|a(y)|\}.
```

Consequently

```math
\boxed{
\min_{[y]}\{b(y)-|a(y)|\}\le q_m.
}
\tag{R37.W4}
```

This is the direct scalar conclusion of the one-block deletion argument: it
concerns an unrestricted projective word, not a favorable child word.

There is a useful center bound.  Since

```math
f_{S,y}(x)+f_{S,y}(-x)=2x^TA[S^c]x,
```

evaluating this identity at a maximizer and at a minimizer of `f` gives

```math
|a(y)|\le Q(A[S^c]).
```

Combining with (R37.W4),

```math
\boxed{
\min_{[y]}b_S(y)\le q_m+Q(A[S^c]).
}
\tag{R37.W5}
```

Again, the minimizing label is unrestricted.  The bound is of bulk order
and supplies no favorable support.

## 3. Both signs give an exact self-exposure-or-migration dichotomy

Now let `[y]` be a chosen child-favorable projective word and orient its
signed child energy positively:

```math
E_y=\sigma_y y^TA[S]y>0.
```

Relabel the two response orientations by `sigma_y`; equivalently replace
`a(y)` below by `sigma_y a(y)`.  Then the parent constraint is
`|a(y)+E_y|<=b(y)`, while every expression symmetric in `|a+-q_m|` is
unchanged.

Switch an exact order-`m` minimizer so that this same `y` is a positive
ground.  In the two full signings obtained from its internal blocks `+C_y`
and `-C_y`, the maximum score **within the same projective fiber** is

```math
P_\pm(y)=q-b(y)+|a(y)\pm q_m|.
\tag{R37.W6}
```

Define the two-sign self-exposure shortfall

```math
\varepsilon_y
=\max\{[q-P_+(y)]_+,[q-P_-(y)]_+\}.
\tag{R37.W7}
```

If either term is positive, exact minimality says only that the corresponding
replacement norm reaches at least `q` on **another** projective fiber.  This
is witness migration, not favorable support.

From (R37.W6)--(R37.W7),

```math
b(y)\le\min\{|a(y)+q_m|,|a(y)-q_m|\}+\varepsilon_y
=\big||a(y)|-q_m\big|+\varepsilon_y.
\tag{R37.W8}
```

Together with the parent constraint `|a+E_y|<=b` and (R37.W2), this gives
the exact dichotomy

```math
\boxed{
\varepsilon_y\ge E_y-q_m
\quad\text{or}\quad
2b_S(y)\le q_n-q_m+\varepsilon_y.
}
\tag{R37.W9}
```

Indeed, if `epsilon_y<E_y-q_m`, then (R37.W3) and (R37.W8) rule out
`|a(y)|<=q_m`.  Thus `b<=|a|-q_m+epsilon`, and (R37.W2) yields (R37.W9).

This explicitly survives the `A_6/A_9` wall: no migrated witness is called
favorable.  It also shows why defeating migration would still be
insufficient.  Even `epsilon_y=0` gives only

```math
b_S(y)\le\frac{q_n-q_m}{2}=\Theta(n^{3/2})
```

at fixed density, whereas (10.1038) requires
`b_S(y)=O(n(d+1))=o(n^(3/2))` in the project regime.  Moreover, at the full
favorable allowance one can have `E_y<=q_m`, in which case the first branch
of (R37.W9) is vacuous from the start.

## 4. Exact finite falsification of natural favorable transfers

For every selector, the checker minimizes `b_S(y)` over exact child-ground
projective labels.  The histograms are

| parent, `m` | `q_n` | `q_m` | histogram of best child-ground `b_S` |
|:--|--:|--:|:--|
| `A_6,3` | 10 | 6 | `2^20` |
| `A_8,5` | 20 | 8 | `4^16,8^16,10^8,12^16` |
| `A_9,6` | 24 | 10 | `8^14,12^16,16^50,20^4` |

The attractive pointwise proposal

```math
\min_{y\in F_S^{\rm gr}}b_S(y)\le q_n-q_m
\tag{R37.W10}
```

holds on `A_6,A_8` but fails on `54` of the `84` order-six selectors of
`A_9`; its maximum excess is `6`.  Even its uniform-selector average fails
on `A_9`, by the exact amount

```math
\mathbb E_S\min_{y\in F_S^{\rm gr}}b_S(y)-(q_9-q_6)=\frac2{21}.
\tag{R37.W11}
```

More decisively, the favorable version of the proved unrestricted bound
(R37.W5) is false.  For four `A_9` selectors, including

```text
S={0,1,2,3,5,7},
```

every exact child ground has `b_S>=20`, while

```math
q_6+Q(A_9[S^c])=10+6=16.
```

Thus the label supplied by block deletion is necessarily outside the exact
child-ground fiber in these cases.  This is an exact width-level instance of
witness migration, not merely a failure of one minimax proof.

The checker also verifies (R37.W9) for every exact child ground.  The simpler
finite inequality

```math
\min_{y\in F_S^{\rm gr}}b_S(y)
\le q_n-q_m+Q(A[S])-q_m
\tag{R37.W12}
```

happens to hold on all three matrices, but no proof is known and its right
side remains `Theta(n^(3/2))`; it is recorded only as finite evidence.

## 5. Precise successor

One-block internal replacement should not be pursued as a standalone proof
of the width term in (10.1038).  It is off scale even under perfect
same-fiber exposure, and unrestricted deletion witnesses fail favorable
support exactly on `A_9`.

A viable successor must add a mechanism which cancels the bulk
`q_n-q_m` **before** estimating width.  Concretely, it would have to be a
joint anchored-family theorem or a multi-block replacement law producing
row-qualified favorable labels with

```math
\mathbb E_{S\in\mathcal G}b_S(y^S)=O(n(d+1)),
\qquad
\log\beta^{-1}=O(TL_0),
```

while charging every migrated replacement witness to a cross-selector
collision or to the already-centered recurrence gap.  Merely proving a
positive favorable mass in each separate replacement game is not enough:
the unavoidable one-block bulk term is already too large.
