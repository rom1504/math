# Independent full audit: geodesic synchronization beyond scalar BLR

**Status.**  Fresh adversarial reconstruction.  Theorems GS.1 and
GS.5--GS.11 in `phase3_geodesic_synchronization.md` survive with their
stated constants and scopes.  This audit deliberately treats the scalar
BLR lemma as already checked and concentrates on the shortest-word,
all-future, extension, counting, and bent-family arguments.

The new verifier
[`verify_phase3_geodesic_sync_full_audit.py`](../experiments/verify_phase3_geodesic_sync_full_audit.py)
independently checks the complete selector formula through `h=4` and the
coupled-bent construction at both `k=4` and `k=8`.  The latter is important:
it checks the first nontrivial block extension rather than only the displayed
four-dimensional seed.

Two expository repairs are recommended before surface promotion.  Neither
changes a theorem.

1. In the binary specialization following GS.7, the general theorem requires
   `A subseteq G\K`, whereas the text takes `A=S\B`.  This is valid because
   GS.2 applied to a singleton first proves `S cap W=B`; that one-line
   implication should be stated.
2. In GS.11, either choose the canonical section `L_0(q)=(0,q)` or interpret
   `(z,q)` relative to the chosen `L_0`.  With an arbitrary section
   `L_0(q)=(a(q),q)`, the displayed residual in the diameter proof is
   `z+a(q)+f(q)`, not literally `z+f(q)`.  Linearity makes the same proof go
   through.

## 1. GS.1: exact shortest-word criterion

Let `A subseteq S` represent `t=sum B`, put `I=A cap B` and `R=A\B`.
Projection forces `sum pi(R)=0`, and in `B` coordinates

```math
\sum R=t+\sum I=\sum_{b_i\notin I}b_i.
```

Thus cycle contraction gives `D-|I|<=|R|`, hence `|A|>=D`.  Conversely, a
projected-zero `R` of kernel weight `r_B` gives the representation

```math
t=\sum R+\sum_{i\notin\operatorname{supp}_B(\sum R)}b_i
```

of length at most `|R|+D-r_B`; shortestness forces `r_B<=|R|`.  Repeated
letters cause no gap in the binary group because they cancel in pairs.
This proves the advertised **if and only if**, and it uses shortestness of
the chosen word rather than diametrality.

A small consequence used implicitly later is

```math
S\cap W=B.                                      \tag{GA.1}
```

Indeed, if `s in (S\B) cap W`, the one-element cycle gives `|s|_B<=1`.
Since `s` is nonzero, it must equal a member of `B`, a contradiction.

## 2. GS.5: complete fibres and arbitrary raw futures

Fix any word over `S union T` and let `R` be its letters from `S\B` after
assigning overlaps with `T` to either side.  Put `q=sum pi(R)`.

* If `q=0`, GS.2 replaces `R` by at most `|R|` basis letters.
* If `q ne 0`, the cycle `R triangle {s_q}` has at most `|R|+1` letters.
  Therefore `sum R+s_q` has `B`-weight at most `|R|+1`.  The synchronized
  error `s_q+L(q)` has weight at most nine.  One `L(q)` and at most
  `|R|+10` basis letters replace `R`, an overhead of at most eleven.

Nothing here assumes that the future is disjoint from, spans with, or has
any structural relation to the source.  This proves the first inequality
for every raw `T` and every root.

Conversely, if a word contains `r` generators `L(q_i)`, linearity replaces
their sum by zero when `sum q_i=0`, and otherwise by one `L(q)`.  In the
second case `r>=1`, and replacing that generator by `s_q` plus at most nine
basis letters increases length by at most `10-r<=9`.  Hence the two
directional constants are respectively `11` and `9` exactly as printed.

The fibre estimate used upstream is also sound: two elements in one
nonzero fibre form a two-letter quotient cycle, so every full-fibre member
is within Hamming distance two of the chosen representative and within
eleven of `L(q)`.  The response theorem itself is stronger than a covering-
radius comparison: it controls the complete rooted word profile after every
raw append.

## 3. GS.6: cover size and telescoping scope

For fixed `B`, the set of linear sections is a torsor for
`Hom(Q,W)`, and therefore has exactly `2^(Dk)` elements.  The projection
label identifies each generator `L(q)`, so no extra truth table is hidden in
the centre `S_L`.  If the chart varies, an ordered basis and the values of a
section on a quotient basis use at most

```math
wD+wk=w^2
```

binary coordinates.  These facts verify the cover count.

For `m` sources, replace them one at a time.  At the `i`th replacement, the
union of the external future with all other raw or already replaced sources
is a legitimate arbitrary raw future in GS.5.  Pointwise triangle
inequality yields `11m`.  This is a bounded-composition **telescoping**
statement.  It does not define a product on summaries and does not claim an
error controlled only by the final union, so it does not overstate exact
closure.

## 4. GS.7: finite-abelian conventions and the factor `2h`

The correct general-group convention is an undirected Cayley word metric:
`A=-A`, `B=-B`, and source words may repeat letters.  The sign-compatible
choice

```math
s_{-u}=-s_u
```

is essential.  In particular, when `u=-u`, it requires an order-two lift;
this need not exist in an arbitrary extension and is properly an explicit
hypothesis.  It is automatic over an elementary abelian binary group.

If the `r` original source letters sum to `a` and a shortest quotient word
of length `m<=h` has selected lift sum `c`, then the word formed from the
original letters and the inverses of the selected lifts is a projected-zero
source word of length `r+m`.  Cycle contraction represents `a-c` by at most
`r+m` kernel letters.  Thus `a` is represented by

```math
m+(r+m)=r+2m<=r+2h
```

letters in the transversal-plus-kernel support.  Common future and kernel
letters are untouched.  This proves the full all-context upper bound.  The
reverse inequality is monotonicity because the transversal is a subset of
the full source.

In the binary specialization, (GA.1) justifies taking `A=S\B`: it really is
disjoint from the kernel.  Repeated binary words reduce to subsets, and
their reduced cardinality can only decrease, so GS.1 supplies the word-form
cycle hypothesis required by GS.7.

## 5. GS.8: exact selector cube

The construction is a direct product of `h` independent three-coordinate
blocks.  In block `i`, a queried `c_i'` costs one when `i in J`; otherwise it
is uniquely obtained as `c_i+e_{2i-1}+e_{2i}` at cost three.  A block not in
the query costs zero.  Consequently

```math
\ell_{S_J}(x_P)=|P|+2|P\setminus J|.
```

Subtracting two profiles gives

```math
2\bigl(|P\cap(K\setminus J)|-|P\cap(J\setminus K)|\bigr).
```

Selecting one of the two directed differences proves the exact supremum

```math
2\max\{|J\setminus K|,|K\setminus J|\}.
```

The audit verifier exhausts all `J,K,P` through `h=4` (340 rooted-query and
340 profile-pair checks).  Thus `2h` is both the theorem's upper constant
and an attained response gap, not merely a loose proof charge.

## 6. GS.9: partial projection and both constants

The revised theorem correctly assumes `Q` nontrivial, hence `h>=1`; when
`Q=0`, (GA.1) gives `S=B` and exact recovery.  This separation is needed to
avoid the meaningless displayed reverse bound `10h-1=-1` at `h=0`.

Choose shortest selected-representative words `C_q`, of length at most
`h`, and let `g(q)=sum C_q`.  The symmetric difference of
`C_x,C_y,C_{x+y}` is a quotient cycle of cardinality at most `3h`.
The vector BLR theorem therefore gives

```math
\max_q|g(q)+L(q)|_B\le9h.
```

For a raw source subword `R` with nonzero quotient total `q`, comparison
with `C_q` costs at most `|R|+h` kernel letters.  Adding the synchronization
error and one graph generator gives total length

```math
1+|R|+h+9h=|R|+10h+1.
```

In the reverse direction, a nonempty collection of graph generators first
collapses to one `L(q)`, which is replaced by at most `h` selected source
letters and `9h` kernel letters.  Relative to at least one original graph
letter, the overhead is at most

```math
h+9h-1=10h-1.
```

Zero quotient totals delete exactly.  Hence both constants in GS.25 are
correct, and the proof remains uniform over raw futures.  When `h=1`, they
reduce to the `11` and `9` constants of GS.5.

## 7. GS.11: coupled bent pairs and diametrality

The three displayed alternating forms are nonsingular and obey
`beta_3=beta_1+beta_2`; block-diagonal repetition preserves both facts.  For
a projected-zero set `R`:

* odd cardinality gives exactly one one in each paired output, hence weight
  three, and a nonempty odd zero-sum set has at least three elements;
* even cardinality at least six has output weight at most six;
* a four-element zero-sum set is an affine two-flat.  Its three pair bits
  are `beta_i(u,v)`, and the relation among the forms prevents all three
  from being one, so its total output weight is at most four.

This exhausts all possible cycles and proves GS.2.

For diametrality, take the canonical splitting for clarity.  A target over
nonzero `q` uses its single quotient generator and at most `D-1` basis
letters unless the residual is the all-one vector `t`.  In the exceptional
case, split `q=q_1+q_2` with both summands nonzero.  The three paired defects
have total weight three, so the two quotient generators plus `D-3` basis
letters use `D-1` steps.  Every point is therefore within `D`, while GS.1
puts `t` at distance exactly `D`.

Finally, each quadratic coordinate is bent.  On all of `Q`, its nearest
linear distance is `2^(k-1)-2^(k/2-1)`.  On `Q\{0}`, the `c+phi` coordinate
is the complement and has nearest distance one smaller.  The six coordinate
choices of a vector-valued linear map are independent, so the exact total is

```math
3\left(2^k-2^{k/2}-1\right).
```

The independent verifier obtains `33` at `k=4` and `717` at `k=8`, exactly
the formula; it checks every nondegenerate additive triangle at both sizes
and computes full Cayley diameters six.  This supports the block-family
argument, while the proof above supplies all larger multiples of four.

## 8. Numbering, information scope, and final judgment

The theorem-like headings occur exactly once and consecutively as GS.1
through GS.11.  Equation tags occur exactly once and consecutively as GS.1
through GS.31.  Sharing the `GS` prefix between theorem and equation
numbering is typographically conventional here and creates no ambiguous
cross-reference in the current prose.

The information scope is also consistent:

* GS.5 and GS.9 are all-raw-context response comparisons;
* GS.6 is a cover plus a finite telescoping law;
* GS.7 strips full fibres against arbitrary raw contexts;
* none of them claims an exact summary homomorphism or an error independent
  of the number of already summarized fragments.

No mathematical flaw was found in the audited claims.  The substantive law
that survives is sharper than scalar stability alone: global shortest-word
compatibility synchronizes an entire quotient-indexed family, while the
quotient Cayley diameter exactly governs how microscopic fibre information
can amplify into future rooted responses.  The selector construction proves
that the `h` dependence cannot be removed under these hypotheses.
