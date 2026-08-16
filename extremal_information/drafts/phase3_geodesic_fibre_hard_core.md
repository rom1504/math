# Geodesic fibres and an exponentially thinner syndrome hard core

**Status.** Independent theorem audit.  All statements below are proved.  The
structural part was checked for every spanning support through `w=4` and on
deterministic random samples at `w=5,6` by
[`verify_phase3_geodesic_fibre_bound.py`](../experiments/verify_phase3_geodesic_fibre_bound.py).
The computation also records the necessary `D=2` endpoint exception.

This strengthens the coarse additive-combinatorial count in
`phase3_hard_core_rees_compression.md`.  A support of large Cayley diameter
must be sparse not only globally: relative to every diametral geodesic, its
intersection with each affine fibre is a diameter-two Hamming anticode.  The
resulting exact count gives a depth-stable syndrome quotient with
`2^{(1-2 epsilon)w+O(log w)}` **bits** at response error `epsilon w`.

## 1. Setup

Let `G=F_2^w`, let `S subseteq G\{0}` span `G`, and put

```math
\ell_S(x)=\min\{|A|:A\subseteq S,\ \sum_{a\in A}a=x\},
\qquad
D(S)=\max_{x\in G}\ell_S(x).                    \tag{GF.1}
```

Binary repetitions cancel, so subsets suffice.  This is the Cayley diameter
of `S`, and it is the covering radius/response observable of the associated
full-rank syndrome fragment.

Choose `t` with `ell_S(t)=D=D(S)` and a shortest representation

```math
t=b_1+\cdots+b_D,
\qquad B=\{b_1,\ldots,b_D\}\subseteq S.          \tag{GF.2}
```

Write `W=span(B)`.  Once `B` is ordered, identify `W` with `F_2^D` by its
`B`-coordinates.  Every affine `W`-coset is then a torsor for this Hamming
cube; pairwise coordinate distances do not depend on the choice of its
origin.

## 2. The geodesic-fibre theorem

### Theorem GF.1 (diametral geodesic fibres)

For every choice (GF.2):

1. `B` is linearly independent, so `dim W=D`;
2. `S cap W=B`;
3. in every nonzero coset `q+W`, the coordinate image of
   `S cap(q+W)` has pairwise Hamming distance at most two;
4. if `D>=3`, every such fibre has at most `D+1` elements.

The restriction in item 4 is necessary.  At `D=2`, a fibre can be the full
two-cube of size four.

#### Proof

If a nonempty subfamily of `B` summed to zero, deleting it from (GF.2) would
give a representation of `t` with fewer than `D` generators (and if it were
all of `B`, it would say `t=0`).  Thus `B` is independent.

Take `s in S cap W` and write uniquely

```math
s=\sum_{i\in I}b_i.
```

The set `I` is nonempty because `0 notin S`.  If `|I|>=2`, then

```math
t=s+\sum_{i\notin I}b_i
```

uses at most `1+D-|I|<D` generators, a contradiction.  Hence `|I|=1`,
which proves `S cap W=B`.

Now take distinct `s,s'` in the same nonzero `W`-coset and write

```math
s+s'=\sum_{i\in I}b_i.
```

Both are outside `W`, hence are distinct from every member of `B`.  If
`|I|>=3`, then

```math
t=s+s'+\sum_{i\notin I}b_i
```

uses at most `2+D-|I|<D` generators.  Therefore `|I|<=2`, proving item 3.

It remains to use the elementary diameter-two anticode bound.  Translate a
family `F subseteq F_2^D` of pairwise distance at most two so that it contains
zero.  Every member then has Hamming weight at most two.  Regard the
weight-two members as edges of a graph.  They are pairwise intersecting, so
they form either a star or a triangle.  If there are no edges, there are at
most `1+D` sets.  With one edge there are at most that edge, its two singleton
endpoints, and zero.  With at least two star edges, the only compatible
singleton is their common centre, giving at most

```math
1+(D-1)+1=D+1.
```

A triangle allows no singleton and has four sets including zero.  This is at
most `D+1` exactly when `D>=3`. `square`

The endpoint example is explicit.  In `F_2^3` take

```math
B=\{e_1,e_2\},
\qquad
S=B\cup(e_3+\operatorname{span}B).
```

Then `D(S)=2`, while the displayed nonzero fibre has size four.

## 3. Counting the hard core

The proof above also controls the number of possible fibres.  Let `a_D` be
the number of subsets of `F_2^D` with pairwise Hamming distance at most two.
For `D>=3`,

```math
a_D
\le
1+2^D\left(2^D+D2^{D+1}+{D\choose3}\right)
\le 4D\,2^{2D}.                                  \tag{GF.3}
```

The initial one accounts for the empty fibre.  For a nonempty fibre, choose
a translating member in at most `2^D` ways.  After translation
the weight-two graph is empty, a star, or a triangle.  The empty case has at
most `2^D` singleton subfamilies.  For the star case, choose its centre, an
arbitrary **nonempty** incident-edge subfamily, and at most four choices for
the compatible singleton subfamily, giving the deliberately loose term
`D2^{D+1}`.  A nonstar
pairwise-intersecting edge family is a full triangle.  This proves (GF.3);
overcounting translated families is harmless.

### Theorem GF.2 (enumeration of large-diameter supports)

Let `N_(w,R)` be the number of spanning supports
`S subseteq F_2^w\{0}` with `D(S)>=R`, where `4<=R<=w`.  Then

```math
\log_2 N_{w,R}
\le
w^2
+2^{w-R}(2w+2+\log_2 w)
+\log_2 w.                                       \tag{GF.4}
```

In particular,

```math
\log_2 N_{w,R}=O\!\left(w^2+w2^{w-R}\right).     \tag{GF.5}
```

#### Proof

Fix a diameter `D>=R`.  Choose and order a geodesic basis `B`; there are at
most `(2^w)^D=2^{wD}` choices.  The zero `W`-fibre is exactly `B`.  There are
`2^{w-D}-1` nonzero fibres, and Theorem GF.1 plus (GF.3) gives at most `a_D`
possibilities for each.  Thus the number at diameter `D` is at most

```math
2^{wD}a_D^{\,2^{w-D}-1}.                        \tag{GF.6}
```

Sum (GF.6) over the at most `w` values `R<=D<=w`, take logarithms, and use

```math
wD\le w^2,
\qquad
2^{w-D}\le2^{w-R},
\qquad
\log_2 a_D\le2D+2+\log_2D
                 \le2w+2+\log_2w.
```

This proves (GF.4). `square`

No divisibility, parity, or typical-support hypothesis enters this count.
It applies to every full-rank binary syndrome support.

## 4. A depth-stable approximate response algebra

For an integer `4<=R<=w`, define

```math
q_R(S)=
\begin{cases}
\bot,&D(S)<R,\\
S,&D(S)\ge R.
\end{cases}                                      \tag{GF.7}
```

Give these states the product

```math
\bot\odot z=\bot,
\qquad
S\odot T=q_R(S\cup T)                            \tag{GF.8}
```

on retained spanning supports.  This is well defined because Cayley diameter is
antitone under support union: if either hidden support has diameter below
`R`, every later union does too.

### Corollary GF.3 (geodesic hard-core quotient)

Equations (GF.7)--(GF.8) define a commutative associative idempotent response
algebra with at most `1+N_(w,R)` states.  Decode a retained state by its exact
diameter and decode `bot` by `R/2`.  After any number of support unions, the
decoded covering radius differs from the truth by at most

```math
R/2-1.                                           \tag{GF.9}
```

Consequently, for every fixed `0<epsilon<1/2` and every integer

```math
w\ge \max\{1/\epsilon,\,2/(1-2\epsilon)\},       \tag{GF.10}
```

taking

```math
R=\lfloor2\epsilon w+2\rfloor                   \tag{GF.11}
```

gives uniform all-future error at most `epsilon w` with worst-case message
length

```math
2^{(1-2\epsilon)w+O(\log w)}+O(w^2)
```

bits.  In particular this is `o(2^w)` bits.  The more conservative choice
`R asymp epsilon w` already gives the advertised
`2^{(1-epsilon)w+O(log w)}` bound with smaller-than-requested error.

#### Proof

The displayed lower bound on `w` guarantees `4<=R<=w`; it is only a
convenient sufficient threshold, not an optimized endpoint.  The set
`{S:D(S)<R}` is an absorbing union ideal, which proves closure and
associativity of (GF.8).  Every collapsed, still-spanning support has integer
diameter in `[1,R-1]`; the midpoint `R/2` has maximum error `R/2-1` on this
interval.  Retained states are exact.  Theorem GF.2 bounds the number of
retained states.  Under (GF.11), (GF.9) is at most `epsilon w` and

```math
w-R\le(1-2\epsilon)w-1.
```

Substitution in (GF.4) proves the bit bound. `square`

The same decoder answers an arbitrary appended support `U`, even when `U`
does not itself span: apply (GF.7) to the spanning union `S union U`.
If the stored state of `S` is `bot`, antitonicity says that union is still
collapsed; if it is retained, `S` is known exactly.  The binary closed-algebra
statement above concerns the natural source class of spanning fragments.

This quotient is fundamentally different from a one-shot landmark sketch.
Its state count is larger, but its error never accumulates: the detailed
support is retained until it enters the low-response ideal and is then
forgotten irreversibly.  It is also noncircular.  State membership is decided
by the present Cayley diameter, while the counting theorem—not knowledge of
future queries—shows that the uncollapsed hard core is succinct.

## 5. What this proves and what it does not

The result establishes the positive branch of the unrestricted syndrome
dichotomy in a composition-stable form.  It gives a genuine general-theory
mechanism:

> a monotone extremal deficit can admit strong future-response compression
> when its high-deficit level sets have small structural entropy, even if its
> exact response quotient has exponential ambient description size.

The new ingredient is the geodesic-fibre theorem: high response forces a
low-entropy affine-fibre geometry.  The Rees collapse itself is classical.

The exponent is only an upper bound.  The current Grassmannian construction
gives a quadratic-bit lower bound at fixed macroscopic distortion, leaving a
large gap between `Omega(w^2)` and
`2^{(1-2epsilon)w+O(log w)}` bits.  Closing that gap is a separate metric-
entropy problem; it is no longer necessary to establish subexponential
compression.
