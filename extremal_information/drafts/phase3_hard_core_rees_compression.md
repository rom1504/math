# Hard-core/Rees compression for future extremal responses

**Status.** Theorem draft by the response-law adversary.  The abstract Rees
construction and matroid closure lattice are classical.  The project-level
consequence proved here is quantitative: a response ball about an absorbing
terminal state is exactly a forward-stable region that can be collapsed with
uniform error for **every future context**.  For nonnegative antitone
deficits, midpoint decoding identifies the largest possible terminal cell
exactly.  The matroid application gives a third strict
composable quotient outside the code and fixed-rank mean-field examples.

## 1. Response ideals

Let `(Q,star,1,0)` be a commutative monoid with an absorbing element
`0 star q=0`.  Let `F:Q->R` be an observable and define its complete future
response metric

```math
d_F(p,q)=\sup_{c\in Q}|F(p\star c)-F(q\star c)|.   \tag{HRC.1}
```

As in Proposition CRL.1 of `phase3_contextual_response_law.md`, this is a
translation-contractive pseudometric.  For `eta>=0`, put

```math
I_\eta=\{q:d_F(q,0)\le\eta\},
\qquad
H_\eta=Q\setminus I_\eta.                         \tag{HRC.2}
```

The first set is the response-easy ideal; the second is the response hard
core.

### Theorem HRC.1 (quantitative Rees response quotient)

For every `eta>=0`:

1. `I_eta` is a semigroup ideal: `q in I_eta` implies
   `q star u in I_eta` for every `u`.
2. Collapse all of `I_eta` to one absorbing symbol and retain the exact
   response-congruence classes in `H_eta`.  The resulting Rees quotient is a
   closed composition state.
3. Decode a collapsed state using the terminal response
   `c -> F(0 star c)=F(0)`.  Every future answer has uniform error at most
   `eta`.
4. Conversely, if any uniformly `eta`-accurate summary gives `q` and `0`
   the same message, then

   ```math
   d_F(q,0)\le2\eta.                               \tag{HRC.3}
   ```

Thus `I_eta` is always safe to terminal-collapse, while no point outside
`I_(2eta)` can be terminal-collapsed at error `eta`.

#### Proof

Translation contraction and absorption give

```math
d_F(q\star u,0)=d_F(q\star u,0\star u)\le d_F(q,0),
```

proving the ideal property.  The ideal is saturated under exact response
equivalence, since `d_F(q,q')=0` implies
`d_F(q,0)=d_F(q',0)`.  Thus it descends to the exact response quotient.
A Rees quotient by a semigroup ideal is a monoid: any product entering the
ideal remains there.  For `q in I_eta`,
the definition of `d_F` says directly that

```math
|F(q\star c)-F(0\star c)|\le\eta
```

for every future `c`.  This proves the upper guarantee.  If `q,0` share a
summary message, their two true responses are each within `eta` of the same
decoded answer at every context.  The triangle inequality and the supremum
over contexts prove (HRC.3). `square`

This is an exact composition theorem for the summary symbols, not merely a
one-shot approximation of `F(q)`.  Once a partial product reaches the easy
ideal, all its descendants remain collapsed, so repeated composition need
retain detailed states only on `H_eta`.

The factor two is unavoidable at this level.  Two scalar response vectors at
distance `2eta` can both be decoded by their midpoint with error `eta`, while
decoding the collapsed class specifically by the terminal vector protects
only the radius-`eta` ball.

### Corollary HRC.2 (optimal terminal cell for an antitone deficit)

Suppose `F:Q->[0,infinity)` satisfies

```math
F(q\star u)\le F(q),
\qquad
F(0)=0.                                           \tag{HRC.4}
```

Then

```math
d_F(q,0)=F(q),                                    \tag{HRC.5}
```

and therefore

```math
I_\eta=\{q:F(q)\le\eta\}.                        \tag{HRC.6}
```

More sharply, for a target response error `eta`, the larger set

```math
J_\eta=\{q:F(q)\le2\eta\}                        \tag{HRC.7}
```

is an ideal.  Collapse it to one absorbing state and decode that state by the
constant `eta`.  This gives uniform error at most `eta` after every future
composition.  Conversely, `J_eta` is the largest possible summary cell
containing the terminal state under uniform error `eta`.

#### Proof

Every future value `F(q star c)` lies between zero and `F(q)`, while the
identity context attains `F(q)`, proving (HRC.5)--(HRC.6).  If `q in J_eta`,
all its future values lie in `[0,2eta]`, so their distance from the midpoint
`eta` is at most `eta`; antitonicity makes `J_eta` an ideal.  If an
`eta`-accurate summary identifies `q` with the terminal state, the decoder
triangle inequality gives `d_F(q,0)=F(q)<=2eta`, hence `q in J_eta`.
`square`

This is the useful noncircular case: membership in the easy ideal is checked
from the present deficit, rather than by enumerating the full future-response
map.

## 2. Complexity consequence and its limit

Let `Q_F` denote the exact syntactic response quotient.  If it is finite,
Theorem HRC.1 gives an `eta`-accurate closed algebra with at most

```math
1+|H_\eta/\!\equiv_F|
```

states.  More generally, the same statement holds with cardinality replaced
by a description scheme for the hard-core classes.  On the other hand, every
`2eta`-packing in `(Q,d_F)` requires distinct messages under uniform error
`eta`, independently of the Rees construction.

For a nonnegative antitone deficit, midpoint decoding improves this count by
replacing `H_eta` with `{q:F(q)>2eta}`.  The displayed counts are criteria,
not automatic upper bounds.  The
universality theorem CRL.2 can put an arbitrary response table entirely in
the hard core.  Rees compression is substantial only when model structure
forces most composed states toward a terminal response or makes the hard
core itself succinct.

This separates two mechanisms which were previously mixed together:

1. **algebraic closure:** response equivalence supplies a monoid state;
2. **terminal absorption:** a response ball about zero may be forgotten
   forever once reached.

The first is a syntactic-semigroup fact.  The second is the extremal-response
consequence: the same error guarantee survives arbitrary adversarial future
composition.

## 3. Non-code application: matroid residual rank

Let `mathcal M=(E,r)` be a finite matroid of rank `R=r(E)`.  Compose subsets
by union and query the residual-rank landscape

```math
F(S)=R-r(S),
\qquad
R_S(T)=R-r(S\cup T).                              \tag{HRC.8}
```

### Proposition HRC.3 (flat quotient and exact response metric)

Two subsets have identical responses to every future union if and only if
they have the same matroid closure.  Hence the exact response quotient is the
join-semilattice of flats, with

```math
X\vee Y=\operatorname{cl}(X\cup Y).               \tag{HRC.9}
```

On this quotient the complete future-response metric is

```math
d_F(X,Y)=\max\{r(X\vee Y)-r(X),\ r(X\vee Y)-r(Y)\}.                \tag{HRC.10}
```

#### Proof

Equal closures have equal rank after every common union.  For two flats
`X,Y`, the contexts `X` and `Y` attain respectively the two directed join
increments in (HRC.10).  Conversely, if for some context `T` one has
`r(X union T)>=r(Y union T)`, submodularity gives

```math
r(X\cup T)-r(Y\cup T)
\le r((X\vee Y)\cup T)-r(Y\cup T)
\le r(X\vee Y)-r(Y).
```

The opposite ordering is bounded by the other directed increment.  This
proves (HRC.10).  Distinct flats have a positive directed join increment, so
they cannot be response-equivalent. `square`

This quotient is often strict: parallel elements, redundant subsets, and all
spanning subsets are identified without retaining their membership lists.

### Theorem HRC.4 (optimal truncated flat algebra)

Fix an integer threshold `0<=k<=R`.  Define

```math
Q_k=
\{X:X\text{ is a flat and }R-r(X)>k\}\cup\{\star_k\}.             \tag{HRC.11}
```

Compose two retained flats by taking their join (HRC.9), unless the resulting
codimension is at most `k`, in which case output `star_k`; make `star_k`
absorbing.  Then `Q_k` is a closed associative composition algebra.  It
answers

```math
R-r(S\cup T)
```

for every future `T` with uniform additive error at most `k/2`: retained
flats give the exact answer until their join enters the collapsed region, and
`star_k` is decoded as `k/2`.

This is the largest possible terminal cell at that error: no uniformly
`k/2`-accurate summary can identify a flat of codimension greater than `k`
with the spanning terminal state.  For an arbitrary target error `eta`, take
`k=floor(2eta)` and decode the collapsed state by `eta` (or by `k/2`).

#### Proof

Residual rank is nonnegative, antitone under union, and zero on the absorbing
flat `E`.  Corollary HRC.2 identifies the optimal terminal cell at error
`k/2` with the flats of codimension at most `k`, and gives both the algebra
and the converse. `square`

For a binary vector matroid whose ground set contains every nonzero vector of
`F_2^w`, flats are precisely linear subspaces.  A subspace has an `O(w^2)`-
bit basis description, whereas an arbitrary ground-set subset has `2^w-1`
membership bits.  Thus even the **exact** all-future residual-rank state is a
strict sub-landscape quotient; (HRC.11) further collapses every subspace of
codimension at most `k` into one absorbing approximate state.  Repeated union
can only move upward in the subspace lattice, so detailed information is
permanently discarded upon reaching that region.

### Corollary HRC.5 (macroscopic response rate for projective matroids)

For the preceding binary vector matroid and every fixed `epsilon<1/4`, the
deterministic message complexity of answering every future residual-rank
query within uniform error `epsilon*w` is `Theta(w^2)` bits.

#### Proof

The exact flat state gives the `O(w^2)` upper bound.  For even `w=2d`, use
the `d`-dimensional subspaces.  Equation (HRC.10) becomes the injection
metric

```math
d_F(X,Y)=d-\dim(X\cap Y).                          \tag{HRC.12}
```

There are at least `2^(d^2)` such subspaces.  More generally, for
`1<=s<=d`, a greedy constant-dimension packing with pairwise injection
distance at least `s` has size at least

```math
{2^{d^2-s(2d-s)}\over16d}.                       \tag{HRC.13}
```

Indeed, the exact number at injection distance `j` from one subspace is
`2^(j^2){d bracket j}_2^2`; using
`{d bracket j}_2<=4*2^(j(d-j))` bounds the ball of radius below `s` by
`16d*2^(s(2d-s))`.  Divide the total population by this ball bound.

Now fix `epsilon<1/4`, choose `2epsilon<gamma<1/2`, and take
`s=ceil(gamma*w)`.  Equation (HRC.12) gives response separation at least
`s>2epsilon*w`, so distinct packing points cannot share one uniformly
`epsilon*w`-accurate message.  Taking logarithms in (HRC.13), the necessary
number of bits is at least

```math
\left((1/2-\gamma)^2-o(1)\right)w^2              \tag{HRC.14}
```

along even widths.  At odd width, take a fixed line `U`, a complementary
even-dimensional hyperplane `H`, and flats `U direct-sum X` with `X` ranging
over the preceding middle-dimensional packing inside `H`.  Their injection
distances are unchanged, so the same asymptotic order follows. `square`

Thus this model has both a strict exact composable quotient and a sharp
quadratic-order macroscopic response complexity.  The result is a response-
rate theorem, not merely the observation that matroid closure exists.

The same construction applies to residual coverage, rank functions, and
other nonnegative monotone deficits, but the matroid example is stronger
than a renaming: its exact quotient and product are the familiar closure and
join operations, and the response theorem explains precisely why truncating
the flat lattice remains valid against every future extension.

## 4. Code coordinate

For a spanning syndrome support `S subset F_2^w minus {0}`, the shifted
deficit

```math
F(S)=\rho(S)-1                                    \tag{HRC.15}
```

is nonnegative, antitone under union, and zero at the support containing all
nonzero syndromes.  Work on the upward-closed semigroup of spanning supports,
with the empty support allowed as a future identity context (equivalently,
adjoin a formal identity only for testing the present response).  The proof
of Corollary HRC.2 then applies verbatim.  Therefore all supports of radius at
most `k+1` form a response ideal and may be collapsed, with the shifted
deficit decoded as `k/2`, at all-future error `k/2`.  This terminal cell is
maximal at that error.  This is much
weaker quantitatively than the landmark quotient: it retains every
large-radius support and gives no count of that hard core.  Its role is to
separate terminal absorption from the landmark theorem's genuinely stronger
geometric compression of profiles that remain far from terminal.

## 5. What is classical and what is learned here

**Classical ingredients.** Rees quotients collapse semigroup ideals; future-
context equivalence is a syntactic congruence; matroid subsets quotient by
closure and flats compose by join.

**Extremal-response theorem.** Translation contraction makes every terminal
response ball an ideal, so its Rees quotient has a uniform, adversarial
all-future distortion bound.  For a nonnegative antitone deficit, midpoint
decoding makes the full `2eta` sublevel set the largest possible terminal
cell at error `eta`, eliminating the general factor-two gap.  This supplies
a feature-algebra stopping rule: preserve the hard core exactly or by a
separate geometric sketch, and stop retaining features once an absorbing
response interval is reached.

This principle is generative but limited.  It produces a third strict
composable quotient and distinguishes absorption from landmark compression.
It does not explain the size of a model's hard core.  The strongest next
question is therefore:

> Under which algebraic or geometric hypotheses does the number (or metric
> entropy) of response-hard-core classes grow sub-landscape under repeated
> composition?

The syndrome landmark theorem answers this through a Lipschitz chart even
inside the hard core; the arbitrary-table semilattice proves that some such
additional hypothesis is necessary.
