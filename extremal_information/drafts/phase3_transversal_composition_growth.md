# Dense transversals: a bounded-depth quotient and a sharp mixed-cycle obstruction

**Status.**  The statements below are proved.  The independent finite checks
in
[`verify_phase3_transversal_composition.py`](../experiments/verify_phase3_transversal_composition.py)
test the uniform synchronization refinement, the affine-state profile bound,
the mixed-cycle identity, and the equal-depth obstruction.

This note answers the composition question left open by the geodesic
synchronization theorem.  A single complete, cycle-contracting transversal
is within constant all-context distance of a linear graph.  A bounded union
of exact linear graphs has a closed affine feature state.  These two facts
give an `O(ell)`-error decoder for `ell` general sources.  For exact linear
sources the sharper affine-circuit-rank law is independent of `ell`, but
its error must grow with the rank: mixed cycles, absent inside every atomic
source, are created by union.  Already for a one-dimensional quotient,
`Theta(D)` linear sources with the same affine state and rank can have
covering radii separated by `D/2`.

## 1. Setup

Let

```math
G=W\oplus Q,\qquad W=\mathbb F_2^D,\qquad Q=\mathbb F_2^k,
```

and fix the coordinate basis
`B={(e_i,0):1<=i<=D}` of `W`.  Write `|.|_B` for its Hamming
weight.  For a map `f:Q to W` with `f(0)=0`, put

```math
S_f=B\cup\{(f(q),q):q\in Q\setminus\{0\}\}.       \tag{TC.1}
```

We call `f` **cycle-contracting** when

```math
\left|\sum_{q\in R}f(q)\right|_B\le |R|          \tag{TC.2}
```

for every `R subseteq Q\{0}` whose binary sum is zero.  This is precisely
the complete-transversal specialization of the geodesic cycle criterion.

For a spanning support `S subseteq G\{0}`, let `lambda_S(x)` be Cayley word
length and `rho(S)=max_x lambda_S(x)`.  Its complete future-response distance
from `T` is

```math
d_{\rm resp}(S,T)=
\sup_{E\subseteq G\setminus\{0\}}
|\rho(S\cup E)-\rho(T\cup E)|.                  \tag{TC.3}
```

All claims below are actually pointwise statements about `lambda`, and are
therefore stronger than (TC.3).

## 2. A complete transversal has a uniformly close linear centre

The joint BLR theorem gives average distance at most three.  Completeness of
the quotient and the same triangle constraints upgrade this to a uniform
constant.

### Theorem TC.1 (uniform synchronization and all-context replacement)

If `f:Q to W` is cycle-contracting, there is a linear map `L:Q to W` such
that

```math
\max_{q\in Q}|f(q)+L(q)|_B\le8.                 \tag{TC.4}
```

For this same `L`, every raw future `E subseteq G\{0}` and every `x in G`
satisfy

```math
\lambda_{S_f\cup E}(x)
\le\lambda_{S_L\cup E}(x)+8,
\qquad
\lambda_{S_L\cup E}(x)
\le\lambda_{S_f\cup E}(x)+10.                 \tag{TC.5}
```

In particular `d_resp(S_f,S_L)<=10`.  If `f_1,...,f_ell` are
cycle-contracting and `L_i` are the corresponding linear maps, then

```math
\sup_{E,x}
\left|
\lambda_{(\bigcup_iS_{f_i})\cup E}(x)
-\lambda_{(\bigcup_iS_{L_i})\cup E}(x)
\right|
\le10\ell.                                     \tag{TC.6}
```

#### Proof

Put `N=|Q|`.  When `k=0` there is nothing to prove, and when `k=1` every
map with value zero at zero is linear.  Assume `k>=2`.  For

```math
\partial f(a,b)=f(a)+f(b)+f(a+b),
```

the defect is zero on degenerate pairs and has Hamming weight at most three
on the `(N-1)(N-2)` ordered nondegenerate pairs, by (TC.2).  The joint BLR
argument therefore gives a linear `L` such that, with `e=f+L` and
`A=sum_q |e(q)|_B`,

```math
A\le {3(N-1)(N-2)\over N}.                     \tag{TC.7}
```

Fix nonzero `q`.  For every `a notin {0,q}`, the set
`{q,a,q+a}` is a three-element zero-sum cycle, hence

```math
|e(q)|_B\le3+|e(a)|_B+|e(q+a)|_B.
```

Summing over those `N-2` values of `a` gives

```math
N|e(q)|_B
\le3(N-2)+2A
\le3(N-2)+{6(N-1)(N-2)\over N}<9N.            \tag{TC.8}
```

The weight is integral, so it is at most eight.  This proves (TC.4).

We next record a consolidation consequence of (TC.2).  If `R` is any set
of graph generators with quotient sum `q`, then:

* for `q=0`, its `W`-sum costs at most `|R|` basis generators;
* for `q ne 0`, compare `R` with the graph generator over `q`.  Their
  symmetric difference is a zero-sum cycle of size at most `|R|+1`.
  Thus all of `R` can be replaced by that one graph generator and at most
  `|R|+1` basis generators, an overhead of at most two.

Start with a shortest word over `S_L union E`.  All `S_L` graph letters
combine exactly to zero or to one graph letter.  Replacing the latter by
the corresponding `S_f` letter and at most eight basis letters proves the
first inequality in (TC.5).  Conversely, consolidate all `S_f` graph
letters at overhead at most two and then replace the surviving letter by
its `S_L` counterpart at cost at most eight.  This proves the second
inequality.

For (TC.6), replace the sources one at a time.  At each step, the arbitrary
future in (TC.5) is the union of `E` with every other raw or already
replaced source.  Pointwise telescoping proves the claim. `square`

The constants eight and ten are certified bounds, not claimed sharp.  The
important feature is that they do not depend on `D` or `k` and use the
joint Hamming budget before splitting output coordinates.

## 3. The affine-hull state for a bounded union of linear graphs

Let `L_1,...,L_ell in Hom(Q,W)` and put

```math
V=\sum_{i=2}^{\ell}\operatorname{Im}(L_i+L_1). \tag{TC.9}
```

The reference map is well defined modulo `V`.  Write the resulting state as
`[L_1,V]`.  For `(w,q) in W direct-sum Q`, define

```math
\delta_{[L,V]}(w,q)
=d_B(w+L(q),V)+\mathbf1_{q\ne0},               \tag{TC.10}
```

where `d_B(z,V)=min_{v in V}|z+v|_B`.

### Theorem TC.2 (affine quotient algebra and bounded-depth decoder)

Let `U=bigcup_(i=1)^ell S_(L_i)`.  Then

```math
\delta_{[L_1,V]}(x)
\le\lambda_U(x)
\le\delta_{[L_1,V]}(x)+\ell                    \tag{TC.11}
```

for every `x in G`.  Consequently, for every future support `E`,

```math
\max_x(\delta_{[L_1,V]}\mathbin\square\lambda_E)(x)
```

approximates `rho(U union E)` to additive error at most `ell`, where
`square` is min-plus convolution.

These feature states have an exact commutative idempotent update

```math
[L,V]\odot[K,Z]
=
[L,\,V+Z+\operatorname{Im}(L+K)].              \tag{TC.12}
```

It is well defined, and multiplying the atomic states `[L_i,0]` produces
exactly (TC.9).  This is an exact update of the **approximate feature
state**, not an exact response algebra.  A state can be encoded using at
most `D^2+Dk` bits, independent of `ell`.  Combining TC.1 and TC.2, a union of `ell`
cycle-contracting transversals has an all-future decoder from this state
with error at most

```math
11\ell.                                         \tag{TC.13}
```

#### Proof

Modulo `V`, every graph generator `(L_i(q),q)` has the common image
`(L_1(q)+V,q)`.  Projecting any word for `x` therefore gives a word in the
quotient generated by the images of `B` and this one complete graph.  Its
word length is exactly (TC.10), proving the lower bound.

For the upper bound, choose `v in V` attaining
`d_B(w+L_1(q),V)` and write

```math
v=\sum_{i=2}^{\ell}(L_i+L_1)(q_i).
```

Set `q_1=q+sum_(i=2)^ell q_i`.  Omitting zero labels, the at most `ell`
graph generators `(L_i(q_i),q_i)` sum to `(L_1(q)+v,q)`.  Correct the
remaining `W`-part with `d_B(w+L_1(q),V)` basis generators.  This proves
(TC.11).

If `L` is changed by a map with image in `V`, or `K` by one with image in
`Z`, the subspace on the right of (TC.12) does not change.  Thus the product
is well defined.  Commutativity, associativity, and idempotence also follow
directly from the reference-free description

```math
V({\cal L})=\sum_{L,K\in{\cal L}}\operatorname{Im}(L+K),
```

with any member of the family as reference.  Storing a basis of `V` and a
matrix for `L` gives the crude advertised bit bound.

Finally, min-plus convolution and taking a maximum are nonexpansive in
uniform norm.  Apply TC.1 to replace the `ell` raw transversals at cost
`10ell`, then (TC.11) at cost `ell`. `square`

This is a bounded closed state: it stores a subspace and one quotient map,
not a truth table on `Q` or the union of all graph supports.  It is not
automatically a shorter encoding than the raw list when `ell` is small;
its advantage is closure and a size independent of later composition
depth.  Its response decoder remains approximate.  The next sections
identify the sharper rank parameter and show why some growing error is
necessary.

## 4. The exact composition defect is mixed-cycle excess

The zero-test in the geodesic criterion has a quantitative form.

### Theorem TC.3 (mixed-cycle defect identity)

Let `S subseteq G\{0}` contain the independent basis `B`, put
`t=sum_(i=1)^D(e_i,0)`, and let `pi:G to G/W`.  Define

```math
\Delta_B(S)=
\max_{\substack{R\subseteq S\setminus B\\
                 \sum_{s\in R}\pi(s)=0}}
\left(
\left|\sum_{s\in R}s\right|_B-|R|
\right),                                        \tag{TC.14}
```

including the empty set.  Then

```math
\boxed{\lambda_S(t)=D-\Delta_B(S).}             \tag{TC.15}
```

If each `S_i` is cycle-contracting, every positive term in
`Delta_B(bigcup_i S_i)` is necessarily a mixed cycle.  In particular, two
graph generators from different sources with the same nonzero quotient
label `q` contribute

```math
|(L_i+L_j)(q)|_B-2.                             \tag{TC.16}
```

Thus composition creates exactly the obstruction that each source's
internal hypotheses omit.

#### Proof

In any word for `t`, let `R` be the letters outside `B`.  Its projected sum
is zero.  If its `W`-sum has support `J subseteq[D]`, the uniquely required
basis correction is the complement of `J`, of size `D-|J|`.  The word
therefore has length

```math
|R|+D-\left|\sum_{s\in R}s\right|_B.
```

Minimizing over `R` proves (TC.15).  The remaining assertions follow from
the definition and the two-letter cycle. `square`

For `dim Q=1`, write a linear map as `L_a(1)=a`, `a in W`, and let `A`
be the set of maps already unioned.  Then the union is cycle-contracting if
and only if

```math
\operatorname{diam}_{\rm Ham}(A)\le2.          \tag{TC.17}
```

Necessity comes from two-letter cycles.  Conversely, every zero-quotient
cycle has even size; pair its elements arbitrarily and use the triangle
inequality.  This is a concrete closure criterion: pairwise mixed-cycle
control is sufficient in quotient dimension one, but affine span alone is
not.

More quantitatively, in this case

```math
\Delta_B(U_A)=
\max_{\substack{R\subseteq A\\|R|\text{ even}}}
\left(\left|\sum_{a\in R}a\right|_B-|R|\right), \tag{TC.18}
```

so a pair at Hamming distance `r` already forces
`lambda_(U_A)(t)<=D-r+2`.  Also, for `D>=3`, a diameter-two binary
anticode has at most `D+1` points.  Indeed, after translation all words
have weight at most two; the weight-two supports are pairwise-intersecting
edges and hence form a star or a triangle.  Counting the compatible
singletons gives at most `D+1`.  Therefore `D+1` is the largest possible
number of distinct one-quotient linear sources that can compose while
preserving the original geodesic.  Composition beyond this threshold must
create a positive mixed-cycle excess.

## 5. Equal-depth affine-span obstruction

The `O(ell)` error in TC.2 has the correct order under only affine-state
information.

### Theorem TC.4 (same affine state, linearly separated responses)

Let `D>=6` be even, let `Q=F_2`, and put `t=e_1+...+e_D`.  Consider the two
equal-size collections of `D+1` linear maps indexed by

```math
A_0=\{0,e_1,\ldots,e_D\},
\qquad
A_1=\{0,t+e_1,\ldots,t+e_D\}.                  \tag{TC.19}
```

Both affine spans are all of `W`, so both have exactly the same state
`[0,W]` in TC.2.  Nevertheless, for

```math
U_A=B\cup\{(a,1):a\in A\},
```

one has

```math
\rho(U_{A_0})=D,
\qquad
\rho(U_{A_1})={D\over2}.                       \tag{TC.20}
```

Thus even with source count retained, every decoder using only affine
span/rank data has worst-case error at least `D/4`.  Moreover the complete
rooted word profile recovers the exact set:

```math
A=\{a\in W:\lambda_{U_A}(a,1)=1\}.             \tag{TC.21}
```

Consequently arbitrary-depth exact composition of atomic graph maps is the
full support-union semilattice.  After at most `ell` distinct atoms it has
`sum_(j<=ell) binom(2^D,j)` possible states.

Exact recovery in (TC.21) is elementary (and unrestricted all-future radius
queries can likewise test one missing support element).  The new no-go is
not that exact support information is expensive; it is the macroscopic
empty-future separation (TC.20) between equal-depth families after the
affine state has identified them.

The radii in (TC.20) are the **actual** Cayley radii.  They must not be
confused with the affine decoder profile: for the common state `[0,W]`,
(TC.10) is simply `delta(w,q)=1_(q ne 0)` and its maximum is one.  Theorem
TC.2 permits the actual profile to sit as much as `ell=D+1` above this
collapsed lower profile.  The present example proves that a linear portion
of that allowance is unavoidable.

#### Proof

Both collections contain zero.  The `e_i` form a basis.  For even `D`, the
vectors `t+e_i` also form a basis: if their coefficient sum is `s`, every
coordinate equation says `c_i=s`, and then `s=D s=0` in `F_2`.
Hence both affine states are `[0,W]`.

For `A_0`, the target `(t,1)` costs exactly `D`.  Indeed, if a word uses
the zero lift with indicator `z` and basis lifts indexed by `I`, its cost
after the unique basis correction is

```math
|I|+z+|t+\sum_{i\in I}e_i|=D+z\ge D.
```

Every `(w,0)` costs at most `D` using `B`; every `(w,1)` of weight below
`D` uses the zero lift and at most `D-1` basis letters, while `(t,1)` uses
one basis lift and `D-1` basis letters.  Thus the radius is `D`.

Write `D=2m`.  A word over `U_(A_1)` uses a subset `I` of the nonzero
lifts and possibly the zero lift.  If `r=|I|`, its `W`-sum is

```math
(r\bmod2)t+\sum_{i\in I}e_i.                   \tag{TC.22}
```

Choose any `w` of weight `m`.  The displayed vector has weight `r` when
`r` is even and `D-r` when `r` is odd.  The reverse Hamming triangle
inequality shows that the number of lift and correction letters is at least
`m`, for either parity.  Hence the radius is at least `m`.

For the reverse bound, let `s=|w|<D`.  Using no nonzero lift, or using one
`t+e_i` with `i` outside the support of `w`, gives

```math
\lambda_{U_{A_1}}(w,0)\le\min\{s,D-s+1\},
\qquad
\lambda_{U_{A_1}}(w,1)\le\min\{s+1,D-s\}.     \tag{TC.23}
```

Both minima are at most `m`.  At `w=t`, one nonzero lift leaves one basis
coordinate and costs two for quotient one; adding the zero lift costs three
for quotient zero.  Since `m>=3`, these endpoints also cost at most `m`.
This proves (TC.20).  Finally, a length-one word with quotient coordinate
one is exactly one of the lift generators, proving (TC.21). `square`

The obstruction is not a failure of either source's synchronization:
every atomic source is already an exact linear graph.  It is the creation
of mixed cycles under composition.  `A_0` is a diameter-two anticode and
keeps `Delta_B=0`; in `A_1`, the pair `0,t+e_i` has excess `D-3` and gives
`lambda(t)=3`.  The macroscopic radius separation is the global effect of
these new low-cost cross-source directions.

## 6. Circuit rank, rather than source count, controls the quotient error

The preceding `ell`-bound pays every source.  Linear dependencies among the
maps make most of that payment unnecessary.  Let

```math
{\cal A}=L+{\cal U}\subseteq\operatorname{Hom}(Q,W) \tag{TC.24}
```

be the affine hull of a nonempty family of linear maps, where `cal U` is a
linear subspace of map space and

```math
r=\dim {\cal U},\qquad
V({\cal U})=\sum_{A\in{\cal U}}\operatorname{Im}A. \tag{TC.25}
```

The state is the affine subspace `[L;cal U]`, not merely its pointwise image
span.  Define its lower profile by

```math
\delta_{[L;{\cal U}]}(w,q)
=d_B(w+L(q),V({\cal U}))+\mathbf1_{q\ne0}.      \tag{TC.26}
```

### Theorem TC.5 (affine-circuit-rank response law)

Let `F` be any nonempty set of linear maps with affine hull
`L+cal U`, and put `U_F=bigcup_(K in F)S_K`.  Then, independently of
`|F|`,

```math
\delta_{[L;{\cal U}]}(x)
\le\lambda_{U_F}(x)
\le\delta_{[L;{\cal U}]}(x)+r+1.               \tag{TC.27}
```

Therefore the midpoint profile answers every appended-support
covering-radius query to error at most `(r+1)/2`.

Affine-subspace feature states have the exact join law

```math
[L;{\cal U}]\vee[K;{\cal Z}]
=
[L;{\cal U}+{\cal Z}
       +\operatorname{span}_{\mathbb F_2}\{L+K\}]. \tag{TC.28}
```

This is again an exact update of an approximate response state, not an
exact response algebra.  The law is well defined, commutative, associative,
and idempotent.  If
`m=Dk=dim Hom(Q,W)`, a rank-`r` state uses at most `(r+1)m` bits by storing
one point and a basis.  In particular, arbitrary-depth unions whose map
affine rank is `o(D)` have `o(D)` all-context distortion even if the number
of sources is arbitrarily large, **provided the sources are the exact
linear graph supports appearing in this theorem**.

The order in `r` is necessary.  The two families in TC.4 have the same
state of affine rank `r=D`, the same number of sources, and radii separated
by `D/2`.  Hence the minimax error of an affine-subspace-state decoder is at
least `r/4` in the worst case, while (TC.27) gives at most `(r+1)/2`.

#### Proof

Choose an affine basis `L_0,L_1,...,L_r` from `F`, with
`A_j=L_j+L_0` a basis of `cal U`.  Modulo `V(cal U)`, every graph generator
from every member of `F` becomes the common graph of `L_0`.  Projection
therefore gives the lower bound exactly as in TC.2.

Conversely, choose `v in V(cal U)` attaining the distance in (TC.26).
Because the `A_j` form a basis of map space `cal U`, there are
`q_1,...,q_r in Q` such that

```math
v=\sum_{j=1}^r A_j(q_j).                        \tag{TC.29}
```

Put `q_0=q+sum_(j=1)^r q_j`.  After omitting zero labels, the at most
`r+1` graph generators `(L_j(q_j),q_j)` have sum

```math
(L_0(q)+v,q).
```

The remaining kernel residual costs exactly
`d_B(w+L_0(q),V(cal U))` basis letters.  This proves (TC.27).  The true
profile lies pointwise in an interval of width `r+1`, so centering that
interval and using nonexpansiveness of min-plus convolution and maximum
proves the all-future midpoint claim.

Formula (TC.28) is the ordinary affine hull of a union.  Changing the
representatives by members of their direction spaces leaves the displayed
affine subspace unchanged, and the semilattice laws follow from set union
followed by affine closure.  The bit bound is the cost of a point and `r`
basis maps.  TC.4 supplies the lower bound. `square`

The theorem identifies the useful circuit statistic.  The nullity
`|F|-r-1` consists of redundant atoms that do not enlarge the error budget;
only independent affine directions are paid.  It is a strict closed
quotient of the full atom set, but it cannot resolve the weights of
independent mixed cycles once `r` is macroscopic.

This rank law must not be transferred silently to nonlinear sources.  For
arbitrary cycle-contracting transversals, TC.1 first replaces each raw map
by a synchronized linear centre and currently pays `10ell`.  Knowing that
those centres have affine rank `r` then adds only the TC.5 midpoint error,
but does not remove the accumulated synchronization charge.  Thus the
combined certified bound is `10ell+(r+1)/2`, not `O(r)`.

## 7. Pairwise mixed-cycle spectra do not close

Affine rank is coarse but rooted.  A different tempting state is the full
unrooted histogram of Hamming distances between pairs of atomic maps.  The
next example shows that even this weighted pair information misses a
linearly growing response difference.

### Proposition TC.6 (same pair-distance law, separated responses)

Regard the following decimal integers as six-bit binary words:

```math
C_0=\{0,3,9,10,53,54,60,63\},
\qquad
C_1=\{0,3,12,15,48,51,60,63\}.                 \tag{TC.30}
```

Both are binary linear `[6,3]` codes with weight enumerator

```math
1+3z^2+3z^4+z^6,                               \tag{TC.31}
```

but their covering radii are respectively two and three.  For a linear
code `C subseteq F_2^D`, define the union of one-quotient graph atoms

```math
U_C=B\cup\{(c,1):c\in C\}.                     \tag{TC.32}
```

Then

```math
\lambda_{U_C}(z,1)=1+d(z,C),
\qquad
\lambda_{U_C}(z,0)
=\min\{|z|,\,2+d(z,C\setminus\{0\})\},        \tag{TC.33}
```

and consequently

```math
\rho(C)+1\le\rho(U_C)\le\rho(C)+2.             \tag{TC.34}
```

Let `C_i^(oplus s)` be the `s`-fold direct sum.  The two atom families have
identical unrooted pair-distance histograms at every `s`, including every
pair-distance moment, while

```math
\rho(U_{C_1^{\oplus s}})
-\rho(U_{C_0^{\oplus s}})\ge s-1.              \tag{TC.35}
```

Thus no feature state determined only by the complete unrooted two-atom
distance distribution can have sublinear worst-case distortion under
composition.  This does **not** say that the two labeled distance matrices
are isometric; their higher circuit geometry is exactly what differs.

#### Proof

Closure under xor is immediate from the displayed lists (or from the bases
`{3,9,53}` and `{3,12,48}`).  Counting weights gives (TC.31).  Exhaustion
is unnecessary for the radius claim.  Parity-check bases may be chosen so
that the six syndrome columns are respectively

```math
(1,1,7,1,2,4)\qquad\hbox{and}\qquad(1,1,2,2,4,4), \tag{TC.36}
```

where the entries are three-bit integers.  In the first list, every syndrome
is a sum of at most two columns and syndrome `3`, for example, needs two.
In the second, syndrome `7` needs the three distinct basis values `1,2,4`,
while every syndrome needs at most three.  This proves covering radii two
and three.  The verifier cited at the start independently reconstructs
these values from all 64 ambient words.

An odd subset of `C` has xor in `C`, and every `c in C` is already available
as one lift.  Hence the least odd-lift cost followed by basis correction is
`1+d(z,C)`.  A nonempty even subset whose xor is zero is dominated by the
empty subset.  If its xor is a nonzero `c in C`, the two lifts `0,c` realize
the same xor at no greater cost.  This proves (TC.33).  Taking maxima gives the lower bound
in (TC.34).  For the upper bound, if a nearest codeword is zero, the
`|z|` branch is at most `rho(C)`; otherwise the second branch is at most
`rho(C)+2`.

Direct sums multiply weight enumerators and add covering radii.  In a
linear code every difference word occurs in exactly `|C|` ordered pairs,
so the unrooted pair-distance histogram is determined by the weight
enumerator.  Thus the histograms of the two `s`-fold powers agree, whereas
(TC.34) gives

```math
\rho(U_{C_1^{\oplus s}})\ge3s+1,
\qquad
\rho(U_{C_0^{\oplus s}})\le2s+2.
```

Subtracting proves (TC.35). `square`

This is an old covering-radius obstruction placed inside the new graph-map
composition model.  Its theoretical role is precise: pair weights detect
two-letter mixed-cycle severity but forget how those pairs assemble into
higher circuits.  A weighted mixed-cycle state that genuinely closes must
retain some rooted higher-circuit incidence, not just the pair spectrum.

## 8. Research judgment

This yields a sharp three-regime law, with a necessary distinction between
one nonlinear transversal and unions of exact linear graph sources.

1. **One source:** cycle contraction synchronizes it to a linear graph at
   constant all-context error.
2. **Exact linear sources of sublinear affine circuit rank:** the
   affine-subspace algebra is a strict polynomial-size state and gives
   `o(D)` all-context error, even for arbitrarily many sources.  This
   strictly improves the bounded-depth `ell` law in that linear class; it
   does not erase the per-source synchronization cost for nonlinear maps.
3. **Exact linear sources of linear affine circuit rank:** mixed cycles can
   change the extremal response by `Theta(D)` while affine hull, rank, and
   source count remain identical.  Exact composition retains an arbitrary
   subset of atomic maps.
4. **Unrooted pair weights:** even the complete pair-distance histogram
   fails by `Theta(D)` under direct-sum composition.  Higher rooted circuit
   incidence, rather than more pair moments, is the next missing feature.

The positive theorem is not merely convex duality or dynamic programming:
its input is the deterministic cycle-contraction law, its compression step
is joint vector-Hamming synchronization, and its product is an algebraic
quotient of a union of graph subspaces.  The negative theorem identifies
the exact feature that the quotient forgets and proves a linear response
penalty at equal composition depth.

The strongest next theorem is now precise: characterize a weighted
mixed-cycle profile, strictly smaller than the full set of maps, that
refines affine circuit rank when `r=Theta(D)`.  It should interpolate
between the exact `o(D)`-rank quotient of TC.5 and the full support
semilattice.  Any candidate must distinguish the two equal-depth families
in TC.4; affine hull or abstract rank data cannot.
