# Near-minimizer rigidity does not narrow all-order recovery

Date: 2026-08-16.

Status: focused negative audit. No candidate below satisfies all three required
clauses in a convergence-relevant sense. Conference matrices are used only as
test objects and falsifiers; no conference structure is assumed for extremal
signings.

## 1. Test being applied

Use the repository normalization

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j=\frac12x^{\mathsf T}Ax,
\qquad Q(A)=\max_x|H_A(x)|,
\qquad L=\liminf_n\frac{M_n}{n^{3/2}}.
```

A useful rigidity property `P` must pass all three tests below.

1. **Uniform forcing.**  There must be a modulus tending to zero with the
   near-optimality tolerance, or a canonical selection from the proved
   purified near-minimizers, which forces `P`.  Merely passing to a further
   subsequence and naming whatever statistic converges is not rigidity.
2. **Limit closure.**  `P` must survive in the selected action/extremal-envelope
   limit in the universal direction needed for an upper bound.  An empirical
   moment which permits an exceptional Boolean coloring is not closed in this
   sense.
3. **Objective-safe all-order realization.**  From the resulting extremal
   `P`-object one must construct exact signings on upward ratio-dense orders
   with objective at most the extremal value plus `o(1)`, without selecting
   target-order minimizers.  The existence of unrelated all-order examples
   satisfying `P` is not enough.

The third interpretation is compulsory.  If “realization” meant only that
some matrices with `P` exist at every order, a generic consequence of the
operator bound would pass while saying nothing about the selected extremal
envelope.

Exactly three candidates were frozen before the targeted ledger search:

1. two-walk residuals, encompassing spectral flatness, signed cycle
   deviations, and conference identities;
2. spectral/Boolean delocalization and row-field uniform integrability; and
3. switched ground-state row-field laws.

## 2. Candidate I: spectral flatness and cycle deviations

### 2.1 Weak property and exact theorem

The natural scalar defect is

```math
\Delta_4(A)
:=\frac{\|A^2-(n-1)I\|_F^2}{n^3}
=\frac{\operatorname{tr}(A^4)-n(n-1)^2}{n^3}\ge0.       \tag{2.1}
```

It is simultaneously the excess normalized fourth spectral moment, the
aggregate squared two-walk deviation, and an aggregate signed four-cycle
statistic.  Equality before normalization is exactly the conference identity
`A^2=(n-1)I`.

There is a genuine all-order construction for the *weak* property.  Let `C`
be a symmetric conference signing of order `m=n+h`, and let `A` be an
`n`-vertex principal submatrix.  With the corresponding off-diagonal block
`B`,

```math
A^2-(n-1)I=hI-BB^{\mathsf T}.                            \tag{2.2}
```

Since `rank(BB^T)<=h`, `tr(BB^T)=nh`, and
`||BB^T||_op<=m-1`,

```math
\Delta_4(A)=O(h/n+h^2/n^2).                              \tag{2.3}
```

Paley conference orders and the prime number theorem in the relevant
progression permit `h=o(n)`.  Hence exact hollow signs with
`Delta_4=o(1)` exist at all orders.  Principal deletion also gives

```math
\frac{Q(A)}{n^{3/2}}\le\frac12+o(1).                     \tag{2.4}
```

This is an exact positive theorem, but only for the already known comparison
constant `1/2`.

The missing forcing theorem would have to say, for some
`omega(eta)->0`, that a selected purified sequence with

```math
Q(A_n)\le(L+\eta)n^{3/2},\qquad
\|A_n\|_{op}\le C_\eta\sqrt n
```

satisfies

```math
\limsup_n\Delta_4(A_n)\le\omega(\eta).                  \tag{NMR4}
```

No such theorem is known.  The orientation-even `A^2` gain in the ledger
detects `Delta_4`, but its joint-selection/heavy-field loss does not turn
near-minimality into `(NMR4)`.  Also, if `L<1/2`, the conference realizers in
(2.4) do not realize the selected extremal value.  If `L=1/2`, convergence is
already immediate from the known universal upper bound.

### 2.2 Exact falsifier to profile-level closure

Weak fourth-moment flatness can hide a complete Boolean spectral spike.  Take
a symmetric conference matrix of order `n+1` and split off one vertex:

```math
C=\begin{pmatrix}0&b^{\mathsf T}\\ b&A\end{pmatrix},
\qquad b\in\{\pm1\}^n,
\qquad C^2=nI.
```

Exact multiplication gives

```math
Ab=0,
\qquad A^2-(n-1)I=I-bb^{\mathsf T},
\qquad
\Delta_4(A)=\frac{n(n-1)}{n^3}\longrightarrow0.         \tag{2.5}
```

Nevertheless `b` is itself a Boolean coloring, its Boolean spectral measure
is `delta_0`, and, in probability-space `L_2` normalization,

```math
\left\|\frac{A^2-(n-1)I}{n}b\right\|_{L_2}
=\frac{n-1}{n}\longrightarrow1.                          \tag{2.6}
```

Thus (2.1) is an empirical statement, not a universal profile statement.  It
does not exclude even one outer Boolean law far from the conference relation.
Strengthening it to uniform Boolean two-walk control would repair this defect,
but then even one-vertex conference deletion fails the property.  That
strengthening is therefore not stable under the lossless near-order deletion
operation.  Deletion may still be applied after a scalar upper bound is
obtained, but the strengthened property cannot itself be the closed
all-order class supplied by that operation.

Bounded cycle enrichment does not repair the problem.  The archived planting
calculation flips `h=floor(delta n^(3/2))` edges against one selected Boolean
spin.  Its energy rises by `2h`, while for a uniform spin `X`,

```math
\mathbb E\bigl(H_B(X)-H_A(X)\bigr)^2=4h.                 \tag{2.7}
```

Consequently every fixed-replica normalized energy law and every fixed signed
subgraph/cycle density changes by `o(1)`, although the planted maximum changes
at leading order.  Growing the cycle hierarchy to the scale which detects
this spike reopens the signed Eulerian/coset histogram.

### 2.3 Verdict

The weak property is all-order realizable but neither forced nor sufficient
for envelope recovery.  The strong universal property is profile-relevant
but fails the exact deletion test and has no independent realization theorem.
This candidate collides with the archived bounded-cycle and Boolean-resonance
obstructions; it is not a strict reduction.

## 3. Candidate II: delocalization and row-field tails

For a Boolean spin `x`, orient its local fields by

```math
\ell_i=x_i(Ax)_i.
```

If `x` is a positively oriented absolute ground state, one-flip optimality
gives `ell_i>=0`.  A natural strong delocalization condition is square-field
uniform integrability,

```math
\lim_{K\to\infty}\limsup_n
\frac1{n^2}\sum_i\ell_i^2
\mathbf1_{\{|\ell_i|>K\sqrt n\}}=0,                     \tag{3.1}
```

either for ground states or uniformly for all Boolean spins.

### 3.1 Exact near-minimizer falsifier

Condition (3.1) is not forced even by a fixed normalized operator bound.
Start from any bounded-operator purified near-minimizer `B_m`, switch and
negate it so that `1` is a positive ground state, and add one universally
positive vertex:

```math
\widetilde B_m=
\begin{pmatrix}0&\mathbf1^{\mathsf T}\\
                \mathbf1&B_m\end{pmatrix}.              \tag{3.2}
```

Triangle inequality and the all-one witness give the exact identities and
bounds

```math
Q(\widetilde B_m)=Q(B_m)+m,
\qquad
\|\widetilde B_m\|_{op}
\le\|B_m\|_{op}+\sqrt m.                                \tag{3.3}
```

Thus (3.2) has the same limiting normalized objective and still has
`O(sqrt(m))` operator norm.  At its all-one ground state, however, the new
field equals `m`.  For every fixed `K`, its single contribution to (3.1) is

```math
\frac{m^2}{(m+1)^2}\longrightarrow1.                    \tag{3.4}
```

This is a scalable exact falsifier to ground-state square-field rigidity.  It
also explains the topology: the exceptional vertex has empirical mass
`1/(m+1)` and disappears from the weak one-profile, while its square moment
stays macroscopic.  Adding square-moment uniform integrability would make that
moment closed, but (3.2) proves that near-minimality does not force it.

There is a weak theorem, but it is already contained in purification.  If
`||A||_op<=C sqrt(n)`, then uniformly in every Boolean `x`,

```math
\frac1{n^{3/2}}\sum_i |\ell_i|
 \mathbf1_{\{|\ell_i|>K\sqrt n\}}
\le\frac1{Kn^2}\sum_i\ell_i^2
\le\frac{C^2}{K}.                                       \tag{3.5}
```

This linear-tail bound supplies the uniform integrability already used in
action continuity.  It does not control `A^2`-energy or restrict the outer
Boolean profile, so it gives no new realization state.

### 3.2 Exact falsifiers to ordinary spectral delocalization

Even the standard strongest coordinate-isotropy surrogates do not control
cube alignment.  For a symmetric conference signing `C` of order `N`, set

```math
P_+=\frac12\left(I+\frac C{\sqrt{N-1}}\right).
```

This projection has constant diagonal `1/2` and all off-diagonal magnitudes
`1/(2 sqrt(N-1))`.  Nevertheless the archived Paley square-wave theorem gives
an infinite subsequence and Boolean vectors `x` with

```math
\frac{x^{\mathsf T}P_+x}{N}\longrightarrow1.             \tag{3.6}
```

Thus flat eigenbases, small projector entries, and equiangularity do not
exclude exponentially rare Boolean resonance.

There is also an elementary pointwise field spike in every symmetric
conference matrix.  For a root `r`, let

```math
x=e_r+Ce_r\in\{\pm1\}^N.
```

Since `C^2=(N-1)I`, its oriented local fields are exactly

```math
\{N-1,1,\ldots,1\}.                                      \tag{3.7}
```

There is one such coloring rooted at every coordinate.  Hence no common
`o(N)` exceptional vertex set can make (3.1) hold uniformly over all Boolean
colorings, even in the exactly flat test family.  This is a falsifier to a
proposed implication from conference flatness, not an assumption that
near-minimizers are conference matrices.

The only delocalization strong enough to avoid (3.6) is uniform
anti-alignment with all `2^n` cube points.  For a flat two-eigenspace model,
that condition directly bounds the same Boolean maximum being recovered.
Constructing it at new orders is therefore the universal outer-profile
obligation, not an independently realizable spectral property.

### 3.3 Verdict

Strong field or cube delocalization fails uniform forcing; ordinary
coordinate delocalization is constructible but does not control Boolean
resonance; the forced linear-tail version is vacuous beyond the existing
operator bound.  No version passes all three tests, and none is a strict
reduction.

## 4. Candidate III: ground-state row-field laws

For a positively oriented ground state `x`, define its oriented empirical
field law

```math
\lambda_{A,x}
=\frac1n\sum_{i=1}^n
\delta_{x_i(Ax)_i/\sqrt n}.                              \tag{4.1}
```

Under `||A||_op<=C sqrt(n)`, exact identities give

```math
\operatorname{supp}\lambda_{A,x}\subset[0,\infty),
\qquad
\int u\,d\lambda_{A,x}(u)=\frac{2Q(A)}{n^{3/2}},
\qquad
\int u^2\,d\lambda_{A,x}(u)\le C^2.                    \tag{4.2}
```

Consequently purified ground-state laws are subsequentially compact, and the
linear-tail estimate (3.5) preserves their means.  Limit points are extremal
one-profile laws.  This is the only candidate with a clean compactness
statement, but compactness supplies no unique law and no rigidity modulus.

After switching `x` to `1`, (4.1) is just the empirical signed row-sum law.
Approximate realization of an admissible row-sum distribution is a plausible
signed-degree construction.  It enforces only one-flip stability, not global
optimality.

### 4.1 Exact row-law falsifier

The following two order-six signings have exactly the same labelled row-sum
vector `(5,3,3,1,1,1)`:

```math
A=\begin{pmatrix}
0&1&1&1&1&1\\
1&0&-1&1&1&1\\
1&-1&0&1&1&1\\
1&1&1&0&-1&-1\\
1&1&1&-1&0&-1\\
1&1&1&-1&-1&0
\end{pmatrix},
\qquad
B=\begin{pmatrix}
0&1&1&1&1&1\\
1&0&1&-1&1&1\\
1&1&0&1&-1&1\\
1&-1&1&0&1&-1\\
1&1&-1&1&0&-1\\
1&1&1&-1&-1&0
\end{pmatrix}.                                           \tag{4.3}
```

All row sums are positive, so the all-one spin is a strict one-flip local
maximum for both, with energy `7`.  Exact enumeration of the 32 spins modulo
global negation gives

```text
A: -11:1, -9:1, -5:2, -3:4, -1:6, 1:8, 3:6, 5:3, 7:1;
B:  -7:2, -5:4, -3:4, -1:6, 1:6, 3:4, 5:4, 7:2.
```

Hence

```math
Q(A)=11,
\qquad Q(B)=7.                                           \tag{4.4}
```

The complete labelled row vector, not merely its moments or empirical law,
therefore fails to determine whether its locally stable state is an absolute
ground state.

### 4.2 Exact circularity boundary

Suppose an all-order row-law theorem realizes (4.1) at a spin `x_m` and also
asserts

```math
|H_{A_m}(y)|\le H_{A_m}(x_m)
\quad\hbox{for every }y\in\{\pm1\}^m.                   \tag{4.5}
```

Then (4.2) immediately gives the desired scalar upper bound.  But (4.5) is
exactly the full target-order Boolean optimization constraint.  The row-law
construction has not reduced it.

Alternatively, retain the row-field law for every Boolean spin.  For

```math
\lambda_{A,y}^{\rm full}
=\mathcal L_i\left(y_i,\frac{(Ay)_i}{\sqrt n}\right),
```

one has

```math
\frac{H_A(y)}{n^{3/2}}
=\frac12\int uv\,d\lambda_{A,y}^{\rm full}(u,v).         \tag{4.6}
```

Thus the outer set of all these laws determines `Q(A)` exactly.  Recovering
that set in the directed direction is the Boolean part of one-profile
recovery, with the same universal `2^n` quantifier.  Fixed row moments or a
bounded collection of field laws fall under the archived bounded-local-data
obstruction; the complete collection restores the missing profile.

### 4.3 Verdict

One selected ground-state law is compact and locally realizable but does not
certify the ground state.  Adding globality is scalar recovery itself; adding
all Boolean row laws is directed profile recovery in different notation.
This candidate is not a strict reduction.

## 5. Severe conclusion

| candidate | uniformly forced | envelope-closed at useful strength | objective-safe all-order signs | result |
|---|---|---|---|---|
| fourth-moment/cycle flatness | open | weak version misses Boolean spikes | only unrelated `1/2` conference models | fail |
| delocalization/field tails | strong version false by (3.2) | weak `L_1` version is already automatic | cube-uniform version is the outer profile | fail |
| ground-state row law | only subsequential compactness is forced | yes for a selected law | global realization is (4.5) | fail |

No tested property passes all three clauses.  The weak properties are genuine
information quotients but lose the rare Boolean colorings that determine the
objective.  Every strengthening which retains those colorings either fails an
exact near-minimizer/conference-deletion test or reconstructs the directed
Boolean one-profile.  Near-minimizer rigidity therefore supplies no strict
narrowing of `EER`, `AR_min^->`, or sign-near weighted recovery at this
checkpoint.

The result is a stop, not a new route: spectral flatness collides with
arithmetic resonance and bounded-cycle blindness; delocalization collides
with localized square-field mass and cube alignment; row-field laws collide
with the full all-coset extreme tail.
