# Selective neutralization in finite-field dictionaries

**Status.** Proved draft.  This report extracts the reusable content of the
Grassmannian syndrome packing without promoting a tautological response-
packing axiom.  The theorem applies both to appended-fragment covering radius
and to worst-case sparse synthesis under dictionary augmentation.

## 1. Why the completely abstract statement is rejected

For sources `x`, queries `q`, and scalar responses `R_x(q)`, the implication

```math
\forall x\ne y\ \exists q:\quad
|R_x(q)-R_y(q)|>2\eta
\quad\Longrightarrow\quad
\text{every uniform-`eta` summary separates all sources}              \tag{SN.1}
```

is merely the triangle inequality.  Calling (SN.1) “selective
neutralization” would repackage metric entropy rather than add a theorem.

The content below is instead an algebraic response embedding: a moving dense
linear carrier, when appended as one legal dictionary query, converts
subspace injection distance into a scalar extremal-response gap.  Neither the
response separation nor its scale is assumed.

## 2. Sparse-synthesis response

Let `K=F_q`, let `G=K^w`, and let a **projective dictionary** `D` be a finite
set of one-dimensional subspaces of `G` whose span is `G`.  Its synthesis
length and worst-case synthesis radius are

```math
\ell_D(x)=\min\left\{k:
x=c_1v_1+\cdots+c_kv_k,\quad
c_i\in K^*,\ [v_i]\in D\right\},                 \tag{SN.2}
```

```math
\rho_0(D)=\max_{x\in G}\ell_D(x).                \tag{SN.3}
```

Repeated use of one projective atom is unnecessary: its two coefficients can
be combined.  A future dictionary `E` acts by union, and the declared query
response is

```math
R_D(E)=\rho_0(D\cup E).                           \tag{SN.4}
```

This is a standard exact sparse-representation quantity: it asks for the
largest support size needed to synthesize a target after future atoms are
made available.  For `q=2`, a parity-check matrix with column-type support
`D` has covering radius (SN.3), so (SN.4) is also the repository's syndrome
response.

## 3. The carrier theorem

For a subspace `W<=G`, write `P(W)` for its projective points.  Choose a
linear complement `V_W` and a basis `C_W` of `V_W`, regarded as projective
atoms.  Put

```math
D_W=P(W)\cup C_W.                                \tag{SN.5}
```

### Theorem SN.1 (selective neutralization by a moving carrier)

Let `W,W'` have the same dimension `1<=d<=w` in `K^w`, and put

```math
s=d-\dim(W\cap W')                               \tag{SN.6}
```

for their injection distance.  Then

```math
\rho_0(D_W)=w-d+1,                               \tag{SN.7}
```

```math
\rho_0(D_W\cup D_{W'})
\le w-d-s+2,                                     \tag{SN.8}
```

and the legal query `E=D_W` exposes

```math
\boxed{
|R_{D_W}(D_W)-R_{D_{W'}}(D_W)|\ge s-1.}          \tag{SN.9}
```

The complements and their bases may be chosen independently for every
carrier.

#### Proof

In `G=W direct-sum V_W`, the quotient component of a target has a unique
coordinate representation in `C_W`.  It costs its coordinate Hamming weight;
the remaining nonzero `W` component costs exactly one projective atom.
Maximizing gives (SN.7).

Put `L=W+W'`.  It has dimension `d+s`, hence codimension `w-d-s`.  The image
of `C_W` spans `G/L`, so at most `w-d-s` of those atoms match the quotient
of any target.  The residual lies in `W+W'` and is the sum of at most one
vector from each carrier.  This proves (SN.8).  Appending `D_W` to itself
does nothing, while appending it to `D_W'` gives their union.  Subtract
(SN.8) from (SN.7) to obtain (SN.9). `square`

The `-1` is a real endpoint allowance in this argument: a residual in
`W+W'` may require one atom from each carrier, whereas a residual in one
carrier requires one.  It is irrelevant at linear injection distance but is
why unit-separated carriers do not automatically produce a positive gap.

## 4. Information consequence over every finite field

Let

```math
C_q=\prod_{i=1}^{\infty}(1-q^{-i})^{-1}<\infty.  \tag{SN.10}
```

The Gaussian coefficient satisfies

```math
q^{k(n-k)}\le {n\brack k}_q
\le C_q q^{k(n-k)}.                              \tag{SN.11}
```

For fixed `d`-space `W<=K^{2d}`, the number of `d`-spaces at injection
distance exactly `j` is

```math
q^{j^2}{d\brack j}_q^2.                          \tag{SN.12}
```

Here (SN.12) uses `j=d-dim(W cap W')`; it is the usual Grassmann-scheme
intersection count.

### Corollary SN.2 (quadratic response packing)

For `w=2d` and `1<=r<=d`, there is a family of carrier dictionaries with
pairwise future-response distance at least `r-1` and size at least

```math
{q^{d^2-r(2d-r)}\over d C_q^2}.                  \tag{SN.13}
```

Hence, for every fixed `epsilon<1/4`, every deterministic summary answering
all future dictionary-augmentation responses within `epsilon*w` has, for all
sufficiently large even `w`, worst-case description length at least

```math
\left((1/2-\gamma)^2-o(1)\right)w^2\log_2q       \tag{SN.14}
```

for every fixed `2*epsilon<gamma<1/2`.

#### Proof

Greedily pack the Grassmannian while deleting all carriers at injection
distance below `r`.  By (SN.11)--(SN.12), one deleted ball has size at most

```math
d C_q^2 q^{r(2d-r)},
```

whereas the whole Grassmannian has at least `q^{d^2}` points.  This proves
(SN.13).  Theorem SN.1 supplies response separation `r-1`.  Sources sharing
one decoded summary have response distance at most `2*epsilon*w`; put
`r=ceil(gamma*w)` and take logarithms. `square`

For `q=2`, `C_2<4`, so (SN.13) specializes to the constants in the binary
Grassmannian report.

## 5. Microscopic versus macroscopic exposure

The theorem isolates a sharp scale condition for this architecture.

* A carrier pair at injection distance `s` can be separated by the displayed
  scalar response only to `s-1`.
* Uniform distortion `eta` is information-forcing precisely when the chosen
  carrier code has `s-1>2*eta`; below this threshold the common decoded value
  is not ruled out.
* Therefore bounded-distance carrier changes are **microscopic** at error
  `epsilon*w`, even though they can be detected exactly on the integer
  lattice.  A constant-relative-distance Grassmann code is required for
  macroscopic response information.

This distinction is not a separate-payment argument.  The query appends the
entire carrier at once, the quotient dimension falls before the maximum over
targets is taken, and one scalar radius records the coherent reduction.  The
number of appended atoms is exponential in `d`, but the response gap is only
the `d`-dimensional injection distance.  Thus atom count itself is not the
amplification parameter; **independent quotient directions neutralized** is.

## 6. Two model readings

1. **Linear codes.**  For `q=2`, dictionary atoms are nonzero syndrome column
   types, (SN.3) is coset-leader covering radius, and union is parity-check
   fragment concatenation.  Corollary SN.2 gives the quadratic-bit packing
   for arbitrary syndrome supports.

2. **Sparse representation.**  For arbitrary finite `q`, `D` is a synthesis
   dictionary, `ell_D(x)` is exact `ell_0` representation cost, and a query
   supplies future atoms.  The theorem says that any data structure which
   predicts worst-case sparse synthesis after every dictionary augmentation
   to additive error `epsilon*w` requires `Omega(w^2 log q)` bits on a
   concrete family.  No codeword or Hamming ambient space is used in this
   formulation.

The two readings are algebraically equivalent, but operationally distinct:
one answers future code covering radii, the other answers robust dictionary-
completion queries.  This is the intended kind of cross-model validation for
the response framework: the theorem concerns the query algebra rather than
the vocabulary of either model.

## 7. What did and did not abstract

What survives is a three-part mechanism:

1. a source contains a low-cost carrier and an expensive quotient;
2. a future query can append another carrier in the same interface;
3. carrier join converts a geometric distance into a loss of expensive
   quotient directions before the extremum is evaluated.

What does **not** survive as a substantive theorem is the bare statement
“separated response maps require different summaries.”  That is generic
metric packing.  The reusable mathematical task in a new model is to prove
an analogue of (SN.8) from its composition algebra.  If the analogue must be
assumed, selective neutralization has explained nothing.

### Proposition SN.3 (one-carrier ceiling)

The complete class of dictionaries (SN.5), allowing every carrier dimension,
every carrier, every complement, and every complement basis, has at most

```math
q^{O(w^2)}                                       \tag{SN.15}
```

members.  Consequently no response-packing argument confined to this class,
under any future-query family and at any distortion, can force more than
`O(w^2 log q)` description bits.

#### Proof

There are at most `w+1` choices of dimension and at most
`C_q q^{d(w-d)}` carriers of dimension `d`.  After fixing a carrier, an
ordered list of at most `w` ambient vectors over `K` overcounts all possible
complement bases by at most `q^{w^2}`.  Summing over dimensions proves
(SN.15); the logarithm of the source cardinality bounds every packing.
`square`

Thus the quadratic lower bound is of the correct order for the whole
one-carrier architecture.  Exponential-in-`w` macroscopic response
complexity, if true for arbitrary supports, requires either many carriers
whose union does not collapse the quotient or a qualitatively different
nonlinear feature.  This is a genuine stopping statement for the present
construction, not a claim that arbitrary syndrome supports have only
quadratic response complexity.

## 8. Bounded multi-carriers and flags collapse

The obvious next attempt is to store several dense carriers or a flag.  For
the same query interface this does not evade the ceiling.

Let `m>=1`, let `B` be any background dictionary, and let
`W_1,...,W_m<=G`.  Assume the dictionary `D` below spans `G`, and put

```math
D=B\cup\bigcup_{i=1}^m P(W_i),
\qquad
L=W_1+\cdots+W_m,
\qquad
\bar D=B\cup P(L).                               \tag{SN.16}
```

### Theorem SN.4 (carrier-span compression)

For every future dictionary `E` and every target `x`,

```math
0\le
\ell_{D\cup E}(x)-\ell_{\bar D\cup E}(x)
\le m-1.                                         \tag{SN.17}
```

Consequently

```math
0\le R_D(E)-R_{\bar D}(E)\le m-1                \tag{SN.18}
```

uniformly over **all** future augmentations.  In particular, `m=o(w)` dense
carriers contain no more normalized future-response information than their
single span carrier, up to `o(w)` error.

#### Proof

Since `D subseteq bar D`, the left inequality in (SN.17) is immediate.  In a
shortest representation using `bar D union E`, combine all atoms drawn from
`P(L)` into their sum.  Thus at most one such atom is needed.  If it is
nonzero, write it as

```math
z=z_1+\cdots+z_m,
\qquad z_i\in W_i.
```

Delete zero terms and replace the one `P(L)` atom by at most `m` atoms, one
from each `P(W_i)`.  Every background and future atom is left unchanged, so
the representation length grows by at most `m-1`.  This proves (SN.17).
Taking maxima over `x` proves (SN.18). `square`

This is stronger than a source-counting observation: it gives a uniform
query-by-query approximation and permits an arbitrary, even adversarial,
background `B` and future context `E`.

### Corollary SN.5 (bounded-carrier and flag ceiling)

Fix `m`.  The class of states built from at most `m` dense linear carriers,
with complement bases and a bounded number of labels drawn from fixed or
`q^{O_m(w^2)}`-sized alphabets, has at most
`q^{O_m(w^2)}` members and hence at most `O_m(w^2 log q)` response bits under
every query family.  Moreover, by Theorem SN.4 its dense-carrier part is
uniformly within `m-1` of the state containing only their span.

A nested flag is even more degenerate: if

```math
W_1\subseteq\cdots\subseteq W_m,
```

then `union_i P(W_i)=P(W_m)` exactly.  Allowing a complete labeled flag does
not create superquadratic positional entropy either, because every flag is
the sequence of initial spans of an ordered basis and there are fewer than
`q^{w^2}` ordered bases.

#### Proof

For each carrier there are `w+1` dimension choices and at most
`C_q q^{w^2/4}` subspaces.  A bounded number of complement bases and labels
adds at most `q^{O_m(w^2)}` possibilities.  Take logarithms.  The exact flag
collapse and ordered-basis count prove the remaining claims. `square`

### Resource required beyond this architecture

Neither a flag nor any bounded (indeed, any `o(w)`) collection of dense
linear carriers can produce a new macroscopic feature: Theorem SN.4 replaces
it by one span at sublinear error.  A packing with exponentially many bits in
`w` would therefore need at least one of the following genuinely new
resources:

1. a number of independently positioned carriers growing so rapidly that
   their positional source entropy itself is exponential in `w`, together
   with a query that selects a macroscopic subset **without** their union
   collapsing to its span;
2. carriers that are nonlinear and cannot synthesize their joint span with
   one atom per carrier; or
3. a response whose future composition retains rooted multiplicity or
   incidence information erased by projective support union.

Item 1 cannot be achieved merely by listing `m` subspaces unless
`m` is itself exponential up to polynomial factors, since an `m`-tuple has
only `O(mw^2 log q)` positional bits.  This identifies the exact obstruction
to iterating the Grassmannian construction rather than leaving “try more
carriers” as an open-ended variant.
