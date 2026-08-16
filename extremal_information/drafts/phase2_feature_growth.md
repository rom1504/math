# Phase 2 feature-algebra growth: extremal transfer kernels

## Verdict

There is a second exact closure class beyond fixed-rank mean-field models:
finite-width landscapes whose future environment meets a component only
through a finite separator.  The correct reusable extremal state is the
component's **boundary response kernel**.  For a two-terminal component this
kernel is a max-plus matrix.

This is related to transfer matrices and dynamic programming, but the result
below is stronger than the algorithmic observation that a transfer table is
sufficient:

1. the kernel is the coarsest exact statistic for all endpoint contexts;
2. the endpoint-response metric is *isometric* to entrywise sup distance of
   kernels;
3. serial gluing is max-plus matrix multiplication and is nonexpansive;
4. arbitrary binary kernels give a sharp `Q^2`-bit lower bound, so the
   separator dependence cannot be improved for the universal class; and
5. the upper response roof is precisely an optimal-transport extension of
   the kernel.

Thus this model supplies both a genuine sub-landscape closure theorem and a
scalable lower bound on feature-algebra growth.

The result does **not** say that arbitrary bounded-width ground-state dynamic
programming is new.  Its contribution to the present theory is the exact
operational minimality, distortion geometry, and compositional information
lower bound.

## 1. The endpoint-query experiment

Let `B` be a finite boundary alphabet of cardinality `Q`.  A finite
two-terminal landscape is

```math
\mathcal C=(\Omega,H,\ell,r),
```

where `H:Omega -> R` and `ell,r:Omega -> B`.  Assume for now that every
endpoint fibre

```math
\Omega_{ab}=\{x:\ell(x)=a,\ r(x)=b\}
```

is nonempty.  Define the boundary response kernel

```math
K_{\mathcal C}(a,b)=\max_{x\in\Omega_{ab}}H(x).       \tag{FG.1}
```

The declared query apparatus adds arbitrary left and right endpoint rewards
`f,g:B -> R` and observes

```math
V_{\mathcal C}(f,g)
=\max_{x\in\Omega}
 \{H(x)+f(\ell(x))+g(r(x))\}.                         \tag{FG.2}
```

This is deliberately a restricted experiment: it does not query internal
occupancies, multiplicities, or overlaps.

For two kernels on the same boundary alphabet put

```math
K\odot L(a,c)=\max_{b\in B}\{K(a,b)+L(b,c)\}.         \tag{FG.3}
```

This is max-plus matrix multiplication.

## 2. Exact quotient, isometry, and roof identification

### Theorem FG.1 (the extremal transfer-kernel theorem)

For finite two-terminal landscapes with finite kernels, the following hold.

**(a) Exact response quotient.**

```math
V_{\mathcal C}(f,g)
=\max_{a,b\in B}\{K_{\mathcal C}(a,b)+f(a)+g(b)\}.
                                                               \tag{FG.4}
```

Consequently `K_C` determines every declared response.

**(b) Operational minimality and isometry.**  Define

```math
d_{\rm end}(\mathcal C,\mathcal D)
=\sup_{f,g\in\mathbb R^B}
 |V_{\mathcal C}(f,g)-V_{\mathcal D}(f,g)|.
```

Then

```math
d_{\rm end}(\mathcal C,\mathcal D)
=\|K_{\mathcal C}-K_{\mathcal D}\|_\infty.          \tag{FG.5}
```

In particular, any deterministic summary that answers all endpoint queries
exactly determines the entire kernel.  Equality of kernels is the coarsest
exact response equivalence, up to a one-to-one recoding.

**(c) Identification with the upper response roof.**  Give a state the
feature

```math
\phi(x)=(e_{\ell(x)},e_{r(x)})\in\mathbb R^Q\times\mathbb R^Q.
```

For `p,q` in the probability simplex `Delta_B`, its upper response roof is

```math
\widehat H_\phi(p,q)
=\max_{\substack{\pi\in\Delta(B\times B)\\
                  \pi_1=p,\ \pi_2=q}}
  \sum_{a,b}\pi(a,b)K_{\mathcal C}(a,b).             \tag{FG.6}
```

Thus the roof is the optimal-transport concave extension of the boundary
kernel.  Conversely,

```math
\widehat H_\phi(e_a,e_b)=K_{\mathcal C}(a,b),        \tag{FG.7}
```

so the roof and kernel contain exactly the same information.

#### Proof

For fixed endpoints `(a,b)`, an endpoint reward is constant on the fibre
`Omega_ab`; maximizing there gives (FG.4).

The standard maximum inequality applied to (FG.4) gives

```math
d_{\rm end}(\mathcal C,\mathcal D)
\le \|K_{\mathcal C}-K_{\mathcal D}\|_\infty.
```

For the reverse inequality, fix `(a_0,b_0)`.  Put `f(a_0)=g(b_0)=0` and put
`f(a)=g(b)=-T` off the selected coordinates.  For `T` larger than the
oscillations of both finite kernels, `(a_0,b_0)` is the maximizing endpoint
pair for both landscapes.  Their response difference is then exactly

```math
|K_{\mathcal C}(a_0,b_0)-K_{\mathcal D}(a_0,b_0)|.
```

Choose an entry attaining the sup norm.  This proves (FG.5) and minimality.

For (FG.6), any mixture of landscape states induces a joint endpoint law
`pi`, and its conditional expected energy on the fibre `(a,b)` is at most
`K_C(a,b)`.  Conversely, for each `(a,b)` choose a fibre maximizer and mix
these maximizers according to any coupling `pi`.  The prescribed feature
means are exactly the two marginals of `pi`, proving equality.  At the vertex
`(e_a,e_b)`, the only possible joint endpoint law is the point mass at
`(a,b)`, giving (FG.7). `square`

### Why the roof geometry matters

The visible feature vector has ambient dimension only `2Q`, but its roof can
carry `Q^2` independently variable vertex heights.  Feature dimension alone
therefore does not measure extremal information.  The metric entropy of the
exposed roof does.  The lower bound in Section 4 makes this distinction
sharp.

Equation (FG.6) also identifies the interior roof: it is not an arbitrary
concave interpolation, but an optimal-transport value with cost `K`.  This is
one concrete case in which a standard model's transfer state is exactly a
response roof rather than merely analogous to one.

## 3. Closure under repeated composition

Serially glue `C=(Omega,H,ell,r)` and
`D=(Xi,G,ell',r')` by identifying their adjacent endpoint states.  The parent
state space is the fibre product

```math
\Omega\mathbin{\times_B}\Xi
=\{(x,y):r(x)=\ell'(y)\},
```

with energy `H(x)+G(y)` and exposed endpoints `(ell(x),r'(y))`.

### Theorem FG.2 (exact max-plus closure and data processing)

The serial composite satisfies

```math
K_{\mathcal C\circ\mathcal D}
=K_{\mathcal C}\odot K_{\mathcal D}.                \tag{FG.8}
```

The product is associative.  Moreover,

```math
\|K\odot L-K'\odot L'\|_\infty
\le \|K-K'\|_\infty+\|L-L'\|_\infty.              \tag{FG.9}
```

Consequently, if `t` component kernels have respective sup errors
`epsilon_1,...,epsilon_t`, every endpoint response of their serial composite
has error at most

```math
\epsilon_1+\cdots+\epsilon_t.                       \tag{FG.10}
```

#### Proof

Condition on the shared endpoint `b` and maximize independently in the two
fibres.  This gives (FG.8).  Associativity follows either from associativity
of gluing or by expanding both bracketings as a maximum over all intermediate
boundary states.

For every `(a,c,b)`, the two summands in the primed and unprimed products
differ in total by at most the right-hand side of (FG.9).  Taking a maximum
cannot enlarge a uniform error.  Now use Theorem FG.1(b) and induction to get
(FG.10). `square`

### Corollary FG.3 (bounded-separator closure)

Let a factor landscape be split into an interior component and an environment
that share only a separator assignment `s in S`, with additive energy across
the separator.  Every local factor is allocated to exactly one component, so
no separator factor is double counted.  Then the conditional maximum profile

```math
m_{\mathcal C}(s)=\max\{H(x):\text{boundary}(x)=s\}  \tag{FG.11}
```

is the coarsest exact summary for arbitrary separator potentials.  Gluing
uses addition of conditional profiles followed by maximization over eliminated
separator variables.  An infeasible boundary assignment is represented by
`-infinity` (or by an explicit feasibility bit), and all maxima are restricted
to feasible assignments.

For a `q`-ary separator of width `w`, this profile has `q^w` values.  A
reusable strip segment with a left and a right width-`w` boundary has a
`q^w` by `q^w` kernel, hence `q^{2w}` values.  These sizes are independent of
the number of internal variables.

If local energies are integers of magnitude at most `W` and there are `N`
local factors, exact storage takes

```math
O(q^w\log(1+NW))
```

bits for a one-sided separator profile and

```math
O(q^{2w}\log(1+NW))
```

bits for a reusable two-sided strip kernel, plus one feasibility bit per entry
when infeasible fibres occur.  The full landscape table can have exponentially
many state entries even though, under the displayed integer-range assumption,
it has only `O(NW)` distinct numerical energy values.

The proof is Theorem FG.1 with `B=S` in the one-sided case and repeated
application of Theorem FG.2 in a tree or path decomposition, allocating every
local factor exactly once. `square`

This is a deterministic finite synchronization statement: conditional on
the separator assignment, the internal optimizer and every external context
decouple exactly.  No overlap identities or randomness are required.  It is
also a precise stopping criterion: if no small separator screens the future
query, this theorem supplies no compression.

## 4. Sharp separator information lower bound

For the universal class of arbitrary boundary kernels, the exponential
dependence on separator width is not an artifact of the proof.  A restricted
local constraint language may generate a much smaller kernel semigroup.

### Theorem FG.4 (sharp metric entropy of universal transfer kernels)

For each binary matrix `A in {0,1}^{B times B}`, define the component

```math
\Omega_A=B\times B,
\qquad H_A(a,b)=A_{ab},
\qquad \ell(a,b)=a,
\qquad r(a,b)=b.
```

Then its kernel is `A`.  Distinct matrices have endpoint-response distance
exactly one.  Therefore any deterministic summary that answers all endpoint
queries on this class with uniform additive error `epsilon<1/2` needs at
least

```math
\log_2 2^{Q^2}=Q^2                              \tag{FG.12}
```

bits.  Storing the binary kernel attains this bound exactly.

More generally, for kernels in `[-B_0,B_0]^(Q times Q)`, the response metric
is entrywise sup distance, so in the nontrivial regime
`0<epsilon<=B_0/2` its logarithmic covering
complexity is

```math
\Theta\!\left(Q^2\log {B_0\over\epsilon}\right),   \tag{FG.13}
```

with universal constant-factor changes from the choice of packing or
covering convention.

The binary packing persists for arbitrarily long nearest-neighbor chains.
Fix `C>1` and let

```math
I_C(a,b)=
\begin{cases}
0,&a=b,\\
-C,&a\ne b.
\end{cases}
```

For every binary kernel `A`,

```math
A\odot I_C=A.                                       \tag{FG.14}
```

Hence the length-`L` chain with first transition kernel `A` and all remaining
transition kernels `I_C` has endpoint kernel `A` for every `L`.  Its full
path space has `Q^(L+1)` states, while the `2^(Q^2)`-element hard packing and
the `Q^2`-bit lower bound remain unchanged.

#### Proof

The first assertion follows directly from (FG.1).  Theorem FG.1(b) says the
response distance is sup distance, which is one for two distinct binary
matrices.  Two such matrices cannot share a message under error less than
one half, since their decoded response functions would be at distance less
than one by the triangle inequality.  Hence all `2^(Q^2)` matrices require
different messages.  The continuous estimate is the elementary grid
packing/covering estimate for a `Q^2`-dimensional cube in sup norm. `square`

For the persistence claim, the term with intermediate state `b=c` in
`(A odot I_C)(a,c)` equals `A_ac`.  Every term with `b ne c` is at most
`1-C<0`, while `A_ac>=0`; therefore the maximum is `A_ac`. `square`

For `B=[q]^w`, Theorem FG.4 gives a `q^{2w}`-bit lower bound for universal
binary two-sided kernels.  For a one-sided separator, the identical argument
with binary vectors gives `q^w` bits.  Thus the familiar exponential-in-width
barrier is an extremal-information lower bound, not merely a limitation of a
particular elimination algorithm.

The lower bound uses a universal factor on the exposed boundary variables.
Restricted local-factor families may have smaller kernel classes; in that
case the correct question is the metric entropy of the actually realizable
kernel semigroup, not the ambient cube.

## 5. External applications

### 5.1 Finite-state spin chains

Let `B` have `Q` spin states and consider an inhomogeneous nearest-neighbor
chain

```math
H(s_0,\ldots,s_L)=\sum_{i=1}^L W_i(s_{i-1},s_i).
```

The full landscape has `Q^(L+1)` states.  A one-edge component has kernel
`W_i`, and an interval has the exact reusable state

```math
K_{[1,L]}=W_1\odot W_2\odot\cdots\odot W_L.         \tag{FG.15}
```

It contains only `Q^2` values, independent of `L`, answers every pair of
endpoint-field queries, and composes without bracket dependence.  For an
Ising chain `Q=2`, four values replace the exponential energy table.  The
same construction applies to finite-range chains after taking a boundary
state to be a word of the required memory length.

This is genuinely outside fixed-rank Curie--Weiss: interactions are sparse,
local, and may be completely inhomogeneous.  Compression comes from an exact
separator/Markov property, not a global low-rank order parameter.

### 5.2 Trellis decoding and weighted automata

A path through a finite-state trellis is a landscape state; its score is the
sum of branch metrics.  A segment's best score for every start/end state is
exactly `K`.  Concatenating trellis segments is (FG.8).  Therefore the
standard max-product/Viterbi transfer matrix is the minimal exact statistic
for a segment that must remain reusable under arbitrary initial and terminal
metrics.

Classical sources establish the algorithmic max-product/distributive-law
background:

- G. D. Forney, “The Viterbi Algorithm,” *Proceedings of the IEEE* 61
  (1973), 268--278, DOI
  [10.1109/PROC.1973.9030](https://doi.org/10.1109/PROC.1973.9030).
- S. M. Aji and R. J. McEliece, “The Generalized Distributive Law,”
  *IEEE Transactions on Information Theory* 46 (2000), 325--343,
  DOI [10.1109/18.825794](https://doi.org/10.1109/18.825794), with an
  [author-hosted copy](https://authors.library.caltech.edu/records/sw1pm-bwj40).
- F. R. Kschischang, B. J. Frey, and H.-A. Loeliger, “Factor Graphs and the
  Sum-Product Algorithm,” *IEEE Transactions on Information Theory* 47
  (2001), 498--519, DOI
  [10.1109/18.910572](https://doi.org/10.1109/18.910572).

Those works provide the semiring/message-passing architecture.  The exact
response isometry, operational minimality for all endpoint contexts, and
the sharp query-information lower bound above are the additions relevant to
this program.

## 6. What the state remembers and forgets

The kernel remembers:

- every conditional ground-state energy at the declared separator;
- every response to arbitrary endpoint rewards;
- exactly the information needed for any further serial gluing through that
  separator.

It forgets:

- the identity and multiplicity of internal optimizers;
- second-best and near-optimal energies;
- internal overlaps and local observables;
- finite-temperature weight; and
- responses to a future edge that bypasses the separator and touches an
  internal variable.

The forgetting can be enormous.  Each endpoint fibre may contain
exponentially many states with arbitrary energies below its fibre maximum;
all such modifications leave `K` and every declared future response
unchanged.

## 7. Falsifiers and scope boundaries

The theorem should be rejected outside its declared experiment in any of the
following situations.

1. **An internal query is permitted.**  With a singleton boundary, one
   component may have states of energies `0,-1` and another only the energy
   `0` state.  They have the same kernel, but a field coupled to the second
   internal state separates them.
2. **The environment bypasses the separator.**  A cross interaction depending
   on an internal spin is not a function of the endpoint kernel.
3. **Counts or positive temperature matter.**  Maxima discard fibre
   multiplicities.  A log-sum-exp transfer kernel, not (FG.1), is then needed.
4. **Endpoint fields are bounded too weakly to pin endpoint pairs.**  The
   exact isometry (FG.5) uses arbitrary fields.  For a bounded query set the
   operational quotient may be strictly coarser.
5. **Width grows with volume.**  If `w=Theta(N)`, the `q^(2w)` state is again
   exponential and this closure is not a useful compression.
6. **The local-factor class is restricted.**  The universal lower bound may
   overstate the complexity of its realizable transfer semigroup.  One must
   prove that the hard kernel packing lies in the model class.

These are falsifiable structural boundaries, not technical caveats.

## 8. Director checkpoint

### Is this more than repackaged convex duality or dynamic programming?

Partly classical, but yes at the theorem level.

- Sufficiency and max-plus composition are the familiar transfer-matrix
  mechanism.
- The response-roof identification (FG.6) says exactly which convex object
  the transfer matrix represents: an optimal-transport roof over endpoint
  marginals.
- The isometry (FG.5) turns endpoint experiments into entrywise kernel
  geometry.
- Operational minimality and Theorem FG.4 prove that exponential separator
  dependence is unavoidable for the universal composable experiment.
- Theorems FG.1--FG.2 distinguish the state needed to answer one ground-state
  query from the state needed to remain reusable under every future context.

This produces something the framework could not previously formulate: a
sharp law saying that **extremal feature-algebra complexity is controlled by
the information capacity of the interaction separator**, with exact closure
at fixed separator width and a matching exponential lower bound when width
grows.

### Strongest next theorem suggested by this result

Let `K_{n,w}` be the semigroup of boundary kernels actually realizable by a
specified structured model of width `w` (linear codes on trellises, a fixed
finite CSP language, or bounded-rank transition factors).  Determine its
response-metric entropy:

```math
\log\operatorname{Cov}(K_{n,w},\|\cdot\|_\infty,\epsilon).
```

The universal ambient value is of order `q^(2w)`, but algebraic restrictions
may force a dramatically smaller tropical semigroup.  Proving such a bound,
together with closure under `odot`, would be a nontrivial structured
feature-algebra theorem rather than generic separator dynamic programming.
