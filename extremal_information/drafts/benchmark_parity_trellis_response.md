# Solution-hidden benchmark: parity trellis response

**Status.**  The operational quotient, composition law, raw-interface caveat,
and approximate response-rate bound were frozen on 2026-08-17 before the
literature comparison in Section 8.  The finite verifier is
[`verify_parity_trellis_response.py`](../experiments/verify_parity_trellis_response.py).

## 1. Declared continuation experiment

All vector spaces are over `F_2`.  First consider a width-`w` boundary
space `X=F_2^w`.  Let the boundary assignments reachable from a fragment be
the affine space

```math
R=x_0+U\subseteq X.                                           \tag{PT.1}
```

Let `L<=X^*` be the span of every parity form that a declared future may
read, constrain, or price.  The visible port of `x` is the evaluation map

```math
\pi_L(x)\in L^*,\qquad \pi_L(x)(\ell)=\ell(x).                \tag{PT.2}
```

The future class is assumed to be rich enough to pin every realizable value
of this port.  This is the exact exposure assumption used below.  It is
satisfied, for example, when affine parity constraints on a basis of `L`
are allowed.

Two reachable assignments are **operationally equivalent** if every
declared future has the same feasible completions and the same conditional
optimum from them.  Absolute optimum costs are compared literally; an
additive prefix cost is not silently quotiented out.

### Proposition PT.1 (visible-row-space quotient)

Under the exposure assumption,

```math
x\sim_L y
\quad\Longleftrightarrow\quad
\pi_L(x)=\pi_L(y)
\quad\Longleftrightarrow\quad
x+y\in L^\perp.                                               \tag{PT.3}
```

Consequently the coarsest exact boundary state is the affine quotient

```math
R/(U\cap L^\perp),                                            \tag{PT.4}
```

and its number of states is

```math
q=2^d,\qquad
d=\dim U-\dim(U\cap L^\perp).                                \tag{PT.5}
```

#### Proof

Every allowed future depends on `x` only through `pi_L(x)`, so equal ports
are sufficient.  If the ports differ, the future that pins the port to
`pi_L(x)` accepts `x` and rejects `y`; hence they are not equivalent.  The
kernel of `pi_L|_U` is `U intersect L^perp`, and rank-nullity gives (PT.5).
`square`

This proposition already records the main no-go condition.  If every raw
boundary bit is exposed, then `L=X^*` and `L^perp=0`: distinct reachable raw
assignments cannot be merged.  Syndrome compression is valid only when the
semantic interface is a parity port rather than the full raw assignment.

## 2. Linear-code cut and compatible-future cosets

Let

```math
C\leq F_2^P\oplus F_2^F                                      \tag{PT.6}
```

be a binary linear code split into past and future coordinates.  Define the
past- and future-supported shortened subcodes

```math
C_P=\{p:(p,0)\in C\},\qquad
C_F=\{f:(0,f)\in C\}.                                        \tag{PT.7}
```

For a reachable past word `p in proj_P(C)`, its compatible-future set is

```math
\mathcal F(p)=\{f:(p,f)\in C\}.                              \tag{PT.8}
```

The declared futures may impose arbitrary affine parity constraints and
costs on future variables.  In particular, they may hard-fix a future word.

### Theorem PT.2 (coarsest exact trellis state)

For reachable `p,p'`,

```math
\mathcal F(p)=\mathcal F(p')
\quad\Longleftrightarrow\quad p+p'\in C_P.                   \tag{PT.9}
```

Thus the coarsest exact state of the code cut is

```math
\boxed{\operatorname{proj}_P(C)/C_P
       \ \cong\ C/(C_P\oplus C_F)}.                          \tag{PT.10}
```

Its dimension and size are

```math
r=\dim C-\dim C_P-\dim C_F,
\qquad q=2^r.                                                 \tag{PT.11}
```

#### Proof

If the two future sets agree, choose `f` in their common nonempty set.
Adding `(p,f)` and `(p',f)` gives `(p+p',0) in C`, so `p+p' in C_P`.
Conversely, adding `(p+p',0) in C` to any `(p,f) in C` gives `(p',f) in C`,
so the future sets agree.

Each nonempty `mathcal F(p)` is an affine coset of `C_F`; two such cosets
are equal or disjoint.  If two states in (PT.10) differ, choose a word `f`
in one compatible-future coset and hard-fix the future to `f`.  That context
is feasible from one state and infeasible from the other, proving exact
minimality.  Finally, the kernel of the natural map from `C` to the state
space is `C_P direct-sum C_F`, giving (PT.10)--(PT.11). `square`

The theorem distinguishes the **state label** from its dynamic-programming
metric.  Prefix paths with the same state label may be merged for a minimum
query by retaining the least accumulated cost.  Paths in different labels
cannot be merged under the declared future class.

## 3. Partial-syndrome realization and strict compression

Suppose `C=ker H` with a parity-check matrix split as

```math
H=[H_P\ H_F].                                                  \tag{PT.12}
```

For a compatible pair `(p,f)`,

```math
H_Pp=H_Ff.                                                    \tag{PT.13}
```

The partial-syndrome map

```math
p+C_P\longmapsto \sigma=H_Pp                                \tag{PT.14}
```

is therefore an isomorphism onto

```math
S=\operatorname{im}H_P\cap\operatorname{im}H_F.             \tag{PT.15}
```

Indeed, its kernel is precisely `C_P`, and (PT.13) proves that its image is
the intersection.  Hence

```math
\begin{aligned}
r
 &=\dim(\operatorname{im}H_P\cap\operatorname{im}H_F)\\
 &=\operatorname{rank}H_P+\operatorname{rank}H_F
   -\operatorname{rank}[H_P\ H_F]\\
 &=\dim C-\dim C_P-\dim C_F.                                \tag{PT.16}
\end{aligned}
```

Only reachable **and completable** syndromes survive.  Thus `q=2^r` may be
strictly smaller than the `2^w` assignments of a declared `w`-bit check
port.  For example, in `F_2^4`, take

```math
\operatorname{im}H_P=\langle e_1,e_2,e_3\rangle,
\qquad
\operatorname{im}H_F=\langle e_1,e_2,e_4\rangle.             \tag{PT.17}
```

Their intersection has dimension two, giving four exact states instead of
the naive sixteen check-bit patterns.

## 4. Exact weighted response and composition

For a min-cost prefix with possibly many paths, define one value per exact
state:

```math
\alpha(s)=\min\{\text{prefix cost of a realization ending in }s\}. \tag{PT.18}
```

Infeasibility has value `+infinity`.  A future context `D` induces

```math
\beta_D(s)=\min\{\text{future cost completing }s\},           \tag{PT.19}
```

and the closed response is

```math
Q_\alpha(D)=\min_s\{\alpha(s)+\beta_D(s)\}.                  \tag{PT.20}
```

### Theorem PT.3 (coarsest aggregated response)

Under hard future-word pins, two finite prefix tables have the same response
to every future if and only if they agree coordinatewise.

#### Proof

Equality is sufficient by (PT.20).  For each state `s`, choose a future word
whose syndrome is `s` and hard-fix it.  Its future table is zero at `s` and
`+infinity` elsewhere, so the response is exactly `alpha(s)`.  All
coordinates are therefore operationally exposed. `square`

A block between left states `S` and right states `T` has transfer kernel

```math
K(s,t)=\min\{\text{block cost conditional on }s,t\}.          \tag{PT.21}
```

Serial composition is min-plus matrix multiplication:

```math
(K\star L)(s,u)=\min_t\{K(s,t)+L(t,u)\},
\qquad
\alpha'(t)=\min_s\{\alpha(s)+K(s,t)\}.                      \tag{PT.22}
```

For a binary block whose assignment `b` contributes syndrome `delta(b)`, a
transition obeys `delta(b)=s+t`.  The max-completion version replaces
`min,+infinity` by `max,-infinity`; all quotient and minimality statements
are unchanged.

## 5. Approximate response rate

Assume `q` states are probeable and the weighted fragment family can realize
every table `alpha in [0,B]^q`.  Define

```math
d_{\rm resp}(\alpha,\alpha')
=\sup_D|Q_\alpha(D)-Q_{\alpha'}(D)|.                          \tag{PT.23}
```

The supremum is over contexts with at least one feasible state, so both
closed responses are finite.

### Theorem PT.4 (future-response isometry and coding rate)

Under the declared hard probes,

```math
\boxed{d_{\rm resp}(\alpha,\alpha')
       =\|\alpha-\alpha'\|_\infty}.                          \tag{PT.24}
```

If a summary must answer every future query within additive `epsilon`, then,
for `0<epsilon<B/6`, its worst-case bit length `b_epsilon` satisfies

```math
q\log_2\!\left(\left\lfloor{B\over3\epsilon}\right\rfloor+1\right)
\le b_\epsilon
\le q\left\lceil\log_2\!\left(\left\lceil{B\over\epsilon}\right\rceil+1\right)\right\rceil.
                                                                    \tag{PT.25}
```

In particular,

```math
b_\epsilon=\Theta\!\left(q\log(B/\epsilon)\right)
=\Theta\!\left(2^r\log(B/\epsilon)\right).                  \tag{PT.26}
```

#### Proof

The elementary minimum inequality gives the upper bound in (PT.24), while a
coordinate hard probe attains each coordinate difference.  For the rate
lower bound, use the Cartesian grid of mesh `3 epsilon` in `[0,B]^q`.  Any
two grid tables have response distance at least `3 epsilon`; if they shared
one summary, a common decoded answer accurate to `epsilon` for both would
force their response distance to be at most `2 epsilon`.  Thus all grid
tables need distinct summaries.  Quantizing each coordinate to an
`epsilon`-mesh gives the upper bound. `square`

Negating all costs gives the identical isometry and rate for max-completion
queries.

The table family assumption in Theorem PT.4 is material.  Arbitrary bounded
boundary-state penalties realize it in a weighted parity-CSP.  A strict
Hamming-cost code family may occupy a lower-dimensional subset of the cube;
its sharp lower bound is the packing entropy of that actual response family,
not automatically the right side of (PT.25).

There is also a useful restricted-future pruning statement.  Suppose every
future table is finite and has oscillation at most `L`.  Put
`m=min_s alpha(s)`.  States with `alpha(s)>m+L` never win.  More generally,
for `0<=epsilon<=L`, retaining only states satisfying

```math
\alpha(s)\le m+L-\epsilon                                    \tag{PT.27}
```

changes every response by at most `epsilon`.  If `s_0` minimizes `alpha`
and a removed `s` were optimal, then

```math
\alpha(s_0)+\beta(s_0)
\le m+\beta(s)+L
<\alpha(s)+\beta(s)+\epsilon.                                \tag{PT.28}
```

Hard pins have infinite oscillation, so the full arbitrary-future class
intentionally rules out this pruning.

## 6. Caveats and scope

1. If a future may fix each raw separator variable, unequal raw assignments
   are distinguishable.  A syndrome quotient is then not an operational
   quotient unless additional structure first proves that the response
   factors through that syndrome.
2. The `2^r` count concerns concrete trellis labels.  A prefix containing
   alternatives generally needs the entire `q`-entry metric table
   (PT.18), not one syndrome and one number.
3. Absolute optimum queries retain a common additive offset.  Removing it
   changes the experiment to responses modulo constants; it is a gauge
   choice, not exact equivalence for literal costs.
4. The proof uses linear/affine parity constraints.  Nonlinear future
   predicates may refine the quotient.
5. The finite verifier is a falsification audit of tiny instances.  The
   theorem follows from the algebraic proof, not from finite enumeration.

## 7. Finite verification

The dependency-free verifier represents a parity-check column by an integer
bit mask and exhausts every past and future word.  For each code it checks

- `F(p)=F(p')` iff `p+p' in C_P`;
- equality of this partition with partial-syndrome equality;
- hard-future-word exposure of every state;
- `|proj_P(C)/C_P|=|C/(C_P direct-sum C_F)|`;
- the intersection image and all three dimensions in (PT.16).

It audits four named examples, including strict compression, zero
intersection, and nontrivial shortened subcodes.  It also checks every
parity-check column tuple with one or two past columns, one or two future
columns, and ambient check dimension one or two.

Run from the repository root:

```bash
.venv/bin/python extremal_information/experiments/verify_parity_trellis_response.py
```

## 8. Post-freeze literature comparison

The independently derived quotient specializes exactly to classical
minimal-trellis state spaces.

- Wolf constructs a syndrome trellis with at most `q^(n-k)` states and
  observes that code structure can give fewer:
  [*Efficient maximum likelihood decoding of linear block codes using a
  trellis*](https://doi.org/10.1109/TIT.1978.1055821) (1978).
- Forney and Trott's State Space Theorem identifies the cut state as the
  quotient by past- and future-supported subcodes, matching
  `C/(C_P direct-sum C_F)`:
  [*The dynamics of group codes: state spaces, trellis diagrams, and
  canonical encoders*](https://doi.org/10.1109/18.259635) (1993).
- Forney relates these quotient dimensions to minimal block-code trellis
  complexity:
  [*Dimension/length profiles and trellis complexity of linear block
  codes*](https://doi.org/10.1109/18.340452) (1994).
- Min-plus transfer composition is the optimization-semiring instance of
  generalized distributive-law message passing:
  [Aji--McEliece, *The generalized distributive law*](https://doi.org/10.1109/18.825794)
  (2000), with the code-graph realization in
  [Forney, *Codes on graphs: normal realizations*](https://doi.org/10.1109/18.910573)
  (2001).

The visible-row-space formulation, raw-interface no-go distinction,
future-probe response isometry, rate bound, and bounded-oscillation pruning
were obtained before this comparison.  The checked primary sources confirm
the exact trellis quotient and composition.  No novelty claim is made for
the operational reformulation or the elementary metric-entropy consequence.

## Benchmark verdict

**Pass, independently predicted.**  Future-context equivalence recovers the
classical minimal syndrome/coset state, explains exactly when it is smaller
than all width-`w` assignments, and extends directly to an exact response
metric and a matching approximate storage-rate law under a declared rich
weighted fragment family.
