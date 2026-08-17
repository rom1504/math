# Bridge-query isometry and a sparse incompressibility theorem

Status: main-agent proof draft for independent audit.

This note derives the exact semantic interface of an arbitrary bipartite
coupling and uses it to show that bounded bridge degree alone gives no
extremal compression, even at macroscopic additive error.

## 1. The restricted bridge transform

Let `X,Y` be finite, let `B:X times Y -> R`, and let `h:X->R` be an internal
landscape.  Define its bridge-response table

```math
(P_Bh)(y)=\max_{x\in X}\{h(x)+B(x,y)\}.                       \tag{BQ.1}
```

A future continuation is an arbitrary `g:Y->R`, and gluing returns

```math
Opt_B(h,g)=\max_{x,y}\{h(x)+B(x,y)+g(y)\}
          =\max_y\{(P_Bh)(y)+g(y)\}.                          \tag{BQ.2}
```

### Theorem BQ.1 (bridge-query isometry)

For arbitrary real future continuations,

```math
\sup_g|Opt_B(h,g)-Opt_B(h',g)|
=\|P_Bh-P_Bh'\|_\infty.                                      \tag{BQ.3}
```

If landscapes are identified modulo additive constants, the quotient metric
is

```math
\inf_c\|P_Bh-P_Bh'-c\mathbf1\|_\infty
={1\over2}\operatorname{osc}(P_Bh-P_Bh').                    \tag{BQ.4}
```

#### Proof

The maximum map is one-Lipschitz in sup norm, proving `<=` in (BQ.3).  For
the reverse inequality choose a coordinate `y_0` at which the largest signed
difference is attained and choose `g(y_0)=0`, with every other `g(y)` a
sufficiently large negative number.  Both maxima in (BQ.2) are then forced to
`y_0`.  Apply the same argument after swapping `h,h'` if the largest absolute
difference has the other sign.  Formula (BQ.4) is the standard best constant
approximation of a finite real vector. `square`

This theorem makes three familiar interfaces instances of one law:

- if `B(x,y)=<phi(x),psi(y)>` with `phi` in `R^r`, `P_Bh` is a restriction of
  the upper-roof transform to the bridge query set `psi(Y)`;
- if `B` meets the past only through `w` boundary spins, `P_Bh` is the usual
  `2^w` boundary response table;
- if a symmetry of `h` and `B` has few query orbits, the response table
  descends to those orbits.

The state is not guessed from rank, sparsity, or symmetry.  It is exactly the
realizable image of (BQ.1), and its response rate--distortion is precisely the
metric entropy of that image in sup norm.

## 2. A degree-one bridge with extensive response information

### Theorem BQ.2 (macroscopic sparse-bridge packing)

Fix `0<delta<1/2`.  There are sets

```math
C_n\subseteq\{-1,1\}^n,
\qquad |C_n|\ge 2^{(1-h_2(\delta)-o(1))n},                    \tag{BQ.5}
```

whose pairwise Hamming distances are at least `delta n`, and a family of
landscapes `{h_sigma:sigma in {0,1}^{C_n}}` such that, for the matching
bridge

```math
B(x,y)=\langle x,y\rangle,                                    \tag{BQ.6}
```

their response tables have pairwise distance `delta n`:

```math
\|P_Bh_\sigma-P_Bh_\tau\|_\infty=\delta n
\quad(\sigma\ne\tau).                                        \tag{BQ.7}
```

Consequently, for every `epsilon<delta/2`, an **absolute-score** `epsilon n`-
accurate summary for this class under arbitrary appended landscapes requires
at least

```math
2^{|C_n|}
```

summary states, or at least

```math
|C_n|\ge2^{(1-h_2(\delta)-o(1))n}
```

bits.  The bipartite interaction graph in (BQ.6) is a matching: it has degree
one and treewidth one as an isolated graph.

#### Proof

A greedy Hamming packing gives (BQ.5).  Take `X=Y=C_n`, put

```math
h_\sigma(c)=\delta n\,\sigma(c),                               \tag{BQ.8}
```

and fix a query `y=c`.  The intended state has value

```math
h_\sigma(c)+\langle c,c\rangle=n+\delta n\,\sigma(c).
```

Every `c'!=c` has

```math
\langle c',c\rangle=n-2d_H(c',c)\le n-2\delta n,
```

so even after receiving the largest possible bonus in (BQ.8), it remains at
least `delta n` below `c` when `sigma(c)=0`, and farther below when
`sigma(c)=1`.  Thus `c` is the unique optimizer and

```math
(P_Bh_\sigma)(c)=n+\delta n\,\sigma(c).                        \tag{BQ.9}
```

Distinct binary labels differ at some `c`, proving (BQ.7).  Balls of radius
less than `delta n/2` around the `2^|C_n|` tables are disjoint, and Theorem
BQ.1 transfers this packing to future-response distortion. `square`

If a full-cube state space is desired, extend each `h_sigma` by a sufficiently
negative common value outside `C_n`; the proof is unchanged.

The same exponential-bit lower bound holds projectively.  Restrict the labels
to a constant-weight family, say `|sigma^(-1)(1)|=floor(|C_n|/2)`.  For two
distinct such labels, their response difference has both a `+delta n` and a
`-delta n` coordinate, so (BQ.4) gives projective distance `delta n`.  The
family has `binom(|C_n|,floor(|C_n|/2))=2^(|C_n|-o(|C_n|))`
members.

## 3. What the theorem rules out

The obstruction does not say that every sparse graphical model is hard.
Ordinary sparse dynamic programming relies on a decomposition whose *live
separator* is small.  The matching in Theorem BQ.2 exposes `Theta(n)` past
spins simultaneously, so its live interface is extensive despite degree one.

It does establish three boundaries.

1. Bounded degree, bounded local arity, and small treewidth of the bridge
   graph considered alone do not imply response compression.
2. Low algebraic rank and small live separator are genuinely different
   mechanisms: the identity matching has rank `n` and an extensive exposed
   query set.
3. Any theorem compressing a sparse bridge must also use structure of the
   internal landscape class, a sublinear live boundary, contraction, or a
   restricted future family.  Graph sparsity by itself cannot be the state
   variable.

The landscapes in Theorem BQ.2 are arbitrary, not quadratic Boolean forms.
The theorem is therefore a scalable theory falsifier and bridge benchmark,
not an incompressibility theorem for the motivating signing problem.
