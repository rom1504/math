# Benchmark validation: separator, mean-field, and residual states

**Status.** Working synthesis.  The Ising and weighted-automaton derivations
were produced by a solution-hidden agent before it read the repository.  The
profile calculus and fixed-rank comparison are director deductions.  Claims
marked as benchmark interpretation are not literature-novelty claims.

## 1. The common response-profile calculus

Let `S` be a finite declared interface set.  A partial landscape `P` has a
feature `phi_P:Omega_P -> S` and conditional profile

```math
h_P(s)=\max\{H_P(x):\phi_P(x)=s\},                 \tag{PC.1}
```

with `-infinity` for an empty fibre.  A future context contributes an
arbitrary interface potential `g:S -> R` and asks

```math
V_P(g)=\max_{s\in S}\{h_P(s)+g(s)\}.               \tag{PC.2}
```

The standard pinning argument gives

```math
\sup_g|V_P(g)-V_{P'}(g)|=\|h_P-h_{P'}\|_\infty.   \tag{PC.3}
```

Thus the profile is the coarsest exact state and its approximate contextual
complexity is precisely the covering complexity of the *realizable* profile
class in sup norm.  This last qualification matters: arbitrary contexts can
expose every coordinate without the model being able to realize an arbitrary
table of coordinate values.

Suppose interfaces compose through a finite relation
`R subset S_1 times S_2 times S_3` and a feature-level interaction `c(s,t)`.
Then the parent profile is

```math
(h_1\star h_2)(u)
=\max_{(s,t,u)\in R}\{h_1(s)+h_2(t)+c(s,t)\}.       \tag{PC.4}
```

This is an exact feature algebra whenever the relation and interaction are
associative in the evident cocycle sense.  It yields two familiar but sharply
different growth laws.

* **Separator calculus.**  `S` is the set of assignments on the separator;
  gluing imposes equality and eliminates a shared assignment.  A binary
  width-`w` boundary has `2^w` types (or `2^(w-1)` in flip-invariant
  Max-Cut), independently of interior volume.
* **Aggregate calculus.**  `S_n` is the set of reachable sums of fixed-rank
  site features and `R` enforces `u=s+t`.  If the feature alphabet lies in a
  fixed box of `Z^d`, then `|S_n|=O(n^d)`.

The formal composition law alone does not explain compression.  Compression
comes from the growth of the reachable interface types and from the metric
entropy of realizable profiles on them.  The Max-Cut benchmark gives a
well-conditioned exponential profile packing; the fixed-rank mean-field
benchmark has only polynomially many aggregate types to begin with.

## 2. One-dimensional and finite-width Ising

For a width-`w` prefix with exposed boundary spin `s in {+1,-1}^w`, the
solution-hidden derivation produced

```math
f_P(s)=\max_x E_P(x,s).                              \tag{PC.5}
```

Arbitrary boundary fields select each assignment, so (PC.3) proves exact
contextual minimality.  A two-ended fragment has kernel

```math
K_A(s,t)=\max_xE_A(s,x,t),
\qquad
K_{A\circ B}(s,u)=\max_t\{K_A(s,t)+K_B(t,u)\}.      \tag{PC.6}
```

The framework therefore rediscovers the zero-temperature transfer matrix:
not because a boundary table was postulated, but because contextual fibres
are the coarsest quotient and conditional screening supplies max-plus
composition.

### Proposition PC.1 (exact width-one relative recurrence)

For a chain prefix put `a=f(-1)` and `d=f(+1)-f(-1)`.  Appending a spin with
coupling `Jxy` and field `hy` gives

```math
a'=a-h+\max\{d-J,J\},                               \tag{PC.7}
```

```math
d'=2h+\operatorname {sgn}(J)
       \operatorname {clip}(d,-2|J|,2|J|).          \tag{PC.8}
```

Thus `(a,d)` is an exact absolute state and `d` the exact projective state.
An edge transmits the old relative response only up to its coupling capacity;
this is the one-sided version of the directed bottleneck law.

Every real two-ended `2 by 2` kernel has the Walsh expansion

```math
K(x,y)=c+ux+vy+Jxy.                                 \tag{PC.9}
```

It is therefore realized by one general Ising edge with endpoint fields.
For projective kernels of bounded range `B`, the exact small-scale covering
number is

```math
\Theta((B/\epsilon)^3).                             \tag{PC.10}
```

This is a quantitative sanity check: the response theory finds the three
relative transfer parameters and does not retain the exponential path table.

### Proposition PC.2 (sharp width-`w` profile exponent for general pairwise Ising)

For a general pairwise Ising/QUBO fragment of treewidth `w`, every table
`F:{+1,-1}^w -> R` is realizable as a conditional maximum.  Consequently,
after quotienting a common offset, profiles of oscillation at most `B` have
`epsilon`-covering complexity

```math
\Theta\left(2^w\log {B\over\epsilon}\right)\quad\text{bits}. \tag{PC.11}
```

#### Proof

Subtract `c=min_sF(s)` and put `lambda_a=F(a)-c>=0`.  For each assignment
`a`, add an independent Boolean internal variable `y_a` and the pairwise
energy

```math
y_a\lambda_a\left(\sum_i a_i s_i-(w-1)\right).     \tag{PC.12}
```

At `s=a` maximizing `y_a in {0,1}` contributes `lambda_a`; at every other
assignment its coefficient is nonpositive and contributes zero.  The sum of
the gadgets realizes `F`.  Converting `y_a` to an Ising spin is affine and
keeps all terms unary or pairwise.  Processing one `y_a` at a time gives bags
of size `w+1`.

For the upper bound, anchor one coordinate and quantize the remaining
`2^w-1`.  For the lower bound, put those coordinates on a grid separated by
`5 epsilon`.  The projective contextual metric is

```math
\inf_c\|F-F'-c\boldsymbol1\|_\infty
={1\over2}\{\max(F-F')-\min(F-F')\},                \tag{PC.13}
```

so the grid is separated and has the claimed logarithm. `square`

The lower bound is for general width-`w` pairwise interactions, not a strict
bounded-degree lattice strip.  It complements the pure-Max-Cut lower bound:
arbitrary signed pairwise Ising realizes the whole profile cube directly,
whereas nonnegative Max-Cut needs a structured projective-distance gadget
algebra.

## 3. Fixed-rank mean-field as the same calculus with slow type growth

Let the site feature alphabet be a fixed finite subset `F subset Z^d`.  For a
block of `n` sites retain the conditional profile over its total feature

```math
u=\sum_i\phi(s_i).                                  \tag{PC.14}
```

For a fixed bilinear interaction matrix `J`, merging two blocks obeys

```math
h_{m+n}(w)=\max_{u+v=w}
 \{h_m(u)+h_n(v)+u^TJv\}.                           \tag{PC.15}
```

The bilinear term satisfies the cocycle identity

```math
u^TJv+(u+v)^TJz=v^TJz+u^TJ(v+z),                   \tag{PC.16}
```

so composition is bracket-independent.  The reachable `u` lie in a
fixed-dimensional box of side `O(n)`, giving `O(n^d)` profile coordinates.
This recovers total magnetization for Curie--Weiss and the vector of color
counts for fixed-state mean-field Potts models.

The separator and mean-field states are thus instances of the same response
profile calculus, but not of the same complexity regime:

```math
\begin{array}{c|c|c|c}
\text{model}&\text{interface type}&\text{number of types}&\text{closure}\\
\hline
\text{width-}w\text{ separator}&\text{boundary assignment}&q^w
  &\text{eliminate shared label}\\
\text{rank-}d\text{ mean field}&\text{aggregate in }\mathbb Z^d&O(n^d)
  &\text{max-plus convolution plus cocycle}
\end{array}                                         \tag{PC.17}
```

Repeated composition creates no microscopic information in the second row
because all cross-block interactions factor through the aggregate monoid.
If the rank grows with `n`, this explanation disappears; the existing
Shapley--Folkman counterexample shows that the resulting loss can become
macroscopic.

## 4. A harder quotient: weighted-language residuals

Consider a max-plus weighted automaton with transition matrices `T_a`, initial
row `alpha`, and final column `beta`.  For a prefix `x` and suffix `y`, put

```math
p_x=\alpha\odot T_x,
\qquad h_y=T_y\odot\beta,
\qquad L(xy)=\max_i\{p_x(i)+h_y(i)\}.               \tag{PC.18}
```

Contextual equivalence does **not** automatically retain the raw forward
vector.  Its coarsest exact state is the weighted residual

```math
\rho_x:y\longmapsto L(xy),                          \tag{PC.19}
```

or equivalently the quotient

```math
p\sim p'\quad\Longleftrightarrow\quad
\max_i(p_i+h_y(i))=\max_i(p_i'+h_y(i))
\quad\text{for every reachable suffix }y.           \tag{PC.20}
```

Appending a word `u` acts directly on residuals by

```math
(D_u\rho)(y)=\rho(uy),
\qquad D_vD_u=D_{uv}.                               \tag{PC.21}
```

Thus the contextual quotient composes even when it is strictly smaller than
the presented automaton vector.  If arbitrary terminal metrics are declared,
coordinate pinning makes the forward vector minimal and the projective
`epsilon`-covering exponent is `m-1` for an `m`-state bounded box.  For a
fixed suffix language, inaccessible coordinates and identical future
behaviours collapse automatically under (PC.20).

This is the weighted analogue of a Nerode residual, but the quantitative
consequence is response-metric rather than merely finite index.  Real weights
can generate infinitely many exact residuals even for a fixed finite
automaton; a Boolean finite-state conclusion must not be imported.

## 5. Director assessment

The benchmarks support one coherent law:

> A composable extremal state is a conditional response profile on the
> smallest interface through which all declared futures factor.  Its true
> complexity is the metric entropy of the realizable profile semimodule, and
> its growth is governed jointly by the number of reachable interface types
> and the conditioning of gadget-generated profile directions.

This is more informative than “use dynamic programming,” because it predicts
the projective Max-Cut quotient, distinguishes context exposure from profile
realizability, quantifies approximation, and recognizes when a presented
automaton state is nonminimal.  It remains mostly a unification of known
algebras.  The spectral atom-packing theorem and the pure-Max-Cut application
are the genuinely generative additions in this benchmark phase.
