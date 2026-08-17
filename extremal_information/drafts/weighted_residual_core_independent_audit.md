# Independent audit: bounded-delay residual cores

**Verdict: REPAIR (the scoped theorem is correct).**  The row-residual,
cycle-potential, greatest-core, and terminal-pin calculations are valid after
the scope is made explicit.  Two readings must be ruled out in the canonical
statement:

1. the greatest core characterizes a **single support attached to each fixed
   residual context**, not arbitrary anticipatory-support presentations; and
2. failure of that core does not imply a scalar response gap.  The pumpable
   converse comes from a violating **residual-context cycle**, via the exact
   row eigenprofile identity, not from failure of a chosen support
   certificate.

There is also a naming issue.  The relation

```math
T_c(i,j)\geq \kappa(q,c)
```

is an edgewise scalar-threshold relation.  It is generally *not* the tight
argmax relation for the projective residual equation.  Calling it simply
"tight" or "the residual relation" would make a false target-surjectivity
claim look automatic.

The corrected theorem and checks follow.

## 1. Typed max-plus setup

Let `P` be a finite phase set.  A decorated letter `c` has a source phase
`s(c)`, a target phase `t(c)`, and a finite real max-plus block

```math
T_c\in\mathbb R^{I_{s(c)}\times I_{t(c)}}.
```

Words are legal when their phases compose.  Products use the row convention

```math
(p\otimes T)_j=\max_i\{p_i+T(i,j)\}.
```

Assume that there is a fixed `D>=1` such that every legal `D`-letter product
has max-plus row rank one:

```math
T_v(i,j)=A_v(i)+p_v(j),\qquad \max_j p_v(j)=0.                 \tag{WR.1}
```

All coordinates are finite.  Let `Q` be the set of distinct pairs consisting
of the terminal phase and normalized profile `p_v` occurring in (WR.1).
Enabled decorated letters are retained as part of the context data.

### Theorem WR.1 (finite suffix-residual quotient)

For every enabled `c` at `q`, there are unique

```math
\delta(q,c)\in Q,\qquad \kappa(q,c)\in\mathbb R
```

such that

```math
p_q\otimes T_c=\kappa(q,c)\mathbf 1+p_{\delta(q,c)}.          \tag{WR.2}
```

The transition is well defined even if two `D`-words have the same
normalized terminal profile.  Moreover, if

```math
q_0\xrightarrow{c_1}q_1\xrightarrow{c_2}\cdots
\xrightarrow{c_t}q_t=q_0                                    \tag{WR.3}
```

is a legal residual-context cycle and `w=c_1\cdots c_t`, then

```math
\rho(T_w)=\sum_{a=1}^t\kappa(q_{a-1},c_a).                    \tag{WR.4}
```

#### Verification

Write a legal `D`-word as `v=c_1\cdots c_D` and let `v'` be the
length-`D` suffix of `vc`.  Since `T_{v'}` has the form (WR.1), left
multiplication by `T_{c_1}` changes only its row factor.  On the other hand,
using (WR.1) for `T_v` shows that the normalized row of `T_vT_c` is the
normalization of `p_vT_c`.  Therefore it is exactly `p_{v'}`.  If two words
give the same phase and `p_v`, their products with `T_c` plainly have the
same normalization and scalar normalization toll.  This proves (WR.2) and
well-definedness.

Around (WR.3), iteration gives

```math
p_{q_0}\otimes T_w=K\mathbf1+p_{q_0},
\qquad K=\sum_a\kappa(q_{a-1},c_a).                            \tag{WR.5}
```

A finite max-plus left eigenvector forces its eigenvalue to equal spectral
radius.  Indeed, (WR.5) bounds every directed cycle mean above by `K` after
the profile terms telescope.  Choosing, for every target coordinate, one
predecessor attaining the maximum gives a finite tight-predecessor graph;
it contains a directed cycle attaining mean `K`.  This proves (WR.4).

Short prefixes of length below `D` are a finite transient.  They must either
be included explicitly in the quotient or excluded from claims about the
initialized machine.  They do not alter recurrent cycle slopes.

## 2. Exact residual-cycle excess and the potential criterion

Suppose an advertised scalar upper toll satisfies

```math
T_c(i,j)\leq\lambda_c                                          \tag{WR.6}
```

for every entry.  Since `max p_q=max p_delta=0`, (WR.2) implies
`kappa(q,c)<=lambda_c`.  Put

```math
d(q,c)=\lambda_c-\kappa(q,c)\geq0.                            \tag{WR.7}
```

For nonnegative declared budgets `beta_c`, the following are equivalent:

1. there is a potential `psi:Q->R` with

   ```math
   d(q,c)\leq\beta_c+\psi(\delta(q,c))-\psi(q);               \tag{WR.8}
   ```

2. every directed residual-context cycle `C` satisfies

   ```math
   \sum_{(q,c)\in C}(d(q,c)-\beta_c)\leq0.                    \tag{WR.9}
   ```

This is the standard finite difference-constraints criterion.  Summing
(WR.8) also gives, on any context path,

```math
\sum d(q,c)\leq\sum\beta_c+\psi(q_{\rm end})-\psi(q_{\rm start}).
                                                                    \tag{WR.10}
```

The converse is genuinely observable.  If a context cycle `C` violates
(WR.9), then (WR.4) gives, for its word `w`,

```math
\sum_{c\in w}\lambda_c-\rho(T_w)-\sum_{c\in w}\beta_c
=\sum_{(q,c)\in C}(d(q,c)-\beta_c)>0.                         \tag{WR.11}
```

Repeating `w` multiplies (WR.11) by the repetition count.  This is the
correct max-plus spectral-radius pumpable converse.  It depends only on the
row-residual quotient and should be stated separately from the support-core
test below.

If one is optimizing over successor-support choices rather than fixing them,
ordinary difference constraints apply only **after** a policy has been
fixed.  Choosing a policy is a finite mean-payoff-game/strategy problem, not
one undifferentiated linear program.

## 3. The greatest scalar-threshold support core

For an enabled transition `q --c--> q'`, define

```math
R_{q,c}(S)=\{j\in I_{q'}:\exists i\in S,
                         \ T_c(i,j)\geq\kappa(q,c)\}.          \tag{WR.12}
```

Optionally fix allowed endpoint sets `A_q`; otherwise take `A_q=I_q`.  Set

```math
K_q^0=A_q,
```

and iterate

```math
K_q^{n+1}=K_q^n\cap
 \bigcap_{(r,c):\delta(r,c)=q}R_{r,c}(K_r^n).                 \tag{WR.13}
```

The intersection over no incoming transitions is the full target fibre.
Only reachable contexts required by the declared query should be retained.

### Theorem WR.2 (greatest fixed-context support lift)

Iteration (WR.13) terminates after at most `sum_q |A_q|` strict coordinate
deletions.  Its limit `K^infty` is the greatest family, by componentwise
inclusion, satisfying

```math
K_q\subseteq A_q,
\qquad
K_{\delta(q,c)}\subseteq R_{q,c}(K_q)                         \tag{WR.14}
```

for every enabled transition.  In particular, all required components of
`K^infty` are nonempty iff a nonempty **one-support-per-residual-context**
family satisfying (WR.14) exists.

Given such a family, every finite context path has a genuine raw path whose
`c`-edge has weight at least `kappa(q,c)`: choose a terminal point and lift
backwards through (WR.14).  Repetition of a context cycle and finiteness then
close a raw cycle with the corresponding lower weight.

#### Verification and orientation

The map on the right of (WR.13) is monotone.  Any postfixed family `L`
satisfying (WR.14) lies in `K^0`, and induction gives `L subseteq K^n` for
every `n`; hence it lies in the limit.  The limit itself satisfies (WR.14),
which proves greatestness and the equivalence.

The incoming-image orientation in (WR.13) is essential.  Pruning a source
point merely because it does not continue to every target would instead
enforce a rowwise/source-total simulation, which is a different and often
larger resource.

### The threshold relation is not the row-residual argmax relation

Equation (WR.2) has the automatically target-surjective argmax relation

```math
\widehat R_{q,c}
=\{(i,j):p_q(i)+T_c(i,j)
       =\kappa(q,c)+p_{q'}(j)\}.                               \tag{WR.15}
```

But (WR.12) instead asks `T_c(i,j)>=kappa(q,c)`.  Neither relation contains
the other in general, because of the profile difference
`p_{q'}(j)-p_q(i)`.  Thus the support core is not automatically nonempty
from (WR.2).  Its content is precisely the existence of locally
edgewise-`kappa` witnesses, a stronger path presentation than the residual
eigenprofile identity alone.

## 4. A decisive scope counterexample

Let `I={1,2}`, let `p=(0,-1)`, and use two all-finite max-plus generators

```math
T_a=\begin{pmatrix}0&-1\\-2&-3\end{pmatrix},
\qquad
T_b=\begin{pmatrix}-2&-3\\1&0\end{pmatrix}.                  \tag{WR.16}
```

Both matrices have row-rank-one form `A+p`, so the bounded-delay hypothesis
holds with `D=1`.  Directly,

```math
pT_a=pT_b=p,
\qquad \kappa(a)=\kappa(b)=0.                                 \tag{WR.17}
```

Consequently every word has spectral radius zero.  Nevertheless,

```math
R_a=\{(1,1)\},
\qquad R_b=\{(2,1),(2,2)\}.                                   \tag{WR.18}
```

No nonempty `K subseteq I` satisfies both
`K subseteq R_a(K)` and `K subseteq R_b(K)`.  The greatest one-context core
is empty.

Therefore:

```math
\boxed{\text{empty fixed-context support core}
       \not\Rightarrow\text{scalar response loss}.}           \tag{WR.19}
```

It only says that this stricter locally thresholded witness presentation
does not exist.  Allowing several support states over one residual context,
or using the eigenprofile directly, can still answer the scalar query.

This example also prevents the phrase "nonempty iff sharp witnesses exist"
from being read without its quantifiers.  The exact statement is "nonempty
iff a single `Q`-indexed family satisfying (WR.14) exists."

## 5. Terminal pins and projective minimality

For a finite residual profile `p`, declare arbitrary terminal rewards
`z in R^I` and set

```math
V_p(z)=\max_i\{p_i+z_i\}.
```

Then the projective all-terminal response metric is exactly

```math
\inf_{a\in\mathbb R}\sup_{z\in\mathbb R^I}
 |V_p(z)-V_{p'}(z)-a|
=\frac12\operatorname{osc}(p-p').                             \tag{WR.20}
```

Indeed, for every `z`, the difference lies between the minimum and maximum
coordinates of `p-p'`.  A terminal pin which gives one coordinate reward
zero and all other coordinates reward `-L` approaches either endpoint as
`L->infinity`.  Centering the resulting interval gives half its oscillation.

Thus distinct normalized profiles are the minimal exact rooted residual
states when arbitrary terminal pins and history-dependent scalar prefix
tolls are allowed.  Quantitatively, an `epsilon`-accurate merger can contain
only profiles at pairwise projective distance at most `2epsilon`, yielding
the corresponding packing lower bound.

This does **not** establish minimality for the unrooted scalar spectral query.
That query can identify many terminal profiles.  The de Bruijn examples in
Theorems 17.1q, 17.1r, and 18.5 already show the gap between wordwise scalar
response, anticipatory certificate memory, and rooted residual memory.
Bounded terminal fields likewise give only a truncated version of (WR.20).

## 6. Relation to weighted-automaton determinization

The finite family of normalized rows and update (WR.2) are the max-plus
transpose of classical weighted-subset determinization.  This part should be
credited as weighted residual/determinization algebra, not presented as a
new characterization of all finite weighted lumpability.

Mohri's 1997 results give:

- the twins property as a sufficient condition for tropical weighted
  determinization (Theorem 11); and
- for trim unambiguous transducers, determinizable iff twins (Theorem 12).

The bounded-delay row-rank-one hypothesis is a strong, directly checkable
synchronizing sufficient condition.  It is neither necessary for weighted
determinization nor equivalent to twins.  For example, a deterministic
two-state identity automaton is already deterministic but no power of its
max-plus permutation/identity matrix has row rank one.  Even within finite
all-entry matrices,

```math
T=\begin{pmatrix}0&-C\\-C&0\end{pmatrix},\qquad C>0,
```

is max-plus idempotent and never row rank one, while its projective product
semigroup and scalar response are finite.

Recent unary results are also strictly broader in a different direction.
Almagor--Jecker--Mazowiecki--Orlikowski--Purser--Sinclair-Banks prove that
every `n`-state unary tropical WFA is the pointwise minimum of at most `n`
deterministic unary WFAs, each of quadratic size (their Theorem 3).  This is
an anticipatory lower-envelope decomposition, not a single finite row-
residual quotient, and their word-to-final-output semantics is not the cyclic
max-plus spectral response used here.  It shows in particular that the
bounded-delay condition cannot be advertised as necessary even in the unary
setting.

Primary sources:

- Mehryar Mohri, [*Finite-State Transducers in Language and Speech
  Processing*](https://aclanthology.org/J97-2003/), especially Theorems
  11--12.
- Shaull Almagor et al., [*Representing One Letter Weighted Automata Over
  the Tropical Semiring*](https://arxiv.org/abs/2606.26038), especially
  Theorem 3.

## 7. Exact canonical repairs

The result is ready for canonical use after the following edits.

1. State typed phases/fibres, legal words, finiteness of all profile
   coordinates, and the treatment of prefixes shorter than `D`.
2. Call (WR.12) the **scalar `kappa`-threshold relation**.  Display the true
   row-argmax relation (WR.15) so target-surjectivity is not attributed to
   the wrong object.
3. State greatest-core completeness only for one nonempty support per fixed
   residual context (or for the explicitly declared allowed sets).  Do not
   claim completeness among all multi-support anticipatory carriers.
4. Keep the potential/cycle theorem independent of the support core.  A
   violating residual-context cycle pumps by (WR.4); empty-core failure does
   not pump, as (WR.16)--(WR.19) show.
5. State terminal-profile minimality only for the rooted query class with
   arbitrary coordinate pins, modulo history-dependent scalar offsets.
6. Describe `D`-block row rank one as a strong sufficient finite-context
   hypothesis, not an iff form of twins, determinization, or finite weighted
   support lumpability.

With these repairs, the theorem is a useful, checkable middle-level result:
bounded-delay row synchronization produces a finite rooted response quotient;
the greatest core decides whether that quotient also admits a particularly
economical local extremal-witness presentation; and context-cycle potentials
give an exact, observable accumulated-error law.
